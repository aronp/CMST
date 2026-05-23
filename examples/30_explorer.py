import os
import sys
import threading
import numpy as np
import h5py
import requests
from scipy.signal import spectrogram
from scipy.interpolate import interp1d

import tkinter as tk
from tkinter import messagebox, filedialog, ttk

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class GWExplorerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("ulti-Detector CMST Transient Explorer")
        self.root.geometry("1400x900")
        
        # Core Global State Engine
        self.detectors = {
            'H1': {'raw': None, 'whitened': None, 'Sxx': None, 't': None, 'f': None, 'loaded': False, 'active_corr': tk.BooleanVar(value=True)},
            'L1': {'raw': None, 'whitened': None, 'Sxx': None, 't': None, 'f': None, 'loaded': False, 'active_corr': tk.BooleanVar(value=True)},
            'V1': {'raw': None, 'whitened': None, 'Sxx': None, 't': None, 'f': None, 'loaded': False, 'active_corr': tk.BooleanVar(value=True)}
        }
        self.dt = None
        self.fs = None
        self.total_duration = 0.0
        self.is_dragging = False
        self.drag_start_x = None
        
        # Build Menus & Layout Frames
        self.create_menu()
        self.create_widgets()
        
        # OS window close event
        self.root.protocol("WM_DELETE_WINDOW", self.safe_exit)
        
        # Terminal termination sequence (Ctrl+C)
        self.root.bind_all("<Control-c>", lambda event: self.safe_exit())
        
    def safe_exit(self):
        try:
            plt.close('all')
        except Exception:
            pass
        self.root.quit()
        self.root.destroy()
        sys.exit(0)

    def create_menu(self):
        menubar = tk.Menu(self.root)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Browse & Download from GWOSC...", command=self.open_gwosc_browser)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.safe_exit)
        menubar.add_cascade(label="File", menu=file_menu)
        
        param_menu = tk.Menu(menubar, tearoff=0)
        param_menu.add_command(label="Whiten Segmentation Config...", command=self.open_whiten_dialog)
        param_menu.add_command(label="FFT/Stiffness Configuration...", command=self.open_fft_dialog)
        param_menu.add_command(label="Display Limits & Window Width...", command=self.open_display_dialog)
        menubar.add_cascade(label="Parameters", menu=param_menu)
        
        self.root.config(menu=menubar)

    def create_widgets(self):
        self.whiten_interval = tk.DoubleVar(value=5.0)
        self.nperseg = tk.IntVar(value=256)
        self.nfft = tk.IntVar(value=1024)
        self.overlap_pct = tk.DoubleVar(value=85.0)
        
        self.t_center = tk.DoubleVar(value=0.0)
        self.t_width_seconds = 0.65  # Internal float variable tracking evaluated width seconds
        self.t_width_str = tk.StringVar(value="0.65")  # CHANGED: Entry-bound tracking string for decimals or fractions
        self.f_min = tk.DoubleVar(value=0)
        self.f_max = tk.DoubleVar(value=500)
        self.pct_low = tk.DoubleVar(value=60.0)
        self.pct_high = tk.DoubleVar(value=99.99)
        
        sidebar = ttk.Frame(self.root, padding="10", width=300)
        sidebar.pack(side=tk.LEFT, fill=tk.Y)
        sidebar.pack_propagate(False)
        
        ttk.Label(sidebar, text="Active Detector Loading:", font=('Helvetica', 10, 'bold')).pack(anchor=tk.W, pady=5)
        for det in ['H1', 'L1', 'V1']:
            f = ttk.LabelFrame(sidebar, text=f"Interferometer {det}", padding=5)
            f.pack(fill=tk.X, pady=4)
            ttk.Button(f, text="Load Local HDF5...", command=lambda d=det: self.load_local_detector(d)).pack(fill=tk.X)
            self.detectors[det]['status_lbl'] = ttk.Label(f, text="Status: Empty", foreground="gray")
            self.detectors[det]['status_lbl'].pack(anchor=tk.W)
            
        ttk.Separator(sidebar, orient='horizontal').pack(fill=tk.X, pady=10)
        
        ttk.Label(sidebar, text="Correlation Mask Switches:", font=('Helvetica', 10, 'bold')).pack(anchor=tk.W, pady=5)
        for det in ['H1', 'L1', 'V1']:
            ttk.Checkbutton(sidebar, text=f"Include {det} in Joint Product", variable=self.detectors[det]['active_corr'], command=self.update_all_tabs).pack(anchor=tk.W, pady=2)
            
        self.global_status = ttk.Label(sidebar, text="System: Idle", font=('Helvetica', 9, 'italic'), foreground="blue")
        self.global_status.pack(side=tk.BOTTOM, fill=tk.X, pady=10)

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.tabs = {}
        for tab_name in ['H1', 'L1', 'V1', 'Joint Correlation']:
            frame = ttk.Frame(self.notebook)
            self.notebook.add(frame, text=tab_name)
            
            container = ttk.Frame(frame)
            container.pack(fill=tk.BOTH, expand=True)
            
            pbar = ttk.Progressbar(container, orient=tk.HORIZONTAL, mode='determinate')
            pbar.pack(fill=tk.X, padx=20, pady=10)
            pbar.pack_forget()
            
            fig, ax = plt.subplots(figsize=(8, 6))
            canvas = FigureCanvasTkAgg(fig, master=container)
            canvas_widget = canvas.get_tk_widget()
            canvas_widget.pack(fill=tk.BOTH, expand=True)
            
            canvas_widget.bind("<ButtonPress-1>", self.start_drag)
            canvas_widget.bind("<B1-Motion>", self.drag_motion)
            canvas_widget.bind("<ButtonRelease-1>", self.end_drag)
            
            self.tabs[tab_name] = {'fig': fig, 'ax': ax, 'canvas': canvas, 'pbar': pbar}
            
        self.notebook.bind("<<NotebookTabChanged>>", lambda e: self.update_all_tabs())

    def open_gwosc_browser(self):
        win = tk.Toplevel(self.root)
        win.title("GWOSC API Catalog Browser")
        win.geometry("500x150")
        ttk.Label(win, text="Paste Target REST Endpoint or Dataset Handle URL:").pack(pady=5)
        ent = ttk.Entry(win, width=60)
        ent.pack(pady=5)
        ent.insert(0, "https://gwosc.org/archive/data/O4c1DiscC00_16KHZ_R1/1422917632/V-V1_GWOSC_O4c1DiscC00_16KHZ_R1-1422962688-4096.hdf5")
        
        def download():
            url = ent.get().strip()
            win.destroy()
            det_key = 'V1' if 'V-V1' in url else ('H1' if 'H-H1' in url else 'L1')
            self.global_status.config(text=f"Downloading asset to thread layer...")
            threading.Thread(target=self.download_and_ingest, args=(url, det_key), daemon=True).start()
            
        ttk.Button(win, text="Fetch & Whitening Process", command=download).pack(pady=10)

    def open_whiten_dialog(self):
        win = tk.Toplevel(self.root)
        win.title("Segmentation Window Properties")
        ttk.Label(win, text="Rolling Interval Chunk Length (s):").grid(row=0, column=0, padx=10, pady=10)
        ttk.Entry(win, textvariable=self.whiten_interval).grid(row=0, column=1, padx=10, pady=10)
        ttk.Button(win, text="Save Parameters", command=win.destroy).grid(row=1, columnspan=2, pady=10)

    def open_fft_dialog(self):
        win = tk.Toplevel(self.root)
        win.title("FFT Properties")
        ttk.Label(win, text="NPERSEG Size:").grid(row=0, column=0, padx=10, pady=5)
        ttk.Entry(win, textvariable=self.nperseg).grid(row=0, column=1, padx=10, pady=5)
        ttk.Label(win, text="NFFT (Padding):").grid(row=1, column=0, padx=10, pady=5)
        ttk.Entry(win, textvariable=self.nfft).grid(row=1, column=1, padx=10, pady=5)
        ttk.Label(win, text="Overlap Target (%):").grid(row=2, column=0, padx=10, pady=5)
        ttk.Entry(win, textvariable=self.overlap_pct).grid(row=2, column=1, padx=10, pady=5)
        ttk.Button(win, text="Apply Changes", command=win.destroy).grid(row=3, columnspan=2, pady=10)

    def open_display_dialog(self):
        win = tk.Toplevel(self.root)
        win.title("Display Scale & View Width")
        
        # CHANGED: Added input row fields for the window width to accept decimals or fractions string values
        ttk.Label(win, text="Window Width (s or fraction e.g., 1/4):").grid(row=0, column=0, padx=10, pady=5, sticky=tk.W)
        ttk.Entry(win, textvariable=self.t_width_str).grid(row=0, column=1, padx=10, pady=5)
        
        ttk.Label(win, text="Lower Floor Percentile (%):").grid(row=1, column=0, padx=10, pady=5, sticky=tk.W)
        ttk.Entry(win, textvariable=self.pct_low).grid(row=1, column=1, padx=10, pady=5)
        
        ttk.Label(win, text="Upper Peak Percentile (%):").grid(row=2, column=0, padx=10, pady=5, sticky=tk.W)
        ttk.Entry(win, textvariable=self.pct_high).grid(row=2, column=1, padx=10, pady=5)
        
        def save_and_eval():
            # CHANGED: Custom robust parser inside save loop to digest raw fraction entries safely
            raw_str = self.t_width_str.get().strip()
            try:
                if "/" in raw_str:
                    num, denom = map(float, raw_str.split("/"))
                    val = num / denom
                else:
                    val = float(raw_str)
                
                if val <= 0: raise ValueError()
                self.t_width_seconds = val
                win.destroy()
                self.update_all_tabs()
            except Exception:
                messagebox.showerror("Parsing Error", "Invalid Width entry. Please input a positive float number or explicit fraction pattern like '1/4'.")
                
        ttk.Button(win, text="Apply Changes", command=save_and_eval).grid(row=3, columnspan=2, pady=10)

    def start_drag(self, event):
        self.is_dragging = True
        self.drag_start_x = event.x

    def drag_motion(self, event):
        if not self.is_dragging or self.total_duration == 0:
            return
        dx = event.x - self.drag_start_x
        time_shift = -(dx / 100.0) * self.t_width_seconds
        new_center = self.t_center.get() + time_shift
        new_center = max(self.t_width_seconds/2, min(new_center, self.total_duration - self.t_width_seconds/2))
        self.t_center.set(new_center)
        self.drag_start_x = event.x
        self.update_all_tabs()

    def end_drag(self, event):
        self.is_dragging = False
        self.update_all_tabs()

    def load_local_detector(self, det):
        filename = filedialog.askopenfilename(filetypes=[("HDF5 keys", "*.hdf5 *.h5")])
        if filename:
            self.global_status.config(text=f"Ingesting local array for {det}...")
            threading.Thread(target=self.process_pipeline_worker, args=(filename, det), daemon=True).start()

    def download_and_ingest(self, url, det):
        try:
            local_name = url.split("/")[-1]
            if not os.path.exists(local_name):
                r = requests.get(url)
                with open(local_name, 'wb') as f:
                    f.write(r.content)
            self.process_pipeline_worker(local_name, det)
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Network Crash", str(e)))

    def show_pbar(self, det):
        self.tabs[det]['pbar'].pack(fill=tk.X, padx=20, pady=10, before=self.tabs[det]['canvas'].get_tk_widget())
        self.tabs[det]['pbar']['value'] = 0

    def update_pbar(self, det, val):
        self.tabs[det]['pbar']['value'] = val

    def hide_pbar(self, det):
        self.tabs[det]['pbar'].pack_forget()

    def process_pipeline_worker(self, filepath, det):
        try:
            self.root.after(0, lambda: self.notebook.select(self.notebook.tabs()[['H1', 'L1', 'V1'].index(det)]))
            self.root.after(0, lambda: self.show_pbar(det))
            
            with h5py.File(filepath, 'r') as f:
                key = 'strain/Strain' if 'strain/Strain' in f else 'strain/strain'
                self.detectors[det]['raw'] = f[key][:]
                self.dt = f[key].attrs['Xspacing']
                
                # NEW DETECTOR STATE: Read the true absolute GPS anchor attribute from the file meta
                if 'meta/GPSstart' in f:
                    self.detectors[det]['gps_start'] = f['meta/GPSstart'][()]
                elif 'meta' in f and 'GPSstart' in f['meta'].attrs:
                    self.detectors[det]['gps_start'] = f['meta'].attrs['GPSstart']
                else:
                    self.detectors[det]['gps_start'] = 0.0 # Fallback
            
            self.fs = 1.0 / self.dt
            self.total_duration = len(self.detectors[det]['raw']) * self.dt
            
            self.detectors[det]['whitened'] = self.whiten_by_intervals(self.detectors[det]['raw'], self.fs, self.whiten_interval.get(), det)
            self.compute_complete_spectrogram(det)
            
            self.detectors[det]['loaded'] = True
            
            if self.t_center.get() == 0.0:
                self.root.after(0, lambda: self.t_center.set(self.total_duration / 2.0))
                
            self.root.after(0, lambda: self.detectors[det]['status_lbl'].config(text="Status: Loaded & Cached", foreground="green"))
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, self.update_all_tabs)
        except Exception as e:
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, lambda: messagebox.showerror("Pipeline Failure", str(e)))

    def render_joint_correlation(self, t_start, t_end):
        ax = self.tabs['Joint Correlation']['ax']
        canvas = self.tabs['Joint Correlation']['canvas']
        ax.clear()
        
        active_matrices = [det for det in ['H1', 'L1', 'V1'] if self.detectors[det]['loaded'] and self.detectors[det]['active_corr'].get()]
        if len(active_matrices) < 2:
            ax.text(0.5, 0.5, "Load at least TWO detectors to run true cross-correlation.", ha='center', va='center')
            canvas.draw()
            return
            
        # Use the first active detector as our absolute master clock mapping reference
        ref = active_matrices[0]
        ref_gps = self.detectors[ref]['gps_start']
        
        t_ref = self.detectors[ref]['t']
        f_ref = self.detectors[ref]['f']
        
        t_mask_ref = (t_ref >= t_start) & (t_ref <= t_end)
        f_mask_ref = (f_ref >= self.f_min.get()) & (f_ref <= self.f_max.get())
        
        if not np.any(t_mask_ref) or not np.any(f_mask_ref): return
        
        # Initialize our clean joint product matrix map
        joint_product = np.ones((np.sum(f_mask_ref), np.sum(t_mask_ref)))
        
        # Loop through detectors and interpolate their data onto our master clock mapping reference
        for det in active_matrices:
            det_gps = self.detectors[det]['gps_start']
            time_offset = ref_gps - det_gps # Calculate relative physical displacement window
            
            # Translate current display boundaries into the secondary detector's personal array timeline
            t_det_start = t_start + time_offset
            t_det_end = t_end + time_offset
            
            t_arr = self.detectors[det]['t']
            f_arr = self.detectors[det]['f']
            
            t_mask_det = (t_arr >= t_det_start) & (t_arr <= t_det_end)
            f_mask_det = (f_arr >= self.f_min.get()) & (f_arr <= self.f_max.get())
            
            if not np.any(t_mask_det): continue
            
            # Extract and normalize the sub-slice matrix frame cleanly
            S_sub = self.detectors[det]['Sxx'][np.ix_(f_mask_det, t_mask_det)]
            S_norm = S_sub / np.max(S_sub)
            
            # If the matrix sizing doesn't line up perfectly due to floating sample points, resample it
            if S_norm.shape != joint_product.shape:
                from scipy.interpolate import interp2d
                t_det_actual = t_arr[t_mask_det] - time_offset
                f_det_actual = f_arr[f_mask_det]
                
                # Map onto reference frame grid layout
                interp_func = interp2d(t_det_actual, f_det_actual, S_norm, kind='linear')
                S_norm = interp_func(t_ref[t_mask_ref], f_ref[f_mask_ref])
                
            joint_product *= S_norm
            
        vmin, vmax = np.percentile(joint_product, self.pct_low.get()), np.percentile(joint_product, self.pct_high.get())
        levels = np.linspace(vmin, max(vmax, vmin + 1e-10), 25)
        
        ax.pcolormesh(t_ref[t_mask_ref], f_ref[f_mask_ref], joint_product, shading='gouraud', cmap='magma')
        ax.contourf(t_ref[t_mask_ref], f_ref[f_mask_ref], joint_product, levels=levels, cmap='hot', extend='both')
        ax.set_xlim(t_start, t_end)
        ax.set_ylim(self.f_min.get(), self.f_max.get())
        canvas.draw_idle()
        
    def whiten_by_intervals(self, strain, fs, interval_sec, det):
        N = len(strain)
        chunk_size = int(interval_sec * fs)
        whitened = np.zeros_like(strain)
        t = np.linspace(-1, 1, chunk_size)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            win = np.where(np.abs(t) < 1, np.exp(t**4 / (t**2 - 1)), 0.0)
            
        starts = list(range(0, N, chunk_size // 2))
        total_steps = len(starts)
        
        for idx, start in enumerate(starts):
            end = start + chunk_size
            if end > N: break
            slice_data = strain[start:end] * win
            spec = np.fft.rfft(slice_data)
            mag = np.abs(spec)
            freqs = np.fft.rfftfreq(chunk_size, 1/fs)
            psd_smooth = interp1d(freqs, mag, kind='nearest', fill_value="extrapolate")(freqs)
            
            with np.errstate(divide='ignore', invalid='ignore'):
                whitened_spec = spec / (psd_smooth + 1e-20)
                
            whitened[start:end] += np.fft.irfft(whitened_spec, n=chunk_size)
            
            pct = int((idx / total_steps) * 100)
            if idx % 5 == 0:
                self.root.after(0, lambda p=pct: self.update_pbar(det, p))
                
        return whitened

    def compute_complete_spectrogram(self, det):
        nperseg = self.nperseg.get()
        nfft = self.nfft.get()
        noverlap = int(nperseg * (self.overlap_pct.get() / 100.0))
        
        t_win = np.linspace(-1, 1, nperseg)
        with np.errstate(divide='ignore', invalid='ignore'):
            win_seg = np.where(np.abs(t_win) < 1, np.exp(t_win**4 / (t_win**2 - 1)), 0.0)
        
        f, t, Sxx_power = spectrogram(self.detectors[det]['whitened'], self.fs,
                                      window=win_seg, nperseg=nperseg, nfft=nfft, noverlap=noverlap, scaling='spectrum')
        
        Sxx = np.sqrt(Sxx_power)
        alpha = np.sum(t_win*t_win*win_seg)/np.sum(win_seg)
        laplacian_f = np.zeros_like(Sxx)
        laplacian_f[1:-1, :] = Sxx[2:, :] - 2*Sxx[1:-1, :] + Sxx[0:-2, :]
        Sxx_sharp = Sxx - (alpha / (nfft / nperseg)) * laplacian_f
        
        self.detectors[det]['Sxx'] = np.maximum(np.nan_to_num(Sxx_sharp, nan=1e-20), 1e-20)
        self.detectors[det]['t'] = t
        self.detectors[det]['f'] = f

    def update_all_tabs(self):
        if self.total_duration == 0:
            return
        t_center = self.t_center.get()
        # CHANGED: Replaced self.t_width.get() references to use self.t_width_seconds evaluation engine
        t_width = self.t_width_seconds
        t_start, t_end = max(0.0, t_center - t_width/2), min(self.total_duration, t_center + t_width/2)
        
        for det in ['H1', 'L1', 'V1']:
            if self.detectors[det]['loaded']:
                self.render_canvas_frame(self.tabs[det]['ax'], self.tabs[det]['canvas'], self.detectors[det]['Sxx'], self.detectors[det]['t'], self.detectors[det]['f'], t_start, t_end)
                
        self.render_joint_correlation(t_start, t_end)

    def render_canvas_frame(self, ax, canvas, Sxx, t, f, t_start, t_end):
        ax.clear()
        t_mask = (t >= t_start) & (t <= t_end)
        f_mask = (f >= self.f_min.get()) & (f <= self.f_max.get())
        
        if not np.any(t_mask) or not np.any(f_mask): return
        
        Sxx_sub = Sxx[np.ix_(f_mask, t_mask)]
        t_sub, f_sub = t[t_mask], f[f_mask]
        
        vmin = np.percentile(Sxx_sub, self.pct_low.get())
        vmax = np.percentile(Sxx_sub, self.pct_high.get())
        levels = np.linspace(vmin, max(vmax, vmin + 1e-100), 25)
        
        ax.pcolormesh(t_sub, f_sub, Sxx_sub, shading='gouraud', cmap='viridis')
        ax.contourf(t_sub, f_sub, Sxx_sub, levels=levels, cmap='inferno', extend='both')
        ax.set_xlim(t_start, t_end)
        ax.set_ylim(self.f_min.get(), self.f_max.get())
        canvas.draw_idle()


if __name__ == "__main__":
    root = tk.Tk()
    app = GWExplorerApp(root)
    root.mainloop()
