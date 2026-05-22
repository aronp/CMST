import os
import io
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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class GWExplorerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("CMST Gravitational Wave Spectrogram Explorer")
        self.root.geometry("1300x850")
        
        # Internal State Management
        self.strain_full = None
        self.dt = None
        self.fs = None
        self.total_duration = 0.0
        self.whitened_strain = None
        self.is_updating_scroll = False
        
        self.create_widgets()
        
    def create_widgets(self):
        # ----------------- Control Sidebar Panel -----------------
        sidebar = ttk.Frame(self.root, padding="10", width=320)
        sidebar.pack(side=tk.LEFT, fill=tk.Y, expand=False)
        sidebar.pack_propagate(False)
        
        # File Source Input
        ttk.Label(sidebar, text="File Input / URL Location:", font=('Helvetica', 10, 'bold')).pack(anchor=tk.W, pady=(0,5))
        self.file_entry = ttk.Entry(sidebar, width=35)
        self.file_entry.pack(fill=tk.X, pady=2)
        self.file_entry.insert(0, "https://gwosc.org/archive/data/O4c1DiscC00_16KHZ_R1/1422917632/V-V1_GWOSC_O4c1DiscC00_16KHZ_R1-1422962688-4096.hdf5")
        
        btn_frame = ttk.Frame(sidebar)
        btn_frame.pack(fill=tk.X, pady=5)
        ttk.Button(btn_frame, text="Browse Local...", command=self.browse_file).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        ttk.Button(btn_frame, text="Load / Download", command=self.start_load_thread).pack(side=tk.RIGHT, expand=True, fill=tk.X, padx=2)
        
        self.status_label = ttk.Label(sidebar, text="Status: Ready", foreground="blue", wraplength=280)
        self.status_label.pack(anchor=tk.W, pady=5)
        
        ttk.Separator(sidebar, orient='horizontal').pack(fill=tk.X, pady=10)
        
        # Processing Window Parameters
        ttk.Label(sidebar, text="Processing Specs:", font=('Helvetica', 10, 'bold')).pack(anchor=tk.W, pady=(0,5))
        
        ttk.Label(sidebar, text="NPERSEG (Window Stiffness Size):").pack(anchor=tk.W)
        self.nperseg_entry = ttk.Entry(sidebar)
        self.nperseg_entry.pack(fill=tk.X, pady=2)
        self.nperseg_entry.insert(0, "256")

        ttk.Label(sidebar, text="NFFT (Interpolation Padded Size):").pack(anchor=tk.W)
        self.nfft_entry = ttk.Entry(sidebar)
        self.nfft_entry.pack(fill=tk.X, pady=2)
        self.nfft_entry.insert(0, "1024")
        
        ttk.Label(sidebar, text="Overlap Target (Percentage 0-99%):").pack(anchor=tk.W)
        self.overlap_pct_entry = ttk.Entry(sidebar)
        self.overlap_pct_entry.pack(fill=tk.X, pady=2)
        self.overlap_pct_entry.insert(0, "85")
        
        ttk.Separator(sidebar, orient='horizontal').pack(fill=tk.X, pady=10)
        
        # View Axis Frame Controls
        ttk.Label(sidebar, text="View Windows Control:", font=('Helvetica', 10, 'bold')).pack(anchor=tk.W, pady=(0,5))
        
        ttk.Label(sidebar, text="Time Window Width (seconds):").pack(anchor=tk.W)
        self.t_width_entry = ttk.Entry(sidebar)
        self.t_width_entry.pack(fill=tk.X, pady=2)
        self.t_width_entry.insert(0, "0.65")

        ttk.Label(sidebar, text="Time Window Center (seconds):").pack(anchor=tk.W)
        self.t_center_entry = ttk.Entry(sidebar)
        self.t_center_entry.pack(fill=tk.X, pady=2)
        self.t_center_entry.insert(0, "0.00") # Dynamically overwritten upon successful loading
        
        ttk.Label(sidebar, text="Min Frequency (Hz):").pack(anchor=tk.W)
        self.f_min_entry = ttk.Entry(sidebar)
        self.f_min_entry.pack(fill=tk.X, pady=2)
        self.f_min_entry.insert(0, "10")
        
        ttk.Label(sidebar, text="Max Frequency (Hz):").pack(anchor=tk.W)
        self.f_max_entry = ttk.Entry(sidebar)
        self.f_max_entry.pack(fill=tk.X, pady=2)
        self.f_max_entry.insert(0, "500")
        
        ttk.Button(sidebar, text="Update Plot Window", command=self.update_plot).pack(fill=tk.X, pady=15)
        
        # ----------------- Visual Canvas Panel -----------------
        self.plot_frame = ttk.Frame(self.root, padding="5")
        self.plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        self.fig, self.ax = plt.subplots(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)
        
        # Native Time Navigation Scrollbar Integration
        self.time_scrollbar = tk.Scrollbar(self.plot_frame, orient=tk.HORIZONTAL, command=self.on_scroll)
        self.time_scrollbar.pack(fill=tk.X, side=tk.BOTTOM)
        
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        
    def browse_file(self):
        filename = filedialog.askopenfilename(filetypes=[("HDF5 files", "*.hdf5 *.h5")])
        if filename:
            self.file_entry.delete(0, tk.END)
            self.file_entry.insert(0, filename)
            
    def start_load_thread(self):
        target = self.file_entry.get().strip()
        if not target:
            messagebox.showerror("Error", "Please provide a valid file location path or URL.")
            return
        self.status_label.config(text="Status: Loading file asset...", foreground="orange")
        threading.Thread(target=self.load_data, args=(target,), daemon=True).start()
        
    def load_data(self, target):
        try:
            if target.startswith("http://") or target.startswith("https://"):
                local_filename = target.split("/")[-1]
                if not os.path.exists(local_filename):
                    self.root.after(0, lambda: self.status_label.config(text="Status: Downloading stream target...", foreground="orange"))
                    r_file = requests.get(target)
                    with open(local_filename, 'wb') as f_out:
                        f_out.write(r_file.content)
                target = local_filename
            
            self.root.after(0, lambda: self.status_label.config(text="Status: Reading HDF5 arrays...", foreground="orange"))
            
            with h5py.File(target, 'r') as f:
                if 'strain/Strain' in f:
                    self.strain_full = f['strain/Strain'][:]
                    self.dt = f['strain/Strain'].attrs['Xspacing']
                elif 'strain/strain' in f:
                    self.strain_full = f['strain/strain'][:]
                    self.dt = f['strain/strain'].attrs['Xspacing']
                else:
                    possible_keys = list(f.keys())
                    raise KeyError(f"Could not find standard strain path. Available top keys are: {possible_keys}")
                
            self.fs = 1.0 / self.dt
            self.total_duration = len(self.strain_full) * self.dt
            
            self.root.after(0, lambda: self.status_label.config(text="Status: Whitening whole sample...", foreground="orange"))
            self.whitened_strain = self.whiten_whole_sample(self.strain_full, self.dt, self.cmst_window_func)
            
            self.root.after(0, self.signal_load_success)
            
        except Exception as e:
            self.root.after(0, lambda err=e: self.signal_load_failure(err))

    def signal_load_success(self):
        self.status_label.config(text=f"Status: Loaded. Size: {self.total_duration:.1f}s", foreground="green")
        
        # CHANGED: Calculate exact halfway mark configuration of file
        midpoint_time = self.total_duration / 2.0
        
        # Automatically insert the safe midpoint calculation value into user control text entry box
        self.t_center_entry.delete(0, tk.END)
        self.t_center_entry.insert(0, f"{midpoint_time:.2f}")
        
        self.update_plot()

    def signal_load_failure(self, err):
        self.status_label.config(text="Status: Processing crash.", foreground="red")
        messagebox.showerror("Background Thread Error", f"The loader failed!\n\nError type: {type(err).__name__}\nDetails: {str(err)}")
        
    @staticmethod
    def cmst_window_func(N):
        t = np.linspace(-1, 1, N)
        with np.errstate(divide='ignore', invalid='ignore'):
            w = np.exp(t**4 / (t**2 - 1))
        return np.where(np.abs(t) < 1, w, 0.0)

    @staticmethod
    def calculate_alpha(window):
        sum_w = np.sum(window)
        t = np.linspace(-1, 1, len(window))
        return np.sum(t*t*window)/sum_w

    def whiten_whole_sample(self, strain, dt, window_func):
        N = len(strain)
        spec = np.fft.rfft(strain * window_func(N))
        freqs = np.fft.rfftfreq(N, dt)
        mag = np.abs(spec)
        psd_smooth = interp1d(freqs, mag, kind='nearest', fill_value="extrapolate")(freqs)
        whitened_spec = spec / (psd_smooth + 1e-20)
        return np.fft.irfft(whitened_spec, n=N)

    def on_scroll(self, *args):
        if self.whitened_strain is None or self.is_updating_scroll:
            return
        
        try:
            t_center = float(self.t_center_entry.get())
            t_width = float(self.t_width_entry.get())
            
            if args[0] == 'moveto':
                frac = float(args[1])
                new_center = (frac * self.total_duration) + (t_width / 2.0)
            elif args[0] == 'scroll':
                steps = int(args[1])
                unit = args[2]
                
                if unit == 'pages' or unit == 'units':
                    new_center = t_center + (steps * 0.5 * t_width)
                else:
                    return
            else:
                return

            # CHANGED: Tighten boundaries to prevent scrolling out past edge parameters
            new_center = max(t_width / 2.0, min(new_center, self.total_duration - t_width / 2.0))
            
            self.t_center_entry.delete(0, tk.END)
            self.t_center_entry.insert(0, f"{new_center:.3f}")
            
            self.update_plot(update_scrollbar=True)
            
        except Exception as e:
            pass

    def update_plot(self, update_scrollbar=True):
        if self.whitened_strain is None:
            self.ax.clear()
            self.ax.text(0.5, 0.5, "Load data asset to initialize visualization.", ha='center', va='center')
            self.canvas.draw()
            return
            
        try:
            nperseg = int(self.nperseg_entry.get())
            nfft = int(self.nfft_entry.get())
            overlap_pct = float(self.overlap_pct_entry.get())
            
            if overlap_pct >= 100.0:
                overlap_pct = 99.0
                self.overlap_pct_entry.delete(0, tk.END)
                self.overlap_pct_entry.insert(0, "99")
            elif overlap_pct < 0.0:
                overlap_pct = 0.0
                self.overlap_pct_entry.delete(0, tk.END)
                self.overlap_pct_entry.insert(0, "0")
                
            noverlap = int(np.floor(nperseg * (overlap_pct / 100.0)))
            if noverlap >= nperseg:
                noverlap = nperseg - 1
                
            t_center = float(self.t_center_entry.get())
            t_width = float(self.t_width_entry.get())
            
            # CHANGED: Force dynamic input cleaning to confirm requested value falls inside physical bounds of file duration
            if t_center < 0.0 or t_center > self.total_duration:
                t_center = self.total_duration / 2.0
                self.t_center_entry.delete(0, tk.END)
                self.t_center_entry.insert(0, f"{t_center:.2f}")

            t_start = max(0.0, t_center - t_width / 2.0)
            t_end = min(self.total_duration, t_center + t_width / 2.0)
            
            f_min = float(self.f_min_entry.get())
            f_max = float(self.f_max_entry.get())
            
            i_start = int(t_start * self.fs)
            i_end = int(t_end * self.fs)
            slice_strain = self.whitened_strain[i_start:i_end]
            
            if len(slice_strain) < nperseg:
                return 
                
            win_seg = self.cmst_window_func(nperseg)
            alpha = self.calculate_alpha(win_seg)
            
            f, t_spec, Sxx_power = spectrogram(slice_strain, self.fs,
                                               window=win_seg,
                                               nperseg=nperseg,
                                               nfft=nfft,
                                               noverlap=noverlap,
                                               scaling='spectrum')
            
            Sxx = np.sqrt(Sxx_power)
            
            laplacian_f = np.zeros_like(Sxx)
            laplacian_f[1:-1, :] = Sxx[2:, :] - 2*Sxx[1:-1, :] + Sxx[0:-2, :]
            Sxx_sharp = Sxx - (alpha / (nfft / nperseg)) * laplacian_f
            
            Sxx_sharp = np.nan_to_num(Sxx_sharp, nan=1e-20, posinf=1e-20, neginf=1e-20)
            Sxx_sharp = np.maximum(Sxx_sharp, 1e-20)
            
            t_spec_global = t_spec + t_start
            
            self.ax.clear()
            
            low_p = 10
            vmax = np.percentile(Sxx_sharp, 99.99)
            vmin = np.percentile(Sxx_sharp, low_p)
            
            if np.isnan(vmin) or np.isnan(vmax) or vmin >= vmax:
                vmin, vmax = 1e-20, 1e-5
                
            levels = np.linspace(vmin, vmax, 30)
            
            self.ax.pcolormesh(t_spec_global, f, Sxx_sharp, shading='gouraud', cmap='viridis')
            contour_filled = self.ax.contourf(t_spec_global, f, Sxx_sharp, levels=levels, cmap='inferno', extend='both')
            self.ax.contour(t_spec_global, f, Sxx_sharp, levels=levels, colors='white', linewidths=0.5, alpha=0.3)
            
            self.ax.set_xlim(t_start, t_end)
            self.ax.set_ylim(f_min, f_max)
            
            self.ax.set_title("Interferometer Data Analysis: CMST-Sharpened Matrix")
            self.ax.set_ylabel("Frequency (Hz)")
            self.ax.set_xlabel("Time (s)")
            self.ax.grid(True, alpha=0.2, color='white')
            
            if update_scrollbar:
                self.is_updating_scroll = True
                left_frac = t_start / self.total_duration
                right_frac = t_end / self.total_duration
                self.time_scrollbar.set(left_frac, right_frac)
                self.is_updating_scroll = False

            self.canvas.draw()
            
        except Exception as e:
            messagebox.showerror("Plot Error", f"Failed to compute window frame: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = GWExplorerApp(root)
    root.mainloop()
