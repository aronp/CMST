import os
import sys
import threading
import numpy as np
import h5py
import requests
from scipy import signal
from scipy.signal import spectrogram
from scipy.interpolate import interp1d

import tkinter as tk
from tkinter import messagebox, filedialog, ttk

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from gwosc.locate import get_event_urls

class GWExplorerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Multi-Detector CMST Transient Explorer")
        self.root.geometry("1450x950")

        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            
            if sys.platform.startswith('win'):
                ico_path = os.path.join(script_dir, 'explorer_icon.ico')
                if os.path.exists(ico_path):
                    self.root.iconbitmap(ico_path)
            else:
                png_path = os.path.join(script_dir, 'explorer_icon.png')
                if os.path.exists(png_path):
                    self.icon_img = tk.PhotoImage(file=png_path)
                    self.root.iconphoto(True, self.icon_img)
                else:
                    print(f"Could not find icon at: {png_path}")
        except Exception as e:
            print(f"Icon missing or failed to load: {e}")
            
        # ---------------------------------------        
        # Core Global State Engine
        self.detectors = {
            'H1': {'raw': None, 'whitened': None, 'Sxx': None, 't': None, 'f': None, 'loaded': False, 'active_corr': tk.BooleanVar(value=True)},
            'L1': {'raw': None, 'whitened': None, 'Sxx': None, 't': None, 'f': None, 'loaded': False, 'active_corr': tk.BooleanVar(value=True)},
            'V1': {'raw': None, 'whitened': None, 'Sxx': None, 't': None, 'f': None, 'loaded': False, 'active_corr': tk.BooleanVar(value=False)}
        }
        self.dt = None
        self.fs = None
        self.total_duration = 0.0
        self.is_dragging = False
        self.drag_start_x = None
        
        # Active file metadata cache trackers
        self.current_file_duration = "N/A"
        self.current_file_rate = "N/A"
        self.current_event_offset = "N/A"
        self.current_filenames = {'H1': 'N/A', 'L1': 'N/A', 'V1': 'N/A'}
        
        # Playback Control Variables
        self.is_playing = False
        self.play_direction = 1  
        self.play_speed = 1.0    
        self.base_fps = 10       
        
        # Active file metadata cache trackers
        self.current_event_name = "N/A"
        self.current_file_duration = "N/A"
        
        # Build Menus & Layout Frames
        self.create_menu()
        self.create_widgets()
        
        self.root.protocol("WM_DELETE_WINDOW", self.safe_exit)
        self.root.bind_all("<Control-c>", lambda event: self.safe_exit())
        
        # Check for the specific H1 file and auto-load if present
        startup_file = 'H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5'
        if os.path.exists(startup_file):
            # 1. Explicitly set metadata to ensure centering logic works
            self.current_event_name = "GW150914"
            # Hardcoded GPS merger time for this specific file to match your metadata needs
            self.cached_target_gps = "1126259462" 
            
            # 2. Trigger the load
            # We call the worker directly. Note: process_pipeline_worker 
            # already handles its own threading for file I/O.
            self.root.after(1000, lambda: self.process_pipeline_worker(startup_file, 'H1'))            
            
            
    def safe_exit(self):
        self.is_playing = False
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
        file_menu.add_command(label="GWOSC Event Query Catalog...", command=self.open_gwosc_catalog_browser)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.safe_exit)
        menubar.add_cascade(label="File", menu=file_menu)
        
        param_menu = tk.Menu(menubar, tearoff=0)
        param_menu.add_command(label="Whiten Segmentation Config...", command=self.open_whiten_dialog)
        param_menu.add_command(label="FFT/Stiffness Configuration...", command=self.open_fft_dialog)
        param_menu.add_command(label="Display Limits & Window Width...", command=self.open_display_dialog)
        menubar.add_cascade(label="Parameters", menu=param_menu)
        
        self.root.config(menu=menubar)

    def copy_to_clipboard(self, text, name):
        if text and text not in ["N/A", "Unknown", "Unknown / Local"]:
            # Clean up the string slightly for the offset so it is purely numerical
            clean_text = str(text).replace("s in file", "").strip()
            self.root.clipboard_clear()
            self.root.clipboard_append(clean_text)
            self.global_status.config(text=f"Copied {name} to clipboard: {clean_text}")

    def create_widgets(self):
        # 1. Pack the global status bar FIRST so it claims the entire bottom edge
        self.global_status = ttk.Label(self.root, text="System: Idle", font=('Helvetica', 9, 'italic'), foreground="blue", padding=5)
        self.global_status.pack(side=tk.BOTTOM, fill=tk.X)

        self.whiten_interval = tk.DoubleVar(value=5.0)
        self.nperseg = tk.IntVar(value=256)
        self.nfft = tk.IntVar(value=1024)
        self.overlap_pct = tk.DoubleVar(value=85.0)
        
        self.t_center = tk.DoubleVar(value=0.0)
        self.t_width_seconds = 1.00
        self.t_width_str = tk.StringVar(value="0.65")
        self.f_min = tk.DoubleVar(value=0)
        self.f_max = tk.DoubleVar(value=500)
        self.pct_low = tk.DoubleVar(value=50.0)
        self.pct_high = tk.DoubleVar(value=99.99)
        
        sidebar = ttk.Frame(self.root, padding="10", width=340)
        sidebar.pack(side=tk.LEFT, fill=tk.Y)
        sidebar.pack_propagate(False)
        
        # Global Stream Information Monitor Desk
        meta_frame = ttk.LabelFrame(sidebar, text="Active Stream Hardware Specs", padding=8)
        meta_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.lbl_meta_len = ttk.Label(meta_frame, text="Sample Length: N/A", font=('Courier', 9, 'bold'))
        self.lbl_meta_len.pack(anchor=tk.W, pady=2)
        self.lbl_meta_rate = ttk.Label(meta_frame, text="Data Rate:     N/A", font=('Courier', 9, 'bold'))
        self.lbl_meta_rate.pack(anchor=tk.W, pady=2)
        
        # Clickable Event Offset
        self.lbl_meta_evt = ttk.Label(meta_frame, text="Event Offset:  N/A", font=('Courier', 9, 'bold'), foreground="darkorange", cursor="hand2")
        self.lbl_meta_evt.pack(anchor=tk.W, pady=2)
        self.lbl_meta_evt.bind("<Button-1>", lambda e: self.copy_to_clipboard(self.current_event_offset, "Event Offset"))

        ttk.Separator(meta_frame, orient='horizontal').pack(fill=tk.X, pady=5)
        
        # Clickable Event Identifier
        self.lbl_meta_name = ttk.Label(meta_frame, text="Event ID:      N/A", font=('Courier', 9, 'bold'), foreground="purple", cursor="hand2")
        self.lbl_meta_name.pack(anchor=tk.W, pady=2)
        self.lbl_meta_name.bind("<Button-1>", lambda e: self.copy_to_clipboard(self.current_event_name, "Event ID"))
        
        # Clickable File Names
        self.file_labels = {}
        for det in ['H1', 'L1', 'V1']:
            lbl = ttk.Label(meta_frame, text=f"{det} File: N/A", font=('Courier', 8), cursor="hand2", foreground="blue")
            lbl.pack(anchor=tk.W, pady=1)
            lbl.bind("<Button-1>", lambda e, d=det: self.copy_to_clipboard(self.current_filenames[d], f"{d} File Name"))
            self.file_labels[det] = lbl

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

        right_panel = ttk.Frame(self.root, padding="5")
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.notebook = ttk.Notebook(right_panel)
        self.notebook.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        self.tabs = {}
        for tab_name in ['H1', 'L1', 'V1', 'Joint Correlation']:
            frame = ttk.Frame(self.notebook)
            self.notebook.add(frame, text=tab_name)
            
            container = ttk.Frame(frame)
            container.pack(fill=tk.BOTH, expand=True)
            
            pbar = ttk.Progressbar(container, orient=tk.HORIZONTAL, mode='determinate')
            pbar.pack(fill=tk.X, padx=20, pady=5)
            pbar.pack_forget()
            
            fig, ax = plt.subplots(figsize=(8, 5))
            canvas = FigureCanvasTkAgg(fig, master=container)
            canvas_widget = canvas.get_tk_widget()
            canvas_widget.pack(fill=tk.BOTH, expand=True)
            
            canvas_widget.bind("<ButtonPress-1>", self.start_drag)
            canvas_widget.bind("<B1-Motion>", self.drag_motion)
            canvas_widget.bind("<ButtonRelease-1>", self.end_drag)
            
            self.tabs[tab_name] = {'fig': fig, 'ax': ax, 'canvas': canvas, 'pbar': pbar}
            
        self.notebook.bind("<<NotebookTabChanged>>", lambda e: self.update_all_tabs())

        # --- Cross-Correlation Control Panel ---
        self.corr_frame = ttk.LabelFrame(right_panel, text="Interferometer Time Delay Analysis", padding=10)
#        self.corr_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)

        self.btn_correlate = ttk.Button(self.corr_frame, text="Correlate Current Frame", command=self.trigger_frame_correlation)
        self.btn_correlate.pack(side=tk.LEFT, padx=5)

        self.lbl_offset_result = ttk.Label(self.corr_frame, text="Time Offset: Waiting for input...", font=('Courier', 10, 'bold'), foreground="blue")
        self.lbl_offset_result.pack(side=tk.LEFT, padx=15)
        # ---------------------------------------
        # ---------------------------------------

        # Integrated Control Console Layout
        controls_bar = ttk.LabelFrame(right_panel, text="Playback Control Desk & Global Positioning Timeline", padding="10")
        controls_bar.pack(side=tk.BOTTOM, fill=tk.X, pady=5)

        slider_frame = ttk.Frame(controls_bar)
        slider_frame.pack(fill=tk.X, pady=(0, 5))
        
        self.time_lbl = ttk.Label(slider_frame, text="Current Window Center: 0.00s", font=('Helvetica', 9, 'bold'))
        self.time_lbl.pack(side=tk.LEFT, padx=5)
        
        self.timeline_slider = ttk.Scale(slider_frame, from_=0, to=100, orient=tk.HORIZONTAL, command=self.on_slider_move)
        self.timeline_slider.pack(side=tk.RIGHT, fill=tk.X, expand=True, padx=5)

        btn_layout = ttk.Frame(controls_bar)
        btn_layout.pack(fill=tk.X, pady=(5, 0))

        ttk.Button(btn_layout, text="◀◀ Rev 2x", command=lambda: self.set_play(direction=-1, speed=2.0)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_layout, text="◀ Rev 1x", command=lambda: self.set_play(direction=-1, speed=1.0)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_layout, text="◂ Rev 0.5x", command=lambda: self.set_play(direction=-1, speed=0.5)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_layout, text=" Frame Rev (1/5)", command=lambda: self.step_frame(direction=-1)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_layout, text="█ Stop", command=self.stop_play).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_layout, text="Frame Fwd (1/5) ", command=lambda: self.step_frame(direction=1)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_layout, text="Play 0.5x ▸", command=lambda: self.set_play(direction=1, speed=0.5)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_layout, text="Play 1x ▶", command=lambda: self.set_play(direction=1, speed=1.0)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_layout, text="Play 2x ▶▶", command=lambda: self.set_play(direction=1, speed=2.0)).pack(side=tk.LEFT, padx=2)

        # Quick Jump Time Coordinate Entry Field
        jump_frame = ttk.Frame(btn_layout)
        jump_frame.pack(side=tk.RIGHT, padx=10)
        
        ttk.Label(jump_frame, text="Jump to Time (s):").pack(side=tk.LEFT, padx=2)
        self.jump_var = tk.StringVar()
        jump_ent = ttk.Entry(jump_frame, textvariable=self.jump_var, width=8)
        jump_ent.pack(side=tk.LEFT, padx=2)
        jump_ent.bind("<Return>", lambda e: self.execute_time_jump())
        ttk.Button(jump_frame, text="Go", command=self.execute_time_jump, width=4).pack(side=tk.LEFT, padx=2)

    def execute_time_jump(self):
        if self.total_duration == 0: return
        try:
            target = float(self.jump_var.get().strip())
            self.stop_play()
            self.clamp_and_set_center(target)
        except ValueError:
            messagebox.showerror("Format Error", "Please provide a valid numeric float coordinate value.")

    # --- Live GWOSC Catalog Window with Meta Reporting ---

    def pull_exact_detector_suite(self, event_name, duration_val, rate_val, gps_merger, files_dict):
    
        try:
            self.current_file_duration = f"{duration_val} Seconds"
            self.current_event_name = event_name
            
            # Standardize the target sampling matrix definition
            target_hz = 16384 if "16" in rate_val else 4096
            self.current_file_rate = f"{target_hz} Hz"
            self.cached_target_gps = gps_merger
            
            if not files_dict:
                self.root.after(0, lambda: messagebox.showwarning("Empty Matrix", f"No verified URLs found in the cache map for {event_name}."))
                return
                
            # Loop exclusively through the URLs confirmed natively by the API cache
            for det_key, data in files_dict.items():
                exact_url = data['url']
                threading.Thread(target=self.download_and_ingest, args=(exact_url, det_key), daemon=True).start()
                    
        except Exception as e:
            self.root.after(0, lambda err=e: messagebox.showerror("Pipeline Loader Error", str(err)))

    # --- Transportation/Playback Engine Mechanics ---
    def set_play(self, direction, speed):
        if self.total_duration == 0: return
        self.play_direction = direction
        self.play_speed = speed
        if not self.is_playing:
            self.is_playing = True
            self.play_step()

    def stop_play(self):
        self.is_playing = False

    def step_frame(self, direction):
        if self.total_duration == 0: return
        self.stop_play()
        shift = direction * (self.t_width_seconds / 5.0)
        new_center = self.t_center.get() + shift
        self.clamp_and_set_center(new_center)

    def play_step(self):
        if not self.is_playing: return
        step_increment = self.play_direction * (1.0 / self.base_fps) * self.play_speed
        new_center = self.t_center.get() + step_increment
        
        if new_center <= (self.t_width_seconds / 2.0) or new_center >= (self.total_duration - self.t_width_seconds / 2.0):
            self.is_playing = False
            return
            
        self.clamp_and_set_center(new_center)
        self.root.after(int(1000 / self.base_fps), self.play_step)

    # --- Accelerated & Cached GWOSC Catalog Window ---
    def open_gwosc_catalog_browser(self):
        import json
        import os
        import re
        import threading
        from concurrent.futures import ThreadPoolExecutor, as_completed
        import tkinter as tk
        from tkinter import ttk, messagebox

        CACHE_FILE = "gwosc_catalog_cache.json"

        win = tk.Toplevel(self.root)
        win.title("Verified GWOSC Stream Catalog (Exact Filenames)")
        win.geometry("1100x580")
        
        filter_frame = ttk.LabelFrame(win, text="Data Selection Criteria Filter", padding=10)
        filter_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(filter_frame, text="Substring Search Filter (e.g., 'GW150914', 'GW170817'):").pack(anchor=tk.W)
        search_var = tk.StringVar()
        search_entry = ttk.Entry(filter_frame, textvariable=search_var)
        search_entry.pack(fill=tk.X, pady=4)

        ttk.Label(filter_frame, text="4096 second samples will be very slow").pack(anchor=tk.W)

        results_frame = ttk.LabelFrame(win, text="Available Telemetry Configurations", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        columns = ('event_id', 'duration', 'rate', 'filenames')
        self.tree = ttk.Treeview(results_frame, columns=columns, show='headings', selectmode='browse')
        
        self.tree.heading('event_id', text='Event ID')
        self.tree.heading('duration', text='Duration (s)')
        self.tree.heading('rate', text='Rate')
        self.tree.heading('filenames', text='Exact Filenames')
        
        self.tree.column('event_id', width=120, anchor=tk.CENTER)
        self.tree.column('duration', width=90, anchor=tk.CENTER)
        self.tree.column('rate', width=80, anchor=tk.CENTER)
        self.tree.column('filenames', width=650, anchor=tk.W)
        
        self.tree.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
        
        scr = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.tree.yview)
        scr.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.config(yscrollcommand=scr.set)

        progress_frame = ttk.Frame(win)
        progress_frame.pack(fill=tk.X, padx=10, pady=(0, 5))
        
        cat_pbar = ttk.Progressbar(progress_frame, orient=tk.HORIZONTAL, mode='determinate')
        cat_pbar.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        cat_lbl = ttk.Label(progress_frame, text="Cache Status: Waiting...", font=('Helvetica', 9, 'italic'), foreground="blue")
        cat_lbl.pack(side=tk.RIGHT, padx=10)

        self.master_records = []

        def query_api(force_download=False):
            for i in self.tree.get_children(): self.tree.delete(i)
            self.master_records.clear()
            
            if not force_download and os.path.exists(CACHE_FILE):
                cat_lbl.config(text="Loading local disk cache...")
                try:
                    with open(CACHE_FILE, 'r') as f:
                        cached_data = json.load(f)
                    self.master_records = cached_data.get('records', [])
                    apply_filter()
                    cat_pbar.config(value=100)
                    cat_lbl.config(text=f"Loaded {len(self.master_records)} variants from cache.")
                    return
                except Exception:
                    cat_lbl.config(text="Cache broken, reverting to API download...", foreground="red")

            cat_lbl.config(text="Initializing API connection...", foreground="blue")
            cat_pbar.config(value=0)
            threading.Thread(target=download_and_compile_catalog, daemon=True).start()

        def process_single_event(name):
            """Worker target for parallel thread pooling execution."""
            from gwosc.locate import get_event_urls
            if not isinstance(name, str) or len(name) < 2: 
                return []
            
            all_urls = []
            
            # THE FIX: Explicitly request both 4kHz and 16kHz arrays from the server
            for target_hz in [4096, 16384]:
                try:
                    # Append the returned URLs to our master list for this event
                    urls = get_event_urls(name, sample_rate=target_hz)
                    all_urls.extend(urls)
                except Exception:
                    # If an event doesn't have a 16kHz release, the API throws an error. 
                    # We just silently catch it and move on.
                    pass
                    
            if not all_urls:
                return []
                
            variants = {}
            
            for url in all_urls:
                if not url.endswith('.hdf5'): continue
                filename = url.split("/")[-1]
                fname_upper = filename.upper()
                
                det = 'H1' if 'H-H1' in fname_upper else ('L1' if 'L-L1' in fname_upper else ('V1' if 'V-V1' in fname_upper else 'Unknown'))
                if det == 'Unknown': continue
                
                rate = "16kHz" if ("16K" in fname_upper or "_16_" in fname_upper) else "4kHz"
                
                try:
                    duration_str = re.search(r'-(\d+)\.HDF5$', fname_upper).group(1)
                except AttributeError:
                    duration_str = "Unknown"
                    
                key = (rate, duration_str)
                if key not in variants:
                    variants[key] = {}
                
                variants[key][det] = {'filename': filename, 'url': url}
            
            event_records = []
            for (rate, duration), files_dict in variants.items():
                filenames_str = " | ".join([f"{d}: {files_dict[d]['filename']}" for d in sorted(files_dict.keys())])
                event_records.append({
                    'event': name,
                    'duration': duration,
                    'rate': rate,
                    'files_dict': files_dict,
                    'filenames_str': filenames_str
                })
            return event_records

        def download_and_compile_catalog():
            try:
                from gwosc.datasets import find_datasets
                
                self.root.after(0, lambda: cat_lbl.config(text="Fetching master event signature list..."))
                all_events = find_datasets(type='event')
                if not all_events:
                    raise ValueError("No distinct event signatures returned from client query module.")
                
                compiled_records = []
                total_events = len(all_events)
                completed_count = 0
                
                # EXECUTION SPEED UP: Spin up 20 parallel network threads
                with ThreadPoolExecutor(max_workers=20) as executor:
                    future_to_event = {executor.submit(process_single_event, name): name for name in all_events}
                    
                    for future in as_completed(future_to_event):
                        name = future_to_event[future]
                        completed_count += 1
                        
                        # Dynamic thread-safe UI updates
                        pct = int((completed_count / total_events) * 100)
                        self.root.after(0, lambda p=pct, n=name, c=completed_count, t=total_events: [
                            cat_pbar.config(value=p),
                            cat_lbl.config(text=f"Processed {n} ({c}/{t})...")
                        ])
                        
                        result = future.result()
                        if result:
                            compiled_records.extend(result)
                
                self.master_records = compiled_records
                
                cache_structure = {'records': self.master_records}
                with open(CACHE_FILE, 'w') as f:
                    json.dump(cache_structure, f, indent=2)
                    
                self.root.after(0, apply_filter)
                self.root.after(0, lambda: [
                    cat_pbar.config(value=100),
                    cat_lbl.config(text=f"Complete! Indexed {len(self.master_records)} valid tracks.", foreground="green")
                ])
                
            except Exception as e:
                self.root.after(0, lambda err=e: [
                    messagebox.showerror("API Connection Error", f"Failed to parse catalog: {str(err)}"),
                    cat_lbl.config(text="Download failed.", foreground="red")
                ])
        def apply_filter(*args):
            # THE FIX: Ensure the widget still exists before trying to modify it
            if not self.tree.winfo_exists():
                return
                
            for i in self.tree.get_children(): self.tree.delete(i)
            query_string = search_var.get().strip().upper()
            
            for item in self.master_records:
                if not query_string or query_string in item['event'].upper():
                    self.tree.insert('', tk.END, values=(item['event'], item['duration'], item['rate'], item['filenames_str']))

        search_var.trace_add("write", apply_filter)
        
        ttk.Button(filter_frame, text="Force Clear Cache & Re-Download Server Registry", 
                   command=lambda: query_api(force_download=True)).pack(fill=tk.X, pady=5)

        def load_selected_record():
            sel = self.tree.selection()
            if not sel: return
            row_values = self.tree.item(sel[0])['values']
            chosen_event = str(row_values[0])
            chosen_duration = str(row_values[1])
            chosen_rate = str(row_values[2])
            
            record = next((r for r in self.master_records if str(r['event']) == chosen_event and str(r['duration']) == chosen_duration and str(r['rate']) == chosen_rate), None)
            
            if not record: return
            
            gps_merger = "N/A"
            try:
                from gwosc.datasets import event_gps
                gps_merger = event_gps(chosen_event)
            except Exception:
                pass
            
            win.destroy()
            self.global_status.config(text=f"Initiating direct mapped download for {chosen_event}...")
            
            threading.Thread(target=self.pull_exact_detector_suite, 
                             args=(chosen_event, record['duration'], record['rate'], gps_merger, record['files_dict']), 
                             daemon=True).start()

        btn = ttk.Button(win, text="Download Selected Stream Matrix", command=load_selected_record)
        btn.pack(fill=tk.X, padx=10, pady=10)
        
        self.tree.bind("<Double-1>", lambda event: load_selected_record())

        query_api(force_download=False)    
    

    def on_slider_move(self, val):
        if self.total_duration == 0 or self.is_dragging: return
        self.clamp_and_set_center(float(val), update_slider_ui=False)

    def clamp_and_set_center(self, target_center, update_slider_ui=True):
        min_c = self.t_width_seconds / 2.0
        max_c = self.total_duration - min_c
        clamped = max(min_c, min(target_center, max_c))
        
        self.t_center.set(clamped)
        self.time_lbl.config(text=f"Current Window Center: {clamped:.2f}s")
        
        if update_slider_ui:
            self.timeline_slider.set(clamped)
            
        self.update_all_tabs()

    # --- Config Panels ---
    def open_whiten_dialog(self):
        win = tk.Toplevel(self.root)
        win.title("Segmentation Window Properties")
        
        ttk.Label(win, text="Rolling Interval Chunk Length (s):").grid(row=0, column=0, padx=10, pady=10)
        ttk.Entry(win, textvariable=self.whiten_interval).grid(row=0, column=1, padx=10, pady=10)
        
        def save_and_rewhiten():
            win.destroy()
            active_detectors = [det for det in ['H1', 'L1', 'V1'] if self.detectors[det]['loaded']]
            
            if not active_detectors:
                self.global_status.config(text="Interval saved. (No active data loaded to re-whiten).")
                return
                
            self.global_status.config(text="Re-whitening active cache streams...")
            
            for det in active_detectors:
                threading.Thread(target=self.reprocess_active_cache_worker, args=(det,), daemon=True).start()

        ttk.Button(win, text="Save & Re-Whiten Data", command=save_and_rewhiten).grid(row=1, columnspan=2, pady=10)

    def reprocess_active_cache_worker(self, det):
        try:
            self.root.after(0, lambda: self.show_pbar(det))
            self.root.after(0, lambda: self.global_status.config(text=f"Re-whitening {det} data vectors..."))
            
            self.detectors[det]['whitened'] = self.whiten_by_intervals(
                self.detectors[det]['raw'], self.fs, self.whiten_interval.get(), det
            )
            
            self.root.after(0, lambda: self.update_pbar(det, 75))
            self.root.after(0, lambda: self.global_status.config(text=f"Re-generating 2D Spectrogram maps for {det}..."))
            
            self.compute_complete_spectrogram(det)
            
            self.root.after(0, lambda: self.update_pbar(det, 100))
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, self.update_all_tabs)
            self.root.after(0, lambda: self.global_status.config(text="System: Idle"))
        except Exception as e:
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, lambda: messagebox.showerror("Re-Whitening Error", str(e)))
            
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
        
        ttk.Label(win, text="Window Width (s or fraction e.g., 1/4):").grid(row=0, column=0, padx=10, pady=5, sticky=tk.W)
        ttk.Entry(win, textvariable=self.t_width_str).grid(row=0, column=1, padx=10, pady=5)
        
        ttk.Label(win, text="Lower Floor Percentile (%):").grid(row=1, column=0, padx=10, pady=5, sticky=tk.W)
        ttk.Entry(win, textvariable=self.pct_low).grid(row=1, column=1, padx=10, pady=5)
        
        ttk.Label(win, text="Upper Peak Percentile (%):").grid(row=2, column=0, padx=10, pady=5, sticky=tk.W)
        ttk.Entry(win, textvariable=self.pct_high).grid(row=2, column=1, padx=10, pady=5)
        
        ttk.Label(win, text="Lower Frequency Limit (f_min Hz):").grid(row=3, column=0, padx=10, pady=5, sticky=tk.W)
        ttk.Entry(win, textvariable=self.f_min).grid(row=3, column=1, padx=10, pady=5)
        
        ttk.Label(win, text="Upper Frequency Limit (f_max Hz):").grid(row=4, column=0, padx=10, pady=5, sticky=tk.W)
        ttk.Entry(win, textvariable=self.f_max).grid(row=4, column=1, padx=10, pady=5)
        
        def save_and_eval():
            raw_str = self.t_width_str.get().strip()
            try:
                if "/" in raw_str:
                    num, denom = map(float, raw_str.split("/"))
                    val = num / denom
                else:
                    val = float(raw_str)
                if val <= 0: raise ValueError()
                self.t_width_seconds = val
                
                float(self.f_min.get())
                float(self.f_max.get())
                
                win.destroy()
                self.update_all_tabs()
            except Exception:
                messagebox.showerror("Parsing Error", "Invalid input format detected. Please verify numerical configurations.")
                
        ttk.Button(win, text="Apply Changes", command=save_and_eval).grid(row=5, columnspan=2, pady=10)
        
    def render_canvas_frame(self, ax, canvas, Sxx, t, f, t_start, t_end):
        ax.clear()
        t_mask = (t >= t_start) & (t <= t_end)
        f_mask = (f >= self.f_min.get()) & (f <= self.f_max.get())
        
        if not np.any(t_mask) or not np.any(f_mask): 
            return
            
        Sxx_sub = Sxx[np.ix_(f_mask, t_mask)]
        vmin = np.percentile(Sxx_sub, self.pct_low.get())
        vmax = np.percentile(Sxx_sub, self.pct_high.get())

        # DEFENSIVE CHECK: Ensure range is significant
        if vmax <= vmin + 1e-12:
            vmax = vmin + 1e-9 # Force a small, detectable range
            
        levels = np.linspace(vmin, vmax, 25)
        
        # Use a try-except to catch potential contouring errors gracefully
        try:
            ax.pcolormesh(t[t_mask], f[f_mask], Sxx_sub, shading='gouraud', cmap='viridis')
            ax.contourf(t[t_mask], f[f_mask], Sxx_sub, levels=levels, cmap='inferno', extend='both')
        except ValueError:
            # Fallback if contouring still fails (e.g., constant Z data)
            ax.pcolormesh(t[t_mask], f[f_mask], Sxx_sub, shading='gouraud', cmap='viridis')
            
        ax.set_xlim(t_start, t_end)
        ax.set_ylim(self.f_min.get(), self.f_max.get())
        canvas.draw_idle()
        
    def start_drag(self, event):
        self.stop_play()
        self.is_dragging = True
        self.drag_start_x = event.x

    def drag_motion(self, event):
        if not self.is_dragging or self.total_duration == 0: 
            return
        
        # 1. Identify active tab axes
        active_tab_name = self.notebook.tab(self.notebook.select(), "text")
        if active_tab_name not in self.tabs: return
        ax = self.tabs[active_tab_name]['ax']
        
        # 2. Convert pixel coordinates (event.x) to data coordinates
        # inverted() translates display -> data
        inv = ax.transData.inverted()
        current_data_x, _ = inv.transform((event.x, event.y))
        start_data_x, _ = inv.transform((self.drag_start_x, 0))
        
        # 3. Calculate the shift
        data_shift = start_data_x - current_data_x
        new_center = self.t_center.get() + data_shift
        
        # 4. Apply
        self.clamp_and_set_center(new_center)
        
        # Update start position for the next drag event
        self.drag_start_x = event.x     
           
    def end_drag(self, event):
        self.is_dragging = False
        self.update_all_tabs()

    def load_local_detector(self, det):
        filename = filedialog.askopenfilename(filetypes=[("HDF5 keys", "*.hdf5 *.h5")])
        if filename:
            self.current_file_duration = "Local Readout"
            self.current_file_rate = "Extracting..."
            self.cached_target_gps = "N/A"
            self.global_status.config(text=f"Ingesting local array for {det}...")
            threading.Thread(target=self.process_pipeline_worker, args=(filename, det), daemon=True).start()
            
    def download_and_ingest(self, url, det):
        try:
            local_name = url.split("/")[-1]
            
            if not os.path.exists(local_name):
                # 1. Alert the UI that network contact has started
                self.root.after(0, lambda: self.global_status.config(text=f"Contacting GWOSC server for {det}..."))
                
                # 2. Use stream=True to prevent RAM explosion and add a timeout safety net
                with requests.get(url, stream=True, timeout=20) as r:
                    r.raise_for_status()
                    
                    # Safely calculate file size for the UI
                    total_size = int(r.headers.get('content-length', 0))
                    size_mb = total_size / (1024 * 1024)
                    
                    self.root.after(0, lambda: self.global_status.config(text=f"Downloading {det} stream ({size_mb:.1f} MB)..."))
                    
                    # 3. Stream in 8KB chunks directly to the hard drive
                    with open(local_name, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=8192): 
                            if chunk: 
                                f.write(chunk)
                                
            self.root.after(0, lambda: self.global_status.config(text=f"Download complete. Ingesting {det} array..."))
            
            # 4. Safely hand off to the processing worker
            self.process_pipeline_worker(local_name, det)
            
        except requests.exceptions.Timeout:
            self.root.after(0, lambda: self.global_status.config(text="System: Idle (Download Timeout)"))
            self.root.after(0, lambda: messagebox.showerror("Network Timeout", f"The server took too long to respond for {det}. It might be under heavy load."))
        except Exception as e:
            self.root.after(0, lambda: self.global_status.config(text="System: Idle (Download Failed)"))
            self.root.after(0, lambda: messagebox.showerror("Network Crash", f"Failed to fetch {det}:\n{str(e)}"))

    def show_pbar(self, det):
        self.tabs[det]['pbar'].pack(fill=tk.X, padx=20, pady=5, before=self.tabs[det]['canvas'].get_tk_widget())
        self.tabs[det]['pbar']['value'] = 0

    def update_pbar(self, det, val):
        self.tabs[det]['pbar']['value'] = val

    def hide_pbar(self, det):
        self.tabs[det]['pbar'].pack_forget()

    # --- Worker Core Processing Chain Layout ---
    def process_pipeline_worker(self, filepath, det):
        try:
            self.root.after(0, lambda: self.notebook.select(self.notebook.tabs()[['H1', 'L1', 'V1'].index(det)]))
            self.root.after(0, lambda: self.show_pbar(det))
            
            # Extract filename dynamically and update the corresponding clickable label
            filename_only = os.path.basename(filepath)
            self.current_filenames[det] = filename_only
            self.root.after(0, lambda d=det, fn=filename_only: self.file_labels[d].config(text=f"{d} File: {fn}"))

            with h5py.File(filepath, 'r') as f:
                key = 'strain/Strain' if 'strain/Strain' in f else 'strain/strain'
                self.detectors[det]['raw'] = f[key][:]
                self.dt = f[key].attrs['Xspacing']
                
                gps_start = f['meta/GPSstart'][()] if 'meta/GPSstart' in f else f['meta']['GPSstart'][()]
            
            self.fs = 1.0 / self.dt
            self.total_duration = len(self.detectors[det]['raw']) * self.dt
            
            if self.current_file_duration == "Local Readout":
                self.current_file_duration = f"{self.total_duration:.1f} Seconds"
                self.current_file_rate = f"{int(self.fs)} Hz"
                
            try:
                if self.cached_target_gps != "N/A" and float(self.cached_target_gps) > 0:
                    offset_calculated = float(self.cached_target_gps) - float(gps_start)
                    self.current_event_offset = f"{offset_calculated:.2f}s in file"
                else:
                    self.current_event_offset = "Unknown / Local"
            except Exception:
                self.current_event_offset = "Unknown"
                
            self.root.after(0, lambda: self.lbl_meta_len.config(text=f"Sample Length: {self.current_file_duration}"))
            self.root.after(0, lambda: self.lbl_meta_rate.config(text=f"Data Rate:     {self.current_file_rate}"))
            self.root.after(0, lambda: self.lbl_meta_evt.config(text=f"Event Offset:  {self.current_event_offset}"))
            self.root.after(0, lambda: self.lbl_meta_name.config(text=f"Event ID:      {self.current_event_name}"))
            
            self.detectors[det]['whitened'] = self.whiten_by_intervals(self.detectors[det]['raw'], self.fs, self.whiten_interval.get(), det)
            
            self.root.after(0, lambda: self.update_pbar(det, 75))
            self.root.after(0, lambda: self.global_status.config(text=f"Executing 2D Spectrogram & Laplacian Matrix loops for {det}..."))
            
            self.compute_complete_spectrogram(det)
            
            self.root.after(0, lambda: self.update_pbar(det, 100))
            self.detectors[det]['loaded'] = True

            # 1. Determine target time: Event offset if valid, else center of file
            target_time = self.total_duration / 2.0
            try:
                if self.cached_target_gps != "N/A" and float(self.cached_target_gps) > 0:
                    potential_offset = float(self.cached_target_gps) - float(gps_start)
                    # Validate: Ensure time is within valid bounds [0, total_duration]
                    if 0 <= potential_offset <= self.total_duration:
                        target_time = potential_offset
            except Exception:
                pass

            # 2. Update Slider and Jump exactly once
            if self.timeline_slider.cget('to') != self.total_duration:
                self.root.after(0, lambda: self.timeline_slider.config(to=self.total_duration))
            self.root.after(100, lambda: self.clamp_and_set_center(target_time))
            
            # 3. Final UI cleanup
            self.root.after(0, lambda: self.detectors[det]['status_lbl'].config(text="Status: Loaded & Cached", foreground="green"))
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, self.update_all_tabs)
            self.root.after(0, lambda: self.global_status.config(text="System: Idle"))

        except Exception as e:
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, lambda: messagebox.showerror("Pipeline Failure", str(e)))
    def whiten_by_intervals(self, strain, fs, interval_sec, det):
        N = len(strain)
        chunk_size = int(interval_sec * fs)
        whitened = np.zeros_like(strain)
        window_sum = np.zeros_like(strain)  
        
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
                
            clean_chunk = np.fft.irfft(whitened_spec, n=chunk_size) * win
            
            whitened[start:end] += clean_chunk
            window_sum[start:end] += win * win  
            
            pct = int((idx / total_steps) * 70)
            if idx % 5 == 0:
                self.root.after(0, lambda p=pct: self.update_pbar(det, p))
                
        with np.errstate(divide='ignore', invalid='ignore'):
            whitened = np.where(window_sum > 1e-10, whitened / window_sum, 0.0)
            
        return whitened
        
    def trigger_frame_correlation(self):
        self.lbl_offset_result.config(text="Calculating...", foreground="orange")
        self.root.update_idletasks()
        
        try:
            # 1. Access whitened strain data
            h1_data = self.detectors['H1']['whitened']
            l1_data = self.detectors['L1']['whitened']
            
            # 2. Slice to the UI window (Crucial to match the spectrogram view)
            t_center = self.t_center.get()
            t_width = self.t_width_seconds
            t_start, t_end = max(0.0, t_center - t_width/2), min(self.total_duration, t_center + t_width/2)
            
            idx_min = int(t_start * self.fs)
            idx_max = int(t_end * self.fs)
            
            h1_slice = h1_data[idx_min:idx_max]
            l1_slice = l1_data[idx_min:idx_max]
            
            # 3. Compute FFTs on the fly
            h1_fft = np.fft.rfft(h1_slice)
            l1_fft = np.fft.rfft(l1_slice)
            freqs = np.fft.rfftfreq(len(h1_slice), 1/self.fs)
            
            # 4. Call the alignment method
            offset_ms, lead_text = self.calculate_delay_from_existing_ffts(
                h1_fft, 
                l1_fft, 
                freqs, 
                max(0.1, self.f_min.get()), 
                self.f_max.get()
            )
            
            self.lbl_offset_result.config(
                text=f"Time Offset: {offset_ms:.2f} ms ({lead_text})", 
                foreground="green"
            )
            
        except Exception as e:
            self.lbl_offset_result.config(text=f"Error: {str(e)}", foreground="red")
            print(f"Correlation Error: {e}")

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
        if self.total_duration == 0: return
        t_center = self.t_center.get()
        t_width = self.t_width_seconds
        t_start, t_end = max(0.0, t_center - t_width/2), min(self.total_duration, t_center + t_width/2)
        for det in ['H1', 'L1', 'V1']:
            if self.detectors[det]['loaded']:
                self.render_canvas_frame(self.tabs[det]['ax'], self.tabs[det]['canvas'], self.detectors[det]['Sxx'], self.detectors[det]['t'], self.detectors[det]['f'], t_start, t_end)
        self.render_joint_correlation(t_start, t_end)

    def render_joint_correlation(self, t_start, t_end):
        ax = self.tabs['Joint Correlation']['ax']
        canvas = self.tabs['Joint Correlation']['canvas']
        ax.clear()
        
        # 1. Identify valid matrices
        active_matrices = [det for det in ['H1', 'L1', 'V1'] 
                           if self.detectors[det]['loaded'] and self.detectors[det]['active_corr'].get()]
        if not active_matrices: return
        
        # 2. Use the FIRST detector as the reference clock
        ref = active_matrices[0]
        t, f = self.detectors[ref]['t'], self.detectors[ref]['f']
        
        # 3. Create initial mask based on time window
        t_mask = (t >= t_start) & (t <= t_end)
        f_mask = (f >= self.f_min.get()) & (f <= self.f_max.get())
        if not np.any(t_mask) or not np.any(f_mask): return
        
        # 4. Initialize the product using the reference detector's slice
        joint_product = None
        
        for det in active_matrices:
            # SAFETY: Ensure we don't exceed the bounds of THIS specific detector's matrix
            Sxx = self.detectors[det]['Sxx']
            
            # Create a localized mask that respects this matrix's shape
            # Sxx.shape[1] is the time dimension
            valid_t_mask = t_mask.copy()
            if np.sum(valid_t_mask) > Sxx.shape[1]:
                # If our mask is too wide, truncate it to the matrix size
                valid_t_mask = np.zeros(len(t), dtype=bool)
                valid_t_mask[:Sxx.shape[1]] = True
                valid_t_mask &= (t >= t_start) & (t <= t_end)
            
            # Slice with safe indices
            try:
                S_sub = Sxx[np.ix_(f_mask, valid_t_mask)]
                
                # Normalize and multiply
                norm_S = S_sub / (np.max(S_sub) + 1e-20)
                if joint_product is None:
                    joint_product = norm_S
                else:
                    # Match shapes if they differ slightly due to sample rate variations
                    min_h = min(joint_product.shape[0], norm_S.shape[0])
                    min_w = min(joint_product.shape[1], norm_S.shape[1])
                    joint_product = joint_product[:min_h, :min_w] * norm_S[:min_h, :min_w]
            except IndexError:
                continue # Skip this detector if it's incompatible
        
        if joint_product is None: return

        vmin, vmax = np.percentile(joint_product, self.pct_low.get()), np.percentile(joint_product, self.pct_high.get())
        levels = np.linspace(vmin, max(vmax, vmin + 1e-10), 25)
        ax.pcolormesh(t[t_mask], f[f_mask], joint_product, shading='gouraud', cmap='magma')
        ax.contourf(t[t_mask], f[f_mask], joint_product, levels=levels, cmap='hot', extend='both')
        ax.set_xlim(t_start, t_end)
        ax.set_ylim(self.f_min.get(), self.f_max.get())
        canvas.draw_idle()

    def calculate_frame_correlation(self, h1_strain, l1_strain, time_array, sample_rate, current_t_min, current_t_max, f_min, f_max):
        # 1. Slice to the UI window
        idx_min = np.searchsorted(time_array, current_t_min)
        idx_max = np.searchsorted(time_array, current_t_max)
        h1 = h1_strain[idx_min:idx_max]
        l1 = l1_strain[idx_min:idx_max]
    
        # 2. Apply bandpass filter to the merger band
        # This isolates the signal from seismic/instrumental noise
        sos = signal.butter(4, [max(0.1, f_min), f_max], btype='bandpass', fs=sample_rate, output='sos')
        h1_filt = signal.sosfiltfilt(sos, h1)
        l1_filt = signal.sosfiltfilt(sos, l1)
    
        # 3. Compute Instantaneous Power: P(t) = s(t)^2
        h1_power = h1_filt**2
        l1_power = l1_filt**2
    
        # 4. Smooth the power (Moving Average)
        window_size = int(0.01 * sample_rate) 
        h1_smooth = np.convolve(h1_power, np.ones(window_size)/window_size, mode='same')
        l1_smooth = np.convolve(l1_power, np.ones(window_size)/window_size, mode='same')
    
        # 5. Cross-correlate the filtered, smoothed power time series
        correlation = signal.correlate(h1_smooth, l1_smooth, mode='full')
        lags = signal.correlation_lags(len(h1_smooth), len(l1_smooth), mode='full')
    
        # 6. Find the peak
        peak_idx = np.argmax(correlation)
        tau = lags[peak_idx] / sample_rate
    
        return abs(tau * 1000), ("L1 leads H1" if tau > 0 else "H1 leads L1")        


    def calculate_delay_from_existing_ffts(self, H1_fft, L1_fft, freqs, f_min, f_max):
        mask = (freqs >= max(0.1, f_min)) & (freqs <= f_max)
        cross_spectrum = H1_fft[mask] * np.conj(L1_fft[mask])
        phase = np.unwrap(np.angle(cross_spectrum))
        slope, _ = np.polyfit(freqs[mask], phase, 1)
        tau = slope / (2 * np.pi)
        return abs(tau * 1000), ("L1 leads H1" if tau > 0 else "H1 leads L1")
            
if __name__ == "__main__":
    root = tk.Tk(className="GWExplorer")
    app = GWExplorerApp(root)
    root.mainloop()    
    
