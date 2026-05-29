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
import matplotlib.patches
import json
import webbrowser
from datetime import datetime, timezone, timedelta
from urllib.parse import urlencode

C_LIGHT = 299792458.0

DETECTOR_ECEF = {
    "H1": np.array([-2161414.926, -3834695.178, 4600350.226]),
    "L1": np.array([-74276.044, -5496283.719, 3224257.017]),
    "V1": np.array([4546374.099, 842989.697, 4378576.963]),
}


class RangeSlider(tk.Canvas):
    """Simple two-handle horizontal range slider for Tkinter."""

    def __init__(
            self,
            master,
            from_,
            to,
            low_var,
            high_var,
            command=None,
            width=190,
            height=42,
            min_gap=0.0,
            decimals=1,
            **kwargs
    ):
        super().__init__(
            master,
            width=width,
            height=height,
            highlightthickness=0,
            **kwargs
        )
        self.from_ = float(from_)
        self.to = float(to)
        self.low_var = low_var
        self.high_var = high_var
        self.command = command
        self.slider_width = int(width)
        self.slider_height = int(height)
        self.pad = 14
        self.handle_radius = 6
        self.min_gap = float(min_gap)
        self.decimals = int(decimals)
        self.active_handle = None

        self.bind("<Configure>", lambda event: self.redraw())
        self.bind("<Button-1>", self._on_press)
        self.bind("<B1-Motion>", self._on_drag)
        self.bind("<ButtonRelease-1>", self._on_release)

        try:
            self.low_var.trace_add("write", lambda *_: self.redraw())
            self.high_var.trace_add("write", lambda *_: self.redraw())
        except Exception:
            pass

        self.redraw()

    def _track_bounds(self):
        width = max(self.slider_width, self.winfo_width())
        return self.pad, max(self.pad + 1, width - self.pad)

    def _clip(self, value):
        return max(self.from_, min(self.to, float(value)))

    def _round(self, value):
        if self.decimals <= 0:
            return round(value)
        return round(value, self.decimals)

    def _value_to_x(self, value):
        left, right = self._track_bounds()
        frac = (self._clip(value) - self.from_) / (self.to - self.from_)
        return left + frac * (right - left)

    def _x_to_value(self, x):
        left, right = self._track_bounds()
        frac = (float(x) - left) / (right - left)
        value = self.from_ + max(0.0, min(1.0, frac)) * (self.to - self.from_)
        return self._round(value)

    def _current_values(self):
        try:
            low = self._clip(float(self.low_var.get()))
        except Exception:
            low = self.from_
        try:
            high = self._clip(float(self.high_var.get()))
        except Exception:
            high = self.to

        if high < low:
            low, high = high, low
        if high - low < self.min_gap:
            high = min(self.to, low + self.min_gap)
            low = max(self.from_, high - self.min_gap)
        return low, high

    def redraw(self):
        self.delete("all")
        low, high = self._current_values()
        left, right = self._track_bounds()
        y = max(10, self.winfo_height() // 2)

        low_x = self._value_to_x(low)
        high_x = self._value_to_x(high)

        self.create_line(left, y, right, y, width=4)
        self.create_line(low_x, y, high_x, y, width=6)

        for x in (low_x, high_x):
            self.create_oval(
                x - self.handle_radius,
                y - self.handle_radius,
                x + self.handle_radius,
                y + self.handle_radius,
                width=2,
                fill=self.cget("background") or "SystemButtonFace"
            )

    def _on_press(self, event):
        low, high = self._current_values()
        low_x = self._value_to_x(low)
        high_x = self._value_to_x(high)
        self.active_handle = "low" if abs(event.x - low_x) <= abs(event.x - high_x) else "high"
        self._set_from_x(event.x)

    def _on_drag(self, event):
        self._set_from_x(event.x)

    def _on_release(self, event):
        self._set_from_x(event.x)
        self.active_handle = None

    def _set_from_x(self, x):
        if self.active_handle is None:
            return

        value = self._clip(self._x_to_value(x))
        low, high = self._current_values()

        if self.active_handle == "low":
            value = min(value, high - self.min_gap)
            value = self._clip(value)
            self.low_var.set(value)
        else:
            value = max(value, low + self.min_gap)
            value = self._clip(value)
            self.high_var.set(value)

        self.redraw()
        if self.command is not None:
            self.command()


class GWExplorerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Multi-Detector CMST Transient Explorer")
        self.root.geometry("1450x950")

        self.display_defaults_path = "display_defaults.json"
        self.recent_files_path = "recent_files.json"
        self.recent_files = self.load_recent_files()

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

        # Core Global State Engine
        self.detectors = {
            'H1': {'raw': None, 'whitened': None, 'Sxx': None, 't': None, 'f': None, 'loaded': False,
                   'active_corr': tk.BooleanVar(value=True)},
            'L1': {'raw': None, 'whitened': None, 'Sxx': None, 't': None, 'f': None, 'loaded': False,
                   'active_corr': tk.BooleanVar(value=True)},
            'V1': {'raw': None, 'whitened': None, 'Sxx': None, 't': None, 'f': None, 'loaded': False,
                   'active_corr': tk.BooleanVar(value=False)}
        }

        self.dt = None
        self.fs = None
        self.total_duration = 0.0
        self.is_dragging = False
        self.drag_start_x = None

        self.current_file_duration = "N/A"
        self.current_file_rate = "N/A"
        self.current_event_offset = "N/A"
        self.current_filenames = {'H1': 'N/A', 'L1': 'N/A', 'V1': 'N/A'}

        # Playback Control Variables
        self.is_playing = False
        self.play_direction = 1
        self.play_speed = 1.0
        self.base_fps = 10
        self.current_event_name = "N/A"

        self.create_menu()
        self.create_widgets()

        self.root.protocol("WM_DELETE_WINDOW", self.safe_exit)
        self.root.bind_all("<Control-c>", lambda event: self.safe_exit())

        # Check for the specific H1 file and auto-load if present
        startup_files = {
            "H1": "H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5",
            "L1": "L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5",
        }

        available_startup_files = {
            det: filename
            for det, filename in startup_files.items()
            if os.path.exists(filename)
        }

        if available_startup_files:
            self.current_event_name = "GW150914"
            self.current_file_duration = "32 Seconds"
            self.current_file_rate = "16384 Hz"
            self.cached_target_gps = "1126259462"

            def load_startup_suite():
                for det, filename in available_startup_files.items():
                    threading.Thread(
                        target=self.process_pipeline_worker,
                        args=(filename, det),
                        kwargs={"add_to_recent": False},
                        daemon=True
                    ).start()

            self.root.after(1000, load_startup_suite)

    def load_display_defaults(self):
        defaults = {
            "pct_low": 50.0,
            "pct_high": 99.99,
            "f_min": 0.0,
            "f_max": 500.0,
            "t_width_seconds": 0.5,
            "show_grid": False,
            "show_contours": False,
            "contour_count": 25,
            "colormap_name": "inferno",
        }

        try:
            if os.path.exists(self.display_defaults_path):
                with open(self.display_defaults_path, "r") as f:
                    saved = json.load(f)
                defaults.update(saved)
        except Exception as e:
            print("Display defaults load error:", e)

        return defaults

    def save_display_defaults(self):
        self.display_defaults = self.get_current_display_settings()
        try:
            with open(self.display_defaults_path, "w") as f:
                json.dump(self.display_defaults, f, indent=2)
            self.global_status.config(text="Saved display settings as default.")
        except Exception as e:
            messagebox.showerror("Save Error", f"Could not save display defaults:\n{e}")

    def load_recent_files(self):
        try:
            if os.path.exists(self.recent_files_path):
                with open(self.recent_files_path, "r") as f:
                    return json.load(f)
        except Exception:
            pass
        return []

    def save_recent_files(self):
        try:
            with open(self.recent_files_path, "w") as f:
                json.dump(self.recent_files[:20], f, indent=2)
        except Exception as e:
            print("Recent save error:", e)

    def add_recent_catalog_entry(self, event_name, duration_val, rate_val, gps_merger, files_dict):
        entry = {
            "source": "catalog",
            "event": event_name,
            "duration": str(duration_val),
            "rate": str(rate_val),
            "gps": str(gps_merger),
            "files_dict": files_dict
        }

        self.recent_files = [
            x for x in self.recent_files
            if not (
                    x.get("source") == "catalog"
                    and x.get("event") == entry["event"]
                    and x.get("duration") == entry["duration"]
                    and x.get("rate") == entry["rate"]
            )
        ]

        self.recent_files.insert(0, entry)
        self.save_recent_files()
        self.refresh_recent_menu()

    def add_recent_file(self, filepath, det):
        entry = {
            "source": "local",
            "path": filepath,
            "det": det
        }

        self.recent_files = [
            x for x in self.recent_files
            if not (
                    x.get("source") == "local"
                    and x.get("path") == filepath
                    and x.get("det") == det
            )
        ]

        self.recent_files.insert(0, entry)
        self.save_recent_files()
        self.refresh_recent_menu()

    def refresh_recent_menu(self):
        self.recent_menu.delete(0, tk.END)

        if not self.recent_files:
            self.recent_menu.add_command(label="(Empty)", state=tk.DISABLED)
            return

        for item in self.recent_files[:20]:
            if item.get("source") == "catalog":
                label = f"{item['event']} | {item['rate']} | {item['duration']}s"
                self.recent_menu.add_command(
                    label=label,
                    command=lambda x=item: self.load_recent_catalog_entry(x)
                )
            else:
                path = item.get("path", "")
                det = item.get("det", "?")
                label = f"{det}: {os.path.basename(path)}"
                self.recent_menu.add_command(
                    label=label,
                    command=lambda p=path, d=det: self.load_recent_direct(p, d)
                )

    def load_recent_catalog_entry(self, item):
        self.global_status.config(
            text=f"Loading recent catalog entry {item['event']}..."
        )
        threading.Thread(
            target=self.pull_exact_detector_suite,
            args=(
                item["event"],
                item["duration"],
                item["rate"],
                item["gps"],
                item["files_dict"]
            ),
            daemon=True
        ).start()

    def load_recent_direct(self, filepath, det):
        if not os.path.exists(filepath):
            messagebox.showerror("Missing File", f"Could not locate:\n{filepath}")
            return
        self.global_status.config(text=f"Loading recent {det} file...")
        threading.Thread(target=self.process_pipeline_worker, args=(filepath, det), daemon=True).start()

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
        self.recent_menu = tk.Menu(file_menu, tearoff=0)

        file_menu.add_command(label="GWOSC Event Query Catalog...", command=self.open_gwosc_catalog_browser)
        file_menu.add_cascade(label="Recent Files", menu=self.recent_menu)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.safe_exit)
        menubar.add_cascade(label="File", menu=file_menu)

        param_menu = tk.Menu(menubar, tearoff=0)
        param_menu.add_command(label="Whiten Segmentation Config...", command=self.open_whiten_dialog)
        param_menu.add_command(label="FFT/Stiffness Configuration...", command=self.open_fft_dialog)
        param_menu.add_command(label="Display Limits & Window Width...", command=self.open_display_dialog)
        menubar.add_cascade(label="Parameters", menu=param_menu)

        self.root.config(menu=menubar)
        self.refresh_recent_menu()

    def copy_to_clipboard(self, text, name):
        if text and text not in ["N/A", "Unknown", "Unknown / Local"]:
            clean_text = str(text).replace("s in file", "").strip()
            self.root.clipboard_clear()
            self.root.clipboard_append(clean_text)
            self.global_status.config(text=f"Copied {name} to clipboard: {clean_text}")

    def create_widgets(self):
        self.global_status = ttk.Label(self.root, text="System: Idle", font=('Helvetica', 9, 'italic'),
                                       foreground="blue", padding=5)
        self.global_status.pack(side=tk.BOTTOM, fill=tk.X)

        self._init_ui_variables()
        self.display_defaults = self.load_display_defaults()
        self._apply_display_defaults()

        self.sidebar = ttk.Frame(self.root, padding="10", width=340)
        self.sidebar.pack(side=tk.LEFT, fill=tk.Y)
        self.sidebar.pack_propagate(False)

        self.right_panel = ttk.Frame(self.root, padding="5")
        self.right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.chart_area = ttk.Frame(self.right_panel)
        self.chart_area.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self._build_sidebar()
        self._build_notebook_tabs(self.chart_area)
        self.build_display_side_controls(self.chart_area)
        self._build_correlation_panel()
        self._build_playback_controls()

    def _init_ui_variables(self):
        self.whiten_interval = tk.DoubleVar(value=5.0)
        self.nperseg = tk.IntVar(value=256)
        self.nfft = tk.IntVar(value=1024)
        self.overlap_pct = tk.DoubleVar(value=85.0)
        self.t_center = tk.DoubleVar(value=0.0)
        self.t_width_seconds = 0.2
        self.t_width_str = tk.StringVar(value="0.65")
        self.f_min = tk.DoubleVar(value=0)
        self.f_max = tk.DoubleVar(value=500)
        self.pct_low = tk.DoubleVar(value=50.0)
        self.pct_high = tk.DoubleVar(value=99.99)
        self.time_window_log = tk.DoubleVar(value=np.log10(self.t_width_seconds))
        self.show_grid = tk.BooleanVar(value=False)
        self.show_contours = tk.BooleanVar(value=False)
        self.contour_count = tk.IntVar(value=25)
        self.colormap_name = tk.StringVar(value="inferno")

        self.detector_offsets_ms = {"H1": tk.DoubleVar(value=0.0), "L1": tk.DoubleVar(value=0.0),
                                    "V1": tk.DoubleVar(value=0.0)}
        self.signal_flips = {"H1": tk.BooleanVar(value=False), "L1": tk.BooleanVar(value=False),
                             "V1": tk.BooleanVar(value=False)}

        self.display_control_update_job = None
        self.display_controls_recenter = False
        self.file_labels = {}
        self.tabs = {}
        self.offset_spinboxes = {}

        self.sky_result_text = tk.StringVar(value="Sky: Waiting...")
        self.sky_link_text = tk.StringVar(value="N/A")
        self.sky_link_url = None
        self.last_sky_ra_deg = None
        self.last_sky_dec_deg = None
        self.sky_recompute_job = None
        self.jump_var = tk.StringVar()

    def _apply_display_defaults(self):
        self.pct_low.set(self.display_defaults["pct_low"])
        self.pct_high.set(self.display_defaults["pct_high"])
        self.f_min.set(self.display_defaults["f_min"])
        self.f_max.set(self.display_defaults["f_max"])
        self.t_width_seconds = self.display_defaults["t_width_seconds"]
        self.t_width_str.set(f"{self.t_width_seconds:.3f}")
        self.time_window_log.set(np.log10(self.t_width_seconds))
        self.show_grid.set(self.display_defaults["show_grid"])
        self.show_contours.set(self.display_defaults["show_contours"])
        self.contour_count.set(self.display_defaults["contour_count"])
        self.colormap_name.set(self.display_defaults["colormap_name"])

    def _build_sidebar(self):
        meta_frame = ttk.LabelFrame(self.sidebar, text="Active Stream Hardware Specs", padding=8)
        meta_frame.pack(fill=tk.X, pady=(0, 10))

        self.lbl_meta_len = ttk.Label(meta_frame, text="Sample Length: N/A", font=('Courier', 9, 'bold'))
        self.lbl_meta_len.pack(anchor=tk.W, pady=2)
        self.lbl_meta_rate = ttk.Label(meta_frame, text="Data Rate:     N/A", font=('Courier', 9, 'bold'))
        self.lbl_meta_rate.pack(anchor=tk.W, pady=2)

        self.lbl_meta_evt = ttk.Label(meta_frame, text="Event Offset:  N/A", font=('Courier', 9, 'bold'),
                                      foreground="darkorange", cursor="hand2")
        self.lbl_meta_evt.pack(anchor=tk.W, pady=2)
        self.lbl_meta_evt.bind("<Button-1>",
                               lambda e: self.copy_to_clipboard(self.current_event_offset, "Event Offset"))

        ttk.Separator(meta_frame, orient='horizontal').pack(fill=tk.X, pady=5)

        self.lbl_meta_name = ttk.Label(meta_frame, text="Event ID:      N/A", font=('Courier', 9, 'bold'),
                                       foreground="purple", cursor="hand2")
        self.lbl_meta_name.pack(anchor=tk.W, pady=2)
        self.lbl_meta_name.bind("<Button-1>", lambda e: self.copy_to_clipboard(self.current_event_name, "Event ID"))

        for det in ['H1', 'L1', 'V1']:
            lbl = ttk.Label(meta_frame, text=f"{det} File: N/A", font=('Courier', 8), cursor="hand2", foreground="blue")
            lbl.pack(anchor=tk.W, pady=1)
            lbl.bind("<Button-1>", lambda e, d=det: self.copy_to_clipboard(self.current_filenames[d], f"{d} File Name"))
            self.file_labels[det] = lbl

        ttk.Label(self.sidebar, text="Active Detector Loading:", font=('Helvetica', 10, 'bold')).pack(anchor=tk.W,
                                                                                                      pady=5)
        for det in ['H1', 'L1', 'V1']:
            f = ttk.LabelFrame(self.sidebar, text=f"Interferometer {det}", padding=5)
            f.pack(fill=tk.X, pady=4)
            ttk.Button(f, text="Load Local HDF5...", command=lambda d=det: self.load_local_detector(d)).pack(fill=tk.X)
            self.detectors[det]['status_lbl'] = ttk.Label(f, text="Status: Empty", foreground="gray")
            self.detectors[det]['status_lbl'].pack(anchor=tk.W)

        ttk.Separator(self.sidebar, orient='horizontal').pack(fill=tk.X, pady=10)

        ttk.Label(self.sidebar, text="Correlation Mask Switches:", font=('Helvetica', 10, 'bold')).pack(anchor=tk.W,
                                                                                                        pady=5)
        for det in ['H1', 'L1', 'V1']:
            ttk.Checkbutton(self.sidebar, text=f"Include {det} in Joint Product",
                            variable=self.detectors[det]['active_corr'],
                            command=self.update_all_tabs).pack(anchor=tk.W, pady=2)

    def _build_notebook_tabs(self, parent):
        self.notebook = ttk.Notebook(parent)
        self.notebook.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

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

            canvas_widget.bind("<Shift-ButtonPress-1>", self.start_zoom_drag)
            canvas_widget.bind("<Shift-B1-Motion>", self.zoom_drag_motion)
            canvas_widget.bind("<Shift-ButtonRelease-1>", self.end_zoom_drag)
            canvas_widget.bind("<ButtonPress-1>", self.start_drag)
            canvas_widget.bind("<B1-Motion>", self.drag_motion)
            canvas_widget.bind("<ButtonRelease-1>", self.end_drag)

            self.tabs[tab_name] = {'fig': fig, 'ax': ax, 'canvas': canvas, 'pbar': pbar}

        wave_frame = ttk.Frame(self.notebook)
        self.notebook.add(wave_frame, text="Whitened Waveforms")

        wave_fig, wave_axs = plt.subplots(4, 1, figsize=(8, 8), sharex=True)
        wave_fig.subplots_adjust(left=0.08, right=0.95, bottom=0.10, top=0.95, hspace=0.3)

        wave_canvas = FigureCanvasTkAgg(wave_fig, master=wave_frame)
        wave_canvas_widget = wave_canvas.get_tk_widget()
        wave_canvas_widget.pack(fill=tk.BOTH, expand=True)

        self.tabs["Whitened Waveforms"] = {'fig': wave_fig, 'ax': wave_axs, 'canvas': wave_canvas, 'pbar': None}
        self.notebook.bind("<<NotebookTabChanged>>", lambda e: self.update_all_tabs())

    def _build_correlation_panel(self):
        self.corr_frame = ttk.LabelFrame(self.right_panel, text="Interferometer Time Delay Analysis", padding=10)
        self.corr_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)

        top_row = ttk.Frame(self.corr_frame)
        top_row.pack(fill=tk.X, pady=(0, 6))

        self.btn_correlate = ttk.Button(top_row, text="Correlate Current Frame", command=self.trigger_frame_correlation)
        self.btn_correlate.pack(side=tk.LEFT, padx=5)

        self.lbl_offset_result = ttk.Label(top_row, textvariable=self.sky_result_text, font=('Courier', 10, 'bold'),
                                           foreground="blue", cursor="hand2")
        self.lbl_offset_result.pack(side=tk.LEFT, padx=15)
        self.lbl_offset_result.bind("<Button-1>",
                                    lambda e: self.copy_to_clipboard(self.sky_result_text.get(), "Sky Location"))

        bottom_row = ttk.Frame(self.corr_frame)
        bottom_row.pack(fill=tk.X)

        ttk.Label(bottom_row, text="Offsets ms:").pack(side=tk.LEFT, padx=(5, 4))

        for det in ["L1", "V1"]:
            ttk.Label(bottom_row, text=f"{det}:").pack(side=tk.LEFT, padx=(8, 2))
            spin = ttk.Spinbox(bottom_row, from_=-30.0, to=30.0, increment=0.1, format="%.1f", width=6,
                               textvariable=self.detector_offsets_ms[det], command=self.on_offset_spinbox_changed)
            spin.pack(side=tk.LEFT)
            spin.bind("<KeyRelease>", self.on_offset_spinbox_changed)
            spin.bind("<Return>", self.on_offset_spinbox_changed)
            spin.bind("<FocusOut>", self.on_offset_spinbox_changed)
            self.offset_spinboxes[det] = spin

        ttk.Button(bottom_row, text="Sky Location", command=self.estimate_sky_location).pack(side=tk.LEFT, padx=10)

        self.lbl_sky_link = ttk.Label(bottom_row, textvariable=self.sky_link_text, font=("Courier", 10, "bold"),
                                      foreground="purple", cursor="hand2")
        self.lbl_sky_link.pack(side=tk.LEFT, padx=10)
        self.lbl_sky_link.bind("<Button-1>", lambda e: self.open_sky_link())

        flip_row = ttk.Frame(self.corr_frame)
        flip_row.pack(fill=tk.X, pady=(5, 0))
        ttk.Label(flip_row, text="Invert Polarity (Bottom Chart):").pack(side=tk.LEFT, padx=(5, 4))
        for det in ["H1", "L1", "V1"]:
            ttk.Checkbutton(flip_row, text=det, variable=self.signal_flips[det], command=self.update_all_tabs).pack(
                side=tk.LEFT, padx=5)

    def _build_playback_controls(self):
        controls_bar = ttk.LabelFrame(self.right_panel, text="Playback Control Desk & Global Positioning Timeline",
                                      padding="10")
        controls_bar.pack(side=tk.BOTTOM, fill=tk.X, pady=5)

        slider_frame = ttk.Frame(controls_bar)
        slider_frame.pack(fill=tk.X, pady=(0, 5))

        self.time_lbl = ttk.Label(slider_frame, text="Current Window Center: 0.00s", font=('Helvetica', 9, 'bold'))
        self.time_lbl.pack(side=tk.LEFT, padx=5)

        self.timeline_slider = ttk.Scale(slider_frame, from_=0, to=100, orient=tk.HORIZONTAL,
                                         command=self.on_slider_move)
        self.timeline_slider.pack(side=tk.RIGHT, fill=tk.X, expand=True, padx=5)

        btn_layout = ttk.Frame(controls_bar)
        btn_layout.pack(fill=tk.X, pady=(5, 0))

        ttk.Button(btn_layout, text="◀◀ 2x", command=lambda: self.set_play(direction=-1, speed=2.0)).pack(side=tk.LEFT,
                                                                                                          padx=2)
        ttk.Button(btn_layout, text="◀ 1x", command=lambda: self.set_play(direction=-1, speed=1.0)).pack(side=tk.LEFT,
                                                                                                         padx=2)
        ttk.Button(btn_layout, text="◂ 0.5x", command=lambda: self.set_play(direction=-1, speed=0.5)).pack(side=tk.LEFT,
                                                                                                           padx=2)
        ttk.Button(btn_layout, text="◂ Step", command=lambda: self.step_frame(direction=-1)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_layout, text="█", command=self.stop_play).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_layout, text="▸ Step", command=lambda: self.step_frame(direction=1)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_layout, text="0.5x ▸", command=lambda: self.set_play(direction=1, speed=0.5)).pack(side=tk.LEFT,
                                                                                                          padx=2)
        ttk.Button(btn_layout, text="1x ▶", command=lambda: self.set_play(direction=1, speed=1.0)).pack(side=tk.LEFT,
                                                                                                        padx=2)
        ttk.Button(btn_layout, text="2x ▶▶", command=lambda: self.set_play(direction=1, speed=2.0)).pack(side=tk.LEFT,
                                                                                                         padx=2)

        jump_frame = ttk.Frame(btn_layout)
        jump_frame.pack(side=tk.RIGHT, padx=10)
        ttk.Button(btn_layout, text="Reset Zoom", command=self.reset_view).pack(side=tk.LEFT, padx=10)
        ttk.Label(jump_frame, text="Jump to Time (s):").pack(side=tk.LEFT, padx=2)
        jump_ent = ttk.Entry(jump_frame, textvariable=self.jump_var, width=8)
        jump_ent.pack(side=tk.LEFT, padx=2)
        jump_ent.bind("<Return>", lambda e: self.execute_time_jump())
        ttk.Button(jump_frame, text="Go", command=self.execute_time_jump, width=4).pack(side=tk.LEFT, padx=2)

    def get_current_display_settings(self):
        return {
            "pct_low": float(self.pct_low.get()),
            "pct_high": float(self.pct_high.get()),
            "f_min": float(self.f_min.get()),
            "f_max": float(self.f_max.get()),
            "t_width_seconds": float(self.t_width_seconds),
            "show_grid": bool(self.show_grid.get()),
            "show_contours": bool(self.show_contours.get()),
            "contour_count": int(self.contour_count.get()),
            "colormap_name": str(self.colormap_name.get()),
        }

    def build_display_side_controls(self, parent):
        panel = ttk.LabelFrame(parent, text="Display Controls", padding=8)
        panel.pack(side=tk.RIGHT, fill=tk.Y, padx=(8, 0))
        panel.pack_propagate(False)
        panel.configure(width=230)

        ttk.Label(panel, text="Intensity percentiles").pack(anchor=tk.W)
        self.pct_range_label = ttk.Label(panel, text="")
        self.pct_range_label.pack(anchor=tk.W, pady=(0, 2))
        self.pct_slider = RangeSlider(
            panel, 0, 100, self.pct_low, self.pct_high,
            command=self.on_display_control_changed, width=200, min_gap=0.1, decimals=2
        )
        self.pct_slider.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(panel, text="Time window, seconds").pack(anchor=tk.W)
        self.time_window_label = ttk.Label(panel, text="")
        self.time_window_label.pack(anchor=tk.W, pady=(0, 2))
        self.time_window_slider = ttk.Scale(
            panel, from_=-1.0, to=1.0, orient=tk.HORIZONTAL,
            variable=self.time_window_log, command=self.on_time_window_slider_changed
        )
        self.time_window_slider.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(panel, text="Frequency range, Hz").pack(anchor=tk.W)
        self.freq_range_label = ttk.Label(panel, text="")
        self.freq_range_label.pack(anchor=tk.W, pady=(0, 2))
        self.freq_slider = RangeSlider(
            panel, 0, 1000, self.f_min, self.f_max,
            command=self.on_display_control_changed, width=200, min_gap=1.0, decimals=0
        )
        self.freq_slider.pack(fill=tk.X, pady=(0, 10))

        ttk.Checkbutton(panel, text="Grid", variable=self.show_grid, command=self.on_display_control_changed).pack(
            anchor=tk.W, pady=(2, 2))

        contour_row = ttk.Frame(panel)
        contour_row.pack(fill=tk.X, pady=(2, 4))
        ttk.Checkbutton(contour_row, text="Contour lines", variable=self.show_contours,
                        command=self.on_display_control_changed).pack(side=tk.LEFT)

        contour_count_row = ttk.Frame(panel)
        contour_count_row.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(contour_count_row, text="Line count:").pack(side=tk.LEFT)
        contour_entry = ttk.Entry(contour_count_row, textvariable=self.contour_count, width=6)
        contour_entry.pack(side=tk.LEFT, padx=(6, 0))
        contour_entry.bind("<Return>", self.on_display_control_changed)
        contour_entry.bind("<FocusOut>", self.on_display_control_changed)

        ttk.Label(panel, text="Colour scheme").pack(anchor=tk.W)
        self.colormap_combo = ttk.Combobox(
            panel, textvariable=self.colormap_name,
            values=("viridis", "inferno", "magma", "plasma", "cividis", "turbo", "gray", "hot"),
            state="readonly", width=16
        )
        self.colormap_combo.pack(fill=tk.X, pady=(2, 10))
        self.colormap_combo.bind("<<ComboboxSelected>>", self.on_display_control_changed)

        ttk.Button(panel, text="Reset display", command=self.reset_display_controls).pack(fill=tk.X, pady=(8, 0))
        ttk.Button(panel, text="Save as default", command=self.save_display_defaults).pack(fill=tk.X, pady=(4, 0))
        self.update_display_control_labels()

    def get_percentile_limits(self):
        try:
            low = float(self.pct_low.get())
            high = float(self.pct_high.get())
        except Exception:
            low, high = 50.0, 99.99

        low = max(0.0, min(100.0, low))
        high = max(0.0, min(100.0, high))

        if high <= low:
            high = min(100.0, low + 0.1)
            low = max(0.0, high - 0.1)
        return low, high

    def get_frequency_limits(self):
        try:
            low = float(self.f_min.get())
            high = float(self.f_max.get())
        except Exception:
            low, high = 0.0, 500.0

        low = max(0.0, min(1000.0, low))
        high = max(0.0, min(1000.0, high))

        if high <= low:
            high = min(1000.0, low + 1.0)
            low = max(0.0, high - 1.0)
        return low, high

    def get_contour_level_count(self):
        try:
            count = int(self.contour_count.get())
        except Exception:
            count = 25
        return max(2, min(100, count))

    def get_colormap_name(self):
        try:
            cmap = str(self.colormap_name.get()).strip()
        except Exception:
            cmap = "viridis"
        return cmap or "viridis"

    def apply_chart_grid(self, ax):
        show = bool(self.show_grid.get())
        if show:
            ax.set_axisbelow(False)
            ax.grid(True, which="major", linewidth=0.6, alpha=0.85)
        else:
            ax.grid(False)

    def update_display_control_labels(self):
        if not hasattr(self, "pct_range_label"):
            return
        pct_low, pct_high = self.get_percentile_limits()
        f_low, f_high = self.get_frequency_limits()
        self.pct_range_label.config(text=f"{pct_low:.2f}% to {pct_high:.2f}%")
        self.freq_range_label.config(text=f"{f_low:.0f} Hz to {f_high:.0f} Hz")
        self.time_window_label.config(text=f"{self.t_width_seconds:.3f} s")
        if hasattr(self, "pct_slider"):
            self.pct_slider.redraw()
        if hasattr(self, "freq_slider"):
            self.freq_slider.redraw()

    def schedule_display_update(self, recenter=False):
        self.update_display_control_labels()
        if recenter:
            self.display_controls_recenter = True
        if getattr(self, "display_control_update_job", None) is not None:
            try:
                self.root.after_cancel(self.display_control_update_job)
            except Exception:
                pass
        self.display_control_update_job = self.root.after(75, self.apply_display_control_update)

    def apply_display_control_update(self):
        self.display_control_update_job = None
        if self.total_duration > 0:
            if self.display_controls_recenter:
                self.display_controls_recenter = False
                self.clamp_and_set_center(self.t_center.get())
            else:
                self.update_all_tabs()
        else:
            self.display_controls_recenter = False

    def on_display_control_changed(self, event=None):
        self.schedule_display_update(recenter=False)
        return None

    def on_time_window_slider_changed(self, event=None):
        try:
            log_value = float(self.time_window_log.get())
            width_seconds = 10.0 ** log_value
        except Exception:
            width_seconds = 1.0

        width_seconds = max(0.1, min(10.0, width_seconds))
        self.t_width_seconds = width_seconds
        self.t_width_str.set(f"{width_seconds:.3f}")
        self.schedule_display_update(recenter=True)
        return None

    def sync_display_controls(self):
        try:
            width_seconds = max(0.1, min(10.0, float(self.t_width_seconds)))
            self.time_window_log.set(np.log10(width_seconds))
        except Exception:
            self.time_window_log.set(0.0)
        self.update_display_control_labels()

    def reset_display_controls(self):
        defaults = getattr(self, "display_defaults", None)
        if defaults is None:
            defaults = {
                "pct_low": 50.0, "pct_high": 99.99, "f_min": 0.0, "f_max": 500.0,
                "t_width_seconds": 0.5, "show_grid": False, "show_contours": False,
                "contour_count": 25, "colormap_name": "inferno",
            }

        self.pct_low.set(defaults["pct_low"])
        self.pct_high.set(defaults["pct_high"])
        self.f_min.set(defaults["f_min"])
        self.f_max.set(defaults["f_max"])
        self.t_width_seconds = defaults["t_width_seconds"]
        self.t_width_str.set(f"{self.t_width_seconds:.3f}")
        self.time_window_log.set(np.log10(self.t_width_seconds))
        self.show_grid.set(defaults["show_grid"])
        self.show_contours.set(defaults["show_contours"])
        self.contour_count.set(defaults["contour_count"])
        self.colormap_name.set(defaults["colormap_name"])
        self.schedule_display_update(recenter=True)

    def start_zoom_drag(self, event):
        self.zoom_dragging = True
        active_tab_name = self.notebook.tab(self.notebook.select(), "text")
        if active_tab_name not in self.tabs:
            return

        tab = self.tabs[active_tab_name]
        ax = tab["ax"]
        canvas = tab["canvas"]
        canvas_height = canvas.get_tk_widget().winfo_height()
        inv = ax.transData.inverted()
        self.zoom_start_data = inv.transform((event.x, canvas_height - event.y))

    def zoom_drag_motion(self, event):
        if not getattr(self, "zoom_dragging", False):
            return
        active_tab_name = self.notebook.tab(self.notebook.select(), "text")
        if active_tab_name not in self.tabs:
            return

        tab = self.tabs[active_tab_name]
        ax = tab["ax"]
        canvas = tab["canvas"]

        canvas_height = canvas.get_tk_widget().winfo_height()
        inv = ax.transData.inverted()

        x0, y0 = self.zoom_start_data
        x1, y1 = inv.transform((event.x, canvas_height - event.y))
        xmin, xmax = sorted([x0, x1])
        ymin, ymax = sorted([y0, y1])

        if hasattr(self, "zoom_rect") and self.zoom_rect is not None:
            try:
                self.zoom_rect.remove()
            except Exception:
                pass

        self.zoom_rect = matplotlib.patches.Rectangle(
            (xmin, ymin), xmax - xmin, ymax - ymin,
            fill=False, linewidth=1.5, linestyle="--", edgecolor="white")

        ax.add_patch(self.zoom_rect)
        canvas.draw_idle()

    def end_zoom_drag(self, event):
        if not getattr(self, "zoom_dragging", False):
            return
        self.zoom_dragging = False
        active_tab_name = self.notebook.tab(self.notebook.select(), "text")
        if active_tab_name not in self.tabs:
            return

        tab = self.tabs[active_tab_name]
        ax = tab["ax"]
        canvas = tab["canvas"]
        canvas_height = canvas.get_tk_widget().winfo_height()
        inv = ax.transData.inverted()

        x0, y0 = self.zoom_start_data
        x1, y1 = inv.transform((event.x, canvas_height - event.y))
        xmin, xmax = sorted([x0, x1])
        ymin, ymax = sorted([y0, y1])

        if hasattr(self, "zoom_rect") and self.zoom_rect is not None:
            try:
                self.zoom_rect.remove()
            except Exception:
                pass
            self.zoom_rect = None

        if abs(xmax - xmin) < 0.01 or abs(ymax - ymin) < 1.0:
            self.update_all_tabs()
            return

        self.t_width_seconds = xmax - xmin
        self.t_width_str.set(f"{self.t_width_seconds:.4f}")
        self.t_center.set((xmin + xmax) / 2.0)
        self.f_min.set(max(0.0, ymin))
        self.f_max.set(max(self.f_min.get() + 1.0, ymax))
        self.sync_display_controls()
        self.clamp_and_set_center(self.t_center.get())

    def reset_view(self):
        self.stop_play()
        self.t_width_seconds = 0.5
        self.t_width_str.set("0.5")
        self.f_min.set(0)
        self.f_max.set(500)
        self.pct_low.set(50.0)
        self.pct_high.set(99.99)
        self.sync_display_controls()

        if self.total_duration > 0:
            self.clamp_and_set_center(self.t_center.get())
        else:
            self.update_all_tabs()

    def execute_time_jump(self):
        if self.total_duration == 0: return
        try:
            target = float(self.jump_var.get().strip())
            self.stop_play()
            self.clamp_and_set_center(target)
        except ValueError:
            messagebox.showerror("Format Error", "Please provide a valid numeric float coordinate value.")

    def pull_exact_detector_suite(self, event_name, duration_val, rate_val, gps_merger, files_dict):
        try:
            self.root.after(0, self.clear_sky_solution)
            for det in ["H1", "L1", "V1"]:
                self.root.after(0, lambda d=det: self.clear_detector(d))

            self.current_file_duration = f"{duration_val} Seconds"
            self.current_event_name = event_name

            target_hz = 16384 if "16" in rate_val else 4096
            self.current_file_rate = f"{target_hz} Hz"
            self.cached_target_gps = gps_merger

            if not files_dict:
                self.root.after(0, lambda: messagebox.showwarning(
                    "Empty Matrix", f"No verified URLs found in the cache map for {event_name}."
                ))
                return

            for det_key, data in files_dict.items():
                exact_url = data["url"]
                threading.Thread(
                    target=self.download_and_ingest,
                    args=(exact_url, det_key, False),
                    daemon=True).start()

        except Exception as e:
            self.root.after(0, lambda err=e: messagebox.showerror("Pipeline Loader Error", str(err)))

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
        if not self.is_playing:
            return
        step_increment = self.play_direction * self.t_width_seconds * 0.06 * self.play_speed
        new_center = self.t_center.get() + step_increment

        if (new_center <= self.t_width_seconds / 2.0 or new_center >= self.total_duration - self.t_width_seconds / 2.0):
            self.is_playing = False
            return

        self.clamp_and_set_center(new_center)
        self.root.after(int(1000 / self.base_fps), self.play_step)

    def open_sky_link(self):
        if not self.sky_link_url:
            self.recompute_sky_from_offsets(show_errors=False)
        if self.sky_link_url:
            webbrowser.open(self.sky_link_url)

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

        cat_lbl = ttk.Label(progress_frame, text="Cache Status: Waiting...", font=('Helvetica', 9, 'italic'),
                            foreground="blue")
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
            from gwosc.locate import get_event_urls
            if not isinstance(name, str) or len(name) < 2:
                return []

            all_urls = []
            for target_hz in [4096, 16384]:
                try:
                    urls = get_event_urls(name, sample_rate=target_hz)
                    all_urls.extend(urls)
                except Exception:
                    pass

            if not all_urls:
                return []

            variants = {}
            for url in all_urls:
                if not url.endswith('.hdf5'): continue
                filename = url.split("/")[-1]
                fname_upper = filename.upper()

                det = 'H1' if 'H-H1' in fname_upper else (
                    'L1' if 'L-L1' in fname_upper else ('V1' if 'V-V1' in fname_upper else 'Unknown'))
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

                with ThreadPoolExecutor(max_workers=20) as executor:
                    future_to_event = {executor.submit(process_single_event, name): name for name in all_events}

                    for future in as_completed(future_to_event):
                        name = future_to_event[future]
                        completed_count += 1

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
                    cat_lbl.config(text=f"Complete! Indexed {len(self.master_records)} valid tracks.",
                                   foreground="green")
                ])

            except Exception as e:
                self.root.after(0, lambda err=e: [
                    messagebox.showerror("API Connection Error", f"Failed to parse catalog: {str(err)}"),
                    cat_lbl.config(text="Download failed.", foreground="red")
                ])

        def apply_filter(*args):
            if not self.tree.winfo_exists():
                return
            for i in self.tree.get_children(): self.tree.delete(i)
            query_string = search_var.get().strip().upper()

            for item in self.master_records:
                if not query_string or query_string in item['event'].upper():
                    self.tree.insert('', tk.END,
                                     values=(item['event'], item['duration'], item['rate'], item['filenames_str']))

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

            record = next((r for r in self.master_records if str(r['event']) == chosen_event
                           and str(r['duration']) == chosen_duration and str(r['rate']) == chosen_rate), None)

            if not record: return

            gps_merger = "N/A"
            try:
                from gwosc.datasets import event_gps
                gps_merger = event_gps(chosen_event)
            except Exception:
                pass

            win.destroy()
            self.global_status.config(text=f"Initiating direct mapped download for {chosen_event}...")

            self.add_recent_catalog_entry(
                chosen_event,
                record["duration"],
                record["rate"],
                gps_merger,
                record["files_dict"]
            )

            threading.Thread(target=self.pull_exact_detector_suite,
                             args=(chosen_event, record['duration'], record['rate'], gps_merger, record['files_dict']),
                             daemon=True).start()

        btn = ttk.Button(win, text="Download Selected Stream Matrix", command=load_selected_record)
        btn.pack(fill=tk.X, padx=10, pady=10)
        self.tree.bind("<Double-1>", lambda event: load_selected_record())
        query_api(force_download=False)

    def clear_detector(self, det):
        self.detectors[det]["raw"] = None
        self.detectors[det]["whitened"] = None
        self.detectors[det]["Sxx"] = None
        self.detectors[det]["t"] = None
        self.detectors[det]["f"] = None
        self.detectors[det]["loaded"] = False

        self.current_filenames[det] = "N/A"

        if hasattr(self, "file_labels"):
            self.file_labels[det].config(text=f"{det} File: N/A")

        if "status_lbl" in self.detectors[det]:
            self.detectors[det]["status_lbl"].config(
                text="Status: Empty",
                foreground="gray"
            )

        if det in self.tabs:
            ax = self.tabs[det]["ax"]
            canvas = self.tabs[det]["canvas"]
            ax.clear()
            canvas.draw_idle()

    def clear_sky_solution(self):
        self.sky_link_url = None
        self.last_sky_ra_deg = None
        self.last_sky_dec_deg = None
        if hasattr(self, "sky_link_text"):
            self.sky_link_text.set("N/A")
        if hasattr(self, "sky_result_text"):
            self.sky_result_text.set("Sky: Waiting...")

    def on_slider_move(self, val):
        if self.total_duration == 0 or self.is_dragging: return
        self.clamp_and_set_center(float(val), update_slider_ui=False)

    def clamp_and_set_center(self, target_center, update_slider_ui=True):
        if self.total_duration <= 0:
            return

        if self.total_duration <= self.t_width_seconds:
            clamped = self.total_duration / 2.0
        else:
            min_c = self.t_width_seconds / 2.0
            max_c = self.total_duration - min_c
            clamped = max(min_c, min(target_center, max_c))

        self.t_center.set(clamped)
        self.time_lbl.config(text=f"Current Window Center: {clamped:.2f}s")

        if update_slider_ui:
            self.timeline_slider.set(clamped)

        self.update_all_tabs()

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
            self.root.after(0,
                            lambda: self.global_status.config(text=f"Re-generating 2D Spectrogram maps for {det}..."))

            self.compute_complete_spectrogram(det)

            self.root.after(0, lambda: self.update_pbar(det, 100))
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, self.update_all_tabs)
            self.root.after(0, lambda: self.global_status.config(text="System: Idle"))
        except Exception as e:
            err = str(e)
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, lambda err=err: messagebox.showerror("Re-Whitening Error", err))

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

        ttk.Label(win, text="Window Width (s or fraction e.g., 1/4):").grid(row=0, column=0, padx=10, pady=5,
                                                                            sticky=tk.W)
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
                self.sync_display_controls()
                if self.total_duration > 0:
                    self.clamp_and_set_center(self.t_center.get())
                else:
                    self.update_all_tabs()
            except Exception:
                messagebox.showerror("Parsing Error",
                                     "Invalid input format detected. Please verify numerical configurations.")

        ttk.Button(win, text="Apply Changes", command=save_and_eval).grid(row=5, columnspan=2, pady=10)

    def render_canvas_frame(self, ax, canvas, Sxx, t, f, t_start, t_end):
        ax.clear()
        f_low, f_high = self.get_frequency_limits()
        pct_low, pct_high = self.get_percentile_limits()

        t_mask = (t >= t_start) & (t <= t_end)
        f_mask = (f >= f_low) & (f <= f_high)

        if not np.any(t_mask) or not np.any(f_mask):
            canvas.draw_idle()
            return

        Sxx_sub = Sxx[np.ix_(f_mask, t_mask)]
        vmin = np.percentile(Sxx_sub, pct_low)
        vmax = np.percentile(Sxx_sub, pct_high)

        if vmax <= vmin + 1e-12:
            vmax = vmin + 1e-9

        contour_count = self.get_contour_level_count()
        levels = np.linspace(vmin, vmax, contour_count)
        cmap = self.get_colormap_name()

        try:
            ax.pcolormesh(
                t[t_mask],
                f[f_mask],
                Sxx_sub,
                shading='gouraud',
                cmap=cmap,
                vmin=vmin,
                vmax=vmax
            )

            if self.show_contours.get() and Sxx_sub.shape[0] >= 2 and Sxx_sub.shape[1] >= 2:
                ax.contour(
                    t[t_mask],
                    f[f_mask],
                    Sxx_sub,
                    levels=levels,
                    colors="white",
                    linewidths=0.25, alpha=0.5,
                    zorder=2
                )
        except Exception:
            ax.pcolormesh(
                t[t_mask],
                f[f_mask],
                Sxx_sub,
                shading='gouraud',
                cmap='viridis'
            )

        ax.set_xlim(t_start, t_end)
        ax.set_ylim(f_low, f_high)
        self.apply_chart_grid(ax)
        canvas.draw_idle()

    def start_drag(self, event):
        self.stop_play()
        self.is_dragging = True
        self.drag_start_x = event.x

    def drag_motion(self, event):
        if not self.is_dragging or self.total_duration == 0:
            return

        active_tab_name = self.notebook.tab(self.notebook.select(), "text")
        if active_tab_name not in self.tabs: return
        ax = self.tabs[active_tab_name]['ax']

        inv = ax.transData.inverted()
        current_data_x, _ = inv.transform((event.x, event.y))
        start_data_x, _ = inv.transform((self.drag_start_x, 0))

        data_shift = start_data_x - current_data_x
        new_center = self.t_center.get() + data_shift

        self.clamp_and_set_center(new_center)
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
            threading.Thread(
                target=self.process_pipeline_worker,
                args=(filename, det),
                kwargs={"add_to_recent": True},
                daemon=True).start()

    def download_and_ingest(self, url, det, add_to_recent=False):
        try:
            local_name = url.split("/")[-1]

            if not os.path.exists(local_name):
                self.root.after(0, lambda: self.global_status.config(text=f"Contacting GWOSC server for {det}..."))
                with requests.get(url, stream=True, timeout=20) as r:
                    r.raise_for_status()
                    total_size = int(r.headers.get('content-length', 0))
                    size_mb = total_size / (1024 * 1024)

                    self.root.after(0, lambda: self.global_status.config(
                        text=f"Downloading {det} stream ({size_mb:.1f} MB)..."))
                    with open(local_name, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)

            self.root.after(0, lambda: self.global_status.config(text=f"Download complete. Ingesting {det} array..."))

            self.process_pipeline_worker(
                local_name,
                det,
                add_to_recent=add_to_recent
            )
        except requests.exceptions.Timeout:
            self.root.after(0, lambda: self.global_status.config(text="System: Idle (Download Timeout)"))
            self.root.after(0, lambda: messagebox.showerror("Network Timeout",
                                                            f"The server took too long to respond for {det}. It might be under heavy load."))
        except Exception as e:
            err = str(e)
            self.root.after(0, lambda: self.global_status.config(text="System: Idle (Download Failed)"))
            self.root.after(0, lambda err=err, det=det: messagebox.showerror("Network Crash",
                                                                             f"Failed to fetch {det}:\n{err}"))

    def show_pbar(self, det):
        self.tabs[det]['pbar'].pack(fill=tk.X, padx=20, pady=5, before=self.tabs[det]['canvas'].get_tk_widget())
        self.tabs[det]['pbar']['value'] = 0

    def update_pbar(self, det, val):
        self.tabs[det]['pbar']['value'] = val

    def hide_pbar(self, det):
        self.tabs[det]['pbar'].pack_forget()

    def process_pipeline_worker(self, filepath, det, add_to_recent=True):
        try:
            self.root.after(0, lambda: self.notebook.select(self.notebook.tabs()[['H1', 'L1', 'V1'].index(det)]))
            self.root.after(0, lambda: self.show_pbar(det))

            filename_only = os.path.basename(filepath)
            self.current_filenames[det] = filename_only

            if add_to_recent:
                self.root.after(0, lambda fp=filepath, d=det: self.add_recent_file(fp, d))

            if hasattr(self, "detector_offsets_ms"):
                for d in ["L1", "V1"]:
                    self.detector_offsets_ms[d].set(0)

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

            self.detectors[det]['whitened'] = self.whiten_by_intervals(self.detectors[det]['raw'], self.fs,
                                                                       self.whiten_interval.get(), det)

            self.root.after(0, lambda: self.update_pbar(det, 75))
            self.root.after(0, lambda: self.global_status.config(
                text=f"Executing 2D Spectrogram & Laplacian Matrix loops for {det}..."))

            self.compute_complete_spectrogram(det)

            self.root.after(0, lambda: self.update_pbar(det, 100))
            self.detectors[det]['loaded'] = True

            target_time = self.total_duration / 2.0
            try:
                if self.cached_target_gps != "N/A" and float(self.cached_target_gps) > 0:
                    potential_offset = float(self.cached_target_gps) - float(gps_start)
                    if 0 <= potential_offset <= self.total_duration:
                        target_time = round(potential_offset, 2)
                        if self.current_event_name == "GW150914":
                            target_time = 15.4
            except Exception:
                pass

            if self.timeline_slider.cget('to') != self.total_duration:
                self.root.after(0, lambda: self.timeline_slider.config(to=self.total_duration))
            self.root.after(100, lambda: self.clamp_and_set_center(target_time))

            self.root.after(0, lambda: self.detectors[det]['status_lbl'].config(text="Status: Loaded & Cached",
                                                                                foreground="green"))
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, self.update_all_tabs)
            self.root.after(0, lambda: self.global_status.config(text="System: Idle"))

        except Exception as e:
            err = str(e)
            self.root.after(0, lambda: self.hide_pbar(det))
            self.root.after(0, lambda err=err: messagebox.showerror("Pipeline Failure", err))

    def whiten_by_intervals(self, strain, fs, interval_sec, det):
        N = len(strain)
        chunk_size = int(interval_sec * fs)
        whitened = np.zeros_like(strain)
        window_sum = np.zeros_like(strain)

        t = np.linspace(-1, 1, chunk_size)
        with np.errstate(divide='ignore', invalid='ignore'):
            win = np.where(np.abs(t) < 1, np.exp(t ** 4 / (t ** 2 - 1)), 0.0)

        starts = list(range(0, N, chunk_size // 2))
        total_steps = len(starts)

        for idx, start in enumerate(starts):
            end = start + chunk_size
            if end > N: break

            slice_data = strain[start:end] * win
            spec = np.fft.rfft(slice_data)
            mag = np.abs(spec)
            freqs = np.fft.rfftfreq(chunk_size, 1 / fs)

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

    def get_detector_offset_ms(self, det):
        try:
            if hasattr(self, "offset_spinboxes") and det in self.offset_spinboxes:
                return float(self.offset_spinboxes[det].get())
            return float(self.detector_offsets_ms[det].get())
        except Exception:
            return 0.0

    def on_offset_spinbox_changed(self, event=None):
        if self.total_duration > 0:
            self.update_all_tabs()
        self.schedule_sky_recompute()
        return None

    def schedule_sky_recompute(self, event=None):
        if getattr(self, "sky_recompute_job", None) is not None:
            try:
                self.root.after_cancel(self.sky_recompute_job)
            except Exception:
                pass
        self.sky_recompute_job = self.root.after(200, lambda: self.recompute_sky_from_offsets(show_errors=False))

    def trigger_frame_correlation(self):
        self.lbl_offset_result.config(text="Calculating...", foreground="orange")
        self.root.update_idletasks()

        try:
            if not self.detectors["H1"]["loaded"]:
                raise ValueError("H1 must be loaded as the reference detector.")

            h1_data = self.detectors["H1"]["whitened"]
            t_center = self.t_center.get()
            t_width = self.t_width_seconds
            t_start = max(0.0, t_center - t_width / 2)
            t_end = min(self.total_duration, t_center + t_width / 2)

            idx_min = int(t_start * self.fs)
            idx_max = int(t_end * self.fs)
            h1_slice = h1_data[idx_min:idx_max]
            results = []

            for det in ["L1", "V1"]:
                if not self.detectors[det]["loaded"]:
                    continue

                other_data = self.detectors[det]["whitened"]
                other_slice = other_data[idx_min:idx_max]

                offset_ms, lead_text = self.calculate_delay_xcorr(
                    h1_slice, other_slice, self.fs,
                    max(35.0, self.f_min.get()), self.f_max.get(), max_delay_ms=30.0
                )

                signed_ms = -offset_ms if "leads H1" in lead_text else offset_ms
                rounded_ms = round(signed_ms, 1)

                if hasattr(self, "detector_offsets_ms"):
                    self.detector_offsets_ms[det].set(rounded_ms)

                results.append(f"{det}: {rounded_ms:+.1f} ms")

            if not results:
                raise ValueError("No L1 or V1 detector loaded.")

            self.lbl_offset_result.config(text="Offsets: " + " | ".join(results), foreground="green")
            self.update_all_tabs()
            self.schedule_sky_recompute()

        except Exception as e:
            self.lbl_offset_result.config(text=f"Error: {str(e)}", foreground="red")
            print(f"Correlation Error: {e}")

    def gps_to_gmst(self, gps_time):
        jd = 2444244.5 + float(gps_time) / 86400.0
        d = jd - 2451545.0
        gmst_deg = 280.46061837 + 360.98564736629 * d
        return np.deg2rad(gmst_deg % 360.0)

    def radec_to_ecef_unit(self, ra, dec, gmst):
        x = np.cos(dec) * np.cos(ra)
        y = np.cos(dec) * np.sin(ra)
        z = np.sin(dec)

        x_ecef = np.cos(gmst) * x + np.sin(gmst) * y
        y_ecef = -np.sin(gmst) * x + np.cos(gmst) * y
        z_ecef = z
        return np.array([x_ecef, y_ecef, z_ecef])

    def estimate_sky_location(self):
        self.recompute_sky_from_offsets(show_errors=True)

    def recompute_sky_from_offsets(self, show_errors=False):
        self.sky_recompute_job = None
        try:
            if not (self.detectors["H1"]["loaded"] and self.detectors["L1"]["loaded"] and self.detectors["V1"][
                "loaded"]):
                if show_errors:
                    raise ValueError("H1, L1 and V1 must all be loaded.")
                return None

            gps_time = float(self.cached_target_gps)
            observed = {
                "L1": self.get_detector_offset_ms("L1"),
                "V1": self.get_detector_offset_ms("V1"),
            }
            gmst = self.gps_to_gmst(gps_time)

            best_err = np.inf
            best_ra = None
            best_dec = None
            best_model = None

            for dec_deg in np.linspace(-90, 90, 181):
                dec = np.deg2rad(dec_deg)
                for ra_deg in np.linspace(0, 360, 361, endpoint=False):
                    ra = np.deg2rad(ra_deg)
                    n = self.radec_to_ecef_unit(ra, dec, gmst)
                    model = {}

                    for det in ["L1", "V1"]:
                        baseline = DETECTOR_ECEF[det] - DETECTOR_ECEF["H1"]
                        dt = -np.dot(n, baseline) / C_LIGHT
                        model[det] = dt * 1000.0

                    err = (model["L1"] - observed["L1"]) ** 2 + (model["V1"] - observed["V1"]) ** 2
                    if err < best_err:
                        best_err = err
                        best_ra = ra_deg
                        best_dec = dec_deg
                        best_model = dict(model)

            if best_ra is None or best_dec is None or best_model is None:
                raise ValueError("Could not solve sky position from current offsets.")

            self.last_sky_ra_deg = float(best_ra)
            self.last_sky_dec_deg = float(best_dec)
            self.sky_link_url = self.build_sky_url(self.last_sky_ra_deg, self.last_sky_dec_deg)

            result_text = (
                f"RA={best_ra:.1f} deg, DEC={best_dec:+.1f} deg, "
                f"L1={best_model['L1']:+.1f} ms, V1={best_model['V1']:+.1f} ms"
            )

            self.sky_result_text.set(result_text)
            if hasattr(self, "sky_link_text"):
                self.sky_link_text.set("Open Link")
            if hasattr(self, "lbl_offset_result"):
                self.lbl_offset_result.config(foreground="purple")
            return self.last_sky_ra_deg, self.last_sky_dec_deg

        except Exception as e:
            if show_errors:
                self.lbl_offset_result.config(text=f"Sky Error: {str(e)}", foreground="red")
                print(f"Sky Location Error: {e}")
            return None

    def build_sky_url(self, ra_deg, dec_deg):
        try:
            event_dt = self.gps_to_utc_datetime(float(self.cached_target_gps))
        except Exception:
            event_dt = datetime.now(timezone.utc)

        params = {
            "no_cookie": 1, "latitude": 51.51, "longitude": -0.13, "timezone": "0.00",
            "year": event_dt.year, "month": event_dt.month, "day": event_dt.day,
            "hour": event_dt.hour, "min": event_dt.minute, "PLlimitmag": 2, "zoom": 5,
            "ra": f"{ra_deg / 15.0:.5f}", "dec": f"{dec_deg:.5f}",
        }
        return "https://in-the-sky.org/skymap.php?" + urlencode(params)

    def compute_complete_spectrogram(self, det):
        nperseg = self.nperseg.get()
        nfft = self.nfft.get()
        noverlap = int(nperseg * (self.overlap_pct.get() / 100.0))
        t_win = np.linspace(-1, 1, nperseg)
        with np.errstate(divide='ignore', invalid='ignore'):
            win_seg = np.where(np.abs(t_win) < 1, np.exp(t_win ** 4 / (t_win ** 2 - 1)), 0.0)
        f, t, Sxx_power = spectrogram(self.detectors[det]['whitened'], self.fs,
                                      window=win_seg, nperseg=nperseg, nfft=nfft, noverlap=noverlap, scaling='spectrum')
        Sxx = np.sqrt(Sxx_power)
        alpha = np.sum(t_win * t_win * win_seg) / np.sum(win_seg)
        laplacian_f = np.zeros_like(Sxx)
        laplacian_f[1:-1, :] = Sxx[2:, :] - 2 * Sxx[1:-1, :] + Sxx[0:-2, :]
        Sxx_sharp = Sxx - (alpha / (nfft / nperseg)) * laplacian_f
        self.detectors[det]['Sxx'] = np.maximum(np.nan_to_num(Sxx_sharp, nan=1e-20), 1e-20)
        self.detectors[det]['t'] = t
        self.detectors[det]['f'] = f

    def update_all_tabs(self):
        if self.total_duration == 0: return
        t_center = self.t_center.get()
        t_width = self.t_width_seconds
        t_start, t_end = max(0.0, t_center - t_width / 2), min(self.total_duration, t_center + t_width / 2)
        for det in ['H1', 'L1', 'V1']:
            if self.detectors[det]['loaded']:
                self.render_canvas_frame(self.tabs[det]['ax'], self.tabs[det]['canvas'], self.detectors[det]['Sxx'],
                                         self.detectors[det]['t'],
                                         self.detectors[det]['f'], t_start, t_end)
        self.render_joint_correlation(t_start, t_end)
        self.render_whitened_waveforms(t_start, t_end)

    def render_joint_correlation(self, t_start, t_end):
        ax = self.tabs['Joint Correlation']['ax']
        canvas = self.tabs['Joint Correlation']['canvas']
        ax.clear()

        active_matrices = [det for det in ['H1', 'L1', 'V1'] if
                           self.detectors[det]['loaded'] and self.detectors[det]['active_corr'].get()]

        if not active_matrices:
            canvas.draw_idle()
            return

        f_low, f_high = self.get_frequency_limits()
        pct_low, pct_high = self.get_percentile_limits()

        ref = active_matrices[0]
        ref_t = self.detectors[ref]['t']
        ref_f = self.detectors[ref]['f']

        t_mask = (ref_t >= t_start) & (ref_t <= t_end)
        f_mask = (ref_f >= f_low) & (ref_f <= f_high)

        if not np.any(t_mask) or not np.any(f_mask):
            canvas.draw_idle()
            return

        ref_t_window = ref_t[t_mask]
        ref_f_window = ref_f[f_mask]

        joint_product = None

        for det in active_matrices:
            Sxx = self.detectors[det]['Sxx']
            det_t = self.detectors[det]['t']
            det_f = self.detectors[det]['f']

            offset_sec = self.get_detector_offset_ms(det) / 1000.0
            sample_t = ref_t_window - offset_sec

            det_f_mask = (det_f >= f_low) & (det_f <= f_high)
            S_freq = Sxx[det_f_mask, :]

            if S_freq.shape[0] == 0:
                continue

            shifted_rows = []
            for row in S_freq:
                shifted = np.interp(sample_t, det_t, row, left=0.0, right=0.0)
                shifted_rows.append(shifted)

            S_shifted = np.array(shifted_rows)
            norm_S = S_shifted / (np.max(S_shifted) + 1e-20)

            if joint_product is None:
                joint_product = norm_S
            else:
                min_h = min(joint_product.shape[0], norm_S.shape[0])
                min_w = min(joint_product.shape[1], norm_S.shape[1])
                joint_product = joint_product[:min_h, :min_w] * norm_S[:min_h, :min_w]

        if joint_product is None:
            canvas.draw_idle()
            return

        vmin = np.percentile(joint_product, pct_low)
        vmax = np.percentile(joint_product, pct_high)
        if vmax <= vmin + 1e-12:
            vmax = vmin + 1e-9

        contour_count = self.get_contour_level_count()
        levels = np.linspace(vmin, vmax, contour_count)
        cmap = self.get_colormap_name()

        ax.pcolormesh(
            ref_t_window, ref_f_window[:joint_product.shape[0]], joint_product,
            shading='gouraud', cmap=cmap, vmin=vmin, vmax=vmax
        )

        if self.show_contours.get() and joint_product.shape[0] >= 2 and joint_product.shape[1] >= 2:
            try:
                ax.contour(
                    ref_t_window, ref_f_window[:joint_product.shape[0]], joint_product,
                    levels=levels, colors="white", linewidths=0.25, alpha=0.5, zorder=2
                )
            except Exception:
                pass

        ax.set_xlim(t_start, t_end)
        ax.set_ylim(f_low, f_high)
        self.apply_chart_grid(ax)

        title_parts = [f"{det} {self.get_detector_offset_ms(det):+.1f} ms" for det in active_matrices if det != "H1"]
        ax.set_title("Joint Correlation " + " | ".join(title_parts))
        canvas.draw_idle()

    def render_whitened_waveforms(self, t_start, t_end):
        if "Whitened Waveforms" not in self.tabs or not hasattr(self, 'fs') or self.fs is None:
            return

        axs = self.tabs["Whitened Waveforms"]['ax']
        canvas = self.tabs["Whitened Waveforms"]['canvas']

        # Clear existing lines
        for ax in axs:
            ax.clear()

        # 1. Get frequency limits from the UI controls
        f_low, f_high = self.get_frequency_limits()

        # Ensure valid bandpass frequencies (avoiding 0 Hz and staying below Nyquist)
        nyquist = self.fs / 2.0
        low = max(1.0, float(f_low))
        high = min(float(f_high), nyquist - 1.0)

        # Create the filter if the range is valid
        apply_filter = False
        if high > low:
            apply_filter = True
            sos = signal.butter(4, [low, high], btype='bandpass', fs=self.fs, output='sos')

        colors = {'H1': 'red', 'L1': 'green', 'V1': 'purple'}
        detectors = ['H1', 'L1', 'V1']
        combined_ax = axs[3]

        for i, det in enumerate(detectors):
            if self.detectors[det]['loaded'] and self.detectors[det]['whitened'] is not None:
                # Calculate exact slice indices
                idx_start = int(max(0, t_start * self.fs))
                idx_end = int(min(len(self.detectors[det]['whitened']), t_end * self.fs))

                if idx_start >= idx_end:
                    continue

                # Generate exact time array and get the data slice
                t_arr = np.arange(idx_start, idx_end) / self.fs
                data = self.detectors[det]['whitened'][idx_start:idx_end]

                # 2. Apply the zero-phase bandpass filter
                if apply_filter and len(data) > 33:
                    try:
                        data = signal.sosfiltfilt(sos, data)
                    except ValueError:
                        pass

                        # 3. Normalize the data strictly to a [-1, 1] scale based on the max absolute peak in the window
                max_amp = np.max(np.abs(data))
                if max_amp > 1e-20:
                    data = data / max_amp

                # 4. Plot individual chart
                axs[i].plot(t_arr, data, color=colors[det], label=det, linewidth=0.8)
                axs[i].set_ylabel("Norm")
                axs[i].set_ylim(-1.1, 1.1)
                axs[i].grid(True, alpha=0.5)
                axs[i].legend(loc="upper left")

                # 5. Plot combined chart with offset alignment AND potential polarity flip
                offset_sec = self.get_detector_offset_ms(det) / 1000.0
                shifted_t_arr = t_arr - offset_sec

                # Check if the user has requested to invert this specific signal
                is_flipped = getattr(self, 'signal_flips', {}).get(det, tk.BooleanVar(value=False)).get()
                plot_data = -data if is_flipped else data
                flip_text = " [Inverted]" if is_flipped else ""

                combined_ax.plot(
                    shifted_t_arr,
                    plot_data,
                    color=colors[det],
                    label=f"{det} ({offset_sec * 1000:+.1f}ms){flip_text}",
                    linewidth=0.6,
                    alpha=0.8
                )

        # Format the combined axis
        filter_text = f"Aligned Norm (BP: {low:.0f}-{high:.0f} Hz)" if apply_filter else "Aligned Norm"
        combined_ax.set_ylabel(filter_text)
        combined_ax.set_ylim(-1.1, 1.1)
        combined_ax.set_xlabel("Time (s)")
        combined_ax.grid(True, alpha=0.5)
        combined_ax.set_xlim(t_start, t_end)

        # Avoid empty legend warnings
        handles, labels = combined_ax.get_legend_handles_labels()
        if handles:
            combined_ax.legend(loc="upper left")

        canvas.draw_idle()
        
    def calculate_delay_xcorr(self, h1_slice, l1_slice, fs, f_min, f_max, max_delay_ms=15.0):
        h1 = np.asarray(h1_slice, dtype=float)
        l1 = np.asarray(l1_slice, dtype=float)

        n = min(len(h1), len(l1))
        h1 = h1[:n]
        l1 = l1[:n]

        if n < 32:
            raise ValueError("Frame too short")

        low = max(20.0, float(f_min))
        high = min(float(f_max), fs / 2.0 - 1.0)

        if high <= low:
            raise ValueError("Invalid bandpass range")

        sos = signal.butter(4, [low, high], btype="bandpass", fs=fs, output="sos")

        h1 = signal.sosfiltfilt(sos, h1)
        l1 = signal.sosfiltfilt(sos, l1)

        h1 -= np.mean(h1)
        l1 -= np.mean(l1)

        win = signal.windows.tukey(n, alpha=0.25)
        h1 *= win
        l1 *= win

        h1_norm = np.linalg.norm(h1)
        l1_norm = np.linalg.norm(l1)

        if h1_norm < 1e-20 or l1_norm < 1e-20:
            raise ValueError("Signal energy too small")

        h1 /= h1_norm
        l1 /= l1_norm

        corr = signal.correlate(l1, h1, mode="full", method="fft")
        lags = signal.correlation_lags(len(l1), len(h1), mode="full")

        max_lag = int((max_delay_ms / 1000.0) * fs)
        mask = np.abs(lags) <= max_lag

        corr_m = corr[mask]
        lags_m = lags[mask]

        peak = np.argmax(np.abs(corr_m))
        lag_samples = float(lags_m[peak])

        if 0 < peak < len(corr_m) - 1:
            y0, y1, y2 = corr_m[peak - 1], corr_m[peak], corr_m[peak + 1]
            denom = y0 - 2 * y1 + y2
            if abs(denom) > 1e-20:
                lag_samples += 0.5 * (y0 - y2) / denom

        tau_ms = 1000.0 * lag_samples / fs
        lead_text = "other leads H1" if tau_ms < 0 else "other lags H1"

        return abs(tau_ms), lead_text

    def gps_to_utc_datetime(self, gps_time):
        gps_epoch = datetime(1980, 1, 6, tzinfo=timezone.utc)
        leap_seconds = 18
        return gps_epoch + timedelta(seconds=float(gps_time) - leap_seconds)


if __name__ == "__main__":
    root = tk.Tk(className="GWExplorer")
    app = GWExplorerApp(root)
    root.mainloop()
