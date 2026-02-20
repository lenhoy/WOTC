import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox

def plot_tau_comparison(file_paths=None):
    # --- Configuration Defaults ---
    DEFAULT_TITLE = "Thrust Sum Comparison"
    DEFAULT_X_LABEL = "Time [s]"
    DEFAULT_Y_LABEL = "Force [N]"
    
    # Fonts
    FONT_SIZE = 12
    plt.rcParams.update({'font.size': FONT_SIZE})

    SHOW_GRID = True
    
    # --- 1. Handle Input (Variable Number of Files) ---
    selected_files = []
    
    # Check if files provided via CLI
    if file_paths and len(file_paths) > 0:
        selected_files = file_paths
    else:
        # Create hidden root window for dialogs
        root = tk.Tk()
        root.withdraw() 
        
        # Ask how many series
        num_series_input = simpledialog.askinteger("Input", "How many series do you want to compare?", minvalue=1, maxvalue=20)
        
        if not num_series_input:
            print("No number entered or cancelled.")
            return

        print(f"Selecting {num_series_input} files...")
        
        for i in range(num_series_input):
            f = filedialog.askopenfilename(
                title=f"Select Data File for Series {i+1}",
                filetypes=[("MAT files", "*_py.mat")]
            )
            if f:
                selected_files.append(f)
            else:
                print(f"Selection cancelled for Series {i+1}.")
                break
        
    if not selected_files:
        print("No files selected.")
        return

    num_series = len(selected_files)
    print(f"Selected {num_series} files for comparison:")
    for f in selected_files:
        print(f" - {f}")

    # --- 2. Pre-Plot Configuration (Popup) ---
    # Create hidden root window if not already created (CLI case might not have it)
    if 'root' not in locals() or root is None:
        root = tk.Tk()
        root.withdraw()

    # Default legends
    default_legends = [f"Series {i+1}" for i in range(num_series)]
    
    config = {
        "title": DEFAULT_TITLE,
        "xlabel": DEFAULT_X_LABEL,
        "ylabel": DEFAULT_Y_LABEL,
        "legends": default_legends
    }
    
    def get_config(num_files):
        dialog = tk.Toplevel(root)
        dialog.title("Plot Configuration")
        
        # Dynamic height based on number of files
        height = 250 + (num_files * 35)
        dialog.geometry(f"400x{height}")
        
        # Title & Axis Labels
        tk.Label(dialog, text="Figure Title:").pack(pady=2)
        e_title = tk.Entry(dialog, width=50)
        e_title.insert(0, DEFAULT_TITLE)
        e_title.pack(pady=2)
        
        tk.Label(dialog, text="X Label:").pack(pady=2)
        e_xlabel = tk.Entry(dialog, width=50)
        e_xlabel.insert(0, DEFAULT_X_LABEL)
        e_xlabel.pack(pady=2)
        
        tk.Label(dialog, text="Y Label:").pack(pady=2)
        e_ylabel = tk.Entry(dialog, width=50)
        e_ylabel.insert(0, DEFAULT_Y_LABEL)
        e_ylabel.pack(pady=2)
        
        # Legend Labels
        tk.Label(dialog, text="Legends:").pack(pady=5)
        
        legend_entries = []
        
        # Scrollable frame could be better for MANY files, but let's stick to simple for now
        # assuming realistic use case of < 10 files.
        
        frame_legends = tk.Frame(dialog)
        frame_legends.pack(pady=2, fill="both", expand=True)
        
        # Add scrollbar if needed or just simple pack? 
        # For simplicity, just pack. If n is huge, valid concern, but YAGNI for now.
        
        for i in range(num_files):
            f_row = tk.Frame(frame_legends)
            f_row.pack(fill="x", padx=10, pady=1)
            tk.Label(f_row, text=f"S{i+1}:", width=5).pack(side="left")
            e = tk.Entry(f_row)
            e.insert(0, default_legends[i])
            e.pack(side="left", fill="x", expand=True)
            legend_entries.append(e)
        
        def on_submit():
            config["title"] = e_title.get()
            config["xlabel"] = e_xlabel.get()
            config["ylabel"] = e_ylabel.get()
            config["legends"] = [e.get() for e in legend_entries]
            dialog.destroy()
            
        tk.Button(dialog, text="Continue", command=on_submit).pack(pady=10)
        
        root.wait_window(dialog)
        
    get_config(num_series)
    root.destroy()

    # --- 3. Load Data & Extract Signal ---
    series_data = [] # List of (time, signal) tuples

    for file_path in selected_files:
        if not os.path.exists(file_path):
            print(f"Error: File not found: {file_path}")
            series_data.append(None)
            continue
            
        try:
            mat = scipy.io.loadmat(file_path, squeeze_me=True, struct_as_record=False, simplify_cells=True)
        except Exception as e:
            print(f"Error loading .mat file {file_path}: {e}")
            series_data.append(None)
            continue

        if 'data' not in mat:
            print(f"Error: 'data' struct not found in {file_path}.")
            series_data.append(None)
            continue
            
        data = mat['data']
        
        # Handle struct/dict differences
        if not isinstance(data, dict):
            try:
                field_names = data.dtype.names
                data = {name: data[name] for name in field_names}
            except:
                pass

        # Extract tau_thr_sum
        def get_signal(name, source):
            if name in source:
                val = source[name]
                if isinstance(val, dict) and 'data' in val: return val['time'], val['data']
                if hasattr(val, 'data'): return val.time, val.data
            return None, None

        t_sig, sig_vals = get_signal('tau_thr_sum', data)
        
        if sig_vals is None:
            print(f"Warning: 'tau_thr_sum' not found in {os.path.basename(file_path)}. Skipping this trace.")
            series_data.append(None)
            continue
            
        # Handle 2D cases
        if sig_vals.ndim == 2:
             if sig_vals.shape[0] == 1: sig_vals = sig_vals.flatten()
             elif sig_vals.shape[1] == 1: sig_vals = sig_vals.flatten()
        
        series_data.append((t_sig, sig_vals))


    # --- 4. Plotting ---
    # Larger Figure for visibility
    fig, ax = plt.subplots(figsize=(10, 6), dpi=100)
    
    # Extended colors for many series
    colors = ['b', 'r', 'g', 'm', 'c', 'y', 'k', 'orange', 'purple', 'brown']
    linestyles = ['-', '--', '-.', ':']
    
    for i, data_tuple in enumerate(series_data):
        if data_tuple is None:
            continue
            
        t, y = data_tuple
        lbl = config["legends"][i]
        col = colors[i % len(colors)]
        ls = linestyles[i % len(linestyles)]
        
        ax.plot(t, y, color=col, linestyle=ls, linewidth=2, label=lbl)

    # Apply Styles
    ax.set_title(config["title"], fontsize=FONT_SIZE + 4)
    ax.set_xlabel(config["xlabel"], fontsize=FONT_SIZE)
    ax.set_ylabel(config["ylabel"], fontsize=FONT_SIZE)
    
    if SHOW_GRID: ax.grid(True)
    ax.legend(fontsize=FONT_SIZE)
    
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare tau_thr_sum from 3 .mat files.")
    parser.add_argument("files", nargs='*', help="Paths to 3 _py.mat files")
    args = parser.parse_args()
    
    plot_tau_comparison(args.files)
