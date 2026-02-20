import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox

def plot_speed(file_path=None):
    # --- Configuration Defaults ---
    DEFAULT_TITLE = "Velocity"
    DEFAULT_X_LABEL = "Time [s]"
    DEFAULT_Y_LABEL = "Velocity [m/s]"
    
    # Fonts
    FONT_SIZE = 15
    plt.rcParams.update({'font.size': FONT_SIZE})

    # Line Styles
    U_COLOR = 'g'
    U_LINESTYLE = '-'
    U_LINEWIDTH = 1.5
    U_LABEL = r'$U$ (Actual)'
    
    UD_COLOR = 'm'
    UD_LINESTYLE = '--'
    UD_LINEWIDTH = 1.5
    UD_LABEL = r'$U_d$ (Desired)'
    
    SHOW_GRID = True
    
    # --- 1. Pre-Plot Configuration (Popup) ---
    root = tk.Tk()
    root.withdraw() 
    
    config = {
        "title": DEFAULT_TITLE,
        "xlabel": DEFAULT_X_LABEL,
        "ylabel": DEFAULT_Y_LABEL
    }
    
    def get_config():
        dialog = tk.Toplevel(root)
        dialog.title("Plot Configuration")
        dialog.geometry("300x250")
        
        tk.Label(dialog, text="Figure Title:").pack(pady=5)
        e_title = tk.Entry(dialog, width=40)
        e_title.insert(0, DEFAULT_TITLE)
        e_title.pack(pady=5)
        
        tk.Label(dialog, text="X Label:").pack(pady=5)
        e_xlabel = tk.Entry(dialog, width=40)
        e_xlabel.insert(0, DEFAULT_X_LABEL)
        e_xlabel.pack(pady=5)
        
        tk.Label(dialog, text="Y Label:").pack(pady=5)
        e_ylabel = tk.Entry(dialog, width=40)
        e_ylabel.insert(0, DEFAULT_Y_LABEL)
        e_ylabel.pack(pady=5)
        
        def on_submit():
            config["title"] = e_title.get()
            config["xlabel"] = e_xlabel.get()
            config["ylabel"] = e_ylabel.get()
            dialog.destroy()
            
        tk.Button(dialog, text="Continue", command=on_submit).pack(pady=20)
        
        root.wait_window(dialog)
        
    get_config()

    # --- 2. Handle Input ---
    if not file_path:
        file_path = filedialog.askopenfilename(
            title="Select Converted Python Data File (_py.mat)",
            filetypes=[("MAT files", "*_py.mat")]
        )
        
    root.destroy()
        
    if not file_path:
        print("No file selected.")
        return

    print(f"Loading: {file_path}")
    
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return

    # --- 3. Load Data ---
    try:
        mat = scipy.io.loadmat(file_path, squeeze_me=True, struct_as_record=False, simplify_cells=True)
    except Exception as e:
        print(f"Error loading .mat file: {e}")
        return

    if 'data' not in mat:
        print("Error: Could not find 'data' struct. Run convert_for_python.m first.")
        return
        
    data = mat['data']
    
    if not isinstance(data, dict):
        try:
            field_names = data.dtype.names
            data = {name: data[name] for name in field_names}
        except:
            pass

    # --- 4. Extract Data ---
    def get_signal(name, source):
        if name in source:
            val = source[name]
            if isinstance(val, dict) and 'data' in val: return val['time'], val['data']
            if hasattr(val, 'data'): return val.time, val.data
        return None, None

    # Get Speed (U)
    # Check for 'nu' (index 0) or 'U'
    time_u = None
    u_vals = None
    
    # Priority 1: 'U' (Total Speed)
    t_temp, u_temp = get_signal('U', data)
    if u_temp is not None:
        u_vals = u_temp
        time_u = t_temp

    # Priority 2: 'nu' (Surge Velocity) - only if U is missing
    if u_vals is None and 'nu' in data:
        t_nu, nu_data = get_signal('nu', data)
        if nu_data is not None:
            if nu_data.ndim == 2:
                if nu_data.shape[1] != 3 and nu_data.shape[0] == 3:
                     nu_data = nu_data.T
            u_vals = nu_data[:, 0]
            time_u = t_nu

    # Priority 3: 'u' (lowercase)
    if u_vals is None:
        t_temp, u_temp = get_signal('u', data)
        if u_temp is not None:
             u_vals = u_temp
             time_u = t_temp

    if u_vals is None:
        print("Error: Speed signal 'U' or 'nu' not found.")
        return

    # Get Desired Speed (Ud)
    time_ud, ud_vals = get_signal('Ud', data)
    if ud_vals is None: time_ud, ud_vals = get_signal('U_d', data)
    
    if ud_vals is None:
         print("Warning: Desired Speed 'Ud' not found. Plotting only actual speed.")

    # --- 5. Plotting ---
    # Square Figure (400x400 pixels approx)
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)

    # Plot Lines
    ax.plot(time_u, u_vals, color=U_COLOR, linestyle=U_LINESTYLE, 
            linewidth=U_LINEWIDTH, label=U_LABEL)
    
    if ud_vals is not None:
        ax.plot(time_ud, ud_vals, color=UD_COLOR, linestyle=UD_LINESTYLE, 
                linewidth=UD_LINEWIDTH, label=UD_LABEL)

    # Apply Styles from Config
    ax.set_title(config["title"], fontsize=FONT_SIZE + 2)
    ax.set_xlabel(config["xlabel"], fontsize=FONT_SIZE)
    ax.set_ylabel(config["ylabel"], fontsize=FONT_SIZE)
    if SHOW_GRID: ax.grid(True)
    ax.legend(fontsize=FONT_SIZE)
    
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Speed (U vs Ud) from converted MAT file.")
    parser.add_argument("file", nargs='?', help="Path to _py.mat file")
    args = parser.parse_args()
    
    plot_speed(args.file)
