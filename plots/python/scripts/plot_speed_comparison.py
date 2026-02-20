import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox

def plot_speed_comparison(file_paths=None):
    # --- Configuration Defaults ---
    DEFAULT_TITLE = "Velocity Comparison"
    DEFAULT_X_LABEL = "Time [s]"
    DEFAULT_Y_LABEL = "Velocity [m/s]"
    
    # Fonts and Styles (Matching plot_speed.py)
    FONT_SIZE = 15
    plt.rcParams.update({'font.size': FONT_SIZE})

    # Line Styles
    # Actual speed will cycle colors
    # Desired speed will be single series
    UD_COLOR = 'm'
    UD_LINESTYLE = '--'
    UD_LINEWIDTH = 1.5
    UD_LABEL = r'$U_d$ (Desired)'
    
    U_LINEWIDTH = 1.5

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
        num_series_input = simpledialog.askinteger("Input", "How many runs do you want to compare?", minvalue=1, maxvalue=20)
        
        if not num_series_input:
            print("No number entered or cancelled.")
            try:
                root.destroy()
            except:
                pass
            return

        print(f"Selecting {num_series_input} files...")
        
        for i in range(num_series_input):
            f = filedialog.askopenfilename(
                title=f"Select Data File for Run {i+1}",
                filetypes=[("MAT files", "*_py.mat")]
            )
            if f:
                selected_files.append(f)
            else:
                print(f"Selection cancelled for Run {i+1}.")
                break
        
    if not selected_files:
        print("No files selected.")
        try:
            root.destroy()
        except:
            pass
        return

    num_series = len(selected_files)
    print(f"Selected {num_series} files for comparison:")
    for f in selected_files:
        print(f" - {f}")

    # --- 2. Pre-Plot Configuration (Popup) ---
    # Create hidden root window if not already created
    if 'root' not in locals() or root is None:
        try:
            root = tk.Tk()
            root.withdraw()
        except:
            pass # might already exist

    # Default legends
    default_legends = [f"Run {i+1}" for i in range(num_series)]
    
    config = {
        "title": DEFAULT_TITLE,
        "xlabel": DEFAULT_X_LABEL,
        "ylabel": DEFAULT_Y_LABEL,
        "legends": default_legends,
        "visible": [True] * num_series
    }
    
    def get_config(num_files):
        try:
            dialog = tk.Toplevel(root)
        except:
            dialog = tk.Toplevel()
            
        dialog.title("Plot Configuration")
        
        # Dynamic height based on number of files
        height = 300 + (num_files * 35)
        dialog.geometry(f"500x{height}")
        
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
        
        # Series Configuration
        tk.Label(dialog, text="Series Configuration:", font=("bold")).pack(pady=10)
        
        series_widgets = []
        
        frame_legends = tk.Frame(dialog)
        frame_legends.pack(pady=2, fill="both", expand=True)
        
        for i in range(num_files):
            f_row = tk.Frame(frame_legends)
            f_row.pack(fill="x", padx=10, pady=2)
            
            # Show/Hide Checkbox
            var_show = tk.BooleanVar(value=True)
            chk = tk.Checkbutton(f_row, variable=var_show, text=f"Show")
            chk.pack(side="left")
            
            # Label
            tk.Label(f_row, text=f"Run {i+1} Name:", width=10).pack(side="left")
            
            # Entry
            e = tk.Entry(f_row)
            # Default name could be filename base if desired, but sticking to "Run i" as per default
            fname = os.path.basename(selected_files[i])
            e.insert(0, f"Run {i+1}: {fname[:15]}..." if len(fname)>15 else fname)
            e.pack(side="left", fill="x", expand=True)
            
            series_widgets.append({"entry": e, "var": var_show})
        
        def on_submit():
            config["title"] = e_title.get()
            config["xlabel"] = e_xlabel.get()
            config["ylabel"] = e_ylabel.get()
            config["legends"] = [w["entry"].get() for w in series_widgets]
            config["visible"] = [w["var"].get() for w in series_widgets]
            dialog.destroy()
            
        tk.Button(dialog, text="Plot", command=on_submit, height=2, width=15).pack(pady=20)
        
        # Wait for dialog
        dialog.wait_window()
        
    get_config(num_series)
    try:
        root.destroy()
    except:
        pass

    # --- 3. Load Data & Extract Signals ---
    
    series_data = [] # List of dicts: {'u': (t, y), 'label': str, 'visible': bool}
    
    # We only need one Ud (Desired Speed) series, assume same for all or take from first valid
    common_ud = None # (t, u_d)

    for i, file_path in enumerate(selected_files):
        if not config["visible"][i]:
            series_data.append(None)
            continue
            
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
                
        # Helper to extract signal
        def get_signal(name, source):
            if name in source:
                val = source[name]
                if isinstance(val, dict) and 'data' in val: return val['time'], val['data']
                if hasattr(val, 'data'): return val.time, val.data
            return None, None

        # Extract U (Actual Speed)
        time_u = None
        u_vals = None
        
        # Priority 1: 'U'
        t_temp, u_temp = get_signal('U', data)
        if u_temp is not None:
            u_vals = u_temp
            time_u = t_temp
            
        # Priority 2: 'nu' (Surge Velocity)
        if u_vals is None and 'nu' in data:
            t_nu, nu_data = get_signal('nu', data)
            if nu_data is not None:
                if nu_data.ndim == 2:
                    if nu_data.shape[1] != 3 and nu_data.shape[0] == 3:
                         nu_data = nu_data.T
                u_vals = nu_data[:, 0]
                time_u = t_nu
                
        # Priority 3: 'u'
        if u_vals is None:
            t_temp, u_temp = get_signal('u', data)
            if u_temp is not None:
                 u_vals = u_temp
                 time_u = t_temp
                 
        if u_vals is None:
            print(f"Warning: Speed signal not found in {os.path.basename(file_path)}")
            series_data.append(None)
            continue

        # Extract Ud (Desired Speed) - only if we haven't found one yet
        if common_ud is None:
            time_ud, ud_vals = get_signal('Ud', data)
            if ud_vals is None: time_ud, ud_vals = get_signal('U_d', data)
            
            if ud_vals is not None:
                common_ud = (time_ud, ud_vals)
        
        series_data.append({
            'u': (time_u, u_vals),
            'label': config["legends"][i],
            'visible': True
        })


    # --- 4. Plotting ---
    if not any(s is not None for s in series_data):
        print("No data to plot.")
        return

    # Square Figure matching plot_speed.py "sizing"
    # plot_speed.py uses figsize=(4, 4)
    # Since this is a comparison, it might get crowded, but user asked for "same sizing".
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
    
    # Extended colors for many series (excluding magenta 'm' used for Ud if possible to avoid confusion, 
    # but 'm' is distinct enough if line style differs)
    colors = ['b', 'r', 'g', 'c', 'y', 'k', 'orange', 'purple', 'brown']
    
    # Plot Actual Speeds
    plotted_count = 0
    for i, item in enumerate(series_data):
        if item is None:
            continue
            
        col = colors[i % len(colors)]
        label_base = item['label']
        
        t_u, u_val = item['u']
        ax.plot(t_u, u_val, color=col, linestyle='-', linewidth=U_LINEWIDTH, label=label_base)
        
        plotted_count += 1

    if plotted_count == 0:
        print("No series visible.")
        return

    # Plot Desired Speed (Once)
    if common_ud is not None:
        t_ud, ud_vals = common_ud
        ax.plot(t_ud, ud_vals, color=UD_COLOR, linestyle=UD_LINESTYLE, 
                linewidth=UD_LINEWIDTH, label=UD_LABEL)
    else:
        print("Warning: Desired Speed not found in any of the selected files.")

    # Apply Styles
    ax.set_title(config["title"], fontsize=FONT_SIZE + 2)
    ax.set_xlabel(config["xlabel"], fontsize=FONT_SIZE)
    ax.set_ylabel(config["ylabel"], fontsize=FONT_SIZE)
    
    if SHOW_GRID: ax.grid(True)
    
    ax.legend(fontsize=FONT_SIZE)
    
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Speed (U) from multiple .mat files.")
    parser.add_argument("files", nargs='*', help="Paths to _py.mat files")
    args = parser.parse_args()
    
    plot_speed_comparison(args.files)
