import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox

def plot_heading(file_path=None):
    # --- Configuration Defaults ---
    DEFAULT_TITLE = "Heading Variation"
    DEFAULT_X_LABEL = "Time [s]"
    DEFAULT_Y_LABEL = "Heading [deg]"
    
    # Fonts
    FONT_SIZE = 12
    plt.rcParams.update({'font.size': FONT_SIZE})

    # Line Styles
    ETA_COLOR = 'b'
    ETA_LINESTYLE = '-'
    ETA_LINEWIDTH = 1.5
    ETA_LABEL_DEF = r'$\psi$ (Actual)'
    
    ETAD_COLOR = 'r'
    ETAD_LINESTYLE = '--'
    ETAD_LINEWIDTH = 1.5
    ETAD_LABEL_DEF = r'$\psi_d$ (Desired)'
    
    BETA_COLOR = 'g'
    BETA_LINESTYLE = ':'
    BETA_LINEWIDTH = 1.5
    BETA_LABEL_DEF = r'$\beta_c$ (Current)'
    
    SHOW_GRID = True
    
    # --- 1. Pre-Plot Configuration (Popup) ---
    # Create hidden root window
    root = tk.Tk()
    root.withdraw() 
    
    # Prompt for user configuration
    config = {
        "title": DEFAULT_TITLE,
        "xlabel": DEFAULT_X_LABEL,
        "ylabel": DEFAULT_Y_LABEL,
        "show_psid": True,
        "show_beta": False,
        "label_psi": ETA_LABEL_DEF,
        "label_psid": ETAD_LABEL_DEF,
        "label_beta": BETA_LABEL_DEF
    }
    
    def get_config():
        dialog = tk.Toplevel(root)
        dialog.title("Plot Configuration")
        dialog.geometry("400x500")
        
        # Frame: General Labels
        lf_gen = tk.LabelFrame(dialog, text="General Settings")
        lf_gen.pack(pady=5, padx=10, fill="x")
        
        tk.Label(lf_gen, text="Figure Title:").grid(row=0, column=0, sticky="e", padx=5, pady=2)
        e_title = tk.Entry(lf_gen, width=30)
        e_title.insert(0, DEFAULT_TITLE)
        e_title.grid(row=0, column=1, padx=5, pady=2)
        
        tk.Label(lf_gen, text="X Label:").grid(row=1, column=0, sticky="e", padx=5, pady=2)
        e_xlabel = tk.Entry(lf_gen, width=30)
        e_xlabel.insert(0, DEFAULT_X_LABEL)
        e_xlabel.grid(row=1, column=1, padx=5, pady=2)
        
        tk.Label(lf_gen, text="Y Label:").grid(row=2, column=0, sticky="e", padx=5, pady=2)
        e_ylabel = tk.Entry(lf_gen, width=30)
        e_ylabel.insert(0, DEFAULT_Y_LABEL)
        e_ylabel.grid(row=2, column=1, padx=5, pady=2)
        
        # Frame: Series Labels
        lf_labels = tk.LabelFrame(dialog, text="Legend Labels")
        lf_labels.pack(pady=5, padx=10, fill="x")
        
        tk.Label(lf_labels, text="Psi Label:").grid(row=0, column=0, sticky="e", padx=5, pady=2)
        e_l_psi = tk.Entry(lf_labels, width=30)
        e_l_psi.insert(0, ETA_LABEL_DEF)
        e_l_psi.grid(row=0, column=1, padx=5, pady=2)

        tk.Label(lf_labels, text="Psi_d Label:").grid(row=1, column=0, sticky="e", padx=5, pady=2)
        e_l_psid = tk.Entry(lf_labels, width=30)
        e_l_psid.insert(0, ETAD_LABEL_DEF)
        e_l_psid.grid(row=1, column=1, padx=5, pady=2)

        tk.Label(lf_labels, text="Beta Label:").grid(row=2, column=0, sticky="e", padx=5, pady=2)
        e_l_beta = tk.Entry(lf_labels, width=30)
        e_l_beta.insert(0, BETA_LABEL_DEF)
        e_l_beta.grid(row=2, column=1, padx=5, pady=2)

        # Frame: Visibility
        lf_vis = tk.LabelFrame(dialog, text="Visibility")
        lf_vis.pack(pady=5, padx=10, fill="x")
        
        var_psid = tk.BooleanVar(value=True)
        tk.Checkbutton(lf_vis, text="Show Desired Heading (Psi_d)", variable=var_psid).pack(anchor="w", padx=5)
        
        var_beta = tk.BooleanVar(value=False)
        tk.Checkbutton(lf_vis, text="Show Current Beta (Beta_c)", variable=var_beta).pack(anchor="w", padx=5)

        # Frame: Legend Settings
        lf_leg = tk.LabelFrame(dialog, text="Legend Settings")
        lf_leg.pack(pady=5, padx=10, fill="x")
        
        tk.Label(lf_leg, text="Location:").grid(row=0, column=0, sticky="e", padx=5, pady=2)
        
        legend_locs = [
            "best", "upper right", "upper left", "lower left", "lower right", 
            "right", "center left", "center right", "lower center", "upper center", "center"
        ]
        var_loc = tk.StringVar(value="best")
        opt_loc = tk.OptionMenu(lf_leg, var_loc, *legend_locs)
        opt_loc.config(width=15)
        opt_loc.grid(row=0, column=1, padx=5, pady=2, sticky="w")
        
        var_draggable = tk.BooleanVar(value=True)
        tk.Checkbutton(lf_leg, text="Draggable Legend", variable=var_draggable).grid(row=1, column=0, columnspan=2, sticky="w", padx=5)
        
        def on_submit():
            config["title"] = e_title.get()
            config["xlabel"] = e_xlabel.get()
            config["ylabel"] = e_ylabel.get()
            config["label_psi"] = e_l_psi.get()
            config["label_psid"] = e_l_psid.get()
            config["label_beta"] = e_l_beta.get()
            config["show_psid"] = var_psid.get()
            config["show_beta"] = var_beta.get()
            config["legend_loc"] = var_loc.get()
            config["legend_draggable"] = var_draggable.get()
            dialog.destroy()
            
        tk.Button(dialog, text="Continue", command=on_submit).pack(pady=15)
        
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

    t_eta, eta_vals = get_signal('eta', data)
    if eta_vals is None:
        print("Error: 'eta' signal not found.")
        return
        
    if eta_vals.ndim == 2:
        if eta_vals.shape[1] != 3 and eta_vals.shape[0] == 3:
            eta_vals = eta_vals.T
    
    psi = eta_vals[:, 2]

    t_etad, etad_vals = get_signal('etad', data)
    psi_d = None
    if etad_vals is not None:
        if etad_vals.ndim == 2:
            if etad_vals.shape[1] != 3 and etad_vals.shape[0] == 3:
                etad_vals = etad_vals.T
        psi_d = etad_vals[:, 2]
    else:
        print("Warning: 'etad' signal not found. Plotting only actual heading.")

    # Helper: Wrap to +/- PI
    def wrap_to_pi(angle):
        return (angle + np.pi) % (2 * np.pi) - np.pi

    psi_deg = np.rad2deg(wrap_to_pi(psi))
    if psi_d is not None:
        psi_d_deg = np.rad2deg(wrap_to_pi(psi_d))
        
    # Extract Beta
    beta_deg = None
    t_beta = None
    if config["show_beta"]:
        possible_keys = ['beta_current', 's_beta_current_', 'beta_c', 'beta']
        beta_vals = None
        for k in possible_keys:
            t_b, val = get_signal(k, data)
            if val is not None:
                beta_vals = val
                t_beta = t_b
                print(f"Found beta signal: {k}")
                break
        
        if beta_vals is not None:
             # Transform: beta_plot = beta_raw - pi ("Toward" to "From")
             beta_plot = beta_vals - np.pi
             beta_deg = np.rad2deg(wrap_to_pi(beta_plot))
        else:
            print("Warning: Beta signal requested but not found.")

    # --- 5. Plotting ---
    # Square Figure (400x400 pixels approx)
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
    
    # Plot Lines
    ax.plot(t_eta, psi_deg, color=ETA_COLOR, linestyle=ETA_LINESTYLE, 
            linewidth=ETA_LINEWIDTH, label=config["label_psi"])
    
    if psi_d is not None and config["show_psid"]:
        ax.plot(t_etad, psi_d_deg, color=ETAD_COLOR, linestyle=ETAD_LINESTYLE, 
                linewidth=ETAD_LINEWIDTH, label=config["label_psid"])

    if beta_deg is not None and config["show_beta"]:
         # Ensure time vector matches if different (usually same for simulink logsout, but good to be safe)
         # If t_beta and t_eta are same size/values just plot, otherwise handle mismatch if strictly needed.
         # For this specific context, we assume alignment or simply plot vs its own time.
         ax.plot(t_beta, beta_deg, color=BETA_COLOR, linestyle=BETA_LINESTYLE,
                 linewidth=BETA_LINEWIDTH, label=config["label_beta"])

    # Apply Styles from Config
    ax.set_title(config["title"], fontsize=FONT_SIZE + 2) # Slightly larger title
    ax.set_xlabel(config["xlabel"], fontsize=FONT_SIZE)
    ax.set_ylabel(config["ylabel"], fontsize=FONT_SIZE)
    if SHOW_GRID: ax.grid(True)
    leg = ax.legend(fontsize=FONT_SIZE, loc=config.get("legend_loc", "best"))
    if config.get("legend_draggable", True):
        leg.set_draggable(True)
    
    # Adjust ticks size
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Heading (Psi vs Psi_d) from converted MAT file.")
    parser.add_argument("file", nargs='?', help="Path to _py.mat file")
    args = parser.parse_args()
    
    plot_heading(args.file)
