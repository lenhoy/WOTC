import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

def plot_selection(file_path=None):
    # --- Configuration Defaults ---
    DEFAULT_TITLE = "Selected Signals"
    DEFAULT_X_LABEL = "Time [s]"
    DEFAULT_Y_LABEL = "Value"
    
    # Fonts
    FONT_SIZE = 12
    plt.rcParams.update({'font.size': FONT_SIZE})
    
    SHOW_GRID = True
    
    # --- 1. Handle Input ---
    if not file_path:
        root = tk.Tk()
        root.withdraw() 
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

    # --- 2. Load Data ---
    try:
        # consistent loading with other scripts
        mat = scipy.io.loadmat(file_path, squeeze_me=True, struct_as_record=False, simplify_cells=True)
    except Exception as e:
        print(f"Error loading .mat file: {e}")
        return

    if 'data' not in mat:
        print("Error: Could not find 'data' struct.")
        return
        
    data = mat['data']
    
    # Ensure data is a dict
    if not isinstance(data, dict):
        try:
            field_names = data.dtype.names
            data = {name: data[name] for name in field_names}
        except:
            pass

    # --- 3. Signal Discovery ---
    # We want to find all items that have 'time' and 'data'
    
    found_signals = {} # "Display Name" -> (Time Array, Value Array)

    def process_node(node, path_prefix=""):
        if isinstance(node, dict):
            # Check if this node itself is a signal
            if 'time' in node and 'data' in node:
                t = node['time']
                v = node['data']
                
                # Check validity
                if isinstance(t, (np.ndarray, list)) and isinstance(v, (np.ndarray, list)):
                    v = np.array(v)
                    t = np.array(t)
                    
                    # Handle Dimensions
                    if v.ndim == 1:
                        found_signals[path_prefix] = (t, v)
                    elif v.ndim == 2:
                        # Handle (3, N) vs (N, 3)
                        # Assume time dimension matches t length
                        if v.shape[0] == len(t) and v.shape[1] < 20: # Likely (N, Cols)
                            for i in range(v.shape[1]):
                                found_signals[f"{path_prefix} (Col {i+1})"] = (t, v[:, i])
                        elif v.shape[1] == len(t) and v.shape[0] < 20: # Likely (Cols, N)
                            for i in range(v.shape[0]):
                                found_signals[f"{path_prefix} (Col {i+1})"] = (t, v[i, :])
                        else:
                             found_signals[f"{path_prefix} (Multi-dim)"] = (t, v) # Fallback
            
            # Recurse into children
            for key, child in node.items():
                new_prefix = f"{path_prefix}.{key}" if path_prefix else key
                # simple cycle prevention or depth check could be added if needed
                process_node(child, new_prefix)

    process_node(data)
    
    if not found_signals:
        print("No plot-able signals found in data.")
        return

    sorted_names = sorted(found_signals.keys())

    # --- 4. Selection UI ---
    selection_root = tk.Tk()
    selection_root.title("Select Signals to Plot")
    selection_root.geometry("400x600")
    
    # DataVars
    check_vars = {}
    
    # Layout
    lbl = tk.Label(selection_root, text="Available Signals:", font=("Arial", 12, "bold"))
    lbl.pack(pady=5)
    
    # Scrollable Frame
    container = tk.Frame(selection_root)
    container.pack(fill="both", expand=True, padx=10, pady=5)
    
    canvas = tk.Canvas(container)
    scrollbar = tk.Scrollbar(container, orient="vertical", command=canvas.yview)
    scrollable_frame = tk.Frame(canvas)

    scrollable_frame.bind(
        "<Configure>",
        lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
    )

    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)

    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")
    
    # Populate Checkboxes
    for name in sorted_names:
        var = tk.BooleanVar()
        chk = tk.Checkbutton(scrollable_frame, text=name, variable=var, anchor='w')
        chk.pack(fill='x', expand=True)
        check_vars[name] = var

    selected_signals = []

    def on_plot():
        for name, var in check_vars.items():
            if var.get():
                selected_signals.append(name)
        selection_root.destroy()
        
    btn_plot = tk.Button(selection_root, text="Plot Selected", command=on_plot, bg="#dddddd", height=2)
    btn_plot.pack(fill="x", padx=10, pady=10)
    
    selection_root.mainloop()

    if not selected_signals:
        print("No signals selected.")
        return

    # --- 5. Interactive Plotting with Controls ---
    from matplotlib.widgets import TextBox, Button, CheckButtons
    
    # Main Plot Figure
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100, layout='constrained')
    fig.canvas.manager.set_window_title("Selected Signals Plot")

    # Values for interactive update
    current_config = {
        'title': DEFAULT_TITLE,
        'xlabel': DEFAULT_X_LABEL,
        'ylabel': DEFAULT_Y_LABEL,
        'grid': SHOW_GRID
    }
    
    # Store lines and visibility
    plot_lines = {}
    line_visibility = {name: True for name in selected_signals}
    display_names = {name: name for name in selected_signals}

    # Initial Plot
    for name in selected_signals:
        t, v = found_signals[name]
        line, = ax.plot(t, v, linewidth=1.5, label=display_names[name])
        plot_lines[name] = line
        
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE)
    
    def update_plot(val=None):
        # Update Titles/Labels
        ax.set_title(current_config['title'], fontsize=FONT_SIZE + 2)
        ax.set_xlabel(current_config['xlabel'], fontsize=FONT_SIZE)
        ax.set_ylabel(current_config['ylabel'], fontsize=FONT_SIZE)
        ax.grid(current_config['grid'])
        
        # Update Visibility
        for name, line in plot_lines.items():
            line.set_visible(line_visibility[name])
            
        # Update Legend
        handles = []
        labels = []
        for name, line in plot_lines.items():
            if line.get_visible():
                handles.append(line)
                labels.append(display_names[name])
                
        if handles:
            ax.legend(handles, labels, fontsize=FONT_SIZE, loc='best')
        else:
            if ax.get_legend(): ax.get_legend().remove()
            
        fig.canvas.draw_idle()

    # --- Controls Window ---
    fig_ctrl = plt.figure(figsize=(4, 6))
    fig_ctrl.canvas.manager.set_window_title("Controls")
    
    # Header
    fig_ctrl.text(0.5, 0.95, "Plot Settings", ha='center', fontsize=12, weight='bold')
    
    # 1. Title/Labels
    ax_tit = fig_ctrl.add_axes([0.3, 0.85, 0.6, 0.05])
    tb_title = TextBox(ax_tit, "Title: ", initial=current_config['title'])
    
    ax_xl = fig_ctrl.add_axes([0.3, 0.78, 0.6, 0.05])
    tb_xlabel = TextBox(ax_xl, "X Label: ", initial=current_config['xlabel'])
    
    ax_yl = fig_ctrl.add_axes([0.3, 0.71, 0.6, 0.05])
    tb_ylabel = TextBox(ax_yl, "Y Label: ", initial=current_config['ylabel'])
    
    def on_text_submit(text):
        current_config['title'] = tb_title.text
        current_config['xlabel'] = tb_xlabel.text
        current_config['ylabel'] = tb_ylabel.text
        update_plot()
        
    tb_title.on_submit(on_text_submit)
    tb_xlabel.on_submit(on_text_submit)
    tb_ylabel.on_submit(on_text_submit)
    
    # 2. Visibility Toggles
    fig_ctrl.text(0.5, 0.60, "Visible Signals", ha='center', fontsize=10, weight='bold')
    
    # Calculate height needed for checkboxes
    num_sigs = len(selected_signals)
    chk_height = min(0.35, num_sigs * 0.05)
    ax_chk = fig_ctrl.add_axes([0.1, 0.60 - chk_height - 0.05, 0.8, chk_height])
    
    chk = CheckButtons(ax_chk, selected_signals, actives=[True]*num_sigs)
    
    def on_check(label):
        if label in line_visibility:
            line_visibility[label] = not line_visibility[label]
            update_plot()
            
    chk.on_clicked(on_check)
    
    # 3. Rename Button
    ax_ren = fig_ctrl.add_axes([0.1, 0.15, 0.8, 0.05])
    b_ren = Button(ax_ren, 'Rename Series')
    
    def open_rename_dialog(event):
        rename_win = tk.Tk()
        rename_win.title("Rename Series")
        rename_win.geometry("400x400")
        
        entries = {}
        
        tk.Label(rename_win, text="Edit Legend Labels:", font=("Arial", 10, "bold")).pack(pady=10)
        
        # Scrollable container
        container = tk.Frame(rename_win)
        container.pack(fill="both", expand=True, padx=10)
        canvas = tk.Canvas(container)
        scrollbar = tk.Scrollbar(container, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Define a function to limit the size of the window based on content
        # ... logic skipped for brevity ...

        for name in selected_signals:
            row = tk.Frame(scrollable_frame)
            row.pack(fill='x', pady=2)
            tk.Label(row, text=f"{name}: ", width=20, anchor='e').pack(side='left')
            e = tk.Entry(row, width=25)
            e.insert(0, display_names[name])
            e.pack(side='left', padx=5)
            entries[name] = e
            
        def apply_renames():
            for name, entry in entries.items():
                display_names[name] = entry.get()
            rename_win.destroy()
            update_plot()
            
        tk.Button(rename_win, text="Apply", command=apply_renames, bg="#dddddd").pack(pady=10)
        rename_win.mainloop()
        
    b_ren.on_clicked(open_rename_dialog)
    
    # 4. Save Button
    ax_save = fig_ctrl.add_axes([0.1, 0.05, 0.8, 0.08])
    b_save = Button(ax_save, 'Save Plot')
    
    def save_plot(event):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.dirname(script_dir)
        output_dir = os.path.join(base_dir, 'output')
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        script_name = os.path.splitext(os.path.basename(__file__))[0]
        input_name = os.path.splitext(os.path.basename(file_path))[0]
        if input_name.endswith('_py'): input_name = input_name[:-3]
        
        # Use title as part of filename if clean
        clean_title = "".join(x for x in current_config['title'] if x.isalnum() or x in " -_").strip()
        out_filename = f"{script_name}_{input_name}_{clean_title}.png"
        out_path = os.path.join(output_dir, out_filename)
        
        fig.savefig(out_path, dpi=300)
        print(f"Plot saved to: {out_path}")
        
    b_save.on_clicked(save_plot)

    # Initial Draw
    update_plot()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interactive generic signal plotter.")
    parser.add_argument("file", nargs='?', help="Path to _py.mat file")
    args = parser.parse_args()
    
    plot_selection(args.file)
