import scipy.io
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import argparse
from tkinter import filedialog, Tk

def plot_trajectory(file_path=None):
    # --- 1. Handle Input ---
    if not file_path:
        root = Tk()
        root.withdraw() # Hide small window
        file_path = filedialog.askopenfilename(
            title="Select Converted Python Data File (_py.mat)",
            filetypes=[("MAT files", "*_py.mat")]
        )
        root.destroy()
        
    if not file_path:
        print("No file selected.")
        return

    print(f"Loading: {file_path}")
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return

    # --- 2. Load Data ---
    # struct_as_record=False ensures structs are loaded as objects/dictionaries 
    # simplify_cells=True helps flatten the hierarchy
    try:
        mat = scipy.io.loadmat(file_path, squeeze_me=True, struct_as_record=False, simplify_cells=True)
    except Exception as e:
        print(f"Error loading .mat file: {e}")
        return

    # The data struct is inside 'data' key
    if 'data' not in mat:
        print("Error: Could not find 'data' struct in file. Run convert_for_python.m first.")
        return
        
    data = mat['data']
    
    # Extract keys (signal names)
    # Depending on how it loaded, 'data' might be a dictionary or object
    # With simplify_cells=True, it should be a dictionary-like structure
    if not isinstance(data, dict):
        # Try to convert to dict if it's a structured array
        try:
            # If it's a void object (scipy behavior for structs sometimes)
            field_names = data.dtype.names
            data_dict = {name: data[name] for name in field_names}
            data = data_dict
        except:
             # If simplify_cells worked well, data is a dict
             pass

    # Extract Eta (Vessel State)
    if 'eta' in data:
        eta_struct = data['eta']
        t = eta_struct['time']
        eta_vals = eta_struct['data'] # (Time, 3) -> [N, E, Psi]
        
        # Determine index mapping (Scipy loading can sometimes transpose)
        # We expect (Time, 3) based on our conversion script
        if eta_vals.shape[1] != 3 and eta_vals.shape[0] == 3:
             eta_vals = eta_vals.T
             
        y_east = eta_vals[:, 1]
        x_north = eta_vals[:, 0]
        psi = eta_vals[:, 2]
    else:
        print("Error: 'eta' signal not found in data.")
        return

 

    # --- 3. Configuration ---
    ship_size_scale = 100
    plot_interval_s = 120 # Plot every 60 seconds

    # Style Configuration
    # Vessel Path
    eta_line_color = 'b'         
    eta_line_style = '-'         
    eta_line_width = 1.5         
    
    # Background Current Arrows
    current_arrow_color = 'lightblue'
    current_arrow_alpha = 0.5
    current_arrow_width = 0.005  
    current_arrow_headwidth = 4
    current_arrow_scale = 20 # Scale factor relative to plot width (smaller = longer)
    
    # Vessel Shape (Triangle)
    vessel_face_color = 'navajowhite'
    vessel_edge_color = 'darkorange'
    vessel_line_width = 2
    
    # Create Indices for plotting
    # Assuming t is valid
    t_start = np.floor(t[0])
    t_end = np.floor(t[-1])
    plot_times = np.arange(t_start, t_end + 1, plot_interval_s)
    
    plot_indices = []
    for pt in plot_times:
        idx = (np.abs(t - pt)).argmin()
        plot_indices.append(idx)
    
    # Ensure unique and include end
    plot_indices = sorted(list(set(plot_indices)))
    if len(t) - 1 not in plot_indices:
        plot_indices.append(len(t) - 1)


    # --- 4. Plotting ---
    # Create figure and axis with space at the bottom for widgets
    fig, ax = plt.subplots(figsize=(10, 9))
    # Create figure and axis with space at the bottom for widgets
    fig, ax = plt.subplots(figsize=(10, 9))
    plt.subplots_adjust(bottom=0.20) # Compact margin

    # --- 4a. Plot Current/Wind Vector Field (Background) ---
    # Try to find current data: 's_V_current_' and 's_beta_current_'
    # Names might vary slightly depending on cleaning, check typical candidates
    vc_keys = ['s_V_current_', 'V_current', 'vc']
    beta_keys = ['s_beta_current_', 'beta_current', 'beta_c']
    
    vc_data = None
    beta_data = None
    
    # Helper to find key in data
    def get_signal_data(keys, source):
        for k in keys:
            if k in source:
                # Handle struct/dict differences
                val = source[k]
                if hasattr(val, 'data'): return val.data
                if isinstance(val, dict) and 'data' in val: return val['data']
                # If it's a direct array/value
                return val
        return None

    vc_vals = get_signal_data(vc_keys, data)
    beta_vals = get_signal_data(beta_keys, data)

    if vc_vals is not None and beta_vals is not None:
        # Calculate Mean Current
        # Assuming values might be arrays or scalars
        if np.ndim(vc_vals) > 0: vc_mag = np.mean(vc_vals)
        else: vc_mag = vc_vals
        
        if np.ndim(beta_vals) > 0: beta_dir = np.mean(beta_vals)
        else: beta_dir = beta_vals
        
        print(f"Plotting Background Current: Vc={vc_mag:.2f} m/s, Beta={beta_dir:.2f} rad")
        
        print(f"Plotting Background Current: Vc={vc_mag:.2f} m/s, Beta={beta_dir:.2f} rad")
        
        # --- NEW: Screen-Relative Background Arrows ---
        # Create a grid in axes coordinates (0,1)
        # This ensures arrows cover the whole visible area regardless of zoom.
        n_arr_x = 15
        n_arr_y = 15
        
        grid_x = np.linspace(0.05, 0.95, n_arr_x)
        grid_y = np.linspace(0.05, 0.95, n_arr_y)
        X_grid, Y_grid = np.meshgrid(grid_x, grid_y)
        
        # Velocity Components (Direction)
        # We want the arrows to point in the physical direction (North/East).
        # Since plot X=East, Y=North, and screen X/Y align with these:
        # U_plot (East) corresponds to Screen X direction
        # V_plot (North) corresponds to Screen Y direction
        
        U_plot = vc_mag * np.sin(beta_dir) # East component
        V_plot = vc_mag * np.cos(beta_dir) # North component

        # Normalize vectors for constant visual length
        magnitude = np.sqrt(U_plot**2 + V_plot**2)
        if magnitude > 0:
            U_norm = U_plot / magnitude
            V_norm = V_plot / magnitude
        else:
            U_norm = 0
            V_norm = 0
        
        # Broadcast to grid
        U_field = np.full_like(X_grid, U_norm)
        V_field = np.full_like(Y_grid, V_norm)
        
        # Plot using transform=ax.transAxes
        # pivoting='mid' centers the arrow.
        # scale_units='width', scale=current_arrow_scale makes length relative to width.
        ax.quiver(X_grid, Y_grid, U_field, V_field, 
                  transform=ax.transAxes, 
                  color=current_arrow_color, alpha=current_arrow_alpha, 
                  pivot='mid',
                  units='width', width=current_arrow_width, headwidth=current_arrow_headwidth, 
                  scale=current_arrow_scale, scale_units='width',
                  zorder=0, label='Current')


    
    # Plot Trajectory Line
    # Plot Trajectory Line
    vessel_line, = ax.plot(y_east, x_north, color=eta_line_color, linestyle=eta_line_style, linewidth=eta_line_width, label='Vessel Path')
    
    # --- Interactive Vessel Plotting ---
    
    # Container for vessel patches to clear easily
    vessel_patches = []

    def update_vessels(val=None):
        # Clear existing patches
        for p in vessel_patches:
            p.remove()
        vessel_patches.clear()
        
        # Get values from sliders
        curr_interval = s_interval.val
        curr_scale = s_scale.val
        curr_alpha_max = s_alpha_max.val
        curr_alpha_min = s_alpha_min.val
        
        # Recalculate Indices
        t_s = np.floor(t[0])
        t_e = np.floor(t[-1])
        p_times = np.arange(t_s, t_e + 1, curr_interval)
        
        p_indices = []
        for pt in p_times:
            idx = (np.abs(t - pt)).argmin()
            p_indices.append(idx)
        p_indices = sorted(list(set(p_indices)))
        if len(t) - 1 not in p_indices:
            p_indices.append(len(t) - 1)
            
        # Recalculate Geometry
        # Use current axis span to determine size relative to 'screen range'
        # b_len represents a fraction of the visible width
        xlim = ax.get_xlim()
        view_span = xlim[1] - xlim[0]
        
        # If view_span is tiny (e.g. initial 0-1), scale might be off, but usually data sets limits.
        scale_factor = curr_scale / 1000.0
        b_len = view_span * scale_factor
        s_wid = b_len / 2
        
        b_verts = np.array([
            [b_len/2, 0],             # Bow
            [-b_len/2, s_wid/2],      # Rear Starboard
            [-b_len/2, -s_wid/2]      # Rear Port
        ])
        
        # Plot
        for i, idx in enumerate(p_indices):
            n_p = x_north[idx]
            e_p = y_east[idx]
            h_p = psi[idx]
            
            # Rotate
            R = np.array([
                [np.cos(h_p), -np.sin(h_p)],
                [np.sin(h_p), np.cos(h_p)]
            ])
            
            r_verts = (R @ b_verts.T).T 
            
            # Plot (X=East, Y=North)
            f_verts = np.column_stack((
                e_p + r_verts[:, 1], 
                n_p + r_verts[:, 0]
            ))
            
            # Alpha gradient (from min to max)
            if len(p_indices) > 1:
                # Linear interp: min + (max-min) * progress
                progress = i / (len(p_indices) - 1)
                a_val = curr_alpha_min + (curr_alpha_max - curr_alpha_min) * progress
            else:
                a_val = curr_alpha_max
                
            poly = patches.Polygon(f_verts, closed=True, 
                                   facecolor=vessel_face_color, edgecolor=vessel_edge_color, 
                                   linewidth=vessel_line_width, alpha=a_val, zorder=10)
            ax.add_patch(poly)
            vessel_patches.append(poly)
            
        fig.canvas.draw_idle()

    # Styling
    initial_title = f"Vessel Trajectory: {os.path.basename(file_path)}"
    ax.set_title(initial_title, fontsize=14)
    ax.set_xlabel("East [m]", fontsize=12)
    ax.set_ylabel("North [m]", fontsize=12)
    ax.grid(True)
    ax.axis('equal')
    ax.legend()
    
    # --- Widgets ---
    from matplotlib.widgets import Button, TextBox
    
    # Define locations [left, bottom, width, height]
    # Row 1: Title
    ax_text = plt.axes([0.15, 0.11, 0.75, 0.035])
    
    # Row 2: Interval, Scale, Line Width
    ax_interval = plt.axes([0.15, 0.06, 0.15, 0.035])
    ax_scale    = plt.axes([0.45, 0.06, 0.15, 0.035])
    ax_lw       = plt.axes([0.75, 0.06, 0.15, 0.035])
    
    # Row 3: Alpha Min, Alpha Max, Save
    ax_alpha_min = plt.axes([0.15, 0.01, 0.15, 0.035])
    ax_alpha_max = plt.axes([0.45, 0.01, 0.15, 0.035])
    ax_save      = plt.axes([0.70, 0.01, 0.20, 0.035])
    
    # Create Widgets
    text_box = TextBox(ax_text, 'Title: ', initial=initial_title)
    
    tb_interval = TextBox(ax_interval, 'Interval [s]: ', initial=str(plot_interval_s))
    tb_scale    = TextBox(ax_scale,    'Scale (%): ',    initial=str(ship_size_scale)) 
    tb_lw       = TextBox(ax_lw,       'Line Width: ',   initial=str(eta_line_width))
    tb_amin     = TextBox(ax_alpha_min,'Min Alpha: ',    initial="0.3")
    tb_amax     = TextBox(ax_alpha_max,'Max Alpha: ',    initial="1.0")
    
    def submit_update(text):
        # Refresh plot with new values
        update_vessels()
        fig.canvas.draw_idle()

    def update_vessels(text=None):
        # Clear existing patches
        for p in vessel_patches:
            p.remove()
        vessel_patches.clear()
        
        # Parse Values (safe parsing)
        try:
            val_interval = float(tb_interval.text)
            val_scale = float(tb_scale.text)
            val_lw = float(tb_lw.text)
            val_amin = float(tb_amin.text)
            val_amax = float(tb_amax.text)
        except ValueError:
            print("Invalid input in text boxes. Using defaults/last valid.")
            return

        # Update Line Width
        vessel_line.set_linewidth(val_lw)

        # Recalculate Indices
        t_s = np.floor(t[0])
        t_e = np.floor(t[-1])
        # Protect against zero interval
        if val_interval <= 0: val_interval = 10
        
        p_times = np.arange(t_s, t_e + 1, val_interval)
        
        p_indices = []
        for pt in p_times:
            idx = (np.abs(t - pt)).argmin()
            p_indices.append(idx)
        p_indices = sorted(list(set(p_indices)))
        if len(t) - 1 not in p_indices:
            p_indices.append(len(t) - 1)
            
        # Recalculate Geometry
        # Use current axis span to determine size relative to 'screen range'
        xlim = ax.get_xlim()
        view_span = xlim[1] - xlim[0]
        
        scale_factor = val_scale / 1000.0
        b_len = view_span * scale_factor
        s_wid = b_len / 2
        
        b_verts = np.array([
            [b_len/2, 0],             # Bow
            [-b_len/2, s_wid/2],      # Rear Starboard
            [-b_len/2, -s_wid/2]      # Rear Port
        ])
        
        # Plot
        for i, idx in enumerate(p_indices):
            n_p = x_north[idx]
            e_p = y_east[idx]
            h_p = psi[idx]
            
            # Rotate
            R = np.array([
                [np.cos(h_p), -np.sin(h_p)],
                [np.sin(h_p), np.cos(h_p)]
            ])
            
            r_verts = (R @ b_verts.T).T 
            
            # Plot (X=East, Y=North)
            f_verts = np.column_stack((
                e_p + r_verts[:, 1], 
                n_p + r_verts[:, 0]
            ))
            
            # Alpha gradient (from min to max)
            if len(p_indices) > 1:
                # Linear interp: min + (max-min) * progress
                progress = i / (len(p_indices) - 1)
                a_val = val_amin + (val_amax - val_amin) * progress
            else:
                a_val = val_amax
                
            poly = patches.Polygon(f_verts, closed=True, 
                                   facecolor=vessel_face_color, edgecolor=vessel_edge_color, 
                                   linewidth=vessel_line_width, alpha=a_val, zorder=10)
            ax.add_patch(poly)
            vessel_patches.append(poly)
            
        fig.canvas.draw_idle()
    
    def update_title(text):
        ax.set_title(text, fontsize=14)
        fig.canvas.draw_idle()

    # Link Callbacks
    text_box.on_submit(update_title)
    tb_interval.on_submit(submit_update)
    tb_scale.on_submit(submit_update)
    tb_lw.on_submit(submit_update)
    tb_amin.on_submit(submit_update)
    tb_amax.on_submit(submit_update)
    
    # Save Button
    b_save = Button(ax_save, 'Save Plot')
    
    def save_plot(event):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.dirname(script_dir)
        output_dir = os.path.join(base_dir, 'output')
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        script_name = os.path.splitext(os.path.basename(__file__))[0]
        input_name = os.path.splitext(os.path.basename(file_path))[0]
        if input_name.endswith('_py'):
            input_name = input_name[:-3]
            
        # Append params to filename?
        # out_filename = f"{script_name}_{input_name}_int{int(s_interval.val)}_s{int(s_scale.val)}.png"
        
        # Keep simple name or overwrite
        out_filename = f"{script_name}_{input_name}.png"
        out_path = os.path.join(output_dir, out_filename)
        
        plt.savefig(out_path, dpi=300)
        print(f"Plot saved to: {out_path}")
        
    b_save.on_clicked(save_plot)
    
    # Initial Draw
    update_vessels()
    
    # Show plot
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot trajectory from converted MAT file.")
    parser.add_argument("file", nargs='?', help="Path to _py.mat file")
    args = parser.parse_args()
    
    plot_trajectory(args.file)
