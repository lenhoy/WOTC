import scipy.io
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.widgets import Button, TextBox, CheckButtons
import numpy as np
import os
import argparse
from tkinter import filedialog, Tk

def plot_trajectory(file_path=None):
    # Style Configuration
    # Vessel Path
    eta_line_color = 'b'         
    eta_line_width = 1.5         
    
    # Vessel Shape (Triangle)
    vessel_face_color = 'navajowhite'
    vessel_edge_color = 'darkorange'
    vessel_line_width = 2
    
    # Reference (p0)
    p0_line_width = 2
    p0_color = 'grey' # Base color
    p0_marker_style = 'x'
    p0_marker_color = 'black'
    
    # Etad
    etad_color = 'green'
    etad_style = '--'
    etad_width = 1.5
    # etad_marker_style = 'o' # No longer used
    
    # Etad Shape (Triangle)
    etad_face_color = 'lightgreen' 
    etad_edge_color = 'green'
    etad_line_width = 1.5
    
    # Background Current
    current_arrow_color = 'lightblue'
    current_arrow_alpha = 0.5
    current_arrow_width = 0.06  # Inches
    current_arrow_headwidth = 6
    current_arrow_scale = 1.5   # Data units per inch (Lower = Larger arrows)
    
    # Fonts
    FONT_SIZE = 12
    plt.rcParams.update({'font.size': FONT_SIZE})

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
    if not isinstance(data, dict):
        # Try to convert to dict if it's a structured array
        try:
            field_names = data.dtype.names
            data_dict = {name: data[name] for name in field_names}
            data = data_dict
        except:
             pass

    # Helper to find key in data
    def get_signal_struct(keys, source):
        for k in keys:
            if k in source:
                return source[k]
        return None

    # Helper for extracting array data
    def get_signal_data(keys, source):
        for k in keys:
            if k in source:
                val = source[k]
                if hasattr(val, 'data'): return val.data
                if isinstance(val, dict) and 'data' in val: return val['data']
                return val
        return None

    # Extract Eta (Vessel State)
    if 'eta' in data:
        eta_struct = data['eta']
        t = eta_struct['time']
        eta_vals = eta_struct['data'] # (Time, 3) -> [N, E, Psi]
        
        if eta_vals.shape[1] != 3 and eta_vals.shape[0] == 3:
             eta_vals = eta_vals.T
             
        y_east = eta_vals[:, 1]
        x_north = eta_vals[:, 0]
        psi = eta_vals[:, 2]
    else:
        print("Error: 'eta' signal not found in data.")
        return

    # Extract p0 (Reference Trajectory) if available
    p0_available = False
    p0_north = None
    p0_east = None
    # We assume p0 is time-synced with eta if it comes from Simulink logsout
    
    p0_struct = get_signal_struct(['p0', 'ref_pos'], data)
    
    if p0_struct:
        try:
            p0_vals = p0_struct['data']
            if p0_vals.ndim == 2:
                if p0_vals.shape[1] < 2 and p0_vals.shape[0] >= 2:
                     p0_vals = p0_vals.T
                
                p0_north = p0_vals[:, 0]
                p0_east = p0_vals[:, 1]
                p0_available = True
                print("Found 'p0' signal.")
        except Exception as e:
            print(f"Error processing 'p0': {e}")

    # Extract etad (Desired Position) if available
    etad_available = False
    etad_north = None
    etad_east = None
    etad_psi = None
    
    etad_struct = get_signal_struct(['etad', 'pd', 'pos_d'], data)
    
    if etad_struct:
        try:
            etad_vals = etad_struct['data']
            if etad_vals.ndim == 2:
                # Correct orientation if needed
                if etad_vals.shape[1] < 3 and etad_vals.shape[0] >= 3:
                     etad_vals = etad_vals.T
                
                # Assuming [N, E, Psi]
                etad_north = etad_vals[:, 0]
                etad_east = etad_vals[:, 1]
                
                if etad_vals.shape[1] >= 3:
                    etad_psi = etad_vals[:, 2]
                else:
                    etad_psi = np.zeros_like(etad_north) # Default to 0 if no heading
                    
                etad_available = True
                print("Found 'etad' signal.")
        except Exception as e:
            print(f"Error processing 'etad': {e}")
    else:
        pass # Signal 'etad' not found, optional.


    # --- 3. Configuration (Runtime Params) ---
    ship_size_scale = 150
    plot_interval_s = 60 # Plot every X seconds

    # --- 4. Prepare Data ---
    # We will compute plot indices dynamically in update, but do initial setup

    # --- 5. Setup Figures ---
    # Main Plot Figure
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100, layout='constrained')
    fig.canvas.manager.set_window_title("Trajectory Plot")
    
    # Controls Figure
    fig_ctrl = plt.figure(figsize=(5, 7.5))
    fig_ctrl.canvas.manager.set_window_title("Controls")
    
    # Dictionary to store visibility states
    visibility = {
        'etad': True,
        'current': True,
        'p0': True,
        'pose': True,
        'p0_sync': True,
        'etad_sync': True
    }
    
    # Objects storage
    plot_objs = {
        'current_quiver': None,
        'p0_line': None,
        'p0_markers': None,
        'etad_line': None,
        # 'etad_markers': None, # Replaced by etad_patches
        'etad_patches': [], 
        'eta_line': None,
        'vessel_patches': []
    }

    # --- 5a. Plot Background Current ---
    vc_vals = get_signal_data(['s_V_current_', 'V_current', 'vc'], data)
    beta_vals = get_signal_data(['s_beta_current_', 'beta_current', 'beta_c'], data)

    if vc_vals is not None and beta_vals is not None:
        if np.ndim(vc_vals) > 0: vc_mag = np.mean(vc_vals)
        else: vc_mag = vc_vals
        if np.ndim(beta_vals) > 0: beta_dir = np.mean(beta_vals)
        else: beta_dir = beta_vals
        
        # Grid for arrows
        n_arr_x, n_arr_y = 4, 4
        grid_x = np.linspace(0.05, 0.95, n_arr_x)
        grid_y = np.linspace(0.05, 0.95, n_arr_y)
        X_grid, Y_grid = np.meshgrid(grid_x, grid_y)
        
        U_plot = vc_mag * np.sin(beta_dir) 
        V_plot = vc_mag * np.cos(beta_dir)
        mag = np.sqrt(U_plot**2 + V_plot**2)
        if mag > 0:
            U_norm, V_norm = U_plot/mag, V_plot/mag
        else:
            U_norm, V_norm = 0, 0
            
        U_field = np.full_like(X_grid, U_norm)
        V_field = np.full_like(Y_grid, V_norm)
        
        # Initial Plotting using Inches
        # Note: scale=1.5 means 1 data unit = 1/1.5 inches = 0.66 inches long
        plot_objs['current_quiver'] = ax.quiver(
            X_grid, Y_grid, U_field, V_field, 
            transform=ax.transAxes, 
            color=current_arrow_color, alpha=current_arrow_alpha, 
            pivot='mid', units='inches', width=current_arrow_width, 
            headwidth=current_arrow_headwidth, scale=current_arrow_scale, scale_units='inches',
            zorder=1, label='Current'
        )

    # --- 5b. Plot Static Lines (P0, Etad, Eta Path) ---
    
    # P0 Line
    if p0_available:
        # We plot the full line once, we don't need gradient segments for the line itself if using markers
        p0_line, = ax.plot(p0_east, p0_north, color=p0_color, linewidth=p0_line_width, alpha=0.5, label='p0 Reference', zorder=2)
        plot_objs['p0_line'] = p0_line
        
        # Initialize Markers
        p0_markers, = ax.plot([], [], linestyle='None', marker=p0_marker_style, color=p0_marker_color, label='p0 Sync', zorder=3)
        plot_objs['p0_markers'] = p0_markers

    # Etad Line
    if etad_available:
        etad_line, = ax.plot(etad_east, etad_north, color=etad_color, linestyle=etad_style, linewidth=etad_width, label='etad Desired', zorder=5)
        plot_objs['etad_line'] = etad_line
        
        # Etad Markers - REMOVED, replaced by patches in update_plot
        # etad_markers, = ax.plot([], [], linestyle='None', marker=etad_marker_style, color=etad_color, alpha=0.5, markersize=4, label='etad Sync', zorder=6)
        # plot_objs['etad_markers'] = etad_markers

    # Eta Path
    eta_line, = ax.plot(y_east, x_north, color=eta_line_color, linestyle='-', linewidth=eta_line_width, label='Vessel Path', zorder=4)
    plot_objs['eta_line'] = eta_line

    # --- Interactive Update Logic ---
    
    display_coords = True # Flag to avoid multiple regenerations if not needed, but here we redraw cheaply

    def update_plot(val=None):
        # 1. Update Visibility
        if plot_objs['current_quiver']:
            plot_objs['current_quiver'].set_visible(visibility['current'])
        if plot_objs['p0_line']:
            plot_objs['p0_line'].set_visible(visibility['p0'])
        if plot_objs['p0_markers']:
            plot_objs['p0_markers'].set_visible(visibility['p0_sync'])
        if plot_objs['etad_line']:
            plot_objs['etad_line'].set_visible(visibility['etad'])
        
        # 2. Update Parameters
        try:
            val_interval = float(tb_interval.text)
            val_scale = float(tb_scale.text)
            val_amin = float(tb_amin.text)
            val_amax = float(tb_amax.text)
            val_grid_n = int(tb_grid_n.text)
            val_arr_scale = float(tb_arr_scale.text)
            val_arr_width = float(tb_arr_width.text)
            val_lw = float(tb_lw.text)
        except ValueError:
            return

        # Update Line Widths
        if plot_objs['eta_line']: plot_objs['eta_line'].set_linewidth(val_lw)
        if plot_objs['p0_line']: plot_objs['p0_line'].set_linewidth(val_lw)
        if plot_objs['etad_line']: plot_objs['etad_line'].set_linewidth(val_lw)

        if val_interval <= 0: val_interval = 10
        if val_grid_n < 2: val_grid_n = 2
        
        # 3. Regenerate Current Quiver
        if plot_objs['current_quiver'] is not None:
             plot_objs['current_quiver'].remove()
             plot_objs['current_quiver'] = None
             
        if visibility['current'] and vc_vals is not None:
            # Recompute grid
            grid_x = np.linspace(0.05, 0.95, val_grid_n)
            grid_y = np.linspace(0.05, 0.95, val_grid_n)
            X_grid, Y_grid = np.meshgrid(grid_x, grid_y)
            
            U_field = np.full_like(X_grid, U_norm)
            V_field = np.full_like(Y_grid, V_norm)
            
            # Auto-scale width to preserve aspect ratio (Base scale 20)
            new_width = val_arr_width * (current_arrow_scale / val_arr_scale)
            
            plot_objs['current_quiver'] = ax.quiver(
                X_grid, Y_grid, U_field, V_field, 
                transform=ax.transAxes, 
                color=current_arrow_color, alpha=current_arrow_alpha, 
                pivot='mid', units='inches', width=new_width, 
                headwidth=current_arrow_headwidth, scale=val_arr_scale, scale_units='inches',
                zorder=1, label='Current'
            )
             
        
        # Calculate indices
        t_s, t_e = np.floor(t[0]), np.floor(t[-1])
        p_times = np.arange(t_s, t_e + 1, val_interval)
        
        p_indices = []
        for pt in p_times:
            idx = (np.abs(t - pt)).argmin()
            p_indices.append(idx)
        p_indices = sorted(list(set(p_indices)))
        if len(t) - 1 not in p_indices:
            p_indices.append(len(t) - 1)
            
        # Common Geometry
        xlim = ax.get_xlim()
        view_span = xlim[1] - xlim[0]
        scale_factor = val_scale / 1000.0 
        b_len = view_span * scale_factor
        s_wid = b_len / 2
        
        b_verts = np.array([[b_len/2, 0], [-b_len/2, s_wid/2], [-b_len/2, -s_wid/2]])

        # Update Vessel Triangles (Pose)
        for p in plot_objs['vessel_patches']:
            p.remove()
        plot_objs['vessel_patches'].clear()
        
        if visibility['pose']:
            for i, idx in enumerate(p_indices):
                n_p, e_p, h_p = x_north[idx], y_east[idx], psi[idx]
                
                # Rotate
                R = np.array([[np.cos(h_p), -np.sin(h_p)], [np.sin(h_p), np.cos(h_p)]])
                r_verts = (R @ b_verts.T).T 
                f_verts = np.column_stack((e_p + r_verts[:, 1], n_p + r_verts[:, 0]))
                
                # Alpha
                if len(p_indices) > 1:
                    progress = i / (len(p_indices) - 1)
                    a_val = val_amin + (val_amax - val_amin) * progress
                else:
                    a_val = val_amax
                    
                poly = patches.Polygon(f_verts, closed=True, 
                                       facecolor=vessel_face_color, edgecolor=vessel_edge_color, 
                                       linewidth=vessel_line_width, alpha=a_val, zorder=3)
                ax.add_patch(poly)
                plot_objs['vessel_patches'].append(poly)
        
        # Update Etad Triangles
        for p in plot_objs['etad_patches']:
            p.remove()
        plot_objs['etad_patches'].clear()

        if visibility['etad_sync'] and etad_available:
             for i, idx in enumerate(p_indices):
                if idx >= len(etad_north): continue

                n_d, e_d, h_d = etad_north[idx], etad_east[idx], etad_psi[idx]
                
                # Rotate using desired heading
                R = np.array([[np.cos(h_d), -np.sin(h_d)], [np.sin(h_d), np.cos(h_d)]])
                r_verts = (R @ b_verts.T).T 
                f_verts = np.column_stack((e_d + r_verts[:, 1], n_d + r_verts[:, 0]))
                
                # Alpha (Same logic as vessel)
                if len(p_indices) > 1:
                    progress = i / (len(p_indices) - 1)
                    a_val = val_amin + (val_amax - val_amin) * progress
                else:
                    a_val = val_amax
                
                # Use etad colors
                poly = patches.Polygon(f_verts, closed=True, 
                                       facecolor=etad_face_color, edgecolor=etad_edge_color, 
                                       linewidth=etad_line_width, alpha=a_val, zorder=6)
                ax.add_patch(poly)
                plot_objs['etad_patches'].append(poly)


        # Update Synchronized Markers (P0)
        if plot_objs['p0_markers']:
            # Ensure safe indexing
            if len(p_indices) > 0 and len(p0_north) == len(t):
                 m_n = p0_north[p_indices]
                 m_e = p0_east[p_indices]
                 plot_objs['p0_markers'].set_data(m_e, m_n)
            else:
                 plot_objs['p0_markers'].set_data([], [])

        
        # 3. Dynamic Legend Update
        handles = []
        labels = []
        
        def add_if_visible(key, label_override=None):
            obj = plot_objs.get(key)
            if obj and hasattr(obj, 'get_visible') and obj.get_visible():
                handles.append(obj)
                labels.append(label_override if label_override else obj.get_label())
        
        # Special handling for patches which are lists
        if visibility['pose'] and plot_objs['vessel_patches']:
             handles.append(plot_objs['vessel_patches'][-1])
             labels.append('Vessel Pose')

        if visibility['etad_sync'] and plot_objs['etad_patches']:
             handles.append(plot_objs['etad_patches'][-1])
             labels.append('etad Pose')

        add_if_visible('current_quiver', 'Current')
        add_if_visible('p0_line', r'$p_0$')
        add_if_visible('etad_line', r'$\eta_d$')
        add_if_visible('eta_line', r'$\eta$')
        add_if_visible('p0_markers', r'$p_0$ Sync')
        
        
        # Re-create legend
        if handles:
            ax.legend(handles, labels, loc='upper right', fontsize=FONT_SIZE)
        
        fig.canvas.draw_idle()


    # --- Controls Window Components ---
    # Compact Layout: 2 Columns
    
    # Header 1: General Settings
    fig_ctrl.text(0.5, 0.96, "General Settings", ha='center', fontsize=10, weight='bold')

    # 1. Title
    ax_tit = fig_ctrl.add_axes([0.1, 0.92, 0.8, 0.04])
    tb_title = TextBox(ax_tit, "Title: ", initial="Position and Heading")
    
    def on_title_submit(text):
        ax.set_title(text, fontsize=FONT_SIZE + 2, pad=20)
        fig.canvas.draw_idle()
    tb_title.on_submit(on_title_submit)
    
    # 2. Plot Parameters
    fig_ctrl.text(0.5, 0.88, "Plot Parameters", ha='center', fontsize=10, weight='bold')
    
    # Row 1: Interval | Scale
    ax_int = fig_ctrl.add_axes([0.30, 0.84, 0.15, 0.04])
    tb_interval = TextBox(ax_int, "Int [s]: ", initial=str(plot_interval_s))
    
    ax_scl = fig_ctrl.add_axes([0.75, 0.84, 0.20, 0.04])
    tb_scale = TextBox(ax_scl, "Scale: ", initial=str(ship_size_scale))

    # Row 2: Line Width | Axis Ratio
    ax_lw = fig_ctrl.add_axes([0.30, 0.79, 0.15, 0.04])
    tb_lw = TextBox(ax_lw, "Line Width: ", initial=str(eta_line_width))
    
    ax_ratio = fig_ctrl.add_axes([0.75, 0.79, 0.20, 0.04])
    tb_ratio = TextBox(ax_ratio, "Y/X Ratio: ", initial="1.0")

    # Row 3: Alpha Min | Alpha Max
    ax_amn = fig_ctrl.add_axes([0.30, 0.74, 0.15, 0.04])
    tb_amin = TextBox(ax_amn, "Min Alpha: ", initial="0.3")
    
    ax_amx = fig_ctrl.add_axes([0.75, 0.74, 0.20, 0.04])
    tb_amax = TextBox(ax_amx, "Max Alpha: ", initial="1.0")
    
    # Row 4: Grid N | Arr Scl | Arr Wid
    fig_ctrl.text(0.5, 0.69, "Current Arrows", ha='center', fontsize=9, weight='bold', color='gray')
    ax_grd = fig_ctrl.add_axes([0.25, 0.65, 0.15, 0.04])
    tb_grid_n = TextBox(ax_grd, "Grid N: ", initial="4")
    
    ax_arr = fig_ctrl.add_axes([0.55, 0.65, 0.15, 0.04])
    tb_arr_scale = TextBox(ax_arr, "Scl: ", initial=str(current_arrow_scale))

    ax_wid = fig_ctrl.add_axes([0.80, 0.65, 0.15, 0.04])
    tb_arr_width = TextBox(ax_wid, "Wid: ", initial=str(current_arrow_width))
    
    def on_param_submit(text):
        update_plot()
        
    def on_ratio_submit(text):
        try:
            r = float(text)
            ax.set_aspect(r)
            fig.canvas.draw_idle()
        except ValueError: pass
    
    tb_ratio.on_submit(on_ratio_submit)
    
    for tb in [tb_interval, tb_scale, tb_lw, tb_amin, tb_amax, tb_grid_n, tb_arr_scale, tb_arr_width]:
        tb.on_submit(on_param_submit)

    # 3. Actions & Visibility
    # Auto Scale
    ax_auto = fig_ctrl.add_axes([0.1, 0.58, 0.8, 0.05])
    b_auto = Button(ax_auto, 'Auto Scale View')
    
    def auto_scale_view(event):
        # Gather all visible points
        all_x = []
        all_y = []
        
        # 1. P0
        if plot_objs['p0_line'] and plot_objs['p0_line'].get_visible():
            all_x.extend(p0_east)
            all_y.extend(p0_north)
            
        # 2. Etad (Using markers or line points if we had full line data, here we have line data)
        if plot_objs['etad_line'] and plot_objs['etad_line'].get_visible() and etad_available:
             all_x.extend(etad_east)
             all_y.extend(etad_north)
             
        # 3. Eta (Vessel Path)
        if plot_objs['eta_line'] and plot_objs['eta_line'].get_visible():
            all_x.extend(y_east)
            all_y.extend(x_north)
            
        if not all_x or not all_y:
            return
            
        min_x, max_x = np.min(all_x), np.max(all_x)
        min_y, max_y = np.min(all_y), np.max(all_y)
        
        # Add padding (5%)
        pad_x = (max_x - min_x) * 0.05
        pad_y = (max_y - min_y) * 0.05
        
        if pad_x == 0: pad_x = 10
        if pad_y == 0: pad_y = 10
        
        ax.set_xlim(min_x - pad_x, max_x + pad_x)
        ax.set_ylim(min_y - pad_y, max_y + pad_y)
        fig.canvas.draw_idle()

    b_auto.on_clicked(auto_scale_view)

    fig_ctrl.text(0.5, 0.52, "Visibility", ha='center', fontsize=10, weight='bold')
    ax_chk = fig_ctrl.add_axes([0.1, 0.25, 0.5, 0.25]) 
    labels = ['etad', 'current', 'p0', 'pose', 'p0_sync', 'etad_sync']
    initial_status = [visibility[l] for l in labels]
    chk = CheckButtons(ax_chk, labels, actives=initial_status)
    def on_check(label):
        visibility[label] = not visibility[label]
        update_plot()
        
    chk.on_clicked(on_check)

    # 4. Save Button
    ax_save = fig_ctrl.add_axes([0.1, 0.10, 0.8, 0.08])
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
            
        out_filename = f"{script_name}_{input_name}.png"
        out_path = os.path.join(output_dir, out_filename)
        
        fig.savefig(out_path, dpi=300)
        print(f"Plot saved to: {out_path}")

    b_save.on_clicked(save_plot)
    
    # 5. Instructions
    fig_ctrl.text(0.5, 0.03, "Press Enter to Apply Changes", ha='center', fontsize=8, color='gray')


    # Initial Plot
    ax.set_title(tb_title.text, fontsize=FONT_SIZE + 2, pad=20)
    ax.set_xlabel("East [m]", fontsize=FONT_SIZE)
    ax.set_ylabel("North [m]", fontsize=FONT_SIZE)
    ax.grid(True, zorder=0)
    ax.axis('equal')
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE)
    
    update_plot()
    
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot trajectory WO from converted MAT file.")
    parser.add_argument("file", nargs='?', help="Path to _py.mat file")
    args = parser.parse_args()
    
    plot_trajectory(args.file)
