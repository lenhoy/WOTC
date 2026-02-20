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
    
    # Vessel Shape (Arrow)
    vessel_face_color = 'lightblue' 
    vessel_edge_color = 'steelblue'
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
    etad_marker_style = 'o'
    
    # Fonts
    FONT_SIZE = 15
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
        # psi = eta_vals[:, 2] # NOTE: Not used for heading in this script
    else:
        print("Error: 'eta' signal not found in data.")
        return

    # Extract p0 (Reference Trajectory) if available
    p0_available = False
    p0_north = None
    p0_east = None
    
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
    
    etad_struct = get_signal_struct(['etad', 'pd', 'pos_d'], data)
    
    if etad_struct:
        try:
            etad_vals = etad_struct['data']
            if etad_vals.ndim == 2:
                if etad_vals.shape[1] < 2 and etad_vals.shape[0] >= 2:
                     etad_vals = etad_vals.T
                
                etad_north = etad_vals[:, 0]
                etad_east = etad_vals[:, 1]
                etad_available = True
                print("Found 'etad' signal.")
        except Exception as e:
            print(f"Error processing 'etad': {e}")

    # --- Load Current Data (Beta) for Heading ---
    beta_vals = get_signal_data(['s_beta_current_', 'beta_current', 'beta_c'], data)
    
    # Ensure beta_vals covers the time range
    if beta_vals is not None:
        if np.ndim(beta_vals) == 0:
            beta_vals = np.full_like(t, beta_vals)
        elif len(beta_vals) != len(t):
             # Resize/Interpolate if needed (Naive resize for now as typically logsout is synced)
             if len(beta_vals) > len(t): beta_vals = beta_vals[:len(t)]
             else: beta_vals = np.resize(beta_vals, len(t))
    else:
        beta_vals = np.zeros_like(t) # Default to 0 if missing


    # --- 3. Configuration (Runtime Params) ---
    ship_size_scale = 150
    plot_interval_s = 60 # Plot every X seconds

    # --- 5. Setup Figures ---
    # Main Plot Figure
    fig, ax = plt.subplots(figsize=(6, 5))
    fig.canvas.manager.set_window_title("Trajectory Plot")
    plt.subplots_adjust(left=0.183) 
    
    # Controls Figure
    fig_ctrl = plt.figure(figsize=(5, 6)) # Reduced height
    fig_ctrl.canvas.manager.set_window_title("Controls")
    
    # Dictionary to store visibility states
    visibility = {
        'etad': True,
        'p0': True,
        'pose': True,
        'p0_sync': True,
        'etad_sync': True
    }
    
    # Objects storage
    plot_objs = {
        'p0_line': None,
        'p0_markers': None,
        'etad_line': None,
        'etad_markers': None,
        'eta_line': None,
        'vessel_patches': []
    }

    # --- 5b. Plot Static Lines (P0, Etad, Eta Path) ---
    
    # P0 Line
    if p0_available:
        p0_line, = ax.plot(p0_east, p0_north, color=p0_color, linewidth=p0_line_width, alpha=0.5, label='p0 Reference', zorder=2)
        plot_objs['p0_line'] = p0_line
        
        p0_markers, = ax.plot([], [], linestyle='None', marker=p0_marker_style, color=p0_marker_color, label='p0 Sync', zorder=3)
        plot_objs['p0_markers'] = p0_markers

    # Etad Line
    if etad_available:
        etad_line, = ax.plot(etad_east, etad_north, color=etad_color, linestyle=etad_style, linewidth=etad_width, label='etad Desired', zorder=2)
        plot_objs['etad_line'] = etad_line
        
        etad_markers, = ax.plot([], [], linestyle='None', marker=etad_marker_style, color=etad_color, alpha=0.5, markersize=4, label='etad Sync', zorder=3)
        plot_objs['etad_markers'] = etad_markers

    # Eta Path
    eta_line, = ax.plot(y_east, x_north, color=eta_line_color, linestyle='-', linewidth=eta_line_width, label='Vessel Path', zorder=2)
    plot_objs['eta_line'] = eta_line

    # --- Interactive Update Logic ---
    
    def update_plot(val=None):
        # 1. Update Visibility
        if plot_objs['p0_line']:
            plot_objs['p0_line'].set_visible(visibility['p0'])
        if plot_objs['p0_markers']:
            plot_objs['p0_markers'].set_visible(visibility['p0_sync'])
        if plot_objs['etad_line']:
            plot_objs['etad_line'].set_visible(visibility['etad'])
        if plot_objs['etad_markers']:
            plot_objs['etad_markers'].set_visible(visibility['etad_sync'])
        
        # 2. Update Parameters
        try:
            val_interval = float(tb_interval.text)
            val_scale = float(tb_scale.text)
            val_amin = float(tb_amin.text)
            val_amax = float(tb_amax.text)
            val_lw = float(tb_lw.text)
        except ValueError:
            return

        # Update Line Widths
        if plot_objs['eta_line']: plot_objs['eta_line'].set_linewidth(val_lw)
        if plot_objs['p0_line']: plot_objs['p0_line'].set_linewidth(val_lw)
        if plot_objs['etad_line']: plot_objs['etad_line'].set_linewidth(val_lw)

        if val_interval <= 0: val_interval = 10
        
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
            
        # Update Vessel Triangles (Pose)
        # Clear old patches
        for p in plot_objs['vessel_patches']:
            p.remove()
        plot_objs['vessel_patches'].clear()
        
        if visibility['pose']:
            # Geometry
            xlim = ax.get_xlim()
            view_span = xlim[1] - xlim[0]
            scale_factor = val_scale / 1000.0 
            b_len = view_span * scale_factor
            s_wid = b_len / 2
            
            # Arrow Shape Dimensions
            head_len = b_len * 0.4
            shaft_wid = s_wid * 0.3
            
            b_verts = np.array([
                [b_len/2, 0],                           # Tip
                [b_len/2 - head_len, s_wid/2],          # Head Top
                [b_len/2 - head_len, shaft_wid/2],      # Shaft Top Join
                [-b_len/2, shaft_wid/2],                # Shaft Tail Top
                [-b_len/2, -shaft_wid/2],               # Shaft Tail Bottom
                [b_len/2 - head_len, -shaft_wid/2],     # Shaft Bottom Join
                [b_len/2 - head_len, -s_wid/2]          # Head Bottom
            ])
            
            for i, idx in enumerate(p_indices):
                n_p, e_p = x_north[idx], y_east[idx]
                h_p = beta_vals[idx] # USING CURRENT DIRECTION
                
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
                                       linewidth=vessel_line_width, alpha=a_val, zorder=10)
                ax.add_patch(poly)
                plot_objs['vessel_patches'].append(poly)
        
        # Update Synchronized Markers (P0 / Etad)
        if plot_objs['p0_markers']:
            # Ensure safe indexing
            if len(p_indices) > 0 and len(p0_north) == len(t):
                 m_n = p0_north[p_indices]
                 m_e = p0_east[p_indices]
                 plot_objs['p0_markers'].set_data(m_e, m_n)
            else:
                 plot_objs['p0_markers'].set_data([], [])

        if plot_objs['etad_markers']:
             if etad_available and len(etad_north) == len(t) and len(p_indices) > 0:
                 m_n = etad_north[p_indices]
                 m_e = etad_east[p_indices]
                 plot_objs['etad_markers'].set_data(m_e, m_n)
             else:
                  plot_objs['etad_markers'].set_data([], [])
        
        # 3. Dynamic Legend Update
        handles = []
        labels = []
        
        def add_if_visible(key, label_override=None):
            obj = plot_objs.get(key)
            if obj and obj.get_visible():
                handles.append(obj)
                labels.append(label_override if label_override else obj.get_label())

        add_if_visible('p0_line', r'$p_0$')
        add_if_visible('etad_line', r'$\eta_d$')
        add_if_visible('eta_line', r'$\eta$')
        
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
    tb_title = TextBox(ax_tit, "Title: ", initial="Current Direction at Vessel")
    
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
    
    # Removed Row 4: Grid N | Arr Scl | Arr Wid (Background arrows removed)
    
    def on_param_submit(text):
        update_plot()
        
    def on_ratio_submit(text):
        try:
            r = float(text)
            ax.set_aspect(r)
            fig.canvas.draw_idle()
        except ValueError: pass
    
    tb_ratio.on_submit(on_ratio_submit)
    
    for tb in [tb_interval, tb_scale, tb_lw, tb_amin, tb_amax]:
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
            
        # 2. Etad
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
    labels = ['etad', 'p0', 'pose', 'p0_sync', 'etad_sync']
    initial_status = [visibility[l] for l in labels]
    chk = CheckButtons(ax_chk, labels, actives=initial_status)
    def on_check(label):
        visibility[label] = not visibility[label]
        update_plot()
        
    chk.on_clicked(on_check)

    # 4. Save Button
    ax_save = fig_ctrl.add_axes([0.1, 0.10, 0.8, 0.08]) # Adjusted Y prop since controls are fewer
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
    parser = argparse.ArgumentParser(description="Plot trajectory varying current alt2.")
    parser.add_argument("file", nargs='?', help="Path to _py.mat file")
    args = parser.parse_args()
    
    plot_trajectory(args.file)
