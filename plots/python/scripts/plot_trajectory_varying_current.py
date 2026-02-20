import scipy.io
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
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
    p0_color = 'grey'
    p0_marker_style = 'x'
    p0_marker_color = 'black'
    
    # Etad
    etad_color = 'green'
    etad_style = '--'
    etad_width = 1.5
    etad_marker_style = 'o'
    
    # Background Current (Arrows)
    current_arrow_color = 'cyan'
    current_arrow_alpha_base = 0.6
    current_arrow_width_scale = 0.005 # Width relative to plot
    
    # Fonts
    FONT_SIZE = 12
    plt.rcParams.update({'font.size': FONT_SIZE})

    # --- 1. Handle Input ---
    if not file_path:
        root = Tk()
        root.withdraw() 
        file_path = filedialog.askopenfilename(
            title="Select Converted Python Data File (_py.mat)",
            filetypes=[("MAT files", "*_py.mat")]
        )
        root.destroy()
        
    if not file_path or not isinstance(file_path, str):
        print("No file selected.")
        return

    print(f"Loading: {file_path}")
    
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return

    # --- 2. Load Data ---
    try:
        mat = scipy.io.loadmat(file_path, squeeze_me=True, struct_as_record=False, simplify_cells=True)
    except Exception as e:
        print(f"Error loading .mat file: {e}")
        return

    if 'data' not in mat:
        print("Error: Could not find 'data' struct in file.")
        return
        
    data = mat['data']
    
    if not isinstance(data, dict):
        try:
            field_names = data.dtype.names
            data = {name: data[name] for name in field_names}
        except: pass

    def get_signal_data(keys, source):
        for k in keys:
            if k in source:
                val = source[k]
                if hasattr(val, 'data'): return val.data
                if isinstance(val, dict) and 'data' in val: return val['data']
                return val
        return None

    # Extract Eta
    if 'eta' in data:
        eta_struct = data['eta']
        t = eta_struct['time']
        eta_vals = eta_struct['data']
        if eta_vals.shape[1] != 3 and eta_vals.shape[0] == 3: eta_vals = eta_vals.T
        y_east = eta_vals[:, 1]
        x_north = eta_vals[:, 0]
        psi = eta_vals[:, 2] # Radians
    else:
        print("Error: 'eta' signal not found in data.")
        return

    # Extract p0
    p0_available = False
    p0_north, p0_east = None, None
    p0_struct = None
    for k in ['p0', 'ref_pos']:
        if k in data: p0_struct = data[k]; break
            
    if p0_struct:
        try:
            p0_vals = p0_struct['data']
            if p0_vals.ndim == 2:
                if p0_vals.shape[1] < 2 and p0_vals.shape[0] >= 2: p0_vals = p0_vals.T
                p0_north, p0_east = p0_vals[:, 0], p0_vals[:, 1]
                p0_available = True
        except: pass

    # Extract etad
    etad_available = False
    etad_north, etad_east = None, None
    etad_struct = None
    for k in ['etad', 'pd', 'pos_d']:
        if k in data: etad_struct = data[k]; break
            
    if etad_struct:
        try:
            etad_vals = etad_struct['data']
            if etad_vals.ndim == 2:
                if etad_vals.shape[1] < 2 and etad_vals.shape[0] >= 2: etad_vals = etad_vals.T
                etad_north, etad_east = etad_vals[:, 0], etad_vals[:, 1]
                etad_available = True
        except: pass

    # Extract Current
    # We need time-varying current
    vc_vals = get_signal_data(['s_V_current_', 'V_current', 'vc'], data)
    beta_vals = get_signal_data(['s_beta_current_', 'beta_current', 'beta_c'], data)
    
    # Ensure they are arrays matching time t, or interpolate
    if vc_vals is not None and beta_vals is not None:
        if np.ndim(vc_vals) == 0: vc_vals = np.full_like(t, vc_vals)
        if np.ndim(beta_vals) == 0: beta_vals = np.full_like(t, beta_vals)
        
        if len(vc_vals) != len(t):
             if len(vc_vals) > len(t): vc_vals = vc_vals[:len(t)]
             else: vc_vals = np.resize(vc_vals, len(t))
             
        if len(beta_vals) != len(t):
             if len(beta_vals) > len(t): beta_vals = beta_vals[:len(t)]
             else: beta_vals = np.resize(beta_vals, len(t))
    else:
        # Defaults
        vc_vals = np.zeros_like(t)
        beta_vals = np.zeros_like(t)


    # --- 3. Configuration ---
    ship_size_scale = 150
    plot_interval_s = 60
    current_circle_radius_scale = 1.5 # Radius relative to ship size

    # --- 4. Figures ---
    # Setup GridSpec for Main Figure: Top (Trajectory), Bottom (Heading)
    fig = plt.figure(figsize=(8, 10))
    fig.canvas.manager.set_window_title("Trajectory Time-Varying Current")
    
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.3)
    ax = fig.add_subplot(gs[0])
    ax_head = fig.add_subplot(gs[1])
    
    # Controls Figure
    # Widen controls window for descriptive labels
    fig_ctrl = plt.figure(figsize=(5.5, 7.5)) # Increased size
    fig_ctrl.canvas.manager.set_window_title("Controls")
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    # --- State Storage ---
    visibility = {
        'etad': True,
        'current': True,
        'p0': True,
        'pose': True,
        'p0_sync': True,
        'etad_sync': True
    }
    
    plot_objs = {
        'p0_line': None,
        'p0_markers': None,
        'etad_line': None,
        'etad_markers': None,
        'eta_line': None,
        'vessel_patches': [],
        'current_arrows': [],
        'rose_line': None,
        'rose_scatter': None
    }

    # --- Initial Static Plotting ---
    
    # 1. Trajectory Plot (Static Parts)
    if p0_available:
        l, = ax.plot(p0_east, p0_north, color=p0_color, linewidth=p0_line_width, alpha=0.5, label='p0 Ref', zorder=2)
        plot_objs['p0_line'] = l
        m, = ax.plot([], [], linestyle='None', marker=p0_marker_style, color=p0_marker_color, label='p0 Sync', zorder=3)
        plot_objs['p0_markers'] = m

    if etad_available:
        l, = ax.plot(etad_east, etad_north, color=etad_color, linestyle=etad_style, linewidth=etad_width, label='etad Desired', zorder=2)
        plot_objs['etad_line'] = l
        m, = ax.plot([], [], linestyle='None', marker=etad_marker_style, color=etad_color, alpha=0.5, markersize=4, label='etad Sync', zorder=3)
        plot_objs['etad_markers'] = m

    l, = ax.plot(y_east, x_north, color=eta_line_color, linestyle='-', linewidth=eta_line_width, label='Vessel Path', zorder=2)
    plot_objs['eta_line'] = l
    
    ax.set_xlabel("East [m]")
    ax.set_ylabel("North [m]")
    ax.axis('equal')
    ax.grid(True)
    ax.set_title("Trajectory & Current")

    # 2. Current Direction Subplot (Beta)
    beta_deg = np.degrees(beta_vals)
    beta_unwrapped = np.degrees(np.unwrap(beta_vals))
    
    ax_head.plot(t, beta_unwrapped, color='cyan', linewidth=1.5, label=r'$\beta_c$ (Current Dir)')
    
    ax_head.set_xlabel("Time [s]")
    ax_head.set_ylabel("Beta [deg]")
    ax_head.grid(True)
    ax_head.legend(loc='upper right', fontsize=10)
    ax_head.set_title("Current Direction History")

    # 3. Compass Rose Init
    # Removed as per user request

    # --- Update Logic ---
    
    def update_plot(val=None):
        # 1. Parameters
        try:
            val_interval = float(tb_interval.text)
            val_scale = float(tb_scale.text)
            val_amin = float(tb_amin.text)
            val_amax = float(tb_amax.text)
            val_arr_scale = float(tb_arr_scale.text) 
            val_arrow_len = float(tb_arr_len.text)   
            val_lw = float(tb_lw.text)
            val_arr_width = float(tb_arr_width.text) # New: Independent Arrow Width
        except ValueError: return

        if val_interval <= 0: val_interval = 10
        
        # Update Line Widths
        if plot_objs['eta_line']: plot_objs['eta_line'].set_linewidth(val_lw)
        if plot_objs['p0_line']: plot_objs['p0_line'].set_linewidth(val_lw)
        if plot_objs['etad_line']: plot_objs['etad_line'].set_linewidth(val_lw)
        
        # 2. Indices
        t_s, t_e = np.floor(t[0]), np.floor(t[-1])
        p_times = np.arange(t_s, t_e + 1, val_interval)
        p_indices = [(np.abs(t - pt)).argmin() for pt in p_times]
        p_indices = sorted(list(set(p_indices)))
        
        # 3. Update Vessel Ghosts & Arrows
        for p in plot_objs['vessel_patches']: p.remove()
        for a in plot_objs['current_arrows']: a.remove()
        plot_objs['vessel_patches'].clear()
        plot_objs['current_arrows'].clear()
        
        if visibility['pose']:
            # Ship Dimensions
            xlim = ax.get_xlim()
            
            L_ship = val_scale # e.g. 100 meters
            W_ship = L_ship * 0.5 # Aspect ratio
            
            b_verts = np.array([[L_ship/2, 0], [-L_ship/2, W_ship/2], [-L_ship/2, -W_ship/2]])
            
            # Circle Radius for Current Arrow
            R_circle = L_ship * val_arr_scale
            
            for i, idx in enumerate(p_indices):
                # Fade
                if len(p_indices) > 1:
                    progress = i / (len(p_indices) - 1)
                    a_val = val_amin + (val_amax - val_amin) * progress
                else: a_val = val_amax
                
                n_p, e_p, h_p = x_north[idx], y_east[idx], psi[idx]
                
                # Vessel
                R = np.array([[np.cos(h_p), -np.sin(h_p)], [np.sin(h_p), np.cos(h_p)]])
                r_verts = (R @ b_verts.T).T
                f_verts = np.column_stack((e_p + r_verts[:, 0], n_p + r_verts[:, 1]))
                
                poly = patches.Polygon(f_verts, closed=True, facecolor=vessel_face_color, 
                                       edgecolor=vessel_edge_color, linewidth=vessel_line_width, alpha=a_val, zorder=10)
                ax.add_patch(poly)
                plot_objs['vessel_patches'].append(poly)
                
                # Current Arrow
                if visibility['current']:
                    vc = vc_vals[idx]
                    beta = beta_vals[idx] 
                    
                    tail_ang = beta + np.pi 
                    t_x = e_p + R_circle * np.cos(tail_ang)
                    t_y = n_p + R_circle * np.sin(tail_ang)
                    
                    arrow_len = vc * val_arrow_len 
                    dx = arrow_len * np.cos(beta)
                    dy = arrow_len * np.sin(beta)
                    
                    # Independent Width Scaling
                    width = val_arr_width
                    head_width = width * 3
                    
                    arr = patches.FancyArrow(t_x, t_y, dx, dy, 
                                             width=width, 
                                             length_includes_head=True, 
                                             head_width=head_width,
                                             color=current_arrow_color, alpha=a_val, zorder=11)
                    ax.add_patch(arr)
                    plot_objs['current_arrows'].append(arr)
                    
        # 4. Sync Markers
        if plot_objs['p0_markers']:
            if visibility['p0_sync'] and len(p_indices)>0:
                plot_objs['p0_markers'].set_data(p0_east[p_indices], p0_north[p_indices])
                plot_objs['p0_markers'].set_visible(True)
            else: plot_objs['p0_markers'].set_visible(False)
            plot_objs['p0_line'].set_visible(visibility['p0'])

        if plot_objs['etad_markers']:
            if visibility['etad_sync'] and len(p_indices)>0:
                 plot_objs['etad_markers'].set_data(etad_east[p_indices], etad_north[p_indices])
                 plot_objs['etad_markers'].set_visible(True)
            else: plot_objs['etad_markers'].set_visible(False)
            plot_objs['etad_line'].set_visible(visibility['etad'])
            
        fig.canvas.draw_idle()

    # --- Controls Layout ---
    # Title
    ax_tit = fig_ctrl.add_axes([0.1, 0.94, 0.8, 0.04])
    tb_title = TextBox(ax_tit, "Plot Title: ", initial="Trajectory & Current")
    def on_title(start): 
        ax.set_title(start, fontsize=14)
        fig.canvas.draw_idle()
    tb_title.on_submit(on_title)
    
    # Params
    fig_ctrl.text(0.5, 0.90, "General Parameters", ha='center', weight='bold')
    
    # Row 1
    ax_int = fig_ctrl.add_axes([0.45, 0.86, 0.2, 0.04])
    tb_interval = TextBox(ax_int, "Plot Interval [s]: ", initial=str(plot_interval_s), label_pad=0.1)
    
    ax_scl = fig_ctrl.add_axes([0.45, 0.81, 0.2, 0.04])
    tb_scale = TextBox(ax_scl, "Ship Length [m]: ", initial=str(10.0), label_pad=0.1)

    # Row 2
    ax_lw = fig_ctrl.add_axes([0.45, 0.76, 0.2, 0.04])
    tb_lw = TextBox(ax_lw, "Line Width: ", initial="1.5", label_pad=0.1)
    
    ax_asp = fig_ctrl.add_axes([0.45, 0.71, 0.2, 0.04])
    tb_ratio = TextBox(ax_asp, "Aspect Ratio (Y/X): ", initial="1.0", label_pad=0.1)

    # Row 3 (Opacity)
    fig_ctrl.text(0.5, 0.67, "Opacity (Fading)", ha='center', weight='bold')
    
    ax_amn = fig_ctrl.add_axes([0.35, 0.63, 0.15, 0.04])
    tb_amin = TextBox(ax_amn, "Min Opacity: ", initial="0.3", label_pad=0.1)
    
    ax_amx = fig_ctrl.add_axes([0.75, 0.63, 0.15, 0.04])
    tb_amax = TextBox(ax_amx, "Max Opacity: ", initial="1.0", label_pad=0.1)
    
    # Row 4 (Current)
    fig_ctrl.text(0.5, 0.58, "Current Arrow Config", ha='center', weight='bold')
    
    ax_rad = fig_ctrl.add_axes([0.55, 0.54, 0.2, 0.04])
    tb_arr_scale = TextBox(ax_rad, "Arrow Dist (xLen): ", initial=str(current_circle_radius_scale), label_pad=0.1)
    
    ax_len = fig_ctrl.add_axes([0.55, 0.49, 0.2, 0.04])
    tb_arr_len = TextBox(ax_len, "Arrow Mag Scale: ", initial="50.0", label_pad=0.1)
    
    ax_wid = fig_ctrl.add_axes([0.55, 0.44, 0.2, 0.04])
    tb_arr_width = TextBox(ax_wid, "Arrow Thickness: ", initial=str(ship_size_scale * 0.005), label_pad=0.1) # Approx default

    # Plot Update Callback
    def submit_all(text): update_plot()
    
    # Aspect Ratio Callback
    def on_ratio(text):
        try:
            r = float(text)
            ax.set_aspect(r)
            fig.canvas.draw_idle()
        except: pass

    tb_ratio.on_submit(on_ratio)
    for tb in [tb_interval, tb_scale, tb_lw, tb_amin, tb_amax, tb_arr_scale, tb_arr_len, tb_arr_width]:
        tb.on_submit(submit_all)

    # Actions: Auto Scale
    ax_auto = fig_ctrl.add_axes([0.1, 0.35, 0.8, 0.06])
    b_auto = Button(ax_auto, "Auto Scale View")
    
    def auto_scale_view(event):
        all_x, all_y = [], []
        
        if plot_objs['p0_line'] and plot_objs['p0_line'].get_visible() and p0_available:
            all_x.extend(p0_east); all_y.extend(p0_north)
            
        if plot_objs['etad_line'] and plot_objs['etad_line'].get_visible() and etad_available:
             all_x.extend(etad_east); all_y.extend(etad_north)
             
        if plot_objs['eta_line'] and plot_objs['eta_line'].get_visible():
            all_x.extend(y_east); all_y.extend(x_north)
            
        if not all_x: return
            
        min_x, max_x = np.min(all_x), np.max(all_x)
        min_y, max_y = np.min(all_y), np.max(all_y)
        pad_x = (max_x - min_x) * 0.05
        pad_y = (max_y - min_y) * 0.05
        if pad_x==0: pad_x=10
        if pad_y==0: pad_y=10
        
        ax.set_xlim(min_x - pad_x, max_x + pad_x)
        ax.set_ylim(min_y - pad_y, max_y + pad_y)
        fig.canvas.draw_idle()
        
    b_auto.on_clicked(auto_scale_view)

    # Visibility
    fig_ctrl.text(0.5, 0.30, "Visibility", ha='center', weight='bold')
    ax_chk = fig_ctrl.add_axes([0.1, 0.10, 0.5, 0.20])
    labels = ['etad', 'current', 'p0', 'pose', 'p0_sync', 'etad_sync']
    chk = CheckButtons(ax_chk, labels, actives=[visibility[l] for l in labels])
    
    def on_chk(lbl):
        visibility[lbl] = not visibility[lbl]
        update_plot()
    chk.on_clicked(on_chk)
    
    # Save
    ax_sav = fig_ctrl.add_axes([0.1, 0.02, 0.8, 0.08])
    b_sav = Button(ax_sav, "Save Plot")
    def save(e):
        out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'output')
        if not os.path.exists(out_dir): os.makedirs(out_dir)
        fig.savefig(os.path.join(out_dir, f"traj_curr_{os.path.basename(file_path[:-7])}.png"), dpi=300)
        print("Saved.")
    b_sav.on_clicked(save)

    update_plot()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", nargs='?')
    args = parser.parse_args()
    plot_trajectory(args.file)
