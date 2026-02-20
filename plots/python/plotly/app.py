import os
import sys
import dash
from dash import dcc, html, Input, Output, State, callback_context
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import base64
import io

# Add conversion folder to path to import loader
current_dir = os.path.dirname(os.path.abspath(__file__))
conversion_dir = os.path.join(current_dir, '../conversion')
sys.path.append(conversion_dir)

try:
    from load_mat_data import get_sim_data
except ImportError:
    # Fallback mock for development if file missing
    print("Warning: Could not import load_mat_data. Ensure '../conversion/load_mat_data.py' exists.")
    def get_sim_data(path): return {}

# Initialize App
app = dash.Dash(__name__, title="SimPlotter", update_title=None)
server = app.server

# Global Data Cache (Simple solution for local single-user app)
DATA_CACHE = {}

# --- Helper Functions ---
def get_available_files(root_dir):
    mat_files = []
    # limit search to reasonable depth or specific folder strings
    # For now, let's search specifically in 'output data' if we can find it
    # The structure seems to be: .../plots/python/plotly -> .../plots/python -> .../plots -> .../Matlab_simulasjon
    # Usually output data is in .../Matlab_simulasjon/output data
    
    # Let's try to find "output data" relative to this script
    base_proj_dir = os.path.abspath(os.path.join(current_dir, '../../../../')) # Adjust based on depth
    output_data_dir = os.path.join(base_proj_dir, 'output data')
    
    if not os.path.exists(output_data_dir):
        # Fallback to just scanning around
        output_data_dir = os.path.abspath(os.path.join(current_dir, '../../../'))

    for root, dirs, files in os.walk(output_data_dir):
        for file in files:
            if file.endswith('_py.mat'):
                full_path = os.path.join(root, file)
                # Make path relative for display if possible
                rel_path = os.path.relpath(full_path, output_data_dir)
                mat_files.append({'label': rel_path, 'value': full_path})
    
    # Sort by value (path) effectively sorts by date usually if folders are dated
    mat_files.sort(key=lambda x: x['value'], reverse=True) 
    return mat_files

def create_ghost_shape(x, y, psi, scale=1.0, color='rgba(100, 100, 100, 0.5)'):
    """
    Creates a vessel triangle shape.
    Shape: Simple isosceles triangle pointing at Psi.
    """
    # Base shape (Ned: x=North, y=East) -> but Plotly is X=East, Y=North usually
    # Let's align with Plotly coordinates: X=Horizontal, Y=Vertical.
    # If the data is NED:
    # Plot X = East (y_ned)
    # Plot Y = North (x_ned)
    # Heading Psi is compass (0=North, 90=East) using NED geometry.
    # To plot this in standard math (0=East, 90=North), we transform:
    # Math Angle = 90 - Psi (degrees)
    
    # Or we use the rotation matrix approach from the script
    # R = [[cos(psi), -sin(psi)], [sin(psi), cos(psi)]]
    # Vertices (Body): Bow=(L/2, 0), RearR=(-L/2, W/2), RearL=(-L/2, -W/2)
    # Note on size: Arbitrary units, scaled by 'scale'.
    
    length = 50 * scale
    width = 25 * scale
    
    # Vertices in body frame (nose pointing East for 0 angle in math, but let's stick to the script logic)
    # The script logic:
    # Bow: (base_length/2, 0)
    # Rear R: (-base_length/2, ship_width/2) --> This is technically 'Port' if Y is East? 
    # Let's assume Body X is forward.
    
    v_body = np.array([
        [length/2, 0],   # Nose
        [-length/2, width/2], # Rear Right
        [-length/2, -width/2], # Rear Left
    ])
    
    # Rotation Matrix for NED (Psi is angle from North (X) towards East (Y))
    # North = X_plot, East = Y_plot? NO. 
    # Usually: Plot X axis is East, Plot Y axis is North.
    # So:
    # Pos = [East, North]
    # Body X (Forward) maps to North * cos(psi) + East * sin(psi)
    # Wait, let's keep it simple.
    # Psi = 0 -> Pointing North (Up, +Y plot)
    # Psi = 90 -> Pointing East (Right, +X plot)
    
    # Rotation angle for Plotly (standard math, 0=Right):
    # If Psi=0 (North) -> Math=90.
    # If Psi=90 (East) -> Math=0.
    # Angle = 90 - degrees(psi)
    
    # Using manual rotation to get vertices relative to (0,0) then translate
    # Math Angle
    theta = np.radians(90) - psi # Psi in radians
    
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))
    
    v_rot = v_body.dot(R.T) # (N, 2)
    
    # Translate
    v_final = v_rot + np.array([x, y])
    
    path = f"M {v_final[0,0]},{v_final[0,1]} L {v_final[1,0]},{v_final[1,1]} L {v_final[2,0]},{v_final[2,1]} Z"
    
    return dict(
        type="path",
        path=path,
        fillcolor=color,
        line=dict(color=color, width=1),
        layer="above"
    )

# --- Layout ---
app.layout = html.Div(style={'display': 'flex', 'height': '100vh', 'fontFamily': 'sans-serif', 'backgroundColor': '#f4f4f4'}, children=[
    # Sidebar
    html.Div(className='sidebar', style={
        'width': '350px', 
        'padding': '20px', 
        'backgroundColor': '#2c3e50', 
        'color': '#ecf0f1',
        'display': 'flex',
        'flexDirection': 'column',
        'overflowY': 'auto',
        'boxShadow': '2px 0 5px rgba(0,0,0,0.1)'
    }, children=[
        html.H2("SimPlotter", style={'textAlign': 'center', 'marginBottom': '20px', 'color': '#3498db'}),
        
        # 1. File Selection
        html.Label("1. Select Simulation File:", style={'fontWeight': 'bold'}),
        dcc.Dropdown(
            id='file-dropdown',
            options=get_available_files('.'), 
            placeholder="Select a _py.mat file...",
            style={'color': '#000', 'marginBottom': '15px'}
        ),
        html.Button('Refresh Files', id='refresh-files-btn', n_clicks=0, style={'width': '100%', 'marginBottom': '20px', 'padding': '5px'}),
        
        html.Hr(style={'borderColor': '#7f8c8d'}),
        
        # 2. Data Mapping
        html.Div(id='data-controls', style={'display': 'none'}, children=[
            html.Label("2. Axis Mapping:", style={'fontWeight': 'bold'}),
            
            html.Label("X-Axis (East):", style={'fontSize': '12px'}),
            dcc.Dropdown(id='x-axis-col', style={'color': '#000', 'marginBottom': '5px'}),
            
            html.Label("Y-Axis (North):", style={'fontSize': '12px'}),
            dcc.Dropdown(id='y-axis-col', style={'color': '#000', 'marginBottom': '5px'}),
            
            html.Label("Psi (Heading) [rad]:", style={'fontSize': '12px'}),
            dcc.Dropdown(id='psi-col', style={'color': '#000', 'marginBottom': '15px'}),
            
            html.Hr(style={'borderColor': '#7f8c8d'}),
            
            # 3. Plot Settings
            html.Label("3. Plot Settings:", style={'fontWeight': 'bold'}),
            
            html.Label("Ghost Interval [s]:", style={'fontSize': '12px'}),
            dcc.Slider(id='ghost-interval', min=10, max=600, step=10, value=60, 
                       marks={60: '60s', 300: '5m', 600: '10m'}, tooltip={'placement': 'bottom', 'always_visible': True}),
            
            html.Label("Ship Scale:", style={'fontSize': '12px', 'marginTop': '10px'}),
            dcc.Slider(id='ship-scale', min=0.1, max=5, step=0.1, value=1.0, 
                       marks={1: '1x', 5: '5x'}, tooltip={'placement': 'bottom', 'always_visible': True}),
            
            html.Div([
                dcc.Checklist(
                    id='plot-options',
                    options=[
                        {'label': ' Show Grid', 'value': 'grid'},
                        {'label': ' Equal Axis', 'value': 'equal'},
                        {'label': ' Show Ghosts', 'value': 'ghosts'},
                    ],
                    value=['grid', 'equal', 'ghosts'],
                    labelStyle={'display': 'block', 'cursor': 'pointer'}
                )
            ], style={'marginTop': '15px'}),
             
            html.Label("Line Color:", style={'fontSize': '12px', 'marginTop': '10px'}),
            dcc.Input(id='line-color', type='text', value='#3498db', style={'width': '100%'}),
            
            html.Hr(style={'borderColor': '#7f8c8d'}),
            
            # 4. Export
            html.Label("4. Export:", style={'fontWeight': 'bold'}),
             dcc.Input(id='export-filename', type='text', placeholder='plot_export', value='my_plot', style={'width': '100%', 'marginBottom': '5px'}),
             html.Button('Download Image', id='download-btn', style={'width': '100%', 'padding': '10px', 'backgroundColor': '#27ae60', 'color': 'white', 'border': 'none', 'cursor': 'pointer'}),
             dcc.Download(id="download-image")
        ]),
        
        html.Div(id='status-msg', style={'marginTop': 'auto', 'fontSize': '12px', 'color': '#bdc3c7'})
    ]),
    
    # Main Content
    html.Div(style={'flex': '1', 'padding': '20px', 'display': 'flex', 'flexDirection': 'column'}, children=[
        dcc.Graph(
            id='main-graph', 
            style={'height': '100%', 'width': '100%'},
            config={
                'scrollZoom': True, 
                'displayModeBar': True,
                'toImageButtonOptions': {
                    'format': 'png', # one of png, svg, jpeg, webp
                    'filename': 'custom_image',
                    'height': 1000,
                    'width': 1200,
                    'scale': 2 # Multiply title/legend/axis/canvas sizes by this factor
                }
            }
        )
    ])
])

# --- Callbacks ---

@app.callback(
    [Output('file-dropdown', 'options'),
     Output('status-msg', 'children')],
    [Input('refresh-files-btn', 'n_clicks')]
)
def refresh_file_list(n):
    files = get_available_files('.')
    msg = f"Found {len(files)} files."
    return files, msg

@app.callback(
    [Output('data-controls', 'style'),
     Output('x-axis-col', 'options'),
     Output('y-axis-col', 'options'),
     Output('psi-col', 'options'),
     Output('x-axis-col', 'value'),
     Output('y-axis-col', 'value'),
     Output('psi-col', 'value')],
    [Input('file-dropdown', 'value')]
)
def load_and_populate_columns(file_path):
    if not file_path:
        return {'display': 'none'}, [], [], [], None, None, None
        
    # Load Data
    try:
        # Use our loaded script
        # Check cache
        if file_path in DATA_CACHE:
            dfs = DATA_CACHE[file_path]
        else:
            dfs = get_sim_data(file_path)
            DATA_CACHE[file_path] = dfs
            
        # Extract available columns across all dataframes
        # We flatten the structure: "SignalName.ColumnName"
        options = []
        for sig_name, df in dfs.items():
            for col in df.columns:
                options.append({'label': f"{sig_name} : {col}", 'value': f"{sig_name}|{col}"})
        
        # Smart Defaults
        def find_default(keywords):
            for opt in options:
                val = opt['value'].lower()
                for k in keywords:
                    if k in val:
                        return opt['value']
            return None
        
        # Default X (East), Y (North), Psi
        # Looking for 'eta_2' or 'East' or 'y'
        default_x = find_default(['eta_2', 'east', 'y_pos'])
        default_y = find_default(['eta_1', 'north', 'x_pos'])
        default_psi = find_default(['eta_3', 'psi', 'yaw'])
        
        return {'display': 'block'}, options, options, options, default_x, default_y, default_psi
        
    except Exception as e:
        print(f"Error loading file: {e}")
        return {'display': 'none'}, [], [], [], None, None, None

@app.callback(
    Output('main-graph', 'figure'),
    [Input('x-axis-col', 'value'),
     Input('y-axis-col', 'value'),
     Input('psi-col', 'value'),
     Input('ghost-interval', 'value'),
     Input('ship-scale', 'value'),
     Input('plot-options', 'value'),
     Input('line-color', 'value'),
     Input('file-dropdown', 'value')]
)
def update_graph(x_val, y_val, psi_val, ghost_interval, ship_scale, options, line_color, file_path):
    if not file_path or not x_val or not y_val or file_path not in DATA_CACHE:
        return go.Figure()
        
    dfs = DATA_CACHE[file_path]
    
    # helper to fetch series
    def get_series(selection):
        sig, col = selection.split('|')
        return dfs[sig][col]
        
    try:
        x_series = get_series(x_val)
        y_series = get_series(y_val)
        
        # Align timestamps (Assuming they are mostly aligned if from same sim)
        # Ideally we use the time index of one and reindex the others if needed
        # For now, assumes exact alignment or same DataFrame origin which is likely
        
        fig = go.Figure()
        
        # 1. Main Trajectory
        fig.add_trace(go.Scatter(
            x=x_series,
            y=y_series,
            mode='lines',
            name='Trajectory',
            line=dict(color=line_color, width=2)
        ))
        
        # 2. Ghosts
        if 'ghosts' in options and psi_val:
            psi_series = get_series(psi_val)
            
            # Resample for ghosts
            # Get time range
            t = x_series.index 
            t_start = np.floor(t.min())
            t_end = np.floor(t.max())
            
            ghost_times = np.arange(t_start, t_end, ghost_interval)
            
            # Find indices
            shapes = []
            for gt in ghost_times:
                idx = (np.abs(t - gt)).argmin()
                
                # Get Values using positional index
                xi = x_series.iloc[idx]
                yi = y_series.iloc[idx]
                psii = psi_series.iloc[idx]
                
                # Create Shape
                shape = create_ghost_shape(xi, yi, psii, scale=ship_scale)
                shapes.append(shape)
                
            fig.update_layout(shapes=shapes)
            
        # 3. Layout Update
        fig.update_layout(
            title=f"Simulation Analysis: {os.path.basename(file_path)}",
            xaxis_title=x_val.split('|')[1],
            yaxis_title=y_val.split('|')[1],
            template="plotly_white",
            margin=dict(l=40, r=40, t=50, b=40),
            hovermode='closest'
        )
        
        if 'grid' in options:
            fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
            fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
        else:
             fig.update_xaxes(showgrid=False)
             fig.update_yaxes(showgrid=False)
             
        if 'equal' in options:
            fig.update_yaxes(scaleanchor="x", scaleratio=1)
            
    except Exception as e:
        print(f"Plot Logic Error: {e}")
        return go.Figure()
        
    return fig

# Manual Download Handler (Uses dcc.Download)
# Note: Plotly's camera button is client-side. 
# If we want a specific server-side high-res render pipeline we need kaleido.
# For now, we rely on the client-side 'toImageButtonOptions' configured in the Graph component.
# But I added a button that could ostensibly do more if needed.
# Let's make the button just trigger a client-side download via JS or advanced callback if the user requested "clean ways".
# The build-in camera button is usually cleanest for SVG/High-Res PNG.
# I will leave the button as a placeholder or remove it if I can't serve a better file easily without kaleido.
# Actually, let's keep it simple: the built-in ModeBar button is configured to give high-res. 

if __name__ == '__main__':
    app.run(debug=True, port=8050)
