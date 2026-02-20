import scipy.io as spio
import pandas as pd
import numpy as np
import os

def loadmat(filename):
    """
    This function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects.
    
    Source: https://scipy-cookbook.readthedocs.io/items/Reading_mat_files.html
    """
    try:
        data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None
    return _check_keys(data)

def _check_keys(dict):
    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def load_sim_data(file_path):
    """
    Loads a .mat file (specifically those converted by convert_for_python.m)
    and returns a dictionary of Pandas DataFrames.
    
    Args:
        file_path (str): Path to the .mat file.
        
    Returns:
        dict: A dictionary where keys are signal names and values are pd.DataFrames
              with a 'Time' index.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    print(f"Loading: {file_path} ...")
    raw_data = loadmat(file_path)
    
    if raw_data is None:
        return {}
        
    dfs = {}
    
    # Identify the container for signals.
    # convert_for_python.m saves signals inside a 'data' struct.
    # If using formatted file:
    source_data = None
    if 'data' in raw_data:
        source_data = raw_data['data']
    elif 'logsout' in raw_data:
        # Fallback if someone tried to load a raw file that happened to be a struct
        source_data = raw_data['logsout']
    else:
        # Maybe the file *is* the dictionary of signals (flat save)
        source_data = raw_data

    # Iterate through signals
    for name, signal in source_data.items():
        # Clean metadata keys like __header__, __version__, etc.
        if name.startswith('__'):
            continue
            
        # Check if it looks like a signal (has 'data' and 'time')
        if isinstance(signal, dict) and 'data' in signal and 'time' in signal:
            try:
                data_vals = signal['data']
                time_vals = signal['time']
                
                # Check consistency
                # 1D array -> Series or DataFrame
                # 2D array -> DataFrame
                
                # Make sure time is 1D
                if hasattr(time_vals, 'ndim') and time_vals.ndim > 1:
                    time_vals = time_vals.flatten()
                
                # Handle Data
                if isinstance(data_vals, (np.ndarray, list)):
                    data_arr = np.array(data_vals)
                    
                    # If data is 2D (N x Dim) and N matches Time
                    if data_arr.ndim == 2 and data_arr.shape[0] == len(time_vals):
                        columns = [f"{name}_{i+1}" for i in range(data_arr.shape[1])]
                        df = pd.DataFrame(data_arr, index=time_vals, columns=columns)
                    
                    # If data is 1D and matches Time
                    elif data_arr.ndim == 1 and len(data_arr) == len(time_vals):
                         df = pd.DataFrame(data_arr, index=time_vals, columns=[name])
                         
                    # Transposed case? (Dim x N)
                    elif data_arr.ndim == 2 and data_arr.shape[1] == len(time_vals):
                         data_arr = data_arr.T
                         columns = [f"{name}_{i+1}" for i in range(data_arr.shape[1])]
                         df = pd.DataFrame(data_arr, index=time_vals, columns=columns)
                    
                    else:
                        # Mismatched lengths, skip
                        # print(f"Skipping {name}: Dimension mismatch. Time {len(time_vals)}, Data {data_arr.shape}")
                        continue
                        
                    df.index.name = 'Time'
                    dfs[name] = df
                    
            except Exception as e:
                print(f"Failed to convert {name} to DataFrame: {e}")
                
    return dfs

def get_sim_data(raw_file_path):
    """
    End-to-end loader that accepts the RAW .mat file path (typical simulation output).
    It automatically looks for the converted '_py.mat' version.
    
    Rules:
    1. If raw_file_path doesn't exist, raise Error.
    2. Construct py_file_path = raw_file_path + '_py.mat' (ignoring extension replacement for safety, or smart replace).
    3. Check if py_file_path exists.
    4. If not, raise specific error instructing user to run 'convert_for_python.m'.
    5. If exists, check timestamps. If raw is newer than py, warn user.
    6. Load using load_sim_data.
    """
    if not os.path.exists(raw_file_path):
        raise FileNotFoundError(f"Raw simulation file not found: {raw_file_path}")
    
    # Construct expected converted filename
    # Strategy: /path/to/file.mat -> /path/to/file_py.mat
    folder, filename = os.path.split(raw_file_path)
    if filename.lower().endswith('.mat'):
        basename = filename[:-4]
    else:
        basename = filename
        
    py_filename = f"{basename}_py.mat"
    py_file_path = os.path.join(folder, py_filename)
    
    if not os.path.exists(py_file_path):
        raise FileNotFoundError(
            f"Converted file missing: {py_file_path}\n"
            f"Please run the MATLAB script 'plots/python/conversion/convert_for_python.m' \n"
            f"on your raw file: '{filename}' first."
        )
        
    # Check timestamps
    raw_mtime = os.path.getmtime(raw_file_path)
    py_mtime = os.path.getmtime(py_file_path)
    
    if raw_mtime > py_mtime:
        print("WARNING: The raw simulation file is NEWER than the converted python file.")
        print(f"Raw: {filename}, Converted: {py_filename}")
        print("You should re-run 'convert_for_python.m' to ensure you are plotting latest data.")
        
    return load_sim_data(py_file_path)


if __name__ == "__main__":
    # Test block
    from tkinter import filedialog, Tk
    
    print("Select a .mat file to test loader...")
    root = Tk()
    root.withdraw()
    path = filedialog.askopenfilename(filetypes=[("MAT files", "*.mat")])
    root.destroy()
    
    if path:
        # Try the smart loader if it looks like a raw file (doesn't end in _py.mat)
        if not path.endswith('_py.mat'):
            try:
                print(f"Attempting smart load of: {path}")
                results = get_sim_data(path)
                print("End-to-end load successful.")
            except Exception as e:
                print(f"End-to-end load failed: {e}")
                print("Falling back to direct load...")
                results = load_sim_data(path)
        else:
             results = load_sim_data(path)

        if results:
            print(f"\nLoaded {len(results)} signals:")
            for k, v in results.items():
                print(f" - {k}: {v.shape} (Time: {v.index.min():.2f} to {v.index.max():.2f})")
                if k == 'eta':
                    print("   Head of eta:\n", v.head())
    else:
        print("No file selected.")
