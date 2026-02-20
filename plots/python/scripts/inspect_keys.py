import scipy.io
import numpy as np

file_path = r'c:\Users\lenna\SynologyDrive\Ntnu\Master\Matlab_simulasjon\output data\20260116_111613\20260116_111613_py.mat'

try:
    mat = scipy.io.loadmat(file_path, squeeze_me=True, struct_as_record=False, simplify_cells=True)
    if 'data' in mat:
        data = mat['data']
        # Handle if data is a struct object or dict
        if not isinstance(data, dict):
             try:
                 names = data.dtype.names
                 data = {n: data[n] for n in names}
             except:
                 pass

        if isinstance(data, dict):
            if 'nu' in data:
                nu = data['nu']
                if hasattr(nu, 'data'): val = nu.data
                elif isinstance(nu, dict) and 'data' in nu: val = nu['data']
                else: val = nu
                
                if np.ndim(val) == 2 and val.shape[1] == 3:
                    u_surge = val[:, 0]
                    print(f"nu[:,0] (Surge): Min={np.min(u_surge)}, Max={np.max(u_surge)}")
                else:
                    print(f"nu shape {val.shape} not expected 3 columns")
            
            if 'U' in data:
                U = data['U']
                if hasattr(U, 'data'): val = U.data
                elif isinstance(U, dict) and 'data' in U: val = U['data']
                else: val = U
                print(f"U (Total Speed): Min={np.min(val)}, Max={np.max(val)}")
                     
        else:
            print(f"Data is of type {type(data)}")
    else:
        print("No 'data' key found.")

except Exception as e:
    print(f"Error: {e}")
