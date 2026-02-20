
import sys
import os
import pprint

# Add current directory to path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from load_mat_data import loadmat

test_file = r"c:\Users\lenna\SynologyDrive\Ntnu\Master\Matlab_simulasjon\output data\20260111_135543\20260111_135543.mat"

print(f"Inspecting: {test_file}")
data = loadmat(test_file)

if data:
    print("\nTop Level Keys:")
    keys = list(data.keys())
    print(keys)
    
    # Check 'None' key
    if 'None' in data:
        print("\nInside 'None':")
        val = data['None']
        print(f"Type: {type(val)}")
        print(val)
        
    # Check function workspace
    if '__function_workspace__' in data:
        print("\nInside '__function_workspace__':")
        fw = data['__function_workspace__']
        print(f"Type: {type(fw)}")
        # print(fw) # Might be huge binary
        
        # Try to decode if it's a mat_struct
        if hasattr(fw, '_fieldnames'):
             print(f"Fields: {fw._fieldnames}")
    
else:
    print("Could not load data.")
