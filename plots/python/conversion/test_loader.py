
import sys
import os

# Add current directory to path so we can import load_mat_data
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from load_mat_data import load_sim_data

# Path to a raw .mat file
test_file = r"c:\Users\lenna\SynologyDrive\Ntnu\Master\Matlab_simulasjon\output data\20260111_135543\20260111_135543.mat"

if not os.path.exists(test_file):
    print(f"Test file not found: {test_file}")
    sys.exit(1)

print(f"Testing loader on: {test_file}")

try:
    results = load_sim_data(test_file)
    
    if not results:
        print("Loader returned empty dictionary.")
    else:
        print(f"Success! Loaded {len(results)} signals.")
        for name, df in results.items():
            print(f"Signal: {name}")
            print(f"  Shape: {df.shape}")
            print(f"  Time Range: {df.index.min():.2f} to {df.index.max():.2f}")
            print(f"  Columns: {df.columns.tolist()}")
            
except Exception as e:
    print(f"Loader failed with error: {e}")
    import traceback
    traceback.print_exc()
