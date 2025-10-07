import uproot
import h5py
import numpy as np
import re

root_file_path = "ORCA6_433kton_v_0_5.root"
hdf5_file_path = "ORCA6_433kton_v_0_5.h5"

# Open ROOT file
with uproot.open(root_file_path) as root_file, h5py.File(hdf5_file_path, "w") as h5_file:
    
    for key, obj in root_file.items():
        key = re.sub(r";\d+$", "", key)
        if isinstance(obj, uproot.TTree):
            group = h5_file.create_group(key)
            for branch_name, array in obj.arrays(library="np").items():
                group.create_dataset(branch_name, data=array)
    
        # Handle TAxis objects (histogram axis)
        elif isinstance(obj, uproot.models.TH.Model_TAxis_v10):
            # Extract bin edges and axis info
            axis_group = h5_file.create_group(key)
            axis_group.create_dataset("edges", data=obj.edges())
            axis_group.create_dataset("centers", data=obj.centers())

        else:
            print("Failed to convert object ", key, obj, " of type ", type(obj) )

print(f"\n✅ Converted '{root_file_path}' to '{hdf5_file_path}'")