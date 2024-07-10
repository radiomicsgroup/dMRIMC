import bpy
import numpy as np
import os

sub = "mouse_example"

# Set working directory
try:
    os.chdir(f"../playgrounds/{sub}/")
except FileNotFoundError:
    return "[Errno 2] No such file or directory"

stl_folder = "./cells/stls"
ply_folder = "./cells/plys"
if not os.path.exists(ply_folder):
    os.makedirs(ply_folder)
files = os.listdir(stl_folder)

# Check if the dimensions.txt exists, if yes delete it
# This is needed because we are appending to the opened txt ('at'),
# so when running the script again the dimensions.txt is not recreated
# we are just adding lines at end
if os.path.exists(f"{ply_folder}/dimensions.txt"):
    os.remove(f"{ply_folder}/dimensions.txt")

for i in files:
    name = i
    bpy.ops.import_mesh.stl(filepath=f"{stl_folder}/{i}")
    # Select object
    obj = bpy.data.objects[0]
    bpy.context.view_layer.objects.active = obj
    selection = obj.select_get()
    obj.select_set(True)
    bpy.ops.export_mesh.ply(filepath=f"{ply_folder}/{name[:-4]}.ply",
                            check_existing=True,
                            use_ascii=True,
                            use_selection=True,
                            use_colors=False,
                            use_normals=False)
    x, y, z = obj.dimensions
    z_new = 0.150
    vol = x * y * z_new
    vol2 = np.round(vol * 10**9, decimals=6)  # μm^3
    f = open(f"{ply_folder}/dimensions.txt", 'at')
    f.write(f"{name[:-4]},{x},{y},{z},{vol2}\n")
    f.close()
    bpy.ops.object.delete()
