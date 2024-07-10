import bpy
import os

sub = "mouse_example"

# Set working directory
try:
    os.chdir(f"../playgrounds/{sub}/")
except FileNotFoundError:
    return "[Errno 2] No such file or directory"

folder = "."
stl_folder = "."

name = f"{sub}_ALL_STRUCTURES.svg"

# Import svg, convert to mesh, extrude as a percent of total size, export, delete
# SVG
bpy.ops.import_curve.svg(filepath=f"{folder}/{name}")
obj = bpy.data.objects[0]
bpy.context.view_layer.objects.active = obj
selection = obj.select_get()
obj.select_set(True)

# Scale all dimensions
scale_factor = 0.778

bpy.context.object.scale[0] = scale_factor
bpy.context.object.scale[1] = scale_factor
bpy.context.object.scale[2] = scale_factor
# Extrude
# Convert to mesh
bpy.ops.object.convert(target="MESH")
bpy.ops.object.mode_set(mode='EDIT')
# Select all vertices
bpy.ops.mesh.select_all(action='SELECT')
# Get dimensions
x, y, z = bpy.data.objects[0].dimensions
z_new = 0.150
bpy.ops.mesh.extrude_region_move(
    MESH_OT_extrude_region={
        "use_normal_flip": False,
        "use_dissolve_ortho_edges": False,
        "mirror": False},
    TRANSFORM_OT_translate={"value": (0, 0, z_new),
                            "orient_axis_ortho": 'X',
                            "orient_type": 'NORMAL',
                            #    "orient_matrix":((0.895266, 0.445533, -0), (-0.445533, 0.895266, 0), (0, 0, 1)), "orient_matrix_type":'NORMAL',
                            "constraint_axis": (False, False, True),
                            "mirror": False,
                            "use_proportional_edit": False,
                            "proportional_edit_falloff": 'SMOOTH',
                            "proportional_size": 1,
                            "use_proportional_connected": False,
                            "use_proportional_projected": False,
                            "snap": False,
                            "snap_elements": {'INCREMENT'},
                            "use_snap_project": False,
                            "snap_target": 'CLOSEST',
                            "use_snap_self": False,
                            "use_snap_edit": False,
                            "use_snap_nonedit": False,
                            "use_snap_selectable": False,
                            "snap_point": (0, 0, 0),
                            "snap_align": False,
                            "snap_normal": (0, 0, 0),
                            "gpencil_strokes": False,
                            "cursor_transform": False,
                            "texture_space": False,
                            "remove_on_cancel": False,
                            "view2d_edge_pan": False,
                            "release_confirm": False,
                            "use_accurate": False,
                            "use_automerge_and_split": False})

# Select all vertices again to subdivide them
bpy.ops.mesh.select_all(action='SELECT')
# Subdivide
bpy.ops.mesh.subdivide()
# Convert faces to triangles for .ply
bpy.ops.mesh.quads_convert_to_tris(quad_method='BEAUTY', ngon_method='BEAUTY')
# Export as stl
bpy.ops.export_mesh.stl(filepath=f"{stl_folder}/{name[:-4]}.stl",
                        check_existing=True,
                        ascii=True,
                        use_selection=True)
bpy.ops.object.mode_set(mode='OBJECT')
bpy.ops.object.delete()
