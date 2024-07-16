from pyntcloud import PyntCloud
import glob as gb
import pickle as pk
from scipy.spatial import ConvexHull

sub_directories = gb.glob("../playgrounds/*/")
# This should be <cells> or <feature_name> if you have something else
# and it must correspond to the folder name that you used
geom_type = "cells"
all_areas = {}
for i in sub_directories:
    files = gb.glob(i + f"{geom_type}/plys/*.ply")
    if files != []:
        name = files[0].split("/")[2]
        areas = {}  # mm
        vols = []   # mm
        for j in files:
            fname = j.split("/")[-1]
            point_cloud = PyntCloud.from_file(j)
            convex_hull_id = point_cloud.add_structure("convex_hull")
            convex_hull = point_cloud.structures[convex_hull_id]

            # Convex hull
            hull = ConvexHull(convex_hull.points[convex_hull.vertices])
            volume = hull.volume
            vols.append(volume)
            area = volume / 0.15
            areas.update({f"{fname}": area * 1e6})   # to μm^2
        all_areas.update({f"{name}": areas})

# Save the area dict for future use
hfile = open(
    f"../SIGNAL_SYNTHESIS/metrics_and_cell_data/{geom_type}_areas.bin", 'wb')
pk.dump(all_areas, hfile, protocol=5)
hfile.close()

# Compare by hand the dictionaries to make sure there is no "leakage" between
# them


def compare_dictionaries(dict1, dict2):
    lfile = open(f"{dict1}", 'rb')
    lfile2 = open(f"{dict2}", 'rb')
    areas1 = pk.load(lfile)
    areas2 = pk.load(lfile2)

    common_entries = set(areas1.keys()) & set(areas2.keys())
    common_entries_values = []

    for key in common_entries:
        if areas1[key] == areas2[key]:
            common_entries_values.append((key, areas1[key]))

    return common_entries_values
