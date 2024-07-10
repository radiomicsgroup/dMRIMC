import cv2
import glob as gb
import pickle as pk

sub_directories = gb.glob("../playgrounds/*/")
ext_areas = {}
for i in sub_directories:
    name = i.split("/")[2]
    img = cv2.imread(i + f'{name}_ALL_STRUCTURES.png', 0)
    _, r = cv2.threshold(img, 75, 255, cv2.THRESH_BINARY)
    n_white = cv2.countNonZero(r)

    fraction = (n_white / img.size)
    # Dimensions
    fs = open(i + "dimensions.txt")
    lines = fs.read().splitlines()
    fs.close()
    ll = lines[0].split(",")
    dim_x = float(ll[1])
    dim_y = float(ll[2])
    img_area = dim_x * dim_y * 1e6  # μm2
    area = fraction * img_area

    ext_areas.update({f"{name}": area})

# Save the area dict for future use
hfile = open(
    f"../SIGNAL_SYNTHESIS/metrics_and_cell_data/extracellular_areas.bin", 'wb')
pk.dump(ext_areas, hfile, protocol=5)
hfile.close()
