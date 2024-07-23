from vol_EXTRA_func_CUSTOM import vol_w_extra
import sys
import os

"""
Example usage:
# python aggregate_all_extra.py <protocol_name> <protocol_name>
"""

substrates = os.listdir("../playgrounds/")
sequence = sys.argv[1]
protocol = sys.argv[2] # e.g PREDICT, has to match the folder name
dest_fold = f"aggregated_signals/complete_extra_signals_{protocol}"
if not os.path.exists(dest_fold):
    os.makedirs(dest_fold)
for i in substrates:
    print(f"working on {i} for {sequence}")
    vol_w_extra(i, sequence, dest_fold, protocol)
