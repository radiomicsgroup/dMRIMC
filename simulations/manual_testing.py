import sys
import os
from post_simulation_tests import sim_tests

# For quick prototyping
sim_path = sys.argv[1]
sim_type = sys.argv[2]
try:
    os.chdir(sim_path)
except:
    ValueError("Simulation path unreachable or doesn't exist")

print("#####################")
print("Running the sim tests")
print("#####################")
sim_tests(sim_type)