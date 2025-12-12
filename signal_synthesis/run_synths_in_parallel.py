import sys
from standalone_synth_PGSE_for_parallel_running import synth_run_bash_commands_in_parallel_with_timeout

substrates = [
    "mouse_example"
]

Nspins = 20000
tdur = 0.110 # in sec
tsteps = 2370

seq = sys.argv[1]

commands = []
for i in substrates:
    proc = ["python", f"standalone_synth_PGSE_for_parallel_running.py", f"--Nspins", f"{Nspins}", f"--Ntime", f"{tsteps}", f"--Tdur", f"{tdur}", f"--substrate", f"{i}", f"--SEQUENCE_NAME", f"{seq}"]
    commands.append(proc)

cores = int(sys.argv[2])
synth_run_bash_commands_in_parallel_with_timeout(commands, n_parallel=cores, t=2000000)