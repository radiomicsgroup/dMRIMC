import subprocess
import time
from collections import deque
import numpy as np

# https://stackoverflow.com/questions/69009892/python-run-multiple-subprocesses-in-parallel-but-allow-to-retry-when-a-bash-co
def run_bash_commands_in_parallel(commands, max_tries, n_parallel):
    """
    Run a list of bash commands in parallel with maximum number of processes
    """
    # Log file for errors
    f = open(f"errors.txt",'wt')

    # we use a tuple of (command, tries) to store the information
    # how often a command was already tried.
    waiting = deque([(command, 1) for command in commands])
    running = deque()
    while len(waiting) > 0 or len(running) > 0:
        print(f"Running: {len(running)}, Waiting: {len(waiting)}",end='\r')
        # if less than n_parallel jobs are running and we have waiting jobs,
        # start new jobs
        while len(waiting) > 0 and len(running) < n_parallel:
            command, tries = waiting.popleft()
            try:
                running.append((subprocess.Popen(command,stdout=subprocess.DEVNULL), command, tries))
                print(f"Started working on {command[1].split('/')[1]}")
            except OSError:
                print(f"{command[1].split('/')[1]} simulation failed")
        # poll running commands
        for _ in range(len(running)):
            process, command, tries = running.popleft()
            ret = process.poll()
            if ret is None:
                running.append((process, command, tries))
            # retry errored jobs
            elif ret != 0:
                if tries < max_tries:
                    waiting.append((command, tries  + 1))
                else:
                    print(f"Simulation of {command[1].split('/')[1]} errored after {max_tries} tries")
                    f.write(f"{command[1]}\n")
            else:
                print(f"Simulation of {command[1].split('/')[1]} finished successfully")
        # sleep a bit to reduce CPU usage
        time.sleep(0.5)
    print('All tasks done')
    f.close()

def run_bash_commands_in_parallel_with_timeout_with_retry(commands, max_tries, n_parallel, t):
    """
    Run a list of bash commands in parallel with maximum number of processes
    and kill the process if it takes more than {t} minutes. If the
    process fails it will be retried {max_tries} times
    """
    # Log files for errors and timeouts
    error_file = open("errors.txt", 'wt')
    timeout_file = open("timed_out.txt", 'wt')
    waiting = deque([(command, 1) for command in commands])
    running = deque()
    while len(waiting) > 0 or len(running) > 0:
        print(f"Running: {len(running)}, Waiting: {len(waiting)}", end='\r')

        while len(waiting) > 0 and len(running) < n_parallel:
            command, tries = waiting.popleft()
            try:
                process = subprocess.Popen(command, stdout=subprocess.DEVNULL)
                running.append((process, command, tries, time.time()))
                print(f"Started working on {command[1].split('/')[1]}")
            except OSError:
                print(f"{command[1].split('/')[1]} simulation failed to start")
                error_file.write(f"{command[1]}\n")

        for _ in range(len(running)):
            process, command, tries, start_time = running.popleft()
            ret = process.poll()
            if ret is None:
                if time.time() - start_time > t:
                    process.kill()
                    timeout_file.write(f"{command}\n")
                    print(f"{command[1].split('/')[1]} timed out")
                else:
                    running.append((process, command, tries, start_time))
            elif ret != 0:
                if tries < max_tries:
                    waiting.append((command, tries + 1))
                else:
                    print(f"Simulation of {command[1].split('/')[1]} errored after {max_tries} tries")
                    error_file.write(f"{command[1]}\n")
            else:
                print(f"Simulation of {command[1].split('/')[1]} finished successfully")

        time.sleep(0.5)

    print('All tasks done')
    error_file.close()
    timeout_file.close()


def run_bash_commands_in_parallel_with_timeout(commands, n_parallel, t, retry=False):
    """
    Run a list of bash commands in parallel with maximum number of processes
    and kill the process if it takes more than t seconds.
    """
    # Log files for errors and timeouts
    if retry:
        error_file = open("errors_retry.txt", 'wt')
        timeout_file = open("timed_out_retry.txt", 'wt')
    else:
        error_file = open("errors.txt", 'wt')
        timeout_file = open("timed_out.txt", 'wt')

    waiting = deque([(command, 1) for command in commands])
    running = deque()
    cnt_errs = 0
    cnt_timeouts = 0
    cnt_completed = 0
    while len(waiting) > 0 or len(running) > 0:
        print(f"Running: {len(running)}, Waiting: {len(waiting)}, Completed: {cnt_completed}, Errors: {cnt_errs}, Timeouts: {cnt_timeouts}", end='\r')
        while len(waiting) > 0 and len(running) < n_parallel:
            command, tries = waiting.popleft()
            try:
                process = subprocess.Popen(command, stdout=subprocess.DEVNULL)
                running.append((process, command, tries, time.time()))
                print(f"Started working on {command[1].split('/')[1]}")
            except OSError:
                print(f"{command[1].split('/')[1]} simulation failed to start")
                error_file.write(f"{command[1]}\n")
                cnt_errs += 1

        for _ in range(len(running)):
            process, command, tries, start_time = running.popleft()
            ret = process.poll()
            if ret is None:
                if time.time() - start_time > t:
                    process.kill()
                    timeout_file.write(f"{command[1]}\n")
                    print(f"{command[1].split('/')[1]} timed out")
                    cnt_timeouts += 1
                else:
                    running.append((process, command, tries, start_time))
            elif ret != 0:
                print(f"Simulation of {command[1].split('/')[1]} returned an error")
                error_file.write(f"{command[1]}\n")
                cnt_errs += 1
            else:
                exec_time = np.round((time.time() - start_time)/60, decimals=4)
                print(f"Simulation of {command[1].split('/')[1]} finished successfully after {exec_time} minutes")
                cnt_completed += 1
        time.sleep(0.5)

    print('All tasks done')
    error_file.close()
    timeout_file.close()