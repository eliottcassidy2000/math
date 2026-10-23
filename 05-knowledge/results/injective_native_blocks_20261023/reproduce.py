"""Run every frozen checker normally and with -O; require identical outputs."""
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
jobs = [('electrical/check.py',[]),('electrical/check_larger_image.py',[]),
        ('search/check.py',[]),('search/larger_check.py',[]),
        ('search/larger_check.py',['--sixteen']),('assembly/check.py',[]),
        ('dyadic/check.py',[]),('audit.py',[]),('audit_larger.py',[])]
for script,extra in jobs:
    outputs = []
    for flags in (['-B'], ['-O','-B']):
        run = subprocess.run([sys.executable,*flags,str(HERE/script),*extra],capture_output=True,text=True)
        if run.returncode:
            raise RuntimeError(f'{script} {flags}:\n{run.stdout}\n{run.stderr}')
        outputs.append(run.stdout)
    if outputs[0] != outputs[1]:
        raise RuntimeError(f'Optimization changed output: {script}')
    print(f'PASS normal/-O: {script} {" ".join(extra)}',flush=True)
print('All frozen receipts reproduced; normal and optimized outputs agree.')
