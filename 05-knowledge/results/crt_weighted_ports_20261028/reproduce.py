#!/usr/bin/env python3
"""Read-only exact replay of the frozen receipts, in normal and optimized mode."""
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
NAMES=['crt_probe','weighted_ports','linear_repair','energy_and_clocks','audit_alternate']

def main():
    total=0
    for name in NAMES:
        expected=json.loads((ROOT/(name+'.json')).read_text())
        for optimized in [False,True]:
            command=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(ROOT/(name+'.py'))]
            run=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,check=True)
            got=json.loads(run.stdout)
            if got!=expected:raise RuntimeError(f'{name}: frozen receipt mismatch, optimized={optimized}')
        total+=expected['active_checks']
        print(f'{name}: normal / optimized / frozen agree',flush=True)
    print(json.dumps(dict(programs=len(NAMES),checks_per_mode=total,status='PASS'),sort_keys=True))

if __name__=='__main__':main()
