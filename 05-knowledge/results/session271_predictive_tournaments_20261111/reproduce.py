#!/usr/bin/env python3
"""Run scoped controls, optionally regenerate stochastic prediction tables.
No exact A000568 cache or complete inherited census is checked.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--full',action='store_true',help='also rerun the stochastic and target prediction probes')
    a=ap.parse_args()
    cmds=[['forecast/validate.py'],['structures/atlas.py'],['chanpak/controls.py'],
          ['ensemble/predict.py','--output',str(HERE/'ensemble/predictions.json')]]
    if a.full:
        cmds += [['forecast/predict.py','--output',str(HERE/'forecast/forecasts.jsonl')],
                 ['forecast/compare_methods.py'],['chanpak/predict.py'],['chanpak/target_probe.py']]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    out=[]
    for rel,*args in cmds:
        cmd=[sys.executable,'-B',str(HERE/rel),*args]
        p=subprocess.run(cmd,cwd=ROOT,env=env,capture_output=True,text=True)
        out.append({'script':rel,'returncode':p.returncode,'stdout':p.stdout[-1600:],'stderr':p.stderr[-1600:]})
        print(rel, 'PASS' if p.returncode==0 else 'FAIL',flush=True)
    result={'status':'PASS' if all(x['returncode']==0 for x in out) else 'FAIL',
            'scope':'new deliverable controls; no complete inherited audits or exact A000568 cache checks',
            'full_prediction_regeneration':a.full,'commands':out}
    (HERE/'FINAL_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
