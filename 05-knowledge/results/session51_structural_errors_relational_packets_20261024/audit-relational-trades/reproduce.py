#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys, json
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,'-B',str(P/'verify_manifest.py')],cwd=P,check=True)
scripts=['audit.py','structural_audit.py','native/verify_native.py','native/replay_filtration.py']
for script in scripts:
    outputs=[]
    for optimized in (False,True):
        command=[sys.executable]+(['-O'] if optimized else [])+['-B',str(P/script)]
        out=subprocess.check_output(command,cwd=P,text=True)
        print(out,end='',flush=True);outputs.append(out)
    if outputs[0]!=outputs[1]:raise RuntimeError('normal/optimized mismatch: '+script)
print(json.dumps({'status':'PASS','normal_optimized_equal':True,'scripts':len(scripts)},sort_keys=True))
