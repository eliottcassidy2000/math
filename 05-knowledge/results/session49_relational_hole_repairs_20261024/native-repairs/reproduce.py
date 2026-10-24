#!/usr/bin/env python3
from pathlib import Path
import json,subprocess,sys
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,'-B',str(P/'verify_manifest.py')],cwd=P,check=True)
for script in ('check.py','inherited/verify_native.py','inherited/replay_filtration.py'):
 outputs=[]
 for optimized in (False,True):
  cmd=[sys.executable]+(['-O'] if optimized else [])+['-B',str(P/script)]
  out=subprocess.check_output(cmd,cwd=P,text=True);print(out,end='',flush=True);outputs.append(out)
 if outputs[0]!=outputs[1]:raise RuntimeError('Normal/-O mismatch: '+script)
print(json.dumps({'status':'PASS','normal_optimized_equal':True,'scripts':3},sort_keys=True))
