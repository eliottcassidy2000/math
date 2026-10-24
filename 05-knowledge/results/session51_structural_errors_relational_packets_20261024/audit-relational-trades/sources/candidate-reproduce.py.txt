#!/usr/bin/env python3
from pathlib import Path
import json,subprocess,sys
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,'-B',str(P/'verify_manifest.py')],cwd=P,check=True)
scripts=['check.py','structural.py','inherited/receiver-thresholds/check.py','inherited/receiver-thresholds/inherited/native-repairs/check.py','inherited/receiver-thresholds/inherited/native-repairs/inherited/verify_native.py','inherited/receiver-thresholds/inherited/native-repairs/inherited/replay_filtration.py']
for script in scripts:
 outputs=[]
 for optimized in (False,True):
  command=[sys.executable]+(['-O'] if optimized else [])+['-B',str(P/script)]
  out=subprocess.check_output(command,cwd=P,text=True);outputs.append(out);print(out,end='',flush=True)
 if outputs[0]!=outputs[1]:raise RuntimeError('normal/-O mismatch: '+script)
print(json.dumps({'status':'PASS','normal_optimized_equal':True,'scripts':len(scripts)},sort_keys=True))
