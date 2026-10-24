#!/usr/bin/env python3
from pathlib import Path
import json,subprocess,sys
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,'-B',str(P/'verify_manifest.py')],cwd=P,check=True)
a=subprocess.check_output([sys.executable,'-B',str(P/'check.py')],cwd=P,text=True)
b=subprocess.check_output([sys.executable,'-O','-B',str(P/'check.py')],cwd=P,text=True)
if a!=b:raise RuntimeError('normal/-O mismatch')
print(a,end='');print(json.dumps({'status':'PASS','normal_optimized_equal':True},sort_keys=True))
