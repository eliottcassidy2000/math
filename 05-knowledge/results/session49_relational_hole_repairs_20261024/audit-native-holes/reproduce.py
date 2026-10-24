#!/usr/bin/env python3
from pathlib import Path
import json,subprocess,sys
P=Path(__file__).resolve().parent
results={}
for script in ('audit.py','replay_inherited_filtration.py','compare_producer.py'):
    normal=subprocess.run([sys.executable,str(P/script)],cwd=P,check=True,capture_output=True,text=True)
    optimized=subprocess.run([sys.executable,'-O',str(P/script)],cwd=P,check=True,capture_output=True,text=True)
    if normal.stdout!=optimized.stdout:raise RuntimeError('normal/-O outputs differ for '+script)
    results[script]=json.loads(normal.stdout)
print(json.dumps({'status':'PASS','normal_optimized_equal':True,'results':results},sort_keys=True))
