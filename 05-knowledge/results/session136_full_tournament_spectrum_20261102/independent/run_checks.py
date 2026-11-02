#!/usr/bin/env python3
"""Replay independent checks normally/optimized; compare frozen receipts."""
import hashlib,json,os,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent
jobs=[('audit.py',''),('local_check.py','local.'),('grammar_check.py','grammar.'),('width5_check.py','width5.'),('obstruction_check.py','obstruction.'),('ear_check.py','ear.'),('tiny_extensions.py','tiny_extensions.')]
result=[]
for script,prefix in jobs:
 outputs=[]
 for optimized in (False,True):
  argv=[sys.executable,'-B']+(['-O']if optimized else[])+[str(root/script)]
  p=subprocess.run(argv,cwd=root,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},capture_output=True,check=True)
  outputs.append(p.stdout)
  wanted=root/(prefix+('optimized.json'if optimized else'normal.json'))
  if p.stdout!=wanted.read_bytes():raise RuntimeError('Frozen receipt differs: '+wanted.name)
 if outputs[0]!=outputs[1]:raise RuntimeError('Normal/optimized differ: '+script)
 result.append({'script':script,'status':'PASS','output_sha256':hashlib.sha256(outputs[0]).hexdigest()})
print(json.dumps({'status':'PASS','checks':result},indent=2,sort_keys=True))
