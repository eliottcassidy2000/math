#!/usr/bin/env python3
"""Portable replay plus deliberately damaged variants. Does not mutate sources."""
import hashlib,json,pathlib,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent
src=(ROOT/'check_figurate.py').read_text()
mutations=[
 ('wrong_Cassini_cross_term','return x*x+x*y-y*y','return x*x-x*y-y*y'),
 ('wrong_triangular_index','return n*(n+1)//2','return n*(n-1)//2'),
 ('counter_presence_instead_of_parity','return sum((e%2)<<i for i,e in enumerate(v))','return sum(int(e>0)<<i for i,e in enumerate(v))'),
 ('replace_15_over2_by_15_over14','(15,2),(1,7),(55,1)','(15,14),(1,7),(55,1)'),
 ('omit_converse_quotient','min(canonical[m],canonical[full^m])','canonical[m]'),
 ('break_recurrence_window','return 1,n,n+1,2*n+1','return 1,n,n+1,2*n+2'),
]
receipt={'status':'PASS','replays':[],'mutations':[]}
with tempfile.TemporaryDirectory(prefix='s129-figurate-replay-') as directory:
 p=pathlib.Path(directory)/'relocated.py';p.write_text(src)
 for opts in ([],['-O']):
  r=subprocess.run([sys.executable,*opts,str(p)],cwd=directory,capture_output=True,text=True)
  if r.returncode:raise RuntimeError(r.stderr)
  data=json.loads(r.stdout)
  if r.stdout!=(ROOT/'normal.json').read_text():raise RuntimeError('relocated output changed')
  receipt['replays'].append({'optimized':bool(opts),'status':data['status'],'explicit_checks':data['explicit_checks'],'stdout_sha256':hashlib.sha256(r.stdout.encode()).hexdigest()})
 for name,old,new in mutations:
  if src.count(old)!=1:raise RuntimeError('mutation target not unique: '+name)
  p.write_text(src.replace(old,new));r=subprocess.run([sys.executable,'-O',str(p)],cwd=directory,capture_output=True,text=True)
  if r.returncode==0:raise RuntimeError('mutation survived: '+name)
  receipt['mutations'].append({'name':name,'status':'REJECTED','exit_code':r.returncode,'last_error':r.stderr.strip().splitlines()[-1]})
print(json.dumps(receipt,sort_keys=True,indent=2))
