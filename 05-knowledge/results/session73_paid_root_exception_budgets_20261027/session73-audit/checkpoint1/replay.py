from pathlib import Path
import subprocess,json,concurrent.futures,hashlib
BASE=Path(__file__).parent
jobs=[('ports','verify.py','receipt.json'),('moduli','check.py','receipt.json'),('arithmetic','check.py','normal.json'),('root','check_compiler_budget.py','compiler.normal.json'),('forcing','check_core.py','core.normal.json')]
def run(job):
 lane,script,native=job;p=BASE/'reviewed'/lane
 rows=[]
 for mode in [[],['-O']]:
  tag='optimized' if mode else 'normal';r=subprocess.run(['python3','-B',*mode,script],cwd=p,capture_output=True)
  (p/('replay-'+tag+'.json')).write_bytes(r.stdout)
  if r.returncode:raise RuntimeError(lane+' '+tag+' '+r.stderr.decode())
  if r.stdout!=(p/native).read_bytes():raise RuntimeError(lane+' receipt changed')
  rows.append({'mode':tag,'exit':r.returncode,'matches_frozen_receipt':True,'bytes':len(r.stdout),'sha256':hashlib.sha256(r.stdout).hexdigest()})
 return {'lane':lane,'replays':rows}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(run,jobs))
controls=[];p=BASE/'reviewed'/'arithmetic'
for mode in [[],['-O']]:
 for flag in ['target-only','drop-guard','omit-cube','count-is-weight']:
  r=subprocess.run(['python3','-B',*mode,'check.py','--negative-control',flag],cwd=p,capture_output=True,text=True)
  if r.returncode==0:raise RuntimeError('hostile passed '+flag)
  controls.append({'mode':'optimized' if mode else 'normal','control':flag,'exit':r.returncode,'stderr_last_line':r.stderr.strip().splitlines()[-1]})
# Exact source snapshots used for this review must still match their producers.
for x in json.loads((BASE/'REVIEWED_INPUTS.json').read_text()):
 b=Path(x['producer_path']).read_bytes()
 if len(b)!=x['bytes'] or hashlib.sha256(b).hexdigest()!=x['sha256']:raise RuntimeError('producer changed '+x['producer_path'])
result={'status':'PASS','producer_replays':rows,'negative_controls':controls,'producer_bytes_unchanged_at_review':True}
(BASE/'REPLAY.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
