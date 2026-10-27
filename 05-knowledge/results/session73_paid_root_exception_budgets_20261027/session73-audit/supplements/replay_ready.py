from pathlib import Path
import subprocess,json,hashlib,concurrent.futures
BASE=Path(__file__).parent
jobs=[('arithmetic','check.py','normal.json'),('scaling','verify_scaling.py','receipt.json'),('root-scalar','check_scalar_budget.py','scalar.normal.json')]
def run(job):
 lane,script,native=job;p=BASE/'reviewed'/lane;out=[]
 for mode in [[],['-O']]:
  tag='optimized' if mode else 'normal';r=subprocess.run(['python3','-B',*mode,script],cwd=p,capture_output=True)
  (p/('replay-'+tag+'.json')).write_bytes(r.stdout)
  if r.returncode:raise RuntimeError(lane+' failed '+r.stderr.decode())
  if r.stdout!=(p/native).read_bytes():raise RuntimeError(lane+' receipt mismatch')
  out.append({'mode':tag,'exit':0,'bytes':len(r.stdout),'sha256':hashlib.sha256(r.stdout).hexdigest(),'equal_native':True})
 return {'lane':lane,'runs':out}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,jobs))
controls=[]
for mode in [[],['-O']]:
 r=subprocess.run(['python3','-B',*mode,'check.py','--negative-control'],cwd=BASE/'reviewed/arithmetic',capture_output=True,text=True)
 if not r.returncode:raise RuntimeError('equal prime hostile passed')
 controls.append({'mode':mode,'exit':r.returncode,'last_error':r.stderr.strip().splitlines()[-1]})
result={'status':'PASS','replays':rows,'negative_controls':controls}
(BASE/'READY_REPLAY.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
