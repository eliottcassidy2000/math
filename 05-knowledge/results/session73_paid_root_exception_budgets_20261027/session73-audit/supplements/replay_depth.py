from pathlib import Path
import subprocess,hashlib,json,concurrent.futures
B=Path(__file__).parent;expected=(B/'reviewed/depth/receipt.json').read_bytes()
def run(opt):
 name='optimized' if opt else 'normal';p=B/('depth-exec-'+name)
 r=subprocess.run(['python3','-B',*(['-O'] if opt else []),'search_depth.py'],cwd=p,capture_output=True)
 (p/'stdout.json').write_bytes(r.stdout)
 if r.returncode:raise RuntimeError(name+' failed '+r.stderr.decode())
 result=(p/'search_result.json').read_bytes()
 if result!=expected or r.stdout!=expected:raise RuntimeError(name+' differs')
 return {'mode':name,'exit':0,'bytes':len(result),'sha256':hashlib.sha256(result).hexdigest(),'equals_native':True}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:out=list(pool.map(run,[False,True]))
p=B/'reviewed/depth';ver=[]
for mode in [[],['-O']]:
 r=subprocess.run(['python3','-B',*mode,'verify_depth_receipt.py'],cwd=p,capture_output=True)
 if r.returncode:raise RuntimeError(r.stderr.decode())
 (p/('replay-verify-'+('optimized' if mode else 'normal')+'.json')).write_bytes(r.stdout)
 ver.append({'mode':mode,'exit':0,'bytes':len(r.stdout),'sha256':hashlib.sha256(r.stdout).hexdigest()})
result={'status':'PASS','full_search_runs':out,'independent_receipt_verifier':ver}
(B/'DEPTH_REPLAY.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
