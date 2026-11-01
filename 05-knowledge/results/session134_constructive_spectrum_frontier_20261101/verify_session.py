#!/usr/bin/env python3
"""Relocate and replay every delivered lane, with normal and optimized Python."""
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,shutil,subprocess,sys,tempfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def manifest_check(root):
 entries=0;manifests=0
 for p in sorted(root.rglob('SHA256SUMS')):
  manifests+=1
  for line in p.read_text().splitlines():
   sha,name=line.split('  ',1);q=p.parent/name
   need(q.is_file(),'missing manifest target '+name)
   need(hashlib.sha256(q.read_bytes()).hexdigest()==sha,'hash mismatch '+str(q.relative_to(root)));entries+=1
 for p in sorted(root.rglob('*MANIFEST.json')):
  manifests+=1
  meta=json.loads(p.read_text());files=meta.get('files',meta.get('sha256'))
  pairs=files.items() if isinstance(files,dict) else [(row['name'],row) for row in files]
  for name,meta in pairs:
   q=p.parent/name;sha=meta if isinstance(meta,str) else meta['sha256']
   need(q.is_file(),'missing JSON manifest target '+name)
   need(hashlib.sha256(q.read_bytes()).hexdigest()==sha,'JSON hash mismatch '+name);entries+=1
 return {'manifests':manifests,'entries':entries}

JOBS=[('run_entropy','check.py',[],None,'normal.json'),
 ('intervals','verify_arc_budget.py',[],None,'normal.json'),
 ('two_apex','check_word_recovery.py',[],None,'VERIFIED.normal.json'),
 ('constructor','verify.py',[],None,'normal.json'),
 ('descent','run_checks.py',[],None,'verification.json'),
 ('radix_probe','verify.py',[],None,'normal.json'),
 ('independent','run_checks.py',[],None,None)]

def one(job,optimized):
 lane,script,args,generated,expected=job
 with tempfile.TemporaryDirectory(prefix='s134-replay-') as tmp:
  p=Path(tmp)/'packet';shutil.copytree(ROOT,p,ignore=shutil.ignore_patterns('__pycache__','.git','fresh.json'))
  env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
  cmd=[sys.executable]+(['-O'] if optimized else [])+[str(p/lane/script)]+args
  t=time.monotonic();r=subprocess.run(cmd,cwd=p/lane,env=env,capture_output=True,timeout=1800)
  need(r.returncode==0,f'{lane} failed: '+r.stderr.decode()[-3000:]+r.stdout.decode()[-3000:])
  actual=(p/lane/generated).read_bytes() if generated else r.stdout
  if expected is not None:
   want=(ROOT/lane/expected).read_bytes()
   need(actual==want,f'{lane}: replay differs from frozen receipt')
  else:need(json.loads(actual)['status']=='PASS','wrapper not PASS')
  return {'lane':lane,'script':script,'mode':'optimized' if optimized else 'normal','status':'PASS',
   'receipt_sha256':hashlib.sha256(actual).hexdigest(),'seconds':round(time.monotonic()-t,3)}

def hostile_checks():
 with tempfile.TemporaryDirectory(prefix='s134-negative-') as tmp:
  p=Path(tmp)/'packet';shutil.copytree(ROOT,p,ignore=shutil.ignore_patterns('__pycache__','.git'))
  target=p/'constructor/migration.py';old=target.read_bytes()
  for mode in ['corrupt_source','remove_source','stale_receipt']:
   if mode=='corrupt_source':target.write_bytes(old+b'\n# integrity-test mutation\n')
   elif mode=='remove_source':target.unlink()
   else:
    target.write_bytes(old);(p/'run_entropy/normal.json').write_text('{}\n')
   try:manifest_check(p)
   except RuntimeError:pass
   else:raise RuntimeError('negative control survived '+mode)
   if not target.exists() or mode!='stale_receipt':target.write_bytes(old)
 return ['corrupt_source','remove_source','stale_receipt']

def main():
 manifests=manifest_check(ROOT);negative=hostile_checks();jobs=[(j,o) for j in JOBS for o in (False,True)]
 with ThreadPoolExecutor(max_workers=2) as pool:
  futures=[pool.submit(one,j,o) for j,o in jobs];results=[]
  for f in futures:
   row=f.result();results.append(row);print(row['lane'],row['mode'],'PASS',flush=True)
 need(manifest_check(ROOT)==manifests,'manifest set changed')
 for j in JOBS:
  same=[x['receipt_sha256'] for x in results if x['lane']==j[0] and x['script']==j[1]]
  need(len(set(same))==1,'normal/optimized receipt mismatch '+j[0])
 out={'status':'PASS','manifest_verification':manifests,'negative_controls':negative,
      'isolated_principal_replays':results,'scope':'Seven principal verifiers across seven lanes, each normal and optimized, relocated independently; inner verifiers specify finite execution ranges and proof dependencies.'}
 (ROOT/'SESSION_REPLAY.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'PASS','principal_replays':len(results),**manifests,'negative_controls':len(negative)}))
if __name__=='__main__':main()
