#!/usr/bin/env python3
"""Check frozen producer hashes and intentional rejection controls offline."""
from pathlib import Path
from hashlib import sha1,sha256
import json,sys,tempfile,subprocess,shutil
BASE=Path(__file__).resolve().parent.parent

def need(c,m):
 if not c:raise RuntimeError(m)

def main():
 f=BASE/'frontier';j=BASE/'jacobian_bridge';out=[]
 for line in (f/'SHA256SUMS').read_text().splitlines():
  digest,name=line.split(maxsplit=1)
  need(sha256((f/name).read_bytes()).hexdigest()==digest,'frontier payload hash '+name)
 manifest=json.loads((j/'MANIFEST.json').read_text())
 for record in manifest['files']:
  data=(j/record['path']).read_bytes()
  need((len(data),sha256(data).hexdigest())==(record['bytes'],record['sha256']),'jacobian payload hash')
 source=json.loads((j/'source_manifest.json').read_text())
 for record in source['sources']:
  data=(j/'sources'/Path(record['path']).name).read_bytes()
  need((len(data),sha256(data).hexdigest(),sha1(f'blob {len(data)}\0'.encode()+data).hexdigest())==(record['bytes'],record['sha256'],record['git_blob']),'jacobian source bytes')
 for mode in ([],['-O']):
  code=(f/'verify_frontier.py').read_text().replace('5: [9, 11, 13, 15]','5: [9, 11, 13, 15, 21]')
  need(code!=(f/'verify_frontier.py').read_text(),'mutant not injected')
  r=subprocess.run([sys.executable,*mode,'-c',code],capture_output=True,text=True)
  need(r.returncode!=0 and 'RuntimeError: small spectrum' in r.stderr,'false21 not rejected')
  out.append({'mode':'optimized' if mode else 'normal','control':'false21 spectrum','expected_failure':True})
  with tempfile.TemporaryDirectory(prefix='session142-audit-') as tmp:
   d=Path(tmp)/'frontier';shutil.copytree(f,d)
   records=json.loads((d/'repository_sources.json').read_text());records[0]['sha256']='0'*64
   (d/'repository_sources.json').write_text(json.dumps(records))
   r=subprocess.run([sys.executable,*mode,str(d/'verify_sources.py')],capture_output=True,text=True)
   need(r.returncode!=0 and 'source identity mismatch' in r.stderr,'corrupt source hash not rejected')
   out.append({'mode':'optimized' if mode else 'normal','control':'corrupt source digest','expected_failure':True})
 print(json.dumps({'status':'PASS','frontier_payloads':len((f/'SHA256SUMS').read_text().splitlines()),'jacobian_payloads':len(manifest['files']),'jacobian_source_blobs':len(source['sources']),'negative_controls':out},indent=2,sort_keys=True))
if __name__=='__main__':main()
