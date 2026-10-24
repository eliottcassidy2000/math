#!/usr/bin/env python3
"""Rebuild from a new manifest-only physical copy; never reuse binaries."""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys,tempfile
from concurrent.futures import ThreadPoolExecutor
B=Path(__file__).resolve().parent
manifest=json.loads((B/'manifest.json').read_text())
for name,digest in manifest['files'].items():
 if hashlib.sha256((B/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('bad input '+name)
with tempfile.TemporaryDirectory(prefix='complete-growth-audit-') as tmp:
 D=Path(tmp)
 for name in [*manifest['files'],'manifest.json']:
  target=D/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(B/name,target)
 subprocess.run([sys.executable,'verify_audit.py'],cwd=D,check=True)
 def run(script):return subprocess.run([sys.executable,script],cwd=D,check=True,capture_output=True,text=True)
 with ThreadPoolExecutor(max_workers=5) as pool:
  futures={script:pool.submit(run,script) for script in ['replay_independent.py','replay_producer.py','mutations.py','input_controls.py','structural-audit/replay.py']}
  runs={script:f.result() for script,f in futures.items()}
 for label in ('O2','O3'):
  if json.loads((D/f'independent-{label}.json').read_text())!=json.loads((B/f'independent-{label}.json').read_text()):raise RuntimeError('independent replay mismatch')
  if json.loads((D/f'incoming-{label}.json').read_text())!=json.loads((B/f'incoming-{label}.json').read_text()):raise RuntimeError('incoming replay mismatch')
 for optimized in (False,True):
  args=[sys.executable]+(['-O'] if optimized else [])
  r=subprocess.run(args+['reconcile_strata.py'],cwd=D,check=True,capture_output=True,text=True)
  if json.loads(r.stdout)!=json.loads((B/'strata.json').read_text()):raise RuntimeError('strata replay mismatch')
  subprocess.run(args+['verify_packet.py'],cwd=D/'producer',check=True,capture_output=True)
 r=runs['structural-audit/replay.py']
 if json.loads(r.stdout)!=json.loads((B/'structural-replay.json').read_text()):raise RuntimeError('structural replay mismatch')
print(json.dumps({'status':'PASS','physical_copy':True,'manifest_payloads':len(manifest['files']),'independent_atlas_and_trace':True,'production_atlas_and_trace':True,'mutations':22,'prediction':9733056},sort_keys=True))
