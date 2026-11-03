#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent
m=json.loads((root/'manifest.json').read_text())
for x in m['files']:
 b=(root/x['path']).read_bytes()
 if len(b)!=x['bytes'] or hashlib.sha256(b).hexdigest()!=x['sha256']:raise RuntimeError('hash mismatch '+x['path'])
expected=(root/'result.json').read_bytes();tests=[]
for flags in ([],['-O']):
 cmd=[sys.executable,*flags,str(root/'check_checkpoint1.py')]
 actual=subprocess.check_output(cmd)
 if actual!=expected:raise RuntimeError('replay mismatch')
 bad=subprocess.run(cmd+['--mutate-spectrum'],capture_output=True,text=True)
 if bad.returncode!=1 or bad.stdout or 'RuntimeError: nine distinct adjacency spectra' not in bad.stderr:raise RuntimeError('negative control not rejected')
 tests.append({'optimized':bool(flags),'positive':'PASS','negative_spectrum':'REJECTED'})
print(json.dumps({'status':'PASS','checks_per_replay':json.loads(expected)['checks'],'tests':tests},indent=2))
