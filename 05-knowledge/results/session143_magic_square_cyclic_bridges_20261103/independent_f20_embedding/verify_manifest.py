#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent
m=json.loads((root/'manifest.json').read_text())
for row in m['files']:
 b=(root/row['path']).read_bytes()
 if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:raise RuntimeError('hash mismatch '+row['path'])
tests=[]
for flags in ([],['-O']):
 cmd=[sys.executable,*flags,str(root/'check_f20.py')]
 output=subprocess.check_output(cmd)
 if output!=(root/'result.json').read_bytes():raise RuntimeError('baseline mismatch')
 tests.append({'optimized':bool(flags),'case':'baseline','result':'PASS'})
 for mutation in ['centralizer','character','permutation-copies']:
  z=subprocess.run(cmd+['--mutate',mutation],capture_output=True,text=True)
  if z.returncode!=1 or z.stdout or 'RuntimeError:' not in z.stderr:raise RuntimeError('mutation failed '+mutation)
  tests.append({'optimized':bool(flags),'case':mutation,'result':'REJECTED'})
print(json.dumps({'status':'PASS','checks_per_baseline':json.loads(output)['checks'],'tests':tests},indent=2))
