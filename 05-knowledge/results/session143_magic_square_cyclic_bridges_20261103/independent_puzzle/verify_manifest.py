#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent
manifest=json.loads((root/'manifest.json').read_text())
for item in manifest['files']:
 data=(root/item['path']).read_bytes()
 if len(data)!=item['bytes'] or hashlib.sha256(data).hexdigest()!=item['sha256']:raise RuntimeError('hash mismatch '+item['path'])
results=[]
for flags in ([],['-O']):
 cmd=[sys.executable,*flags,str(root/'check_puzzle.py')]
 output=subprocess.check_output(cmd)
 if output!=(root/'result.json').read_bytes():raise RuntimeError('positive replay mismatch')
 results.append({'optimized':bool(flags),'case':'baseline','status':'PASS'})
 for mutation in ['single-swap','inverse-h','reciprocal-adjacency']:
  bad=subprocess.run(cmd+['--mutate',mutation],capture_output=True,text=True)
  if bad.returncode!=1 or bad.stdout or 'RuntimeError:' not in bad.stderr:raise RuntimeError('negative control failed '+mutation)
  results.append({'optimized':bool(flags),'case':mutation,'status':'REJECTED'})
print(json.dumps({'status':'PASS','checks_per_positive':json.loads(output)['checks'],'tests':results},indent=2))
