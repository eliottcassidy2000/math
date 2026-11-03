#!/usr/bin/env python3
"""Portable hash, positive replay, and hostile regression verification."""
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent
m=json.loads((root/'manifest.json').read_text())
for row in m['files']:
 b=(root/row['path']).read_bytes()
 if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:raise RuntimeError('hash mismatch: '+row['path'])
tests=[]
for flags in ([],['-O']):
 for script,result,mutations in [('check_independent.py','result.json',['spectrum','moment','octonion']),('check_full_cubic.py','full_cubic_result.json',['ray'])]:
  cmd=[sys.executable,*flags,str(root/script)]
  actual=subprocess.check_output(cmd)
  if actual!=(root/result).read_bytes():raise RuntimeError('positive replay mismatch '+script)
  tests.append({'script':script,'optimized':bool(flags),'positive':'PASS'})
  for mutation in mutations:
   bad=subprocess.run(cmd+['--mutate-'+mutation],capture_output=True,text=True)
   if bad.returncode!=1 or bad.stdout or 'RuntimeError:' not in bad.stderr:raise RuntimeError('mutation failed protocol '+mutation)
   tests.append({'script':script,'optimized':bool(flags),'mutation':mutation,'result':'REJECTED'})
print(json.dumps({'status':'PASS','tests':tests,'finite_check_count':json.loads((root/'result.json').read_text())['checks'],'full_cubic_systems':945},indent=2))
