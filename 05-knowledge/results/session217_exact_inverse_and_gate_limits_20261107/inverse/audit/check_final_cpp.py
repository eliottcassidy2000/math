"""Run after compiling producer inverse.cpp with UBSan into /tmp/session217_inverse_audit_ubsan."""
import runpy,subprocess,json,pathlib
D=pathlib.Path(__file__).resolve().parent
root=D.parents[1]
a=runpy.run_path(str(D/'audit.py'))
base=list(range(1,3001))
large=[4295971,9948407,1000000007,1000000000039,392,1568,6272]
raw=subprocess.check_output(['/tmp/session217_inverse_audit_ubsan',str(D.parent/'CORE_STATES.txt')]+list(map(str,base+large)),text=True)
rows=list(map(json.loads,raw.splitlines()))
for d in rows:
 n=d['N']
 if n<=3000:
  if (d['status']=='FOUND')!=bool(a['large_tab'].get(n)):raise Exception(('cpp support',n))
 if d['status']=='FOUND':
  recipe=(d['core'],d['P'],d['R'],d['b'])
  if a['val'](*recipe)!=n:raise Exception(('cpp witness',n))
  if n<=10000 and recipe not in a['large_tab'].get(n,set()):raise Exception(('cpp recipe',n))
(D/'DFS_RECEIPT.json').write_text(json.dumps({'status':'PASS','compiler':'g++ -std=c++17 -O1 -g -fsanitize=undefined -fno-sanitize-recover=all','complete_support_controls':[1,3000],'additional_targets':large,'ubsan_findings':0,'large_witnesses_verified_by':'independent Python arbitrary-precision formula'},indent=2)+'\n')
print((D/'DFS_RECEIPT.json').read_text())
