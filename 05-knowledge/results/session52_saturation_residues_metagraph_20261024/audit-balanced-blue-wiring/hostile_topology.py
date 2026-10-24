#!/usr/bin/env python3
import tempfile,subprocess,sys,json,shutil
from pathlib import Path
P=Path(__file__).resolve().parent;base=json.loads((P/'topology-certificate.json').read_text());results=[]
for name in ['duplicate_ear','omit_ear','degenerate_triangle','invalid_ear_endpoint']:
 with tempfile.TemporaryDirectory() as td:
  q=Path(td);shutil.copy2(P/'verify_topology.py',q/'verify_topology.py');c=json.loads(json.dumps(base))
  if name=='duplicate_ear':c['two_edge_ears'][1]=c['two_edge_ears'][0]
  elif name=='omit_ear':c['two_edge_ears'].pop()
  elif name=='degenerate_triangle':c['triangle'][1]=c['triangle'][0]
  else:c['two_edge_ears'][0][0]=-1
  (q/'topology-certificate.json').write_text(json.dumps(c))
  r=subprocess.run([sys.executable,*(['-O'] if sys.flags.optimize else []),str(q/'verify_topology.py'),str(P/'producer/normal.json')],capture_output=True,text=True)
  if not r.returncode:raise RuntimeError('accepted '+name)
  results.append({'mutation':name,'rejected':True,'failure':r.stderr.strip().splitlines()[-1]})
print(json.dumps({'status':'PASS','optimized':bool(sys.flags.optimize),'mutations':results},indent=2))
