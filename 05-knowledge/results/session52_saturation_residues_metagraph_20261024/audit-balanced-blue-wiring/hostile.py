#!/usr/bin/env python3
import tempfile,subprocess,sys,json,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent
P=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'producer'
base=json.loads((P/'normal.json').read_text())
def mutate(d,kind):
 if kind=='missing_port':d['ports'].pop()
 elif kind=='source_path':d['ports'][0][2][0],d['ports'][0][2][1]=d['ports'][0][2][1],d['ports'][0][2][0]
 elif kind=='complement':d['ports'][0][3]^=1
 elif kind=='target_order':d['ports'][0][6][0],d['ports'][0][6][1]=d['ports'][0][6][1],d['ports'][0][6][0]
 elif kind=='edge_multiplicity':d['adjacency_multiplicities'][0][2]+=1
 elif kind=='omitted_SC_node':d['nodes'].pop()
 elif kind=='false_topology':d['articulation_nodes']=[0]
 elif kind=='target_error':d['boundary_target_classes'][0]['error']=1
out=[]
for kind in ['missing_port','source_path','complement','target_order','edge_multiplicity','omitted_SC_node','false_topology','target_error']:
 with tempfile.TemporaryDirectory() as td:
  dest=Path(td);shutil.copytree(P/'sources',dest/'sources');shutil.copy2(P/'SOURCE_PINS.json',dest/'SOURCE_PINS.json');shutil.copy2(P/'check.py',dest/'check.py');d=json.loads(json.dumps(base));mutate(d,kind);(dest/'normal.json').write_text(json.dumps(d))
  p=subprocess.run([sys.executable,*(['-O'] if sys.flags.optimize else []),str(HERE/'audit.py'),str(dest)],capture_output=True,text=True)
  if p.returncode==0:raise RuntimeError('hostile mutation accepted: '+kind)
  out.append({'mutation':kind,'rejected':True,'failure':p.stderr.strip().splitlines()[-1]})
print(json.dumps({'status':'PASS','optimized':bool(sys.flags.optimize),'mutations':out},indent=2))
