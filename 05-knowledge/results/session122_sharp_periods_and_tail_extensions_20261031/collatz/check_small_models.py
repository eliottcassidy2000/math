#!/usr/bin/env python3
"""Executable scheduling and independent-language regression, E5..E7."""
import hashlib,json,pathlib,subprocess,sys,tempfile
HERE=pathlib.Path(__file__).resolve().parent

def run(cmd):return subprocess.run([str(x) for x in cmd],check=True,capture_output=True,text=True).stdout

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
rows=[]
with tempfile.TemporaryDirectory(prefix='collatz-small-') as td:
 t=pathlib.Path(td)
 flags=['-O3','-std=c++17']+(['-DNDEBUG'] if sys.flags.optimize else [])
 for src,name in [('integer_probe_intrusive.cpp','producer'),('fifo_control.cpp','fifo'),('near_tight_corridor.cpp','corridor')]:
  run(['c++',*flags,HERE/src,'-o',t/name])
 for E in (5,6,7):
  p=t/f'pot{E}.bin';f=t/f'fifo{E}.bin'
  prod=run([t/'producer',E,p,2]);fifo=run([t/'fifo',E,f,2])
  if sha(p)!=sha(f):raise ValueError(('scheduler mismatch',E))
  cpp=json.loads(run([t/'corridor',E,p,2,2,E+3]))
  pycmd=[sys.executable,'-B']+(['-O'] if sys.flags.optimize else [])
  ref=json.loads(run([*pycmd,HERE/'near_tight_corridor.py',E,p,2,2,1000000,E+3]))
  for key in ['nodes','edges','entries_checked','long_entry_controls','retained_entries','records']:
   if cpp[key]!=ref[key]:raise ValueError(('language mismatch',E,key))
  rows.append(dict(E=E,potential_sha256=sha(p),graph_nodes=cpp['nodes'],graph_edges=cpp['edges'],records=cpp['records']))
print(json.dumps(dict(status='PASS',scope='E5,E6,E7 independent scheduling and language controls; full E10 audit is separate',models=rows),sort_keys=True,indent=2))
