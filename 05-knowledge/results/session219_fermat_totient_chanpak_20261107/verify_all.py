#!/usr/bin/env python3
"""Portable release replay; requires Python3 and a C++17 compiler."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
R=Path(__file__).resolve().parent
commands=[]
def run(args,stdin=None):
 p=subprocess.run(args,cwd=R,input=stdin,text=True,capture_output=True)
 commands.append({'command':' '.join(map(str,args)),'returncode':p.returncode})
 if p.returncode:raise RuntimeError(p.stdout+p.stderr)
 print('PASS', ' '.join(map(str,args)),flush=True)
 return p.stdout
for f in ['totient/verify.py','ford/verify.py','fermat/verify_fermat_grammar.py','fermat/graph_witnesses.py',
          'audit/verify_dyadic_independent.py','audit/verify_supplier_independent.py','audit/verify_mersenne_independent.py',
          'audit/verify_smooth_degree_independent.py','audit/verify_phase_independent.py','run-audit/verify_witness_independent.py']:
 run([sys.executable,str(R/f)])
run([sys.executable,str(R/'run-audit/verify_graphs_independent.py'),str(R/'fermat/graph_witnesses.json')])
for f in ['ford/verify.py','fermat/verify_fermat_grammar.py']:
 normal=run([sys.executable,str(R/f)]);optimized=run([sys.executable,'-O',str(R/f)])
 if normal!=optimized:raise RuntimeError('Optimization output mismatch: '+f)
with tempfile.TemporaryDirectory(prefix='session219-') as tmp:
 for n,source in enumerate(['fermat/probe_runs.cpp','fermat/check_literal_dp.cpp','fermat/check_full_36_dp.cpp',
                            'audit/verify_run_boundary_independent.cpp','run-audit/verify_runs_independent.cpp']):
  binary=str(Path(tmp)/('check'+str(n)))
  run(['g++','-std=c++17','-O2',str(R/source),'-o',binary])
  data=(R/'fermat/literal_dp_inputs.txt').read_text() if 'check_literal_dp' in source else None
  out=run([binary],data)
  if source=='fermat/check_full_36_dp.cpp' and out!=(R/'fermat/full_36_dp.out').read_text():raise RuntimeError('Full graph output mismatch')
print(json.dumps({'status':'PASS','processes':len(commands),'limits':'Exact universes and all-length proof scopes in individual reports; no global spectrum claim.'}))
