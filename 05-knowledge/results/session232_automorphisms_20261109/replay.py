#!/usr/bin/env python3
"""Replay published producers and independent checks from an isolated copy."""
import ast,concurrent.futures,hashlib,json,pathlib,shutil,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent
JOBS=[('subgroups/verify.py','subgroups/normal.out'),
 ('subgroups/verify_strata.py','subgroups/strata_normal.out'),
 ('orbit_transfer/verify.py','orbit_transfer/CHECK.json'),
 ('packing_mass/check.py','packing_mass/RESULTS.json'),
 ('boundary_graphs/verify.py','boundary_graphs/results.json'),
 ('packing_audit/independent_check.py','packing_audit/INDEPENDENT.json'),
 ('graph_audit/audit.py','graph_audit/independent.json'),
 ('subgroup_transfer_audit/independent.py','subgroup_transfer_audit/independent.normal.json')]
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(data):return hashlib.sha256(data).hexdigest()
def one(job):
 script,expected=job
 with tempfile.TemporaryDirectory(prefix='session232-replay-') as tmp:
  d=pathlib.Path(tmp)/'packet';shutil.copytree(ROOT,d,ignore=shutil.ignore_patterns('__pycache__'))
  source=d/script
  require(not any(isinstance(x,ast.Assert) for x in ast.walk(ast.parse(source.read_text()))),'removable assert: '+script)
  outs=[]
  for optimized in [False,True]:
   run=subprocess.run([sys.executable,'-B']+(['-O'] if optimized else [])+[str(source)],cwd=source.parent,capture_output=True)
   require(run.returncode==0,f'{script} failed: {run.stderr.decode()}')
   require(not run.stderr,f'{script} stderr: {run.stderr.decode()}')
   outs.append(run.stdout)
  require(outs[0]==outs[1],script+' normal/-O mismatch')
  require(outs[0]==(ROOT/expected).read_bytes(),script+' frozen receipt mismatch')
  return {'script':script,'source_sha256':digest((ROOT/script).read_bytes()),'receipt':expected,'stdout_sha256':digest(outs[0]),'normal_optimized_equal':True,'stored_receipt_equal':True}
def main():
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(one,JOBS))
 print(json.dumps({'status':'PASS','runs':16,'isolated_directory_copies':True,'jobs':records,'scope':'Replays session232 controls, not inherited THM-4812 historical certificates.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
