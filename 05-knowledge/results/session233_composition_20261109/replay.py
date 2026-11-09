"""Replay the assembled mathematical packet in normal and optimized Python."""
from pathlib import Path
import subprocess,sys,json,hashlib
P=Path(__file__).resolve().parent
TASKS=[('graft','verify.py','--fail'),('graft','verify_weighted.py','--fail'),
 ('graft_audit','audit.py','--sentinel'),('graft_audit','audit_supplement.py','--sentinel'),
 ('endpoint','verify_endpoint.py','--fail'),('endpoint','response.py','--fail'),
 ('endpoint_audit','audit.py','--fail'),('transducer','bridge_blocks.py','--fail'),
 ('transducer','check_polynomial.py','--fail'),('transducer_audit','audit.py','--failure-sentinel'),
 ('.','check_legacy_flip.py','--fail')]
rows=[]
for folder,script,sentinel in TASKS:
 outputs=[]
 for optimize in [False,True]:
  cmd=[sys.executable,'-B']+(['-O'] if optimize else [])+[script]
  r=subprocess.run(cmd,cwd=P/folder,capture_output=True)
  if r.returncode:raise RuntimeError((folder,script,r.stderr.decode()))
  outputs.append(r.stdout)
  fail=subprocess.run(cmd+[sentinel],cwd=P/folder,capture_output=True)
  if fail.returncode==0:raise RuntimeError(('inactive failure control',folder,script,sentinel))
 if outputs[0]!=outputs[1]:raise RuntimeError(('optimization mismatch',folder,script))
 rows.append({'folder':folder,'script':script,'normal':'PASS','optimized':'PASS','stdout_sha256':hashlib.sha256(outputs[0]).hexdigest(),'failure_control':sentinel,'failure_active_both_modes':True})
 print('PASS',folder,script,flush=True)
result={'status':'PASS','normal_and_optimized_runs':2*len(rows),'failure_control_runs':2*len(rows),'tests':rows,'scope':'Compact certificate and audit replay, not rediscovery or exhaustive enumeration of large constructed graphs.'}
(P/'REPLAY.json').write_text(json.dumps(result,indent=2)+'\n')
