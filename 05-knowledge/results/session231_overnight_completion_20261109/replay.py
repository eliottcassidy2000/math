"""Replay compact exact receipts in a temporary copy. No tracked files are changed.
Run: python3 -I -B replay.py [--full-physical]
The optional final check uses about 3.7 GB RAM for 25-vertex subset tables.
"""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
REPO=next(p for p in HERE.parents if (p/'.git').exists())
CASES=[
 ('descent/verify.py','descent/CHECK.json'),
 ('descent_audit/check.py','descent_audit/CHECK.json'),
 ('descent/round2_verify.py','descent/ROUND2_CHECK.json'),
 ('descent_audit/round2_check.py','descent_audit/ROUND2_CHECK.json'),
 ('fibers/verify.py','fibers/CHECK.json'),
 ('fibers/verify_repairs.py','fibers/REPAIRS.json'),
 ('fibers/verify_family.py','fibers/FAMILY_CHECK.json'),
 ('fibers/verify_unit_family.py','fibers/UNIT_FAMILY_CHECK.json'),
 ('fibers/verify_short12_unit.py','fibers/SHORT12_CHECK.json'),
 ('fibers/verify_short12_independent.py','fibers/SHORT12_INDEPENDENT.json'),
 ('fibers/baseline_and_gcd.py','fibers/BASELINE_GCD.json'),
 ('fibers/verify_contexts.py','fibers/CONTEXT_CHECK.json'),
 ('fibers/check_shorter_loops.py','fibers/SHORTER_LOOP_CHECK.json'),
 ('fiber_audit/unit_audit.py','fiber_audit/UNIT_AUDIT.json'),
 ('fiber_audit/short12_audit.py','fiber_audit/SHORT12_AUDIT.json'),
 ('fiber_audit/minimal_walks.py','fiber_audit/MINIMAL_WALKS.json'),
 ('fiber_audit/placement_audit.py','fiber_audit/PLACEMENT_AUDIT.json'),
 ('fiber_audit/diameter_identity.py','fiber_audit/DIAMETER_IDENTITY.json'),
 ('source/verify.py','source/CHECK.json'),
 ('source_audit/independent.py','source_audit/RECEIPT.json'),
 ('source/ROUND2/verify.py','source/ROUND2/CHECK.json'),
 ('source_audit/round2_independent.py','source_audit/ROUND2_RECEIPT.json'),
 ('dyadic/verify.py','dyadic/CHECK.json'),
 ('dyadic_audit/check.py','dyadic_audit/CHECK.json'),
 ('ribbon/verify.py','ribbon/VERIFICATION.json'),
 ('ribbon/parity_probe.py','ribbon/PARITY_PROBE.json'),
 ('ribbon_audit/check.py','ribbon_audit/RESULTS.jsonl'),
 ('ribbon_audit/check_kernel.py','ribbon_audit/KERNEL.json'),
 ('ribbon_audit/check_active_gates.py','ribbon_audit/GATES.json'),
 ('placement/verify.py','placement/CHECK.json'),
 ('placement/check_selector.py','placement/SELECTOR_CHECK.json'),
 ('placement_audit/audit.py','placement_audit/CHECK.json'),
 ('rerooted/verify.py','rerooted/VERIFY.json'),
]
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def same(a,b):
 if a==b:return True
 try:return json.loads(a)==json.loads(b)
 except (ValueError,UnicodeDecodeError):return False
rows=[]
with tempfile.TemporaryDirectory(prefix='session231-replay-') as tmp:
 work=Path(tmp)/'packet';shutil.copytree(HERE,work,ignore=shutil.ignore_patterns('__pycache__','*.bin'))
 env=dict(os.environ,MATH_REPO=str(REPO));subprocess.run(['g++','-std=c++17','-O2',str(work/'source/inverse.cpp'),'-o',str(work/'source/inverse')],check=True)
 cases=list(CASES)
 if '--full-physical' in sys.argv:
  subprocess.run(['g++','-std=c++17','-O3',str(work/'rerooted_audit/physical.cpp'),'-o',str(work/'rerooted_audit/physical')],check=True)
  cases.append(('rerooted_audit/check.py','rerooted_audit/CHECK.json'))
 for script,receipt in cases:
  outputs=[]
  for optimized in (False,True):
   args=[sys.executable,'-I','-B']+(['-O'] if optimized else [])+[str(work/script)]
   result=subprocess.run(args,cwd=(work/script).parent,env=env,check=True,capture_output=True)
   require(same(result.stdout,(HERE/receipt).read_bytes()),'receipt mismatch: '+script)
   outputs.append(result.stdout)
  require(outputs[0]==outputs[1],'normal/optimized mismatch: '+script)
  rows.append(dict(script=script,receipt=receipt,status='PASS',stdout_sha256=hashlib.sha256(outputs[0]).hexdigest()))
  print('PASS '+script,file=sys.stderr,flush=True)
print(json.dumps(dict(status='PASS',cases=len(rows),normal_and_optimized=True,full_physical='--full-physical' in sys.argv,scope='Exact compact receipts; large census regeneration is documented separately.',checks=rows),indent=2))
