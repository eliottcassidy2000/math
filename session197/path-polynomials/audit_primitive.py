#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,sys
PACKETS=[
 ('verify_primitive_cancellation.py','primitive_cancellation_normal.stdout','primitive_cancellation_optimized.stdout'),
 ('verify_primitive_target_inverse.py','primitive_target_normal.stdout','primitive_target_optimized.stdout'),
 ('verify_primitive_leaf_rigidity.py','primitive_leaf_normal.stdout','primitive_leaf_optimized.stdout'),
]
def run():
 here=Path(__file__).resolve().parent;report=[]
 for script,norm,opt in PACKETS:
  a=subprocess.run([sys.executable,str(here/script)],cwd=here,check=True,capture_output=True).stdout
  b=subprocess.run([sys.executable,'-O',str(here/script)],cwd=here,check=True,capture_output=True).stdout
  if not (a==b==(here/norm).read_bytes()==(here/opt).read_bytes()):raise RuntimeError('Receipt mismatch: '+script)
  report.append({'script':script,'status':'PASS','receipt_sha256':hashlib.sha256(a).hexdigest()})
 print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__':run()
