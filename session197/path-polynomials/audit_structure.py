#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,sys
PACKETS=[
 ('verify_minimum_cover_polynomial.py','minimum_cover_normal.stdout','minimum_cover_optimized.stdout'),
 ('verify_endpoint_reciprocity.py','endpoint_reciprocity_normal.stdout','endpoint_reciprocity_optimized.stdout'),
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
