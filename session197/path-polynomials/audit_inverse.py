#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,sys
PACKETS=[
 ('verify_ratio_lattice.py','ratio_lattice_normal.stdout','ratio_lattice_optimized.stdout'),
 ('verify_target_ray_inverse.py','target_ray_normal.stdout','target_ray_optimized.stdout'),
 ('verify_two_generator_classification.py','two_generator_normal.stdout','two_generator_optimized.stdout'),
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
