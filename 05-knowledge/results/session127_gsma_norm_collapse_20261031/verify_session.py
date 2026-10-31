#!/usr/bin/env python3
"""Verify hashes and replay exact Session127 checks after relocation.
Requires Python3 and SymPy. No network or repository checkout is needed.
"""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent
JOBS=[('tribonacci/check.py','tribonacci/receipt.json',None),
 ('gsma/check_unit_kernel.py','gsma/unit_kernel.normal.json',None),
 ('lucas47/verify.py','lucas47/normal.json',None),
 ('lucas47/independent.py','lucas47/independent.normal.json',None),
 ('uniform_norm/verify_uniform_norm.py','uniform_norm/verify_uniform_norm.normal.txt','uniform_norm/verify_uniform_norm.output.json'),
 ('volume/verify_volume_bridge.py',None,'volume/verification.json')]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifests(root):
 count=0;files=0
 for m in sorted(root.glob('*/SHA256SUMS')):
  count+=1
  for line in m.read_text().splitlines():
   expected,name=line.split('  ',1);p=m.parent/name
   if digest(p)!=expected:raise RuntimeError('hash mismatch '+str(p.relative_to(root)))
   files+=1
 return count,files

def main():
 counts=manifests(ROOT);runs=[]
 with tempfile.TemporaryDirectory(prefix='session127-replay-') as tmp:
  dest=Path(tmp)/'packet';shutil.copytree(ROOT,dest,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
  for script,expected_stdout,expected_file in JOBS:
   outputs=[]
   for optimized in [False,True]:
    args=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(dest/script)]
    p=subprocess.run(args,cwd=tmp,capture_output=True,check=True)
    outputs.append(p.stdout)
    if expected_stdout and p.stdout!=(ROOT/expected_stdout).read_bytes():raise RuntimeError('stdout mismatch '+script)
    if expected_file and digest(dest/expected_file)!=digest(ROOT/expected_file):raise RuntimeError('generated mismatch '+script)
   if outputs[0]!=outputs[1]:raise RuntimeError('normal/optimized mismatch '+script)
   runs.append({'script':script,'normal':'PASS','optimized':'PASS','relocated':True,'stdout_equal':True,'generated_bytes_verified':expected_file})
  manifests(dest)
 report={'status':'PASS','manifest_count':counts[0],'manifest_entries':counts[1],'runs':runs,'dependencies':'Python3 and SymPy for the volume checker','scope':'Exact packet checks and finite controls. Cited global theorems remain mathematical dependencies; no Lean build or global quartic/quintic classification is claimed.'}
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
