#!/usr/bin/env python3
"""Verify hashes and replay exact Session128 checks after relocation.
Requires Python3 and SymPy. No network or repository checkout is needed.
"""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent
JOBS=[('paper_audit/verify_persistent_bias.py','paper_audit/persistent_bias.json',None),
 ('feedback/verify_feedback.py','feedback/primary.normal.json',None),
 ('feedback/independent_feedback.py','feedback/independent.normal.json',None),
 ('triadic/check_triadic.py','triadic/normal.json',None),
 ('native_audit/check_native_controls.py','native_audit/normal.json',None),
 ('pisot/verify.py','pisot/normal.json',None),
 ('resolvent_norm/check_resolvent.py','resolvent_norm/normal.json',None),
 ('nonunit_feedback/verify.py','nonunit_feedback/normal.json',None),
 ('spectral_audit/verify_independent.py','spectral_audit/normal.json',None)]

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
 with tempfile.TemporaryDirectory(prefix='session128-replay-') as tmp:
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
 report={'status':'PASS','manifest_count':counts[0],'manifest_entries':counts[1],'runs':runs,'dependencies':'Python3 and SymPy for the primary native-feedback checker','scope':'Exact packet checks and finite controls. Cited global theorems remain mathematical dependencies; no Lean build, general nonlinear convergence, tiling or pure-spectrum theorem is claimed.'}
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
