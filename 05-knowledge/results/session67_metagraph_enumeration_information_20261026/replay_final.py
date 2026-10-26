#!/usr/bin/env python3
"""Replay frozen final lanes in temporary copies; standard library, no installs."""
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parent
jobs=[
 ('session67-information','verify_information.py',[],'verification.out'),
 ('session67-information','verify_core_information.py',[],'core_information_verification.out'),
 ('session67-information','verify_growth_information.py',[],'growth_information_verification.out'),
 ('session67-metagraph','verify_core_families.py',['--max-order','40','--class-order','7'],'receipt.json'),
 ('session67-recursion','verify_recursion.py',['--max-order','8'],'receipt-order8.json'),
 ('session67-recursion','audit_hostile.py',[],'hostile-receipt.json'),
 ('session67-cycle-audit','verify.py',[],'receipt.json'),
 ('session67-cycle-audit','audit_core_triangle_portable.py',['--metagraph-dir','../session67-metagraph'],'CORE_TRIANGLE_PORTABLE_RECEIPT.json'),
 ('session67-cycle-audit','audit_recurrence_cross_portable.py',['--recursion-dir','../session67-recursion'],'RECURRENCE_CROSS_PORTABLE_RECEIPT.json'),
 ('session67-root','verify_exports.py',[],'export_check.json'),
 ('session67-incoming','independent_review.py',[],'review.normal.out'),
 ('session67-incoming','score_fibre_check.py',[],'score-fibre.normal.out'),
 ('session67-incoming-root','check_dual_modules.py',[],'receipt.json')]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path);a=ap.parse_args();rows=[]
with tempfile.TemporaryDirectory(prefix='session67-frozen-replay-') as tmp:
 for optimized in (False,True):
  p=Path(tmp)/str(optimized);p.mkdir()
  for lane in {j[0] for j in jobs}:shutil.copytree(root/lane,p/lane)
  for lane,script,args,receipt in jobs:
   command=[sys.executable]+(['-O'] if optimized else [])+[script]+args
   result=subprocess.run(command,cwd=p/lane,env=dict(os.environ,PATH='/nonexistent',PYTHONDONTWRITEBYTECODE='1'),capture_output=True)
   if result.returncode or result.stdout!=(root/lane/receipt).read_bytes():raise RuntimeError(f'replay failed: {lane}/{script}, optimized={optimized}, {result.stderr!r}')
   rows.append({'lane':lane,'script':script,'optimized':optimized,'status':'PASS','stdout_sha256':hashlib.sha256(result.stdout).hexdigest(),'source_sha256':digest(root/lane/script),'PATH':'/nonexistent'})
  if a.repo:
   command=[sys.executable]+(['-O'] if optimized else [])+['audit_core_triangle_portable.py','--metagraph-dir','../session67-metagraph','--repo',str(a.repo.resolve())]
   result=subprocess.run(command,cwd=p/'session67-cycle-audit',env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True)
   if result.returncode or result.stdout!=(root/'session67-cycle-audit/CORE_TRIANGLE_PORTABLE_COMMIT_RECEIPT.json').read_bytes():raise RuntimeError('named-commit audit failed: '+repr(result.stderr))
   rows.append({'lane':'session67-cycle-audit','script':'audit_core_triangle_portable.py','optimized':optimized,'status':'PASS','named_commit_verified':True,'stdout_sha256':hashlib.sha256(result.stdout).hexdigest()})
print(json.dumps({'status':'PASS','runs':rows,'all_stdout_exact':True,'inputs_written_only_in_temporary_copies':True},indent=2))
