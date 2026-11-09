#!/usr/bin/env python3
"""Portable exact replay; no third-party packages, no assertions."""
import hashlib,json,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CASES=[('recursive','verify.py','CHECK'),('recursive_audit','audit.py','CHECK'),('recursive_complete','verify.py','CHECK'),('recursive_complete','ternary.py','TERNARY'),('ternary_audit','audit.py','CHECK'),('overlap','verify.py','CHECK'),('overlap_audit','audit.py','CHECK'),('inverse','verify.py','CHECK'),('inverse_audit','check.py','CHECK'),('dynamic','verify.py','CHECK'),('dynamic_audit','check.py','CHECK'),('rich','verify.py','CHECK'),('rich_audit','check.py','CHECK'),('rich','extend_triples.py','EXTEND_TRIPLES'),('rich','extend_unsaturated.py','EXTEND_UNSATURATED'),('extension_audit','check.py','CHECK')]
def main():
 rows=[]
 for directory,script,receipt in CASES:
  for optimized in (False,True):
   source=ROOT/directory/script
   expected=ROOT/directory/(receipt+('.optimized' if optimized else '')+'.json')
   args=[sys.executable,'-I','-B']+(['-O'] if optimized else [])+[str(source)]
   p=subprocess.run(args,cwd=source.parent,capture_output=True)
   if p.returncode:raise RuntimeError(f'{directory}/{script} failed: {p.stderr.decode()}')
   if p.stdout!=expected.read_bytes():raise RuntimeError(f'{directory}/{script} receipt mismatch, optimized={optimized}')
   rows.append({'script':str(source.relative_to(ROOT)),'optimized':optimized,'script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'receipt_sha256':hashlib.sha256(p.stdout).hexdigest()})
 for directory,source_name,receipt,stderr_name in [('dynamic','search.cpp','SEARCH.txt','SEARCH.stats'),('dynamic_audit','census.cpp','CENSUS.txt',None)]:
  source=ROOT/directory/source_name
  for optimized in (False,True):
   with tempfile.TemporaryDirectory(prefix='session229-replay-') as tmp:
    exe=Path(tmp)/'census'
    c=subprocess.run(['g++','-std=c++17','-O3' if optimized else '-O0',str(source),'-o',str(exe)],capture_output=True)
    if c.returncode:raise RuntimeError(c.stderr.decode())
    p=subprocess.run([str(exe)],capture_output=True)
    if p.returncode or p.stdout!=(ROOT/directory/receipt).read_bytes():raise RuntimeError(f'{directory} C++ census mismatch')
    if stderr_name and p.stderr!=(ROOT/directory/stderr_name).read_bytes():raise RuntimeError(f'{directory} C++ statistics mismatch')
    rows.append({'script':str(source.relative_to(ROOT)),'optimized':optimized,'script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'receipt_sha256':hashlib.sha256(p.stdout).hexdigest()})
 print(json.dumps({'status':'PASS' ,'canonical_runs':len(rows),'runs':rows},indent=2,sort_keys=True))
if __name__=='__main__':main()
