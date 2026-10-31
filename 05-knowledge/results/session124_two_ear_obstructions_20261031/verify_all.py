#!/usr/bin/env python3
"""Portable exact replay. Default is read-only and checks bundled receipts."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys
ROOT=Path(__file__).resolve().parent
JOBS=[
 ('transfer','code/check.py',[],False),
 ('sole_outneighbor','code/check_no_go_identities.py',[],False),
 ('spectrum','code/spectral.py',[],False),
 ('ternary','code/verify_ternary_certificates.py',[],True),
 ('radix27','code/verify_radix27_farkas.py',[],True),
 ('radix81_units','code/verify_periodic_word_certificates.py',['4'],True),
 ('radix9_occupancy','code/verify_periodic_reachable_certificates.py',['2'],True),
 ('radix81_occupancy','code/verify_periodic_reachable_certificates.py',['4'],True),
 ('modular_mode','code/verify_alternating_mod3.py',[],True),
 ('independent_audit','audit/independent_audit.py',[],False),
 ('extended_audit','audit/extended_audit.py',[],False),
]
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--stdlib-only',action='store_true',help='run six independent integer/Fraction proof checkers without SymPy')
 p.add_argument('--write-receipts',action='store_true',help='replace bundled receipt files after successful matched normal/optimized runs')
 args=p.parse_args();results=[]
 for name,script,extra,stdlib in JOBS:
  if args.stdlib_only and not stdlib:continue
  outputs=[]
  for mode in ('normal','optimized'):
   cmd=[sys.executable,'-B']+(['-O'] if mode=='optimized' else [])+[str(ROOT/script)]+extra
   proc=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
   if proc.returncode:raise RuntimeError(f'{name}/{mode} failed: {proc.stderr.decode()}')
   if proc.stderr:raise RuntimeError(f'{name}/{mode} unexpected stderr: {proc.stderr.decode()}')
   outputs.append(proc.stdout)
  if outputs[0]!=outputs[1]:raise RuntimeError(f'{name}: normal/-O mismatch')
  for mode,output in zip(('normal','optimized'),outputs):
   target=ROOT/'receipts'/f'{name}.{mode}.out'
   if args.write_receipts:target.write_bytes(output)
   elif not target.is_file() or target.read_bytes()!=output:raise RuntimeError(f'{name}/{mode}: bundled receipt mismatch')
  results.append(dict(check=name,sha256=hashlib.sha256(outputs[0]).hexdigest(),normal_equals_optimized=True))
 print(json.dumps(dict(checks=results,total_checks=len(results),stdlib_only=args.stdlib_only),sort_keys=True))
if __name__=='__main__':main()
