#!/usr/bin/env python3
"""Read-only, source-relative standard-library replay of the curated proof packet."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent
JOBS=[('farkas','code/verify_farkas.py'),('occupancy_farkas','code/verify_occupancy_farkas.py'),('complement_identity','code/verify_complement_identity.py'),('independent_audit','audit/audit.py'),('spectral_example','audit/spectral_and_example.py')]

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--write-receipts',action='store_true')
 args=parser.parse_args();results=[]
 if not args.write_receipts:
  listed=set((ROOT/'WHITELIST.txt').read_text().splitlines())
  covered=set()
  for line in (ROOT/'SHA256SUMS').read_text().splitlines():
   digest,name=line.split('  ',1)
   if name not in listed or name in covered:raise RuntimeError('invalid manifest coverage')
   if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise RuntimeError(f'manifest mismatch: {name}')
   covered.add(name)
  if covered!=listed-{'SHA256SUMS'}:raise RuntimeError('incomplete manifest')
 for name,script in JOBS:
  outputs=[]
  for mode in ('normal','optimized'):
   command=[sys.executable,'-B']+(['-O'] if mode=='optimized' else [])+[str(ROOT/script)]
   proc=subprocess.run(command,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
   if proc.returncode or proc.stderr:raise RuntimeError(f'{name}/{mode} failed: {proc.stderr.decode()}')
   outputs.append(proc.stdout)
  if outputs[0]!=outputs[1]:raise RuntimeError(f'{name}: normal/-O mismatch')
  for mode,output in zip(('normal','optimized'),outputs):
   target=ROOT/'receipts'/f'{name}.{mode}.out'
   if args.write_receipts:target.write_bytes(output)
   elif not target.is_file() or target.read_bytes()!=output:raise RuntimeError(f'{name}/{mode}: receipt mismatch')
  results.append(dict(check=name,sha256=hashlib.sha256(outputs[0]).hexdigest(),normal_equals_optimized=True))
 print(json.dumps(dict(status='passed',standard_library_only=True,total_checks=len(results),checks=results),sort_keys=True))
if __name__=='__main__':main()
