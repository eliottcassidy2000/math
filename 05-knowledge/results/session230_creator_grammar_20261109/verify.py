#!/usr/bin/env python3
"""Rebuild all finite checks in an isolated temporary packet. Standard library only."""
import argparse, hashlib, json, pathlib, shutil, subprocess, sys, tempfile
P=pathlib.Path(__file__).resolve().parent
CASES=[('audit/check.py','audit/check.json'),('audit/check_witnesses.py','audit/check_witnesses.json'),('cut/check.py','cut/RESULTS.json'),('cut/check_transitive.py','cut/TRANSITIVE.json'),('cut/audit/check.py','cut/audit/RESULTS.json'),('primes/probe.py','primes/probe.out'),('primes/audit/check.py','primes/audit/check.out'),('primes/audit/check_factors.py','primes/audit/check_factors.out')]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--optimized',action='store_true');args=ap.parse_args()
 with tempfile.TemporaryDirectory(prefix='session230-replay-') as td:
  q=pathlib.Path(td)/'packet';shutil.copytree(P,q,ignore=shutil.ignore_patterns('__pycache__','classes_*.txt','classgen'))
  subprocess.run(['gcc','-O3','classgen.c','-o','classgen'],cwd=q/'primes',check=True,capture_output=True)
  gen=subprocess.run(['./classgen','8','classes'],cwd=q/'primes',check=True,capture_output=True)
  rows=[]
  for script,receipt in CASES:
   command=[sys.executable]+(['-O'] if args.optimized else [])+[str(q/script)]
   r=subprocess.run(command,cwd=(q/script).parent,check=True,capture_output=True)
   expected=(P/receipt).read_bytes()
   if r.stdout!=expected:raise RuntimeError('Receipt mismatch: '+script+'\n'+r.stderr.decode(errors='replace'))
   rows.append({'script':script,'status':'PASS','receipt_sha256':hashlib.sha256(expected).hexdigest()})
  print(json.dumps({'status':'PASS','checks':rows,'class_generator_sha256':hashlib.sha256((P/'primes/classgen.c').read_bytes()).hexdigest(),'generated_class_counts':[1,1,2,4,12,56,456,6880],'scope':'All new proof controls and independent order-eight orbit certificate; order-eleven witness soundness only, not full large-census replay.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
