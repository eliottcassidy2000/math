#!/usr/bin/env python3
"""Portable stdlib replay: fresh receipts in isolated normal and -O copies.

Run `python3 verify_all.py > verification_receipt.json` from any directory.
No repository metadata, external checkout, packages, or network are needed.
The default suite deliberately omits duplicate large optional censuses.
Only temporary files are written; the caller may redirect the summary.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
# id, producer, receipt, receipt is JSON stdout rather than a written file
CASES = [
 ('actual_language', 'languages/check_actual.py', 'languages/actual_receipt.json', False),
 ('vertex_dag', 'triangularity/check_avoiding_dag.py', 'triangularity/avoiding_dag_receipt.json', False),
 ('all_k_sparse', 'triangularity/check_all_k_sparse.py', 'triangularity/all_k_sparse_receipt.json', False),
 ('last_birth', 'bridge/check_last_birth_certificate.py', 'bridge/last_birth_receipt.json', False),
 ('audit1', 'audit/check_independent.py', 'audit/receipt.json', False),
 ('audit1_all_k_sink', 'audit/check_all_k_dag_and_sink.py', 'audit/all_k_dag_sink_receipt.json', False),
 ('audit2', 'audit2/check_independent.py', 'audit2/receipt.json', False),
 ('necklaces', 'necklaces/check_chord_necklaces.py', 'necklaces/results.json', True),
 ('audit1_final_additions', 'audit/check_final_additions.py', 'audit/final_additions_receipt.json', False),
]
SUPPORT = [
 'languages/source_rotation_witness.json',
 'singularity/actual_grammar.py',
 'singularity/native_kernel.py',
 'singularity/pinned/tournament_results.json',
 'bridge/star_probe.py',
]

def require(ok, message):
 if not ok: raise RuntimeError(message)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--only', nargs='+', choices=[x[0] for x in CASES], help='Optional selected subset; default is all key checks')
 parser.add_argument('--jobs', type=int, default=2, help='Concurrent subprocesses (default 2)')
 args=parser.parse_args()
 require(args.jobs>=1,'--jobs must be positive')
 selected=[x for x in CASES if args.only is None or x[0] in args.only]
 # Include imported producers and input receipts too, but never clone a checkout.
 manifest=sorted(set(SUPPORT+[p for case in CASES for p in case[1:3]]))
 before={p:sha(ROOT/p) for p in manifest}
 expected={case[0]:json.loads((ROOT/case[2]).read_text()) for case in selected}
 for p in manifest:
  if p.endswith('.py'):
   # Optimized replay must not silently disable verification assertions.
   import ast
   require(not any(isinstance(x,ast.Assert) for x in ast.walk(ast.parse((ROOT/p).read_text()))), 'assert-based verifier: '+p)
 env=os.environ.copy();env.pop('PYTHONPATH',None)
 env['PYTHONDONTWRITEBYTECODE']='1';env['PYTHONHASHSEED']='0'
 def run(task):
  case,optimized=task;name,script,receipt,stdout_json=case
  start=time.monotonic()
  # One copy per subprocess avoids cross-producer receipt or import contamination.
  with tempfile.TemporaryDirectory(prefix='session103-verify-') as td:
   dest=Path(td)
   for rel in manifest:
    target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,target)
   target=dest/receipt
   target.unlink() # mandatory freshness gate: no inherited successful output
   require(not target.exists(), 'receipt removal failed: '+receipt)
   command=[sys.executable]+(['-O'] if optimized else [])+[str(dest/script)]
   proc=subprocess.run(command,cwd=dest,env=env,text=True,capture_output=True)
   require(proc.returncode==0, f'{name} optimized={optimized} exited {proc.returncode}\n{proc.stdout}\n{proc.stderr}')
   if stdout_json:target.write_text(proc.stdout)
   require(target.is_file(),f'{name} did not regenerate {receipt}')
   fresh=json.loads(target.read_text())
   require(fresh==expected[name],f'{name} optimized={optimized}: fresh JSON differs from baseline {receipt}')
   print(f'PASS {name} {"-O" if optimized else "normal"}',file=sys.stderr,flush=True)
   return {'check':name,'optimized':optimized,'status':'PASS','fresh_json_equal':True,'receipt':receipt,'baseline_sha256':before[receipt],'elapsed_seconds':round(time.monotonic()-start,3)}
 tasks=[(case,opt) for case in selected for opt in (False,True)]
 with ThreadPoolExecutor(max_workers=args.jobs) as pool:results=list(pool.map(run,tasks))
 after={p:sha(ROOT/p) for p in manifest}
 require(before==after,'Original input/script/receipt changed during replay')
 out={'status':'PASS','scope':'Key proof checks; finite corroboration, not a replacement for all-order proofs','checks':len(selected),'runs':len(results),'normal_and_optimized':True,'fresh_receipts':True,'original_manifest_unchanged':True,'manifest_sha256':before,'results':results}
 print(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
