#!/usr/bin/env python3
"""Replay every independent producer/audit in fresh isolated normal/-O runs."""
import ast,hashlib,json,os,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CASES=[
 ('exterior/check_exterior.py','exterior/results.json'),
 ('descent/check_descent.py','descent/checks.json'),
 ('positive/check_positive.py','positive/results.json'),
 ('summands/check_summands.py','summands/verification.json'),
 ('summands/check_three_chains.py','summands/three_chains_verification.json'),
 ('audit_algebra/check_independent.py','audit_algebra/receipt.normal.json'),
 ('audit_descent/check_independent.py','audit_descent/receipt.json'),
 ('audit_incoming/check_common_rank.py','audit_incoming/receipt.json')]
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def normalized(d):
 d=dict(d);d.pop('optimized',None);return d
def main():
 before={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for c in CASES for p in c}
 rows=[]
 for script,receipt in CASES:
  src=(ROOT/script).read_text();need(not any(isinstance(n,ast.Assert) for n in ast.walk(ast.parse(src))),script+' disables checks under -O')
  expected=json.loads((ROOT/receipt).read_text());answers=[]
  for opt in [False,True]:
   with tempfile.TemporaryDirectory(prefix='session104-') as td:
    target=Path(td)/Path(script).name;target.write_text(src)
    env=os.environ.copy();env.pop('PYTHONPATH',None);env['PYTHONHASHSEED']='0';env['PYTHONDONTWRITEBYTECODE']='1'
    proc=subprocess.run([sys.executable]+(['-O'] if opt else [])+[str(target)],cwd=td,env=env,text=True,capture_output=True)
    need(proc.returncode==0,script+' failed: '+proc.stderr)
    answer=json.loads(proc.stdout)
    if 'optimized' in answer:need(answer['optimized']==opt,'mode flag')
    need(normalized(answer)==normalized(expected),script+' disagrees with receipt')
    answers.append(normalized(answer));rows.append({'script':script,'optimized':opt,'fresh_stdout':True,'baseline_equal':True})
    print('PASS',script,'-O' if opt else 'normal',file=sys.stderr,flush=True)
  need(answers[0]==answers[1],script+' differs across optimization')
 after={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for c in CASES for p in c}
 need(before==after,'Original scripts/receipts changed during replay')
 print(json.dumps({'status':'PASS','checkers':len(CASES),'fresh_runs':len(rows),'original_files_unchanged':True,'input_hashes':before,'runs':rows},indent=2,sort_keys=True))
if __name__=='__main__':main()
