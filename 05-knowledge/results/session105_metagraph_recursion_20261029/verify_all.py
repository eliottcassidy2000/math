#!/usr/bin/env python3
"""Fresh isolated normal/-O packet replay; preserve mathematical, not timing, data."""
import ast,hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CASES=[
 ('enumeration/verify.py','enumeration/verification.json'),
 ('generator/check_preparation.py','generator/results.json'),
 ('metagraph/check_weighted_transport.py','metagraph/checks.json'),
 ('updown/check_updown.py','updown/check_normal.json'),
 ('harmonic/check_harmonic.py','harmonic/checks.json'),
 ('updown/compressed_smith.py','updown/smith_atlas_200.json'),
 ('audit_native/check_independent.py','audit_native/receipt.json'),
 ('audit_recurrence/check_independent.py','audit_recurrence/normal.json')]
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def normalize(x,timing=False):
 if isinstance(x,list):return [normalize(v,timing) for v in x]
 if isinstance(x,dict):return {k:normalize(v,timing) for k,v in x.items() if k not in ['python_optimized','optimized'] and not(timing and (k in ['timings','median_seconds','seconds'] or k.endswith('_seconds')))}
 return x
def main():
 files=json.loads((ROOT/'FILES.json').read_text())
 files=[p for p in files if p not in ['verification.json','MANIFEST.json']]
 before={p:digest(ROOT/p) for p in files}
 records=[]
 for script,receipt in CASES:
  need(not any(isinstance(n,ast.Assert) for n in ast.walk(ast.parse((ROOT/script).read_text()))),'assert check disabled by -O: '+script)
  timing=script.startswith('generator/')
  expected=normalize(json.loads((ROOT/receipt).read_text()),timing)
  answers=[]
  for optimized in [False,True]:
   with tempfile.TemporaryDirectory(prefix='tournament-recursion-') as td:
    dest=Path(td)
    for p in files:
     target=dest/p;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/p,target)
    # Every selected producer emits new stdout, so no old receipt is read as evidence.
    (dest/receipt).unlink()
    env=os.environ.copy();env.pop('PYTHONPATH',None);env['PYTHONDONTWRITEBYTECODE']='1';env['PYTHONHASHSEED']='0';env['PYTHONINTMAXSTRDIGITS']='0'
    result=subprocess.run([sys.executable]+(['-O'] if optimized else [])+[str(dest/script)],cwd=dest,env=env,text=True,capture_output=True)
    need(result.returncode==0,script+' failed: '+result.stderr+'\n'+result.stdout[:1000])
    fresh=json.loads(result.stdout)
    for key in ['python_optimized','optimized']:
     if key in fresh:need(fresh[key]==optimized,script+' wrong optimization flag')
    fresh=normalize(fresh,timing)
    need(fresh==expected,script+' differs from committed mathematical receipt')
    answers.append(fresh)
    records.append({'script':script,'optimized':optimized,'status':'PASS','fresh_stdout':True,'baseline_equal':True,'timing_fields_excluded':timing,'mathematical_json_sha256':hashlib.sha256(json.dumps(fresh,sort_keys=True).encode()).hexdigest()})
    print('PASS',script,'-O' if optimized else 'normal',file=sys.stderr,flush=True)
  need(answers[0]==answers[1],script+' mathematical modes disagree')
 need(before=={p:digest(ROOT/p) for p in files},'packet inputs changed during replay')
 print(json.dumps({'status':'PASS','checkers':len(CASES),'runs':len(records),'input_sha256':before,'isolated_temp_copies':True,'original_inputs_unchanged':True,'results':records},indent=2,sort_keys=True))
if __name__=='__main__':main()
