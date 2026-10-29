"""Verify declared bytes, isolated normal/-O replays, and regenerated atlas."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

JOBS=[('recurrence/check_recurrence.py','recurrence/checks.json'),
      ('recurrence/check_upper_cutoff.py','recurrence/upper_normal.json'),
      ('spectral/check_spectral.py','spectral/checks.json'),
      ('observability/check_observability.py','observability/checks.json'),
      ('descent/check_descent_transfer.py','descent/checks.normal.json'),
      ('reverse/check_reverse.py','reverse/checks.json'),
      ('audit_transfer/check_independent.py','audit_transfer/checks.normal.json'),
      ('audit_reverse/check_independent.py','audit_reverse/normal.json')]

def require(ok,message):
    if not ok: raise RuntimeError(message)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output');args=parser.parse_args()
    root=Path(__file__).resolve().parent
    manifest=json.loads((root/'manifest.json').read_text())
    for rel,digest in manifest['sha256'].items():
        require(hashlib.sha256((root/rel).read_bytes()).hexdigest()==digest,'manifest mismatch: '+rel)
    results=[]
    with tempfile.TemporaryDirectory(prefix='collatz-transfer-replay-') as td:
        copy=Path(td)/'packet';shutil.copytree(root,copy,ignore=shutil.ignore_patterns('__pycache__','*.local.json','published_replay','SESSION_STATE.md'))
        for script,receipt in JOBS:
            expected=json.loads((copy/receipt).read_text()); outputs=[]
            for opt in (False,True):
                command=[sys.executable]+(['-O'] if opt else [])+[str(copy/script)]
                atlas=Path(td)/('atlas_opt.json' if opt else 'atlas_normal.json')
                if script.startswith('descent/'):
                    command+=['--atlas',str(atlas)]
                r=subprocess.run(command,cwd=copy,text=True,capture_output=True,timeout=300)
                require(r.returncode==0,script+' failed: '+r.stderr[-2000:])
                actual=json.loads(r.stdout);require(actual==expected,'receipt mismatch: '+script)
                outputs.append(actual)
                if script.startswith('descent/'):
                    require(atlas.read_bytes()==(copy/'descent/atlas16.json').read_bytes(),'regenerated atlas differs')
            require(outputs[0]==outputs[1],'optimized mismatch: '+script)
            results.append({'script':script,'normal':'PASS','optimized':'PASS','receipt_sha256':hashlib.sha256((root/receipt).read_bytes()).hexdigest()})
    result={'status':'PASS','manifest_entries':len(manifest['sha256']),'isolated_runs':2*len(JOBS),'atlas_regenerated_both_modes':True,'jobs':results,
            'scope':'Byte checks and finite exact replays; all-order claims depend on proofs and audits; no Collatz solution or CI assertion.'}
    data=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(data)
    print(data,end='')
if __name__=='__main__':main()
