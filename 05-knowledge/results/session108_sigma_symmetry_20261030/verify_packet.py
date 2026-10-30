"""Verify research packet bytes and isolated exact normal/optimized replays."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
JOBS=[('codes/verify_codes.py','codes/verification.json'),('arithmetic/verify.py','arithmetic/certificate.json'),('operations/verify_operations.py','operations/normal.json'),('tournaments/verify_palette.py','tournaments/verification.json'),('tournaments/verify_relations.py','tournaments/relations.verification.json'),('audit_arithmetic_codes/check.py','audit_arithmetic_codes/normal.json'),('audit_operations/check_independent.py','audit_operations/normal.json'),('history/check_covering_boundary.py','history/verification.json')]
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output');args=ap.parse_args();root=Path(__file__).resolve().parent
    m=json.loads((root/'manifest.json').read_text())
    for rel,digest in m['sha256'].items():need(hashlib.sha256((root/rel).read_bytes()).hexdigest()==digest,'manifest mismatch: '+rel)
    rows=[]
    with tempfile.TemporaryDirectory(prefix='sigma-symmetry-replay-') as td:
        cp=Path(td)/'packet';shutil.copytree(root,cp,ignore=shutil.ignore_patterns('__pycache__','*.local.json','SESSION_STATE.md'))
        for script,receipt in JOBS:
            expected=json.loads((cp/receipt).read_text());outs=[]
            for opt in (False,True):
                cmd=[sys.executable]+(['-O'] if opt else [])+[str(cp/script)]
                r=subprocess.run(cmd,cwd=(cp/script).parent,text=True,capture_output=True,timeout=300)
                need(r.returncode==0,script+' failed: '+r.stderr[-2000:]);v=json.loads(r.stdout);need(v==expected,script+' receipt mismatch');outs.append(v)
            need(outs[0]==outs[1],script+' optimized mismatch')
            rows.append({'script':script,'normal':'PASS','optimized':'PASS','receipt_sha256':hashlib.sha256((root/receipt).read_bytes()).hexdigest()})
    out={'status':'PASS','manifest_entries':len(m['sha256']),'relative_links':len(m['relative_links']),'isolated_runs':2*len(JOBS),'jobs':rows,'scope':'Exact finite replays; universal claims rely on audited proofs and qualified cited classification. No global sigma-square classification, Collatz proof, or CI assertion.'}
    text=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(text)
    print(text,end='')
if __name__=='__main__':main()
