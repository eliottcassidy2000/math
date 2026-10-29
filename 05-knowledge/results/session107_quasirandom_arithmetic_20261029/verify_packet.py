"""Manifest and isolated exact normal/optimized verification."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
JOBS=[('characters/verify_characters.py','characters/verification.json'),('conditioning/check_conditioning.py','conditioning/normal.json'),('groups/check_groups.py','groups/check_groups.json'),('orders/check_orders.py','orders/verification.json'),('recursions/verify_recursions.py','recursions/checks.json'),('audit_arithmetic/check_independent.py','audit_arithmetic/normal.json'),('audit_structure/check_independent.py','audit_structure/verification.json'),('audit_structure/check_fourier.py','audit_structure/fourier_verification.json'),('incoming/check_concurrent_extensions.py','incoming/concurrent_extensions_receipt.json'),('audit_structure/check_concurrent.py','audit_structure/concurrent_verification.json'),('colored/verify_colored.py','colored/verification.json'),('audit_structure/check_colored.py','audit_structure/colored_verification.json')]
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output');args=ap.parse_args();root=Path(__file__).resolve().parent
    m=json.loads((root/'manifest.json').read_text())
    for rel,digest in m['sha256'].items():need(hashlib.sha256((root/rel).read_bytes()).hexdigest()==digest,'manifest '+rel)
    for rel,digest in json.loads((root/'audit_structure/source_hashes.json').read_text()).items():need(hashlib.sha256((root/rel).read_bytes()).hexdigest()==digest,'audit input '+rel)
    results=[]
    with tempfile.TemporaryDirectory(prefix='quasirandom-replay-') as td:
        cp=Path(td)/'packet';shutil.copytree(root,cp,ignore=shutil.ignore_patterns('__pycache__','*.local.json','SESSION_STATE.md'))
        for script,receipt in JOBS:
            expected=json.loads((cp/receipt).read_text());outputs=[]
            for opt in (False,True):
                cmd=[sys.executable]+(['-O'] if opt else [])+[str(cp/script)]
                r=subprocess.run(cmd,cwd=(cp/script).parent,text=True,capture_output=True,timeout=300)
                need(r.returncode==0,script+' failed '+r.stderr[-2000:]);got=json.loads(r.stdout);need(got==expected,script+' receipt mismatch');outputs.append(got)
            need(outputs[0]==outputs[1],script+' optimization mismatch')
            results.append({'script':script,'normal':'PASS','optimized':'PASS','receipt_sha256':hashlib.sha256((root/receipt).read_bytes()).hexdigest()})
    out={'status':'PASS','manifest_entries':len(m['sha256']),'relative_links':len(m['relative_links']),'isolated_runs':2*len(JOBS),'jobs':results,'scope':'Exact finite replays and hash verification; infinite claims rely on audited proofs and qualified cited inputs. No CI or Collatz convergence assertion.'}
    text=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(text)
    print(text,end='')
if __name__=='__main__':main()
