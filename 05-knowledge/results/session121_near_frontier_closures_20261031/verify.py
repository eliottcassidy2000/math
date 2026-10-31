"""Complete public replay, including regenerated LRC evidence. Standard library only."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent

def sha(b):return hashlib.sha256(b).hexdigest()
def need(ok,label):
    if not ok:raise RuntimeError(label)
def pipeline(optimized,manifest):
    flags=['-B']+(['-O'] if optimized else [])
    def run(name,*args):
        result=subprocess.run([sys.executable,*flags,str(ROOT/name),*map(str,args)],capture_output=True,check=True)
        return result.stdout
    reports=[]
    def checked(name,*args):
        out=run(name,*args);json.loads(out)
        digest=sha(out)
        need(digest==manifest['expected_stdout_sha256'][name],'Receipt mismatch: '+name)
        reports.append({'script':name,'stdout_sha256':digest,'status':'PASS'})
    for name in manifest['light_checks']:checked(name)
    compact_path=ROOT/'lrc/triple_compact_receipt.json'
    compact=json.loads(compact_path.read_text())
    with tempfile.TemporaryDirectory(prefix='session121-replay-') as folder:
        full=Path(folder)/'full.json'; regenerated=Path(folder)/'compact.json'
        run('lrc/triple_certificate.py','--output',full)
        raw=full.read_bytes()
        need(sha(raw)==compact['complete_ledger_sha256'],'Regenerated full ledger hash')
        need(len(raw)==compact['complete_ledger_bytes'],'Regenerated full ledger length')
        run('lrc/compact_receipt.py',full,regenerated)
        need(regenerated.read_bytes()==compact_path.read_bytes(),'Regenerated compact receipt')
        checked('audit_lrc/check.py',full)
        checked('audit_lrc/compact_binding.py',full,compact_path)
        checked('audit_lrc/mutation_controls.py',full)
    return {'checks':reports,'full_ledger_sha256':compact['complete_ledger_sha256'],'compact_receipt_sha256':sha(compact_path.read_bytes()),'full_lrc_replay':'PASS','core_count':286,'second_joint_branches':574336,'finite_leaves':67585}
def main():
    manifest=json.loads((ROOT/'manifest.json').read_text())
    for name,digest in manifest['sha256'].items():need(sha((ROOT/name).read_bytes())==digest,'Sealed file mismatch: '+name)
    with ThreadPoolExecutor(max_workers=2) as pool:
        normal=pool.submit(pipeline,False,manifest);optimized=pool.submit(pipeline,True,manifest)
        a=normal.result();b=optimized.result()
    need(a==b,'Normal/optimized replay disagreement')
    print(json.dumps({'status':'PASS','normal_and_optimized':'identical',**a,'scope':'Complete declared finite replays plus separate analytic proofs; inherited E10 certificate and general Presburger compiler are not run'},indent=2,sort_keys=True))
if __name__=='__main__':main()
