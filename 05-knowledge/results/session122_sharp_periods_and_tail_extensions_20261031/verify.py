"""Complete published finite replays, standard library, normal and optimized modes."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,subprocess,sys,tempfile
from materialize_receipts import receipt_bytes
ROOT=Path(__file__).resolve().parent

def sha(raw):return hashlib.sha256(raw).hexdigest()
def need(ok,why):
    if not ok:raise RuntimeError(why)
def pipeline(opt,m):
    flags=['-B']+(['-O'] if opt else [])
    def run(script,*args):
        return subprocess.run([sys.executable,*flags,str(ROOT/script),*map(str,args)],capture_output=True,check=True).stdout
    def check(script,receipt,args):
        out=run(script,*args);json.loads(out)
        need(out==(ROOT/receipt).read_bytes(),'output receipt mismatch: '+script)
        return dict(script=script,stdout_sha256=sha(out),status='PASS')
    checks=[]
    for c in m['checks']:checks.append(check(c['script'],c['receipt'],[ROOT/a for a in c['args']]))
    replays=[]
    with tempfile.TemporaryDirectory(prefix='session122-replay-') as tmp:
        for r in m['lrc_replays']:
            out=Path(tmp)/('tails%d.json'%r['tails'])
            run(r['producer'],'--output',out)
            need(out.read_bytes()==receipt_bytes(r['receipt']),'regenerated full-tree commitment mismatch')
            checks.append(check(r['auditor'],r['audit_receipt'],[out]))
            checks.append(check(r['controls'],r['controls_receipt'],[out]))
            data=json.loads(out.read_bytes())
            replays.append(dict(tails=r['tails'],core_count=data['core_count'],receipt_sha256=sha(out.read_bytes()),status='PASS'))
    return dict(checks=checks,lrc_replays=replays)
def main():
    m=json.loads((ROOT/'manifest.json').read_text())
    for name,h in m['sha256'].items():need(sha((ROOT/name).read_bytes())==h,'sealed file mismatch: '+name)
    with ThreadPoolExecutor(max_workers=2) as pool:
        a=pool.submit(pipeline,False,m);b=pool.submit(pipeline,True,m);x=a.result();y=b.result()
    need(x==y,'normal/-O disagreement')
    print(json.dumps(dict(status='PASS',normal_and_optimized='identical',**x,scope='All declared new Python finite certificates regenerated and independently replayed. Analytic proofs are separately reviewed. Inherited large Collatz certificates and the prior Session121 LRC certificate are not rerun by this command; no Lean or human refereeing claim.'),sort_keys=True,indent=2))
if __name__=='__main__':main()
