#!/usr/bin/env python3
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
EXPECTED='e1a91122cbc42f9de3edfc8e650959a1d79f55158e796c007b8dc34fd31432d9'
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inventory(root):return {str(p.relative_to(root)):digest(p) for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
def run(script,opt=False,args=(),success=True):
    p=subprocess.run([sys.executable]+(['-O'] if opt else [])+[str(script),*map(str,args)],cwd='/tmp',capture_output=True,text=True)
    need((p.returncode==0)==success,'unexpected subprocess: '+str(script)+'\n'+p.stderr[-1500:]);return p

def main():
    ap=argparse.ArgumentParser();ap.add_argument('candidate',type=Path);arg=ap.parse_args();candidate=arg.candidate.resolve()
    need(digest(candidate/'MANIFEST.json')==EXPECTED,'exact frozen candidate');before=inventory(candidate)
    result={'status':'PASS','candidate_manifest_sha256':EXPECTED,'candidate_inventory':before,'checks':[]}
    for opt in (False,True):run(candidate/'verify_manifest.py',opt)
    with tempfile.TemporaryDirectory(prefix='session51-score-audit-') as td:
        p=Path(td);audit=p/'actual-audit-copy';shutil.copytree(HERE,audit,ignore=shutil.ignore_patterns('__pycache__'))
        need(not audit.is_symlink(),'audit actual copy')
        for opt in (False,True):need(run(audit/'independent.py',opt).stdout==(HERE/'normal.json').read_text(),'independent actual-copy normal/-O byte replay')
        result['checks'].append('independent actual-copy normal/-O replay from /tmp byte-identical')
        producer=p/'actual-candidate-copy';shutil.copytree(candidate,producer,ignore=shutil.ignore_patterns('__pycache__'))
        for opt in (False,True):
            run(producer/'verify_manifest.py',opt);out=p/('candidate-O.json' if opt else 'candidate-N.json');run(producer/'check.py',opt,('--output',out))
            need(out.read_bytes()==(candidate/'normal.json').read_bytes(),'producer actual-copy byte replay')
        result['checks'].append('candidate strict manifest and actual-copy normal/-O replay verified')
        inp=audit/'sources/PATH_TRANSPORT_AND_FOREST_RULE.md';old=inp.read_bytes();inp.write_bytes(old+b' ')
        for opt in (False,True):need('source SHA256 mismatch' in run(audit/'independent.py',opt,success=False).stderr,'source pin rejection')
        inp.write_bytes(old);result['checks'].append('changed inherited transport source rejected normally and under -O')
        script=audit/'independent.py';old=script.read_text()
        for name,needle,replacement,diagnostic in [
            ('omit black same-score endpoint value','need(jump in {-1,-k-1,k-1}', 'need(jump in {-k-1,k-1}', 'actual unrestricted necessary triple list'),
            ('incorrect interval gap coefficient','ll-hi==(p-2)*(k-1)>0','ll-hi==(p-1)*(k-1)>0','strict complete interval separation'),
            ('promote formal value to nonnegative attainment','c-k-1==-2 and u>=0','c-k-1>=0 and u>=0','unattained negative formal value hostile'),
            ('replace C by identity','T[i][j] if abs(i-j)==1 else 1-T[i][j]','T[i][j] if abs(i-j)==1 else T[i][j]','actual full endpoint degree transport')]:
            need(old.count(needle)==1,'unique mutation anchor: '+name);script.write_text(old.replace(needle,replacement))
            for opt in (False,True):need(diagnostic in run(script,opt,success=False).stderr,'semantic rejection: '+name)
            script.write_text(old);result['checks'].append('semantic hostile rejected normally and under -O: '+name)
        theorem=producer/'THEOREM.md';old=theorem.read_bytes();theorem.write_bytes(old+b'\nMutation.\n')
        for opt in (False,True):need('SHA256 mismatch' in run(producer/'verify_manifest.py',opt,success=False).stderr,'candidate theorem mutation rejected')
        result['checks'].append('candidate theorem corruption rejected normally and under -O')
    need(inventory(candidate)==before,'original frozen candidate unchanged')
    result['checks'].append('all original candidate files unchanged; no publication')
    result['independent_receipt_sha256']=digest(HERE/'normal.json');result['candidate_receipt_sha256']=digest(candidate/'normal.json')
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
