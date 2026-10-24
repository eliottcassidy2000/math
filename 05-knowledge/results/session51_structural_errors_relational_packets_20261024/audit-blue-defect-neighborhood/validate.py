#!/usr/bin/env python3
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
EXPECTED='141cbcb8e2696238c2db7db40e4bcd9570438d791a4c405e662427671db8f6da'
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
    with tempfile.TemporaryDirectory(prefix='session51-neighbor-audit-') as td:
        p=Path(td);audit=p/'actual-audit-copy';shutil.copytree(HERE,audit,ignore=shutil.ignore_patterns('__pycache__'))
        need(not audit.is_symlink(),'audit is actual copy')
        for opt in (False,True):need(run(audit/'independent.py',opt).stdout==(HERE/'normal.json').read_text(),'independent actual-copy normal/-O byte replay')
        result['checks'].append('independent actual-copy normal/-O replay from /tmp is byte-identical')
        producer=p/'actual-candidate-copy';shutil.copytree(candidate,producer,ignore=shutil.ignore_patterns('__pycache__'))
        for opt in (False,True):
            run(producer/'verify_manifest.py',opt);out=p/('candidate-O.json' if opt else 'candidate-N.json');run(producer/'check.py',opt,('--output',out))
            need(out.read_bytes()==(candidate/'normal.json').read_bytes(),'producer actual-copy byte replay')
        result['checks'].append('candidate strict manifest and actual-copy normal/-O replay verified')
        # Integrity hostile.
        inp=audit/'sources/independent-wiring-receipt.json';old=inp.read_bytes();inp.write_bytes(old+b' ')
        for opt in (False,True):need('pinned source corruption' in run(audit/'independent.py',opt,success=False).stderr,'source pin rejection')
        inp.write_bytes(old);result['checks'].append('changed source rejected normally and under -O')
        # Semantic hostile: erroneously collapse full permutation search to identity.
        script=audit/'independent.py';old=script.read_text();needle='tuple(permutations(c))';need(old.count(needle)==1,'unique exhaustive-order mutation')
        script.write_text(old.replace(needle,'(c,)'))
        for opt in (False,True):need('1031 distinct boundary classes' in run(script,opt,success=False).stderr,'exhaustive canonicalization mutation rejection')
        script.write_text(old);result['checks'].append('identity-only class grouping bug rejected normally and under -O')
        # The next hostile deliberately repins a disposable mutated receipt to
        # exercise semantic checks rather than simply trigger hash rejection.
        receipt=audit/'sources/candidate-receipt.json';original_receipt=receipt.read_text();pins=audit/'SOURCE_PINS.json';original_pins=pins.read_text()
        bad=json.loads(original_receipt);bad['boundary_edge_witnesses'][0]['canonical_vertex_order'].reverse();receipt.write_text(json.dumps(bad))
        pd=json.loads(original_pins);pd['sources/candidate-receipt.json']=digest(receipt);pins.write_text(json.dumps(pd))
        for opt in (False,True):need('producer actual vertex-coordinate arc-isomorphism witness' in run(script,opt,success=False).stderr,'wrong actual coordinate witness rejection')
        receipt.write_text(original_receipt);pins.write_text(original_pins);result['checks'].append('repinned wrong actual canonical vertex witness rejected semantically normally and under -O')
        theorem=producer/'THEOREM.md';old=theorem.read_bytes();theorem.write_bytes(old+b'\nMutation.\n')
        for opt in (False,True):need('SHA256 mismatch' in run(producer/'verify_manifest.py',opt,success=False).stderr,'candidate theorem corruption rejection')
        result['checks'].append('candidate theorem corruption rejected normally and under -O')
    need(inventory(candidate)==before,'original frozen candidate remained unchanged')
    result['checks'].append('all original candidate files unchanged; no publication')
    result['independent_receipt_sha256']=digest(HERE/'normal.json');result['candidate_receipt_sha256']=digest(candidate/'normal.json')
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
