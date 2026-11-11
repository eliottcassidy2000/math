"""Deterministic normal/optimized replay of the new integrated proof controls."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys

BASE=Path(__file__).resolve().parent
CHECKS=[
 ('audit/retained_reset.py','audit/retained_reset.json'),
 ('audit/template_controls.py','audit/template_controls.json'),
 ('audit/weak_factor.py','audit/weak_factor.json'),
 ('audit/ready_rebase.py','audit/ready_rebase.json'),
 ('audit/atlas_instructions.py','audit/atlas_instructions.json'),
 ('defect/cycle_residual.py','defect/cycle_residual.json'),
 ('defect/pentagon_inflation.py','defect/pentagon_inflation.json'),
 ('arithmetic/dyadic_defect.py','arithmetic/dyadic_defect.json'),
 ('arithmetic/prefix_interface.py','arithmetic/prefix_interface.json'),
 ('arithmetic/tail_squares.py','arithmetic/tail_squares.json'),
 ('arithmetic/chart_inverse.py','arithmetic/chart_inverse.json'),
 ('arithmetic/minimum_squares.py','arithmetic/minimum_squares.json'),
 ('metagraph/minimal_defect.py','metagraph/minimal_defect.json'),
 ('metagraph/response_families.py','metagraph/response_families.json'),
 ('reset/uniform_return.py','reset/uniform_return.json'),
 ('reset/composition.py','reset/composition.json'),
 ('reset/ready_rebase.py','reset/ready_rebase.json'),
 ('reset/iterable_return.py','reset/iterable_return.json'),
 ('reset/atlas_returns.py','reset/atlas_returns.json'),
 ('reset/atlas_two_birth.py','reset/atlas_two_birth.json'),
]


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--only',nargs='+',help='Replay selected scripts and retain earlier checked entries')
    args=parser.parse_args()
    previous={}
    if args.only and (BASE/'REPLAY.json').exists():
        previous={x['script']:x for x in json.loads((BASE/'REPLAY.json').read_text())['checks']}
    known={s for s,_ in CHECKS}
    if args.only and not set(args.only)<=known:raise RuntimeError('Unknown replay selection')
    for script,receipt in CHECKS:
        if args.only and script not in args.only:continue
        path=BASE/script
        before=digest(BASE/receipt)
        hashes=[]
        for mode in ([],['-O']):
            run=subprocess.run([sys.executable,'-B',*mode,str(path)],cwd=path.parent,
                               stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
            if run.returncode:
                raise RuntimeError(script+' failed: '+run.stdout+'\n'+run.stderr)
            hashes.append(digest(BASE/receipt))
        if hashes[0]!=hashes[1] or hashes[0]!=before:
            raise RuntimeError('Receipt changed or normal/-O mismatch: '+receipt)
        previous[script]=dict(script=script,receipt=receipt,sha256=hashes[0],normal_optimized='MATCH')
        print('PASS '+script,flush=True)
    results=[previous[s] for s,_ in CHECKS if s in previous]
    report=dict(status='PASS',method='Normal and optimized subprocess replay, matching existing normalized receipts',checks=results)
    (BASE/'REPLAY.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('PASS: all'+str(len(results))+' integrated checks match normal/optimized and saved evidence',flush=True)


if __name__=='__main__':main()
