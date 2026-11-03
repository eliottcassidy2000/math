"""Reproduce the scoped proofs and compare normal/optimized exact receipts."""
from pathlib import Path
import hashlib,json,subprocess,sys,platform
import sympy
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def need(ok,why):
    if not ok:raise RuntimeError(why)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    pins=json.loads((HERE/'source_pins.json').read_text(encoding='utf8'))
    for relative,expected in pins['sha256'].items():
        need(digest(ROOT/relative)==expected,('inherited source changed: inspect before repinning',relative))
    inherited=json.loads((ROOT/'05-knowledge/results/mixed_returns_20261103/controllers.json').read_text())
    need(inherited['hostile81_27']['groups']==[81,27] and inherited['one_birth55_orientations']==768,'pinned hostile scope')
    scripts={'verify.py':'verification.json','geometry.py':'geometry.json',
             'jump_carrier.py':'jump_carrier.json','recurrences.py':'recurrences.json'}
    modes={}
    for mode in ('normal','optimized'):
        receipts={}
        for script,output in scripts.items():
            command=[sys.executable,'-X','utf8','-B']+(['-O'] if mode=='optimized' else [])+[str(HERE/script)]
            result=subprocess.run(command,cwd=ROOT,text=True,encoding='utf8',capture_output=True,check=True,timeout=180)
            print(mode,script,result.stdout.strip(),flush=True)
            receipts[output]=digest(HERE/output)
        modes[mode]=receipts
    need(modes['normal']==modes['optimized'],'normal/optimized receipts differ')
    report={'status':'PASS','python':platform.python_version(),'sympy':sympy.__version__,
            'source_pin_count':len(pins['sha256']),'same_agent_independent_implementations':True,
            'receipt_sha256_by_mode':modes,
            'inherited_controls':{'source':'mixed_returns_20261103/controllers.json',
                                  'recomputed_here':False,'hostile81_27':inherited['hostile81_27'],
                                  'one_birth55_orientations':inherited['one_birth55_orientations']}}
    (HERE/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('PASS: all scoped certificates; identical normal/optimized receipts; pinned source controls')
if __name__=='__main__':main()
