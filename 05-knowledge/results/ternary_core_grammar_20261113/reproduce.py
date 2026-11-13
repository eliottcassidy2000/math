"""Reproduce the research receipts; --remaining skips the already replayed long census."""
from pathlib import Path
import sys,subprocess,json,hashlib
HERE=Path(__file__).resolve().parent
JOBS=[('triangle_probe',['7']),('local_obstruction',[]),('peeling',[]),('atomic_ears',[]),('neutral_entry',[]),('audit',[])]

def main():
    remaining='--remaining' in sys.argv
    for name,args in JOBS:
        if remaining and name=='triangle_probe':continue
        modes=[True] if remaining and name!='audit' else [False,True]
        for opt in modes:
            cmd=[sys.executable]+(['-O'] if opt else [])+['-B','-X','utf8',str(HERE/(name+'.py'))]+args
            r=subprocess.run(cmd,capture_output=True,text=True,encoding='utf8')
            if r.returncode:print(r.stdout,r.stderr);raise RuntimeError(name)
            print(name,'optimized' if opt else 'normal','PASS',flush=True)
    pairs=[]
    for name,args in JOBS:
        x=(HERE/(name+'.json.gz')).read_bytes();y=(HERE/(name+'.optimized.json.gz')).read_bytes()
        if x!=y:raise RuntimeError('normal/optimized mismatch '+name)
        pairs.append(dict(name=name,sha256=hashlib.sha256(x).hexdigest()))
    (HERE/'reproduction.json').write_text(json.dumps(dict(status='PASS',scope='Six deterministic normal/optimized receipt pairs; independent controls are in AUDIT.md',pairs=pairs),indent=2)+'\n',encoding='utf8')
    print('ALL SIX RECEIPT PAIRS IDENTICAL')

if __name__=='__main__':main()
