"""Run every local verifier normally and with -O, require byte-identical receipts."""
from pathlib import Path
import subprocess,sys,json,hashlib
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ENGINES=['fibre_probe','inflation_probe','degree_probe','repair23','prefix_selector','context_transport','direction_probe','transport_words','arithmetic_probe','audit']
def main():
    pairs=[]
    for name in ENGINES:
        for optimized in (False,True):
            cmd=[sys.executable,'-B']+(['-O'] if optimized else [])+['-X','utf8',str(HERE/(name+'.py'))]
            subprocess.run(cmd,cwd=ROOT,check=True)
        a=HERE/(name+'.json.gz');b=HERE/(name+'.optimized.json.gz')
        if a.read_bytes()!=b.read_bytes():raise RuntimeError('Receipt mismatch: '+name)
        pairs.append(dict(engine=name,sha256=hashlib.sha256(a.read_bytes()).hexdigest()))
    inherited=['ternary_core_grammar_20261113/triangle_probe.py','ternary_core_grammar_20261113/atomic_ears.py','ternary_core_grammar_20261113/peeling.py','multicore_metagraph_20261113/rooted_growth.py','multicore_metagraph_20261113/rooted_growth.json.gz','modular_centers_collatz_pell_20261111/content_checks.py','modular_centers_collatz_pell_20261111/tournament_checks.py']
    sources={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))}
    for rel in inherited:
        p=HERE.parent/rel;sources[str(p.relative_to(ROOT)).replace('\\','/')]=hashlib.sha256(p.read_bytes()).hexdigest()
    import numpy,numba
    (HERE/'REPLAY.json').write_text(json.dumps(dict(status='PASS',python=sys.version,numpy=numpy.__version__,numba=numba.__version__,pairs=pairs,sources=sources,independent_agent_review=False),indent=2)+'\n',encoding='utf-8')
    print('PASS: ten normal/optimized receipt pairs; all explicit inherited executable dependencies pinned')
if __name__=='__main__':main()
