"""Reproduce the complete packet; optional optimization-only replay of existing receipts."""
from pathlib import Path
import subprocess,sys,hashlib,json
HERE=Path(__file__).resolve().parent
JOBS=[('rooted_growth',['8']),('cycle_index',['30']),('family_search',['--core-order','8','--fixed','12']),('next_target',['18']),('clock_probe',[]),('divisor_entry',['15']),('family_audit',[]),('blue_incidence',['7']),('component_growth',[]),('growth_audit',[]),('inflation_descent',[])]

def run():
    pairs=[]
    for name,args in JOBS:
        for flags in ([['-O']] if '--optimized-only' in sys.argv else [[],['-O']]):
            result=subprocess.run([sys.executable,'-B','-X','utf8',*flags,str(HERE/(name+'.py')),*args],capture_output=True,text=True,encoding='utf8')
            if result.returncode:
                print(result.stdout);print(result.stderr);raise RuntimeError(name+' failed')
        normal=HERE/(name+'.json.gz');optimized=HERE/(name+'.optimized.json.gz')
        if normal.read_bytes()!=optimized.read_bytes():raise RuntimeError(name+' replay mismatch')
        pairs.append(dict(name=name,sha256=hashlib.sha256(normal.read_bytes()).hexdigest()))
        print('PAIRED',name,flush=True)
    (HERE/'reproduction.json').write_text(json.dumps(dict(status='PASS',scope='Eleven deterministic normal and optimization-mode receipts agree byte for byte. This tests guard robustness; mathematical independence is specified in AUDIT.md.',pairs=pairs),indent=2)+'\n',encoding='utf8')
if __name__=='__main__':run()
