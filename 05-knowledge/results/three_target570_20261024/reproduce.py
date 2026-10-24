"""Reproduce exact certificates; no silent regeneration of expected results."""
from pathlib import Path
import subprocess,sys,json,hashlib
HERE=Path(__file__).resolve().parent
def main():
    manifest=json.loads((HERE/'MANIFEST.json').read_text())
    for p,h in manifest.items():
        if hashlib.sha256((HERE/p).read_bytes()).hexdigest()!=h:raise RuntimeError('manifest mismatch '+p)
    runs=[('menus',False,'menus.json'),('bank',False,'bank.json')]
    runs += [(p,opt,p+('.optimized' if opt else '.normal')+'.json')
             for p in ('check','guarded','weighted','twelve') for opt in (False,True)]
    for name,opt,expected in runs:
        cmd=[sys.executable]+(['-O'] if opt else [])+['-B',str(HERE/(name+'.py'))]
        result=json.loads(subprocess.check_output(cmd,cwd=HERE))
        if result!=json.loads((HERE/expected).read_text()):raise RuntimeError('reproduction mismatch '+expected)
        print('PASS',expected,flush=True)
    print('PASS manifest and all exact replays',flush=True)
if __name__=='__main__':main()
