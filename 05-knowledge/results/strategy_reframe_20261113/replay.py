"""Pin inputs and reproduce all exact controls in normal and optimized Python."""
from pathlib import Path
import hashlib,importlib,json,subprocess,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SCRIPTS=['growth_probe','local_probe','transfer_probe']
RECEIPTS=['growth','local','transfer']

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def seal():
    sys.path.insert(0,str(HERE))
    for name in SCRIPTS:importlib.import_module(name)
    paths=set(HERE.glob('*.py'))
    for module in list(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path:
            p=Path(path).resolve()
            if p.is_relative_to(ROOT) and p.suffix=='.py':paths.add(p)
    paths.add(HERE.parent/'endpoint_context_recursion_20261113/renewal_catalog.json.gz')
    pins={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    (HERE/'REPLAY.json').write_text(json.dumps(dict(status='sealed executable dependencies and input catalog',pins=pins),indent=2)+'\n',encoding='utf-8')
    print('SEALED',len(pins),flush=True)

def run():
    if '--seal' in sys.argv:seal();return
    pins=json.loads((HERE/'REPLAY.json').read_text(encoding='utf-8'))['pins']
    for path,digest in pins.items():need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'changed input '+path)
    for mode in ([],['-O']):
        for name in SCRIPTS:
            result=subprocess.run([sys.executable,*mode,'-X','utf8',str(HERE/(name+'.py'))],cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
            if result.returncode:print(result.stdout);print(result.stderr);raise RuntimeError(name)
            print('PASS',name,'optimized' if mode else 'normal',flush=True)
    for name in RECEIPTS:need((HERE/(name+'.json')).read_bytes()==(HERE/(name+'.optimized.json')).read_bytes(),'mode mismatch '+name)
    print('BYTE-IDENTICAL',len(RECEIPTS),'receipt pairs',flush=True)

if __name__=='__main__':run()
