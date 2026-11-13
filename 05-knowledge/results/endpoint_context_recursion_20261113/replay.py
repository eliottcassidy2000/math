"""Reproduce this packet in either interpreter mode and compare exact receipts."""
from pathlib import Path
import subprocess,sys,json,gzip,hashlib,importlib
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SCRIPTS=['context_head','typed_glue','feedback_probe','reroot_probe','port_compiler','port_search','renewal_verify','renewal_catalog','language','collatz_graphs','walk_controls','audit']
RECEIPTS=SCRIPTS+['context_catalog']


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def seal():
    sys.path.insert(0,str(HERE))
    for s in SCRIPTS:importlib.import_module(s)
    paths={p for p in HERE.glob('*.py')}
    for module in list(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path:
            p=Path(path).resolve()
            if p.is_relative_to(ROOT) and p.suffix=='.py':paths.add(p)
    paths.update([HERE.parent/'multicore_metagraph_20261113/rooted_growth.json.gz',HERE.parent/'surplus_selector_20261113/transport_words.json.gz'])
    pins={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    (HERE/'REPLAY.json').write_text(json.dumps(dict(status='sealed inputs and executable dependencies',scripts=SCRIPTS,pins=pins),indent=2)+'\n',encoding='utf8')
    print('SEALED',len(pins),'files')


def main():
    if '--seal' in sys.argv:seal();return
    manifest=json.loads((HERE/'REPLAY.json').read_text(encoding='utf8'))
    for path,digest in manifest['pins'].items():need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'changed input '+path)
    optimized='--optimized' in sys.argv
    for s in SCRIPTS:
        cmd=[sys.executable]+(['-O'] if optimized else [])+['-X','utf8',str(HERE/(s+'.py'))]
        result=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
        if result.returncode:print(result.stdout);print(result.stderr);raise RuntimeError('failed '+s)
        print('PASS',s,flush=True)
    if optimized:
        for name in RECEIPTS:
            a=(HERE/(name+'.json.gz')).read_bytes();b=(HERE/(name+'.optimized.json.gz')).read_bytes()
            need(a==b,'normal/optimized receipt mismatch '+name)
        print('BYTE-IDENTICAL',len(RECEIPTS),'receipt pairs',flush=True)


if __name__=='__main__':main()
