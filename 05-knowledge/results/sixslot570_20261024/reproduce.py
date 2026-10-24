#!/usr/bin/env python3
"""Reproduce exact consequences; compare math rather than wall-clock timings."""
from pathlib import Path
import json,subprocess,sys,tempfile,hashlib
HERE=Path(__file__).resolve().parent
def need(v,m):
    if not v:raise RuntimeError(m)
def main():
    if '--check' in sys.argv:
        M=json.loads((HERE/'MANIFEST.json').read_text())
        for p,h in M['files'].items():need(hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h,'file changed '+p)
    # Python normal/-O outputs must agree byte-for-byte.
    for script in ('verify','shelter'):
        a=subprocess.check_output([sys.executable,str(HERE/f'{script}.py')],cwd=HERE)
        b=subprocess.check_output([sys.executable,'-O',str(HERE/f'{script}.py')],cwd=HERE)
        need(a==b==(HERE/f'{script}.normal.json').read_bytes()==(HERE/f'{script}.optimized.json').read_bytes(),script+' replay mismatch')
        print(script,'PASS',flush=True)
    for script,filename in [('seek570','seek570.json'),('contexts570','contexts570.json'),('bridges','bridges.json')]:
        expected=(HERE/filename).read_bytes()
        for flags in ([],['-O']):
            subprocess.run([sys.executable,*flags,str(HERE/f'{script}.py')],cwd=HERE,check=True,stdout=subprocess.DEVNULL)
            need((HERE/filename).read_bytes()==expected,script+' deterministic replay mismatch')
        print(script,'PASS',flush=True)
    for name,start,finish in [('route22',12864,69272576),('route18_full28',22416516,134320219)]:
        for flags in ([],['-O']):
            raw=subprocess.check_output([sys.executable,*flags,str(HERE/'route.py'),'--start',str(start),'--finish',str(finish)],cwd=HERE)
            need(raw==(HERE/f'{name}.json').read_bytes(),'actual route replay '+name)
        print(name,'PASS',flush=True)
    # C++ implementations use different state generation and adjacency predicates.
    with tempfile.TemporaryDirectory(prefix='sixslot-replay-') as d:
        d=Path(d)
        for name,opt in [('fibres','-O3'),('audit_fibres','-O2')]:
            exe=d/name;subprocess.run(['clang++',opt,'-std=c++20',str(HERE/f'{name}.cpp'),'-o',str(exe)],check=True)
            for b in (6,8):
                args=[str(exe),str(HERE/('labels.txt' if name=='fibres' else 'audit_labels.txt'))]
                if name=='fibres':args+=['0' if b==6 else '4','14' if b==6 else '5']
                args+=[str(b)];raw=subprocess.check_output(args)
                new=[json.loads(s) for s in raw.splitlines()];old=[json.loads(s) for s in (HERE/f'{name}{b}.jsonl').read_text().splitlines()]
                need(len(new)==len(old),'census grade count')
                for a,z in zip(new,old):
                    a.pop('seconds',None);z.pop('seconds',None);need(a==z,'census exact replay '+name)
                print(name,b,'PASS',flush=True)
    print('all exact native/census checks PASS')
if __name__=='__main__':main()
