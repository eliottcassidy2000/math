import importlib.util
from pathlib import Path

OLD=Path(__file__).resolve().parents[2]/'trijection_hourglass_ports_20261025/trijection/explore.py'
spec=importlib.util.spec_from_file_location('old',OLD);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
Parity=old.Parity;edges=old.edges;edge=old.edge

def source(n,h,guards=()):
    d=Parity(n);p=tuple(range(n));d.anti(p)
    for i in range(n-1):d.force(i,i+1)
    for v in range(1,n):d.force(v,0,int(v in {h,n-3,n-2,n-1}))
    for u,v in guards:d.force(u,v)
    return d

def status(n,h,words,guards=()):
    base=source(n,h,guards)
    if not base.ok:return None
    flips=set();blue=[];illegal=[]
    for q in words:
        for u,v in zip(q,q[1:]):
            d=source(n,h,guards)
            d.force(u,v,int(edge(u,v) in flips))
            if d.ok:illegal.append((u,v))
        d=source(n,h,guards);d.anti(q,flips);blue.append(d.ok)
        flips ^= set(base.pairs)^edges(q)
    return illegal,blue,flips==(set(base.pairs)^edges(tuple(range(n))))

def single(n,h):
    d=n-1-h
    prefix=tuple(range(1,d))+tuple(range(d+1,h))
    end=tuple(range(h,n))
    return ((0,)+prefix+end+(d,),prefix+(0,)+end+(d,),tuple(range(1,h))+(0,)+end)

def triple(n,h):
    prefix=(1,2)+tuple(range(6,h));end=tuple(range(h,n))
    return ((0,)+prefix+end+(4,5,3),prefix+(0,)+end+(4,5,3),tuple(range(1,h))+(0,)+end)

if __name__=='__main__':
    for n in range(11,46):
        h=n-5
        print(n,'single',status(n,h,single(n,h),[(2,4),(3,5)]),
              'triple',status(n,h,triple(n,h),[(2,4),(5,3),(2,6)]))
