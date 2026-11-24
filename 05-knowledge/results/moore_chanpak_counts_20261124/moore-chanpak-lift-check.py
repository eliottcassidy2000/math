from pathlib import Path
from runpy import run_path
from math import isqrt
F=run_path(str(Path(__file__).with_name('moore-chanpak-fork.py')))
rows=F['rows']; flip=F['flip']; heights=F['heights']; pars=F['tail_params']; data=F['tail_data']

def apex(r,w):
 return (*[row|((1-((w>>i)&1))<<len(r)) for i,row in enumerate(r)],w)

def combine(r,A):
 x=len(r);m=len(A)
 return apex(tuple(e|(((1<<m)-1)<<x) for e in r)+tuple(e<<x for e in A),1240)

def check(A):
 hx,f=heights(A)
 h=[1]+[sum(w) for w in f[1:]]
 q=sum(a*h[-1-i] for i,a in enumerate(h))
 out=tuple(heights(combine(r,A))[0] for r in (r0,r1))
 assert out==(245121*q-74075*hx,245121*q-74073*hx)
 print('H,Q=',hx,q,'heights',out)

base=rows((0,1,4,3,2));r0=apex(base,194);r1=apex(flip(base,3,4),194)
assert [pars(data(r),1240) for r in (r0,r1)]==[(245121,-74075),(245121,-74073)]
for k in (1,2,3):
 check(tuple(((1<<k)-1)^((1<<(i+1))-1) for i in range(k)))
check((2,4,1))
for k in (1,6):
 p=245121*2**k-74075
 assert all(p%d for d in range(2,isqrt(p)+1))
print('trial-divided prime points k=1,6')
