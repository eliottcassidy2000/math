#!/usr/bin/env python3
from independent_carrier import packet
import json

def ck(p,s):
    if not p:raise RuntimeError(s)
def base_count(r):
    bits=tuple(1-((r>>i)&1)for i in range(r.bit_length()))
    if not bits:return 1
    return int(bits[-1]==0)+int(bits[0]==1)+sum(bits[i]==0 and bits[i+1]==1 for i in range(len(bits)-1))+sum(2**(i-j-1)for i in range(len(bits))for j in range(i)if bits[j]==1 and bits[i]==0)
def defect(r,a):
    k=(r>>a)&3;t=r>>(a+2);p=1 if a==0 else(r>>(a-1))&1
    f=3*t+2*p-(t%2)
    return 0 if k in(0,3)else f if k==1 else-f-2

def triangle_core(r,a):
    pre=r%2**a;k=(r>>a)&3;suf=r>>(a+2);b=suf.bit_length();n=a+3+b
    D=[sum(1<<j for j in range(i+1,n))for i in range(n)]
    # replace transitive triple by0->1->2->0
    D[a]&=~(1<<(a+2));D[a+2]|=1<<a
    inc=pre|(((1<<k)-1)<<a)|(suf<<(a+3))
    T=[sum(1<<(i+1)for i in range(n)if not(inc>>i&1))]
    T +=[(row<<1)|((inc>>i)&1)for i,row in enumerate(D)]
    return tuple(T)
literal=0
for r in range(64):
    B=base_count(r)
    for a in range(r.bit_length()+2):
        p=packet(triangle_core(r,a))
        ck(p['H']==3*B+defect(r,a),'literal binary dimer')
        ck(p['h']==3 and p['V']==3*(r+1),'literal receiver')
        literal+=1
for r in range(1<<18):
    ds={defect(r,a)for a in range(r.bit_length()+2)}
    ck(not any(d+2 in ds and d+4 in ds for d in ds),'no3 control')
print(json.dumps(dict(status='PASS',literal_dimer_controls=literal,all_binary_defect_fibers=1<<18),indent=2,sort_keys=True))
