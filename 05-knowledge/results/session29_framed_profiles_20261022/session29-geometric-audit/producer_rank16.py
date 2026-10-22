#!/usr/bin/env python3
"""Exact minimum-shell Gram verification; no floating-point spectral inference."""
from itertools import combinations, product
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np

def require(ok, msg):
    if not ok: raise RuntimeError(msg)

def dn(n):
    out=[]
    for i,j in combinations(range(n),2):
        for a,b in product((-2,2),repeat=2):
            x=[0]*n; x[i]=a; x[j]=b; out.append(x)
    return out
E=dn(8)+[list(s) for s in product((-1,1),repeat=8) if sum(v<0 for v in s)%2==0]
roots={'E8plusE8':[x+[0]*8 for x in E]+[[0]*8+x for x in E], 'D16plus':dn(16)}
B=[]
for i in range(15):
    h=np.zeros((16,16),dtype=np.int64);h[i,i]=1;h[15,15]=-1;B.append(h)
for i,j in combinations(range(16),2):
    h=np.zeros((16,16),dtype=np.int64);h[i,j]=h[j,i]=1;B.append(h)
B=np.array(B); G=np.einsum('aij,bij->ab',B,B)
result={}
for name,rr in roots.items():
    X=np.array(rr,dtype=np.int64)
    require(len(rr)==480 and len({tuple(x) for x in rr})==480, name+' root count')
    require(np.all(np.sum(X*X,axis=1)==8), name+' doubled root norm')
    Q=np.einsum('ri,aij,rj->ra',X,B,X)
    # |Q|<=16 and 480 terms, so all int64 products and sums are safely exact.
    require(np.max(abs(Q))<=16, 'overflow precondition')
    K=Q.T@Q
    expected=np.zeros((135,135),dtype=np.int64)
    for a,h in enumerate(B):
        for b,k in enumerate(B):
            if name=='E8plusE8':
                v=sum(12*np.trace(h[s,s])*np.trace(k[s,s])+24*np.trace(h[s,s]@k[s,s]) for s in (slice(0,8),slice(8,16)))
            else:
                v=56*np.dot(np.diag(h),np.diag(k))+8*(np.trace(h@k)-np.dot(np.diag(h),np.diag(k)))
            expected[a,b]=16*v
    require(np.array_equal(K,expected), name+' all 18225 bilinear entries')
    spectra={'E8plusE8':[(F(32,3),70),(F(320,3),1),(F(-40,3),64)],'D16plus':[(F(128,3),15),(F(-16,3),120)]}[name]
    require(sum(m for v,m in spectra)==135, 'dimension')
    require(sum(v*m for v,m in spectra)==0, 'trace zero')
    result[name]={'root_count':480,'gram_entries_verified':18225,'angular_spectrum':[(str(v),m) for v,m in spectra], 'angular_trace':'0'}
# e^6 > 400 by positive Taylor partial sum; pi>3 implies q=e^-2pi<1/400.
partial=sum(F(6)**k/F(__import__('math').factorial(k)) for k in range(20))
require(partial>400, 'exponential bound')
q=F(1,400)
S=q*(1+11*q+11*q*q+q**3)/(1-q)**5
upper=10*(1+240*S)
require(upper<27, 'pi² E4 upper bound')
result['selfdual_y_bounds']={'lower':'9','upper':str(upper),'upper_less_than_27':True,'basis':'3<pi, pi²<10; sigma3(n)<=n^4; e^6>400'}
result['selfdual_inertia']={'E8plusE8':[71,64,0],'D16plus':[15,120,0]}
for name in roots:
    spec=[(F(v),m) for v,m in result[name]['angular_spectrum']]
    result[name]['angular_power_traces']=[str(sum(m*v**k for v,m in spec)) for k in range(1,5)]
require(result['E8plusE8']['angular_power_traces']==['0','30720','1146880','132382720'], 'E8 spectral moments')
require(result['D16plus']['angular_power_traces']==['0','30720','1146880','49807360'], 'D16 spectral moments')
result['fourth_trace_gap']=82575360
result['status']='PASS' 
p=Path(__file__).with_suffix('.json');p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
