#!/usr/bin/env python3
"""Independent root-Gram audit; require-based exact checks survive -O."""
import itertools, json, math
from fractions import Fraction as F
import numpy as np

def require(c,m):
    if not c: raise RuntimeError(m)
def roots_D(n):
    r=[]
    for i,j in itertools.combinations(range(n),2):
        for a,b in itertools.product([-2,2],repeat=2):
            x=[0]*n; x[i]=a;x[j]=b;r.append(x)
    return r
E=roots_D(8)+[list(x) for x in itertools.product([-1,1],repeat=8) if math.prod(x)==1]
roots=[np.array([x+[0]*8 for x in E]+[[0]*8+x for x in E],dtype=np.int64), np.array(roots_D(16),dtype=np.int64)]
res={}
for name,X,expected in zip(['E8square','D16plus'],roots,[[0,30720,1146880,132382720],[0,30720,1146880,49807360]]):
    require(X.shape==(480,16),'root shape')
    require(np.array_equal(X.T@X,240*np.eye(16,dtype=np.int64)),'isotropic second moment')
    # C=16 times Gram operator of traceless root dyadics.
    C=(X@X.T)**2-4
    require(np.array_equal(C,C.T),'symmetric')
    # Absolute entry bound controls every product/sum through degree4:
    # 480**3*60**4 < 2**63, so even naive products cannot overflow int64.
    require(480**3*60**4 < 2**63,'integer safety')
    require(np.max(abs(C))<=60,'entry bound')
    P=np.eye(480,dtype=np.int64); tm=[F(135)]
    for k in range(1,5):
        P=P@C;tm.append(F(sum(int(P[i,i]) for i in range(480)),16**k))
    c=F(40,3)
    ta=[sum(F(math.comb(k,j))*(-c)**(k-j)*tm[j] for j in range(k+1)) for k in range(1,5)]
    require(ta==expected,'independent root-Gram trace mismatch')
    res[name]={'M_power_traces':list(map(str,tm[1:])), 'A_power_traces':list(map(str,ta))}
# Sharp tracezero variance bounds, no square roots needed.
for d in range(2,15):
    for seed in range(1,35):
        v=[F(((seed+3)*(j+7)**2)%29-14) for j in range(d)]
        mean=sum(v)/d;v=[x-mean for x in v];a=-min(v);var=sum(x*x for x in v)
        require(a*a<=F(d-1,d)*var and a*a>=var/(d*(d-1)),'variance inequality')
    for v in [[F(-(d-1))]+[F(1)]*(d-1),[F(d-1)]+[F(-1)]*(d-1)]:
        a=-min(v);var=sum(x*x for x in v)
        require(a*a in [F(d-1,d)*var,var/(d*(d-1))],'sharpness')
# Symmetric Ramanujan evaluation with exact symbolic pi^2*E4=y.
import sympy as s
y=s.symbols('y',positive=True)
# F''/F=18+10*y/9; F'/F=-4; Delta/F=E4/1728, so b/F=y/1728.
r=s.factor(((18+s.Rational(10,9)*y)/144-s.Rational(4,16))*1728/y)
require(s.simplify(r-(s.Rational(40,3)-216/y))==0,'selfdual ratio')
# Exact tetrahedral harmonic directional evaluation using unnormalized edges:
tetra=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
total=F(0)
for a,b in itertools.combinations(tetra,2):
    v=[a[i]-b[i] for i in range(3)];norm=sum(t*t for t in v)
    total+=F(v[0]**4-6*v[0]**2*v[1]**2+v[1]**4,norm**2)
require(total==-1,'tetra direction sum')
res.update({'tetra_direction_sum':str(total),'selfdual_ratio':str(r),'variance_dimensions_tested':[2,14],'status':'PASS'})
print(json.dumps(res,indent=2))
