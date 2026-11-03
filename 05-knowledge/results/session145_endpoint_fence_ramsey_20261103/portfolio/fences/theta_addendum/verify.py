#!/usr/bin/env python3
"""Independent polynomial identities for the scoped theta optimum."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys
checks=0
def require(ok,label):
    global checks
    checks+=1
    if not ok: raise RuntimeError(label)
def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a): c[i]+=x
    for i,x in enumerate(b): c[i]+=x
    return c
def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def scale(a,s): return [s*x for x in a]
def derivative(a): return [i*a[i] for i in range(1,len(a))]
def qmul(a,b): return (a[0]*b[0]+13*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def qadd(a,b): return (a[0]+b[0],a[1]+b[1])
def qscale(a,c): return (a[0]*c,a[1]*c)
def qeval(poly,z):
    out=(F(0),F(0))
    for c in reversed(poly): out=qadd(qmul(out,z),(c,F(0)))
    return out
L=[F(5),F(-20),F(16)]; C=[F(1),F(-12),F(16)]
B=add([F(5)],scale(mul(L,C),-1))
require(B==[F(0),F(80),F(-336),F(512),F(-256)],'area polynomial')
identity=add(add(mul([0,2,-2],derivative(B)),scale(B,-1)),mul([0,4,-4],mul(C,derivative(L))))
require(all(x==0 for x in identity),'derivative identity all coefficients')
wrong=add(add(mul([0,2,-2],derivative(B)),scale(B,-1)),mul([0,-4,4],mul(C,derivative(L))))
require(any(x!=0 for x in wrong),'wrong derivative sign rejected')
z=(F(5,8),F(-1,8))
require(qeval(L,z)==(F(2),F(0)),'terminal distance two')
require(qeval([3,-20,16],z)==(F(0),F(0)),'quadratic radius root')
# In Q(sqrt13), inverse z = (10/3)+(2/3)sqrt13.
iz=(F(10,3),F(2,3))
require(qmul(z,iz)==(F(1),F(0)),'inverse exact')
ratio=qadd(iz,(F(-1),F(0)))
bz=qeval(B,z); area2=qscale(qmul(ratio,qmul(bz,bz)),F(1,16))
require(area2==(F(89,16),F(13,8)),'cyclic-cap squared area')
# Check the branch and a host outside the positive-derivative domain.
require(F(18,5)**2<13<F(37,10)**2,'rational radical bracket')
require(0<(F(5)-F(37,10))/8<(F(5)-F(18,5))/8<1,'admissible minus root')
require((F(5)+F(18,5))/8>1,'plus root impossible as sine squared')
for zq in [F(7,40)+F(j,1000) for j in range(166)]:
    lq=sum(c*zq**i for i,c in enumerate(L)); cq=sum(c*zq**i for i,c in enumerate(C))
    require(0<lq<2 and cq<0 and 32*zq-20<0,'interior monotonicity probes')
require(sum(c*F(1,100)**i for i,c in enumerate(C))>0,'derivative sign reverses outside scope')
p=Path(__file__); proof=p.with_name('PROOF.md')
result={'status':'PASS','checked_conditions':checks,'derivative_identity_coefficients':[str(x) for x in identity],
 'squared_cap_area':{'rational':str(area2[0]),'sqrt13':str(area2[1])},
 'scope':'theta(5,5,2), shared two-edge separator, straight terminal chord inside outer polygon',
 'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'proof_sha256':hashlib.sha256(proof.read_bytes()).hexdigest()}
out=Path(sys.argv[1] if len(sys.argv)>1 else 'verification.json'); out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','checked_conditions':checks,'receipt':str(out)}))
