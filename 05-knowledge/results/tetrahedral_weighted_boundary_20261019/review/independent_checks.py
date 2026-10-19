#!/usr/bin/env python3
"""Independent checks for the tetrahedral weighted-boundary audit.
No imports from the submitted verifier. Does not modify the source packet.
"""
from pathlib import Path
from itertools import permutations, product
from collections import Counter
import hashlib,json
import sympy as s

checks=[]
def test(name, condition):
    if not condition:
        raise RuntimeError(name)
    checks.append(name)
def zero(name, expression):
    test(name,s.cancel(expression)==0)

h,e,a,b,c,u,v,t,z=s.symbols('h e a b c u v t z')
U=9*h**4*e**2-8*h*e
V=-27*h**6*e**3+36*h**3*e**2-8*e
zero('polynomial Jacobian',s.det(s.Matrix([U,V]).jacobian([h,e]))-64*e)
f=27*a**4-18*u*a**2-8*v*a-u**2
uX=a**2+8*b; vX=a**3-20*a*b-8*c; eX=a**3-2*a*b+c
zero('monogenic differential exactly 64 e',s.diff(f,a).subs({u:uX,v:vX},simultaneous=True)-64*eX)
zero('cone relation from quartic',f.subs({u:uX,v:vX},simultaneous=True)-64*(a*c-b*b))
zero('cusp factorization',f.subs({u:z*z,v:z**3},simultaneous=True)-(a-z)*(3*a+z)**3)
zero('cone chart e valuation',eX.subs({b:a*t,c:a*t*t},simultaneous=True)-a*(a-t)**2)
zero('cone chart h valuation',(a*a-b).subs(b,a*t)-a*(a-t))
zero('cone chart delta valuation',(uX**3-vX**2).subs({b:a*t,c:a*t*t},simultaneous=True)-64*a*a*t*(a-t)**3)
zero('L0 chart simple b-a squared zero',(b-a*a).subs(a,b*b/c)-b*(1-b**3/c**2))
zero('L0 chart a has order two',a.subs(a,b*b/c)-b*b/c)
zero('L0 chart e nonzero at generic L0',eX.subs({a:0,b:0},simultaneous=True)-c)
zero('inverse h relation',27*(U**3-V**2)*h**4+64*V*h-64*U)
zero('inverse discriminant pullback',U**3-V**2-64*e*e*(h**3*e-1))
# Independent matrix-trace calculation in the quartic algebra K[h].
D=u**3-v**2; q=64*v/(27*D); r=-64*u/(27*D)
C=s.Matrix([[0,0,0,-r],[1,0,0,-q],[0,1,0,0],[0,0,1,0]])
for n,wanted in enumerate([0,0,-64*v/(9*D),256*u/(27*D)],1):
    zero('companion matrix trace '+str(n),s.trace(C**n)-wanted)
# A second route to the cyclic cubic resolvent.
p=-2*u/s.Integer(3); q=-8*v/s.Integer(27); r=-u*u/s.Integer(27)
res=t**3-p*t**2-4*r*t+4*p*r-q*q
zero('classical pair resolvent is depressed cubic',res-(t+2*u/9)**3-64*(u**3-v**2)/729)
zero('quartic discriminant',s.discriminant(f,a)+80621568*D**2)
# Singular loci checked through the defining derivatives.
cone=b*b-a*c
cone_grad=s.groebner([s.diff(cone,x) for x in (a,b,c)],a,b,c)
test('cone singular ideal is origin',tuple(g.as_expr() for g in cone_grad.polys)==(a,b,c))
A,B,C0=s.symbols('A B C0'); d4=C0*C0-B**3+64*A**3
crit=s.groebner([s.diff(d4,x) for x in (A,B,C0)],A,B,C0)
test('cubic resolvent singular locus is only origin',tuple(g.as_expr() for g in crit.polys)==(A*A,B*B,C0))
# Independent permutation enumeration: all order-three subgroups are conjugate.
identity=tuple(range(4))
def mul(x,y):return tuple(x[y[i]] for i in range(4))
def inv(x):return tuple(x.index(i) for i in range(4))
def parity(p):return sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2
G={p for p in permutations(range(4)) if parity(p)==0}
orders={}
for g in G:
    cur=identity
    for n in range(1,13):
        cur=mul(cur,g)
        if cur==identity:orders[g]=n;break
test('A4 element orders independently enumerated',Counter(orders.values())==Counter({1:1,2:3,3:8}))
C3={frozenset((identity,g,mul(g,g))) for g in G if orders[g]==3}
test('four order-three subgroups',len(C3)==4)
H=next(iter(C3)); conj={frozenset(mul(mul(g,x),inv(g)) for x in H) for g in G}
test('all degree-four stabilizers conjugate',conj==C3)
# Every surjection from abstract Q8 to V4 kills center, and choosing
# images of i,j independent gives precisely six maps with this kernel.
V4=list(product(range(2),repeat=2))
xor=lambda x,y:(x[0]^y[0],x[1]^y[1])
linear_maps=[]
# Q8 = { +/- 1, +/- i, +/- j, +/- k }; quotient classes are 0,i,j,i+j.
qclasses=[(0,0),(0,0),(1,0),(1,0),(0,1),(0,1),(1,1),(1,1)]
for i,j in product(V4,repeat=2):
    images=[xor(tuple(k*x for x in i),tuple(l*x for x in j)) for k,l in qclasses]
    if len(set(images))==4:linear_maps.append(images)
test('six Q8 to V4 epimorphisms',len(linear_maps)==6)
test('unique epimorphism kernel',all([k for k,x in enumerate(m) if x==(0,0)]==[0,1] for m in linear_maps))
# Reproducibility and source identities.
root=Path(__file__).resolve().parent; packet=root.parent
source_names=('A4_WEIGHTED_BOUNDARY.md','README.md','SHA256SUMS','SOURCE-DEPENDENCIES.json','verification-results.json','verify.optimized.out','verify.out','verify_a4_weighted_boundary.py')
source_hashes={name:hashlib.sha256((packet/name).read_bytes()).hexdigest() for name in source_names}
ledger=json.loads((packet/'SOURCE-DEPENDENCIES.json').read_text())
for dep in ledger['corpus']:
    actual=hashlib.sha256((packet.parents[2]/dep['path']).read_bytes()).hexdigest()
    test('dependency hash '+dep['id'],actual==dep['sha256'])
out={'status':'PASS','checks':len(checks),'check_names':checks,'source_hashes':source_hashes,'scope':'Independent finite algebra, local valuation formulas, subgroup census, and dependency integrity. Geometry is separately audited in REVIEW.md.'}
(root/'independent-results.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(f'PASS: {len(checks)} independent checks')
for name in checks:print('  '+name)
