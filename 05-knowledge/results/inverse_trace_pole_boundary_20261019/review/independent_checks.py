#!/usr/bin/env python3
"""Independent exact audit of the trace-pole lemmas. No Keller realization claim.

Uses binomial sums of explicitly specified root packets and logarithmic
coefficients of characteristic polynomials, rather than the submitted
Newton/companion implementation. Includes ramified and repeated observables.
"""
from pathlib import Path
import hashlib, itertools, json, math
import sympy as S

HERE = Path(__file__).resolve().parent
DRAFT = HERE.parent
t, Z, X, Y, w, z, u, v = S.symbols('t Z X Y w z u v')
checks = 0

def ck(b, msg):
    global checks
    checks += 1
    if not bool(b):
        raise RuntimeError(msg)

def eq(a,b,msg):
    ck(S.cancel(a-b)==0,msg)

def val(f):
    f=S.cancel(f)
    if not f: return S.oo
    n,d=f.as_numer_denom()
    return min(mon[0] for mon,c in S.Poly(n,t).terms())-min(mon[0] for mon,c in S.Poly(d,t).terms())

# Roots of (Z-c)^m - C*t^(-a)*U have trace moments obtained by
# direct root-of-unity summation of the binomial expansion.
def packet_sum(m,a,C,c,U,k):
    return S.expand(m*sum(S.binomial(k,m*j)*c**(k-m*j)*(C*t**(-a)*U)**j for j in range(k//m+1)))

def packet_poly(m,a,C,c,U):
    return S.expand((Z-c)**m-C*t**(-a)*U)

# An independent polynomial path extracts coefficients of log(product(1-beta*s)).
def logarithmic_moments(poly, n):
    d=S.degree(poly,Z)
    coeffs=S.Poly(poly,Z).all_coeffs()
    A={i:coeffs[i] for i in range(1,len(coeffs)) if coeffs[i]!=0}
    # Truncated polynomial dictionaries in an auxiliary formal variable.
    power={0:S.Integer(1)}
    logc={}
    for j in range(1,n+1):
        nxt={}
        for i,ai in power.items():
            for k,ak in A.items():
                if i+k<=n: nxt[i+k]=nxt.get(i+k,0)+ai*ak
        power=nxt
        for k,c in power.items(): logc[k]=logc.get(k,0)+S.Rational((-1)**(j+1),j)*c
    return [None]+[S.cancel(-k*logc.get(k,0)) for k in range(1,n+1)]

# Exhaust finite leading-residue multisets independently of polynomial recurrences.
# Include repeats, signs, and imaginary residues, testing r rather than only d-1.
residue_cases=0
for r in range(1,7):
    for residues in itertools.combinations_with_replacement([S.Integer(1),-S.Integer(1),S.I,-S.I],r):
        sums=[S.expand(sum(c**k for c in residues)) for k in range(1,r+1)]
        ck(any(s != 0 for s in sums),'nonzero leading packet detected by r')
        if all(s == 0 for s in sums[:-1]):
            poly=S.expand(S.prod(Z-c for c in residues))
            ck(all(S.expand(poly).coeff(Z,k)==0 for k in range(1,r)),'maximal delay is Z^r-c')
            ck(len(set(residues))==r,'maximal delay has distinct residues')
        residue_cases+=1

cases=[
    ('one escaping root',[(1,3,2,1,1+t)],[0,1]),
    ('fractional cubic',[(3,2,1,t,1+t)],[2]),
    ('mixed slopes',[(2,3,1,1,1+t),(3,2,2,0,1+t*t)],[0]),
    ('merged quartic packet',[(2,1,1,0,1),(2,1,-1,0,1)],[0]),
    ('repeated nonprimitive observable',[(2,3,1,0,1),(2,3,1,0,1)],[0,0]),
    ('regular translated roots',[(1,0,2,t,1),(2,0,1,2,1+t)],[0]),
]
for m in range(1,8):
    for a in range(1,5):
        cases.append((f'extremal m={m}, a={a}',[(m,a,1,0,1+t)],[1]))
report=[]
for name,packets,finite in cases:
    m=sum(p[0] for p in packets if p[1]>0)
    degree=sum(p[0] for p in packets)+len(finite)
    lam=max([S.Rational(a,m0) for m0,a,C,c,U in packets if a>0] or [S.Integer(0)])
    r=sum(m0 for m0,a,C,c,U in packets if a>0 and S.Rational(a,m0)==lam)
    sums=[None]+[S.cancel(sum(packet_sum(*p,k) for p in packets)+sum(b**k for b in finite)) for k in range(1,degree+1)]
    values=[max(S.Integer(0),-S.sympify(val(sums[k]))/S.Integer(k)) for k in range(1,degree+1)]
    if m:
        ck(max(values[:r])==lam,name+' r-depth exact order')
        ck(max(values[:m])==lam,name+' m-depth exact order')
        ck(max(values[:degree-1])==lam,name+' d-1 exact order')
    else:
        ck(all(v==0 for v in values),name+' all bounded')
    poly=S.expand(S.prod(packet_poly(*p) for p in packets)*S.prod(Z-b for b in finite))
    alt=logarithmic_moments(poly,degree)
    for k in range(1,degree+1): eq(alt[k],sums[k],name+f' independent logarithm k={k}')
    first=next((k for k in range(1,degree+1) if val(sums[k])<0),None)
    report.append({'name':name,'degree':degree,'escaping_values':m,'most_singular_values':r,'largest_pole':str(lam),'first_polar_moment':first})

# Exact finite unramified projection-collision control: both roots are integral,
# but x projections collide and discriminant has positive even order.
for order in range(1,7):
    roots=[0,t**order,1+t,2+t*t]
    disc=S.prod((roots[i]-roots[j])**2 for i in range(4) for j in range(i+1,4))
    ck(val(disc)==2*order,'unramified projection collision even zero')

# Generic d with exactly one retained finite root: discriminant by resultant
# and its valuation, independent of the trinomial d-root hostile family.
for m in range(1,7):
    for a in range(1,4):
        poly=(Z-2)*((Z-t)**m-t**(-a)*(1+t))
        derivative=S.diff(poly,Z)
        resultant=S.resultant(poly,derivative,Z)
        ck(val(resultant)==-(m+1)*a,'retained-root discriminant pole')
        perm=[(j-a)%m for j in range(m)]
        visited=set(); lengths=[]
        for j in range(m):
            if j in visited: continue
            cursor=j; length=0
            while cursor not in visited:
                visited.add(cursor); length+=1; cursor=perm[cursor]
            lengths.append(length)
        g=math.gcd(a,m)
        ck(sorted(lengths)==[m//g]*g,'inertia cycle type')
        inversions=sum(perm[i]>perm[j] for i in range(m) for j in range(i+1,m))
        ck(inversions%2==((m+1)*a)%2,'direct permutation/discriminant parity')

# The no-survivor hostile costs one full moment.
for d in range(2,9):
    direct=[packet_sum(d,1,1,0,1,k) for k in range(1,d+1)]
    ck(all(s==0 for s in direct[:-1]),'no-survivor first d-1 blindness')
    eq(direct[-1],d/t,'no-survivor depth d')

# Four successive source boundary blowup charts, with volume order n-2.
for n in range(0,5):
    det=S.det(S.Matrix([1/z,z**n*w]).jacobian([z,w]))
    eq(det,-z**(n-2),f'boundary chart after {n} blowups')
eq(S.det(S.Matrix([z**3,-w/3]).jacobian([z,w])),-z*z,'triple different order two')
eq(S.det(S.Matrix([X**(-3),-X**4*Y/3]).jacobian([X,Y])),1,'rational model Jacobian one')
# Full rational f4 cotangent model in the same boundary chart, including units.
Q=z**3/(1+z**4)
P=w*(1+z**4)**2/(z**4-3)
eq(S.det(S.Matrix([Q,P]).jacobian([z,w])),-z*z,'f4 exact volume in boundary chart')
ck(S.limit(Q/z**3,z,0)==1,'f4 target parameter unit')
ck(S.limit(P,z,0)==-w/3,'f4 transverse coordinate limit')

source_names=('README.md','TRACE_POLE_BOUNDARY.md','SOURCE-DEPENDENCIES.json','verify_trace_pole_boundary.py','verification-results.json','verify.out','verify.optimized.out')
hashes={name:hashlib.sha256((DRAFT/name).read_bytes()).hexdigest() for name in source_names}
out={'status':'PASS','checks':checks,'leading_residue_multisets':residue_cases,'packet_cases':report,'scope':'Exact finite audit controls, not a realization or a general Jacobian-conjecture proof.','draft_sha256':hashes,'sympy':S.__version__}
(HERE/'independent-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(f'PASS: {checks} independent exact checks; {residue_cases} leading-residue multisets; {len(cases)} valued packet cases.')
print('Includes repeated observables, mixed slopes, rational pole orders, even projection collisions, inertia, and four blowup charts.')
print('No global Keller realization, unconditional trace regularity, or Jacobian-conjecture closure is inferred.')
