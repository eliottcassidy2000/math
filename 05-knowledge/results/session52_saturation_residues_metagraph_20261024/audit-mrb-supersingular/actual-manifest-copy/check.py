#!/usr/bin/env python3
"""Exact finite tests and directed-interval MRB certificate. No network.
All mathematical tests remain active under python -O. Requires mpmath.
"""
from collections import Counter
from fractions import Fraction
from math import comb, gcd, isqrt
from pathlib import Path
import json
from mpmath import iv

ROOT=Path(__file__).resolve().parent
CHECKS=0

def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok: raise RuntimeError(message)

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def legendre(a,p):
    a %= p
    if not a: return 0
    b=pow(a,(p-1)//2,p)
    require(b in (1,p-1),'Euler criterion not a sign')
    return 1 if b==1 else -1

def forms(D):
    require(D<0 and D%4 in (0,1),'invalid form discriminant')
    ans=[]
    for a in range(1,isqrt((-D)//3)+1):
        for b in range(-a,a+1):
            if (b*b-D)%(4*a): continue
            c=(b*b-D)//(4*a)
            if a>c or gcd(gcd(a,b),c)>1: continue
            if (abs(b)==a or a==c) and b<0: continue
            ans.append([a,b,c])
    return ans

def rational_ss(p):
    """One model per j; count squares by direct residue membership, not Euler."""
    squares=Counter(y*y%p for y in range(p))
    js=[]
    for j in range(p):
        if j==0: a,b=0,1
        elif j==1728%p: a,b=1,0
        else:
            k=j*pow((1728-j)%p,-1,p)%p
            a,b=3*k%p,2*k%p
        require((4*a*a*a+27*b*b)%p!=0,'singular j model')
        count=1+sum(squares[(x*x*x+a*x+b)%p] for x in range(p))
        if count==p+1: js.append(j)
    return js

O=[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]
units=[1,2,4,7,8,11,13,14]
primes=[p for p in range(2,1000) if prime(p)]
rows=[]
for r in units:
    positives=[p for p in O if p%15==r]
    negative=next(p for p in primes if p%15==r and p not in O)
    ordinary=next(p for p in primes if p>5 and p%15==r and p%8 in (1,3))
    ss=next(p for p in primes if p>5 and p%15==r and p%8 in (5,7))
    require(positives and prime(negative),'mod15 witness missing')
    rows.append({'r':r,'ogg':positives,'non_ogg':negative,
                 'ordinary_E_prime':ordinary,'supersingular_E_prime':ss})
require(sorted({p%15 for p in O if gcd(p,15)==1})==units,'unit coverage')
require(Counter(p%15 for p in O)==Counter({1:1,2:3,3:1,4:1,5:1,7:1,8:1,11:3,13:1,14:2}),'Ogg residue multiplicities')

crt=[]
for r in units:
    for s in (1,3,5,7):
        xs=[x for x in range(120) if x%15==r and x%8==s]
        require(len(xs)==1,'CRT uniqueness')
        x=xs[0]
        p=next(p for p in primes if p%120==x)
        require(gcd(x,120)==1 and prime(p),'CRT prime witness')
        crt.append({'mod15':r,'mod8':s,'mod120':x,'prime_witness':p,
                    'chi_minus8':1 if s in (1,3) else -1})
require(len({v['mod120'] for v in crt})==32,'CRT complete')

# Q(sqrt(2)) exact elementary calculation, represented by rational pairs.
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def mul(x,y): return (x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(n,x): return (n*x[0],n*x[1])
def power(x,n):
    z=(1,0)
    for _ in range(n): z=mul(z,x)
    return z
lam=(3,-2)
require(add(add(power(lam,2),scale(-6,lam)),(1,0))==(0,0),'lambda polynomial')
u=add((1,0),scale(-1,lam))
num=scale(256,power(add(u,power(lam,2)),3))
den=mul(power(lam,2),power(u,2))
require(num==scale(8000,den),'j(lambda)=8000')
b2,b4,b6,b8=16,4,0,-4
c4=b2*b2-24*b4
Delta=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
require((c4,Delta)==(160,512) and c4**3==8000*Delta,'E invariants')

Erows=[]
for p in primes:
    if p==2: continue
    trace=-sum(legendre(x*x*x+4*x*x+2*x,p) for x in range(p))
    minus_trace=-sum(legendre(x*x*x-4*x*x+2*x,p) for x in range(p))
    require(minus_trace==legendre(-1,p)*trace,'inherited Bolza twist trace relation')
    character=legendre(-2,p)
    require(character==(1 if p%8 in (1,3) else -1),'CM character residue table')
    require((trace==0)==(character==-1),'CM finite point-count law')
    # A second independent point-count route for all small p and key hostiles.
    if p<=107:
        sq=Counter(y*y%p for y in range(p))
        count=1+sum(sq[(x*x*x+4*x*x+2*x)%p] for x in range(p))
        require(count==p+1-trace,'independent E count')
    Erows.append({'p':p,'trace':trace,'points':p+1-trace,'j_mod_p':8000%p,
                  'supersingular':trace==0,'ogg':p in O,'Bolza_minus_trace':minus_trace})
by_p={x['p']:x for x in Erows}
require(by_p[17]['points']==24 and not by_p[17]['supersingular'] and by_p[17]['ogg'],'17 hostile')
require(by_p[37]['points']==38 and by_p[37]['supersingular'] and not by_p[37]['ogg'],'37 hostile')
require(by_p[17]['p']%15==by_p[47]['p']%15 and by_p[17]['supersingular']!=by_p[47]['supersingular'],'17/47 same residue opposite CM')

ogg_cert=[]
for p in primes:
    if not 5<=p<=211: continue
    h4=forms(-4*p)
    h1=forms(-p) if p%4==3 else []
    degree=p//12+{1:0,5:1,7:1,11:2}[p%12]
    numerator=2*degree-len(h4)-len(h1)
    require(numerator>=0 and numerator%4==0,'nonrational pair count')
    require((numerator==0)==(p in O),'finite Ogg class number census')
    row={'p':p,'degree':degree,'forms_minus4p':h4,'forms_minusp':h1,
         'rational_count':(len(h4)+len(h1))//2,'g_plus':numerator//4}
    if p<=107:
        js=rational_ss(p)
        require(2*len(js)==len(h4)+len(h1),'independent class number/point count')
        row['rational_supersingular_j']=js
    ogg_cert.append(row)

# The shared two-sign wheel is exactly one screen, not two independent ones.
P=30030
surv_plus=[d for d in range(P) if gcd(P+d,P)==1]
surv_minus=[d for d in range(P) if gcd(P-d,P)==1]
require(surv_plus==surv_minus and len(surv_plus)==5760,'shared primorial wheel')
require(Fraction(len(surv_plus),P)==Fraction(192,1001),'wheel fraction')
require(prime(P-19) and P+19==151*199,'two-sign counterexample')
require(P-23==37*811 and P+23==41*733,'prime offset both endpoints composite')

# Exact formal generating function and finite identity, all coefficient arithmetic integral.
def mask(n): return int(gcd(n,15)==1)
coeff={n:(-1)**n*mask(n) for n in range(1,601)}
for n,c in coeff.items():
    inclusion=sum(mu*((-1)**(n//d)) for d,mu in [(1,1),(3,-1),(5,-1),(15,1)] if n%d==0)
    require(c==inclusion,'signed divisor inclusion-exclusion')
    expected=(-1)**n*mask(n) if n<=15 else 0
    require(c+(coeff.get(n-15,0))==expected,'(1+z^15)R=P formal coefficients')
    if n<=585: require(coeff[n+15]==-c,'parity flips under translation by 15')
require(sum(coeff[n] for n in range(1,31))==0,'mean-zero signed mask')
require(sum(mask(n) for n in range(1,16))==8,'unsigned mask density')
for period in (1,2,3,5,6,10,15):
    require(any((-1)**(n+period)*mask(n+period)!=(-1)**n*mask(n) for n in range(1,31)),f'false period {period}')
# Hostile control: parity-forgetting is not a period-15 factor of this signed observable.
require(mask(1)==mask(16) and (-1)**1*mask(1)!=(-1)**16*mask(16),'parity loss witness')
# Integer-valued artificial profile tests regrouping without transcendental rounding.
for N in (1,2,3,17,100,500):
    f=lambda n:n*n+3*n+7
    lhs=sum((-1)**n*mask(n)*f(n) for n in range(1,N+1))
    rhs=sum(mu*sum((-1)**k*f(d*k) for k in range(1,N//d+1)) for d,mu in [(1,1),(3,-1),(5,-1),(15,1)])
    require(lhs==rhs,'finite weighted inclusion-exclusion')

# Audit-supplied full-profile Fourier projection, certified in Z[z]/Phi_15.
ramanujan=[sum(d*mu for d,mu in [(1,1),(3,-1),(5,-1),(15,1)] if j%d==0) for j in range(15)]
require(ramanujan==[8,1,1,-2,1,-4,-2,1,1,-2,-4,1,-2,1,1],'Ramanujan15 vector')
phi15=[1,-1,0,1,-1,1,0,-1,1]
for n in range(1,121):
    poly=[0]*15
    for j,c in enumerate(ramanujan): poly[j*n%15]+=c
    for degree in range(14,7,-1):
        lead=poly[degree]
        for k,c in enumerate(phi15): poly[degree-8+k]-=lead*c
    require(poly==[15*mask(n)]+[0]*14,'exact Fourier unit projection')

# Theta parity is coefficientwise, not a numerical fit.
for exponent in range(1,1025):
    n=isqrt(exponent)
    theta=(2+2*(-1)**n) if n*n==exponent else 0
    k=isqrt(exponent//4)
    target=4 if exponent==4*k*k and k>=1 else 0
    require(theta==target,'theta parity coefficient')

# Recomputed interval certificate with analytic Cauchy/Euler tail (Section 6).
iv.dps=130
M,N=256,100
f=lambda n:iv.exp(iv.log(iv.mpf(n))/n)-1
A=lambda n:2/iv.pi*iv.exp(-iv.log(iv.mpf(n))/n)*iv.sin(iv.pi/n)
weights=[sum(comb(j,k)*2**(N-j-1) for j in range(k,N)) for k in range(N)]
tail=sum(((-1)**k*weights[k]*f(M+k) for k in range(N)),iv.mpf(0))/2**N
prefix=sum(((-1)**n*f(n) for n in range(1,M)),iv.mpf(0))
L=iv.log(iv.mpf(3)*M/2)+iv.pi/6
B=2*iv.factorial(N)*iv.exp(2*L/M)*(L/iv.mpf(M)**(N+1)+(L/N+iv.mpf(1)/N**2)/iv.mpf(M)**N)
C=prefix+(-1)**M*tail+iv.mpf([-1,1])*B
D=A(8)-C
q4=iv.exp(-4*iv.pi*iv.sqrt(2))
require(bool(B.b<iv.mpf('2.654e-84')),'analytic remainder')
require(bool(D.a>iv.mpf('1.891807127880073527054533938441983110587504636897727612419589120e-8')),'defect lower bound')
require(bool(D.b<iv.mpf('1.891807127880073527054533938441983110587504636897727612419589122e-8')),'defect upper bound')
require(bool(D.a>0) and bool((q4-D).a>0),'0<D<q^4')
require(bool((D/q4).a>iv.mpf('.9884507438957745199776876885831501151094671392875')) and bool((D/q4).b<iv.mpf('.9884507438957745199776876885831501151094671392876')),'defect ratio')
nearest=[]
for n,r in [(8,C),(105719004,D)]:
    low=r-(A(n)+A(n+1))/2
    high=(A(n-1)+A(n))/2-r
    require(bool(low.a>0) and bool(high.a>0),'nearest amplitude midpoint')
    nearest.append({'index':n,'lower_margin':str(low),'upper_margin':str(high)})

# The same proof for g_d(z)=f(dz) carries an additional factor 1/d.
scaled={1:C}
for d in (3,5,15):
    values=[f(d*n) for n in range(1,M+N)]
    pre=sum(((-1)**n*values[n-1] for n in range(1,M)),iv.mpf(0))
    smooth=sum(((-1)**k*weights[k]*values[M+k-1] for k in range(N)),iv.mpf(0))/2**N
    Ld=iv.log(iv.mpf(3)*d*M/2)+iv.pi/6
    Bd=2*iv.factorial(N)*iv.exp(2*Ld/(d*M))/d*(Ld/iv.mpf(M)**(N+1)+(Ld/N+iv.mpf(1)/N**2)/iv.mpf(M)**N)
    scaled[d]=pre+(-1)**M*smooth+iv.mpf([-1,1])*Bd
unit_C=scaled[1]-scaled[3]-scaled[5]+scaled[15]
nonunit_C=scaled[3]+scaled[5]-scaled[15]
require(bool(unit_C.a>iv.mpf('.53758442393739534102')) and bool(unit_C.b<iv.mpf('.53758442393739534103')),'unit-filter MRB enclosure')
require(bool((unit_C-iv.mpf(8)*C/15).a>iv.mpf('.437')),'density is not an MRB scalar weight')
require(bool(nonunit_C.b<0),'nonunit signed sum is negative')

result={'schema':1,'status':'PASS','checks':CHECKS,
        'scope':{'prime_range':'2 <= p < 1000','CM_point_counts':len(Erows),
                 'Ogg_class_number_range':'5 <= p <= 211','all_j_point_count_range':'5 <= p <= 107',
                 'formal_coefficient_range':'1 <= n <= 600','theta_coefficient_range':'1 <= n <= 1024'},
        'mod15':rows,'ogg_residue_multiplicities':dict(sorted(Counter(p%15 for p in O).items())),
        'CRT120':crt,'CM_curve':{'a_invariants':[0,4,0,2,0],'c4':c4,'discriminant':Delta,'j':8000,'rows':Erows},
        'Ogg_finite_certificate':ogg_cert,
        'shared_wheel':{'modulus':P,'count':len(surv_plus),'fraction':'192/1001'},
        'periodic_filter':{'unit_density':'8/15','unsigned_divergence_leading_coefficient':'4/15','signed_fundamental_period':30,'ramanujan15':ramanujan,'fourier_checked_n':[1,120],'M_d_intervals':{str(k):str(v) for k,v in scaled.items()},'unit_MRB':str(unit_C),'nonunit_MRB':str(nonunit_C),'unit_minus_density_times_MRB':str(unit_C-iv.mpf(8)*C/15)},
        'MRB_intervals':{'C':str(C),'A8':str(A(8)),'D':str(D),'q4':str(q4),'D_over_q4':str(D/q4),'analytic_tail':str(B),'nearest':nearest}}
(ROOT/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print('PASS:',CHECKS,'active checks')
print('Ogg mod15 unit residues:',units)
print('CRT120 classes:',len(crt),'largest least-prime witness:',max(v['prime_witness'] for v in crt))
print('CM point counts:',len(Erows),'odd primes through',max(v['p'] for v in Erows))
print('17 ordinary / Ogg; 37 supersingular / not Ogg: certified')
print('Period30, unit density 8/15, unsigned leading divergence 4/15: exact')
print('MRB D interval:',D)
print('MRB D/q^4 interval:',D/q4)
print('Two nearest-amplitude indices:',[x['index'] for x in nearest])
print('Unit-filtered MRB interval:',unit_C)
print('Nonunit-filtered MRB interval:',nonunit_C)
print('Bolza -1 twist trace relation: 167 odd primes; full Fourier projection: 120 exact cyclotomic checks')
