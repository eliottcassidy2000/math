"""Independent rational Sturm/whole-box referee; no producer import."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json
import sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
HERE=Path(__file__).resolve().parent
STEM='continuing9_20260907_fixed_moment_circuits'
G=0
def gate(ok,label):
    global G
    G+=1
    if not ok: raise RuntimeError(label)
def locate(name):
    for p in (HERE/name,HERE.parent/'05-knowledge/results'/name,
              Path('C:/w/continuing9_20260907_moments')/name):
        if p.is_file(): return p
    raise FileNotFoundError(name)
for suffix,pin in (('.py','a71ee4509bb46644f8e1fcf3a5f6931941363876540ea8565d1addf46e6e8102'),
                   ('.out','5d51e2a2aa33a1dd55d4d3caaa1c64702567ca4787e78a8b5d1e8c4f0486abfe'),
                   ('_certificate.json','291b62a5638f8a057ccf514851cb0b5d100a2a61a3da8eff9151b49040be11b9')):
    raw=locate(STEM+suffix).read_bytes()
    gate(sha256(raw).hexdigest()==pin,'frozen producer '+suffix)
    if suffix=='.out': gate(b'\r' not in raw,'producer raw LF')
data=json.loads(locate(STEM+'_certificate.json').read_text(encoding='utf-8'))
def trim(f):
    f=list(map(Q,f))
    while len(f)>1 and not f[-1]: f.pop()
    return f
def ev(f,x):
    out=Q(0)
    for a in reversed(f): out=out*x+a
    return out
def deriv(f): return trim([i*f[i] for i in range(1,len(f))] or [0])
def rem(f,g):
    f=trim(f); g=trim(g)
    while f!=[0] and len(f)>=len(g):
        n=len(f)-len(g); factor=f[-1]/g[-1]
        for j,a in enumerate(g): f[n+j]-=factor*a
        f=trim(f)
    return f
def sturm(f):
    a=trim(f); b=deriv(a); out=[a,b]
    while b!=[0]:
        c=[-q for q in rem(a,b)]
        if c==[0]: break
        c=[q/abs(c[-1]) for q in c]
        out.append(c); a,b=b,c
    return out
def variation(values):
    signs=[1 if v>0 else -1 for v in values if v]
    return sum(a!=b for a,b in zip(signs,signs[1:]))
def var_at(chain,x):
    return variation([ev(f,x) for f in chain])
def count(chain,a,b):
    va=var_at(chain,a)
    vb=variation([f[-1] for f in chain]) if b is None else var_at(chain,b)
    return va-vb
def reflected(f): return [(-1)**k*a for k,a in enumerate(f)]
def center(d,q): return [Q(comb(d,k))*q**(k*(k-1)//2) for k in range(d+1)]
def ratios(e):
    d=len(e)-1; h=[e[k]/comb(d,k) for k in range(d+1)]
    R=[h[k]*h[k]/(h[k-1]*h[k+1]) for k in range(1,d)]
    C=[R[k]/R[k-1] for k in range(1,len(R))]
    return R,C
def reconstruct(d,q,cs):
    h=[Q(1),Q(1)]; rk=1/q
    for k in range(1,d):
        if k>=2: rk*=cs[k-2]
        h.append(h[-1]*h[-1]/(rk*h[-2]))
    return [comb(d,k)*h[k] for k in range(d+1)]
def refine(f,chain,interval):
    lo,hi=interval
    if ev(f,hi)==0: return hi,hi
    mid=(lo+hi)/2
    return (lo,mid) if count(chain,lo,mid) else (mid,hi)
def isolate_positive(f):
    chain=sturm(f); degree=len(f)-1
    gate(count(chain,Q(0),None)==degree,'fresh exact positive root count')
    upper=Q(1)
    while count(chain,Q(0),upper)<degree: upper*=2
    pending=[(Q(0),upper)]; intervals=[]
    while pending:
        lo,hi=pending.pop(); n=count(chain,lo,hi)
        if n==0: continue
        if n==1:
            for _ in range(34): lo,hi=refine(f,chain,(lo,hi))
            intervals.append((lo,hi)); continue
        mid=(lo+hi)/2; pending.extend(((mid,hi),(lo,mid)))
    return sorted(intervals)

# Fresh induction controls use no supplied roots or producer root solver.
fresh=0
for q in (Q(1,4),Q(3,5),Q(9,10)):
    previous=[(Q(1),Q(1))]
    for d in range(2,8):
        coeff=center(d,q); f=reflected(coeff)
        old=center(d-1,q)
        expected=[Q(0)]*(d+1)
        for k,a in enumerate(old):
            expected[k]+=a; expected[k+1]+=a*q**k
        gate(expected==coeff,'independent finite recurrence coefficients')
        now=isolate_positive(f)
        gate(all(now[i+1][0]>now[i][1]/q for i in range(d-1)),'fresh strict q-separated centers')
        gate(all(now[i][1]<previous[i][0] and now[i+1][0]>previous[i][1]/q for i in range(d-1)),
             'entire induction interval ordering at both boundaries')
        previous=now; fresh+=1

universe=[]; targets=0; brackets=0
for case in data['cases']:
    d=case['d']; q=Q(case['q']); delta=Q(case['delta'])
    coeff=center(d,q); f=reflected(coeff); chain=sturm(f)
    gate(coeff==list(map(Q,case['center_coefficients'])),'exact independently reconstructed center')
    gate(ratios(coeff)==([1/q]*(d-1),[Q(1)]*(d-2)),'fixed-R1 all-one center')
    roots=[tuple(map(Q,pair)) for pair in case['positive_roots_of_P_minus_X']]
    for lo,hi in roots:
        gate(lo>0 and hi>=lo,'positive retained isolator')
        gate(ev(f,lo)==0 and ev(deriv(f),lo)!=0 if lo==hi else ev(f,lo)!=0 and ev(f,hi)!=0 and count(chain,lo,hi)==1,
             'independent Sturm retained root bracket')
        brackets+=1
    gate(all(roots[i+1][0]>roots[i][1]/q for i in range(d-1)),'supplied brackets certify strict q separation')
    samples=list(map(Q,case['samples']))
    vals=[ev(f,x) for x in samples]
    weights=[sum(comb(k,3)*coeff[k]*x**k for k in range(3,d+1)) for x in samples]
    gate(vals==list(map(Q,case['sample_values'])) and weights==list(map(Q,case['weights'])),'literal values and entire coefficient weights')
    bound=min([Q(1,2),Q(1,2*comb(d,3))]+[abs(v)/(4*w) for v,w in zip(vals,weights)])
    gate(bound==Q(case['raw_bound']) and 0<delta<=bound,'exact uniform rational radius')
    gate(all(a<b for a,b in zip([Q(0)]+samples,samples)),'strictly ordered sign addresses')
    # A second whole-box proof: sum exact coefficient intervals, no Bernoulli bound.
    for ell,x in enumerate(samples,1):
        lower=upper=Q(0)
        for k,a in enumerate(coeff):
            E=comb(k,3) if k>=3 else 0
            al=a/(1+delta)**E; ah=a/(1-delta)**E
            if k%2: lower-=ah*x**k; upper-=al*x**k
            else: lower+=al*x**k; upper+=ah*x**k
        gate((lower>0 if ell%2==0 else upper<0),'direct whole-box coefficient envelope sign')
        gate((lower if ell%2==0 else -upper)>=abs(vals[ell-1])/2,'whole-box half-margin independently retained')
    seen=set()
    for row in case['targets']:
        word=tuple(row['word']); seen.add(word)
        cs=[1+delta*w for w in word]
        es=reconstruct(d,q,cs)
        gate(cs==list(map(Q,row['C'])) and es==list(map(Q,row['coefficients'])),'native ratio-recursive target coefficient reconstruction')
        gate(es[:3]==coeff[:3],'first two moments exactly unchanged')
        rr,cc=ratios(es)
        gate(cc==cs and [0 if v==1 else 1 if v>1 else -1 for v in cc]==list(word),'literal target ratios and exact ties')
        ff=reflected(es)
        gate(count(sturm(ff),Q(0),None)==d and len(sturm(ff)[-1])==1,'direct Sturm target simplicity and positivity')
        gate(all(((-1)**ell)*ev(ff,x)>0 for ell,x in enumerate(samples,1)),'literal target sign samples')
        targets+=1
    gate(seen==set(product((-1,0,1),repeat=d-2)),'complete ternary universe without quotient omission')
    universe.append([d,str(q),str(delta),len(seen)])
gate(targets==132 and brackets==25,'complete frozen scope')
anchor=list(map(Q,data['anchor_S13_Q59']))
base=center(5,Q(275,338))
gate(anchor==[a*Q(13,5)**k for k,a in enumerate(base)],'actual anchor scaling')
gate(anchor[1]==13 and anchor[1]**2-2*anchor[2]==59,'literal anchor first two moments')

# Independent exact cubic algebra, endpoint geometry and original hostiles.
n,aa,v,prod=s.symbols('n aa v prod')
poly=n**3-3*n*n+3*(1-v)*n-prod
disc=s.discriminant(poly,n)
gate(s.expand(disc-27*(4*v**3-(prod-(1-3*v))**2))==0,'full cubic discriminant')
LL=(1+aa)**2*(1-2*aa); UU=(1-aa)**2*(1+2*aa); qq=1-aa*aa
gate(s.expand(poly.subs({v:aa*aa,prod:LL})-(n-(1+aa))**2*(n-(1-2*aa)))==0,'lower repeated endpoint')
gate(s.expand(poly.subs({v:aa*aa,prod:UU})-(n-(1-aa))**2*(n-(1+2*aa)))==0,'upper repeated endpoint')
gate(s.cancel(qq**3/UU-1+aa**3*(2+aa)/(1+2*aa))==0,'strict lower circuit boundary below one')
gate(s.cancel(qq**3/LL-1-aa**3*(2-aa)/(1-2*aa))==0,'strict upper circuit boundary above one')
for a0 in (Q(1,5),Q(1,2),Q(3,4)):
    low=(1+a0)**2*(1-2*a0); high=(1-a0)**2*(1+2*a0)
    for position in (Q(-1,4),Q(1,4),Q(1,2),Q(3,4),Q(5,4)):
        pv=max(Q(0),low)+position*(high-max(Q(0),low))
        roots_poly=[-pv,3*(1-a0*a0),Q(-3),Q(1)]
        positive=count(sturm(roots_poly),Q(0),None)
        gate((positive==3)==(max(0,low)<pv<high),'literal cubic full strict image controls')
for key,degree in (('cubic',3),('quartic',4)):
    es=list(map(Q,data['hostiles'][key+'_e']))
    rr,cc=ratios(es)
    gate(rr==list(map(Q,data['hostiles'][key+'_R'])) and cc==list(map(Q,data['hostiles'][key+'_C'])),'original hostile normalized ratios')
    gate(all(r>1 for r in rr),'all hostile Newton inequalities strict')
    gate(count(sturm(reflected(es)),Q(0),None)==degree-2,'exact missing real-root pair despite Newton inequalities')
    if degree==3:
        gate(s.discriminant(sum(s.Rational(v.numerator,v.denominator)*n**(3-k) for k,v in enumerate(es)),n)==-s.Rational(25515,65536),'original cubic hostile discriminant')
gate(Q(338,275)*Q(1,2)<1,'anchored half-amplitude necessary-condition hostile')
print('FRESH CENTERS',fresh,'independent rational Sturm and complete inductive ordering')
print('FROZEN BOXES [d,q,delta,ternary count]',json.dumps(universe,separators=(',',':')))
print('WHOLE BOX: direct coefficient-interval envelopes retain every sign with at least half the center margin')
print('TARGETS',targets,'all literal ratios, moments, simple-positive Sturm counts; root brackets',brackets)
print('ANCHOR S13 Q59: exact radius1/2048, all27 words, no interlacer predicate inferred')
print('CUBIC: complete discriminant boundaries and threshold Q=S^2/2; strict Newton hostiles verified')
print('PASS',G,'always-active exact gates; raw LF')
print('Scope: positive real roots and rational coefficients where inputs rational; no rational-root or integer-root assertion.')
