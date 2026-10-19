#!/usr/bin/env python3
"""Supporting exact controls for the readable prime-pairs research note.
Run in the staged directory; writes only its named JSON/text results there.
No checkout, network, or Git operations. Floating-point outputs are labelled.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import cos, pi, prod, log
from pathlib import Path
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import json

def primes(n):
    s=bytearray(b'\x01')*(n+1); s[:2]=b'\x00\x00'
    for p in range(2,int(n**.5)+1):
        if s[p]: s[p*p:n+1:p]=b'\x00'*((n-p*p)//p+1)
    return [p for p in range(2,n+1) if s[p]]

def gamma(p,t):
    c=2 if t%p==0 else 3 if (9*t*t-1)%p==0 else 4
    return F(p*(p-c),(p-2)**2)

out={'status':'FINITE-EXACT identities; numerical constants labelled separately'}
# General finite cyclic-group theorem; includes nonsymmetric hostile sets.
general=0
for q in range(2,13):
    for bits in range(1<<q):
        B={i for i in range(q) if bits>>i&1}; A=set(range(q))-B
        if not A: continue
        s=len(B)
        tb=sum((a+b)%q in B for a in B for b in B)
        ta=sum((a+b)%q in A for a in A for b in A)
        assert ta==q*q-3*q*s+3*s*s-tb
        if B=={(-b)%q for b in B}:
            corr=[]
            for t in range(q):
                cc=sum((a-t)%q in A for a in A)
                ss=sum((t-a)%q in A for a in A)
                assert cc==ss==q-2*s+sum((b-t)%q in B for b in B)
                corr.append(F(q*cc,(q-s)**2))
            assert sum(corr)==q
            assert sum(corr[t] for t in A)/len(A)==1+F(s**3-q*tb,(q-s)**3)
            eb=sum(sum((b-t)%q in B for b in B)**2 for t in range(q))
            assert sum((c-1)**2 for c in corr)/q==F(q*eb-s**4,(q-s)**4)
            general+=1
out['general_cyclic_symmetric_cases']=general
# Same size and energy, opposite signs of the target-conditioned surplus.
hostiles=[]
for B in [{1,5,8,12},{1,2,11,12}]:
    q=13;s=len(B)
    energy=sum(sum((b-t)%q in B for b in B)**2 for t in range(q))
    triples=sum((a+b)%q in B for a in B for b in B)
    hostiles.append({'B':sorted(B),'energy':energy,'triples':triples,'conditional_mean':str(1+F(s**3-q*triples,(q-s)**3)),'variance':str(F(q*energy-s**4,(q-s)**4))})
assert [h['energy'] for h in hostiles]==[36,36]
assert [h['triples'] for h in hostiles]==[0,6]
out['equal_energy_opposite_surplus']=hostiles
# Every residue of every prime <= 199; direct enumeration independent of formula.
local=0
for p in primes(199):
    if p<5: continue
    a=pow(6,-1,p); A=set(range(p))-{a,(-a)%p}
    gs=[]
    for t in range(p):
        c=sum((t-x)%p in A for x in A)
        assert F(p*c,(p-2)**2)==gamma(p,t)
        gs.append(gamma(p,t)); local+=1
    assert sum(gs)/p==1
    assert sum((g-1)**2 for g in gs)/p==F(6*p-16,(p-2)**4)
    assert sum(gs[t] for t in A)/(p-2)==1+F(8,(p-2)**3)
    # Trigonometric spectrum is a floating check, exact result is proved in note.
    for k in range(1,p):
        real=sum(float(gs[t])*cos(2*pi*k*t/p) for t in range(p))/p
        expected=4*cos(2*pi*k*a/p)**2/(p-2)**2
        assert abs(real-expected)<2e-12
out['prime_residue_cases']=local
# Composite CRT wheel, all t, exact rational moments.
ps=[5,7,11]; W=prod(ps); A={a for a in range(W) if all((36*a*a-1)%p for p in ps)}
counts=[sum((t-a)%W in A for a in A) for t in range(W)]
gs=[F(W*c,len(A)**2) for c in counts]
assert all(gs[t]==prod(gamma(p,t) for p in ps) for t in range(W))
assert sum(gs)/W==1
assert sum(g*g for g in gs)/W==prod(1+F(6*p-16,(p-2)**4) for p in ps)
assert sum(gs[t] for t in A)/len(A)==prod(1+F(8,(p-2)**3) for p in ps)
out['crt_wheel']={'W':W,'allowed':len(A),'min_pair_residues':min(counts),'max_pair_residues':max(counts),'second_moment':str(sum(g*g for g in gs)/W),'conditioned_mean':str(sum(gs[t] for t in A)/len(A))}
# Uniform-in-target periodic interval bound, with wrap-around independently counted.
for m in range(W+1,2*W+50):
    rw=counts[m%W]; actual=sum(a%W in A and (m-a)%W in A for a in range(1,m))
    assert abs(F(actual)-F((m-1)*rw,W))<=rw
    assert actual>=((m-1)//W)*rw>0
out['rough_sum_interval_targets']=W+49
# Prime-power identities checked with exact rational arithmetic.
P=primes(500); f=[F(int(n in set(P))) for n in range(501)]; j=[F(0) for _ in f]
for p in P:
    n=p;k=1
    while n<=500:
        j[n]+=F(1,k);n*=p;k+=1
e=[a-b for a,b in zip(j,f)]
for N in range(2,501):
    rf=sum(f[a]*f[N-a] for a in range(1,N))
    rj=sum(j[a]*j[N-a] for a in range(1,N))
    cross=2*sum(f[a]*e[N-a] for a in range(1,N))
    ee=sum(e[a]*e[N-a] for a in range(1,N))
    assert rj-rf==cross+ee
    assert 0<=rj-rf<=2*sum(e[:N])
out['prime_power_at_10']={'weighted':str(sum(j[a]*j[10-a] for a in range(1,10))),'primes':str(sum(f[a]*f[10-a] for a in range(1,10)))}
# Boolean divisor-poset indistinguishability, using integer unnormalized masses.
r=4; Q=1<<r; signs=[(-1)**i.bit_count() for i in range(Q)]
for U in range(Q):
    for V in range(Q):
        totals=[]
        for sigma in [-1,1]:
            totals.append(sum(1+sigma*signs[S]*signs[T] for S in range(Q) for T in range(Q) if S&U==U and T&V==V))
        if U!=Q-1 or V!=Q-1: assert totals[0]==totals[1]
primepair=[]
for sigma in [-1,1]:
    primepair.append(sum(1+sigma*signs[S]*signs[T] for S in range(Q) for T in range(Q) if S.bit_count()==T.bit_count()==1))
assert primepair==[0,2*r*r]
out['boolean_prime_pair_masses']={'r':r,'minus':primepair[0],'plus':primepair[1],'total_mass':Q*Q}
# Constants are finite products, not certified infinite decimal values.
pp=primes(1_000_000); gp=vp=bp=cp=1.0
for p in pp:
    if p>=3: cp*=1-1/(p-1)**2
    if p>=5:
        gp*=1-4/(p-2)**2
        vp*=1+(6*p-16)/(p-2)**4
        bp*=1+8/(p-2)**3
out['floating_products_p_le_1000000']={'G':gp,'variance':vp-1,'internal_target_boost':bp,'twin_singular_constant':12*cp,'four_form_floor':(12*cp)**2*gp}
# Certified B interval: positive finite products with directed Decimal rounding.
# The infinite tail multiplier is <= 1/(1-x), x=4/(z-2)^2, by log(1+u)<=u.
z=1_000_000
with localcontext() as ctx:
    ctx.prec=45;ctx.rounding=ROUND_FLOOR
    low=Decimal(1)
    for p in pp:
        if p>=5:
            d=(p-2)**3
            low=(low*Decimal(d+8))/Decimal(d)
with localcontext() as ctx:
    ctx.prec=45;ctx.rounding=ROUND_CEILING
    high=Decimal(1)
    for p in pp:
        if p>=5:
            d=(p-2)**3
            high=(high*Decimal(d+8))/Decimal(d)
    D=(z-2)**2
    high=(high*Decimal(D))/Decimal(D-4)
assert low<high
out['certified_internal_target_boost']={'prime_cutoff':z,'lower':str(low),'upper':str(high),'method':'45-digit directed rounding and all-integers tail bound'}
# Arithmetic mean experiment for the exact finite correction weight w(m).
# For m<=N all exceptional primes lie <=3N+1, so all correction weights are included.
N=100_000; weights=[1.0]*(N+1)
for p in primes(3*N+1):
    if p<5: continue
    for m in range(p,N+1,p): weights[m]*=(p-2)/(p-4)
    inv=pow(3,-1,p)
    for r0 in (inv,(-inv)%p):
        for m in range(r0,N+1,p):
            if m>0: weights[m]*=(p-3)/(p-4)
out['arithmetic_moment_experiment']={'N':N,'using_G_p_le_1000000':True,'mean':gp*sum(weights[1:])/N,'second_moment':gp*gp*sum(w*w for w in weights[1:])/N,'minimum_at':[m for m in range(1,N+1) if weights[m]==1.0]}
# A small direct prime census is a hostile local/global comparison, not a claimed extension.
N=20_000; prime_set=set(primes(6*N+1)); K=[k for k in range(1,N+1) if 6*k-1 in prime_set and 6*k+1 in prime_set]; ks=set(K)
holes=[]
for m in range(2,N+1):
    if not any(m-a in ks for a in K if a<=m//2): holes.append(m)
out['finite_twin_rank_holes_2_to_20000']=holes
out['finite_twin_rank_701']={'representations':sum(701-a in ks for a in K if a<701),'local_weight':weights[701]}
out['all_checks_passed']=True
base=Path(__file__).resolve().parent
(base/'prime_pair_checks.json').write_text(json.dumps(out,indent=2)+'\n')
(base/'prime_pair_checks.txt').write_text(json.dumps(out,indent=2)+'\nALL CHECKS PASSED\n')
print(json.dumps(out,indent=2)); print('ALL CHECKS PASSED')
