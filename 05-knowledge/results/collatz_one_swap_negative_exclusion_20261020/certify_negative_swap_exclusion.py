#!/usr/bin/env python3
"""Exact negative-clock one-swap certificate, preserving all integer exceptions.
External analytic input: real Matveev bound in Bugeaud--Mignotte--Siksek,
Theorem 9.4. Internal input: the proved negative-clock sqrt(6) content gate.
Every decision is exact integer or rational arithmetic.
"""
from fractions import Fraction as F
from math import factorial,gcd
from hashlib import sha256
from pathlib import Path
import json
CHECKS=0

def need(condition,*detail):
    global CHECKS
    CHECKS+=1
    if not condition:
        raise ArithmeticError(detail)

def log_interval(t,N=220):
    lower=2*sum((t**(2*j+1)/F(2*j+1) for j in range(N)),F(0))
    tail=2*t**(2*N+1)/((2*N+1)*(1-t*t))
    return lower,lower+tail

def carry(word):
    B=0
    for i,bit in enumerate(word):
        B=(3 if bit else 1)*B+bit*2**i
    return B

def actual_swaps(p,q):
    c=[((i+1)*p)//q-(i*p)//q for i in range(q)]
    for i in range(q):
        j=(i+1)%q
        if c[i]==c[j]:continue
        w=c.copy();w[i],w[j]=w[j],w[i]
        yield i,w

L2,U2=log_interval(F(1,3));L3,U3=log_interval(F(1,2))
alpha_lo,alpha_hi=L3/U2,U3/L2
Cupper=F(7,5)*30**5*23*F(7,10)*F(11,10)
need(Cupper<10**9,'Matveev coefficient')
need(U2<F(7,10) and U3<F(11,10),'Matveev height bounds')
L54,U54=log_interval(F(1,9))
need(U2+12*(3*U2+U54)<30,'cutoff logarithm')
need(F(10**12,6)>10**9*31,'cutoff endpoint')
threshold=100
exp_lower=sum((F(100,6)**j/factorial(j) for j in range(60)),F(0))
need(F(4*threshold)/exp_lower<L2,'negative Legendre endpoint')
need(F(2,3)**2<F(1,2),'small logarithm factor')
need(F(4)**3<F(5,2)**5,'sqrt6 plus bound')

lo,hi=alpha_lo,alpha_hi
hprev,h=0,1;kprev,k=1,0
cf=[];rows=[]
for index in range(100):
    digit=lo.numerator//lo.denominator
    need(digit==hi.numerator//hi.denominator,'ambiguous CF digit',index)
    cf.append(digit)
    hprev,h=h,digit*h+hprev
    kprev,k=k,digit*k+kprev
    q,p=h,k
    negative_lower=p*L3-q*U2
    negative_upper=p*U3-q*L2
    if threshold<=p<10**12:
        if negative_upper<0:
            status='above alpha; incompatible with negative clock'
        else:
            need(negative_lower>0,'ambiguous negative sign',q,p)
            cap=min(p//6,256)
            upper=F(2,2**cap)
            need(negative_lower>upper,'unexcluded negative convergent',q,p)
            status='negative logarithmic gap exceeds 2 exp(-p/6)'
        rows.append({'q':q,'p':p,'status':status,
                     'negative_gap_lower':str(negative_lower),
                     'negative_gap_upper':str(negative_upper)})
    if p>=10**12:break
    lo,hi=1/(hi-digit),1/(lo-digit)
else:raise ArithmeticError('missing terminal convergent')

small=[];small_swaps=0;integer_hits=[]
for p0 in range(2,threshold):
    q0=(3**p0).bit_length()-1
    interval_lo=p0*alpha_lo;interval_hi=p0*alpha_hi
    need(interval_lo.numerator//interval_lo.denominator==q0,'critical floor lower',p0)
    need(interval_hi.numerator//interval_hi.denominator==q0,'critical floor upper',p0)
    if gcd(p0,q0)!=1:
        small.append({'p':p0,'q':q0,'status':'noncoprime shape; handled by repetition lemma','swaps':0});continue
    D=2**q0-3**p0
    need(D<0,'negative clock')
    row_count=0;max_g=0;row_hits=[]
    for i,w in actual_swaps(p0,q0):
        B=carry(w)
        ones=[j for j,bit in enumerate(w) if bit]
        B2=sum(2**j*3**(p0-1-rank) for rank,j in enumerate(ones))
        need(B==B2,'carry implementations',p0,q0,i)
        g=gcd(-D,B)
        need(g*g<6**p0,'negative norm gate',p0,q0,i)
        if g==-D:
            hit={'p':p0,'q':q0,'swap_start':i,'word':''.join(map(str,w)),'D':D,'fixed_point':B//D}
            integer_hits.append(hit);row_hits.append(hit)
        else:
            need(g<-D,'proper gcd')
        max_g=max(max_g,g);row_count+=1;small_swaps+=1
    small.append({'p':p0,'q':q0,'status':'integer exceptions retained' if row_hits else 'all actual swaps excluded',
                  'swaps':row_count,'maximum_content':max_g})
expected={(2,3,'101',-7),(2,3,'110',-5)}
need({(v['p'],v['q'],v['word'],v['fixed_point']) for v in integer_hits}==expected,'unexpected integer exceptions')

hostile_count=0;hostile_hits=[]
for q0 in range(2,121):
    for p0 in range(1,q0):
        if gcd(p0,q0)!=1:continue
        D=2**q0-3**p0
        if D>=0:continue
        for i,w in actual_swaps(p0,q0):
            B=carry(w);g=gcd(-D,B)
            need(g*g<6**p0,'all-negative norm hostile',p0,q0,i)
            if g==-D:
                hostile_hits.append((p0,q0,''.join(map(str,w)),B//D))
            hostile_count+=1
need(set(hostile_hits)==expected,'broad hostile census exceptions')
# Actual negative cycles remain controls, not claimed exclusions.
cycles=[]
for start in (-5,-17):
    values=[];bits=[];n=start
    while n not in values:
        values.append(n);bits.append(n%2)
        n=(3*n+1)//2 if n%2 else n//2
    need(n==start,'cycle did not close at start',start)
    p0=sum(bits);q0=len(bits);D=2**q0-3**p0
    need(carry(bits)==start*D,'known negative cycle carry',start)
    c=[((i+1)*p0)//q0-(i*p0)//q0 for i in range(q0)]
    distance=min(sum(a!=b for a,b in zip(c,bits[s:]+bits[:s])) for s in range(q0))
    cycles.append({'start':start,'values':values,'word':''.join(map(str,bits)),'p':p0,'q':q0,'D':D,
                   'minimum_hamming_distance_to_mechanical_over_rotations':distance})
need(cycles[1]['minimum_hamming_distance_to_mechanical_over_rotations']==4,'minus17 hostile geometry')

print(json.dumps({'status':'EXACT negative-clock one-swap certificate, with genuine exceptions preserved',
    'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'checks':CHECKS,
    'rate_denominator':6,'logarithmic_factor':2,'p_cutoff_exclusive':10**12,'legendre_threshold':100,
    'matveev_prefactor_upper':str(Cupper),'log_series_terms':220,
    'continued_fraction_prefix':cf,'terminal_convergent':{'q':q,'p':p},
    'large_convergents':rows,'small_shapes':small,'small_swap_occurrences':small_swaps,
    'small_coprime_shapes':sum(v['swaps']>0 for v in small),'integer_exceptions':integer_hits,
    'broad_hostile_universe':'q<=120; every primitive negative mechanical shape; all cyclic unequal adjacent swaps',
    'broad_hostile_swap_occurrences':hostile_count,'known_negative_cycle_controls':cycles,
    'log_intervals':{'log2_lower':str(L2),'log2_upper':str(U2),'log3_lower':str(L3),'log3_upper':str(U3)}},indent=2))
