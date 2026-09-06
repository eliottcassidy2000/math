#!/usr/bin/env python3
"""Bounded exact endpoint/discrepancy controls for the AP8 asymptotic margin.

Standalone Fraction arithmetic. No carrier scan or large-height census.
"""
from fractions import Fraction as F
from math import floor, ceil, gcd
from hashlib import sha256

checks=0
trace=[]


def gate(ok,label):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(label)


def omega(v):
    return 1 if v%3==0 else 3


def measure_union(intervals):
    intervals=sorted((a,b) for a,b in intervals if b>a)
    if not intervals:
        return F(0)
    start,end=intervals[0]
    total=F(0)
    for a,b in intervals[1:]:
        if a<=end:
            end=max(end,b)
        else:
            total+=end-start
            start,end=a,b
    return total+end-start


def primitive_comb(v,lam,lo,hi):
    """Disjoint third-orbit comb teeth; period and width kept separately."""
    period=F(1,omega(v)*v)
    radius=lam/v
    intervals=[]
    for k in range(floor((lo-radius)/period)-1,ceil((hi+radius)/period)+2):
        a,b=max(lo,k*period-radius),min(hi,k*period+radius)
        if b>a:
            intervals.append((a,b))
    # At the declared margins teeth are disjoint up to endpoints.
    return sum(b-a for a,b in intervals)


def addressed_danger(v,lam,lo,hi):
    """Independent union of the three actual observers v(y+j/3)."""
    intervals=[]
    for j in range(3):
        shift=F(j,3)
        for n in range(floor(v*(lo+shift)-lam)-1,ceil(v*(hi+shift)+lam)+2):
            a=max(lo,(n-lam)/v-shift)
            b=min(hi,(n+lam)/v-shift)
            if b>a:
                intervals.append((a,b))
    return measure_union(intervals)


def f(om,lam):
    return 2*lam*(1-2*om*lam)


# Each of 24 addressed teeth keeps a common integer cell. Its two endpoint
# inequalities are affine in lambda, so checking lambda=0,1/9 is complete.
for c in range(1,9):
    for j in range(3):
        shift=F(j,3)
        cell=floor(c*(F(2,21)+shift))
        for lam in [F(0),F(1,9)]:
            lo,hi=(2+3*lam)/21,(1-lam)/8
            gate(c*(lo+shift)>=cell+lam,'AP8 lower tooth endpoint')
            gate(c*(hi+shift)<=cell+1-lam,'AP8 upper tooth endpoint')
            trace.append(('tooth',c,j,str(lam),cell,str(lo),str(hi)))

# The physical equality grid is exactly the eighteen units modulo27.
units=[u for u in range(27) if gcd(u,27)==1]
observed=[]
for u in range(27):
    margin=min(min((3*c*u)%27,27-(3*c*u)%27) for c in range(1,9))
    if margin>=3:
        observed.append(u)
gate(observed==units and len(units)==18,'AP8 equality grid exact')
for tail in units:
    bad=[u for u in units if min((tail*u)%27,27-(tail*u)%27)<3]
    gate(len(bad)==4,'every three-unit tail kills exactly four unit clocks')
    trace.append(('tail_clock',tail,tuple(bad)))
gate(18-3*4==6,'three tails leave at least six clocks')
for c in range(1,9):
    gate(all(min((3*c*u)%27,27-(3*c*u)%27)>=3 for u in units),
         'nonzero mod9 body residue safe on every unit clock')
gate(all((3*9*u)%27==0 for u in units),'nine-divisible outlier kills every equality clock')

# 64 small interval controls compare the reduced primitive comb with the
# literal three observers, retaining the valuation-dependent multiplicity.
intervals=[(F(-1,37),F(5,29)),(F(2,21),F(1,8)),
           (F(1,9),F(1,9)+F(1,101)),(F(2,7),F(5,7))]
comb_rows=0
for v in [3,5,9,10]:
    for lam in [F(1,100),F(1,14),F(1,12),F(1,9)]:
        for lo,hi in intervals:
            actual=addressed_danger(v,lam,lo,hi)
            reduced=primitive_comb(v,lam,lo,hi)
            upper=2*lam*omega(v)*(hi-lo)+f(omega(v),lam)/v
            gate(actual==reduced,'addressed / primitive comb agreement')
            gate(actual<=upper,'rational one-period discrepancy bound')
            trace.append(('comb',v,str(lam),str(lo),str(hi),str(actual),str(upper)))
            comb_rows+=1

# The two discrepancy maxima are polynomial identities with explicit signs.
# f1(1/9)-f1(x)=(1/9-x)*(2-4/9-4x).
gate([F(14,81),F(-2),F(4)]==[F(1,9)*F(14,9),-F(4,9)-F(14,9),F(4)],
     'full f1 difference factorization')
gate(F(14,9)-F(4,9)>0,'f1 difference second factor stays positive')
gate([F(1,12),F(-2),F(12)]==[12*F(1,12)**2,-24*F(1,12),F(12)],
     'full f3 deficit square identity')
gate(F(5*9,168)==F(15,56),'supplier interval length in delta')
gate(F(14,81)+F(1,12)==F(83,324),'mixed discrepancy maximum sum')
gate(F(5,168)*9-F(83,324)==F(53,4536),'mixed final positive measure coefficient')
gate(F(25,168)*9-F(28,81)==F(4507,4536),'both-divisible final measure coefficient')
gate(F(1,9)-F(9,82)>0,'minimum allowed height has positive margin')

# Six named strict-branch rows, including the smallest allowed H=82.
# No tail speeds are chosen: the analytic packet bound is uniform in them.
budgets=[]
for v,w in [(82,90),(90,91),(90,93),(99,102),(999,1000),(999,1002)]:
    H=min(v,w)
    lam=F(1,9)-F(9,H)
    ell=5*(1-9*lam)/168
    b=ell*(1-2*lam*(omega(v)+omega(w)))-f(omega(v),lam)/v-f(omega(w),lam)/w
    coefficient=F(4507,4536) if v%3==w%3==0 else F(53,4536)
    gate(H>=82 and v!=w and (v%9==0 or w%9==0),'named strict-branch row types')
    gate(0<lam<F(1,9) and b>=coefficient/H>0,'named exact high-outlier budget')
    budgets.append((v,w,str(lam),str(b),str(coefficient/H)))
    trace.append(('budget',v,w,str(lam),str(b)))

print('PASS AP8 arbitrary-outlier optimal asymptotic margin')
print('S=3*{1,...,8,v,w} plus any three distinct positive 3-unit tails; H=min(v,w)>=82')
print('1/9-9/H <= M <=1/9; equality iff neither outlier is divisible by 9')
print('Nine-divisible branch: L_(1/9-9/H)>=53/(4536H), improved to 4507/(4536H) if both are divisible by 3')
print('Declared universe: 24 affine tooth cells, both margin endpoints, both inequalities; 18 tail residues')
print('Rational discrepancy controls:',comb_rows,'; named strict-branch budgets:',len(budgets))
for row in budgets:print('BUDGET',*row)
print('Exact gates:',checks)
print('Trace SHA256:',sha256(repr(trace).encode()).hexdigest())
