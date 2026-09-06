"""Exact denominator23 bank, gap transport, and physical gluing controls."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys
sys.stdout.reconfigure(newline='\n')
N=0
V=(185370716505,189592176609,193905909065,198314972625,202822562745,205114343115)
U=(1,49,109,331,109561,36264691,12003612721)
T=36883259192

def need(x,why):
    global N
    N+=1
    if not x:raise ArithmeticError(why)

def norm(x):
    x%=1
    return min(x,1-x)

def gap(points):
    points=sorted(set(points))
    return max([b-a for a,b in zip(points,points[1:])]+[1+points[0]-points[-1]])

def safe(values,x):return min(norm(v*x) for v in values)

def atlas(a,b):
    d=gcd(a,b);n=(a+b)//d
    if n>356:return False
    k=2
    while k*k<=n:
        e=0
        while n%k==0:n//=k;e+=1
        if e and (k%3!=2 or e>2):return False
        k+=1
    return n==1 or n%3==2

def connected(n,edges):
    seen={0}
    while True:
        old=len(seen)
        for i,j in edges:
            if i in seen or j in seen:seen.update((i,j))
        if len(seen)==old:return len(seen)==n

def main():
    B=tuple(k for k in range(23) if all(14*min(v*k%23,(-v*k)%23)>23 for v in V))
    need(B==(3,5,6,10,11,12,13,17,18,20),'complete primitive phase bank')
    allowed=tuple(r for r in range(23) if all(14*min(r*k%23,(-r*k)%23)>23 for k in B))
    need(allowed==(1,3,5,6,10,11,12,13,17,18,20,22),'complete phase-preserving speed residues')
    for r in allowed:
        for k in B:need(norm(F(r*k,23))>=F(2,23)>F(1,14),'strict bank margin')
    print('V_RESIDUES',tuple(v%23 for v in V),'BANK_NUMERATORS',B,'ALLOWED_SPEED_RESIDUES',allowed)
    expected={0:23,1:6,2:4,3:10,4:5,5:5,6:6,7:4,8:5,9:7,10:8,11:6,
              12:6,13:8,14:7,15:5,16:4,17:6,18:5,19:5,20:10,21:4,22:6}
    for g in range(23):
        rho=gap(F(g*k,23)%1 for k in B)
        need(rho==F(expected[g],23),'exact residue gap')
        need(rho==gap(F((g+23)*k,23)%1 for k in B),'scale residue covariance')
    print('GAP_NUMERATOR_BY_GMOD23',tuple(expected[g] for g in range(23)))

    # Direct full union, including branches with noncoprime t and g.
    # d=gcd(t,g) is retained: gap=(d/t)rho_(g/d), not rho_g/t.
    cases=0
    for t in range(1,32):
        for g in range(1,47):
            points={F(g*(k+23*j),23*t)%1 for k in B for j in range(t)}
            d=gcd(t,g)
            rho=F(expected[(g//d)%23],23)
            need(gap(points)==F(d,t)*rho,'literal union maxgap and effective order')
            need(len(points)==(t//d)*(1 if (g//d)%23==0 else len(B)),'literal grid-union cardinality')
            cases+=1
    print('LITERAL_UNION_BANK',cases,'small(t,g) pairs, coprime and noncoprime')

    L=max(U)
    need(28*L>3*T,'one arbitrary-phase universal supplier gate misses')
    need(56*L<23*T,'ten-phase gate strict improvement')
    need(gap(F(k,23) for k in B)==F(6,23),'g1 exact maximal gap')
    x=F(11,23)
    eta=(T*x)%1
    need(eta==F(13,23) and int(eta*23) in B,'incoming phase is bank lift')
    j=T*x-eta
    need(j.denominator==1 and 0<=j<T,'incoming branch label')
    row=tuple(T*v for v in V)+U
    need(safe(row,x)==F(3,23),'literal incoming physical phase')
    print('INCOMING', 't',T,'L',L,'g1_Lcap',23*T//56,'strict_gate_slack',23*T-56*L,
          'phase',x,'eta',eta,'branch',j,'clearance',safe(row,x))
    passing=tuple(g for g in range(1,23) if 28*expected[g]*L<=69*T)
    need(passing==(1,2,4,5,6,7,8,9,11,12,14,15,16,17,18,19,21,22),'incoming aspect-ratio residue classes')
    print('ASPECT_RATIO_GOOD_GMOD23',passing,'uniform_nonzero_Lcap',69*T//280)
    # The incoming U has a stronger known phase; do not claim every old
    # single-phase method misses this already-safe example.
    need(all(u%2 for u in U) and safe(U,F(1,2))==F(1,2),'all-odd stronger phase')
    need(F(6*T,7*L)>1,'actual all-odd arc activates inherited single-phase gluing')

    # A different, mixed-parity nonunit seven-body. This control is a direct
    # physical row, not an asserted actual decoder equality entry.
    U2=(2,3,4,5,6,7,L-2)
    need(safe(U2,F(1,8))==F(1,8) and T%8==0,'exact seven-body supplier phase')
    phase=F(1,8)+F(3,23*T)
    row2=tuple(T*v for v in V)+U2
    clearance=F(1,8)-F(3*(L-2),23*T)
    need(safe(row2,phase)==clearance>F(1,14),'new mixed-parity physical row')
    print('MIXED_NONUNIT_DIRECT_CONTROL','U',U2,'phase',phase,'clearance',clearance)

    # Mixed parity and no unit in EITHER primitive shape. Compute the actual
    # strict inert atlas and exclude every mixed support in both orientations.
    Vm=(185370716505,189592176609,193905909065,198314972625,201676672560,205114343115)
    Um=(2,3,6,331,109561,36264691,12003612721)
    Q=91**6
    need(all(v%23 in allowed for v in Vm),'mixed six-star retains the bank')
    for shape in (Vm,Um):
        need(gcd(*shape)==1 and min(shape)>1 and {v%2 for v in shape}=={0,1},'primitive mixed unitless shape')
        edges=[(i,j) for i,j in combinations(range(len(shape)),2) if atlas(shape[i],shape[j])]
        need(connected(len(shape),edges),'strict actual atlas connected')
        print('ACTUAL_MIXED_SHAPE',shape,'EDGES',edges)
    rowm=tuple(T*v for v in Vm)+Um
    need(len(set(rowm))==13 and gcd(*rowm)==1 and sum(rowm)<Q*Q,'physical actual row and box')
    for v in Vm:
        for u in Um:need(not atlas(T*v,u),'no actual cross edge')
    coeff=[]
    for v,w in combinations(Vm,2):
        D=T*gcd(v,w)
        for u in Um:
            c=D//gcd(D,u);coeff.append(c)
            need(c>Q,'two small labels force outside coefficient above Q')
    ratios=[]
    for u,w in combinations(Um,2):
        D=gcd(u,w)
        for v in Vm:
            delta=gcd(D,T*v)
            ratio=F((T*v//delta)*D,Q*(u+w));ratios.append(ratio)
            need(ratio>1,'two large labels cannot reach smallest permitted target')
    need(len(coeff)+len(ratios)==231,'both support orientations complete')
    need(min(coeff)==36446911794129044,'distinguished coefficient minimum')
    need(min(ratios)==F(3112005545498710680,3102632035080254131),'amplitude minimum')
    need(safe(Um,F(1,8))==F(1,8),'literal mixed seven-body supplier')
    need(safe(rowm,phase)==F(26509924715,212078740354)>F(1,14),'mixed unitless actual physical certificate')
    need(28*max(Um)>3*T and 56*max(Um)<23*T,'bank closes beyond supplier-only one-phase gate')
    print('ACTUAL_MIXED_ROW','sum',sum(rowm),'coefficient_min',min(coeff),'amplitude_ratio_min',min(ratios),
          'phase',phase,'clearance',safe(rowm,phase),'mixed_supports',len(coeff)+len(ratios))
    # Scalar gcd maxima are displayed only: no claim of a full JSON profile audit.
    maxima=tuple(max(gcd(*c) for c in combinations(rowm,k)) for k in range(12,6,-1))
    need(maxima==(1,1,1,1,2,3),'literal subset gcd maxima')
    print('ACTUAL_MIXED_BODY_GCD_MAXIMA_12_TO_7',maxima)

    # The fixed incoming all-odd star has cheaper parity banks: denominator23
    # is a residue-class certificate, not an optimum over all possible banks.
    odd12=tuple(k for k in range(12) if all(norm(F(r*k,12))>F(1,14) for r in (1,3,5,7,9,11)))
    need(odd12==(1,2,3,5,6,7,9,10,11) and gap(F(k,12) for k in odd12)==F(1,6),
         'all-odd denominator12 stopping control')
    print('ALL_ODD_OVERLAP','denominator12_bank',odd12,'gap',F(1,6),'less_than23gap',F(6,23))

    # Incoming unitless5+8 is governed by a favorable offset, not worst gap.
    V5=(40341259,287243635,542783995,807423715,14321146945)
    U8=(25806,32538,44022,68034,249458,374187,748374,2245122)
    t5=4912751;x5=F(1403643,9825502)
    need(all(v%2 for v in V5) and all(u%7 for u in U8),'incoming parity/seven-unit scopes')
    need((t5*x5)%1==F(1,2) and abs(x5-F(1,7))==F(1,14*t5),'incoming exact favorable offset')
    need(max(U8)<t5 and safe(tuple(t5*v for v in V5)+U8,x5)==F(696962,t5)>F(1,14),
         'incoming literal alignment certificate')
    # For every t not divisible by7 choose k with t*k=3 modulo7. The nearest
    # half-grid point differs from k/7 by exactly1/(14t).
    for ts in range(1,85):
        if ts%7==0:continue
        k=3*pow(ts,-1,7)%7
        j=(ts*k)//7
        xs=F(2*j+1,2*ts)
        need(abs(xs-F(k,7))==F(1,14*ts) and (ts*xs)%1==F(1,2),'general seven-clock offset identity')
    print('INCOMING5+8_ALIGNMENT','t',t5,'phase',x5,'clearance',F(696962,t5),
          'general_scope','odd V, all U seven-units, 7 does not divide t, maxU<t')

    # At equality a grid point may be an endpoint only. Strict V margins
    # allow a small physical perturbation into the U arc.
    t=56
    Lb=23
    Ub=(2,3,4,5,6,7,Lb)
    boundary=F(1,8)+F(3,23*t)
    rowb=tuple(t*v for v in V)+Ub
    need(56*Lb==23*t and safe(rowb,boundary)==F(1,14),'inclusive gate endpoint equality')
    epsilon=min(F(5,322)/(4*t*max(V)),F(1,10**6))
    perturbed=boundary-epsilon
    need(safe(rowb,perturbed)>F(1,14),'strict-margin perturbation repairs equality')
    print('EQUALITY_CONTROL','t',t,'L',Lb,'endpoint_clearance',safe(rowb,boundary),
          'epsilon',epsilon,'perturbed_clearance',safe(rowb,perturbed))

    # Same bank size, different gap; a noncoprime collapse is also material.
    need(expected[2]==4 and expected[3]==10,'cardinality does not determine gap')
    need(expected[0]==23,'23-divisible multiplier collapses the bank')
    need(gap({F(2*(k+23*j),23*2)%1 for k in B for j in range(2)})==F(6,23)
         !=F(expected[2],46),'dropping gcd produces false maxgap')
    print('SCOPE: exact bank/gluing sidecar; no universal finite-bank completeness or optimality; actual entry unnecessary')
    print('PASS',N,'always-active exact gates')

if __name__=='__main__':main()
