#!/usr/bin/env python3
"""Small exact controls for the symbolic height-controlled entry proof."""
from fractions import Fraction as F
from itertools import product
import json, math
checks=0
def check(v,msg):
    global checks
    checks+=1
    if not v:raise RuntimeError(msg)
def branch(w):
    a,b,m=1,0,1
    for bit in w:
        if bit:a,b=3*a,3*b+m
        m*=2
    return a,b,m
W=tuple(map(int,'11011011011011010101101110101010'))
a,b,m=branch(W);r=(-b*pow(a,-1,m))%m
check((a,b,m,r)==(3486784401,15810005461,4294967296,2962948091),'W')
check(b<20*(m-a) and r>20,'all-lift W descent')
check((a*r+b)//m==2405410912,'anchor endpoint')
q20=F(10**20,9**9*11**11)
check(0<q20<1,'Chernoff q<1')
check(q20<F(999,1000)**20,'rational q bound')
check(3**5<2**8,'rational slope estimate')
check(F(15,16)-F(91,100)==F(11,400),'height exponent margin')
check(73*F(11,400)>=2,'height threshold')
check(F(1,2**38)<F(3,4000),'exception exponent separation')
# Direct finite parity words independently authenticate both moment formulas.
for j in range(13):
    sumA=F(0);sumB=F(0);count=0
    residues=set()
    for w in product((0,1),repeat=j):
        aa,bb,mm=branch(w)
        sumA+=F(aa,mm);sumB+=F(bb,mm);count+=1
        residue=(-bb*pow(aa,-1,mm))%mm if mm>1 else 0
        check(residue not in residues,'finite parity bijection')
        residues.add(residue)
        n=residue
        actual=[]
        for _ in range(j):
            actual.append(n%2)
            n=(3*n+1)//2 if n%2 else n//2
        check(tuple(actual)==w,'direct finite parity replay')
    check(sumA==count,'mean multiplier')
    check(sumB==F(j,4)*count,'mean carry')
# Analytic mean-carry recurrence at larger lengths, no word enumeration.
EB=F(0)
for j in range(257):
    check(EB==F(j,4),'symbolic carry recurrence')
    EB=F(1,2)*(EB/2+(3*EB+1)/2)
# Independent tiny-pattern controls for the disjoint-window counting mechanism.
for L in range(1,5):
    for k in range(1,4):
        target=(1,)+(0,)*(L-1)
        count=0
        for w in product((0,1),repeat=L*k):
            count+=all(w[t*L:(t+1)*L]!=target for t in range(k))
        check(count==(2**L-1)**k,'exact disjoint block count')
# Floors, ceilings, all candidate words wholly inside the read parity window.
for mm in range(512,20001):
    N=7*mm//8;J=(3*mm+3)//4;k=(N-J)//32
    check(0<J<=N<=mm and k>=1,'admissible time window')
    check(J+32*k<=N,'last inspected block completes')
    check(F(k)>=F(mm,256)-2,'floor loss bound')
    check(F(11*J,20)<=F(11*N,20),'ordered count thresholds')
# Approximate displays are not used in any certificate decision.
delta=-math.log1p(-2**-32)/math.log(2)/256
q=float(q20)**(1/20)
print(json.dumps({'status':'PASS','explicit_checks':checks,
 'scope':'symbolic all-size proof plus finite exact controls, not orbit or source census',
 'word':''.join(map(str,W)),'length':len(W),'residue':r,'modulus':m,
 'height_exponent':'15/16','time_factor':'7/8','delta_display':delta,
 'chernoff_q20':str(q20),'chernoff_q_display':q,
 'dyadic_parameter_control':'all integers m in [512,20000]',
 'independent_moment_and_parity_control':'all parity words through length 12'},indent=2,sort_keys=True))
