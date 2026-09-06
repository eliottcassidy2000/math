#!/usr/bin/env python3
"""All-A2 alpha completion: independent coefficient maps and exact controls.

The theorem is proved in the companion report; finite sign banks are controls.
No repository mathematical producer is imported. SymPy only isolates roots.
"""
from math import comb, gcd
from fractions import Fraction as F
import hashlib, json, sys
import sympy as S
sys.stdout.reconfigure(newline='\n')
CHECKS=0
T=S.symbols('t')
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(label)
def choose(n,k): return comb(n,k) if 0<=k<=n else 0
def mul(a,b):
    c={}
    for i,x in a.items():
        for j,y in b.items(): c[i+j]=c.get(i+j,0)+x*y
    return {i:x for i,x in c.items() if x}
def add(*aa):
    c={}
    for a in aa:
        for i,x in a.items(): c[i]=c.get(i,0)+x
    return {i:x for i,x in c.items() if x}
def shift(a,k,scale=1): return {i+k:scale*x for i,x in a.items() if x}
def star(a,b): return {i:x*b[i] for i,x in a.items() if b.get(i,0)}
def val(a,t): return sum(F(x)*t**i for i,x in a.items())
def beta(B,x,y):
    d=B-2
    return {j:choose(x+y-2*j,x+d*j) for j in range(-(x//d),y//B+1)}
def lifted(G,m,t):
    q=max(G)
    return mul({i:comb(m,i) for i in range(m+1)},
               {2*(q-i):F(c)*t**i for i,c in G.items()})
def derivative(a): return {i-1:i*x for i,x in a.items() if i}
def coeff_product(a,b,k): return sum(x*b.get(k-i,0) for i,x in a.items())
def polylist(a): return [F(a.get(i,0)) for i in range(max(a)+1)]
def mod_reduce(a,p):
    a=list(a)
    while len(a)>=len(p):
        x=a[-1]/p[-1]; d=len(a)-len(p)
        for j,c in enumerate(p): a[d+j]-=x*c
        while a and a[-1]==0:a.pop()
    return a+[F(0)]*(len(p)-1-len(a))
def root_remainder(a,p):
    n=len(p)-1
    powers={-1:[-p[j+1]/p[0] for j in range(n)],0:[F(1)]+[F(0)]*(n-1)}
    for j in range(1,max(a)+1):powers[j]=mod_reduce([F(0)]+powers[j-1],p)
    return [sum(F(c)*powers[j][i] for j,c in a.items()) for i in range(n)]
def interval(a,lo,hi):
    low=high=F(0)
    for c in reversed(a):
        pp=(low*lo,low*hi,high*lo,high*hi)
        low,high=min(pp)+c,max(pp)+c
    return low,high
def root_signs(P,responses):
    p=polylist(P); rr=[root_remainder(a,p) for a in responses]
    expr=sum(c*T**j for j,c in P.items())
    for digits in (5,10,20,40,80):
        ii=S.polys.polytools.intervals(expr,eps=S.Rational(1,10**digits))
        if all(interval(r,F(lo),F(hi))[1]<0 for (lo,hi),_ in ii for r in rr):
            need(len(ii)==max(P) and all(m==1 and hi<0 for (lo,hi),m in ii),'complete simple negative root bank')
            for (lo,hi),m in ii:
                for r in rr:need(interval(r,F(lo),F(hi))[1]<0,'exact response negativity')
            return len(ii)
    raise RuntimeError('sign bank unresolved')
def main():
    maps=0; phases=(F(-1),F(-2),F(-1,3))
    for B in range(3,9):
      for h in range(1,5):
       for x in (0,1,2,5):
        for r in (0,B-1):
         for z in (0,1):
            y=B*h+r; m=x+y+z; L=x//(B-2); q=L+h; k=2*h+z
            b=beta(B,x,y); G=shift(b,L)
            a={j:choose(m,z+2*j) for j in range((m-z)//2+1)}
            aa={j:choose(2*m,2*z+2*j) for j in range(-z,m-z+1)}
            P=star(a,b); hit=star(aa,mul(b,b))
            need(min(G)==0 and max(G)==q and G[0]>0 and G[q]>0,'full shifted support')
            need(0<k<m+2*q,'strict coefficient interior')
            for t in phases:
                H=lifted(G,m,t)
                need(H[0]==G[q]*t**q and max(H)==m+2*q,'exact carrier ends')
                need(H.get(k,0)==t**L*val(P,t),'first same-carrier coefficient')
                need(coeff_product(H,H,2*k)==t**(2*L)*val(hit,t),'completed alpha square coefficient')
            maps+=1
    coupled=0
    for h in range(1,7):
      for x in (1,2,5):
        q=x+h; g=x+3*h+1
        b,c,d=(beta(3,x,3*h-i) for i in range(3))
        G,C,D=(shift(v,x) for v in (b,c,d))
        need(add(mul(b,b),shift(mul(c,d),1,2))==beta(3,2*x,6*h),'full beta midpoint split')
        aa={j:choose(2*g,2+2*j) for j in range(-1,g)}
        skip=star(aa,shift(mul(c,d),1,2))
        for t in phases:
            KB={2*(q-i):F(v)*t**i for i,v in G.items()}
            KC={2*(q-1-i):F(v)*t**i for i,v in C.items()}
            KD={2*(q-1-i):F(v)*t**i for i,v in D.items()}
            rhs=shift(add(shift(KC,0,q+2),shift(derivative(KC),1)),1,F(2,3))
            need(derivative(KB)==rhs,'first contiguous transformed carrier identity')
            need(add(shift(KD,0,q+1),shift(derivative(KD),1))==add(shift(KC,0,2),shift(derivative(KC),1,F(3,2))),'second contiguous transformed carrier identity')
            common={i:comb(g,i) for i in range(g+1)}
            HC,HD=mul(common,KC),mul(common,KD)
            need(val(skip,t)==2*t**(1-2*x)*coeff_product(HC,HD,4*h),'single mixed-product residual coefficient')
        coupled+=1
    rows=[];roots=0
    for h in range(5,11):
      for x in (1,3,4,7,13,31):
        g=x+3*h+1
        if gcd(g,6*h+3)!=1:continue
        b,c,d=(beta(3,x,3*h-i) for i in range(3))
        O={j:choose(g,2*j+1) for j in range((g-1)//2+1)}
        E={j:choose(g,2*j) for j in range(g//2+1)}
        OO,EE=mul(O,O),shift(mul(E,E),-1)
        bb,cd=mul(b,b),shift(mul(c,d),1,2)
        P=star(O,b)
        groups=[star(EE,bb),star(OO,cd),star(EE,cd)]
        hit=star(add(OO,EE),bb); skip=add(groups[1],groups[2])
        count=root_signs(P,groups+[hit,skip])
        roots+=count;rows.append((h,x,g,count))
    b,c,d=(beta(3,1,3-i) for i in range(3))
    O={0:5,1:10,2:1}; E={0:1,1:10,2:5}; t=F(-2)
    need([val(star(a,v),t) for a,v in ((O,b),(O,c),(O,d),(E,b))]==[0,15,10,-16],'different contiguous roots remain invalid')
    # Quadratic pullback of 1+t at positive phase is 1+u^2: not real-rooted.
    need(lifted({0:1,1:1},0,F(1))=={0:F(1),2:F(1)},'negative phase assumption has a real-root boundary')
    print('ANALYTIC_THEOREM A=2 complete-alpha beta-hit response is strictly negative at every first-row zero')
    print('OPEN full beta-skip response; no universal actual doubled-row conclusion')
    print('COEFFICIENT_MAPS',maps,'parameter rows,',len(phases),'exact phases each; nonprimitive and x=0 included')
    print('COUPLED_MAPS',coupled,'A2B3 rows; two contiguous identities and single mixed-product response')
    print('FINITE_SIGNS',len(rows),'primitive indexed rows;',roots,'complete roots;',5*roots,'response signs')
    print('BANK',json.dumps(rows,separators=(',',':')))
    print('SAME_CONSTRAINT_HOSTILE h1g5rho-2: OstarB=0,OstarC=15,OstarD=10,EstarB=-16')
    print('PASS',CHECKS,'explicit gates')
    print('semantic_sha256',hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest())
if __name__=='__main__':main()
