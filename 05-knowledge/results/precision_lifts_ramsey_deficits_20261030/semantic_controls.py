#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import gcd
import json
PTH=Path(__file__).parent
def need(x,m):
    if not x:raise RuntimeError(m)
def height(n,E):
    if n==0:return E
    h=0
    while h<E and n%5==0:h+=1;n//=5
    return h
def lam(e):return F(5,8)*(1-F(1,5**e))
def run():
    words=0;transitions=0
    def dfs(P,B,Q,A):
        nonlocal words,transitions
        for a in range(1,17-A):
            pp,bb,qq=3*P,3*B+Q,Q*2**a;words+=1
            for E in range(1,5):
                q=5**E;z=(B-P+Q)*pow(Q,-1,q)%q;m=P*pow(Q,-1,q)%q
                h=height(z,E);small=5**max(1,h)
                z2=(bb-pp+qq)*pow(qq,-1,q)%q;m2=pp*pow(qq,-1,q)%q
                h2=height(z2,E);newmod=5**max(1,h2)
                need(z2==(3*z+2**a-2)*pow(2**a,-1,q)%q,'z transition')
                need(m2%min(small,newmod)==3*(m%small)*pow(2**a,-1,min(small,newmod))%min(small,newmod),'compatible lift')
                need(min(h2,height(m2%newmod-1,E))==height(gcd(pp-qq,bb),E),'actual charge')
                transitions+=1
            dfs(pp,bb,qq,A+a)
    dfs(1,0,1,0)
    phase_checks=0
    for E in range(1,10):
        q=5**E
        for a in range(2,65):
            h=height((1-pow(2,1-a,q))%q,E);ell=max(1,h);T=2*5**(ell-1)
            for s in range(T):
                short=sum(s//(2*5**(j-1)) for j in range(1,E+1))-lam(E+2)*s
                for b in (0,1,2,5,17,1001):
                    k=b*T+s
                    long=sum(k//(2*5**(j-1)) for j in range(1,E+1))-lam(E+2)*k
                    need(long<=short,'least paid phase');phase_checks+=1
    points=list(product(range(5),repeat=2));phi=lambda x:F((-x[0])%5,5)
    for x in points:
        for y in points:need(int(y==(0,0))+phi(y)-phi(x)<=1,'arbitrary boundary fee')
        y=((x[0]+4)%5,(x[1]+4)%5)
        need(int(y==(0,0))+phi(y)-phi(x)<=F(1,5),'one boundary fee')
    audit=json.loads((PTH/'adaptive9_audit.json').read_text())
    eta=F(-audit['terminal'],audit['D']);beta=F(-audit['max'],audit['D']);gamma=min(eta-2,beta)
    families=[(a,b) for b,aa in ((1,range(4)),(2,range(2)),(3,range(1))) for a in aa]
    need(all(b*gamma>=1-F(a,5) for a,b in families),'seven E12 families')
    w=(2,2,3,1,2,2,1,1,2,1,4);P,B,Q,A,S=1,0,1,0,0;charges=[]
    for a in w:
        before=height(B-P+Q,9);P,B,Q,A=3*P,3*B+Q,Q*2**a,A+a
        c=height(gcd(P-Q,B),9);S+=c;charges.append(c)
    need((B-P+Q)==2*5**9 and before==0 and height(P-Q,9)==1,'artificial lift witness')
    need(lam(11)*A-(S-charges[-1]+9)==eta,'abstract endpoint bound attained')
    # Necessary clock condition for any non-one native12-return.
    q=5**12;pow2=1;pow3=1;seen={};first=None
    for A in range(1,50001):
        pow2=2*pow2%q
        if A>=2:
            pow3=3*pow3%q;seen.setdefault(pow3,A-1)
        if pow2 in seen:first=(A,seen[pow2]);break
    need(first is not None,'clock probe stopping limit')
    delta=gamma-1+F(first[0],2*5**12)
    need(delta>F(1,8),'next-price lower bound')
    next_families=[dict(nonone=b,ones=a) for b in range(1,9) for a in range(5) if b*delta+F(a,5)<1]
    out=dict(status='FINITE-EXACT semantic controls; next clock bound is necessary only',
             words_clock_at_most_16=words,adaptive_transition_checks=transitions,
             initial_phase_checks=phase_checks,boundary_checks=650,
             closed_E12_families=[dict(ones=a,nonone=b) for a,b in families],
             critical_artificial_lift=dict(word=w,P=P,B=B,Q=Q,actual_charges=charges,actual_total=S,
                                          actual_terminal_content=charges[-1],artificial_total=S-charges[-1]+9,price=str(eta)),
             first_multiplier12_return_clock=first[0],corresponding_length=first[1],
             resulting_E13_native_price_floor=str(delta),remaining_E13_families=next_families)
    (PTH/'semantic_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('critical_artificial_lift','remaining_E13_families')},indent=2))
if __name__=='__main__':run()
