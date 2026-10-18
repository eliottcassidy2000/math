#!/usr/bin/env python3
"""Exact whole-polynomial check of the reflected Franel / Apery companion.
Use direct integer Lucas blocks, homogeneous polynomial substitution, and
existing mechanism(). No formal-series truncation is used in the comparisons.
"""
import argparse,json
from pathlib import Path
from ordinary_detector_and_mechanism import mechanism
from sporadic_cm_detector import primes_between


def need(b,msg):
    if not b:raise RuntimeError(msg)
def trim(a):
    while len(a)>1 and not a[-1]:a.pop()
    return a
def add(a,b,p):
    c=a+[0]*max(0,len(b)-len(a))
    for i,v in enumerate(b):c[i]=(c[i]+v)%p
    return trim(c)
def scale(a,b,p):return trim([v*b%p for v in a])
def mul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):c[i+j]=(c[i+j]+u*v)%p
    return trim(c)
def hom(a,N,D,p):
    R=[a[-1]];power=[1]
    for u in a[-2::-1]:
        power=mul(power,D,p);R=add(mul(R,N,p),scale(power,u,p),p)
    return R
def evaluate(a,x,p):
    t=0
    for u in a[::-1]:t=(t*x+u)%p
    return t
def seq(N,cube=False):
    a=[1,5 if cube else 2]
    for n in range(1,N):
        num=((2*n+1)*(17*n*n+17*n+5)*a[-1]-n**3*a[-2]) if cube else ((7*n*n+7*n+2)*a[-1]+8*n*n*a[-2])
        q,r=divmod(num,(n+1)**(3 if cube else 2));need(r==0,'integer recurrence');a.append(q)
    return a
def block(u,p):
    F=[a%p for a in u[:p]]
    v=[]
    for k in range(p):
        q,r=divmod(u[p+k]-u[1]*u[k],p);need(r==0,'Lucas quotient');v.append(q*pow(u[1],-1,p)%p)
    c=v[0];G=[(v[k]-c*F[k])%p for k in range(p)]
    return F,G,c

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bound',type=int,default=397);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    fr=seq(2*a.bound);au=seq(2*a.bound,True);rows=[]
    for p in primes_between(7,a.bound):
        H,G,g0=block(fr,p);A,J,j0=block(au,p)
        eps=1 if p%6==1 else -1
        sigN=[pow(8,-1,p),p-1];D=[1,1];tN=[0,1,p-8]
        GH=hom(G,sigN,D,p);HH=hom(H,sigN,D,p)
        need(trim(H[:])==scale(HH,eps,p),'Franel reflection')
        C=eps*evaluate(G,pow(8,-1,p),p)%p
        B=add(add(G,scale(GH,eps,p),p),scale(H,-C,p),p)
        exact_period=hom(A,tN,D,p);exact_companion=hom(J,tN,D,p)
        need(exact_period==mul(H,H,p),'full period pullback')
        need(exact_companion==mul(H,B,p),'full companion pullback')
        need(B[0]==0,'normalized constant')
        mech=mechanism(au,p,True);need(mech['rad|V'],'existing mechanism')
        rows.append({'p':p,'epsilon':eps,'period_degree':len(trim(A))-1,'companion_degree':len(trim(J))-1,'pullback_degree':len(exact_companion)-1,'Franel_companion_constant_removed':g0,'Apery_companion_constant_removed':j0,'normalizing_c':C,'mechanism':mech})
    out={'status':'FINITE-EXACT whole-polynomial identities; all-prime proof independently audited separately','prime_count':len(rows),'bound':a.bound,'rows':rows}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='rows'})
if __name__=='__main__':main()
