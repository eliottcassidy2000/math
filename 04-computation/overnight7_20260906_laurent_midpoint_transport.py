#!/usr/bin/env python3
"""Exact A2/B3 three-group midpoint transport and named sign controls.

The identities are proved for every h>=1,x>=1 in the report. Finite group
sign controls do not prove the universal remaining signed transport.
"""
from math import comb,factorial,gcd
from fractions import Fraction as F
import hashlib,json,sys
import sympy as S
sys.stdout.reconfigure(newline='\n')
T=S.symbols('tau'); CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)
def convolution(a,b):
    ans={}
    for i,x in a.items():
        for j,y in b.items():ans[i+j]=ans.get(i+j,0)+x*y
    return {i:x for i,x in ans.items() if x}
def star(a,b):return {i:x*b[i] for i,x in a.items() if i in b and x*b[i]}
def add(*args):
    ans={}
    for a in args:
        for i,x in a.items():ans[i]=ans.get(i,0)+x
    return {i:x for i,x in ans.items() if x}
def shift(a,d,scale=1):return {i+d:scale*x for i,x in a.items()}
def beta(x,y):return {j:comb(x+y-2*j,x+j) for j in range(-x,y//3+1)}
def expression(a):return sum(v*T**j for j,v in a.items())
def beta_path_counts(U,V,H):
    # A class is recorded at the unique crossing of X+3Y=H.
    # 1=hit, 2=jump from H-1, 3=jump from H-2.
    grid=[[[0]*4 for _ in range(V+1)] for _ in range(U+1)]
    grid[0][0][0]=1
    for x in range(U+1):
        for y in range(V+1):
            for dx,dy in ((1,0),(0,1)):
                xx,yy=x+dx,y+dy
                if xx>U or yy>V:continue
                before,after=x+3*y,xx+3*yy
                for state,n in enumerate(grid[x][y]):
                    new=state
                    if state==0 and after==H:new=1
                    elif state==0 and before<H<after:new=1+H-before
                    grid[xx][yy][new]+=n
    return grid[U][V]
def interval_value(poly,lo,hi):
    low=high=F(0)
    for coefficient in S.Poly(poly,T).all_coeffs():
        products=(low*lo,low*hi,high*lo,high*hi)
        low,high=min(products)+F(coefficient),max(products)+F(coefficient)
    return low,high
def exact_negative_groups(P,groups):
    inverse=S.invert(T,P,T)
    remainders=[S.rem(S.expand(T*group)*inverse,P,T) for group in groups]
    for digits in (4,8,16,32,64):
        intervals=S.polys.polytools.intervals(P,eps=S.Rational(1,10**digits))
        if all(interval_value(r,F(lo),F(hi))[1]<0
               for (lo,hi),multiplicity in intervals for r in remainders):
            return intervals,remainders
    raise RuntimeError('named grouped response sign not certified')
def main():
    manifest=[]; pathcases=0
    for h,g in ((1,5),(1,7),(1,11),(2,8),(2,11),(2,16),
                (3,11),(3,13),(3,17),(4,14),(4,16),(4,19)):
        x,y=g-3*h-1,3*h
        need(x>=1 and gcd(g,6*h+3)==1,'named source admissibility')
        O={i:comb(g,2*i+1) for i in range((g-1)//2+1)}
        E={i:comb(g,2*i) for i in range(g//2+1)}
        B,C,D=beta(x,y),beta(x,y-1),beta(x,y-2)
        BB=convolution(B,B); CD=shift(convolution(C,D),1)
        doubled=beta(2*x,2*y)
        need(add(BB,shift(CD,0,2))==doubled,'two skipped-level classes give exact beta doubling')
        skip1=convolution(beta(x,y-1),beta(x-1,y+1))
        skip2=convolution(beta(x,y-2),beta(x-1,y+2))
        need(skip1==skip2==CD,'full-support skip class reindexing')
        for j,count in doubled.items():
            path=beta_path_counts(2*y-3*j,2*x+j,y+3*x)
            need(path==[0,BB.get(j,0),CD.get(j,0),CD.get(j,0)],'independent classified path dynamic program')
            need(sum(path)==count,'all doubled beta paths retained')
            pathcases+=1
        OO=convolution(O,O); EE=shift(convolution(E,E),-1)
        actual_alpha={j:comb(2*g,2+2*j) for j in range(-1,g)}
        need(add(OO,EE)==actual_alpha,'exact alpha parity Vandermonde partition')
        virtual=star(OO,BB)
        groups=[star(EE,BB),star(OO,shift(CD,0,2)),star(EE,shift(CD,0,2))]
        actual=star(actual_alpha,doubled)
        raw={j:factorial(2*g)//(factorial(2*x+j)*factorial(2*y-3*j)*factorial(2+2*j))
             for j in range(-1,2*h+1)}
        need(actual==raw==add(virtual,*groups),'three complete actual-minus-virtual response groups')
        P=expression(star(O,B))
        intervals,rems=exact_negative_groups(P,[expression(group) for group in groups])
        need(len(intervals)==h and all(m==1 and hi<0 for (lo,hi),m in intervals),'complete first-root coverage')
        for (lo,hi),m in intervals:
            for r in rems:need(interval_value(r,F(lo),F(hi))[1]<0,'exact grouped negative sign interval')
        # Contiguous operator identities on the same full source coefficients.
        for j,b in B.items():
            N=x+y-2*j
            need(N>=2,'Euler denominator has no support pole')
            need(N*C.get(j,0)==(y-3*j)*b,'first contiguous Euler identity')
            need((N-1)*D.get(j,0)==(y-1-3*j)*C.get(j,0),'second contiguous Euler identity')
        oe,ee=expression(O),expression(E)
        need(S.expand(g*ee-(1+(g-1)*T)*oe-2*T*(1-T)*S.diff(oe,T))==0,'same-g even/odd derivative identity')
        manifest.append({'h':h,'g':g,'root_count':h,'groups':3})
    O={0:5,1:10,2:1}; E={0:1,1:10,2:5}
    B,C,D=beta(1,3),beta(1,2),beta(1,1)
    values=[expression(star(A,b)).subs(T,-2) for A,b in ((O,B),(O,C),(O,D),(E,B))]
    need(values==[0,15,10,-16],'contiguous replacements do not preserve the original root condition')
    print('PROVED_IDENTITY beta(2x,2y)=B^2+2*tau*C*D; alpha_double=O^2+tau^-1*E^2')
    print('THREE_GROUPS (tau^-1*E^2)starB^2; 2O^2star(tau*C*D); 2(tau^-1*E^2)star(tau*C*D)')
    print('UNIFORM_GROUPED_NEGATIVITY OPEN; named controls only',json.dumps(manifest,sort_keys=True))
    print('SAME_CONSTRAINT_HOSTILE h1g5rho-2: OstarB=0,OstarC=15,OstarD=10,EstarB=-16')
    print('PASS',CHECKS,'explicit gates;',pathcases,'independent beta path endpoints')
    print('semantic_sha256',hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest())
if __name__=='__main__':main()
