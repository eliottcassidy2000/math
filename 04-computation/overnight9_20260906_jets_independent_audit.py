#!/usr/bin/env python3
"""Independent exact one-digit Deuring-loss referee; standard library only.

Polynomial reflection uses Horner composition; reciprocal coefficients solve
the inverse defining equation. Full observer controls use local pivots.
No producer, SymPy, or integer Smith engine is imported.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
import hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)
def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def add(a,b):return trim([(a[j] if j<len(a) else 0)+(b[j] if j<len(b) else 0) for j in range(max(len(a),len(b)))])
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)
def power(a,n):
    out=[1]
    for _ in range(n):out=mul(out,a)
    return out
def deriv(a):return [j*a[j] for j in range(1,len(a))] or [0]
def evaluate(a,x,mod=None):
    out=0
    for c in reversed(a):
        out=out*x+c
        if mod:out%=mod
    return out
def reflection(a):
    out=[0]
    for c in reversed(a):out=add(mul(out,[1,-1]),[c])
    return out
def operator(a,k):return add(add(mul([0,1,-1],deriv(deriv(a))),mul([1,-2],deriv(a))),scale(a,k*(k+1)))
def valuation(x,p):
    x=F(x)
    if not x:return float('inf')
    a,b=abs(x.numerator),x.denominator;v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v
def inverse(nodes,m):
    answer=[]
    for i,x in enumerate(nodes):
        q=[1]
        for j,y in enumerate(nodes):
            if i!=j:q=mul(q,[x-y,1])
        a=power(q,m);b=[1/F(a[0])]
        for r in range(1,m):b.append(-sum(a[j]*b[r-j] for j in range(1,r+1))/a[0])
        answer.append(b)
    return answer
def smith(nodes,m,p):
    N=len(nodes)*m
    a=[[F(comb(j,r)*x**(j-r) if j>=r else 0) for j in range(N)] for x in nodes for r in range(m)]
    result=[]
    for k in range(N):
        v,i,j=min((valuation(a[i][j],p),i,j) for i in range(k,N) for j in range(k,N))
        need(v>=0 and v!=float('inf'),'finite p-integral Smith pivot')
        a[k],a[i]=a[i],a[k]
        for row in a:row[k],row[j]=row[j],row[k]
        pivot=a[k][k];result.append(v)
        for i in range(k+1,N):
            ratio=a[i][k]/pivot
            need(valuation(ratio,p)>=0,'p-integral row elimination')
            for j in range(k+1,N):a[i][j]-=ratio*a[k][j]
            a[i][k]=F(0)
        for j in range(k+1,N):a[k][j]=F(0)
    need(result==sorted(result),'ordered actual Smith exponents')
    need(sum(result)==m*m*sum(valuation(x-y,p) for x,y in combinations(nodes,2)),'full confluent determinant')
    return result
def family(p):
    k=(p-1)//2
    # Recover the coefficients from their hypergeometric adjacent ratio.
    P=[comb(2*k,k)]
    for j in range(k):
        next=F(P[-1]*(k+j+1)*(k-j),(j+1)*(2*k-j))
        need(next.denominator==1,'integral coefficient from ratio recurrence')
        P.append(int(next))
    need(P==[comb(k+j,j)*comb(2*k-j,k-j) for j in range(k+1)],'ratio versus factorial numerator')
    R=scale(reflection(P),(-1)**k)
    diff=add(P,scale(R,-1))
    need(all(c%p==0 for c in diff) and len(diff)<=k,'integral lower-degree divided reflection')
    G=[c//p for c in diff]
    need(operator(P,k)==scale(deriv(P),p),'integer P ODE by coefficients')
    need(operator(R,k)==scale(deriv(R),-p),'integer reflected ODE')
    need(operator(G,k)==add(deriv(P),deriv(R)),'integer inhomogeneous divided ODE')
    W=add(mul(P,deriv(G)),scale(mul(deriv(P),G),-1))
    residue=add(add(mul([0,1,-1],W),scale(mul(P,P),-1)),[1])
    need(len(residue)<=p and all(c%p==0 for c in residue),'complete low-degree Wronskian constant')
    H=[comb(k,j)**2 for j in range(k+1)]
    need(all((P[j]-(-1)**k*H[j])%p==0 for j in range(k+1)),'correct Deuring sign and scaling')
    A=scale(H,(-1)**k)
    # Independent exact adjacent polynomial from its two chosen powers.
    B=[0]+[(-1)**(k+1)*comb(k,j-1)*comb(k,j) for j in range(1,k+1)]
    need(scale(B,k+1)==scale(mul([0,1],add(scale(A,k),mul([1,-1],deriv(A)))),-1),'exact adjacent coefficient derivative identity')
    if p<=31:
        for a in range(-1,k+2):
            f=power([0,a,-a-1,1],k)
            need(evaluate(A,a)==f[2*k] and evaluate(B,a)==f[2*k-1],'direct cubic expansion gives both coefficients')
    return P,R,G,H,A,B
def main():
    manifest=[];lifts=0;exactjets=0;matrices=0
    for p in (3,5,7,11,13,17,19,23,31,43,59,83,103):
        k=(p-1)//2;m=k+1
        P,R,G,H,A,B=family(p)
        roots=[a for a in range(2,p) if evaluate(H,a,p)==0]
        for a in roots:
            need(a*(1-a)*evaluate(deriv(P),a,p)*evaluate(G,a,p)%p==1,'root forces derivative and companion units')
            need(evaluate(B,a,p)!=0,'adjacent polynomial unit at every residue root')
        if p<=31:
            for a in range(2,p):
                delta=int(a in roots)
                for digit in range(p):
                    x=a+p*digit
                    need(min(valuation(evaluate(P,x),p),valuation(evaluate(R,x),p))==delta,'all ordinary and root residue lifts have exact joint cap')
                    lifts+=1
            ordinary=next((a for a in range(2,p) if a not in roots),None)
            examples=sorted(set(([roots[0],roots[0]+p*p,-1] if roots else [])+([ordinary] if ordinary else [])))
            for a in examples:
                delta=int(evaluate(H,a,p)==0)
                jets=inverse((0,1,a),m)
                need(jets[0][k]==F(evaluate(P,a),a**p),'exact first-node numerator normalization')
                need(jets[1][k]==F(evaluate(R,a),(1-a)**p),'exact reflected-node normalization')
                third=sum(c*a**(k-j)*(a-1)**j for j,c in enumerate(P))
                need(jets[2][k]==F((-1)**k*third,(a*(a-1))**p),'exact third-node numerator normalization')
                need(min(valuation(row[k],p) for row in jets)==delta,'all three actual top jets joint cap')
                if delta:need(all(valuation(row[k-1],p)==0 for row in jets),'all three actual adjacent jets unit')
                for e in (0,1,2,3):
                    nodes=(0,p**e,a*p**e)
                    cost=max(-valuation(c,p) for row in inverse(nodes,m) for c in row if c)
                    want=0 if e==0 else (3*m-1)*e-delta
                    need(cost==want,'literal rational full-jet maximum including shallow mask')
                    exactjets+=1
                if p<=19 and a<p and a>=2:
                    for e in (1,2):
                        actual=smith((0,p**e,a*p**e),m,p)
                        need(actual[-1]==(3*m-1)*e-delta,'independent full observer largest exponent')
                        matrices+=1
        if p==3:need(P==[2,2] and R==[-4,2] and G==[2],'ternary boundary constants')
        manifest.append((p,k,len(roots)))
    print('POLYNOMIALS',json.dumps(manifest,separators=(',',':')))
    print('ALL_P2_LIFTS',lifts,'ordinary and supersingular admissible residue fibers')
    print('ACTUAL_RECIPROCAL_OBSERVERS',exactjets,'all jet orders and four depths')
    print('FULL_LOCAL_SMITH',matrices,'fresh actual Hasse matrices')
    print('SCOPE exact equilateral largest loss only; e0 separate; no close or full-partition extrapolation')
    print('semantic_sha256',hashlib.sha256(json.dumps(manifest,separators=(',',':')).encode()).hexdigest())
    print('PASS',CHECKS,'explicit gates; no primary mathematical imports')
if __name__=='__main__':main()
