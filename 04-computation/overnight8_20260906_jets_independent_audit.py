#!/usr/bin/env python3
"""Independent standard-library referee: p-local pivots, minors, residue jets.

No producer, SymPy, integer Smith routine or closed binomial inverse-jet
formula is imported. Exact inverse series use polynomial coefficient recursion.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb,factorial
from pathlib import Path
import hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)
def vp(x,p):
    x=F(x)
    if not x:return float('inf')
    a,b=abs(x.numerator),x.denominator;v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c
def power(a,n):
    out=[1]
    for _ in range(n):out=mul(out,a)
    return out
def evaluate(a,x):return sum(c*x**j for j,c in enumerate(a))
def det(a):
    a=[row[:] for row in a];prev=1;sign=1
    for k in range(len(a)-1):
        i=next((i for i in range(k,len(a)) if a[i][k]),None)
        if i is None:return 0
        if i!=k:a[i],a[k]=a[k],a[i];sign=-sign
        v=a[k][k]
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                n=a[i][j]*v-a[i][k]*a[k][j]
                need(n%prev==0,'integer Bareiss division')
                a[i][j]=n//prev
        for i in range(k+1,len(a)):a[i][k]=0
        prev=v
    return sign*a[-1][-1]
def interpolate(values):
    delta=list(map(F,values));out=[F(0)]*len(values);fall=[F(1)]
    for k in range(len(values)):
        for j,c in enumerate(fall):out[j]+=delta[0]*c/factorial(k)
        delta=[delta[j+1]-delta[j] for j in range(len(delta)-1)]
        fall=mul(fall,[-k,1])
    while out and out[-1]==0:out.pop()
    need(all(c.denominator==1 for c in out),'interpolated integral determinant polynomial')
    return list(map(int,out))
def smith(nodes,m,p):
    n=len(nodes)*m
    a=[[F(comb(j,r)*x**(j-r) if j>=r else 0) for j in range(n)] for x in nodes for r in range(m)]
    exponents=[]
    for k in range(n):
        v,i,j=min((vp(a[i][j],p),i,j) for i in range(k,n) for j in range(k,n))
        need(v!=float('inf') and v>=0,'finite p-integral pivot')
        a[k],a[i]=a[i],a[k]
        for row in a:row[k],row[j]=row[j],row[k]
        pivot=a[k][k];exponents.append(v)
        for i in range(k+1,n):
            factor=a[i][k]/pivot
            need(vp(factor,p)>=0,'unimodular p-local row multiplier')
            for j in range(k+1,n):a[i][j]-=factor*a[k][j]
            a[i][k]=F(0)
        for j in range(k+1,n):a[k][j]=F(0)
    need(exponents==sorted(exponents),'ordered local Smith list')
    need(sum(exponents)==m*m*sum(vp(x-y,p) for x,y in combinations(nodes,2)),'confluent determinant valuation')
    return exponents
def inverse_coefficients(nodes,m,p=None):
    answer=[]
    for x in nodes:
        Q=[1]
        for y in nodes:
            if y!=x:Q=mul(Q,[x-y,1])
        A=power(Q,m)
        if p:
            inv=pow(A[0]%p,-1,p);b=[inv]
            for r in range(1,m):b.append(-sum(A[j]*b[r-j] for j in range(1,r+1))*inv%p)
        else:
            b=[1/F(A[0])]
            for r in range(1,m):b.append(-sum(A[j]*b[r-j] for j in range(1,r+1))/A[0])
        answer.append(b)
    return answer
def predicted(nodes,p):
    v=sorted(vp(x-y,p) for x,y in combinations(nodes,2));e=v[0];d=v[-1]-e
    if d:return 11*e+7*d-int(p==5)
    if not e:return 0
    a=[(x-nodes[0])//p**e for x in nodes]
    AP=p==7 and any((2*a[i]-a[j]-a[k])%7==0 for i,j,k in ((0,1,2),(1,0,2),(2,0,1)))
    return 11*e-int(p in (3,5) or AP)
def min_weight():
    aa=[0,1];am=[-1,1]
    factors={
      (3,):[4],(7,):[0,4],(3,7):[40*x for x in mul(aa,am)],
      (2,3,7):[200*x for x in mul(mul(aa,am),[-1,2])],
      (3,6,7):[200*x for x in mul(mul(power(aa,4),am),[-2,1])],
      (2,3,6,7):[2100*x for x in mul(power(aa,4),power(am,4))],
      (1,2,3,6,7):[560*x for x in mul(mul(power(aa,4),power(am,4)),[2,-7,7])],
      (2,3,5,6,7):[560*x for x in mul(mul(power(aa,9),power(am,4)),[7,-7,2])],
      (1,2,3,5,6,7):[1680*x for x in mul(power(aa,9),power(am,9))]}
    seen={};counts=[]
    for r,w in enumerate((1,3,7,12,19,27),1):
        count=0
        for I in combinations(range(8),r):
            for J in combinations(range(4,12),r):
                W=sum(J)-sum(i%4 for i in I);count+=1
                need(W>=w,'whole rank1..6 minor weight bound')
                if W!=w:continue
                need(I in factors and J==tuple(range(4,4+r)),'complete equality weight class')
                near=[i for i in I if i>=4]
                degree=sum(sorted(J,reverse=True)[:len(near)])-sum(i%4 for i in near)
                def value(a):return det([[comb(j,i%4)*(a**(j-i%4) if i>=4 else 1) for j in J] for i in I])
                P=interpolate([value(a) for a in range(degree+1)])
                need(P==factors[I],'independently reconstructed factorization')
                need(evaluate(P,-3)==value(-3),'extra exact determinant check')
                seen[I]=P
        counts.append(count)
    need(counts==[64,784,3136,4900,3136,784] and len(seen)==9,'entire weight universe')
    for a in range(2,7**2):
        if a%7 in (0,1):continue
        for r in range(1,7):need(min(vp(evaluate(P,a),7) for I,P in seen.items() if len(I)==r)==int(r>=4),'joint polynomial floor and attainment')
    print('MINOR_AUDIT 12804 weight classes; 9 independently interpolated exact polynomials')
def literal_controls():
    rows=[]
    for p in (3,5,7,11,13):
        for e,d,u in ((0,1,1),(1,1,-1),(2,1,p+2),(1,3,1),(3,2,-1)):
            nodes=(0,p**e,p**(e+d)*u)
            s=smith(nodes,4,p)
            need(s[-1]==predicted(nodes,p),'fresh close-case largest precision')
            rows.append((p,nodes,s))
        for a in range(2,p):
            for e in (1,2):
                nodes=(0,p**e,a*p**e)
                s=smith(nodes,4,p)
                need(s[-1]==predicted(nodes,p),'complete residue equilateral precision')
                if p==7:
                    k=int(a in (2,4,6))
                    need(s==[0,0,0,0,e,2*e,4*e,5*e+1,7*e,8*e,10*e-1+k,11*e-k],'entire p7 equilateral partition')
                rows.append((p,nodes,s))
    for a,e in ((2,0),(3,0),(9,3),(10,4),(46,3),(47,2)):
        s=smith((0,7**e,a*7**e),4,7);k=int(a%7 in (2,4,6))
        want=[0]*12 if e==0 else [0,0,0,0,e,2*e,4*e,5*e+1,7*e,8*e,10*e-1+k,11*e-k]
        need(s==want,'p7 shallow mask and fresh higher lifts')
        rows.append((7,(0,7**e,a*7**e),s))
    A=smith((0,7,14),4,7);B=smith((0,7,21),4,7)
    need((A[-1],B[-1])==(10,11) and sum(A)==sum(B)==48,'literal odd-prime metric hostile')
    need((sum(min(10,x) for x in A),sum(min(10,x) for x in B))==(48,47),'actual finite observer kernel')
    print('LOCAL_SMITH',len(rows)+2,'full matrices; own p-integral Gaussian elimination')
    print('literal_sha256',hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest())
def residue_controls():
    cases=0;ceiling=0
    for p in (3,5,7,11):
        k=(p-1)//2;m=k+1
        universe=[a for n in range(2,p+1) for a in combinations(range(p),n)] if p<=7 else [(0,1,2),(0,1,3),(0,2,5,7),tuple(range(6)),tuple(range(p))]
        for nodes in universe:
            fn=[1]
            for x in nodes:fn=mul(fn,[-x,1])
            coeff=power(fn,k);packet=coeff[p-1::p]
            need(len(coeff)-1<p*p and len(packet)<len(nodes),'degree and evaluation injectivity bounds')
            jets=inverse_coefficients(nodes,m,p)
            for i,x in enumerate(nodes):
                direct=sum(c*comb(j,p-1)*x**(j-p+1) for j,c in enumerate(coeff) if j>=p-1)%p
                need(direct==evaluate(packet,x)%p,'sparse coefficient translation identity')
                Qi=1
                for y in nodes:
                    if x!=y:Qi=Qi*(x-y)%p
                need(jets[i][k]==direct*pow(Qi,-p,p)%p,'recursive series matches packet coefficient')
            need(any(x%p for x in packet)==any(b[k] for b in jets),'top unit iff nonzero coefficient vector')
            if len(nodes)<=4 and nodes[0]==0 and nodes[-1]<=4:
                actual=smith(tuple(p*x for x in nodes),m,p)[-1]
                need((actual==len(nodes)*m-1)==any(x%p for x in packet),'full observer top-ceiling iff')
                ceiling+=1
            cases+=1
    for p in (3,5,7,11,13,19,23,31):
        k=(p-1)//2;H=[comb(k,j)**2 for j in range(k+1)]
        reflection=[sum(H[j]*comb(j,r)*(-1)**r for j in range(r,k+1))%p for r in range(k+1)]
        need(reflection==[(-1)**k*x%p for x in H] and H==H[::-1],'full S3 generating unit identities')
        roots=[]
        for a in range(2,p):
            f=[0,a,-a-1,1];coeff=power(f,k)
            need(coeff[p-1]%p==(-1)**k*evaluate(H,a)%p,'Deuring coefficient scalar identity')
            jets=inverse_coefficients((0,1,a),k+1,p)
            need((evaluate(H,a)%p==0)==all(b[k]==0 for b in jets),'three-node exact top loss iff')
            if evaluate(H,a)%p==0:roots.append(a)
        if p%4==3:
            need(p-1 in roots,'uniform arithmetic-progression cancellation')
            if p>=7:need(len(roots)<p-2,'both ceiling behaviors occur')
        print('DEURING_CONTROL',p,roots)
    print('PACKETS',cases,'independent polynomial/recursive-series controls;',ceiling,'full observer ceiling checks')
def main():
    min_weight();literal_controls();residue_controls()
    need(hashlib.sha256(Path(__file__).with_name('overnight8_20260906_jets_residue.py').read_bytes()).hexdigest()=='102225c5476d59b3aec02c5b479093d105298000b9411fdb4526f9423efd3b65','pinned primary source')
    print('PASS',CHECKS,'explicit gates; no primary mathematical imports')
if __name__=='__main__':main()
