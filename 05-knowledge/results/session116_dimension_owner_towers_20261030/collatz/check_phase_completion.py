"""Exact five-adic cyclic-subgroup depth and periodic content rates.
All verification gates use explicit exceptions, including under python -O.
"""
from fractions import Fraction as F
from itertools import product
from math import gcd
import json
import random
INF = float('inf')
def require(x, msg):
    if not x: raise ValueError(msg)
def val(x):
    if x == 0: return INF
    x=abs(x); n=0
    while x%5==0: x//=5; n+=1
    return n

def coeff(w):
    P=Q=1; B=0
    for a in w: P,Q,B=3*P,Q*2**a,3*B+Q
    return P,Q,B

def compose(a,b): # chronological a then b
    P,Q,B=a; R,S,C=b
    return P*R,Q*S,R*B+C*Q

def params(W):
    P,Q,B=W; m=P*pow(Q,-1,5)%5
    t=next(t for t in (1,2,4) if pow(m,t,5)==1)
    if t==1: r=min(val(P-Q),val(B))
    else: r=val(P**t-Q**t)
    return t,r

def cutoff(W,H):
    P,Q,B=W; R,S,C=H
    t,r=params(W)
    c=val((P-Q)*C-B*(R-S))
    if t>1:
        x=R*pow(S,-1,5)%5
        if pow(x,t,5)!=1: return 0
        beta=val(R**t-S**t)
        return min(c,beta if beta<r else INF)
    s,h=val(P-Q),val(B)
    alpha=val(R-S) if s<=h else val(C)
    return max(0,min(c-r,alpha if alpha<r else INF))

def mass(t,r,H):
    if H==INF: return F(4*r+1,4*t)
    H=int(H)
    if H<=r: return F(H,t)
    return F(r,t)+F(1,4*t)*(1-F(1,5**(H-r)))

def phase_rate(w,u,v):
    W=coeff(w); H=coeff(v+u)
    t,r=params(W); depth=cutoff(W,H)
    return mass(t,r,depth),depth

def rate(w,u=()):
    return sum((phase_rate(w,u,w[:j])[0] for j in range(1,len(w)+1)),F(0))

def pair(C,q):
    P,Q,B=C; inv=pow(Q,-1,q)
    return P*inv%q,B*inv%q

def orbit(W,q):
    m,b=pair(W,q); x,y=1,0; out=set()
    while (x,y) not in out:
        out.add((x,y)); x,y=m*x%q,(m*y+b)%q
    require((x,y)==(1,0),'orbit failed to close')
    return out

def words(budget,prefix=()):
    if prefix: yield prefix
    for a in range(1,budget+1): yield from words(budget-a,prefix+(a,))

def primitive_root(w):
    return next(w[:d] for d in range(1,len(w)+1) if len(w)%d==0 and w[:d]*(len(w)//d)==w)


def completion_matrix(w):
    return [[phase_rate(w,w[j:],w[:i])[1] for i in range(1,len(w)+1)]
            for j in range(1,len(w)+1)]

def run():
    gates=0; ws=list(words(7)); entries=[()]+list(words(7))
    for w in ws:
        d=len(w); mat=completion_matrix(w); t,r=params(coeff(w))
        rows=[sum((mass(t,r,x) for x in row),F(0)) for row in mat]
        for j in range(d):
            rotation=w[j+1:]+w[:j+1]
            require(rows[j]==rate(rotation),('rotation identity',w,j));gates+=1
            require(mat[j][j]==INF,('diagonal',w,j)); gates+=1
            for i in range(d):
                require(mat[j][i]==mat[i][j],('symmetric',w,i,j)); gates+=1
                if primitive_root(w)==w and i!=j:
                    require(mat[j][i]!=INF,('primitive separation',w,i,j));gates+=1
        for i in range(d):
            for j in range(d):
                for k in range(d):
                    require(mat[i][k]>=min(mat[i][j],mat[j][k]),('ultrametric',w,i,j,k));gates+=1
        for u in entries:
            depths=[phase_rate(w,u,w[:i])[1] for i in range(1,d+1)]
            j=max(range(d),key=lambda j:depths[j]); D=depths[j]
            # Identical successful sets until D; no successes beyond D.
            for e in range(1,5):
                hits=[e<=x for x in depths]; replacement=[e<=x for x in mat[j]]
                require(all(not x or y for x,y in zip(hits,replacement)),('domination',w,u,e));gates+=1
                if e<=D:
                    require(hits==replacement,('prefix agreement',w,u,e));gates+=1
            rr=rate(w,u)
            z=primitive_root(w)
            if z!=w:
                require(rr==(len(w)//len(z))*rate(z,u),('root reduction',w,u));gates+=1
            require(rr<=rows[j]<=max(rows),('rate bound',w,u));gates+=1
            if D!=INF:
                require(rr<rows[j],('strict finite depth',w,u));gates+=1
            if primitive_root(w)==w and rr==max(rows):
                require(D==INF and mat[j][j]==INF,('optimizer infinite',w,u));gates+=1
                require(w[:j+1]+u==w*((j+1+len(u))//d),('word completion',w,u));gates+=1
    # Independent group-orbit check of the formula, without depth arithmetic.
    modular=0
    for w in words(6):
        mat=completion_matrix(w)
        for e in (1,2,3):
            q=5**e; orb=orbit(coeff(w),q)
            for j in range(1,len(w)+1):
                for i in range(1,len(w)+1):
                    direct=pair(coeff(w[:i]+w[j:]),q) in orb
                    require(direct==(e<=mat[j-1][i-1]),('direct matrix',w,e,i,j)); modular+=1
    cases=[]
    for w,expected in [((3,1,2),F(13,16)),((1,2,1,1,2),F(9,20))]:
        rows=[rate(w,w[j:]) for j in range(1,len(w)+1)]
        require(max(rows)==expected and max(rows)>rate(w),('counterexample',w))
        cases.append({'word':w,'coefficients':coeff(w),'parameters':params(coeff(w)),
                      'empty_rate':str(rate(w)),'maximum_rate':str(max(rows)),
                      'winning_phases':[j+1 for j,x in enumerate(rows) if x==max(rows)],
                      'depth_matrix':[['infinity' if x==INF else x for x in row] for row in completion_matrix(w)]})
    receipt={'status':'pass','word_clock_bound':7,'entry_clock_bound':7,
             'controllers':len(ws),'entries_per_controller':len(entries),
             'exact_gates':gates,'independent_modular_gates':modular,'counterexamples':cases,
             'scope':'Finite tests support the separate proof; no orbit convergence claim.'}
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':run()
