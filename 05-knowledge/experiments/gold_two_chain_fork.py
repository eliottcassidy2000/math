#!/usr/bin/env python3
"""Exact controls for the symbolic two-chain/fork Gold classification."""
from fractions import Fraction
from itertools import product
from math import comb, isqrt
from pathlib import Path
import json
from gold_terminal_chain import count, require, deficits, gold
from gold_four_collar import table_pred

def C(n,k):return comb(n,k) if 0<=k<=n else 0

def core(m,n,H=()):
    T=m+n;h=len(H)
    p=[(1<<i)-1 for i in range(m)]
    p += [((1<<j)-1)<<m for j in range(n)]
    p += [(x<<T)|((1<<T)-1) for x in H]
    p += [0,1|(1<<(T+h))]
    return tuple(p),(0,T+h+1),(T-1,m-1),(T-2,m-1)

def triangular(a,b,total):return [[a,0],[b,total-a-b]]

def formula(m,n,h=0):
    T=m+n;N=T+h;K=(N+1)*(N+2)//2
    f0=triangular(K*C(T-1,m-1)-C(T+1,m+1),K*C(T-2,m-1)-C(T,m+1),K*C(T,m)-C(T+2,m+2))
    base=triangular(C(T-2,m-2),C(T-3,m-2),C(T-1,m-1))
    fb=triangular(C(T-1,m-1),C(T-2,m-1),C(T,m))
    scale=lambda t,c:[[c*x for x in row] for row in t]
    return [f0,scale(base,N*(N+1)//2),scale(fb,N+1),scale(base,N)]

def direct_flags(p,ports,e,f):
    F=[];a,b=ports
    for af,bl in ((False,False),(True,False),(False,True),(True,True)):
        pp=list(p)
        if af:
            for i in range(len(p)):
                if i!=a:pp[i]|=1<<a
        if bl:pp[b]|=((1<<len(p))-1)^(1<<b)
        F.append(table_pred(tuple(pp),e,f))
    return F

def side_core(m,n,H,D,W,port):
    base,_,e,f=core(m,n,H);base=list(base[:-2]);N=len(base)
    base += [(x<<N)|(1 if W>>i&1 else 0) for i,x in enumerate(D)]
    return tuple(base),(0,N+port),e,f

def upsets(D):
    for W in range(1<<len(D)):
        if all(not(W>>i&1) or all(not(D[j]>>i&1) or W>>j&1 for j in range(len(D))) for i in range(len(D))):yield W

def margin_formula(m,n):
    T=m+n;d=(T-1)*(T-2)
    return min(Fraction(n-m,T),Fraction((m-1)*(m-2),d),Fraction(m*n+(m-1)*(m-2)-n*n,d))

def margin(F):return min(Fraction(x,sum(sum(t,[]))) for t in F for x in sum(deficits(t),[]))

def main():
    direct=0;symbolic=0;admissible=[]
    caps=[(),(0,),(0,1),(0,0),(0,0,1),(0,1,3)]
    for m,n,H in product(range(2,9),range(2,9),caps):
        p,ports,e,f=core(m,n,H);got=direct_flags(p,ports,e,f)
        expected=[[[count(H)*x for x in row] for row in t] for t in formula(m,n,len(H))]
        require(got==expected,('closed count',m,n,H))
        require(margin(got)==margin_formula(m,n),('sharp margin',m,n,H))
        direct+=1
    for m,n in product(range(2,81),repeat=2):
        F=formula(m,n);q=n*n-m*n-(m-1)*(m-2)
        require(all(gold(t) for t in F)==(m<=n and q<=0),('classification',m,n))
        require(margin(F)==margin_formula(m,n),('cone minimum',m,n))
        for t in F:
            prob=Fraction(t[0][0],sum(sum(t,[])))
            require(Fraction(m-1,m+n-1)<=prob<=Fraction(m,m+n),('balance interval',m,n))
        if m<=n and q<=0 and m<=12:admissible.append([m,n,str(margin_formula(m,n))])
        symbolic+=1
    collars=[]
    for r,s in product((0,1,2,10,30),repeat=2):
        B,ports,e,f=core(4,5);nb=len(B)
        A=tuple((1<<i)-1 for i in range(r));Ccap=tuple((1<<i)-1 for i in range(s))
        p=list(A)+[(x<<r)|((1<<r)-1) for x in B]
        p += [(x<<(r+nb))|((1<<(r+nb))-1) for x in Ccap]
        p[r+ports[0]]&=~((1<<r)-1)
        for j in range(s):p[r+nb+j]&=~(1<<(r+ports[1]))
        got=table_pred(tuple(p),tuple(x+r for x in e),tuple(x+r for x in f))
        F=formula(4,5)
        expected=[[F[0][i][j]+r*F[1][i][j]+s*F[2][i][j]+r*s*F[3][i][j] for j in range(2)] for i in range(2)]
        require(got==expected,'unbounded collar kernel')
        require(margin([got])>=Fraction(1,56),'uniform 1/56')
        collars.append(dict(left=r,right=s,vertices=len(p),margin=str(margin([got]))))
    pell=[];L0,L1=2,1;F0,F1=0,1
    for k in range(1,20):
        L0,L1=L1,L0+L1;F0,F1=F1,F0+F1
        if k%4==3:
            require((L0+6)%5==0,'Pell integrality')
            m=(L0+6)//5;require((F0+m)%2==0,'Pell parity');n=(F0+m)//2
            require(n*n-m*n-(m-1)*(m-2)==0,'Pell boundary')
            pell.append(dict(k=k,m=m,n=n))
    side_checks=0
    side_shapes=[(0,),(0,1),(0,0),(0,0,1),(0,1,3),(0,0,1,3),(0,0,1,2)]
    for (m,n),H,D in product(((2,2),(3,3),(4,4),(4,5),(5,6)),((),(0,0)),side_shapes):
        for W in upsets(D):
            for p in range(len(D)):
                if any(x>>p&1 for x in D):continue
                pp,ports,e,f=side_core(m,n,H,D,W,p);F=direct_flags(pp,ports,e,f)
                require(margin(F)>=margin_formula(m,n),'arbitrary side-poset margin')
                for table in F:
                    prob=Fraction(table[0][0],sum(sum(table,[])))
                    require(Fraction(m-1,m+n-1)<=prob<=Fraction(m,m+n),'arbitrary side-poset balance')
                side_checks+=1
    recurrence=[(2,2)]
    while recurrence[-1][0]<=10000:
        m,n=recurrence[-1];recurrence.append((2*m+3*n-3,3*m+5*n-6))
    direct_boundary=[]
    for m in range(2,10001):
        d=5*m*m-12*m+8;y=isqrt(d)
        if y*y==d and (m+y)%2==0:direct_boundary.append((m,(m+y)//2))
    require(direct_boundary==[v for v in recurrence if v[0]<=10000],'Pell recurrence exhaustion control')
    matmul=lambda A,B:[[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    M=[[-2,1],[-1,1]];B=[[0,1],[-1,1]];Bi=[[1,-1],[1,0]]
    require(matmul(matmul(Bi,M),B)==[[0,-1],[-1,-1]],'budget matrix conjugacy')
    out=dict(status='symbolic classification and sharp cone margin proved; exact controls passed',actual_core_cap_checks=direct,formula_parameter_checks=symbolic,small_admissible_parameters=admissible,eleven_vertex_core=core(4,5)[0],eleven_vertex_flag_tables=formula(4,5),sharp_unbounded_context_margin='1/56',unbounded_chain_collar_checks=collars,pell_boundary=pell,boundary_recurrence_control=direct_boundary,budget_matrix_conjugacy_verified=True,arbitrary_side_poset_checks=side_checks)
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
