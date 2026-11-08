#!/usr/bin/env python3
"""Exact recursive join/singleton-cycle compiler; no external data/dependencies."""
from itertools import combinations, permutations
from math import comb, factorial
from fractions import Fraction
from hashlib import sha256
import json

X=4095
N=2+(X-1).bit_length()-1

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def val(p,t):return sum(c*t**i for i,c in enumerate(p))
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return tuple(out)
def shift(p,s):
    return tuple(sum(p[j]*comb(j,i)*s**(j-i) for j in range(i,len(p))) for i in range(len(p)))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def cyc(a,b):
    u=mul(sub(shift(a,1),a),sub(shift(b,1),b))
    v=mul(sub(a,shift(a,-1)),sub(b,shift(b,-1)))
    return (0,)+tuple(x+y+z for x,y,z in zip(mul(a,b),u,v))
def collar(p):return cyc(p,(0,1))
def tournaments(n):
    pairs=list(combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        r=[0]*n
        for k,(i,j) in enumerate(pairs):
            a,b=(i,j) if bits>>k&1 else (j,i);r[a]|=1<<b
        yield tuple(r)
def join(a,b):
    n=len(a);bm=((1<<len(b))-1)<<n
    return tuple(x|bm for x in a)+tuple(x<<n for x in b)
def cycle(a,b):
    # A -> B -> z -> A.
    n,m=len(a),len(b);z=n+m;bm=((1<<m)-1)<<n
    return tuple(x|bm for x in a)+tuple((x<<n)|(1<<z) for x in b)+((1<<n)-1,)
def hsub(r):
    n=len(r);size=1<<n;dp=[[0]*n for _ in range(size)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,size):
        for i,z in enumerate(dp[mask]):
            if not z:continue
            todo=r[i]&~mask
            while todo:
                bit=todo&-todo;todo-=bit
                dp[mask|bit][bit.bit_length()-1]+=z
    h=list(map(sum,dp));h[0]=1
    return h
def qvalues(r,k):
    h=hsub(r);size=len(h);q=[0]*size;q[0]=1;out=[0]
    for _ in range(k):
        p=[0]*size
        for mask in range(size):
            s=mask
            while True:
                p[mask]+=h[s]*q[mask^s]
                if not s:break
                s=(s-1)&mask
        q=p;out.append(q[-1])
    return out
def permutation_poly(r):
    n=len(r);F=[0]*(n+1)
    for p in permutations(range(n)):
        b=sum(not(r[x]>>y&1) for x,y in zip(p,p[1:]))
        for k in range(b+1,n+1):F[k]+=comb(n-1-b,k-1-b)
    out=[0]*(n+1);fall=(1,)
    for d in range(1,n+1):
        need(F[d]%factorial(d)==0,'unordered cover divisibility')
        fall=mul(fall,(-(d-1),1))
        for i,c in enumerate(fall):out[i]+=c*(F[d]//factorial(d))
    return tuple(out)
def rank(A):
    A=[[Fraction(x) for x in r] for r in A];k=0
    for j in range(len(A[0])):
        at=next((i for i in range(k,len(A)) if A[i][j]),None)
        if at is None:continue
        A[k],A[at]=A[at],A[k];c=A[k][j];A[k]=[v/c for v in A[k]]
        for i in range(len(A)):
            if i!=k:
                c=A[i][j];A[i]=[a-c*b for a,b in zip(A[i],A[k])]
        k+=1
    return k

def compile_profiles(seeds,limit):
    # Values are (minimum cyclic depth, canonical expression).
    bank=[{} for _ in range(limit+1)]
    def put(n,p,d,e):
        if val(p,1)>X:return
        if p not in bank[n] or (d,e)<bank[n][p]:bank[n][p]=(d,e)
    for name,rows,p in seeds:
        if len(rows)<=limit:put(len(rows),p,0,name)
    for n in range(2,limit+1):
        for i in range(1,n//2+1):
            for a,(da,ea) in sorted(bank[i].items()):
                for b,(db,eb) in sorted(bank[n-i].items()):put(n,mul(a,b),max(da,db),'J('+ea+','+eb+')')
        for i in range(1,(n-1)//2+1):
            for a,(da,ea) in sorted(bank[i].items()):
                for b,(db,eb) in sorted(bank[n-1-i].items()):put(n,cyc(a,b),1+max(da,db),'C('+ea+','+eb+')')
    return bank

def graph_expr(expr,seedrows):
    if expr[0]=='s':return seedrows[expr]
    op=expr[0];body=expr[2:-1];depth=0
    for i,c in enumerate(body):
        if c=='(':depth+=1
        if c==')':depth-=1
        if c==',' and depth==0:break
    a=graph_expr(body[:i],seedrows);b=graph_expr(body[i+1:],seedrows)
    return join(a,b) if op=='J' else cycle(a,b)

def main():
    distinct={};seed_checks=0;allsmall=[]
    for n in range(1,6):
        for r in tournaments(n):
            p=permutation_poly(r)
            need(qvalues(r,n)==[val(p,j) for j in range(n+1)],'independent seed profile')
            distinct.setdefault((n,p),r);seed_checks+=1
            if n<=3:allsmall.append((r,p))
    seeds=[('s'+str(i),r,p) for i,((n,p),r) in enumerate(sorted(distinct.items()))]
    seedrows={name:r for name,r,p in seeds}
    binary_checks=0
    for a,p in allsmall:
        for b,q in allsmall:
            r=cycle(a,b);v=cyc(p,q)
            need(qvalues(r,len(r))==[val(v,j) for j in range(len(r)+1)],'literal binary transform')
            need(val(v,1)==2*val(p,1)*val(q,1)+(val(p,2)-val(p,1))*(val(q,2)-val(q,1)),'scalar law')
            need(val(v,1)>=(2**len(a)-1)*(2**len(b)-1)+2>=2**(len(r)-2)+1,'order bound')
            binary_checks+=1
    for n in range(1,5):
        for r in tournaments(n):
            p=permutation_poly(r);v=collar(p);out=cycle(r,(0,))
            need(qvalues(out,len(out))==[val(v,j) for j in range(len(out)+1)],'literal unary transform')
    bank=compile_profiles(seeds,N)
    generators={val(p,1):0 for name,r,p in seeds if val(p,1)<=X}
    for a in range(1,N-1):
        for b in range(a,N-a):
            for p,(dp,ep) in bank[a].items():
                for q,(dq,eq) in bank[b].items():
                    v=2*val(p,1)*val(q,1)+(val(p,2)-val(p,1))*(val(q,2)-val(q,1))
                    if v<=X:generators[v]=min(generators.get(v,999),1+max(dp,dq))
    depth={1:0}
    for v in range(3,X+1,2):
        ds=[generators[v]] if v in generators else []
        for g,dg in generators.items():
            if 1<g<v and v%g==0 and v//g in depth:ds.append(max(dg,depth[v//g]))
        if ds:depth[v]=min(ds)
    # A genuine physical same-(H,Q2), different-Q3 hostile.
    p=(0,0,0,28,0,56,0,14,0,1);q=(0,0,0,36,0,46,0,16,0,1)
    witnesses=[]
    for v in (p,q):
        d,e=bank[9][v];r=graph_expr(e,seedrows)
        need(len(r)==9 and qvalues(r,3)==[val(v,j) for j in range(4)],'hostile realization')
        r1=cycle(r,(0,));r2=cycle(r1,(0,))
        need(hsub(r1)[-1]==4419,'one collar equal')
        need(hsub(r2)[-1]==val(collar(collar(v)),1),'two collar separation')
        witnesses.append({'polynomial':v,'rows':r,'expression':e,'depth':d,'Q1_Q2_Q3':[val(v,j) for j in (1,2,3)],'two_collar_H':hsub(r2)[-1]})
    need(witnesses[0]['two_collar_H']!=witnesses[1]['two_collar_H'],'moment hostile')
    rankchecks=[]
    for k in range(9):
        A=[]
        for depthk in range(k+1):
            row=[]
            for n in range(1,k+2):
                p=(0,)*n+(1,)
                for _ in range(depthk):p=collar(p)
                row.append(val(p,1))
            A.append(row)
        need(rank(A)==k+1,'physical response rank')
        rankchecks.append(k+1)
    rows=[1];momentrows=[rows]
    for k in range(1,9):
        nxt=[0]*(len(rows)+1)
        for j,c in enumerate(rows,1):
            nxt[j]+=j*c;nxt[j-1]+=j*j*c
            if j>=2:nxt[j-2]-=j*c
        rows=nxt;need(rows[-1]==factorial(k),'top moment coefficient');momentrows.append(rows)
    recursive_witnesses=[]
    for target in (241,4001):
        candidates=[(d,n,e,p) for n,B in enumerate(bank) for p,(d,e) in B.items() if val(p,1)==target]
        d,n,e,p=min(candidates);r=graph_expr(e,seedrows)
        need(hsub(r)[-1]==target and d==depth[target],'recursive depth witness')
        recursive_witnesses.append({'H':target,'depth':d,'order':n,'expression':e,'rows':r})
    padded=join((2,4,1),(0,))
    need(len(padded)==4 and hsub(padded)[-1]==3,'outer join violates naive strong order bound')
    need(19 not in depth,'seed-five noncoverage hostile')
    encoded=json.dumps(sorted(depth.items()),separators=(',',':')).encode()
    print(json.dumps({'status':'PASS','cutoff':X,'cyclic_order_bound':N,'labelled_seed_checks':seed_checks,'distinct_seed_profiles':len(seeds),'literal_binary_checks':binary_checks,'literal_unary_checks':75,'profile_counts_by_order':[len(x) for x in bank[1:]],'scalar_count':len(depth),'scalar_depth_counts':{str(k):sum(d==k for d in depth.values()) for k in sorted(set(depth.values()))},'smallest_new_at_depth':{str(k):min(v for v,d in depth.items() if d==k) for k in sorted(set(depth.values()))},'first_missing':sorted(set(range(1,X+1,2))-depth.keys())[:40],'scalar_depth_sha256':sha256(encoded).hexdigest(),'recursive_witnesses':recursive_witnesses,'hostile_witnesses':witnesses,'physical_response_ranks':rankchecks,'moment_rows':momentrows,'seed_representatives':[{'name':name,'rows':r,'polynomial':p} for name,r,p in seeds]},sort_keys=True,indent=2))
if __name__=='__main__':main()
