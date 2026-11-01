#!/usr/bin/env python3
"""Literal independent audit of the authenticated moving-apex map.
All Hamiltonian words and ordered covers are counted explicitly; no producer imports.
State order throughout: h,F,q,g. Graphs are adjacency rows of booleans.
"""
import itertools,json

def require(c,msg):
    if not c: raise RuntimeError(msg)

def tournament(n,code):
    a=[[False]*n for _ in range(n)]; bit=0
    for i in range(n):
        for j in range(i+1,n):
            a[i][j]=bool(code>>bit&1);a[j][i]=not a[i][j];bit+=1
    return a

def add_apex(a,cut):
    n=len(a); b=[row[:]+[not bool(cut>>i&1)] for i,row in enumerate(a)]
    b.append([bool(cut>>i&1) for i in range(n)]+[False]);return b

def hp(a):
    n=len(a)
    return sum(all(a[w[i]][w[i+1]] for i in range(n-1)) for w in itertools.permutations(range(n)))

def marked(a,cut):
    n=len(a);h=q=g=0
    for w in itertools.permutations(range(n)):
        good=[a[w[i]][w[i+1]] for i in range(n-1)]
        if all(good):h+=1;g+=bool(cut>>w[-1]&1)
        for k in range(n):
            if not(cut>>w[k]&1):continue
            if all(good[:max(0,k-1)]) and all(good[k:]):q+=1
    return h,hp(add_apex(a,cut)),q,g

def positive(state):
    h,F,q,g=state
    require((F-h)%2==0 and (q-F+h-g)%2==0,'integrality')
    C=(F-h)//2;L=(q-F+h-g)//2;t=h-g
    require(min(C,L,g,t)>=0,'nonnegative actual state')
    return C,L,g,t

def singleton(a,cut,bit):
    n=len(a); b=[row[:]+[True] for row in a]+[[False]*(n+1)]
    return b,cut|(bit<<n)

def migrate(a,cut):
    # New parent A+x, new apex z beats exactly x.
    return add_apex(a,cut),1<<len(a)

def migrate2(a,cut):
    b,c=singleton(a,cut,0)
    # Parent order A,z2,x; new root z1 beats z2 and x.
    return add_apex(b,c),(1<<len(a))|(1<<(len(a)+1))

def main():
    cases=0; tests=0;by_n={};hostiles=[]
    for n in range(1,5):
        local=0
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code)
            for cut in range(1<<n):
                old=marked(a,cut);h,F,q,g=old;C,L,g2,t=positive(old)
                require(g==g2,'g mismatch')
                r=marked(*migrate(a,cut))
                require(r==(F,F+q+g,q+h,h-g),f'R mismatch {n,code,cut,old,r}')
                require(positive(r)==(C+L+g,0,t,2*C+g),'positive R mismatch')
                r2=marked(*migrate2(a,cut))
                require(r2==(F+q+g,F+3*q+g,F+3*q+2*h,F+g),'R2 mismatch')
                require(positive(r2)==(2*C+2*L+g,C+L+g+t,2*C+2*g+t,2*C+2*L+g),'positive R2 mismatch')
                for bit in (0,1):
                    new=marked(*singleton(a,cut,bit))
                    expected=(h,F+q+g,2*q,0) if bit==0 else (h,F,2*q+h,h)
                    require(new==expected,'singleton mismatch')
                    pos=(2*C+L+g,L,0,g+t) if bit==0 else (C,C+2*L+g,g+t,0)
                    require(positive(new)==pos,'positive singleton mismatch')
                    tests+=2
                # Hostile: dropping the empty suffix Q=(x) removes exactly h.
                require(r[2]-q==h and h>0,'empty suffix invoice missing')
                if n==1:hostiles.append({'parent_code':code,'cut':cut,'old':old,'R':r})
                cases+=1;local+=1;tests+=5
        by_n[n]=local
    print(json.dumps({'status':'PASS','marked_cases':cases,'by_parent_order':by_n,'identity_checks':tests,'maps_checked':['S0','S1','R','R2'],'count_method':'literal vertex permutations and all legal ordered two-path cover cuts','hostiles':hostiles,'scope':'Exhaustive parent orders 1..4, all labelled tournaments and all marked cuts. General proof needs path bijection, not this finite census.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
