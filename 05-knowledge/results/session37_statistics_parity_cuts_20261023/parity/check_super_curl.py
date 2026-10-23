#!/usr/bin/env python3
"""Exact SCC-word, free Lie-super span, and labelled signed-curl controls."""
from collections import Counter
from itertools import permutations, combinations
from fractions import Fraction
from math import factorial, gcd
import json

def need(p,msg):
    if not p:raise RuntimeError(msg)

def oddparts(n,lo=1):
    if n==0:yield ()
    for k in range(lo,n+1,2):
        for q in oddparts(n-k,k):yield (k,)+q

def tournament_types(n):
    ans=Fraction(0)
    for q in oddparts(n):
        mult=Counter(q);power=sum(k//2 for k in q)+sum(gcd(q[i],q[j]) for i in range(len(q)) for j in range(i))
        den=1
        for k,m in mult.items():den*=k**m*factorial(m)
        ans+=Fraction(2**power,den)
    need(ans.denominator==1,'integral Burnside');return int(ans)

def cyclic(w):return min(w[i:]+w[:i] for i in range(len(w)))
def weight(w):return sum(x[0] for x in w)

def primitive(w):
    return not any(len(w)%k==0 and w==w[:k]*(len(w)//k) for k in range(1,len(w)))

def signed_cyclic(w):
    n=weight(w);seen={}
    for i in range(len(w)):
        rotated=w[i:]+w[:i];d=weight(w[:i]);s=(-1)**(d*(n-d))
        if rotated in seen and seen[rotated]!=s:return None,0
        seen[rotated]=s
    c=min(seen);return c,seen[c]

def multiply(f,g):
    out=Counter()
    for u,a in f.items():
        for v,b in g.items():out[u+v]+=a*b
    return {w:c for w,c in out.items() if c}

def bracket(f,g,i,j):
    out=Counter(multiply(f,g));s=(-1)**(i*j)
    for w,c in multiply(g,f).items():out[w]-=s*c
    return {w:c for w,c in out.items() if c}

def insert(f,pivots):
    f={w:Fraction(c) for w,c in f.items() if c}
    while f:
        lead=min(f);coefficient=f[lead]
        if lead not in pivots:
            pivots[lead]={w:c/coefficient for w,c in f.items()};return
        for w,c in pivots[lead].items():
            f[w]=f.get(w,0)-coefficient*c
            if not f[w]:del f[w]

def supertrace(f):
    out=Counter()
    for w,c in f.items():
        k,s=signed_cyclic(w)
        if s:out[k]+=s*c
    return {k:c for k,c in out.items() if c}

def perm_sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def matrix(n,mask):
    M=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        M[i][j]=1 if mask>>k&1 else -1;M[j][i]=-M[i][j]
    return M
def curl(M):
    r=len(M)-1
    return tuple(M[i][j]-M[i][r]+M[j][r] for i,j in combinations(range(r),2))

def main():
    top=8;A=[tournament_types(n) for n in range(top+1)];S=[0]*(top+1)
    for n in range(1,top+1):S[n]=A[n]-sum(S[j]*A[n-j] for j in range(1,n))
    words=[[()]]+[[] for _ in range(top)]
    for n in range(1,top+1):
        words[n]=[((k,i),)+w for k in range(1,n+1) for i in range(S[k]) for w in words[n-k]]
        need(len(words[n])==A[n],'word reconstruction')
    C=[0];P=[0];D=[0];rows=[]
    for n in range(1,top+1):
        necklaces={cyclic(w) for w in words[n]};p=sum(primitive(w) for w in necklaces)
        surviving=sum(signed_cyclic(w)[1]!=0 for w in necklaces)
        killed=[]
        for w in necklaces:
            k=next(k for k in range(1,len(w)+1) if len(w)%k==0 and w==w[:k]*(len(w)//k))
            r=len(w)//k;d=weight(w[:k]);bad=r%2==0 and d%2==1
            need(bad==(signed_cyclic(w)[1]==0),'period sign criterion')
            if bad:killed.append(w)
        C.append(len(necklaces));P.append(p);D.append(surviving)
        oddpart=n
        while oddpart%2==0:oddpart//=2
        need(surviving==C[n]-(C[oddpart] if n%2==0 else 0),'closed signed dimension')
        rows.append(dict(n=n,word_count=A[n],strong_generators=S[n],primitive_necklaces=p,curl_necklaces=C[n],signed_curl_dimension=D[n],killed=len(killed)))
    logA=[Fraction(0)]*(top+1)
    for n in range(1,top+1):
        logA[n]=A[n]-sum(Fraction(j,n)*logA[j]*A[n-j] for j in range(1,n))
        burnside=sum(Fraction(sum(gcd(i,d)==1 for i in range(1,d+1)),d)*((-1)**((d-1)*(n//d)))*logA[n//d] for d in range(1,n+1) if n%d==0)
        need(burnside==D[n],('signed Burnside series',n,burnside,D[n]))
    # Full homogeneous super-Lie span: generate all brackets of lower homogeneous bases.
    lie=[[] for _ in range(8)];lie_rows=[];bracket_checks=0
    for n in range(1,8):
        pivots={}
        for i in range(S[n]):insert({((n,i),):1},pivots)
        for i in range(1,n):
            for f in lie[i]:
                for g in lie[n-i]:
                    b=bracket(f,g,i,n-i)
                    need(not supertrace(b),'bracket survives supercocenter')
                    insert(b,pivots);bracket_checks+=1
        lie[n]=list(pivots.values())
        expected=P[n]+(P[n//2] if n%4==2 else 0)
        need(len(lie[n])==expected,('super Lie dimension',n,len(lie[n]),expected))
        lie_rows.append(dict(n=n,dimension=len(lie[n]),parity=n%2))
    # Direct labelled curl S_n/A_n orbit computation, without the SCC classification.
    actual=[]
    for n in range(2,6):
        perms=list(permutations(range(n)));sgns=[perm_sign(p) for p in perms]
        sn=set();an=set();survivors={}
        for mask in range(1<<(n*(n-1)//2)):
            M=matrix(n,mask);base=curl(M)
            orbit=[curl([[M[p[i]][p[j]] for j in range(n)] for i in range(n)]) for p in perms]
            rep=min(orbit);sn.add(rep);an.add(min(k for k,s in zip(orbit,sgns) if s==1))
            admissible=all(s==1 for k,s in zip(orbit,sgns) if k==base)
            need(rep not in survivors or survivors[rep]==admissible,'orbit sign invariance')
            survivors[rep]=admissible
        need(len(sn)==C[n] and sum(survivors.values())==D[n],'labelled curl signed count')
        need(len(an)==C[n]+D[n],'alternating orbit split')
        actual.append(dict(n=n,symmetric_orbits=len(sn),alternating_orbits=len(an),sign_surviving_orbits=sum(survivors.values())))
    print(json.dumps(dict(status='PASS',word_rows=rows,super_lie_rows=lie_rows,bracket_checks=bracket_checks,direct_curl_rows=actual,scope='all SCC words of vertex weight<=8; complete Lie-super spans through7; all labelled tournaments through5 for direct curl signs'),sort_keys=True,indent=2))

if __name__=='__main__':main()
