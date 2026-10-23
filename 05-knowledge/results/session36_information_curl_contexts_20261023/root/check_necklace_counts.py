#!/usr/bin/env python3
import itertools, json, math
from fractions import Fraction
from collections import Counter
from check_curl_decoder import matrix,key,check,encode

def oddparts(n,lo=1):
    if n==0: yield ()
    for k in range(lo,n+1,2):
        for rest in oddparts(n-k,k): yield (k,)+rest

def tournament_types(n):
    ans=Fraction(0)
    for p in oddparts(n):
        counts=Counter(p)
        e=sum(k//2 for k in p)+sum(math.gcd(p[i],p[j]) for i in range(len(p)) for j in range(i))
        den=math.prod(k**m*math.factorial(m) for k,m in counts.items())
        ans+=Fraction(2**e,den)
    check(ans.denominator==1,'Burnside noninteger')
    return int(ans)

def phi(n): return sum(math.gcd(k,n)==1 for k in range(1,n+1))
def mobius(n):
    sign=1;p=2
    while p*p<=n:
        if n%p==0:
            n//=p;sign=-sign
            if n%p==0:return 0
        p+=1
    return -sign if n>1 else sign

def permute(M,p): return [[M[i][j] for j in p] for i in p]
def canon(M): return min(encode(permute(M,p)) for p in itertools.permutations(range(len(M))))
def scc(M):
    n=len(M);R=[[i==j or M[i][j]==1 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):R[i][j]|=R[i][k] and R[k][j]
    blocks=list({tuple(j for j in range(n) if R[i][j] and R[j][i]) for i in range(n)})
    blocks.sort(key=lambda b:-sum(M[b[0]][v]==1 for v in range(n) if v not in b))
    return blocks

def main():
    N=12;a=[tournament_types(n) for n in range(N+1)]
    ell=[Fraction(0)]*(N+1)
    for n in range(1,N+1): ell[n]=a[n]-sum(Fraction(j,n)*ell[j]*a[n-j] for j in range(1,n))
    c=[0]*(N+1)
    for n in range(1,N+1):
        v=sum(Fraction(phi(d),d)*ell[n//d] for d in range(1,n+1) if n%d==0)
        check(v.denominator==1,'necklace noninteger');c[n]=int(v)
    prim=[0]+[sum(mobius(d)*c[n//d] for d in range(1,n+1) if n%d==0) for n in range(1,N+1)]
    even=[0]*(N+1)
    for n in range(1,N+1):
        v=ell[n]-sum(Fraction(d,n)*even[d] for d in range(1,n) if n%d==0)
        check(v.denominator==1,'Euler noninteger');even[n]=int(v)
    check(prim==even,'primitive vs Euler mismatch')
    check(a[:9]==[1,1,1,2,4,12,56,456,6880],'known A000568 mismatch')
    check(even[1:9]==[1,0,1,2,8,43,398,6413],'recovered connected even sequence mismatch')
    rows=[]
    for n in range(1,6):
        curls=set();types=set();seen=set();stabilizer_checks=0
        for b in range(1<<math.comb(n,2)):
            M=matrix(n,b);t=canon(M)
            if t in seen:continue
            seen.add(t);types.add(t)
            perms=list(itertools.permutations(range(n)))
            K=key(M); orbit=[key(permute(M,p)) for p in perms];curls.add(min(orbit))
            blocks=scc(M);letters=[];aut=1
            for B in blocks:
                S=[[M[i][j] for j in B] for i in B]
                letters.append((len(B),canon(S)))
                aut*=sum(permute(S,p)==S for p in itertools.permutations(range(len(B))))
            r=sum(letters[i:]+letters[:i]==letters for i in range(len(letters)))
            check(sum(x==K for x in orbit)==r*aut,'curl stabilizer formula')
            stabilizer_checks+=1
        check(len(types)==a[n] and len(curls)==c[n],'direct orbit census mismatch')
        rows.append(dict(n=n,tournament_types=len(types),curl_types=len(curls),stabilizers_checked=stabilizer_checks))
    print(json.dumps(dict(status='PASS',A000568=a,curl_types=c[1:],primitive_curl_types=prim[1:],connected_even=even[1:],direct_checks=rows),sort_keys=True,indent=2))
if __name__=='__main__':main()
