#!/usr/bin/env python3
"""Independent determinant expansion vs subset path DP; no assertion statements."""
import itertools, json, random, hashlib
from pathlib import Path

def check(ok, msg):
    if not ok: raise RuntimeError(msg)

def rows(n, code):
    a=[[0]*n for _ in range(n)]
    for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
        a[i][j]=(code>>bit)&1; a[j][i]=1-a[i][j]
    return a

def det(a, names):
    s=0
    for p in itertools.permutations(names):
        if all(a[i][j] for i,j in zip(names,p)):
            inv=sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
            s+=(-1)**inv
    return s

def trial(a):
    n=len(a); m=1<<n
    num=[det(a,[i for i in range(n) if mask>>i&1]) for mask in range(m)]
    den=[v*(-1)**mask.bit_count() for mask,v in enumerate(num)]
    q=[0]*m; q[0]=1
    for mask in range(1,m):
        v=num[mask]; sub=mask
        while sub:
            v-=den[sub]*q[mask^sub]; sub=(sub-1)&mask
        q[mask]=v
    hp=[[0]*n for _ in range(m)]; totals=[1]+[0]*(m-1)
    for mask in range(1,m):
        for j in range(n):
            if not mask>>j&1: continue
            rest=mask^(1<<j)
            hp[mask][j]=sum(hp[rest][i]*a[i][j] for i in range(n)) if rest else 1
        totals[mask]=sum(hp[mask])
    zeta=q[:]
    for j in range(n):
        for mask in range(m):
            if mask>>j&1: zeta[mask]+=zeta[mask^(1<<j)]
    check(zeta==totals, 'square-zero determinant/HP disagreement')
    # Explicit simple-walk enumeration independently checks the endpoint
    # resolvent interpretation, against dynamic programming at full support.
    endpoint=[[0]*n for _ in range(n)]
    for word in itertools.permutations(range(n)):
        if all(a[i][j] for i,j in zip(word,word[1:])):
            endpoint[word[0]][word[-1]]+=1
    check([sum(endpoint[i][j] for i in range(n)) for j in range(n)]==hp[-1],
          'endpoint walk/resolvent mismatch')
    check(all(x>=0 for x in q), 'negative cycle-gas coefficient')
    return totals[-1]

def main():
    counts={}; hashes=hashlib.sha256()
    for n in range(1,6):
        counts[str(n)]=1<<(n*(n-1)//2)
        for code in range(counts[str(n)]):
            h=trial(rows(n,code)); hashes.update(f'{n},{code},{h};'.encode())
    rng=random.Random(291)
    for n in (6,7):
        counts[str(n)+'_sample']=80
        for _ in range(80): trial(rows(n,rng.getrandbits(n*(n-1)//2)))
    checks=0
    for p in (2,3,5):
        for r in (1,2):
            length=p**r
            if length>9: continue
            # Mark each of two path choices; count primitive cyclic words by
            # number of occurrences of the second path. Independent orbit census.
            orbit={}
            for word in itertools.product((0,1),repeat=length):
                shifts={word[k:]+word[:k] for k in range(length)}
                if len(shifts)==length: orbit[min(shifts)]=sum(word)
            actual=[0]*(length+1)
            for weight in orbit.values(): actual[weight]+=1
            from math import comb
            expected=[]
            for j in range(length+1):
                numerator=comb(length,j)-(comb(length//p,j//p) if j%p==0 else 0)
                check(numerator%length==0,'nonintegral necklace')
                expected.append(numerator//length)
            check(actual==expected,'marked necklace mismatch'); checks+=1
    for h in range(1,100,2):
        for p in (3,5,7,11):
            for r in (1,2,3):
                a=h**(p**r)-h**(p**(r-1)); d=p**r
                check(a%d==0 and (a//d)%2==0,'parity obstruction')
    print(json.dumps({'status':'PASS','tournament_cases':counts,
                      'marked_necklace_orbit_cases':checks,
                      'exhaustive_h_digest':hashes.hexdigest(),
                      'scope':'all labelled n<=5; 80 seeded samples each n=6,7; not coverage evidence'},indent=2,sort_keys=True))
if __name__=='__main__': main()
