#!/usr/bin/env python3
from itertools import combinations, permutations
from collections import Counter

def tournament(m):
    n=2*m+2
    R={(2*i,2*i+3) for i in range(m)}
    a=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if (i,j) in R:a[j][i]=True
            else:a[i][j]=True
    return a

def dp_h(a):
    n=len(a); D=[[0]*n for _ in range(1<<n)]
    for i in range(n): D[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,num in enumerate(D[mask]):
            if num:
                for j in range(n):
                    if not mask>>j&1 and a[i][j]: D[mask|(1<<j)][j]+=num
    return sum(D[-1])

def odd_cycles(a):
    n=len(a); out=[]
    for k in range(3,n+1,2):
        for U in combinations(range(n),k):
            for p in permutations(U[1:]):
                seq=(U[0],)+p
                if all(a[seq[i]][seq[(i+1)%k]] for i in range(k)):
                    out.append(seq)
    return out

def confirm(condition,msg):
    if not condition:raise RuntimeError(msg)

def check_cycles(m,a):
    cycles=odd_cycles(a); observed=Counter()
    R={(2*i+3,2*i):i for i in range(m)}
    for seq in cycles:
        edges=list(zip(seq,seq[1:]+seq[:1])); ids=sorted(R[e] for e in edges if e in R)
        confirm(bool(ids),'cycle missing backwards arc')
        lo,hi=ids[0],ids[-1]
        confirm(ids==list(range(lo,hi+1)),'nonconsecutive')
        endpoints={x for i in ids for x in (2*i,2*i+3)}
        cap=set(seq)-endpoints
        confirm(cap in [{2*lo+1},{2*hi+2}],'wrong support')
        observed[(lo,hi,next(iter(cap)))]+=1
    expected=Counter((lo,hi,c) for lo in range(m) for hi in range(lo,m) for c in (2*lo+1,2*hi+2))
    confirm(observed==expected, f'cycle enumeration m={m} mismatch')
    return cycles

def packing_count(cycles,n):
    # Count each collection once, using independent include/skip recursion.
    D={0:1}
    for seq in cycles:
        mask=sum(1<<i for i in seq)
        E=D.copy()
        for used,w in D.items():
            if not used&mask:E[used|mask]=E.get(used|mask,0)+2*w
        D=E
    return sum(D.values())

def main():
    h=[1,5,17]
    for m in range(3,17): h.append(3*h[-1]+h[-2]+h[-3])
    for m in range(0,7):
        a=tournament(m); measured=dp_h(a)
        confirm(measured==h[m],f'path DP m={m}: {measured}!={h[m]}')
        if m<=4:
            cycles=check_cycles(m,a)
            confirm(packing_count(cycles,len(a))==measured,f'packing m={m}')
        if m:
            # Explicit cycle union implies strength; independently test reachability.
            for u in range(len(a)):
                R={u}
                for _ in a:R|={v for x in R for v in range(len(a)) if a[x][v]}
                confirm(len(R)==len(a),'not strong')
    s=(1,5,17%7); start=s; vals=[]
    while True:
        vals.append(s[0]);s=(s[1],s[2],(3*s[2]+s[1]+s[0])%7)
        if s==start:break
        confirm(len(vals)<=24,'period too long')
    confirm(len(vals)==24 and all(vals),'mod 7 claim failed')
    print('H_0..H_10:',h[:11]);print('mod-7 period:',vals)
    print('PASS: DP m<=6; independently enumerated cycles/packings m<=4; strength m<=6; modulo-7 period')
if __name__=='__main__':main()
