#!/usr/bin/env python3
"""Exact independent input-space tests of rooted triangle-curl realization."""
import itertools, json
from collections import defaultdict
from pathlib import Path

def check(x, message):
    if not x: raise RuntimeError(message)

def matrix(n, bits):
    M=[[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(itertools.combinations(range(n),2)):
        M[i][j]=1 if bits>>b&1 else -1; M[j][i]=-M[i][j]
    return M

def key(M):
    r=len(M)-1
    return tuple(M[i][j]-M[i][r]+M[j][r] for i,j in itertools.combinations(range(r),2))

def encode(M):
    return sum((M[i][j]==1)<<b for b,(i,j) in enumerate(itertools.combinations(range(len(M)),2)))

def decode(n,K):
    r=n-1
    if len(K)!=r*(r-1)//2 or any(v not in (-3,-1,1,3) for v in K): return []
    reach=[[i==j for j in range(r)] for i in range(r)]
    plus=set(); minus=set(); entries={}
    for v,(i,j) in zip(K,itertools.combinations(range(r),2)):
        entries[i,j]=v
        if v==3: minus.add(i); plus.add(j)
        elif v==-3: plus.add(i); minus.add(j)
        elif v==1: reach[i][j]=True
        else: reach[j][i]=True
    for k in range(r):
        for i in range(r):
            if reach[i][k]:
                for j in range(r): reach[i][j] |= reach[k][j]
    plus={j for i in plus for j in range(r) if reach[i][j]}
    minus={i for j in minus for i in range(r) if reach[i][j]}
    if plus & minus: return []
    free=set(range(r))-plus-minus; blocks=[]
    while free:
        i=min(free); block={j for j in free if reach[i][j] and reach[j][i]}
        blocks.append(block); free-=block
    snapshot=blocks[:]
    blocks.sort(key=lambda b:-sum(reach[min(b)][min(c)] for c in snapshot))
    # The free quotient is a total order, not an arbitrary partial order.
    for i in range(len(blocks)):
        for j in range(i+1,len(blocks)):
            check(reach[min(blocks[i])][min(blocks[j])] and not reach[min(blocks[j])][min(blocks[i])], 'free quotient not total')
    result=[]
    for cut in range(len(blocks)+1):
        p=plus | set().union(*blocks[cut:]) if blocks[cut:] else plus.copy()
        a=[1 if i in p else -1 for i in range(r)]
        M=[[0]*n for _ in range(n)]
        for i in range(r): M[i][r]=a[i]; M[r][i]=-a[i]
        for (i,j),v in entries.items(): M[i][j]=v+a[i]-a[j]; M[j][i]=-M[i][j]
        check(all(abs(M[i][j])==1 for i,j in itertools.combinations(range(n),2)), 'bad realization')
        check(key(M)==tuple(K), 'curl mismatch')
        result.append(encode(M))
    return sorted(result)

def hp(M):
    n=len(M); d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if d[mask][i]:
                for j in range(n):
                    if not mask>>j&1 and M[i][j]==1: d[mask|1<<j][j]+=d[mask][i]
    return sum(d[-1])

def components(M):
    n=len(M); R=[[i==j or M[i][j]==1 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n): R[i][j] |= R[k][j]
    return len({tuple(j for j in range(n) if R[i][j] and R[j][i]) for i in range(n)})

def main():
    rows=[]
    for n in range(1,7):
        groups=defaultdict(list); H={}; C={}
        for b in range(1<<(n*(n-1)//2)):
            M=matrix(n,b); k=key(M); groups[k].append(b)
            h=hp(M); c=components(M)
            check(k not in H or H[k]==h,'H varies inside fibre')
            check(k not in C or C[k]==c,'SCC count varies')
            H[k]=h; C[k]=c
        for k,vals in groups.items():
            check(decode(n,k)==vals,'decoder not exact')
            check(len(vals)==C[k],'fibre cardinality not SCC count')
        invalid=0
        if n<=5:
            for k in itertools.product((-3,-1,1,3),repeat=(n-1)*(n-2)//2):
                check(decode(n,k)==groups.get(k,[]),'invalid input not rejected')
                invalid+=k not in groups
        rows.append(dict(n=n,tournaments=sum(map(len,groups.values())),curl_classes=len(groups),max_fibre=max(map(len,groups.values())),invalid_inputs_checked=invalid))
    # Pairwise allowed entries alone do not imply realizability.
    hostile=(-3,-3,-3)
    check(abs(hostile[0]+hostile[2]-hostile[1]) in (1,3),'hostile violates local alphabet')
    check(not decode(4,hostile),'hostile unexpectedly realizable')
    payload=dict(status='PASS',universe='All labelled tournaments n=1..6; all alphabet-valid rooted curl arrays n=1..5',rows=rows,unrealizable_n4=list(hostile))
    print(json.dumps(payload,sort_keys=True,indent=2))

if __name__=='__main__': main()
