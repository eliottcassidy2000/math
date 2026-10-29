"""Exact unique-cycle growth experiment. No assert statements are proof gates.

Classical degree-one insertion grammar (Shelor 1977); explicit cycle marks
permit isomorphism canonicalization by rotations alone. Integer two-run
observers are built by subset path DP; FLINT supplies exact determinants.
"""
from pathlib import Path
import argparse
import itertools
import json
import time
import numpy as np
from numba import njit
from flint import fmpz_mat


def require(test, message):
    if not test:
        raise RuntimeError(message)


def bits(A):
    return ''.join(str(int(A[i, j])) for i in range(len(A)) for j in range(i+1, len(A)))


def canonical(A, cycle):
    n = len(A)
    choices = []
    for j in range(n):
        order = cycle[j:] + cycle[:j]
        B = A[np.ix_(order, order)]
        choices.append((bits(B), B))
    return min(choices, key=lambda t:t[0])


def children(A):
    n = len(A)
    top = [v for v in range(n) if sum(A[v]) == n-2]
    require(len(top) in (1,2), 'degree-one grammar scope')
    if len(top) == 1:
        a = top[0]
        ports = [('a', (a-1) % n), ('b', a)]
    else:
        a, b = top
        if not A[a,b]:
            a,b = b,a
        require(b == (a+1) % n, 'maximal vertices adjacent in native cycle')
        ports = [('c',(a-1) % n),('d',a),('e',b)]
    for label, parent in ports:
        B = np.zeros((n+1,n+1), dtype=np.int64)
        B[:n,:n] = A
        B[n,:n] = 1
        B[n,parent] = 0
        B[parent,n] = 1
        cyc = list(range(parent+1)) + [n] + list(range(parent+1,n))
        key, C = canonical(B, cyc)
        yield label, key, C


@njit
def path_table(A):
    n = len(A)
    H = np.zeros((1<<n,n,n), dtype=np.int64)
    for i in range(n):
        H[1<<i,i,i] = 1
    for mask in range(1,1<<n):
        if mask & (mask-1) == 0:
            continue
        for t in range(n):
            if not mask & (1<<t):
                continue
            prev = mask ^ (1<<t)
            for k in range(n):
                if A[k,t] and prev & (1<<k):
                    for s in range(n):
                        H[mask,s,t] += H[prev,s,k]
    return H


@njit
def two_run(H):
    n = H.shape[1]
    full = H.shape[0]-1
    pairs = np.zeros((n*(n-1),2), dtype=np.int64)
    z=0
    for a in range(n):
        for b in range(a+1,n):
            pairs[z,0],pairs[z,1] = a,b
            pairs[z+1,0],pairs[z+1,1] = b,a
            z += 2
    M = np.zeros((z,z), dtype=np.int64)
    for i in range(z):
        a,b = pairs[i]
        for j in range(z):
            c,d = pairs[j]
            for mask in range(1,full):
                if H[mask,a,c] and H[full^mask,b,d]:
                    M[i,j] += H[mask,a,c]*H[full^mask,b,d]
    return M


def exact_stats(A):
    H = path_table(A)
    n = len(A)
    hc = sum(int(H[-1,0,t]) for t in range(1,n) if A[t,0])
    require(hc == 1, 'insertion must preserve exactly one HC')
    M = two_run(H)
    # M has [[X,Y],[Y,X]] blocks. Rational sector change preserves det.
    X,Y = M[::2,::2],M[::2,1::2]
    P,S = X+Y,X-Y
    dp, ds = int(fmpz_mat(P.tolist()).det()), int(fmpz_mat(S.tolist()).det())
    return {'det_plus':dp,'det_minus':ds,'det':dp*ds}, M


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--max-n',type=int,default=11)
    p.add_argument('--output',type=Path,default=Path(__file__).with_name('tournament_results.json'))
    args=p.parse_args()
    require(4 <= args.max_n <= 12, 'int64 audited range 4..12')
    A=np.zeros((4,4),dtype=np.int64)
    for u,v in [(0,1),(1,2),(2,3),(3,0),(0,2),(1,3)]: A[u,v]=1
    key,A=canonical(A,list(range(4)))
    layer={key:('',A)}
    output={'scope':'degree-one insertion grammar; all exact unique-HC classes conditional only on classical completeness theorem', 'layers':{}}
    f0,f1=1,3
    for n in range(4,args.max_n+1):
        expected=f0
        require(len(layer)==expected,'Garey Fibonacci count')
        records=[]
        for key,(word,A) in sorted(layer.items()):
            st,M=exact_stats(A)
            require(st['det'] in (0,1,-1),'new nonsingular non-unimodular unique-cycle witness')
            rec={'bits':key,'word':word,'type':'A' if sum(sum(row)==n-2 for row in A)==1 else 'B',**st}
            records.append(rec)
        uni=sum(abs(t['det'])==1 for t in records)
        output['layers'][str(n)]={'classes':len(layer),'unimodular':uni,'records':records}
        args.output.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
        print(f'n={n} classes={len(layer)} unimodular={uni}',flush=True)
        if n < args.max_n:
            nxt={}
            for word,A in layer.values():
                for letter,key,B in children(A):
                    require(key not in nxt,'grammar must give distinct cycle-normalized classes')
                    nxt[key]=(word+letter,B)
            layer=nxt
            f0,f1=f1,3*f1-f0


if __name__=='__main__': main()
