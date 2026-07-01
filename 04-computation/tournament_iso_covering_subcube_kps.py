#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE MINIMAL-FLIP ISO-COVERING SUBCUBE of the tournament hypercube.

kind-pasteur-2026-07-01-S10. The owner's question: at n=4 all 4 iso classes are reachable by flipping
3 arcs (the naive Hamiltonian-path tiling), but ALSO by flipping only 2 arcs IF the 4 FIXED arcs are in
a certain configuration. Study how this minimal shape changes with n; reframe; pattern-seek.

FRAME. Tournaments on n vertices = the hypercube Q_m, m=C(n,2) (bit p = orientation of arc-pair p).
S_n permutes the coordinates; iso classes = S_n-orbits (count = A000568(n) = 1,1,2,4,12,56,456,...).
An ISO-COVERING SUBCUBE of dimension k = a choice of k FREE arcs D + a fixed orientation of the other
m-k FIXED arcs, whose 2^k tournaments hit EVERY iso class. k_min(n) = the least such k.
 * INFO bound: k >= ceil(log2 A000568(n)).
 * GAUGE reframe: #fixed f = m - k; fixing arcs GAUGE-FIXES the S_n relabeling, so f_max = m - k_min
   ~ log2|S_n| = log2(n!) (the SYMMETRY ENTROPY). The naive Ham-path fixes only n-1 arcs (a weak gauge);
   the optimal gauge fixes ~log2(n!) arcs -- MORE, because most tournaments are rigid (Burnside).
Computed EXACTLY for n=3,4,5; the configuration of the fixed arcs characterized.
"""
import sys, itertools
from math import log2, factorial
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
A000568={3:2,4:4,5:12,6:56,7:456,8:6880}

def setup(n):
    prs=list(itertools.combinations(range(n),2)); idx={p:k for k,p in enumerate(prs)}
    perms=list(itertools.permutations(range(n))); return prs,idx,perms
def beats(t,idx,i,j):
    return (t>>idx[(i,j)])&1 if i<j else 1-((t>>idx[(j,i)])&1)
def canon(t,prs,idx,perms):
    best=None
    for pm in perms:
        t2=0
        for k,(a,b) in enumerate(prs):
            if beats(t,idx,pm[a],pm[b]): t2|=(1<<k)
        if best is None or t2<best: best=t2
    return best
def precompute_canon(n):
    prs,idx,perms=setup(n); m=len(prs); arr=[0]*(1<<m)
    for t in range(1<<m): arr[t]=canon(t,prs,idx,perms)
    return arr,prs,idx,m

def covers(Dpos,fixed_t,arr,nclass):
    seen=set(); k=len(Dpos)
    for s in range(1<<k):
        t=fixed_t
        ss=s; b=0
        while ss:
            if ss&1: t^=(1<<Dpos[b])
            ss>>=1; b+=1
        seen.add(arr[t])
        if len(seen)==nclass: return True
    return False

def kmin_and_configs(n, arr, m, nclass, want_configs=8):
    lb=1
    while (1<<lb)<nclass: lb+=1
    for k in range(lb, m+1):
        sols=[]
        for D in itertools.combinations(range(m), k):
            Fpos=[p for p in range(m) if p not in D]
            fbits=len(Fpos)
            for fo in range(1<<fbits):
                fixed_t=0; x=fo; b=0
                while x:
                    if x&1: fixed_t|=(1<<Fpos[b])
                    x>>=1; b+=1
                if covers(list(D),fixed_t,arr,nclass):
                    sols.append((D,fixed_t));
                    break  # one witness orientation per D is enough to record D as viable
            if len(sols)>=want_configs and k>lb: break
        if sols: return k,lb,sols
    return None,lb,[]

def describe(D, fixed_t, prs, n):
    freearcs=[prs[p] for p in D]
    # free-arc graph structure
    verts=set()
    for (a,b) in freearcs: verts|={a,b}
    deg={}
    for (a,b) in freearcs: deg[a]=deg.get(a,0)+1; deg[b]=deg.get(b,0)+1
    shape="matching" if all(v<=1 for v in deg.values()) else ("star" if max(deg.values())==len(freearcs) else "path/other")
    return freearcs, shape, sorted(deg.items())

print("="*94); print(" k_min(n): the minimal FREE-arc count to hit all iso classes; the FIXED-arc gauge"); print("="*94)
print(f"  {'n':>2} {'m=C(n,2)':>8} {'#classes':>8} {'log2#cls':>8} {'k_min':>6} {'f=m-kmin':>9} {'log2(n!)':>8}")
results={}
for n in [3,4,5]:
    arr,prs,idx,m=precompute_canon(n); nclass=len(set(arr))
    k,lb,sols=kmin_and_configs(n,arr,m,nclass,want_configs=200)
    results[n]=(m,nclass,k,lb,sols,prs)
    print(f"  {n:>2} {m:>8} {nclass:>8} {log2(nclass):>8.2f} {k:>6} {m-k:>9} {log2(factorial(n)):>8.2f}")
print("  (info bound ceil(log2#cls) vs k_min; f=fixed arcs vs log2(n!)=symmetry entropy)")

for n in [3,4,5]:
    m,nclass,k,lb,sols,prs=results[n]
    print("\n"+"="*94); print(f" n={n}: k_min={k} (info-bound lb={lb}; TIGHT={k==lb}); {len(sols)} viable FREE-arc sets"); print("="*94)
    # count how many D-subsets are viable, and characterize free-arc shapes
    shapes={}
    for (D,ft) in sols:
        fa,shape,deg=describe(D,ft,prs,n)
        shapes[shape]=shapes.get(shape,0)+1
    print(f"  viable free-arc-set SHAPES: {shapes}")
    # show 3 sample configs
    for (D,ft) in sols[:3]:
        fa,shape,deg=describe(D,ft,prs,n)
        fixedarcs=[prs[p] for p in range(m) if p not in D]
        # orientation of fixed arcs
        fdir=[]
        for p in range(m):
            if p not in D:
                (a,b)=prs[p]; fdir.append(f"{a}->{b}" if (ft>>p)&1 else f"{b}->{a}")
        print(f"   FREE {fa} ({shape}); FIXED arcs {fixedarcs} oriented {fdir}")

print("\n"+"="*94); print(" PATTERN / REFRAME"); print("="*94)
print("  k_min is the DIMENSION of the smallest axis-aligned subcube meeting every S_n-orbit.")
print("  f_max=m-k_min FIXED arcs GAUGE-FIX the relabeling symmetry; the bound f_max<=log2(n!) is the")
print("  symmetry entropy. Naive Ham-path fixes only n-1 (weak); optimal fixes ~log2(n!) (strong gauge).")
print("DONE.")
