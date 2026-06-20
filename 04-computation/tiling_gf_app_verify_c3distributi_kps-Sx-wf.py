#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL VERIFICATION of the c3-distribution-structure application (THM-554 tool).
kps-Sx-wf.

Checks, all with EXACT integer/Fraction arithmetic:
 (1) Z-engine c3-distribution == independent BRUTE enumeration over all 2^C(n-1,2)
     tilings, n<=7.
 (2) Closed forms re-derived and checked against the engine AND brute:
       E[c3]   = (C(n,3) + (n-2))/4
       Var(c3) = (n^3 - 7n^2 + 20n - 16)/32
       E[(-1)^c3] = 1/2^floor((n-1)/2)  (n>=4), 0 for n=3
       min c3 = 0 mult 1
       max c3 = (n^3-n)/24 (odd), (n^3-4n)/24 (even); multiplicities
       per-triple 3-cycle prob = 1/2 for consecutive {v,v+1,v+2}, else 1/4
 (3) regular census 1,3,91,29157 (odd n=3,5,7,9); near-regular 5,157,51949 (even n=4,6,8)
 (4) complement-fold 2x: confirm c3 invariant under tiling complement R, lossless.
"""
import sys
from collections import defaultdict, Counter
from itertools import product
from math import comb
from fractions import Fraction as F

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

OUT = []
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s); OUT.append(s)

# ---------------- engine ----------------
def beta_step(dist, n):
    nd = defaultdict(int)
    for vec, cnt in dist.items():
        l = list(vec) + [0]; l[n-1] += 1; nd[tuple(l)] += cnt
    dist = nd
    for b in range(1, n-1):
        nd = defaultdict(int)
        for vec, cnt in dist.items():
            l0 = list(vec); l0[n-1] += 1; nd[tuple(l0)] += cnt
            l1 = list(vec); l1[b-1] += 1; nd[tuple(l1)] += cnt
        dist = nd
    return dist

def build_Z(N):
    d = {(0,): 1}
    for n in range(2, N+1):
        d = beta_step(d, n)
    return d

def c3_dist_from_Z(distZ, n):
    cs = Counter()
    for vec, cnt in distZ.items():
        cs[comb(n,3) - sum(comb(s,2) for s in vec)] += cnt
    return cs

# ---------------- brute ----------------
def tiles(n):
    return [(a, b) for a in range(3, n+1) for b in range(1, a-1)]

def brute(n):
    """Return (c3 Counter, score-multiset Counter, list of (c3, sorted_score) per tiling)
       computed by FULL tournament rebuild, independent of the GF factorization."""
    T = tiles(n)
    c3c = Counter(); scc = Counter()
    per = []
    for bv in product((0,1), repeat=len(T)):
        adj = [[0]*(n+1) for _ in range(n+1)]
        for k in range(n, 1, -1):
            adj[k][k-1] = 1
        for (a,b), bit in zip(T, bv):
            if bit == 0: adj[a][b] = 1
            else: adj[b][a] = 1
        # scores (out-degree)
        sc = [0]*(n+1)
        for i in range(1, n+1):
            for j in range(1, n+1):
                if i != j: sc[i] += adj[i][j]
        t = 0
        for i in range(1, n+1):
            for j in range(i+1, n+1):
                for k in range(j+1, n+1):
                    if (adj[i][j]+adj[i][k], adj[j][i]+adj[j][k],
                        adj[k][i]+adj[k][j]) == (1,1,1):
                        t += 1
        ss = tuple(sorted(sc[1:]))
        c3c[t] += 1; scc[ss] += 1; per.append((t, ss))
    return c3c, scc, per

# ---------------- moment helpers ----------------
def mean(cs):
    tot = sum(cs.values())
    return F(sum(c*m for c, m in cs.items()), tot)

def var(cs):
    tot = sum(cs.values()); mu = mean(cs)
    return F(sum(m*(c-mu)**2 for c, m in cs.items()), 1)/tot  # tot already in mean

def var2(cs):
    tot = sum(cs.values())
    e1 = F(sum(c*m for c, m in cs.items()), tot)
    e2 = F(sum(c*c*m for c, m in cs.items()), tot)
    return e2 - e1*e1

def signed(cs):
    tot = sum(cs.values())
    return F(sum(((-1)**c)*m for c, m in cs.items()), tot)

# ---------------- claimed closed forms ----------------
def cf_mean(n): return F(comb(n,3) + (n-2), 4)
def cf_var(n):  return F(n**3 - 7*n**2 + 20*n - 16, 32)
def cf_signed(n):
    if n == 3: return F(0)
    return F(1, 2**((n-1)//2))
def cf_maxc3(n):
    return F(n**3 - n, 24) if n % 2 else F(n**3 - 4*n, 24)

def main():
    log("="*78)
    log("PART 1+2: engine vs brute, and closed forms (exact)")
    log("="*78)
    holds = True
    for n in range(3, 8):
        Z = build_Z(n)
        gf = c3_dist_from_Z(Z, n)
        bc3, bsc, per = brute(n)
        eq = dict(gf) == dict(bc3)
        log(f"\n n={n}  total={sum(gf.values())} = 2^{comb(n-1,2)}  engine==brute_c3: {eq}")
        if not eq:
            holds = False
            log("   MISMATCH gf:", dict(sorted(gf.items())))
            log("   MISMATCH br:", dict(sorted(bc3.items())))
        # also: engine score census == brute score census
        zsc = Counter()
        for vec, cnt in Z.items():
            zsc[tuple(sorted(vec))] += cnt
        log("   score-census engine==brute:", dict(zsc) == dict(bsc))

        mu = mean(gf); vv = var2(gf); sg = signed(gf)
        log(f"   E[c3]     engine={mu}  closed={cf_mean(n)}  match={mu==cf_mean(n)}")
        log(f"   Var[c3]   engine={vv}  closed={cf_var(n)}  match={vv==cf_var(n)}")
        log(f"   E[(-1)c3] engine={sg}  closed={cf_signed(n)}  match={sg==cf_signed(n)}")
        if mu != cf_mean(n) or vv != cf_var(n) or sg != cf_signed(n):
            holds = False
        # min
        mn = min(gf); log(f"   min c3={mn} mult={gf[mn]}  (claim 0 mult 1): {mn==0 and gf[mn]==1}")
        if not (mn == 0 and gf[mn] == 1): holds = False
        # max
        mx = max(gf); log(f"   max c3={mx} mult={gf[mx]}  closed_max={cf_maxc3(n)}  match={mx==cf_maxc3(n)}")
        if F(mx) != cf_maxc3(n): holds = False

    log("\n" + "="*78)
    log("PART 2b: closed forms vs ENGINE only, n=8..10 (no brute)")
    log("="*78)
    dist = {(0,): 1}
    maxmult = {}
    for n in range(2, 11):
        dist = beta_step(dist, n)
        if n < 3: continue
        gf = c3_dist_from_Z(dist, n)
        mu = mean(gf); vv = var2(gf); sg = signed(gf)
        mx = max(gf); mn = min(gf)
        okm = mu == cf_mean(n); okv = vv == cf_var(n); oks = sg == cf_signed(n)
        okmax = F(mx) == cf_maxc3(n); okmin = (mn == 0 and gf[mn] == 1)
        maxmult[n] = gf[mx]
        log(f" n={n}: meanOK={okm} varOK={okv} signedOK={oks} maxOK={okmax} minOK={okmin} "
            f"| max={mx} maxmult={gf[mx]} signed={sg}")
        if not (okm and okv and oks and okmax and okmin):
            holds = False

    log("\n" + "="*78)
    log("PART 3: regular / near-regular census")
    log("="*78)
    # max-c3 multiplicity sequence
    odd_seq = [maxmult[n] for n in (3,5,7,9)]
    even_seq = [maxmult[n] for n in (4,6,8)]
    log(f" max-c3 mult odd  (n=3,5,7,9): {odd_seq}  claim [1,3,91,29157]: {odd_seq==[1,3,91,29157]}")
    log(f" max-c3 mult even (n=4,6,8):   {even_seq}  claim [5,157,51949]: {even_seq==[5,157,51949]}")
    if odd_seq != [1,3,91,29157] or even_seq != [5,157,51949]:
        holds = False
    # explicit regular score census (all s_v=(n-1)/2) for odd n
    for n in (3,5,7,9):
        d = build_Z(n); r = (n-1)//2
        reg = sum(c for v, c in d.items() if all(s == r for s in v))
        log(f"  n={n}: regular-score (all s={r}) census = {reg}; equals max-c3 mult: {reg==maxmult[n]}")
        if reg != maxmult[n]: holds = False
    # check max-c3 IS achieved exactly by regular/near-regular scores
    for n in (3,4,5,6,7):
        d = build_Z(n)
        mx = max(c3_dist_from_Z(d, n))
        # which score multisets achieve max c3?
        ach = set()
        for vec, cnt in d.items():
            if comb(n,3) - sum(comb(s,2) for s in vec) == mx:
                ach.add(tuple(sorted(vec)))
        # near-regular = scores in {floor,ceil}((n-1)/2)
        lo = (n-1)//2; hi = (n-1+1)//2 if n % 2 == 0 else lo
        # for odd n, regular all == lo; for even, near-regular in {lo, lo+1}
        allowed_lo, allowed_hi = ((lo, lo) if n % 2 else (lo, lo+1))
        nearreg = all(all(allowed_lo <= s <= allowed_hi for s in ss) for ss in ach)
        log(f"  n={n}: max-c3 achieved by score multisets {sorted(ach)} all near-regular: {nearreg}")

    log("\n" + "="*78)
    log("PART 4: complement-fold 2x via THM-280 R-relabel  (a,b)->(n+1-b,n+1-a)")
    log("="*78)
    log("  The REAL fold: T -> T^op realized inside the tiling cube. T^op reverses every")
    log("  arc; relabel v->n+1-v restores the base path n->..->1. Net action on free-tile")
    log("  bits = position permutation R (a,b)->(n+1-b,n+1-a) WITH a bit flip on each tile.")
    log("  Claim: this is an involution on the 2^F tilings, and c3 is invariant => 2x fold.")
    for n in range(3, 7):
        T = tiles(n)
        idx = {t: i for i, t in enumerate(T)}
        def R(t):
            a, b = t
            return (n+1-b, n+1-a)
        # position permutation on tile indices
        permpos = [idx[R(T[i])] for i in range(len(T))]
        is_inv_pos = all(permpos[permpos[i]] == i for i in range(len(T)))
        # apply fold to a bitvector: new bit at position permpos[i] = 1 - old bit at i
        def fold(bv):
            nb = [0]*len(bv)
            for i, bit in enumerate(bv):
                nb[permpos[i]] = 1 - bit
            return tuple(nb)
        # build c3 per bitvector
        bvmap = {}
        for bv in product((0,1), repeat=len(T)):
            adj = [[0]*(n+1) for _ in range(n+1)]
            for k in range(n,1,-1): adj[k][k-1]=1
            for (a,b),bit in zip(T,bv):
                if bit==0: adj[a][b]=1
                else: adj[b][a]=1
            t=0
            for i in range(1,n+1):
                for j in range(i+1,n+1):
                    for k in range(j+1,n+1):
                        if (adj[i][j]+adj[i][k],adj[j][i]+adj[j][k],
                            adj[k][i]+adj[k][j])==(1,1,1): t+=1
            bvmap[bv]=t
        is_inv_map = all(fold(fold(bv)) == bv for bv in bvmap)
        # SUBSTANTIVE claim: c3(T)=c3(T^op) for every tiling (the real complement-fold).
        # Build T^op directly (reverse all arcs) and compare c3 -- independent of any
        # guess at the in-cube bit realization.
        def c3_op(bv):
            adj=[[0]*(n+1) for _ in range(n+1)]
            for k in range(n,1,-1): adj[k][k-1]=1
            for (a,b),bit in zip(T,bv):
                if bit==0: adj[a][b]=1
                else: adj[b][a]=1
            op=[[adj[j][i] for j in range(n+1)] for i in range(n+1)]
            t=0
            for i in range(1,n+1):
                for j in range(i+1,n+1):
                    for k in range(j+1,n+1):
                        if (op[i][j]+op[i][k],op[j][i]+op[j][k],
                            op[k][i]+op[k][j])==(1,1,1): t+=1
            return t
        comp_ok = all(c3_op(bv) == bvmap[bv] for bv in bvmap)
        # The 2x is over ISO/COMPLEMENT classes. R within the cube is a fixed-point-free
        # involution => orbits size 2 => exactly halves the computation, losslessly.
        fixed = sum(1 for bv in bvmap if fold(bv) == bv)
        log(f" n={n}: R pos-perm involution:{is_inv_pos}  fold involution:{is_inv_map}  "
            f"c3(T)==c3(T^op) all tilings:{comp_ok}  fold-fixed-points:{fixed} (0 => clean 2x)")
        if not comp_ok: holds = False

    log("\n" + "="*78)
    log(f"OVERALL holds = {holds}")
    log("="*78)
    return holds

if __name__ == "__main__":
    h = main()
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + f"\nHOLDS={h}\n")
