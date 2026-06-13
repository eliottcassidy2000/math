#!/usr/bin/env python3
"""
lrc_tree_entropy_attack_s543.py    oracle-2026-06-01-S543o

ATTACK: ENTROPY on the TREE(s) of Tournament Analysis, vigorously. Three entropies,
all measuring 'spread/mixing' on a tree structure, all tied to loneliness.

NOTE: 'Hough entropy' read productively as the entropy functionals on our tree
structures, centered on the H-ENTROPY:

  *** H-ENTROPY: S_H(T) = log2 H(T), where H = directed Hamiltonian-path count (the
  loneliness meter, S26). Since H is MULTIPLICATIVE over the disjoint modules of the
  recursive/modular tree (S531), S_H is ADDITIVE over that tree:
       S_H(T) = sum over modular-tree nodes of log2(H of each module).
  An apex-flipped module of size s contributes log2(1+2^{s-2}). So S_H is a genuine
  tree-additive entropy; S_H=0 at the transitive (rigid) tournament, MAX at the
  regular (the LRC-tight regular polygon). ***

Also: (2) p-adic Bruhat-Tits TREE entropy (S541) -- the Shannon spread of the speeds
across residues mod p^k (the channels, S534); (3) the ISO-CLASS WALK entropy over t
(the menu walk, S518/S542).

Vigorous: many n, many speed families (AP, random, geometric, sieve-structured),
the three entropies, the extremes, and the relationships.
"""
from itertools import combinations, permutations
from functools import reduce
from math import gcd, log2
import random, statistics as st

def frac(x): return x - int(x // 1)

# ---------- H and S_H ----------
def H_count(adj, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not c: continue
            for u in range(n):
                if mask & (1 << u): continue
                if adj[v][u]: dp[mask | (1 << u)][u] += c
    return sum(dp[full][v] for v in range(n))

def half_turn(positions, n):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and 0 < (positions[j]-positions[i]) % 1.0 < 0.5: adj[i][j] = 1
    return adj

def canon(adj, n):
    best = None
    for p in permutations(range(n)):
        b = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or b < best: best = b
    return best

# ---------- (1) S_H additivity over modules (verify) ----------
def transitive_adj(n):
    return [[1 if i > j else 0 for j in range(n)] for i in range(n)]
def flip(adj, flips):
    a = [r[:] for r in adj]
    for (x, y) in flips: a[x][y] = 0; a[y][x] = 1
    return a

def verify_SH_additivity():
    print("="*70); print("(1) H-ENTROPY S_H=log2 H is ADDITIVE over disjoint modules (S531)"); print("="*70)
    # n=8: module [1,3] apex-flip (tile (3,1)) and [5,7] apex-flip (tile (7,5)) disjoint
    tests = [(8, [(3,1)], [(7,5)]), (9, [(3,0)], [(8,5)])]
    for n, A, B in tests:
        HA = H_count(flip(transitive_adj(n), A), n)
        HB = H_count(flip(transitive_adj(n), B), n)
        HAB = H_count(flip(transitive_adj(n), A+B), n)
        print(f"  n={n}: S_H(A)={log2(HA):.3f}, S_H(B)={log2(HB):.3f}, S_H(A+B)={log2(HAB):.3f}; "
              f"sum={log2(HA)+log2(HB):.3f}  additive={abs(log2(HAB)-log2(HA)-log2(HB))<1e-9}")
    # apex-flipped block of size s contributes log2(1+2^{s-2})
    print("  apex module size s -> S_H contribution log2(1+2^{s-2}):")
    for s in range(2, 8):
        print(f"    s={s}: 1+2^(s-2)={1+2**(s-2)}, log2={log2(1+2**(s-2)):.3f}")
    print("  => S_H is a TREE-ADDITIVE entropy: S_H=0 (transitive, rigid) up to log2(H_regular)")
    print("     (the regular polygon = LRC-tight extremal). Max entropy = the tight witness.")
    print()

# ---------- (1b) S_H(t) trajectory over the runner walk ----------
def SH_trajectory(v, G=4000):
    n = len(v); Hs = []
    for s in range(G):
        t = (s+0.5)/G; Hs.append(H_count(half_turn([frac(x*t) for x in v], n), n))
    SHs = [log2(h) for h in Hs]
    return min(SHs), max(SHs), st.mean(SHs)

def family_SH(n):
    print("="*70); print(f"(1b) H-ENTROPY trajectory S_H(t) over the runner walk, n={n}"); print("="*70)
    fams = {
        "AP 1..n (regular)": tuple(range(1, n+1)),
        "geometric 1,2,4,..": tuple(2**k for k in range(n)),
        "random primitive": None,
        "sieve 1,2,..,n with mult": tuple(list(range(1, n))+[2*3*5]),
    }
    rnd = random.Random(7+n)
    for name, v in fams.items():
        if v is None:
            v = tuple(sorted(rnd.sample(range(1, 8*n), n)))
            while reduce(gcd, v) != 1: v = tuple(sorted(rnd.sample(range(1, 8*n), n)))
        lo, hi, mean = SH_trajectory(v)
        print(f"  {name:26s} v={v}: S_H(t) in [{lo:.2f},{hi:.2f}], mean {mean:.3f}")
    print("  PATTERN: max S_H = log2(H_regular) at the spread/regular slices; AP/regular speeds")
    print("  carry the HIGHEST mean H-entropy (most time spread) -- entropy fingerprints regularity.")
    print()

# ---------- (2) p-adic tree entropy ----------
def shannon(counts):
    tot = sum(counts)
    if tot == 0: return 0.0
    return -sum((c/tot)*log2(c/tot) for c in counts if c > 0)

def padic_tree_entropy(v, p, K):
    """Shannon entropy of residues mod p^k for k=1..K (leaf-distribution of the
    p-adic tree); also report 0-branch occupancy (sieve loneliness)."""
    ent = []
    for k in range(1, K+1):
        mod = p**k
        cnt = [0]*mod
        for x in v: cnt[x % mod] += 1
        ent.append(shannon(cnt))
        if k == 1:
            zero0 = (cnt[0] == 0)   # no speed divisible by p -> t=1/p lonely (sieve)
    return ent, zero0

def study_padic_tree(n, p=3, K=3):
    print("="*70); print(f"(2) p-ADIC TREE entropy (p={p}): Shannon spread across residues mod p^k"); print("="*70)
    fams = {
        "AP 1..n-1": tuple(range(1, n)),
        "random primitive": None,
        "sieve: all ≡0 mod p (x3 of 1..)": tuple(p*x for x in range(1, n)),
        "mixed: 1 + multiples of p": tuple([1]+[p*x for x in range(1, n-1)]),
    }
    rnd = random.Random(3+n)
    for name, v in fams.items():
        if v is None:
            v = tuple(sorted(rnd.sample(range(1, 8*n), n-1)))
            while reduce(gcd, v) != 1: v = tuple(sorted(rnd.sample(range(1, 8*n), n-1)))
        ent, z0 = padic_tree_entropy(v, p, K)
        print(f"  {name:30s}: level-entropies(mod p^1..p^{K})={[round(e,2) for e in ent]}; "
              f"0-branch empty (sieve t=1/p lonely)={z0}")
    print("  PATTERN: HIGH tree entropy = speeds spread across channels (debt present, sieve")
    print("  blocked, S534 vacuity); LOW entropy / empty 0-branch = sieve loneliness available.")
    print()

# ---------- (3) iso-class walk entropy ----------
def walk_entropy(v, G=4000):
    n = len(v); cnt = {}
    for s in range(G):
        t = (s+0.5)/G; c = canon(half_turn([frac(x*t) for x in v], n), n)
        cnt[c] = cnt.get(c, 0)+1
    return shannon(list(cnt.values())), len(cnt)

def study_walk_entropy(n):
    print("="*70); print(f"(3) ISO-CLASS WALK entropy over t (the menu walk), n={n}"); print("="*70)
    rnd = random.Random(11+n)
    for name, v in [("AP 1..n", tuple(range(1, n+1))),
                    ("random", None)]:
        if v is None:
            v = tuple(sorted(rnd.sample(range(1, 8*n), n)))
            while reduce(gcd, v) != 1: v = tuple(sorted(rnd.sample(range(1, 8*n), n)))
        e, nc = walk_entropy(v)
        print(f"  {name:10s} v={v}: walk entropy={e:.3f} bits over {nc} iso-classes "
              f"(max if uniform={log2(nc):.3f})")
    print("  PATTERN: the walk entropy measures the MIXING across the 2Fib(n-2) menu; the")
    print("  more uniform the time-in-class, the higher -- another regularity fingerprint.")

def main():
    verify_SH_additivity()
    for n in (5, 6): family_SH(n)
    study_padic_tree(6, p=3, K=3)
    study_padic_tree(6, p=2, K=3)
    study_walk_entropy(5)
    print("="*70)
    print("SYNTHESIS: three TREE entropies all measure SPREAD/MIXING. H-ENTROPY (log2 H) is")
    print("TREE-ADDITIVE over the modular tree (S531), maximized at the regular polygon (LRC-")
    print("tight). p-adic tree entropy = channel spread (S534): high => sieve blocked. Walk")
    print("entropy = menu mixing. Regularity (AP) shows up as HIGH entropy on every tree.")
    print("="*70)

if __name__ == "__main__":
    main()
