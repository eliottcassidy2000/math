#!/usr/bin/env python3
"""
lrc_mappings_restriction_spectrum_s535.py    oracle-2026-06-01-S535o

A MULTITUDE of mappings LRC -> combinatorial structure space, ranked by how
MEANINGFULLY RESTRICTED the realizable iso-class set is.

META-PRINCIPLE (S519/S529): a plain tournament keeps only the ORDER (sign of
(v_i-v_j)t mod 1 in (0,1/2)); the realizable set is then the whole circular menu
(2*Fib(n-2)) -- not restrictive, the conjecture hides in the WALK. Loneliness is
METRIC (the 1/n threshold). So the realizable iso-class set shrinks as the mapping
retains MORE metric. We build a SPECTRUM of mappings and measure
   R = |realizable iso-classes| / |all iso-classes of the target structure|
for small n, and whether LRC = a membership statement.

RUNGS (increasing metric content):
  M0  order tournament (half-turn): target = tournaments on n vertices.
  M1  NEAR-GRAPH (threshold 1/n): i~j iff circular dist < 1/n. target = graphs.
      LRC <=> the OBSERVER-ISOLATED near-graph is exhibitable. <-- the key rung.
  M2  3-LEVEL metric tournament: each edge in {near(<1/n), mid, far} -> target =
      3-edge-colored complete graphs (huge); realizable tiny.
  M3  resonance/QR tournament (n prime): i->j iff (v_i-v_j) is a QR mod n. static.
  (and a multitude of further abstract mappings stated as hypotheses in the writeup)

We compute the realizable sets and R for M1, M2, M3 at small n, and verify the LRC
membership statement for M1.
"""
from itertools import combinations, permutations, product
from functools import reduce
from math import gcd
import random

# ---------- canonical forms ----------
def canon_graph(edges, n):
    """edges: frozenset of frozenset{i,j}. canonical = min adjacency bitstring over perms."""
    adj = [[0]*n for _ in range(n)]
    for e in edges:
        i, j = tuple(e); adj[i][j] = adj[j][i] = 1
    best = None
    for p in permutations(range(n)):
        bits = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(i+1, n))
        if best is None or bits < best: best = bits
    return best

def canon_tournament(adj, n):
    best = None
    for p in permutations(range(n)):
        bits = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or bits < best: best = bits
    return best

def canon_colored(coledges, n):
    """coledges: dict frozenset{i,j}->color. canonical over perms."""
    best = None
    for p in permutations(range(n)):
        seq = tuple(coledges[frozenset((p[i], p[j]))] for i in range(n) for j in range(i+1, n))
        if best is None or seq < best: best = seq
    return best

# all-graph counts (OEIS A000088): n=3:4,4:11,5:34,6:156,7:1044
ALLGRAPHS = {3:4, 4:11, 5:34, 6:156, 7:1044}
A000568 = {3:2,4:4,5:12,6:56,7:456}

def frac(x): return x - int(x//1)
def cdist(a, b):
    d = abs(a-b) % 1.0
    return min(d, 1-d)

# ---------- M1: near-graph ----------
def near_graphs_of_set(speeds, n, samples=3000):
    """union over t of near-graph edge-sets (circular dist < 1/n) on observer+runners."""
    thr = 1.0/n
    seen = set(); lonely_seen = False
    pts_runners = speeds
    for s in range(samples):
        t = (s+0.5)/samples
        pos = [0.0] + [frac(v*t) for v in pts_runners]
        edges = frozenset(frozenset((i, j)) for i in range(n) for j in range(i+1, n)
                          if cdist(pos[i], pos[j]) < thr)
        seen.add(edges)
        # observer isolated?
        if not any(0 in e for e in edges):
            lonely_seen = True
    return seen, lonely_seen

def study_M1(n, n_sets=120, seedbase=1):
    rnd = random.Random(seedbase*7+n)
    global_classes = set()
    lonely_ok = 0; tot = 0
    per_set_sizes = []
    for _ in range(4000):
        v = tuple(sorted(rnd.sample(range(1, 6*n), n-1)))
        if reduce(gcd, v) != 1: continue
        tot += 1
        if tot > n_sets: break
        edgesets, lonely = near_graphs_of_set(v, n)
        per_set_sizes.append(len(edgesets))
        if lonely: lonely_ok += 1
        for es in edgesets:
            global_classes.add(canon_graph(es, n))
    return global_classes, lonely_ok, tot, per_set_sizes

# ---------- M2: 3-level metric tournament ----------
def metric3_of_set(speeds, n, samples=2000):
    """edge color: 0 near (<1/n), 1 mid (1/n..2/n), 2 far (>=2/n) by circular dist.
    union over t of colored-complete-graph iso-classes."""
    seen = set()
    for s in range(samples):
        t = (s+0.5)/samples
        pos = [0.0] + [frac(v*t) for v in speeds]
        col = {}
        for i in range(n):
            for j in range(i+1, n):
                d = cdist(pos[i], pos[j])
                col[frozenset((i, j))] = 0 if d < 1.0/n else (1 if d < 2.0/n else 2)
        seen.add(canon_colored(col, n))
    return seen

def study_M2(n, n_sets=60):
    rnd = random.Random(99+n); g = set(); tot = 0
    for _ in range(4000):
        v = tuple(sorted(rnd.sample(range(1, 6*n), n-1)))
        if reduce(gcd, v) != 1: continue
        tot += 1
        if tot > n_sets: break
        g |= metric3_of_set(v, n)
    # total 3-colorings of C(n,2) edges (labeled, upper bound on iso-classes)
    import math
    total_labeled = 3**(n*(n-1)//2)
    return g, total_labeled, tot

# ---------- M3: resonance / QR tournament (n prime) ----------
def is_qr(a, p):
    a %= p
    if a == 0: return None
    return pow(a, (p-1)//2, p) == 1

def qr_tournament_of_set(speeds, n):
    """on the n-1 runners: i->j iff (v_i - v_j) is QR mod n. (n prime)"""
    m = len(speeds)
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            q = is_qr(speeds[i]-speeds[j], n)
            adj[i][j] = 1 if q else 0
    return canon_tournament(adj, m)

def study_M3(n, n_sets=400):
    """n prime; tournament on m=n-1 runners."""
    rnd = random.Random(5+n); g = set(); tot = 0
    for _ in range(8000):
        v = tuple(sorted(rnd.sample(range(1, 8*n), n-1)))
        if reduce(gcd, v) != 1: continue
        # require all pairwise differences nonzero mod n (else QR undefined)
        if any((v[i]-v[j]) % n == 0 for i in range(len(v)) for j in range(i+1, len(v))): continue
        tot += 1
        if tot > n_sets: break
        g.add(qr_tournament_of_set(v, n))
    return g, A000568.get(n-1), tot

# ----------------------------------------------------------------------
def main():
    print("="*74)
    print("META: realizable iso-class set shrinks as the mapping retains more METRIC")
    print("="*74)

    print("\n--- M1: NEAR-GRAPH (threshold 1/n) --- LRC <=> observer-isolated is exhibitable")
    for n in (4, 5, 6):
        g, lonely_ok, tot, sizes = study_M1(n, n_sets=(120 if n < 6 else 60))
        avg = sum(sizes)/len(sizes)
        print(f"  n={n}: realizable near-graph iso-classes = {len(g)} of {ALLGRAPHS[n]} all graphs "
              f"(R={len(g)/ALLGRAPHS[n]:.2f}); per-set avg distinct near-graphs={avg:.1f}; "
              f"LRC(observer-isolated exhibited)={lonely_ok}/{tot}")
    print("  => near-graphs are a RESTRICTED family (circular indifference graphs); LRC = the")
    print("     marked observer-isolated near-graph lies in EVERY speed set's realizable union.")

    print("\n--- M2: 3-LEVEL METRIC tournament (near/mid/far) --- massive restriction")
    for n in (4, 5):
        g, total_labeled, tot = study_M2(n, n_sets=(80 if n == 4 else 50))
        print(f"  n={n}: realizable 3-level iso-classes = {len(g)}  (vs {total_labeled} labeled "
              f"3-colorings of {n*(n-1)//2} edges); astronomically restricted")

    print("\n--- M3: RESONANCE/QR tournament (n prime) --- arithmetic restriction")
    for n in (5, 7):
        g, total_iso, tot = study_M3(n)
        print(f"  n={n} (prime): realizable QR-tournament iso-classes on {n-1} runners = {len(g)} "
              f"of A000568({n-1})={total_iso} (R={len(g)/total_iso:.3f}); {tot} sets sampled")
    print("  => QR-difference tournaments hit only a thin arithmetic slice (Paley-type),")
    print("     restricted by reciprocity/difference-set constraints.")

    print("\n" + "="*74)
    print("RANKING (small-n restriction R = realizable/all): M0 order ~ 2Fib/A000568;")
    print("M1 near-graph (metric threshold) strong; M2 3-level near-total collapse;")
    print("M3 QR arithmetic slice. The conjecture lives at the METRIC rungs (M1+).")
    print("="*74)

if __name__ == "__main__":
    main()
