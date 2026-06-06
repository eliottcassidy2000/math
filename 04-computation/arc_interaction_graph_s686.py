#!/usr/bin/env python3
"""arc_interaction_graph_s686.py — the EXACT combinatorial criterion for the
arc-flip Hessian M(g,e), completing HYP-2268's abstract Walsh formula and
correcting s685's single-channel guess.

SETUP. H(T) is always ODD (Redei). For an arc e, delta_T(e)=H(T^e)-H(T) is the
discrete gradient (always EVEN). The Hessian (interaction) is
    M(g,e) = H(T) - H(T^g) - H(T^e) + H(T^{ge}),
the change in delta_T(e) caused by flipping g. HYP-2268: M = 4*Sum_{S>={g,e}} c_S
chi_S over OCF terms S, where each S is a collection of VERTEX-DISJOINT directed
ODD cycles and chi_S the product of its arc variables.

THE SUPPORT IS ORIENTATION-FREE (the key correction over s685). A term S
contributes iff its arc set contains BOTH g and e; the support {S : c_S!=0,
S>={g,e}} is a property of the underlying K_n and the POSITIONS of g,e -- NOT of
T's orientation (the c_S are the global Walsh coefficients; only the signs
chi_S(T) depend on T). So:
    M(g,e) != 0   =>   g,e are jointly coverable by vertex-disjoint ODD cycles
                        in K_n  (= support nonempty).   [NECESSARY condition]
The minimal covers are (A) a single odd cycle through both, or (B) two disjoint
odd cycles each carrying one arc (>=6 vertices). The converse FAILS: a nonempty
support can still cancel (signed sum = 0) for particular T.

CONSEQUENCES (verified below):
  * arcs SHARING a vertex: the 3-cycle on their 3 vertices covers both -> support
    nonempty for ALL n>=3, and empirically M!=0 for EVERY T (P=1) -- the signed
    sum never cancels for sharing arcs.
  * arcs DISJOINT: support nonempty iff n>=5 (a 5-cycle a-b-x-c-d). At n=4 NO odd
    cycle links two opposite K4-edges (only the EVEN 4-cycle) and two disjoint
    3-cycles need 6 vertices -> support EMPTY -> M==0 for all disjoint pairs. This
    is EXACTLY HYP-1165 (disjoint-arc Delta_2=0 at n=4, fails n>=5), now explained
    as orientation-free odd-cycle coverability. For n>=5 support is nonempty but
    M!=0 only for a fraction of T (cancellation), rising 0.52->0.69 as n: 5->6.
  s685's "disjoint <=> a COMMON odd cycle" was the realized n<=5 picture; the
  orientation-free support statement is the correct invariant form.

ARC-INTERACTION GRAPH A(T): vertices = the C(n,2) arcs, edges = {g,e : M(g,e)!=0}.
Tests whether its (labeled-by-score) isomorphism type is a FINER tournament
invariant than H (which collapses many non-isomorphic tournaments).

Session: claude-2026-06-06-S686 (arc-interaction-graph)."""
import sys; sys.stdout.reconfigure(line_buffering=True)
from itertools import product, combinations, permutations
from collections import Counter

def H(n, adj):
    size = 1 << n; dp = [[0]*n for _ in range(size)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(size):
        for v in range(n):
            c = dp[mask][v]
            if c:
                av = adj[v]
                for w in range(n):
                    if not (mask >> w & 1) and av >> w & 1: dp[mask | 1 << w][w] += c
    return sum(dp[size-1])

def E(n): return [(i, j) for i in range(n) for j in range(i+1, n)]

def build(n, bits):
    adj = [0]*n
    for (i, j), b in zip(E(n), bits):
        if b: adj[i] |= 1 << j
        else:  adj[j] |= 1 << i
    return adj

def coverable(n, gi, ei, es):
    """Orientation-FREE support test: can arcs es[gi], es[ei] be jointly covered by
    vertex-disjoint ODD cycles in K_n? (the Walsh support {S>={g,e}} is nonempty).
    This is a property of K_n and the arc positions, NOT of any orientation."""
    g = set(es[gi]); e = set(es[ei])
    if g & e:
        return True            # share a vertex -> the 3-cycle on the 3 used vertices covers both (n>=3)
    # disjoint (4 distinct vertices):
    if n >= 5:
        return True            # one odd 5-cycle a-b-x-c-d (a fifth vertex x) covers both edges
    return False               # n==4: only the EVEN 4-cycle links two opposite K4-edges -> support empty

print("=== PART 1: orientation-free SUPPORT is necessary for M(g,e) != 0 ===")
print("    (M != 0  =>  g,e jointly coverable by disjoint odd cycles in K_n)")
for n in [4, 5, 6]:
    es = E(n); m = len(es)
    share_nz = [0, 0]; disj_nz = [0, 0]
    violations = 0; checked = 0        # M!=0 but support empty  (predicted impossible)
    allbits = product([0, 1], repeat=m)
    cap = 4000 if n == 6 else 10**9    # n=6: sample for speed
    cnt = 0
    for bits in allbits:
        cnt += 1
        if cnt > cap: break
        adj = build(n, bits); H0 = H(n, adj)
        bl = list(bits); Hs = [0]*m
        for k in range(m): bl[k] ^= 1; Hs[k] = H(n, build(n, bl)); bl[k] ^= 1
        for gi in range(m):
            for ei in range(gi+1, m):
                bl[gi] ^= 1; bl[ei] ^= 1; Hge = H(n, build(n, bl)); bl[gi] ^= 1; bl[ei] ^= 1
                M = H0 - Hs[gi] - Hs[ei] + Hge
                nz = (M != 0)
                shares = bool(set(es[gi]) & set(es[ei]))
                if shares: share_nz[0] += nz; share_nz[1] += 1
                else:      disj_nz[0] += nz; disj_nz[1] += 1
                cov = coverable(n, gi, ei, es)
                checked += 1
                if nz and not cov: violations += 1   # M!=0 with empty support: must be 0
    tag = f"(all {cnt-1} T)" if cnt-1 <= cap else f"(first {cap} T)"
    print(f"n={n} {tag}: P(M!=0|share)={share_nz[0]/max(1,share_nz[1]):.3f}  "
          f"P(M!=0|disjoint)={disj_nz[0]/max(1,disj_nz[1]):.3f}")
    print(f"        necessary-condition violations (M!=0 but support empty): "
          f"{violations}/{checked}  "
          f"[support empty <=> disjoint arcs at n=4; matches HYP-1165]")

print("\n=== PART 2: is the arc-interaction graph A(T) a finer invariant than H? ===")
# At n=5: enumerate all tournaments, group by H, then within each H-class ask whether
# the (sorted degree sequence of A(T)) distinguishes non-isomorphic ones.
def arc_interaction_degseq(n, bits):
    es = E(n); m = len(es); adj = build(n, bits); H0 = H(n, adj)
    Hs = [0]*m; bl = list(bits)
    for k in range(m): bl[k] ^= 1; Hs[k] = H(n, build(n, bl)); bl[k] ^= 1
    deg = [0]*m
    for gi in range(m):
        for ei in range(gi+1, m):
            bl[gi] ^= 1; bl[ei] ^= 1; Hge = H(n, build(n, bl)); bl[gi] ^= 1; bl[ei] ^= 1
            if H0 - Hs[gi] - Hs[ei] + Hge != 0:
                deg[gi] += 1; deg[ei] += 1
    return H0, tuple(sorted(deg))

n = 5; es = E(n); m = len(es)
byH = Counter(); byHdeg = Counter()
for bits in product([0, 1], repeat=m):
    H0, ds = arc_interaction_degseq(n, bits)
    byH[H0] += 1
    byHdeg[(H0, ds)] += 1
# Compare resolution: # distinct H-values vs # distinct (H, A-degseq) classes
print(f"n=5: {sum(byH.values())} labeled tournaments")
print(f"  distinct H-values:                 {len(byH)}")
print(f"  distinct (H, A(T)-degree-seq):      {len(byHdeg)}")
print(f"  => A(T) degree sequence {'REFINES' if len(byHdeg) > len(byH) else 'does NOT refine'} H.")
# show which H-values get split
splits = {}
for (h, ds), c in byHdeg.items(): splits.setdefault(h, []).append(c)
for h in sorted(splits):
    if len(splits[h]) > 1:
        print(f"     H={h}: split into {len(splits[h])} A(T)-classes (sizes {sorted(splits[h], reverse=True)})")
