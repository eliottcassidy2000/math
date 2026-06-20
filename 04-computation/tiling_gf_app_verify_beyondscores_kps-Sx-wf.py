#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL VERIFICATION of the "beyond-scores-to-OCF" application of THM-554.
(kps-Sx-wf, 2026-06-20)

Claims under test:
  (C1) E[c3] = (C(n,3) + (n-2))/4                 [PROVED]
  (C2) E[c5](n) = (1/160)n^5 - (1/16)n^4 + (9/32)n^3 - (7/8)n^2 + (147/80)n - 7/4
       claimed exact values n=5..12: 19/16,49/8,315/16,199/4,1727/16,1683/8,6055/16,1279/2
  (C3) E[alpha_2 triangle-pair part](n) = 0,0,0,15/16,93/16,173/8,495/8,2395/16 for n=3..10
       (== brute alpha_2 for n<=7)
  (C4) Full E[alpha_2] at n=8 = 173/8 + 447/32 = 1139/32
  (C5) c5/alpha_2 NOT score-determined (need full 2^F state); c3 IS score-determined.
  (C6) complement-invariance => half-tiling 2x is real and lossless.

Strategy: independent BRUTE enumeration over all 2^F tilings (n<=7, n=8 where feasible),
exact rationals (fractions.Fraction), compared to:
  - the Z-engine score census (for c3),
  - an INDEPENDENT per-subset linearity computation of E[c3], E[c5], E[alpha_2-tri].
Default skeptical. One mismatch => holds=False with (n, values).
"""
import sys, time
from collections import defaultdict, Counter
from itertools import product, combinations
from fractions import Fraction
from math import comb
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

LOG = []
def out(s=""):
    print(s); LOG.append(str(s))

# ---------------------------------------------------------------------------
# Engine (copied from tile_address_score_gf_engine_kps.py)
# ---------------------------------------------------------------------------
def tiles(n):
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]

def beta_step(dist, n):
    nd = defaultdict(int)
    for vec, cnt in dist.items():
        l = list(vec) + [0]; l[n - 1] += 1; nd[tuple(l)] += cnt
    dist = nd
    for b in range(1, n - 1):
        nd = defaultdict(int)
        for vec, cnt in dist.items():
            l0 = list(vec); l0[n - 1] += 1; nd[tuple(l0)] += cnt
            l1 = list(vec); l1[b - 1] += 1; nd[tuple(l1)] += cnt
        dist = nd
    return dist

def build_Z(N):
    d = {(0,): 1}
    for n in range(2, N + 1):
        d = beta_step(d, n)
    return d

# ---------------------------------------------------------------------------
# Build the adjacency of a tiling. Convention: base path n->n-1->...->1,
# tile (a,b) bit0 => a->b, bit1 => b->a.  adj[i][j]=1 means arc i->j.
# ---------------------------------------------------------------------------
def build_adj(n, T, bv):
    adj = [[0] * (n + 1) for _ in range(n + 1)]
    for k in range(n, 1, -1):
        adj[k][k - 1] = 1
    for (a, b), bit in zip(T, bv):
        if bit == 0:
            adj[a][b] = 1
        else:
            adj[b][a] = 1
    return adj

def scores(n, adj):
    return [sum(adj[i][j] for j in range(1, n + 1)) for i in range(1, n + 1)]

# directed k-cycles among a specific vertex subset (count of directed Hamiltonian
# cycles on that subset using the tournament restricted to it). For a tournament
# on k vertices the number of directed cycles spanning all k = number of cyclic
# orderings that are consistent. We count properly: number of directed cycles
# (as cyclic sequences) of length k visiting each of the k vertices once.
def count_directed_kcycles_on_subset(adj, subset):
    # subset: tuple of vertices. Count cyclic permutations v0->v1->...->v_{k-1}->v0
    # with all arcs present. Divide by k (rotations) and treat reversal as distinct
    # direction. A directed cycle is an equivalence class under rotation only.
    k = len(subset)
    if k < 3:
        return 0
    total = 0
    s0 = subset[0]
    rest = subset[1:]
    for perm in _perms(rest):
        seq = (s0,) + perm
        ok = True
        for i in range(k):
            u = seq[i]; w = seq[(i + 1) % k]
            if not adj[u][w]:
                ok = False; break
        if ok:
            total += 1
    return total  # fixing s0 already removes rotation overcount

def _perms(t):
    from itertools import permutations
    return permutations(t)

# c3 = total number of 3-cycles (directed triangles)
def c3_brute(n, adj):
    c = 0
    for trip in combinations(range(1, n + 1), 3):
        c += count_directed_kcycles_on_subset(adj, trip)
    return c

# c5 = total number of directed 5-cycles spanning a 5-subset
def c5_brute(n, adj):
    c = 0
    for sub in combinations(range(1, n + 1), 5):
        c += count_directed_kcycles_on_subset(adj, sub)
    return c

# alpha_2 here = number of DISJOINT PAIRS of directed cycles (vertex-disjoint).
# In repo: H = 1 + 2*alpha_1 + 4*alpha_2 + ...; alpha_1 = # odd directed cycles,
# alpha_2 = # unordered pairs of vertex-disjoint odd directed cycles.
# "triangle-pair part" = pairs of vertex-disjoint 3-cycles.
def alpha2_triangle_pair_brute(n, adj):
    # enumerate 3-cycles as vertex sets with multiplicity (a triangle vertex-set
    # carries 1 directed 3-cycle for a tournament). Count unordered disjoint pairs.
    tri = []  # list of (frozenset, count)
    for trip in combinations(range(1, n + 1), 3):
        cc = count_directed_kcycles_on_subset(adj, trip)
        if cc:
            tri.append((set(trip), cc))
    total = 0
    for i in range(len(tri)):
        for j in range(i + 1, len(tri)):
            if tri[i][0].isdisjoint(tri[j][0]):
                total += tri[i][1] * tri[j][1]
    return total

def alpha2_full_brute(n, adj):
    # all unordered pairs of vertex-disjoint odd directed cycles (length 3 and 5
    # relevant for n<=8). Build list of (vertexset, count) for each odd cycle.
    odd = []
    for L in (3, 5, 7):
        if L > n:
            break
        for sub in combinations(range(1, n + 1), L):
            cc = count_directed_kcycles_on_subset(adj, sub)
            if cc:
                odd.append((set(sub), cc))
    total = 0
    for i in range(len(odd)):
        for j in range(i + 1, len(odd)):
            if odd[i][0].isdisjoint(odd[j][0]):
                total += odd[i][1] * odd[j][1]
    return total

# ---------------------------------------------------------------------------
# BRUTE means over all tilings (exact Fraction)
# ---------------------------------------------------------------------------
def brute_means(n, want_alpha2_full=False):
    T = tiles(n); F = len(T); N = 1 << F
    sum_c3 = 0; sum_c5 = 0; sum_a2tri = 0; sum_a2full = 0
    score_census = Counter()
    c3_dist = Counter()
    for bv in product((0, 1), repeat=F):
        adj = build_adj(n, T, bv)
        sv = scores(n, adj)
        score_census[tuple(sorted(sv))] += 1
        c3 = comb(n, 3) - sum(comb(s, 2) for s in sv)
        # cross-check c3 via direct triangle count
        c3d = c3_brute(n, adj)
        assert c3 == c3d, (n, sv, c3, c3d)
        c3_dist[c3] += 1
        sum_c3 += c3
        sum_c5 += c5_brute(n, adj)
        sum_a2tri += alpha2_triangle_pair_brute(n, adj)
        if want_alpha2_full:
            sum_a2full += alpha2_full_brute(n, adj)
    res = {
        "Ec3": Fraction(sum_c3, N),
        "Ec5": Fraction(sum_c5, N),
        "Ea2tri": Fraction(sum_a2tri, N),
        "score_census": score_census,
        "c3_dist": c3_dist,
        "N": N,
    }
    if want_alpha2_full:
        res["Ea2full"] = Fraction(sum_a2full, N)
    return res

# ---------------------------------------------------------------------------
# INDEPENDENT per-subset-linearity computation of the MEANS.
# For E[c3]: each triple {i<j<k} contributes Pr[it is a directed 3-cycle].
# We compute that probability by enumerating the induced sub-tournament over the
# tiling distribution restricted to the 3 vertices: which arcs among i,j,k are
# base-path arcs (forced) vs free tiles (uniform 1/2 each, independent).
# ---------------------------------------------------------------------------
def induced_arc_status(n, subset):
    """For each ordered pair within subset, return 'fwd' if base-path arc (forced
    high->low consecutive), 'free' if a free tile, etc. Base path: arc v->v-1 for
    consecutive v. Free tile (a,b) with a-b>=2 is uniform. So between u>v in subset:
      if u-v==1: forced arc u->v.
      else: free, P(u->v)=1/2."""
    info = {}
    for (u, v) in combinations(sorted(subset, reverse=True), 2):  # u>v
        if u - v == 1:
            info[(u, v)] = 'forced'   # u->v always
        else:
            info[(u, v)] = 'free'
    return info

def prob_directed_cycle_on_subset(n, subset):
    """Exact probability (over uniform free tiles) that the induced tournament on
    `subset` forms exactly the directed-cycle structure giving a directed cycle.
    We compute E[# directed cycles spanning subset] exactly by averaging over the
    free-tile assignments among the subset (each free pair independent uniform)."""
    sub = tuple(sorted(subset))
    pairs = list(combinations(sorted(sub, reverse=True), 2))  # (u,v), u>v
    free_pairs = [p for p in pairs if p[0] - p[1] != 1]
    forced = {p: 1 for p in pairs if p[0] - p[1] == 1}  # u->v
    total = Fraction(0)
    nfree = len(free_pairs)
    for bits in product((0, 1), repeat=nfree):
        # build adjacency among subset
        adj = {}
        for p in pairs:
            u, v = p
            if p in forced:
                adj[(u, v)] = 1; adj[(v, u)] = 0
            else:
                bit = bits[free_pairs.index(p)]
                # bit0 => a->b i.e. u->v ; bit1 => v->u
                if bit == 0:
                    adj[(u, v)] = 1; adj[(v, u)] = 0
                else:
                    adj[(u, v)] = 0; adj[(v, u)] = 1
        # count directed cycles spanning subset
        adjm = lambda x, y: adj[(x, y)] if (x, y) in adj else 0
        cc = _count_dcycles(sub, adjm)
        total += cc
    return total / (1 << nfree)

def _count_dcycles(sub, adjfn):
    from itertools import permutations
    k = len(sub); s0 = sub[0]; rest = sub[1:]; c = 0
    for perm in permutations(rest):
        seq = (s0,) + perm; ok = True
        for i in range(k):
            if not adjfn(seq[i], seq[(i + 1) % k]):
                ok = False; break
        if ok:
            c += 1
    return c

def Ec3_persubset(n):
    return sum(prob_directed_cycle_on_subset(n, s)
               for s in combinations(range(1, n + 1), 3))

def Ec5_persubset(n):
    return sum(prob_directed_cycle_on_subset(n, s)
               for s in combinations(range(1, n + 1), 5))

def Ea2tri_persubset(n):
    """E[# disjoint pairs of 3-cycles] = sum over unordered disjoint triple-pairs
    of P(both directed 3-cycles). The two triples are vertex-disjoint, but their
    induced free tiles are INDEPENDENT only if no free tile connects them -- but a
    free tile connecting the two triples does NOT affect whether each triple is a
    cycle (a 3-cycle depends only on the 3 internal arcs). So independence across
    the two triples' internal arcs holds; the cross arcs are irrelevant. Hence
    P(both) = P(tri1 cycle)*P(tri2 cycle)."""
    trips = list(combinations(range(1, n + 1), 3))
    pc = {t: prob_directed_cycle_on_subset(n, t) for t in trips}
    total = Fraction(0)
    for i in range(len(trips)):
        for j in range(i + 1, len(trips)):
            if set(trips[i]).isdisjoint(trips[j]):
                total += pc[trips[i]] * pc[trips[j]]
    return total

# ---------------------------------------------------------------------------
# Claimed closed forms
# ---------------------------------------------------------------------------
def Ec3_formula(n):
    return Fraction(comb(n, 3) + (n - 2), 4)

def Ec5_formula(n):
    n = Fraction(n)
    return (Fraction(1,160)*n**5 - Fraction(1,16)*n**4 + Fraction(9,32)*n**3
            - Fraction(7,8)*n**2 + Fraction(147,80)*n - Fraction(7,4))

CLAIM_Ec5 = {5:Fraction(19,16),6:Fraction(49,8),7:Fraction(315,16),8:Fraction(199,4),
             9:Fraction(1727,16),10:Fraction(1683,8),11:Fraction(6055,16),12:Fraction(1279,2)}
CLAIM_Ea2tri = {3:Fraction(0),4:Fraction(0),5:Fraction(0),6:Fraction(15,16),
                7:Fraction(93,16),8:Fraction(173,8),9:Fraction(495,8),10:Fraction(2395,16)}

# ---------------------------------------------------------------------------
def main():
    out("="*78)
    out("ADVERSARIAL VERIFICATION: beyond-scores-to-OCF (THM-554 application)")
    out("="*78)

    # ---- (A) E[c3] formula vs per-subset vs Z-engine vs brute -------------
    out("\n[A] E[c3]: formula (C(n,3)+(n-2))/4  vs  per-subset  vs  Z-engine  vs  brute")
    ok_A = True
    for n in range(3, 8):
        f = Ec3_formula(n)
        ps = Ec3_persubset(n)
        Z = build_Z(n)
        Ntil = sum(Z.values())
        sZ = Fraction(sum((comb(n,3)-sum(comb(s,2) for s in v))*c for v,c in Z.items()), Ntil)
        br = brute_means(n)["Ec3"]
        match = (f == ps == sZ == br)
        ok_A &= match
        out(f"  n={n}: formula={f}  persub={ps}  Zengine={sZ}  brute={br}  -> {'OK' if match else 'MISMATCH'}")
    # extend formula==persubset to higher n (no brute)
    for n in range(8, 11):
        f = Ec3_formula(n); ps = Ec3_persubset(n)
        out(f"  n={n}: formula={f}  persub={ps}  -> {'OK' if f==ps else 'MISMATCH'}")
        ok_A &= (f == ps)

    # ---- (B) E[c5] formula vs claimed values vs per-subset vs brute ------
    out("\n[B] E[c5]: poly formula vs claimed list vs per-subset vs brute")
    ok_B = True; ok_B_claimlist = True
    for n in range(5, 13):
        f = Ec5_formula(n)
        cl = CLAIM_Ec5[n]
        m1 = (f == cl)
        ok_B_claimlist &= m1
        ps = Ec5_persubset(n)
        m2 = (f == ps)
        line = f"  n={n}: formula={f}  claimed={cl}  persub={ps}"
        flags = []
        if not m1: flags.append("FORMULA!=CLAIMLIST")
        if not m2: flags.append("FORMULA!=PERSUB")
        out(line + ("  -> " + (", ".join(flags) if flags else "OK")))
        ok_B &= (m1 and m2)
    out("  brute c5 check:")
    for n in range(5, 8):
        br = brute_means(n)["Ec5"]; f = Ec5_formula(n)
        m = (br == f)
        ok_B &= m
        out(f"    n={n}: brute={br}  formula={f}  -> {'OK' if m else 'MISMATCH'}")

    # ---- (C) E[alpha_2 triangle-pair] -----------------------------------
    out("\n[C] E[alpha_2 triangle-pair]: claimed list vs per-subset vs brute")
    ok_C = True
    for n in range(3, 11):
        cl = CLAIM_Ea2tri[n]
        ps = Ea2tri_persubset(n)
        m1 = (cl == ps)
        out(f"  n={n}: claimed={cl}  persub={ps}  -> {'OK' if m1 else 'MISMATCH'}")
        ok_C &= m1
    out("  brute alpha2-triangle-pair check:")
    for n in range(3, 8):
        br = brute_means(n)["Ea2tri"]; ps = Ea2tri_persubset(n)
        m = (br == ps)
        ok_C &= m
        out(f"    n={n}: brute={br}  persub={ps}  -> {'OK' if m else 'MISMATCH'}")

    # ---- (C4) full alpha_2 at n=8 = 1139/32 -----------------------------
    out("\n[C4] full E[alpha_2] at n=8 (3x3 + 3x5 disjoint cross). Brute (heavy)...")
    t0 = time.time()
    try:
        r8 = brute_means(8, want_alpha2_full=True)
        full = r8["Ea2full"]; tri = r8["Ea2tri"]
        claim_full = Fraction(1139,32); claim_tri = Fraction(173,8); claim_cross = Fraction(447,32)
        out(f"    n=8 brute: Ea2tri={tri} (claim {claim_tri})  Ea2full={full} (claim {claim_full})")
        out(f"    cross = full-tri = {full-tri}  (claim {claim_cross})")
        ok_C4 = (full == claim_full and tri == claim_tri and (full-tri)==claim_cross)
        out(f"    -> {'OK' if ok_C4 else 'MISMATCH'}   ({time.time()-t0:.1f}s)")
    except Exception as e:
        ok_C4 = None
        out(f"    n=8 full brute SKIPPED/failed: {e}")

    # ---- (D) score-determinacy: c3 yes, c5/alpha_2 no -------------------
    out("\n[D] Score-determinacy test: does score vector DETERMINE c3 / c5 / alpha2?")
    for n in (6, 7):
        T = tiles(n); F = len(T)
        by_score_c3 = defaultdict(set)
        by_score_c5 = defaultdict(set)
        by_score_a2 = defaultdict(set)
        for bv in product((0,1), repeat=F):
            adj = build_adj(n, T, bv)
            sv = tuple(sorted(scores(n, adj)))
            by_score_c3[sv].add(comb(n,3)-sum(comb(s,2) for s in scores(n,adj)))
            by_score_c5[sv].add(c5_brute(n, adj))
            by_score_a2[sv].add(alpha2_triangle_pair_brute(n, adj))
        c3_det = all(len(s)==1 for s in by_score_c3.values())
        c5_det = all(len(s)==1 for s in by_score_c5.values())
        a2_det = all(len(s)==1 for s in by_score_a2.values())
        out(f"  n={n}: c3 score-determined={c3_det}  c5 score-determined={c5_det}  "
            f"alpha2tri score-determined={a2_det}")
        # show one counterexample for c5 if not determined
        if not c5_det:
            for sv, vals in by_score_c5.items():
                if len(vals) > 1:
                    out(f"    c5 counterexample at score {sv}: c5 in {sorted(vals)}")
                    break
        if not a2_det:
            for sv, vals in by_score_a2.items():
                if len(vals) > 1:
                    out(f"    alpha2 counterexample at score {sv}: in {sorted(vals)}")
                    break

    # ---- (E) complement halving (2x lossless) ---------------------------
    out("\n[E] Complement-invariance: c3(T)==c3(T^op), c5, alpha2 (per-tiling check)")
    for n in (5, 6):
        T = tiles(n); F = len(T); bad = 0; checked = 0
        for bv in product((0,1), repeat=F):
            adj = build_adj(n, T, bv)
            # complement tournament: reverse ALL arcs
            adjc = [[0]*(n+1) for _ in range(n+1)]
            for i in range(1,n+1):
                for j in range(1,n+1):
                    adjc[i][j] = adj[j][i]
            c3a = c3_brute(n, adj); c3b = c3_brute(n, adjc)
            c5a = c5_brute(n, adj); c5b = c5_brute(n, adjc)
            a2a = alpha2_triangle_pair_brute(n, adj); a2b = alpha2_triangle_pair_brute(n, adjc)
            if (c3a,c5a,a2a) != (c3b,c5b,a2b):
                bad += 1
            checked += 1
        out(f"  n={n}: {checked} tilings, complement-invariance violations={bad} "
            f"-> {'LOSSLESS 2x OK' if bad==0 else 'BROKEN'}")

    out("\n" + "="*78)
    out("SUMMARY")
    out(f"  [A] E[c3] formula==persub==Zengine==brute : {ok_A}")
    out(f"  [B] E[c5] formula==claimlist==persub==brute : {ok_B} (claimlist {ok_B_claimlist})")
    out(f"  [C] E[alpha2-tri] claim==persub==brute      : {ok_C}")
    out(f"  [C4] full alpha2 n=8 == 1139/32             : {ok_C4 if 'ok_C4' in dir() else 'n/a'}")
    out("="*78)

if __name__ == "__main__":
    main()
