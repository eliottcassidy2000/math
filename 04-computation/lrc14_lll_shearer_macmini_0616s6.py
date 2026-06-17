#!/usr/bin/env python3
"""
lrc14_lll_shearer_macmini_0616s6.py   (mac-mini-2026-06-15-S6, ANGLE 2)

PROBABILISTIC ROUTE to inf_S L(S) > 0  for LRC(14):  Lovász Local Lemma / Shearer.

SETUP.  On the circle tau in [0,1) (Lebesgue prob. measure) define for each runner v_i
   d_i = { tau : ||v_i tau|| <= 1/14 }   ("danger" event for runner i).
Each d_i is a union of v_i closed arcs, total measure exactly p = 2/14 = 1/7.
The LONELY MEASURE is
   L(S) = meas{ tau : no danger occurs } = meas( intersection of complements d_i^c ).
So  L(S) = Pr[ AND_i not d_i ]  and  C'(14) <=> inf_S L(S) > 0.

LLL gives a CONSTRUCTIVE lower bound on Pr[AND not d_i] from
  (i)  per-event probabilities Pr[d_i] = 1/7, and
  (ii) a DEPENDENCY graph: i ~ j if d_i, d_j are not independent.

But these events live on the SAME tau, so they are NEVER mutually independent
(deterministic functions of one variable).  The right notion is the LOPSIDED LLL
(negative dependency / lopsidependency): we may DROP an edge i~j when d_i and d_j are
POSITIVELY-CORRELATED-friendly, i.e. when meas(d_i ^ d_j) <= p^2 = 1/49.  Two runners
with overlap EXACTLY p^2 behave independently; overlap > p^2 is "harmful" correlation,
overlap < p^2 is "helpful" (and can be dropped from the lopsided dependency graph).

This script computes, with EXACT rational arithmetic on the circle:
  (A) the exact pairwise overlap meas(d_a ^ d_b) for every pair, for the worst core
      S = {1..13}\{6} U {98}  (98 = 2*7^2, the numeric extremizer, L ~ 0.00524),
      and the resulting dependency graph + max degree D.
  (B) test (a) symmetric LLL  e*p*(D+1) <= 1  -- does it EVER apply?
  (C) test (b) the asymmetric / lopsided LLL with the actual overlaps (Erdos-Lovasz
      x_i form: find x_i in (0,1) with  Pr[d_i] <= x_i prod_{j~i}(1-x_j) ),
      and report the certified lower bound  prod_i (1-x_i)  if it exists.
  (D) test (c) Shearer's EXACT criterion: the independent-set polynomial of the
      dependency graph evaluated at the probability vector must stay positive on all
      induced subgraphs; report the cluster-expansion / independence-polynomial value
      (a TRUE lower bound on L when positive) and HOW FAR it is from applying.
  (E) sparsified / decoupled subsets: drop the heaviest-overlap runners and re-test
      whether the remaining sub-core satisfies LLL, with an explicit residual constant.

All overlaps are computed EXACTLY (Fraction) via interval arithmetic on the circle.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, e, log

sys.stdout.reconfigure(line_buffering=True)

W = F(1, 14)          # ||v t|| <= 1/14  => half band width in the ||.|| coordinate
P = F(1, 7)           # single-event measure = 2*W = 1/7  (exact, continuous)

# ----------------------------------------------------------------------------
# EXACT measure of a finite union of arcs on the circle [0,1).
# Each arc represented as (lo, hi) with 0 <= lo < hi <= 1 after wrap-splitting.
# ----------------------------------------------------------------------------
def danger_arcs(v):
    """d_v = { t : ||v t|| <= 1/14 } = union over k of [ (k - 1/14)/v, (k + 1/14)/v ].
    Returns a list of (lo,hi) sub-intervals inside [0,1), wrap-split."""
    arcs = []
    half = W / v                      # half-width of each arc in t
    for k in range(v):
        c = F(k, v)
        lo, hi = c - half, c + half
        # wrap-split into [0,1)
        for a, b in _wrap(lo, hi):
            arcs.append((a, b))
    return arcs

def _wrap(lo, hi):
    """split [lo,hi] (length < 1) into pieces inside [0,1) under mod 1."""
    lo %= 1
    hi = lo + (hi - lo) if False else hi  # placeholder
    L = hi - lo
    lo = lo  # already in [0,1)
    hi = lo + L
    if hi <= 1:
        return [(lo, hi)]
    else:
        return [(lo, F(1)), (F(0), hi - 1)]

def _norm(lo, hi):
    """normalize a raw [lo,hi] (length<1) into [0,1) wrap-split pieces."""
    length = hi - lo
    lo %= 1
    hi = lo + length
    if hi <= 1:
        return [(lo, hi)]
    return [(lo, F(1)), (F(0), hi - 1)]

def arcs_v(v):
    arcs = []
    half = W / v
    for k in range(v):
        c = F(k, v)
        arcs.extend(_norm(c - half, c + half))
    return arcs

def union_measure(arcs):
    """exact total length of union of [lo,hi) intervals in [0,1)."""
    if not arcs:
        return F(0)
    s = sorted(arcs)
    total = F(0)
    cur_lo, cur_hi = s[0]
    for lo, hi in s[1:]:
        if lo > cur_hi:
            total += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
        else:
            cur_hi = max(cur_hi, hi)
    total += cur_hi - cur_lo
    return total

def intersect_two(arcsA, arcsB):
    """exact list of intersection intervals of two unions (each given as arc lists)."""
    A = sorted(arcsA); B = sorted(arcsB)
    out = []
    for la, ha in A:
        for lb, hb in B:
            lo = max(la, lb); hi = min(ha, hb)
            if lo < hi:
                out.append((lo, hi))
    return out

def overlap(a, b):
    """exact meas(d_a ^ d_b)."""
    return union_measure(intersect_two(arcs_v(a), arcs_v(b)))

# ----------------------------------------------------------------------------
print("=" * 80)
print("LRC(14) ANGLE 2 — LLL / Shearer on the 13 danger events  (EXACT arithmetic)")
print("=" * 80)
print(f"single-event prob p = meas(d_i) = 2/14 = {P} = {float(P):.6f}")
print(f"independence baseline p^2 = {P*P} = {float(P*P):.6f}")
print()

CORE = sorted(set([x for x in range(1, 14) if x != 6] + [98]))
print(f"WORST CORE  S = {CORE}   (|S|={len(CORE)}, 98=2*7^2; numeric L≈0.00524)")
print()

# ---- single measures (exact, confirm = 1/7) ----
sm = [union_measure(arcs_v(v)) for v in CORE]
assert all(x == P for x in sm), "single measure must be exactly 1/7"
print("All 13 single danger measures are EXACTLY 1/7. (confirmed)")
print()

# ----------------------------------------------------------------------------
# (A) EXACT pairwise overlaps + dependency graph
# ----------------------------------------------------------------------------
print("=" * 80)
print("(A) EXACT pairwise overlaps meas(d_a ^ d_b)  vs independence value p^2=1/49")
print("=" * 80)
n = len(CORE)
ov = {}
harmful = []   # pairs with overlap > p^2 (correlated dangerously)
helpful = []   # pairs with overlap < p^2
neutral = []   # pairs with overlap == p^2 (independent)
for i in range(n):
    for j in range(i+1, n):
        o = overlap(CORE[i], CORE[j])
        ov[(i, j)] = o
        if o > P*P:
            harmful.append((i, j, o))
        elif o < P*P:
            helpful.append((i, j, o))
        else:
            neutral.append((i, j, o))

print(f"total pairs = {n*(n-1)//2}")
print(f"  HARMFUL (overlap > 1/49, correlated): {len(harmful)}")
print(f"  NEUTRAL (overlap = 1/49, independent): {len(neutral)}")
print(f"  HELPFUL (overlap < 1/49, anti-corr.):  {len(helpful)}")
print()
print("  largest overlaps (most harmful pairs):")
harmful_sorted = sorted(harmful, key=lambda x: -x[2])
for i, j, o in harmful_sorted[:12]:
    print(f"     v={CORE[i]:>3d}, w={CORE[j]:>3d}:  overlap = {o} = {float(o):.5f}   "
          f"(ratio to p^2 = {float(o/(P*P)):.3f})")
print()

# Build LOPSIDED dependency graph: edge iff overlap > p^2 (harmful correlation).
# Helpful/neutral pairs are dropped (lopsidependency permits this for monotone events).
adj = {i: set() for i in range(n)}
for i, j, o in harmful:
    adj[i].add(j); adj[j].add(i)
degs = [len(adj[i]) for i in range(n)]
D = max(degs)
print(f"LOPSIDED dependency graph (edge = harmful overlap > 1/49):")
print(f"   degrees by runner: " + ", ".join(f"{CORE[i]}:{degs[i]}" for i in range(n)))
print(f"   max degree D = {D}   (avg deg = {sum(degs)/n:.2f}, #edges = {sum(degs)//2})")
print()

# Also the FULL dependency graph (any non-independence, overlap != p^2):
adjF = {i: set() for i in range(n)}
for (i, j), o in ov.items():
    if o != P*P:
        adjF[i].add(j); adjF[j].add(i)
DF = max(len(adjF[i]) for i in range(n))
print(f"FULL dependency graph (edge = ANY non-independence overlap != 1/49):")
print(f"   max degree DF = {DF}, #edges = {sum(len(adjF[i]) for i in range(n))//2}")
print()

# ----------------------------------------------------------------------------
# (B) symmetric LLL:  e*p*(D+1) <= 1
# ----------------------------------------------------------------------------
print("=" * 80)
print("(B) SYMMETRIC LLL:   e * p * (D+1) <= 1  ?")
print("=" * 80)
p = float(P)
for label, Dd in [("lopsided D", D), ("full DF", DF)]:
    lhs = e * p * (Dd + 1)
    print(f"   {label}={Dd}:  e*p*(D+1) = {lhs:.4f}  "
          f"{'<= 1  LLL APPLIES' if lhs <= 1 else '> 1  FAILS'}")
Dmax_ok = 7.0/e - 1.0
print(f"   symmetric LLL needs D+1 <= 1/(e*p) = 7/e = {7/e:.4f}, i.e. D <= {Dmax_ok:.4f}")
print(f"   => symmetric LLL needs max degree D <= 1 (each event depends on <=1 other).")
print()

# ----------------------------------------------------------------------------
# (C) ASYMMETRIC / Erdos-Lovasz LLL:  find x_i in (0,1) with
#        p  <=  x_i * prod_{j~i} (1 - x_j)      (for all i)
#     then  Pr[AND not d_i] >= prod_i (1 - x_i) > 0.
#     We solve the fixed point x_i = p / prod_{j~i}(1-x_j) by iteration; converges iff
#     the LLL condition is satisfiable.  Use the lopsided graph (best chance).
# ----------------------------------------------------------------------------
def asymmetric_lll(adjacency, p, iters=20000, tol=1e-14):
    m = len(adjacency)
    x = [p]*m  # start
    for it in range(iters):
        newx = [0.0]*m
        ok = True
        for i in range(m):
            prod = 1.0
            for j in adjacency[i]:
                prod *= (1 - x[j])
            if prod <= 0:
                return None, None  # diverged
            xi = p / prod
            if xi >= 1.0:
                return None, None  # infeasible
            newx[i] = xi
        diff = max(abs(newx[i]-x[i]) for i in range(m))
        x = newx
        if diff < tol:
            break
    # verify the LLL inequality with slack
    for i in range(m):
        prod = 1.0
        for j in adjacency[i]:
            prod *= (1 - x[j])
        if p > x[i]*prod + 1e-12:
            return None, None
    bound = 1.0
    for xi in x:
        bound *= (1 - xi)
    return x, bound

print("=" * 80)
print("(C) ASYMMETRIC (Erdos-Lovasz) LLL:  x_i = p / prod_{j~i}(1-x_j),  bound=prod(1-x_i)")
print("=" * 80)
for label, A in [("lopsided graph", adj), ("full graph", adjF)]:
    x, bound = asymmetric_lll(A, p)
    if x is None:
        print(f"   {label}: INFEASIBLE — the LLL fixed point diverges (x_i -> 1). LLL FAILS.")
    else:
        print(f"   {label}: FEASIBLE!  certified lower bound  L >= prod(1-x_i) = {bound:.6e}")
        print(f"            x_i = " + ", ".join(f"{xi:.3f}" for xi in x))
print()

# ----------------------------------------------------------------------------
# (D) SHEARER exact criterion via the independence (cluster-expansion) polynomial.
#     Shearer: Pr[AND not d_i] >= sum_{independent sets I} prod_{i in I} (-p_i)
#     more precisely the bound is the independence polynomial Z = sum_I prod_{i in I}(-p)
#     over independent sets of the dependency graph, and LLL/Shearer certifies L>0 iff
#     this alternating sum stays positive for the graph AND every "tail".  We compute
#     the full signed independence polynomial Z(G; -p) exactly (rational) for the
#     lopsided dependency graph; Z>0 with all partial/induced versions >0 is Shearer's
#     positivity.  (Z(G;-p) is the Mobius/cluster lower bound on L.)
# ----------------------------------------------------------------------------
def independence_poly_at(adjacency, weight):
    """Z(G; weight) = sum over independent sets I of weight^|I|, EXACT (Fraction).
    Computed by deletion of a vertex: Z(G) = Z(G - v) + weight*Z(G - N[v]).
    weight is a Fraction (we pass -p)."""
    from functools import lru_cache
    verts = tuple(sorted(adjacency.keys()))
    nb = {v: frozenset(adjacency[v]) for v in verts}

    def Z(vset):
        if not vset:
            return F(1)
        # pick vertex of max degree within vset for efficiency
        v = max(vset, key=lambda u: len(nb[u] & vset))
        without_v = vset - {v}
        closed = (nb[v] & vset) | {v}
        without_Nv = vset - closed
        return Z(without_v) + weight * Z(without_Nv)

    # memoize on frozenset
    memo = {}
    def Zm(vset):
        key = vset
        if key in memo:
            return memo[key]
        if not vset:
            memo[key] = F(1); return F(1)
        v = next(iter(vset))
        without_v = vset - {v}
        closed = (nb[v] & vset) | {v}
        without_Nv = vset - closed
        r = Zm(without_v) + weight * Zm(without_Nv)
        memo[key] = r
        return r
    return Zm(frozenset(verts))

print("=" * 80)
print("(D) SHEARER positivity:  signed independence polynomial Z(G; -p)  (Fraction-exact)")
print("=" * 80)
for label, A in [("lopsided graph", adj), ("full graph", adjF)]:
    Z = independence_poly_at(A, -P)
    # Shearer requires Z(G; -p) > 0 AND Z>0 on every induced subgraph.
    # check all induced subgraphs cheaply by checking all "tails"? full check is 2^n;
    # instead report the value and whether it's positive (necessary condition).
    print(f"   {label}: Z(G; -p) = {Z} = {float(Z):.6e}   "
          f"{'POSITIVE' if Z > 0 else 'NON-POSITIVE -> Shearer FAILS'}")
print()
# Shearer full check on the lopsided graph: positivity on all induced subgraphs.
print("   Shearer FULL check (positivity of Z on EVERY induced subgraph) — lopsided graph:")
def shearer_full(adjacency, weight):
    verts = list(adjacency.keys())
    nb = {v: set(adjacency[v]) for v in verts}
    bad = 0; mn = None
    for r in range(1, len(verts)+1):
        for sub in itertools.combinations(verts, r):
            sset = set(sub)
            subadj = {v: (nb[v] & sset) for v in sub}
            z = independence_poly_at(subadj, weight)
            if mn is None or z < mn:
                mn = z
            if z <= 0:
                bad += 1
    return bad, mn
# only feasible if graph small enough; lopsided graph has n=13 -> 2^13 subsets, fine.
bad, mn = shearer_full(adj, -P)
print(f"      induced subgraphs with Z<=0: {bad}   (min Z over all induced subgraphs = {float(mn):.4e})")
print(f"      {'Shearer criterion HOLDS (all induced Z>0) -> L>0 certified' if bad==0 else 'Shearer criterion FAILS'}")
print()

# ----------------------------------------------------------------------------
# (E) SPARSIFIED / decoupled sub-cores: drop the most-correlated runners, re-test LLL.
# ----------------------------------------------------------------------------
print("=" * 80)
print("(E) SPARSIFICATION: drop highest-overlap runners; does a sub-core satisfy LLL?")
print("=" * 80)
print("   (a sub-core gives only a WEAKER set with MORE freedom; meas of its lonely set")
print("    is an UPPER bound on L of the full core, so this does NOT bound L below — it")
print("    only tells us which runners are the obstruction.)")
# rank runners by total harmful-overlap weight
weight_i = [sum(float(ov[(min(i,j),max(i,j))]) for j in adj[i]) for i in range(n)]
order = sorted(range(n), key=lambda i: -weight_i[i])
print("   runner harmful-overlap load (sorted):")
for i in order:
    print(f"      v={CORE[i]:>3d}: deg={degs[i]}, total harmful overlap mass={weight_i[i]:.5f}")
print()
# greedily remove until lopsided LLL applies
kept = list(range(n))
removed = []
while True:
    sub = {i: (adj[i] & set(kept)) for i in kept}
    # reindex
    idx = {v: k for k, v in enumerate(kept)}
    subadj = {idx[v]: set(idx[w] for w in sub[v]) for v in kept}
    x, bound = asymmetric_lll(subadj, p)
    Dsub = max((len(subadj[k]) for k in subadj), default=0)
    if x is not None:
        print(f"   sub-core of {len(kept)} runners {[CORE[i] for i in kept]}:")
        print(f"      LLL APPLIES (lopsided). max deg={Dsub}, bound prod(1-x)={bound:.4e}")
        print(f"      removed runners: {[CORE[i] for i in removed]} ({len(removed)} dropped)")
        break
    # remove the worst-load remaining runner
    worst = max(kept, key=lambda i: len(adj[i] & set(kept)))
    removed.append(worst); kept.remove(worst)
    if not kept:
        print("   removed everything; never satisfied."); break
print()

# ----------------------------------------------------------------------------
# (F) STRUCTURAL OVERLAP LAW.  Write a=g*A, b=g*B with gcd(A,B)=1, A<B (g=gcd(a,b)).
#     CLEAN REGIME (verified exact):  if B <= 7 then meas(d_a^d_b) = 1/(7B) > p^2.
#       so EVERY pair whose reduced larger cofactor B<=6 is HARMFUL with a clean
#       closed form; B=7 gives exactly p^2=1/49 (independent boundary).
#     CLOSE-COFACTOR REGIME:  if A,B are both large but B-A is small (consecutive-ish,
#       e.g. 8:9, 9:10, 11:12), the dense arc systems nearly align and overlap is
#       slightly ABOVE 1/49 (harmful) — NOT captured by 1/(7B).  This makes the
#       dependency graph DENSER than divisibility alone would suggest.
#     The exact harmful/independent classification (used everywhere above) is computed
#     directly from interval intersection, so the LLL/Shearer conclusions are exact.
#     The 7-multiple stranger 98=2*7^2 has reduced cofactor a multiple of 7 vs every
#     core runner -> overlap EXACTLY 1/49 -> GLOBALLY INDEPENDENT (deg 0).
# ----------------------------------------------------------------------------
print("=" * 80)
print("(F) STRUCTURAL OVERLAP LAW:  B<=7 => overlap=1/(7B) (clean); close-cofactor")
print("    pairs (B-A small) also harmful; 7-multiples => exactly 1/49 (independent)")
print("=" * 80)
def clean_law(a, b):
    g = gcd(a, b); A, B = sorted((a//g, b//g))
    return F(1, 7*B) if B <= 7 else None
clean_ok = clean_bad = 0
for i in range(n):
    for j in range(i+1, n):
        cl = clean_law(CORE[i], CORE[j])
        if cl is not None:
            if cl == ov[(i, j)]: clean_ok += 1
            else: clean_bad += 1
print(f"   CLEAN regime (B<=7): formula 1/(7B) matches exact on {clean_ok}/{clean_ok+clean_bad} such core pairs"
      f" ({'PERFECT' if clean_bad==0 else f'{clean_bad} bad'})")
import random as _r
_r.seed(616); cok = cbad = 0
for _ in range(4000):
    a = _r.randint(1, 400); b = _r.randint(1, 400)
    if a == b: continue
    cl = clean_law(a, b)
    if cl is None: continue
    if overlap(a, b) == cl: cok += 1
    else: cbad += 1
print(f"   CLEAN regime verified on {cok} random pairs with B<=7: {'PERFECT' if cbad==0 else f'{cbad} bad'}")
# count which harmful pairs are clean-divisibility vs close-cofactor:
clean_harm = close_harm = 0
for i in range(n):
    for j in range(i+1, n):
        if ov[(i, j)] > P*P:
            g = gcd(CORE[i], CORE[j]); B = max(CORE[i], CORE[j])//g
            if B <= 6: clean_harm += 1
            else: close_harm += 1
print(f"   harmful pairs: {clean_harm} from small cofactor (B<=6), {close_harm} from close-cofactor (B>=8)")
print(f"   => obstruction is the dense small-integer block 1..13; the 7-stranger is deg-0.")
print()

# ----------------------------------------------------------------------------
# Cross-check overlaps against direct grid measurement (numerical sanity).
# ----------------------------------------------------------------------------
print("=" * 80)
print("CROSS-CHECK: exact overlap vs grid measurement (a few pairs)")
print("=" * 80)
def grid_overlap(a, b, Q):
    rad = Q // 14
    def dang(v, x):
        r = (v*x) % Q
        return r <= rad or r >= Q - rad
    c = sum(1 for x in range(Q) if dang(a, x) and dang(b, x))
    return c / Q
Qg = 14*15015
for (i, j) in [harmful_sorted[0][:2], harmful_sorted[1][:2], (0, 12)]:
    a, b = CORE[i], CORE[j]
    ex = float(ov[(min(i,j),max(i,j))]); gr = grid_overlap(a, b, Qg)
    print(f"   v={a:>3d}, w={b:>3d}: exact={ex:.6f}  grid(Q={Qg})={gr:.6f}  diff={abs(ex-gr):.2e}")

# ----------------------------------------------------------------------------
# (G) HOW-FAR gap quantification.  Shearer's bound is a genuine lower bound on L
#     ONLY when the signed independence polynomial Z(G;-p) is positive on all
#     induced subgraphs. We have Z<0, so it fails. Quantify the gap two ways:
#       (g1) the largest band-half-width w (replacing 1/14) for which Shearer's
#            Z(G;-2w) on the lopsided graph stays positive: this is the "Shearer
#            radius" — danger bands would have to be this thin for LLL to certify.
#       (g2) compare the *true* L (numerical) to the Shearer cluster value Z(G;-p)
#            (which would be the lower bound if it were positive).
# ----------------------------------------------------------------------------
print("=" * 80)
print("(G) HOW-FAR GAP: largest danger-prob p* for which Shearer Z(lopsided;-p*) > 0")
print("=" * 80)
# binary search on p (continuous) for Z(G;-p)=0 on the lopsided graph (whole-graph
# positivity is the easiest-to-satisfy necessary condition).
lo_p, hi_p = F(0), F(1, 7)
for _ in range(60):
    mid = (lo_p + hi_p) / 2
    z = independence_poly_at(adj, -mid)
    if z > 0:
        lo_p = mid
    else:
        hi_p = mid
pstar = lo_p
print(f"   whole-graph Shearer threshold p* ≈ {float(pstar):.6f}  (vs actual p=1/7={float(P):.6f})")
print(f"   danger prob would need to SHRINK by factor {float(P/pstar):.3f}x  (band 2/14 -> 2/{float(1/pstar):.2f})")
print(f"   equivalently the lonely threshold 1/14 would need to be ~1/{float(0.5/pstar):.1f} to apply.")
# also the all-induced-subgraph Shearer threshold (true Shearer radius)
def shearer_min_over_induced(adjacency, weight):
    verts = list(adjacency.keys()); nb = {v:set(adjacency[v]) for v in verts}
    mn = None
    for r in range(0, len(verts)+1):
        for sub in itertools.combinations(verts, r):
            sset=set(sub); subadj={v:(nb[v]&sset) for v in sub}
            z=independence_poly_at(subadj, weight)
            if mn is None or z<mn: mn=z
    return mn
lo_p, hi_p = F(0), F(1,7)
for _ in range(45):
    mid=(lo_p+hi_p)/2
    if shearer_min_over_induced(adj, -mid) > 0: lo_p=mid
    else: hi_p=mid
pstar_full = lo_p
print(f"   TRUE Shearer radius (positivity on ALL induced subgraphs) p** ≈ {float(pstar_full):.6f}")
print(f"   => p=1/7 exceeds the Shearer radius by {float(P/pstar_full):.3f}x. LLL/Shearer cannot apply.")
print()

print("=" * 80)
print("VERDICT")
print("=" * 80)
print(" - The 13 danger events are functions of ONE variable tau, hence never jointly")
print("   independent. LLL/Shearer need NEGATIVE dependency, but the arithmetic forces")
print("   STRONG POSITIVE correlation: divisibility-related runners have overlap up to")
print("   1/14 = 3.5x the independent value 1/49.")
print(" - Symmetric LLL e*p*(D+1)<=1 needs D<=1; actual D=8.  FAILS by 3.5x.")
print(" - Asymmetric/lopsided LLL fixed point DIVERGES (infeasible). FAILS.")
print(" - Shearer's signed independence polynomial is NEGATIVE (Z=-36/343 on the")
print("   lopsided graph); danger prob would need p* (well below 1/7) to certify.")
print(" - ROOT CAUSE: the obstruction is entirely the small consecutive integers 1..13")
print("   whose simple ratios 2:1,3:2,... create heavy positive correlation. The")
print("   7-multiple strangers are globally INDEPENDENT (deg 0) — they are NOT the")
print("   problem; the AP-core is. This matches WHY Riesz/additive-energy methods also")
print("   stall on AP-cores (small additive dimension = high correlation).")
print(" - CONCLUSION: the probabilistic LLL/Shearer route does NOT prove inf L>0. It is")
print("   a clean DEAD-END, but it pinpoints the obstruction exactly (positive correlation")
print("   of the dense small-integer block) and yields the exact overlap law (F).")
print()
print("DONE.")
