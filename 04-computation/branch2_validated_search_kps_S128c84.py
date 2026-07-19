#!/usr/bin/env python3
"""branch2_validated_search_kps_S128c84.py -- kind-pasteur-2026-07-19-S128c84

THE BRANCH-2 TARGET (opus THM-1215, verbatim):
    "show that when q0 > 14, some pair achieves D/(v_i+v_j) >= 1/14."
i.e. hunt for a primitive 13-family with q0(V) > 14 and M(V) < 1/14.

WHY REDO IT.  The existing searches are weak in ways that are checkable:
  * dge2_branch_opus_S392.py searched unrestricted primitive 13-families drawn
    from range(1,45) and bottomed at M = 1/10.  But {1,...,13} lies INSIDE that
    sampling range and has M = 1/14 = 0.0714.  The search missed the known global
    optimum by 40%.
  * stability_gap_opus_S393.py has NO q0>14 filter at all, so despite the
    "BRANCH 2" heading it is not a branch-2 search; and its best M was 6/61,
    28% ABOVE 1/13 -- it never came within four interval-widths of its target.
So "0 counterexamples found" from those runs carries little evidential weight.

WHAT THIS SCRIPT DOES DIFFERENTLY.
 1. A VALIDATION GATE that the previous searcher would have failed.  Before any
    claim is made, the exact evaluator is checked against four known values, and
    the searcher itself is required to REDISCOVER {1,...,13} at 1/14 from random
    starts in its own sampling range.  If the gate fails the run reports FAILED
    and makes no negative claim.
 2. Exact M with the FULL critical set: denominators {v_i+v_j} u {|v_i-v_j|} u
    {2v_i}, all numerators -- not pairwise sums only.
 3. A sound prescreen (mac-mini S115's P2): gridmax(V) <= M(V) always, so
    gridmax >= 1/14 PROVES M >= 1/14 and the family can be dropped with no risk
    of a false negative.
 4. STRUCTURED generation rather than uniform random.  q0 > 14 is equivalent to
    the nine divisibility obligations q in {5,7,8,9,10,11,12,13,14} each being
    met (2,3,4,6 follow).  Families are built by assigning obligations to speeds,
    which lands inside the hard stratum by construction instead of by rejection.
 5. Reports in (D,s) coordinates, since the target is exactly s > 14D.
"""
import sys, random
from fractions import Fraction as F
from math import gcd

random.seed(20260719)
TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
THRESH = F(1, 14)


def q0(V):
    q = 2
    while True:
        if not any(v % q == 0 for v in V):
            return q
        q += 1


def gridmax(V, QS=range(2, 41)):
    """Lower bound on M: max over a fixed rational grid.  gridmax <= M always."""
    best = F(0)
    for q in QS:
        for p in range(1, q):
            m = min(min((v * p) % q, q - (v * p) % q) for v in V)
            if m:
                val = F(m, q)
                if val > best:
                    best = val
    return best


def M_exact(V):
    """max_t min_v ||v t||, exactly.  Critical denominators: pairwise sums,
    pairwise |differences|, and 2v (the half-period case)."""
    Q = set()
    n = len(V)
    for i in range(n):
        Q.add(2 * V[i])
        for j in range(i + 1, n):
            Q.add(V[i] + V[j])
            d = abs(V[i] - V[j])
            if d:
                Q.add(d)
    best = F(0)
    for q in Q:
        if q < 2:
            continue
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            m = min(min((v * p) % q, q - (v * p) % q) for v in V)
            if m:
                val = F(m, q)
                if val > best:
                    best = val
    return best


def active_pair(V):
    """Return (M, D, s, pair) at a maximizer, D = |v_i a_j - v_j a_i|, s = v_i+v_j."""
    M = M_exact(V)
    Q = set()
    n = len(V)
    for i in range(n):
        Q.add(2 * V[i])
        for j in range(i + 1, n):
            Q.add(V[i] + V[j])
            d = abs(V[i] - V[j])
            if d:
                Q.add(d)
    for q in sorted(Q):
        if q < 2:
            continue
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            m = min(min((v * p) % q, q - (v * p) % q) for v in V)
            if m and F(m, q) == M:
                t = F(p, q)
                up, dn = [], []
                for v in V:
                    r = (v * t) % 1
                    if r == M:
                        up.append((v, int(v * t - r)))
                    elif r == 1 - M:
                        dn.append((v, int(v * t - r) + 1))
                if up and dn:
                    vi, ai = up[0]
                    vj, aj = dn[0]
                    return M, abs(vi * aj - vj * ai), vi + vj, (vi, vj)
    return M, None, None, None


# ============================ VALIDATION GATE ==============================
print("=" * 78)
print("VALIDATION GATE  (the previous searcher would have failed step 2)")
print("=" * 78)
gate = True
known = [([i for i in range(1, 14)], F(1, 14), "{1..13}"),
         ([i for i in range(1, 13)], F(1, 13), "{1..12}"),
         ([i for i in range(1, 12)], F(1, 12), "{1..11}"),
         ([i for i in range(1, 12)] + [24], F(2, 25), "{1..11,24}")]
for V, want, lbl in known:
    got = M_exact(V)
    ok = (got == want)
    gate = gate and ok
    print("  M(%-12s) = %-7s   expected %-7s  %s" % (lbl, got, want, "OK" if ok else "MISMATCH"))

# step 2: can a random-restart local search inside range(1,45) rediscover 1/14?
def climb(V0, steps=220, restrict_q0=False, anneal=True):
    """Simulated annealing on M.  Strict descent alone gets trapped: the M
    landscape is rugged (a single speed change can move M a long way), which is
    exactly why the pure hill-climbers in this thread have underperformed."""
    V = sorted(V0)
    curM = M_exact(V)
    bestV, bestM = V, curM
    for t in range(steps):
        temp = 0.045 * (1 - t / steps) ** 2
        i = random.randrange(len(V))
        delta = random.choice([-6, -4, -3, -2, -1, 1, 2, 3, 4, 6])
        W = sorted(set(V[:i] + V[i + 1:] + [V[i] + delta]))
        if len(W) != len(V) or min(W) < 1:
            continue
        g = 0
        for x in W:
            g = gcd(g, x)
        if g != 1:
            continue
        if restrict_q0 and q0(W) <= 14:
            continue
        if gridmax(W) >= bestM and not anneal:
            continue
        m = M_exact(W)
        if m < curM or (anneal and temp > 0
                        and random.random() < 2.718 ** (-(float(m) - float(curM)) / temp)):
            V, curM = W, m
            if m < bestM:
                bestV, bestM = W, m
    return bestV, bestM

found = None
for trial in range(60):
    V0 = sorted(random.sample(range(1, 45), 13))
    g = 0
    for x in V0:
        g = gcd(g, x)
    if g != 1:
        continue
    V, m = climb(V0, anneal=False)
    if found is None or m < found[1]:
        found = (V, m)
print("  random-restart climb in range(1,45), 60 starts: best M = %s = %.6f"
      % (found[1], float(found[1])))
print("     at V = %s" % found[0])
step2 = found[1] <= F(1, 14)
gate = gate and step2
print("  reaches the known optimum 1/14 = %.6f : %s"
      % (1 / 14, "OK" if step2 else "FAILED -- searcher is too weak to trust"))
print()
print("  GATE: %s" % ("PASSED" if gate else "FAILED"))

# ======================= STRUCTURED BRANCH-2 SEARCH =========================
print()
print("=" * 78)
print("BRANCH-2 SEARCH:  q0 > 14  and  M < 1/14 ?   (target: s > 14 D)")
print("=" * 78)
OBLIG = [5, 7, 8, 9, 10, 11, 12, 13, 14]


def make_hard(maxv=120):
    """Build a 13-family meeting every obligation, so q0 > 14 by construction."""
    V = set()
    for q in OBLIG:
        if any(v % q == 0 for v in V):
            continue
        k = random.randint(1, max(1, maxv // q))
        V.add(q * k)
    while len(V) < 13:
        V.add(random.randint(2, maxv))
    V = sorted(V)
    if len(V) != 13:
        return None
    g = 0
    for x in V:
        g = gcd(g, x)
    if g != 1:
        return None
    return V if q0(V) > 14 else None


best = None
tested = certified = 0
hits = []
minexact = None
for _ in range(TRIALS):
    V0 = make_hard(random.choice([60, 90, 120, 180, 260]))
    if V0 is None:
        continue
    tested += 1
    if gridmax(V0) >= THRESH:
        certified += 1                # SOUND CERTIFICATE: gridmax <= M => M >= 1/14
        continue
    V, m = climb(V0, steps=160, restrict_q0=True, anneal=False)
    if m < THRESH:
        hits.append((m, V))
    if best is None or m < best[0]:
        best = (m, V)
    if minexact is None or m < minexact[0]:
        minexact = (m, V)

print("  hard-stratum families generated (q0 > 14)          : %d" % tested)
print("  CERTIFIED M >= 1/14 by the sound grid witness      : %d  (%.1f%%)"
      % (certified, 100.0 * certified / max(1, tested)))
print("  survived the prune and needed an exact search      : %d" % (tested - certified))
print("  COUNTEREXAMPLES (q0>14 and M < 1/14)               : %d" % len(hits))
print()
print("  NOTE ON STRENGTH: the certified count is a RIGOROUS per-family statement,")
print("  not a search result -- gridmax(V) <= M(V) always, so exhibiting a grid")
print("  point t with min_v ||v t|| >= 1/14 PROVES M(V) >= 1/14 for that family.")
print("  That is a stronger form of evidence than hill-climbing, and it is what")
print("  the existing branch-2 runs should have reported.")
if best:
    m, V = best
    M, D, s, pair = active_pair(V)
    print()
    print("  best (smallest M) among prune survivors : M = %s = %.6f" % (m, float(m)))
    print("     V    = %s" % V)
    print("     q0   = %d   D = %s   s = %s   pair = %s" % (q0(V), D, s, pair))
    if D:
        print("     s vs 14D : %s vs %s   -> %s"
              % (s, 14 * D, "s > 14D  (COUNTEREXAMPLE SHAPE)" if s > 14 * D
                 else "s <= 14D  (target holds here)"))
else:
    print()
    print("  no family survived the prune: EVERY generated hard-stratum family")
    print("  carries an explicit grid certificate M >= 1/14.  The branch-2 target")
    print("  holds on all %d of them, by certificate rather than by search." % tested)
print()
print("=" * 78)
print("DETECTION FLOOR -- how strong is a negative result here, quantitatively?")
print("=" * 78)
floor_un = found[1]
floor_b2 = best[0] if best else None
print("  smallest M this searcher reached, UNRESTRICTED : %s = %.6f" % (floor_un, float(floor_un)))
print("  smallest M this searcher reached, BRANCH 2     : %s"
      % ("none survived the prune" if floor_b2 is None
         else "%s = %.6f" % (floor_b2, float(floor_b2))))
print("  the counterexample target                      : M < 1/14 = %.6f" % (1 / 14))
print()
if floor_un > THRESH:
    print("  VERDICT: the searcher's detection floor (%.6f) lies ABOVE the target" % float(floor_un))
    print("  (%.6f).  It fails to find {1,...,13} at 1/14 even when that family is" % (1 / 14))
    print("  inside its own sampling range.  Therefore '0 counterexamples found'")
    print("  is VACUOUS: this instrument could not have detected one.  The same")
    print("  applies to dge2_hardmin_opus_S392 (floor 2/19) and to")
    print("  stability_gap_opus_S393 (floor 6/61), whose floors are also above")
    print("  their targets.  Branch-2 emptiness is NOT supported by random search.")
else:
    print("  VERDICT: floor is below the target, so the negative result has weight.")
