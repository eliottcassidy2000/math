"""
lrc14_refute_pairwise-gap-load_M2diffspeedsigns_kps-S2-wf.py

ADVERSARIAL REFUTATION of the forbidden-class claim for METHOD M2_diff_speed_sign.

THE MAP (exact, copied from lrc14_tourmap_pairwise-gap-load):
  Vertices = the n runners (speeds), sorted ascending.
  Arc i->j  iff  snrm((v_i - v_j) * tau) > 0.
  snrm(x) in (-1/2, 1/2] is the SIGNED normalized residue (antisymmetric except at 1/2).
  A tournament is VALID/CLEAN iff for EVERY pair snrm((v_i-v_j)tau) not in {0, +/-1/2}
  (boundary ties make orientation degenerate). We DROP invalid configs.

THE CLAIM UNDER ATTACK:
  At n = 3, 4, 5, EVERY valid M2 tournament (over all crossing times tau) is
  SELF-COMPLEMENTARY (SC). Hence every NON-SC iso class is FORBIDDEN.
  (Breaks at n=6: witness S=(1,2,3,4,5,7), tau=5/14, non-SC.)

REFUTATION STRATEGY: realize a non-SC iso class at n in {3,4,5} with this EXACT
map over a genuinely broad, EXACT (Fraction) search:
  (A) Optimal lonely tau, all primitive speed sets up to large vmax.
  (B) ALL candidate crossing times tau (THM-524 binding-pair crossings + 1/2),
      which is the complete finite set where M2's signs can change -- so this is
      EXHAUSTIVE over the genuinely distinct M2 tournaments for a given S.
  (C) Off-grid tau: arbitrary rationals p/q between consecutive crossing times,
      AND a dense rational sweep, to catch any sign pattern between crossings.
  (D) Non-unit / non-coprime-to-14 speeds, sporadic sets, covering sets,
      tight-locus sets (AP {1..5}, Goddyn-Wong style), large speeds.
  (E) A STRUCTURAL exhaustion: M2 depends only on the signs of snrm(d*tau) over
      the pairwise differences d. We directly enumerate ALL achievable sign
      vectors for every speed set by walking tau across every crossing of every
      difference speed, capturing EVERY M2 tournament the map can produce.

All decisions use EXACT rationals. No floats anywhere.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations

# ---------------- validated M tool (verbatim) ----------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def snrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else r - 1
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at

# ---------------- tournament utilities ----------------
def score_seq(A):
    n = len(A); return tuple(sorted(sum(A[i]) for i in range(n)))
def count_3cycles(A):
    n = len(A); c = 0
    for i in range(n):
        for j in range(n):
            if i == j or not A[i][j]: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[k][i]: c += 1
    return c // 3
def canon(A):
    n = len(A); best = None
    for p in permutations(range(n)):
        flat = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or flat < best: best = flat
    return best
def complement(A):
    n = len(A); B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j: B[i][j] = A[j][i]
    return B
def is_SC(A): return canon(A) == canon(complement(A))
def H_hampaths(A):
    n = len(A); cnt = 0
    for p in permutations(range(n)):
        ok = all(A[p[i]][p[i+1]] for i in range(n-1))
        if ok: cnt += 1
    return cnt

def is_primitive(S):
    g0 = 0
    for v in S: g0 = gcd(g0, v)
    return g0 == 1

# ---------------- THE M2 MAP (exact) ----------------
def adj_m2(S, tau):
    """Return (A, valid). Vertices = sorted speeds. Arc a->b iff snrm((va-vb)tau)>0.
    valid=False if ANY pair has snrm in {0, +/-1/2} (degenerate boundary tie)."""
    S = sorted(set(S)); n = len(S); A = [[0]*n for _ in range(n)]; valid = True
    for a in range(n):
        for b in range(a+1, n):
            s = snrm((S[a] - S[b]) * tau)
            if s == 0 or s == F(1, 2) or s == F(-1, 2):
                valid = False
            elif s > 0:
                A[a][b] = 1
            else:
                A[b][a] = 1
    return A, valid

# ---------------- free sets for n=3,4,5 ----------------
def free_set(n):
    classes = {}
    pairs = list(combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        A = [[0]*n for _ in range(n)]
        for idx, (a, b) in enumerate(pairs):
            if (bits >> idx) & 1: A[a][b] = 1
            else: A[b][a] = 1
        c = canon(A)
        if c not in classes:
            classes[c] = (score_seq(A), count_3cycles(A), H_hampaths(A), is_SC(A))
    return classes

FREE = {n: free_set(n) for n in (3, 4, 5)}
print("FREE SET sizes (A000568 = 2,4,12) and # non-SC classes:")
for n in (3, 4, 5):
    nonsc = sum(1 for v in FREE[n].values() if not v[3])
    print(f"  n={n}: {len(FREE[n])} classes, {nonsc} of them NON-SC")
print()

# ============================================================
# We hunt for ANY single witness: a valid M2 tournament at n in {3,4,5}
# whose iso class is NON-SC. If found -> claim REFUTED.
# ============================================================
WITNESSES = []
def check(S, tau, source):
    A, valid = adj_m2(S, tau)
    if not valid:
        return False
    if not is_SC(A):
        WITNESSES.append((tuple(sorted(set(S))), tau, score_seq(A), count_3cycles(A),
                          H_hampaths(A), source, canon(A)))
        return True
    return False

# ------------------------------------------------------------
# PASS A: optimal lonely tau, all primitive speed sets, big vmax.
# ------------------------------------------------------------
print("="*70)
print("PASS A: optimal lonely tau, primitive speed sets, n in {3,4,5}")
print("="*70)
A_counts = {}
for n_base, vmax in [(3, 40), (4, 26), (5, 18)]:
    tot = 0; sc = 0
    for S in (c for c in combinations(range(1, vmax+1), n_base) if is_primitive(c)):
        _, tau = M(S)
        A, valid = adj_m2(S, tau)
        if not valid: continue
        tot += 1
        if is_SC(A): sc += 1
        else: check(S, tau, f"A:optimal n={n_base}")
    A_counts[n_base] = (tot, sc, vmax)
    print(f"  n={n_base} vmax={vmax}: valid={tot}, SC={sc}, all SC? {sc==tot}", flush=True)

# ------------------------------------------------------------
# PASS B: ALL candidate crossings (THM-524 complete finite set).
# This is the complete set of times where M2 signs can flip on the lonely-
# relevant grid. Exhaustive over the canonical crossing times for each S.
# ------------------------------------------------------------
print("\n" + "="*70)
print("PASS B: ALL candidate crossings (binding-pair + 1/2), n in {3,4,5}")
print("="*70)
B_counts = {}
for n_base, vmax in [(3, 30), (4, 18), (5, 13)]:
    tot = 0; sc = 0; classes = {}
    for S in (c for c in combinations(range(1, vmax+1), n_base) if is_primitive(c)):
        for tau in cand(S):
            A, valid = adj_m2(S, tau)
            if not valid: continue
            tot += 1
            cc = canon(A); classes.setdefault(cc, is_SC(A))
            if is_SC(A): sc += 1
            else: check(S, tau, f"B:cand n={n_base}")
    nonsc_cls = sum(1 for v in classes.values() if not v)
    B_counts[n_base] = (tot, sc, len(classes), nonsc_cls, vmax)
    print(f"  n={n_base} vmax={vmax}: valid={tot}, SC={sc}, all SC? {sc==tot}; "
          f"distinct classes realized={len(classes)} of {len(FREE[n_base])} "
          f"(non-SC classes realized={nonsc_cls})", flush=True)

# ------------------------------------------------------------
# PASS C: OFF-GRID tau. M2's tournament for speed set S is a step function of
# tau, constant on open intervals between crossings of the DIFFERENCE speeds.
# The full set of break-points for M2 is {k/d : d a difference speed, 0<k<d}
# plus {m/(2d)} where snrm hits 1/2. We enumerate ALL these break-points,
# then sample tau STRICTLY BETWEEN consecutive ones (midpoint, exact Fraction)
# to capture EVERY M2 tournament the map produces -- on AND off grid.
# This is an EXHAUSTIVE enumeration of all M2 tournaments for each S.
# ------------------------------------------------------------
def m2_breakpoints(S):
    """All tau in (0,1) where some snrm((va-vb)tau) crosses 0 or +/-1/2.
    snrm(d*tau) has zeros at tau=k/d and hits +/-1/2 at tau=(2k+1)/(2d)."""
    S = sorted(set(S)); bp = set()
    diffs = set()
    for a in range(len(S)):
        for b in range(a+1, len(S)):
            diffs.add(S[b] - S[a])
    for d in diffs:
        # snrm(d*tau)=0 at tau=k/d ; =+/-1/2 at tau=(2k+1)/(2d)
        k = 0
        while F(k, d) < 1:
            bp.add(F(k, d)); k += 1
        k = 0
        while F(2*k+1, 2*d) < 1:
            bp.add(F(2*k+1, 2*d)); k += 1
    bp.add(F(0)); bp.add(F(1))
    return sorted(bp)

def all_m2_tournaments(S):
    """EXHAUSTIVE: every distinct VALID M2 tournament over tau in (0,1).
    Sample the open midpoint between each consecutive pair of break-points."""
    bp = m2_breakpoints(S)
    out = {}
    for i in range(len(bp) - 1):
        lo, hi = bp[i], bp[i+1]
        mid = (lo + hi) / 2  # exact Fraction, strictly interior -> never degenerate
        A, valid = adj_m2(S, mid)
        if not valid:
            continue
        out[canon(A)] = (mid, is_SC(A), score_seq(A), count_3cycles(A))
    return out

print("\n" + "="*70)
print("PASS C: EXHAUSTIVE off-grid -- every M2 tournament over ALL tau in (0,1)")
print("="*70, flush=True)
C_counts = {}
for n_base, vmax in [(3, 16), (4, 14), (5, 12)]:
    realized = {}; tot_classes = 0; nonsc_found = 0
    for S in (c for c in combinations(range(1, vmax+1), n_base) if is_primitive(c)):
        for cc, (mid, sc, ss, c3) in all_m2_tournaments(S).items():
            if cc not in realized:
                realized[cc] = (sc, ss, c3)
            if not sc:
                nonsc_found += 1
                check(S, mid, f"C:offgrid n={n_base}")
    nonsc_cls = sum(1 for v in realized.values() if not v[0])
    C_counts[n_base] = (len(realized), nonsc_cls, vmax)
    print(f"  n={n_base} vmax={vmax}: distinct M2 classes over ALL tau = "
          f"{len(realized)} of {len(FREE[n_base])}; NON-SC classes = {nonsc_cls}; "
          f"non-SC (S,tau) hits = {nonsc_found}")

# ------------------------------------------------------------
# PASS D: sporadic / non-unit / covering / tight / huge-speed sets.
# Speeds not coprime to 14, big spreads, near-AP, Goddyn-Wong-style, etc.
# For each we scan optimal tau AND all M2 tournaments (off-grid exhaustive).
# ------------------------------------------------------------
print("\n" + "="*70)
print("PASS D: sporadic / non-unit / covering / tight / huge-speed sets")
print("="*70)
sporadic = []
# explicit small n in {3,4,5} hand-picked sets
sporadic += [
    (1, 2, 3), (1, 2, 5), (2, 3, 7), (1, 5, 8), (3, 5, 7), (7, 11, 13),
    (1, 6, 8), (2, 5, 11), (1, 13, 14), (1, 7, 14), (2, 7, 14),
    (1, 2, 3, 4), (1, 2, 3, 5), (1, 3, 5, 7), (2, 3, 5, 7), (1, 5, 7, 11),
    (1, 2, 4, 8), (1, 6, 10, 15), (1, 7, 13, 14), (2, 5, 9, 13),
    (1, 2, 3, 4, 5), (1, 3, 5, 7, 9), (2, 3, 5, 7, 11), (1, 2, 4, 8, 16),
    (1, 5, 7, 11, 13), (1, 6, 8, 13, 14), (3, 5, 7, 11, 13), (1, 2, 3, 4, 14),
    (1, 7, 14, 21, 28),  # all multiples-of-7 flavored
    (5, 9, 11),          # all coprime-to-14 units
    (1, 3, 9),           # geometric
    (1, 4, 9, 16, 25),   # squares
]
tot = 0; sc = 0; nonsc = 0
for S in sporadic:
    Sp = tuple(sorted(set(S)))
    if len(Sp) not in (3, 4, 5):
        continue
    # optimal tau
    _, tau = M(Sp)
    A, valid = adj_m2(Sp, tau)
    if valid:
        tot += 1
        if is_SC(A): sc += 1
        else:
            nonsc += 1; check(Sp, tau, "D:sporadic-optimal")
    # full off-grid exhaustion
    for cc, (mid, scv, ss, c3) in all_m2_tournaments(Sp).items():
        if not scv:
            nonsc += 1; check(Sp, mid, "D:sporadic-offgrid")
print(f"  sporadic sets scanned: {len(sporadic)}; optimal-tau valid={tot}, SC={sc}; "
      f"non-SC hits (incl off-grid)={nonsc}")

# ------------------------------------------------------------
# PASS E: structural sanity -- the reflection identity.
# Prove computationally that T(tau)^op == T(1-tau) for M2 at n in {3,4,5}
# (this is WHY the realized SET is complement-closed). Confirm it does NOT
# by itself force each T to be SC, by exhibiting (if any) tau with T(tau)
# != T(1-tau) up to iso among the valid configs.
# ------------------------------------------------------------
print("\n" + "="*70)
print("PASS E: reflection identity T(tau)^op vs T(1-tau), and SC mechanism")
print("="*70)
refl_ok = True; refl_checked = 0
self_recip = 0  # configs where T(tau) iso T(1-tau) -> consistent with SC
for n_base, vmax in [(3, 20), (4, 14), (5, 11)]:
    for S in (c for c in combinations(range(1, vmax+1), n_base) if is_primitive(c)):
        for tau in cand(S):
            A, va = adj_m2(S, tau)
            B, vb = adj_m2(S, 1 - tau)
            if not (va and vb):
                continue
            refl_checked += 1
            # reflection identity: B should equal complement(A) exactly (same vertex order)
            if canon(B) != canon(complement(A)):
                refl_ok = False
            if canon(A) == canon(B):
                self_recip += 1
print(f"  reflection identity T(1-tau) == T(tau)^op holds on all {refl_checked} "
      f"checked valid configs? {refl_ok}")
print(f"  (configs where T(tau) ~ T(1-tau) as iso classes: {self_recip})")

# ------------------------------------------------------------
# VERDICT
# ------------------------------------------------------------
print("\n" + "="*70)
print("VERDICT")
print("="*70)
if WITNESSES:
    print(f"REFUTED. Found {len(WITNESSES)} non-SC M2 tournament witness(es) at n<=5.")
    for w in WITNESSES[:10]:
        S, tau, ss, c3, H, src, cc = w
        print(f"  S={S} tau={tau} score={ss} c3={c3} H={H}  [{src}]")
else:
    print("NO non-SC M2 tournament found at n=3,4,5 across all passes.")
    print("Claim CONFIRMED within the search bound:")
    print("  PASS A optimal-tau bounds:", A_counts)
    print("  PASS B candidate-crossing bounds:", B_counts)
    print("  PASS C EXHAUSTIVE off-grid (ALL tau) bounds:", C_counts)
    print("  -> PASS C is exhaustive over ALL M2 tournaments for each S in range,")
    print("     so within the speed-range bound NO non-SC class is reachable.")

print("\nDONE.")
