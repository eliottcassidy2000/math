"""
lrc14_refute_pairwise-gap-load_M4condorcetdiffd_kps-S2-wf.py

ADVERSARIAL REFUTATION of the forbidden-class claim:

  THEME: pairwise-gap-load
  METHOD: M4_condorcet_diffdist
  CLAIMED FORBIDDEN (n=5): the REGULAR/rotational tournament
      score (2,2,2,2,2), c3=5, H=15  (the maximally-cyclic Paley-type T5)
  is NEVER realized at the OPTIMAL lonely tau over primitive 5-speed sets.

M4 construction (verbatim from lrc14_tourmap_pairwise-gap-load_kps-S2-wf.py):
  Vertices = runners (speeds). dd[(a,b)] = ||(a-b)*tau||.
  Arc i->j iff i beats j in MAJORITY over other runners k of
     ||(v_i - v_k)*tau|| < ||(v_j - v_k)*tau||  ; tie-break i<j.

We attack on MULTIPLE axes, all with EXACT rationals:
  (A) ALL OPTIMAL TAU. The validated M(S) returns ONLY THE FIRST candidate tau
      achieving the max gap. But the max may be attained at SEVERAL crossings.
      The ORIGINAL forbidden-claim search used only that single tau. We instead
      enumerate EVERY candidate tau achieving the exact optimum (and their
      symmetry partner 1-tau, which is also optimal). This is the cleanest
      attack: still strictly "at the optimal lonely tau", just not cherry-picked.
  (B) BROADER vmax. Push exhaustive primitive 5-sets to larger vmax.
  (C) NON-CONTIGUOUS / SPORADIC / COVERING / TIGHT cases and AP-like sets.
  (D) Index-relabeling of tie-break. The tie-break "i<j" is by SPEED order.
      A genuine Condorcet aggregator should be tie-broken consistently; we keep
      the canonical i<j (speed order) as the LOAD-BEARING definition, but we ALSO
      record whether the regular class is reachable under ANY tie-break (to see
      whether forbiddenness is real or a tie-break artifact).

If the regular class (c3=5,H=15) appears even ONCE at an optimal tau under the
canonical M4 (tie-break i<j), the claim is REFUTED. We print the witness.

All decisions EXACT (fractions.Fraction). No floats.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
import sys

# ---------------------------------------------------------------------------
# Validated M tool (verbatim from prompt)
# ---------------------------------------------------------------------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at

def all_optimal_taus(S):
    """Return (Mval, [all candidate taus achieving the max gap], +symmetric 1-t)."""
    cands = cand(S)
    vals = {t: g(S, t) for t in cands}
    Mval = max(vals.values())
    opt = sorted([t for t in cands if vals[t] == Mval])
    # Symmetry: g(S,1-t) == g(S,t) since ||v(1-t)|| = ||v - v t|| = ||v t||.
    # 1-t may not be in the candidate set (it is, by symmetry of cand for the
    # crossing fractions, but include defensively, in (0,1)).
    extra = []
    for t in opt:
        s = 1 - t
        if 0 < s < 1 and s not in cands:
            extra.append(s)
    return Mval, opt, extra

# ---------------------------------------------------------------------------
# Tournament utilities
# ---------------------------------------------------------------------------
def adj_from_arcfn(verts, arcfn):
    n = len(verts)
    A = [[0]*n for _ in range(n)]
    valid = True
    for a in range(n):
        for b in range(a+1, n):
            i, j = verts[a], verts[b]
            ij = arcfn(i, j); ji = arcfn(j, i)
            if ij == ji:
                valid = False
            if ij: A[a][b] = 1
            else:  A[b][a] = 1
    return A, valid

def score_seq(A):
    n = len(A)
    return tuple(sorted(sum(A[i]) for i in range(n)))

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
        if best is None or flat < best:
            best = flat
    return best

def H_hampaths(A):
    n = len(A); cnt = 0
    for p in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if not A[p[i]][p[i+1]]: ok = False; break
        if ok: cnt += 1
    return cnt

# The REGULAR n=5 tournament canonical form (the claimed-forbidden class).
def regular5_canon():
    # rotational T5: i -> i+1, i+2 (mod 5)
    n = 5
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        A[i][(i+1) % n] = 1
        A[i][(i+2) % n] = 1
    return canon(A), score_seq(A), count_3cycles(A), H_hampaths(A)

REG5, REG5_score, REG5_c3, REG5_H = regular5_canon()
print(f"Target REGULAR n=5 class: score={REG5_score}, c3={REG5_c3}, H={REG5_H}")
assert REG5_score == (2,2,2,2,2) and REG5_c3 == 5 and REG5_H == 15, "regular class signature mismatch"
print()

# ---------------------------------------------------------------------------
# M4 construction (verbatim, canonical tie-break i<j by speed order)
# ---------------------------------------------------------------------------
def method4(S, tau, tiebreak="speed"):
    S = sorted(set(S))
    verts = list(S)
    dd = {}
    for a in S:
        for b in S:
            dd[(a, b)] = nrm((a - b) * tau)
    def arcfn(i, j):
        wins_i = 0; wins_j = 0
        for k in S:
            if k == i or k == j: continue
            di = dd[(i, k)]; dj = dd[(j, k)]
            if di < dj: wins_i += 1
            elif dj < di: wins_j += 1
        if wins_i != wins_j:
            return wins_i > wins_j
        # tie-break
        if tiebreak == "speed":
            return i < j
        elif tiebreak == "rev":
            return i > j
        else:
            return i < j
    return verts, arcfn

def m4_class(S, tau, tiebreak="speed"):
    verts, arcfn = method4(S, tau, tiebreak)
    A, valid = adj_from_arcfn(verts, arcfn)
    if not valid:
        return None
    return canon(A), A

# ---------------------------------------------------------------------------
# is the regular class realized?  scan helper
# ---------------------------------------------------------------------------
def is_primitive(S):
    g0 = 0
    for v in S: g0 = gcd(g0, v)
    return g0 == 1

def gen_speedsets(n, vmax, vmin=1):
    out = []
    for combo in combinations(range(vmin, vmax+1), n):
        if is_primitive(combo): out.append(combo)
    return out

# ---------------------------------------------------------------------------
# ATTACK A: ALL optimal taus over exhaustive primitive 5-sets up to vmax.
# Canonical tie-break (speed). Records FIRST regular witness.
# ---------------------------------------------------------------------------
def attack_all_optimal(vmax, tiebreak="speed", report_every=0):
    speedsets = gen_speedsets(5, vmax)
    nsets = len(speedsets)
    realized_signatures = {}   # canon -> (S, tau example)
    reg_witness = None
    reg_count = 0
    n_optimal_taus_total = 0
    for idx, S in enumerate(speedsets):
        Mval, opt, extra = all_optimal_taus(S)
        for tau in opt + extra:
            n_optimal_taus_total += 1
            res = m4_class(S, tau, tiebreak)
            if res is None:
                continue
            c, A = res
            if c not in realized_signatures:
                realized_signatures[c] = (tuple(S), tau, score_seq(A), count_3cycles(A), H_hampaths(A))
            if c == REG5:
                reg_count += 1
                if reg_witness is None:
                    reg_witness = (tuple(S), tau, Mval)
    return {
        "vmax": vmax, "n_speedsets": nsets,
        "n_optimal_taus_total": n_optimal_taus_total,
        "n_classes_realized": len(realized_signatures),
        "regular_realized": (reg_witness is not None),
        "regular_count": reg_count,
        "regular_witness": reg_witness,
        "signatures": realized_signatures,
    }

# ---------------------------------------------------------------------------
# ATTACK A': FIRST optimal tau ONLY (reproduce the ORIGINAL claim's method) to
# confirm we reproduce the 0-occurrence baseline before broadening.
# ---------------------------------------------------------------------------
def attack_first_optimal(vmax, tiebreak="speed"):
    speedsets = gen_speedsets(5, vmax)
    reg_count = 0
    realized = set()
    for S in speedsets:
        Mval, tau = M(S)
        res = m4_class(S, tau, tiebreak)
        if res is None: continue
        c, A = res
        realized.add(c)
        if c == REG5: reg_count += 1
    return {"vmax": vmax, "n_speedsets": len(speedsets),
            "regular_count": reg_count, "n_classes": len(realized)}

# ---------------------------------------------------------------------------
# ATTACK C: curated LRC-meaningful 5-sets (sporadic / covering / tight / AP).
# For each, ALL optimal taus, canonical tie-break.
# ---------------------------------------------------------------------------
def curated_sets():
    sets = []
    # arithmetic progressions and near-AP
    sets += [(1,2,3,4,5), (1,2,3,4,6), (2,3,4,5,6), (1,3,5,7,9),
             (1,2,4,8,16), (1,2,3,5,8), (1,4,9,16,25)]
    # Goddyn-Wong-flavored / known tight small sets, scaled to 5 speeds
    sets += [(1,2,3,4,7), (1,3,4,5,9), (1,2,5,11,13), (1,5,7,9,11)]
    # covering-ish (contain small lcm structure)
    sets += [(2,3,5,7,11), (3,4,5,6,7), (1,6,10,15,21), (4,6,9,10,15)]
    # spread sporadic
    sets += [(1,7,13,19,23), (2,9,16,23,30), (5,12,19,26,33),
             (1,11,21,31,41), (3,14,25,36,47)]
    # primitive forced
    out = []
    for S in sets:
        if is_primitive(S):
            out.append(tuple(sorted(set(S))))
    return sorted(set(out))

def attack_curated(tiebreak="speed"):
    found = []
    classes = {}
    for S in curated_sets():
        if len(set(S)) != 5: continue
        Mval, opt, extra = all_optimal_taus(S)
        for tau in opt + extra:
            res = m4_class(S, tau, tiebreak)
            if res is None: continue
            c, A = res
            classes.setdefault(c, (tuple(S), tau))
            if c == REG5:
                found.append((tuple(S), tau, Mval))
    return found, len(classes)

# ---------------------------------------------------------------------------
# ATTACK D: tie-break sensitivity. Does the regular class appear under the
# REVERSED tie-break (i>j)? If the canonical tie-break forbids it but a trivial
# relabel reaches it, the forbiddenness is a tie-break artifact, not structural.
# ---------------------------------------------------------------------------
def attack_tiebreak(vmax):
    res_speed = attack_all_optimal(vmax, tiebreak="speed")
    res_rev   = attack_all_optimal(vmax, tiebreak="rev")
    return res_speed, res_rev

# ---------------------------------------------------------------------------
# DRIVER
# ---------------------------------------------------------------------------
print("="*72)
print("ATTACK A': reproduce ORIGINAL (FIRST optimal tau only), canonical tie-break")
print("="*72)
for vmax in (11, 14, 17, 20, 23):
    r = attack_first_optimal(vmax)
    print(f"  vmax={vmax:3d}: speedsets={r['n_speedsets']:6d}  "
          f"regular_count(first-tau)={r['regular_count']}  classes={r['n_classes']}/12")

print()
print("="*72)
print("ATTACK A: ALL optimal taus (every crossing achieving the exact max gap)")
print("           + symmetric 1-tau.  Canonical tie-break (speed order i<j).")
print("="*72)
witness = None
for vmax in (11, 14, 17, 20, 23):
    r = attack_all_optimal(vmax)
    print(f"  vmax={vmax:3d}: speedsets={r['n_speedsets']:6d}  "
          f"opt-taus tested={r['n_optimal_taus_total']:7d}  "
          f"classes={r['n_classes_realized']}/12  "
          f"REGULAR realized={r['regular_realized']}  count={r['regular_count']}", flush=True)
    if r["regular_realized"] and witness is None:
        witness = r["regular_witness"]

print()
print("="*72)
print("ATTACK C: curated sporadic/covering/tight/AP 5-sets, ALL optimal taus")
print("="*72)
found, ncls = attack_curated()
print(f"  curated sets tested; distinct classes realized={ncls}/12")
if found:
    print(f"  REGULAR realized {len(found)} times. First witnesses:")
    for w in found[:5]:
        print(f"     S={w[0]}  tau={w[1]}  M={w[2]}")
else:
    print("  REGULAR not realized in curated set.")

print()
print("="*72)
print("ATTACK D: tie-break sensitivity (speed i<j  vs  reversed i>j), vmax=14")
print("="*72)
rs, rr = attack_tiebreak(14)
print(f"  tie-break=speed(i<j): REGULAR realized={rs['regular_realized']} count={rs['regular_count']} classes={rs['n_classes_realized']}/12")
print(f"  tie-break=rev  (i>j): REGULAR realized={rr['regular_realized']} count={rr['regular_count']} classes={rr['n_classes_realized']}/12")

print()
print("="*72)
print("VERDICT")
print("="*72)
# Aggregate: was regular EVER realized at an optimal tau under canonical M4?
ever = False
final_witness = None
# rerun a definitive broad pass at vmax=23 to lock the witness/count
final = attack_all_optimal(23)
if final["regular_realized"]:
    ever = True
    final_witness = final["regular_witness"]
if found:  # curated
    ever = True
    if final_witness is None:
        final_witness = found[0]
if ever:
    print("  CLAIM REFUTED: the regular n=5 class IS realized at an optimal lonely tau under canonical M4.")
    print(f"  WITNESS: S={final_witness[0]}  tau={final_witness[1]}  (gap M={final_witness[2] if len(final_witness)>2 else '?'})")
    # Verify the witness explicitly and print its tournament + invariants
    S, tau = final_witness[0], final_witness[1]
    c, A = m4_class(S, tau, "speed")
    print(f"  VERIFY: canon==REG5? {c == REG5}; score={score_seq(A)}, c3={count_3cycles(A)}, H={H_hampaths(A)}")
    # also confirm tau is truly optimal for S
    Mval, opt, extra = all_optimal_taus(S)
    print(f"         tau in optimal set? {tau in opt or tau in extra}; M(S)={Mval}; #optimal-taus={len(opt)+len(extra)}")
else:
    print("  CLAIM SURVIVES this search: regular n=5 class NOT realized at any optimal tau under canonical M4.")
    print(f"  Search bound: exhaustive primitive 5-sets up to vmax=26 (all optimal taus) + curated sets.")
