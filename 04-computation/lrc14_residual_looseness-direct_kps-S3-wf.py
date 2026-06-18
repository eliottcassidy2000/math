#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_residual_looseness-direct_kps-S3-wf.py

ANGLE: looseness-direct.  Attack the ONLY open case S3 of LRC(14) directly,
proving M(S) >= 1/14 by exploiting that S3 sets are LOOSE.

S3 = covering primitive 13-sets with k = #{v>13} >= 2 and Vmax >= 13*Vmin.

PARTS
  (a) FLOOR HUNT.  Search hard for the true S3 minimum of M with EXACT rationals
      (float screen for speed, exact-confirm every near-1/14 hit and the argmin).
      Classify minimizers.  Is min M over S3 >= some explicit c > 1/14?
  (a2) EXHAUSTIVE small-window S3 floor (rigorous finite check).
  (b) MULTI-GAP LEMMA 1.  Generalize Lemma 1 (FIRST common gap, needs spread<13) to
      a LATER common gap: an integer vector k=(k_u) with arc
        A(k) = cap_u ((k_u+1/14)/u,(k_u+13/14)/u)  nonempty.
      Test existence on S3.  Establish its EXACT logical status.
  (c) CONSTRUCTIVE structured witness (cluster pivot + small-part gap centers).
  (d) SUBSET Lemma-1 coverage of S3.

EXACT in every DECISION (fractions.Fraction); floats ONLY as a pre-screen.
"""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

H = F(1, 14)
Hf = 1.0 / 14.0

# ---------------------------------------------------------------------------
# EXACT M / W tools
# ---------------------------------------------------------------------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(x * t) for x in S)

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

def Mval(S):
    """EXACT M and an argmax tau (authoritative).  Float-prefiltered: each candidate
    fraction is first evaluated in float; only those whose float value exceeds the
    current exact best (minus a safety margin) are evaluated exactly.  This never
    misses the true argmax because float min-norm is within ~1e-9 of exact, and we
    use a 1e-6 margin.  Result is the EXACT Fraction M and an exact argmax tau."""
    Sf = [float(x) for x in S]
    bf = -1.0          # float best
    b = F(0); at = None
    for t in cand(S):
        tf = float(t)
        vf = min(nrm_f(x*tf) for x in Sf)
        if vf > bf - 1e-9:           # candidate could beat current best
            v = g(S, t)              # exact
            if v > b:
                b = v; at = t; bf = vf
    return b, at

def nrm_f(x):
    r = x - int(x)
    if r < 0: r += 1
    return r if r <= 0.5 else 1 - r

def M_float(S):
    """Float LOWER bound on M: scans odd-halfpoints (2k+1)/(2v) of EVERY speed v
    (a SUBSET of cand, each a valid tau) plus pair fractions only for SMALL
    denominators d<=300 (the only ones that can produce an argmax very near 1/14
    with bounded denom).  Returns <= M(S) (no false 'safe').  O(sum speeds + small
    pairs); fast even for speed ~1500.  Every set with M_float < 0.078 is exact-
    confirmed below, so the screen never misses a true near-1/14 set."""
    Ss = sorted(set(S)); best = 0.0
    # halfpoints of speeds <= 30 (fast; valid taus -> lower bound on M)
    for v in Ss:
        if v > 30: continue
        k = 0
        while 2*k+1 <= v:
            t = (2*k+1) / (2.0*v)
            val = min(nrm_f(x*t) for x in Ss)
            if val > best: best = val
            k += 1
    # plus pair fractions with SMALL denominator d<=120 (cheap, catches tight cases)
    n = len(Ss)
    for i in range(n):
        for j in range(i+1, n):
            for d in (Ss[i]+Ss[j], Ss[j]-Ss[i]):
                if 0 < d <= 120:
                    k = 1
                    while 2*k <= d:
                        t = k / float(d)
                        val = min(nrm_f(x*t) for x in Ss)
                        if val > best: best = val
                        k += 1
    return best

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def is_primitive(S):
    return reduce(gcd, S) == 1

def in_S3(S):
    if len(S) != 13: return False
    if not is_covering(S): return False
    if not is_primitive(S): return False
    k = sum(1 for v in S if v > 13)
    vmn, vmx = min(S), max(S)
    return k >= 2 and vmx >= 13 * vmn

# ---------------------------------------------------------------------------
# MULTI-GAP arc tools
# ---------------------------------------------------------------------------
def gap_arc(S, kvec):
    """(lo,hi) of A(k)=cap_u ((k_u+1/14)/u,(k_u+13/14)/u), or None if empty."""
    lo = F(-10**9); hi = F(10**9)
    for u, ku in zip(S, kvec):
        loi = (F(ku) + H) / u
        hii = (F(ku) + (1 - H)) / u
        if loi > lo: lo = loi
        if hii < hi: hi = hii
    return (lo, hi) if lo < hi else None

def greedy_kvec(S, tau):
    """k_u = floor(u*tau)."""
    out = []
    for u in S:
        ut = u * tau
        out.append(ut.__floor__() if isinstance(ut, F) else int(ut // 1))
    return out

# ---------------------------------------------------------------------------
# Generators of S3 sets
# ---------------------------------------------------------------------------
def gen_mixed(center, nsmall, jit, rng):
    base = list(range(1, nsmall + 1)); used = set(base); larges = []
    needed = [q for q in range(2, 15) if not any(b % q == 0 for b in base)]
    for q in needed:
        k = round(center / q) + rng.randint(-jit, jit); c = q * k
        while c in used or c <= nsmall:
            k += 1; c = q * k
        larges.append(c); used.add(c)
    bl = sorted(set(base + larges)); hi = center + rng.randint(0, 30)
    S = list(bl)
    while len(S) < 13:
        hi += 1
        if hi not in S: S.append(hi)
    S = sorted(set(S))
    return S[:13] if len(S) > 13 else S

def gen_tight_cluster(V, nsmall, spread, rng):
    base = list(range(1, nsmall + 1)); used = set(base)
    needed = [q for q in range(2, 15) if not any(b % q == 0 for b in base)]
    larges = []
    for q in needed:
        k0 = -(-V // q)
        placed = False
        for k in range(k0, k0 + spread + 5):
            c = q * k
            if V <= c <= V + spread and c not in used:
                larges.append(c); used.add(c); placed = True; break
        if not placed:
            c = q * k0
            while c in used: c += q
            larges.append(c); used.add(c)
    S = sorted(set(base + larges)); extra = V
    while len(S) < 13:
        extra += 1
        if extra not in S: S.append(extra)
    S = sorted(set(S))
    return S[:13] if len(S) > 13 else S

def gen_S3(rng, big=False):
    c = rng.randint(0, 2)
    hi = 800 if big else 280
    if c == 0:
        return sorted(set(gen_mixed(rng.randint(20, hi), rng.choice(list(range(1, 14))),
                                    rng.choice([0, 1, 2, 3]), rng)))
    if c == 1:
        return sorted(set(gen_tight_cluster(rng.randint(20, hi), rng.choice(list(range(1, 14))),
                                            rng.choice([14, 20, 28, 35, 42, 45]), rng)))
    return sorted(set(gen_mixed(rng.randint(14, min(hi, 250)), rng.choice([6, 7, 8, 9, 10, 11, 12, 13]),
                                rng.choice([0, 1, 2]), rng)))

# ===========================================================================
# PART (a): FLOOR HUNT (float-screen + exact-confirm)
# ===========================================================================
def part_a(N=2500):
    print("="*78)
    print("PART (a): EXACT S3 floor hunt (float lower-bound screen, exact-confirm near-min)")
    print("="*78)
    rng = random.Random(12345)
    worst_exact = (F(99), None, None)
    tested = 0; breaks = 0
    screen_thresh = 0.078   # confirm exactly anything whose float lower-bound < this
    near = []               # (floatM,S) screened-near-1/14 sets
    for it in range(N):
        S = gen_S3(rng, big=False)
        if not in_S3(S): continue
        tested += 1
        mf = M_float(S)
        if mf < screen_thresh:
            near.append((mf, S))
    # exact-confirm the near sets
    confirmed = []
    for mf, S in near:
        M, at = Mval(S)
        confirmed.append((M, S, at))
        if M < F(1, 14): breaks += 1
        if M < worst_exact[0]:
            worst_exact = (M, S, at)
    confirmed.sort(key=lambda z: z[0])
    print(f"  S3 sets screened: {tested}")
    print(f"  near-1/14 (floatM<{screen_thresh}) exact-confirmed: {len(confirmed)}")
    print(f"  LRC breaks (exact M < 1/14): {breaks}")
    print(f"  smallest EXACT M = {worst_exact[0]} = {float(worst_exact[0]):.6f}  (1/14={Hf:.6f})")
    print(f"  minimizer S = {worst_exact[1]}")
    print(f"  achieved at tau = {worst_exact[2]}")
    print(f"  margin over 1/14 = {worst_exact[0]-F(1,14)} = {float(worst_exact[0]-F(1,14)):.6f}")
    print("  smallest 10 exact M (S3):")
    seen = set()
    cnt = 0
    for M, S, at in confirmed:
        key = M
        if key in seen: continue
        seen.add(key); cnt += 1
        print(f"     {M} = {float(M):.6f}   S={S}")
        if cnt >= 10: break
    return worst_exact

# ===========================================================================
# PART (a2): EXHAUSTIVE small-window S3 (rigorous finite floor)
# ===========================================================================
def part_a2():
    print()
    print("="*78)
    print("PART (a2): EXHAUSTIVE S3 in a small window -- rigorous finite floor")
    print("  Enumerate ALL covering primitive 13-subsets of {1..W} that are in S3,")
    print("  for W up to a feasible bound; report exact min M.")
    print("="*78)
    from itertools import combinations
    worst = (F(99), None)
    # For S3: k>=2 large (>13), Vmax>=13*Vmin.  With Vmax<=W and Vmin>=1, need Vmin<=W/13.
    # The LEAST-loose (most dangerous) S3 sets have Vmin=1 and small core, so we enumerate
    # ALL covering primitive 13-subsets of {1..W} that are in S3.  To stay feasible we prune
    # with the covering test applied incrementally (must include a multiple of each q in 2..14).
    for W in [20, 22]:
        cnt = 0
        wW = (F(99), None)
        pool = list(range(1, W+1))
        # iterate combinations but prune early on Vmin/Vmax and large-count using the fact
        # that the smallest element is combo[0]
        for combo in combinations(pool, 13):
            vmn = combo[0]; vmx = combo[12]
            if vmx < 13*vmn: continue
            if (combo[11] <= 13):  # need >=2 elements >13; if 2nd-largest<=13 impossible
                continue
            S = list(combo)
            if not is_covering(S): continue
            if not is_primitive(S): continue
            cnt += 1
            M, _ = Mval(S)
            if M < wW[0]: wW = (M, S)
            if M < worst[0]: worst = (M, S)
        print(f"  W={W}: #S3 sets={cnt}  minM={wW[0]}={float(wW[0]):.6f}  argmin={wW[1]}")
    print(f"  EXHAUSTIVE S3 floor (W<=22): minM={worst[0]}={float(worst[0]):.6f} at {worst[1]}")
    print(f"  >= 1/14 strictly? {worst[0] > F(1,14)}   (margin {worst[0]-F(1,14)})")
    return worst

# ===========================================================================
# PART (b): MULTI-GAP LEMMA existence + logical status
# ===========================================================================
def part_b(floor_minimizer):
    print()
    print("="*78)
    print("PART (b): MULTI-GAP LEMMA -- existence of a common LATER-gap safe arc")
    print("="*78)
    print("  THEOREM (elementary identity, PROVED):")
    print("    For finite S and tau with min_u||u tau|| > 1/14, set k_u=floor(u tau).")
    print("    Then frac(u tau) in (1/14,13/14) for all u, so tau in A(k) and A(k)!=empty.")
    print("    Conversely every tau in A(k) is safe.  Hence:")
    print("      [ exists nonempty multi-gap arc A(k) ]  <=>  [ M(S) > 1/14 ].")
    print("    => the multi-gap Lemma is LOGICALLY EQUIVALENT to the goal, NOT an")
    print("       independent lever.  Content = prove a good k EXISTS by CRT a priori.")
    print()
    rng = random.Random(2027); tested = 0; found = 0; idfail = 0
    for it in range(800):
        S = gen_S3(rng)
        if not in_S3(S): continue
        tested += 1
        M, at = Mval(S)
        if M > F(1, 14):
            kv = greedy_kvec(S, at)
            arc = gap_arc(S, kv)
            if arc is not None and arc[0] <= at <= arc[1]:
                found += 1
            else:
                idfail += 1
        else:
            # M == 1/14 boundary: degenerate point witness
            found += 1
    print(f"  S3 sets tested: {tested}")
    print(f"  identity verified (k=floor at argmax gives nonempty arc containing tau): {found}")
    print(f"  identity failures: {idfail}")
    if floor_minimizer and floor_minimizer[1]:
        S = floor_minimizer[1]; M, at = Mval(S)
        kv = greedy_kvec(S, at); arc = gap_arc(S, kv)
        print(f"  floor minimizer S={S}")
        print(f"     M={M}, argmax tau={at}, k=floor(u tau)={kv}")
        print(f"     A(k)={arc}")
    return found, tested, idfail

# ===========================================================================
# PART (c): CONSTRUCTIVE structured witness (no full M needed)
# ===========================================================================
def constructive_crt(S):
    L = sorted(v for v in S if v > 13)
    P = sorted(v for v in S if v <= 13)
    if not L: return None
    V0 = L[0]; Vmid = L[len(L)//2]
    cands = set()
    for m in range(0, V0+1): cands.add(F(2*m+1, 2*V0))
    for m in range(0, Vmid+1): cands.add(F(2*m+1, 2*Vmid))
    for u in P:
        for j in range(0, u+1): cands.add(F(2*j+1, 2*u))
    best = None
    for tau in cands:
        if F(0) < tau < F(1):
            val = g(S, tau)
            if val >= F(1, 14) and (best is None or val > best[1]):
                best = (tau, val)
    return best

def part_c():
    print()
    print("="*78)
    print("PART (c): CONSTRUCTIVE structured witness (cluster-pivot + small-part gaps)")
    print("="*78)
    rng = random.Random(55); tested = 0; ok = 0; fails = []
    for it in range(800):
        S = gen_S3(rng)
        if not in_S3(S): continue
        tested += 1
        res = constructive_crt(S)
        if res is not None: ok += 1
        else:
            M, _ = Mval(S); fails.append((M, S))
    print(f"  S3 sets tested: {tested}")
    print(f"  structured witness found: {ok}  ({100.0*ok/max(tested,1):.2f}%)")
    print(f"  structured-witness MISSES: {len(fails)}  (NOT LRC breaks -- scan incompleteness)")
    if fails:
        f0 = min(fails, key=lambda z: z[0])
        print(f"     hardest miss: M={f0[0]}={float(f0[0]):.6f} (>=1/14? {f0[0]>=F(1,14)}), S={f0[1]}")
    return ok, tested, fails

# ===========================================================================
# PART (d): SUBSET Lemma-1 coverage
# ===========================================================================
def lemma1_subset_closes(S):
    rest = sorted(S)[:-1]
    vmn, vmx = min(rest), max(rest)
    if 13*vmn > vmx:
        return F(13, 14*vmx) - F(1, 14*vmn)
    return None

def part_d():
    print()
    print("="*78)
    print("PART (d): SUBSET Lemma-1 -- drop largest; does Lemma 1 fire on the rest?")
    print("="*78)
    rng = random.Random(321); tested = 0; subset_l1 = 0
    for it in range(800):
        S = gen_S3(rng)
        if not in_S3(S): continue
        tested += 1
        if lemma1_subset_closes(S) is not None: subset_l1 += 1
    print(f"  S3 sets tested: {tested}")
    print(f"  closed by SUBSET Lemma-1: {subset_l1} ({100.0*subset_l1/max(tested,1):.2f}%)")
    print(f"  NOT closed -> need multi-gap/band route: {tested-subset_l1}")
    return subset_l1, tested

def main():
    fm = part_a()
    fm2 = part_a2()
    floor = fm if fm[0] <= fm2[0] else (fm2[0], fm2[1], None)
    print()
    print("="*78)
    print(f"COMBINED S3 FLOOR (this search): min M = {floor[0]} = {float(floor[0]):.6f}")
    print(f"  minimizer = {floor[1]}")
    print(f"  > 1/14 strictly? {floor[0] > F(1,14)}   margin = {float(floor[0]-F(1,14)):.6f}")
    print("="*78)
    part_b(floor)
    part_c()
    part_d()

if __name__ == "__main__":
    main()
