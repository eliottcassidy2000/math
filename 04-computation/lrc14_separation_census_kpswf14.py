#!/usr/bin/env python3
"""
kind-pasteur-2026-06-22 (kpswf14) THREAD 2 DECISIVE TEST.

Question: is there a MAGNITUDE-AWARE tournament whose ISO CLASS (here fingerprinted by the
near-complete invariant H = #Hamiltonian paths, plus score-seq + c3) SEPARATES the TIGHT
LRC(14) rows {AP, GW} from the LOOSE q=14 rows -- in particular the apex-twins AP / 12->26 /
12->96 which are identical under any single periodic winding tournament?

We test, over a BANK of tight + loose rows:
  (a) winding T(3/41)                        [periodic, magnitude-blind on apex-twins]
  (c) divisibility-cutoff  i->j iff s_i%s_j < s_j/2   [NON-periodic, magnitude-aware]
  (d) ADDITIVE-GAP tournament i->j iff (s_i - s_j) mod 14 in (0,7) but ENRICHED with magnitude:
      tie-break apex residue-ties by the REAL difference sign  [hybrid]
  (e) RECIPROCAL-WINDING T at t = 1/s for each speed: stacked "self-clock"
We report, per candidate, whether the SET of H-values on tight rows is DISJOINT from the set
on loose rows  ==> clean separation of the class.

Exact arithmetic (Fraction). H via bitmask DP. n=13.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce

# ---------------- iso fingerprint ----------------
def H_paths(A):
    n = len(A); size = 1 << n
    dp = [[0]*n for _ in range(size)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c == 0: continue
            Av = A[v]
            for w in range(n):
                if mask & (1 << w): continue
                if Av[w]: dp[mask | (1 << w)][w] += c
    return sum(dp[size - 1])

def scores(A):
    n = len(A); return tuple(sorted(sum(A[i]) for i in range(n)))
def c3(A):
    n = len(A); cnt = 0
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]: cnt += 1
        if A[a][c] and A[c][b] and A[b][a]: cnt += 1
    return cnt
def is_tournament(A):
    n = len(A)
    return all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(i+1, n))
def fingerprint(A):
    return (scores(A), c3(A), H_paths(A))

# ---------------- exact LRC measure (tight iff M==1/14) ----------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t): return min(nrm(v * t) for v in S)
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
def M_exact(S):
    b = F(0)
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v
    return b
TIGHT = F(1, 14)

# ---------------- candidate tournaments ----------------
def winding(S, t):
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            rel = (F(S[i]-S[j]) * t) % 1
            if 0 < rel < F(1, 2): A[i][j] = 1
            elif rel > F(1, 2): A[i][j] = 0
            else: A[i][j] = 1 if S[i] < S[j] else 0
    return A

def divis_cutoff(S):
    k = len(S)
    raw = [[2*(S[i] % S[j]) < S[j] for j in range(k)] for i in range(k)]
    A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            a, b = raw[i][j], raw[j][i]
            if a and not b: A[i][j] = 1
            elif b and not a: A[j][i] = 1
            else:  # tie -> smaller speed beats
                if S[i] < S[j]: A[i][j] = 1
                else: A[j][i] = 1
    return A

def divis_cutoff_floor(S):
    """variant: i->j iff floor(s_i / s_j) is ODD  (Beatty/continued-fraction flavored)."""
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            qi = S[i] // S[j]; qj = S[j] // S[i]
            ai = (qi % 2 == 1); aj = (qj % 2 == 1)
            # orient by whichever quotient is odd; tie-break by speed
            if ai and not aj: A[i][j] = 1
            elif aj and not ai: A[j][i] = 1
            else:
                if S[i] < S[j]: A[i][j] = 1
                else: A[j][i] = 1
    return A

def additive_hybrid(S):
    """apex residue order, but BREAK the antipodal/equal residue ties by REAL difference sign.
    i->j iff [ (s_i-s_j) mod 14 in (0,7) ] ; if mod-14 diff in {0,7} (apex tie), orient by
    whether the REAL gap |s_i - s_j| is < (max speed)/2 (a magnitude cutoff)."""
    k = len(S); mx = max(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            d = (S[i] - S[j]) % 14
            if 0 < d < 7:
                A[i][j] = 1
            elif 7 < d < 14:
                A[j][i] = 1
            else:  # d in {0,7}: apex tie -> magnitude cutoff
                gap = abs(S[i] - S[j])
                if 2*gap < mx:  # 'close' -> smaller index beats (use speed)
                    if S[i] < S[j]: A[i][j] = 1
                    else: A[j][i] = 1
                else:  # 'far' -> larger speed beats (REVERSED) -> this is the magnitude signal
                    if S[i] > S[j]: A[i][j] = 1
                    else: A[j][i] = 1
    return A

def reciprocal_stack(S):
    """self-clock: arc i->j iff sum over each speed s of [frac((s_i-s_j)/s) in (0,1/2)] majority.
    Magnitude-aware because 1/s scales differ per speed."""
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            votes = 0; tot = 0
            for s in S:
                rel = (F(S[i]-S[j], s)) % 1
                if rel == 0 or rel == F(1, 2): continue
                tot += 1
                if rel < F(1, 2): votes += 1
            if 2*votes > tot: A[i][j] = 1
            elif 2*votes < tot: A[j][i] = 1
            else:
                if S[i] < S[j]: A[i][j] = 1
                else: A[j][i] = 1
    return A

CANDIDATES = {
    "winding T(3/41)   ": lambda S: winding(S, F(3, 41)),
    "divis-cutoff <s/2 ": divis_cutoff,
    "divis floor-odd   ": divis_cutoff_floor,
    "additive-hybrid   ": additive_hybrid,
    "reciprocal-stack  ": reciprocal_stack,
}

# ---------------- the BANK ----------------
def primitive(S): return reduce(gcd, S) == 1

bank = []  # (name, S, is_tight)
# tight anchors
bank.append(("AP", list(range(1, 14)), True))
bank.append(("GW", list(range(1, 12)) + [13, 24], True))
# single-swap family {1..12, 13k}: M = k/(13k+1); tight only iff that == 1/14 i.e. k=1 (=AP).
# (these are LOOSE for k>=2, give a clean loose family)
for k in range(2, 12):
    S = list(range(1, 13)) + [13 * k]
    if primitive(S): bank.append((f"{{1..12,13*{k}}}", S, False))
# the GW petal family {1..11,13, v}: tight for v=24 only; loose otherwise
for v in [16, 20, 26, 28, 36, 48, 52, 96, 120]:
    S = list(range(1, 12)) + [13, v]
    if primitive(S): bank.append((f"12->{v}", S, False))
# AP single replacements (loose unless they reproduce AP/GW) -- a few representative ones
for (rem, ins) in [(7, 14), (8, 16), (9, 18), (10, 20), (11, 22), (6, 19), (5, 18)]:
    S = sorted(set(range(1, 14)) - {rem} | {ins})
    if len(S) == 13 and primitive(S):
        bank.append((f"{rem}->{ins}", S, M_exact(S) == TIGHT))

# label tightness exactly + print bank
print("BANK (exact M, tight iff M=1/14):")
for nm, S, claimed in bank:
    M = M_exact(S)
    t = (M == TIGHT)
    flag = "" if t == claimed else "  <-- LABEL MISMATCH"
    print(f"  {nm:16s} M={str(M):8s} ({float(M):.5f}) tight={t}{flag}")

tights = [(nm, S) for nm, S, _ in bank if M_exact(S) == TIGHT]
looses = [(nm, S) for nm, S, _ in bank if M_exact(S) != TIGHT]
print(f"\n  #tight={len(tights)}  #loose={len(looses)}")

# ---------------- run each candidate, test class separation ----------------
print("\n" + "="*86)
print("SEPARATION TEST: does the invariant H separate {tight} from {loose}?")
print("="*86)
for cname, fn in CANDIDATES.items():
    tight_fps = {}
    loose_fps = {}
    tight_H = set(); loose_H = set()
    bad = False
    for nm, S in tights:
        A = fn(S)
        if not is_tournament(A): bad = True
        fp = fingerprint(A)
        tight_fps[nm] = fp; tight_H.add(fp[2])
    for nm, S in looses:
        A = fn(S)
        if not is_tournament(A): bad = True
        fp = fingerprint(A)
        loose_fps[nm] = fp; loose_H.add(fp[2])
    H_separates = tight_H.isdisjoint(loose_H)
    # also: is H CONSTANT on tight (a clean signature value)?
    print(f"\n[{cname}] valid_tournament={not bad}")
    print(f"   tight H-values: {sorted(tight_H)}")
    print(f"   loose H-values (first 12): {sorted(loose_H)[:12]}{' ...' if len(loose_H)>12 else ''}")
    print(f"   >>> H SEPARATES tight from loose: {H_separates}")
    # collisions: which loose row shares a tight H?
    if not H_separates:
        clash = sorted(tight_H & loose_H)
        bad_loose = [nm for nm, fp in loose_fps.items() if fp[2] in clash]
        bad_tight = [nm for nm, fp in tight_fps.items() if fp[2] in clash]
        print(f"   COLLISION H={clash}: tight={bad_tight} vs loose={bad_loose}")
    # show AP/GW vs apex-twins specifically
    twins = {k: v for k, v in {**tight_fps, **loose_fps}.items()
             if k in ("AP", "GW", "12->26", "12->96", "12->36")}
    for k in ("AP", "GW", "12->26", "12->36", "12->96"):
        if k in twins:
            sc, cc, hh = twins[k]
            print(f"      {k:8s}: c3={cc:4d} H={hh:>10d}")
