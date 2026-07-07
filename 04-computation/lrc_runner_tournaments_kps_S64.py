#!/usr/bin/env python3
r"""
lrc_runner_tournaments_kps_S64.py   (kind-pasteur-2026-07-07-S64, HYP-4927)

RUNNER TOURNAMENTS: turn LRC pair statistics into tournaments via creative binary
cutoffs (owner directive).  The naive full-circle "who leads" rule is EXACTLY fair
(S59/S63 pair-featurelessness) -- every arrow is a coin flip.  The move that breaks
the symmetry is to CONDITION ON THE LOAD-BEARING SET: the density floor lives on
Good_{1/7}(E) = {x : maxgap{frac(v_i x)} > 1/7}, and "who leads on Good" is NOT fair.

RULE 1 -- THE GOOD-SET LEADER TOURNAMENT  T_good(E):
    w_ij := meas{ x in Good_{1/7}(E) : frac(v_i x) < frac(v_j x) }
    arrow  i -> j   iff   w_ij > w_ji   (i leads j more often on the load-bearing set).
  EXACT via order cells: within a cell the fractional-part order of all points is
  fixed, so [i before j] is constant; w_ij = sum over cells with i-before-j of
  meas(cell cap Good).  Ties broken lexicographically (recorded).

RULE 2 -- THE MOD-7 PALEY RULE (CRT / composite-14 tie-in, S63 catalog item 1):
    arrow  i -> j   iff   (v_j - v_i) mod 7  in  QR_7 = {1,2,4}
  (7 = 3 mod 4 => -1 is a non-residue => antisymmetric off the 7|d diagonal;
   tiebreak on 7|d by (v_j - v_i) mod 11 residue, then by index).  Pulls the
   tournament H-machinery onto the mod-7 half of 14 = 2*7.

TARGET THEOREM (both rules): the REVERSAL functor intertwines with COMPLEMENTATION.
  Reversal E -> E* = (max+min) - E is x -> -x on phases (S62); it reverses circular
  order, so w_ij(E*) = w_{sigma(i),sigma(j)}(E) with the order-reversing relabel
  sigma(m) = k+1-m -- i.e. T_good(E*) = sigma . T_good(E)^op.  Likewise QR(-d) =
  QNR(d) flips every mod-7 arrow.  So PALINDROMIC families (E* = E up to translation)
  map to SELF-COMPLEMENTARY tournaments -- fusing the LRC palindromic-extremizer
  phenomenon (S62) with the tournament SC theory (THM-024/SC-maximizer).

OUTPUTS:
  (0) sanity: full-circle w_ij == 1/2 (fair); Good-conditioned w_ij != 1/2.
  (1) T_good on the zoo: adjacency, score sequence, #3-cycles c3, H (Ham-path count),
      regular?/circulant?/SC?  -- exact for modest diameter.
  (2) INTERTWINING: verify T_good(E*) == relabel(T_good(E)^op) exactly; SC-status of
      palindromic families; mirror pair (1^11,4)/(4,1^11) -> complement pair?
  (3) mod-7 Paley rule: same invariants; agreement matrix vs T_good; the AP's
      T_mod7 vs the Paley tournament T_7 blown up.
  (4) THE PUNCHLINE SCAN: does H(T_good(E)) or the SC-status TRACK mu_{1/7}(E)?
      (a tournament invariant read off pair statistics that sees the density floor.)
"""
from fractions import Fraction as F
from itertools import combinations
import random

THETA = F(1, 7)

# ------------------------------------------------------------------ order cells
def cell_data(E, theta=THETA):
    """Yield (order_tuple, good_measure) per order cell of E (positive ints).
    order_tuple = the indices of E sorted by fractional part at the cell (a tie-free
    interior point); good_measure = meas(cell cap {maxgap>theta})."""
    E = list(E); k = len(E)
    diffs = [abs(E[i] - E[j]) for i in range(k) for j in range(i + 1, k)]
    dmax = max(diffs + [max(E)])
    bps = {F(0), F(1)}
    for d in range(1, dmax + 1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    out = []
    for a, b in zip(bps[:-1], bps[1:]):
        mid = (a + b) / 2
        fl = {e: (e * mid).__floor__() for e in range(k)}  # placeholder
        frac = [(E[i] * mid) - (E[i] * mid).__floor__() for i in range(k)]
        order = sorted(range(k), key=lambda i: frac[i])
        # affine gap segments for the maxgap superlevel
        flv = {i: (E[i] * mid).__floor__() for i in range(k)}
        segs = []
        for s in range(k):
            i1 = order[s]; i2 = order[(s + 1) % k]
            if s < k - 1:
                c = E[i2] - E[i1]; b0 = -(flv[i2] - flv[i1])
            else:
                c = E[order[0]] - E[order[-1]]; b0 = -(flv[order[0]] - flv[order[-1]]) + 1
            segs.append((c, F(b0)))
        pts = {a, b}
        for c, b0 in segs:
            if c != 0:
                xc = (theta - b0) / c
                if a < xc < b: pts.add(xc)
        for i in range(len(segs)):
            for j in range(i + 1, len(segs)):
                ci, bi = segs[i]; cj, bj = segs[j]
                if ci != cj:
                    xc = (bj - bi) / (ci - cj)
                    if a < xc < b: pts.add(xc)
        pts = sorted(pts)
        good = F(0)
        for u, v in zip(pts[:-1], pts[1:]):
            m2 = (u + v) / 2
            if max(c * m2 + b0 for c, b0 in segs) > theta:
                good += v - u
        out.append((tuple(order), good, b - a))
    return out

def w_matrix(E, theta=THETA):
    """w[i][j] = meas{Good : frac(v_i x) < frac(v_j x)}; and full-circle version."""
    k = len(E)
    w = [[F(0)] * k for _ in range(k)]
    wfull = [[F(0)] * k for _ in range(k)]
    for order, good, length in cell_data(E, theta):
        pos = {idx: r for r, idx in enumerate(order)}
        for i in range(k):
            for j in range(k):
                if i == j: continue
                if pos[i] < pos[j]:
                    w[i][j] += good
                    wfull[i][j] += length
    return w, wfull

# ------------------------------------------------------------------ tournaments
def T_good(E, theta=THETA):
    w, wfull = w_matrix(E, theta)
    k = len(E); A = [[0] * k for _ in range(k)]
    ties = 0
    for i in range(k):
        for j in range(i + 1, k):
            if w[i][j] > w[j][i]: A[i][j] = 1
            elif w[j][i] > w[i][j]: A[j][i] = 1
            else:
                A[i][j] = 1; ties += 1  # lex tiebreak
    return A, ties, w, wfull

QR7 = {1, 2, 4}
def T_mod7(E):
    k = len(E); A = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            d = (E[j] - E[i]) % 7
            if d in QR7: A[i][j] = 1
            elif d != 0: A[j][i] = 1
            else:
                d11 = (E[j] - E[i]) % 11
                if d11 in {1, 3, 4, 5, 9}: A[i][j] = 1
                else: A[j][i] = 1
    return A

# ------------------------------------------------------------------ invariants
def scores(A):
    k = len(A); return tuple(sorted(sum(A[i]) for i in range(k)))

def c3(A):
    k = len(A); n = 0
    for i in range(k):
        for j in range(k):
            if not A[i][j]: continue
            for l in range(k):
                if A[j][l] and A[l][i]: n += 1
    return n // 3

def ham_paths(A):
    """Number of Hamiltonian paths (directed, following arrows) by subset DP."""
    k = len(A)
    from functools import lru_cache
    # dp[mask][last] = # paths covering mask ending at last
    dp = [[0] * k for _ in range(1 << k)]
    for i in range(k): dp[1 << i][i] = 1
    for mask in range(1 << k):
        for last in range(k):
            c = dp[mask][last]
            if not c: continue
            for nxt in range(k):
                if mask & (1 << nxt): continue
                if A[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    full = (1 << k) - 1
    return sum(dp[full][i] for i in range(k))

def H_ocf(A):
    """H = number of Hamiltonian paths (= Redei); use ham_paths."""
    return ham_paths(A)

def is_regular(A):
    k = len(A); s = [sum(A[i]) for i in range(k)]
    return len(set(s)) == 1

def complement(A):
    k = len(A); return [[A[j][i] for j in range(k)] for i in range(k)]

def relabel(A, sigma):
    k = len(A); return [[A[sigma[i]][sigma[j]] for j in range(k)] for i in range(k)]

def iso_equal(A, B):
    """cheap iso test via sorted (score, local) signatures + brute perm for small k
    fallback; here use exact equality after trying the reversal relabel only."""
    return A == B

# ------------------------------------------------------------------ families
def fam(word, v1=1):
    v = [v1]
    for s in word: v.append(v[-1] + s)
    return v

def steps(E):
    E = sorted(E); return tuple(E[i+1]-E[i] for i in range(len(E)-1))

def mu_of(E, theta=THETA):
    tot = F(0)
    for order, good, length in cell_data(E, theta):
        tot += good
    return tot

print("=" * 98)
print("PART 0 -- the naive full-circle rule is EXACTLY FAIR; conditioning on Good breaks it")
print("=" * 98)
E = [1, 2, 3, 5, 8]   # small test
A, ties, w, wfull = T_good(E)
k = len(E)
full_fair = all(wfull[i][j] == wfull[j][i] for i in range(k) for j in range(i+1, k))
good_fair = all(w[i][j] == w[j][i] for i in range(k) for j in range(i+1, k))
print(f"  E={E}: full-circle w symmetric (fair)? {full_fair}   Good-conditioned w symmetric? {good_fair}")
print(f"  sample Good-conditioned w[0][1]={w[0][1]} vs w[1][0]={w[1][0]}  (meas Good={mu_of(E)})")

print()
print("=" * 98)
print("PART 1 -- T_good on the zoo (exact, modest diameter): scores, c3, H, regular/SC")
print("=" * 98)
zoo = {
    "AP {1..8}": list(range(1, 9)),
    "record-ish {1..7,9}": [1,2,3,4,5,6,7,9],
    "mirror A (1^6,3)-word": fam((1,1,1,1,1,1,3)),
    "mirror B (3,1^6)-word": fam((3,1,1,1,1,1,1)),
    "primes8": [2,3,5,7,11,13,17,19],
    "geom {1,2,4,8,16,32,64,128}": [1,2,4,8,16,32,64,128],
    "GW-ish {1..6,8}": [1,2,3,4,5,6,8],
}
inv = {}
for nm, E in zoo.items():
    A, ties, w, wf = T_good(E)
    sc = scores(A); c = c3(A); H = H_ocf(A); reg = is_regular(A)
    Ac = complement(A)
    # SC via reversal relabel sigma(m)=k-1-m
    kk = len(E); sigma = list(range(kk-1, -1, -1))
    sc_check = iso_equal(A, relabel(Ac, sigma))
    inv[nm] = (steps(E), sc, c, H, reg, sc_check, float(mu_of(E)), ties)
    palin = "PAL" if steps(E) == steps(E)[::-1] else "   "
    print(f"  {nm:26s} {palin} scores={sc} c3={c} H={H} reg={reg} SC(rev)={sc_check} mu={float(mu_of(E)):.4f} ties={ties}")

print()
print("=" * 98)
print("PART 2 -- INTERTWINING: T_good(E*) == relabel(T_good(E)^op)?  (E* = max+min - E)")
print("=" * 98)
for nm, E in zoo.items():
    Estar = sorted(max(E) + min(E) - e for e in E)
    A, *_ = T_good(E)
    As, *_ = T_good(Estar)
    kk = len(E); sigma = list(range(kk-1, -1, -1))
    predicted = relabel(complement(A), sigma)
    ok = (As == predicted)
    palin = steps(E) == steps(E)[::-1]
    tag = "SELF-COMP" if (palin and iso_equal(A, relabel(complement(A), sigma))) else ""
    print(f"  {nm:26s} T_good(E*) == relabel(comp(T_good(E)))? {ok}   {'[palindromic]' if palin else ''} {tag}")

print()
print("=" * 98)
print("PART 3 -- the MOD-7 PALEY rule: invariants + agreement with T_good + reversal intertwine")
print("=" * 98)
for nm, E in zoo.items():
    A7 = T_mod7(E)
    sc = scores(A7); c = c3(A7); H = H_ocf(A7); reg = is_regular(A7)
    kk = len(E); sigma = list(range(kk-1, -1, -1))
    # reversal on mod-7: E* differences negate; QR(-d)=QNR(d) => complement
    Estar = sorted(max(E) + min(E) - e for e in E)
    A7s = T_mod7(Estar)
    intertw = (A7s == relabel(complement(A7), sigma))
    # agreement with T_good
    Ag, *_ = T_good(E)
    agree = sum(1 for i in range(kk) for j in range(kk) if i != j and A7[i][j] == Ag[i][j]) / (kk*(kk-1))
    print(f"  {nm:26s} mod7: scores={sc} c3={c} H={H} reg={reg}  rev->comp:{intertw}  agree(T_good)={agree:.2f}")

print()
print("=" * 98)
print("PART 4 -- does a TOURNAMENT invariant of T_good track the density floor mu_{1/7}?")
print("=" * 98)
rows = sorted(((v[6], nm, v[3], v[2], v[5]) for nm, v in inv.items()))
print(f"  {'mu_1/7':>8} {'family':26s} {'H(T_good)':>10} {'c3':>5} {'SC-rev':>7}")
for mu, nm, H, c, sc in rows:
    print(f"  {mu:8.4f} {nm:26s} {H:10d} {c:5d} {str(sc):>7}")
print()
print("  (reading: if H or c3 rises with mu, T_good's tournament complexity SEES the floor;")
print("   if SC-status clusters the palindromic minimizers, the SC-maximizer analogy holds.)")
print()
print("DONE.")
