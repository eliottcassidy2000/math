#!/usr/bin/env python3
"""
kind-pasteur-2026-06-22 (kpswf14) THREAD 1 -- verify the APEX-BLINDNESS claim.

The difference-winding ROUND tournament T(t) on speeds S:
    arc i->j  iff  frac((s_i - s_j) * t) in (0, 1/2).
(matches HYP-2605 / lrc14_diffwind_local_tournament_kps-S2.py exactly; exact Fraction)

Claim (kps-S40): at the AP apex optimum t* = a/14 (the LRC denominator-14 apex), the
AP winding tournament is the REGULAR rotational tournament on 13 vertices (all scores 6,
H = #Hamiltonian-paths = 3711175, the H-MAXIMIZER = Redei/Paley-type).
Moreover the apex tournament is MAGNITUDE-BLIND: AP {1..13}, GW {1..11,13,24}, and the
loose rows 12->26, 12->36, 12->96 all reduce mod-14-apex to residues {1..13} => SAME regular
tournament => SAME H. So apex-tournament is NECESSARY-only.

H is computed with bitmask DP over subsets (n=13: 2^13 * 13 states, fast & exact).
"""
from fractions import Fraction as F
from itertools import combinations

# ---------- tournaments ----------
def winding_tournament(S, t):
    """A[i][j]=1 iff arc i->j present.  t a Fraction.  Ties (rel in {0,1/2}) flagged."""
    k = len(S)
    A = [[0]*k for _ in range(k)]
    ties = 0
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            rel = (F(S[i]-S[j]) * t) % 1
            if 0 < rel < F(1, 2):
                A[i][j] = 1
            elif rel > F(1, 2):
                A[i][j] = 0
            else:
                # tie: rel==0 (same residue) or rel==1/2 (antipodal). break by speed value.
                ties += 1
                A[i][j] = 1 if S[i] < S[j] else 0
    return A, ties

# ---------- H = number of Hamiltonian paths (bitmask DP) ----------
def H_paths(A):
    """Count directed Hamiltonian paths in tournament A via subset DP. O(2^n * n^2)."""
    n = len(A)
    # dp[mask][v] = # paths covering exactly 'mask', ending at v
    size = 1 << n
    dp = [[0]*n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        # skip empty / singletons handled
        for v in range(n):
            c = row[v]
            if c == 0:
                continue
            # extend by w not in mask with arc v->w
            Av = A[v]
            for w in range(n):
                if mask & (1 << w):
                    continue
                if Av[w]:
                    dp[mask | (1 << w)][w] += c
    full = size - 1
    return sum(dp[full])

# ---------- iso invariants ----------
def scores(A):
    n = len(A)
    return tuple(sorted(sum(A[i]) for i in range(n)))

def c3(A):
    n = len(A)
    cnt = 0
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cnt += 1
        if A[a][c] and A[c][b] and A[b][a]:
            cnt += 1
    return cnt

def score_variance_x_n(A):
    # n * sum (s_i)^2 style fingerprint kept simple: just return multiset of scores already
    pass

# ---------- candidate rows ----------
AP   = list(range(1, 14))                 # {1..13}
GW   = list(range(1, 12)) + [13, 24]      # {1..11, 13, 24}
L26  = list(range(1, 12)) + [13, 26]      # 12 -> 26  (loose; residue-liar q-issue)
L36  = list(range(1, 12)) + [13, 36]      # 12 -> 36  (loose; M=3/41 off apex)
L96  = list(range(1, 12)) + [13, 96]      # 12 -> 96  (loose; large)

def residues_mod(S, q):
    return tuple(sorted(s % q for s in S))

print("=== ROW residues mod 14 (the apex clock) ===")
for name, S in [("AP", AP), ("GW", GW), ("12->26", L26), ("12->36", L36), ("12->96", L96)]:
    print(f"  {name:8s} S={S}")
    print(f"           res mod14 = {residues_mod(S,14)}")

# The apex optimum for AP is t* = a/14 with a coprime to 14.  Use a=1 (any unit gives an
# iso copy by relabeling of the rotation).  We verify all a in (Z/14)^* give SAME iso class.
print("\n=== APEX tournament at t* = a/14, a in (Z/14)^* ===")
units14 = [a for a in range(1, 14) if F(a, 14).denominator == 14]  # gcd(a,14)=1
print("  units mod 14:", units14)

def canon_small(A):
    """canonical form via brute permutations -- ONLY for sanity at tiny n. NOT used at n=13."""
    from itertools import permutations
    n = len(A)
    best = None
    for p in permutations(range(n)):
        key = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or key < best:
            best = key
    return best

for name, S in [("AP", AP), ("GW", GW), ("12->26", L26), ("12->36", L36), ("12->96", L96)]:
    print(f"\n  --- {name} ---")
    Hs = set()
    scs = set()
    c3s = set()
    tie_report = set()
    for a in units14:
        t = F(a, 14)
        A, ties = winding_tournament(S, t)
        Hs.add(H_paths(A))
        scs.add(scores(A))
        c3s.add(c3(A))
        tie_report.add(ties)
    print(f"     H over units      : {sorted(Hs)}")
    print(f"     score-seq over u  : {[s for s in scs]}")
    print(f"     c3 over units     : {sorted(c3s)}")
    print(f"     #ties (rel in 0,1/2) per unit: {sorted(tie_report)}")

print("\n=== Reference: the REGULAR ROTATIONAL tournament on 13 vertices (QR/Paley T_13-like) ===")
# rotational tournament: i->j iff (j-i) mod 13 in a 'symmetric' difference set of size 6.
# Standard regular rotational on odd n: i beats i+1,...,i+(n-1)/2 (the 'consecutive' DRT).
n = 13
RR = [[0]*n for _ in range(n)]
for i in range(n):
    for d in range(1, (n-1)//2 + 1):
        RR[i][(i+d) % n] = 1
print("   consecutive-DRT scores:", scores(RR), " c3:", c3(RR), " H:", H_paths(RR))

# Paley-type QR rotational at p=13: QR mod 13 = {1,3,4,9,10,12}
QR13 = {pow(x, 2, 13) for x in range(1, 13)}
print("   QR mod 13:", sorted(QR13))
PAL = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and (j - i) % 13 in QR13:
            PAL[i][j] = 1
print("   Paley-QR-DRT scores:", scores(PAL), " c3:", c3(PAL), " H:", H_paths(PAL))
