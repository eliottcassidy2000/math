#!/usr/bin/env python3
"""
kind-pasteur-2026-06-22 (kpswf14) THREAD 2 -- MAGNITUDE-AWARE tournaments that try to
SEPARATE the tight rows {AP, GW} from the loose q=14 rows that are apex-identical to AP.

The hard pair: AP {1..13} vs 12->26 {1..11,13,26} vs 12->96 {1..11,13,96}.
All three have residues {1..13} mod 14 => identical apex winding tournament (magnitude-blind).
The magnitude (13 vs 26 vs 96) must show up at a FAREY-NEIGHBOR scale.

Candidates (compute iso invariants for AP, GW, 12->26, 12->36, 12->96):
  (a) winding T(3/41) at Farey neighbor 3/41 of 1/14  (det[[1,3],[14,41]]=-1; mod-41 sees magnitude)
  (b) MULTI-SCALE / STACKED: i->j iff [first Farey-denominator scale a/Q where ||(s_i-s_j)*a/Q||
      crosses 1/2 boundary direction]; tie-break orient by stacked sign.
  (c) DIVISIBILITY-with-cutoff: i->j iff (s_i mod s_j) < s_j/2.
  (c') 2-adic / 7-adic valuation cutoff order.

Iso invariants: score sequence, c3, H (#Ham paths, bitmask DP), c5-ish via trace? We use the
cheap-but-strong triple (scores, c3, H). H is a near-complete iso invariant at n=13.
"""
from fractions import Fraction as F
from itertools import combinations

# ---------- H = number of Hamiltonian paths (bitmask DP) ----------
def H_paths(A):
    n = len(A)
    size = 1 << n
    dp = [[0]*n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c == 0:
                continue
            Av = A[v]
            for w in range(n):
                if mask & (1 << w):
                    continue
                if Av[w]:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[size - 1])

def scores(A):
    n = len(A)
    return tuple(sorted(sum(A[i]) for i in range(n)))

def c3(A):
    n = len(A); cnt = 0
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]: cnt += 1
        if A[a][c] and A[c][b] and A[b][a]: cnt += 1
    return cnt

def is_tournament(A):
    n = len(A)
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] + A[j][i] != 1:
                return False
    return True

# ---------- candidate tournaments ----------
def winding(S, t):
    """arc i->j iff frac((s_i-s_j)t) in (0,1/2). Ties broken by speed value."""
    k = len(S); A = [[0]*k for _ in range(k)]; ties = 0
    for i in range(k):
        for j in range(k):
            if i == j: continue
            rel = (F(S[i]-S[j]) * t) % 1
            if 0 < rel < F(1, 2): A[i][j] = 1
            elif rel > F(1, 2): A[i][j] = 0
            else:
                ties += 1
                A[i][j] = 1 if S[i] < S[j] else 0
    return A, ties

def multiscale_stacked(S, scales):
    """For each ordered pair, walk the list of scales t in order. The FIRST scale at which
    frac((s_i-s_j)t) != 1/2 exactly determines the arc by which half it falls in. If it lands
    exactly on a boundary at every scale, break by speed.
    This 'stacks' resolution: coarse scales decide most pairs, fine scales the residual ties."""
    k = len(S); A = [[0]*k for _ in range(k)]; unresolved = 0
    for i in range(k):
        for j in range(i+1, k):
            decided = False
            for t in scales:
                rel = (F(S[i]-S[j]) * t) % 1
                if rel == 0:
                    continue  # same residue at this scale -> go finer
                if rel < F(1, 2):
                    A[i][j] = 1; A[j][i] = 0; decided = True; break
                elif rel > F(1, 2):
                    A[i][j] = 0; A[j][i] = 1; decided = True; break
                else:
                    continue  # exactly 1/2 -> antipodal at this scale -> go finer
            if not decided:
                unresolved += 1
                if S[i] < S[j]: A[i][j] = 1
                else: A[j][i] = 1
    return A, unresolved

def first_crossing_scale(S, scales):
    """Variant (b'): orient by the INDEX of the first scale where the pair separates, with the
    direction at that scale. The index itself is a magnitude fingerprint."""
    return multiscale_stacked(S, scales)  # same orientation; kept for clarity

def divisibility_cutoff(S):
    """(c) arc i->j iff (s_i mod s_j) < s_j/2.  Need a tournament: must be antisymmetric.
    Note: this is NOT automatically a tournament. We RECORD the raw relation and report
    how many pairs are mutual/empty, then symmetrize by a tie-break to make a tournament."""
    k = len(S)
    raw = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            r = S[i] % S[j]
            raw[i][j] = 1 if 2*r < S[j] else 0  # s_i mod s_j < s_j/2
    # build tournament: for pair (i,j) look at raw[i][j], raw[j][i]
    A = [[0]*k for _ in range(k)]
    mutual = empty = clean = 0
    for i in range(k):
        for j in range(i+1, k):
            a, b = raw[i][j], raw[j][i]
            if a == 1 and b == 0:
                A[i][j] = 1; clean += 1
            elif a == 0 and b == 1:
                A[j][i] = 1; clean += 1
            else:
                # tie (both or neither): break by speed (smaller beats larger -> matches winding tie rule)
                if a == 1 and b == 1: mutual += 1
                else: empty += 1
                if S[i] < S[j]: A[i][j] = 1
                else: A[j][i] = 1
    return A, (clean, mutual, empty)

def padic_cutoff(S, p):
    """(c') order by p-adic valuation then residue; arc i->j iff (v_p(s_i), s_i//p^{v}) precedes.
    This is a TOTAL order => transitive tournament (acyclic) => trivial H=1. Useful as a CONTROL:
    if magnitude lives in valuation, a *windowed* version (cutoff on the quotient) is needed."""
    def key(s):
        v = 0; t = s
        while t % p == 0:
            v += 1; t //= p
        return (v, t % p)  # valuation, then unit residue mod p
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            A[i][j] = 1 if key(S[i]) < key(S[j]) else 0
    return A

# ---------- rows ----------
ROWS = {
    "AP        ": (list(range(1, 14)), "TIGHT"),
    "GW        ": (list(range(1, 12)) + [13, 24], "TIGHT"),
    "12->26    ": (list(range(1, 12)) + [13, 26], "loose M=1/12"),
    "12->36    ": (list(range(1, 12)) + [13, 36], "loose M=3/41"),
    "12->96    ": (list(range(1, 12)) + [13, 96], "loose large"),
}

def report(name, A, extra=""):
    ok = is_tournament(A)
    print(f"  {name} | tourn={ok} | scores={scores(A)} | c3={c3(A):4d} | H={H_paths(A):>10d} {extra}")

# ===================== (a) Farey-neighbor winding T(3/41) =====================
print("="*78)
print("(a) FAREY-NEIGHBOR WINDING  T(3/41)   [det[[1,3],[14,41]] = 1*41-3*14 = -1]")
print("    arc i->j iff frac((s_i-s_j)*3/41) in (0,1/2). Sees magnitude mod 41.")
print("="*78)
t = F(3, 41)
for name, (S, lab) in ROWS.items():
    A, ties = winding(S, t)
    report(name, A, f"ties={ties}  [{lab}]")

# also the pure 1/41 and a few other Farey neighbors / mediants of 1/14
print("\n  -- other scales near 1/14 for contrast --")
for tt in [F(1,41), F(2,27), F(1,13), F(1,14), F(1,15), F(4,55)]:
    print(f"   scale t={tt}:")
    for name, (S, lab) in ROWS.items():
        A, ties = winding(S, tt)
        report(name, A, f"ties={ties}")

# ===================== (b) multi-scale STACKED =====================
print("\n" + "="*78)
print("(b) MULTI-SCALE STACKED winding over Farey scales of 1/14 (coarse->fine)")
print("    scales = [1/14, 3/41, 1/13, 1/41, 1/97, ...]; first separating scale orients arc")
print("="*78)
scale_lists = {
    "14 then 41-neighbors": [F(1,14), F(3,41), F(1,41), F(1,13), F(1,97)],
    "Farey tower of 1/14 ": [F(1,14), F(1,13), F(1,15), F(3,41), F(2,27), F(1,97), F(1,193)],
    "denominators 14..50 ": [F(1,d) for d in range(14, 51)],
}
for sl_name, scales in scale_lists.items():
    print(f"\n  scales = {sl_name}:")
    for name, (S, lab) in ROWS.items():
        A, un = multiscale_stacked(S, scales)
        report(name, A, f"unresolved={un}  [{lab}]")

# ===================== (c) divisibility-with-cutoff =====================
print("\n" + "="*78)
print("(c) DIVISIBILITY-CUTOFF  arc i->j iff (s_i mod s_j) < s_j/2  (tie->smaller speed)")
print("="*78)
for name, (S, lab) in ROWS.items():
    A, (clean, mut, emp) = divisibility_cutoff(S)
    report(name, A, f"clean={clean} mutual={mut} empty={emp}  [{lab}]")

# ===================== (c') p-adic order (CONTROL) =====================
print("\n" + "="*78)
print("(c') p-ADIC ORDER (valuation then unit-residue) -- transitive control, H should be 1")
print("="*78)
for p in (2, 7):
    print(f"  p={p}:")
    for name, (S, lab) in ROWS.items():
        A = padic_cutoff(S, p)
        report(name, A, f"[{lab}]")
