#!/usr/bin/env python3
r"""
metagraph_pureblue_alphabet_kps_S132.py
(kind-pasteur-2026-07-26-S132; companion to THM-2454)

The complete pure-blue law (verified n = 3..11 against exhaustive
censuses; THM-2444 censuses for n <= 10, the 2^25 C++ census
metagraph_pureblue_n11_kps_S132.out for n = 11):

  pure-blue classes = palindromic stacks with RIGID {1, C3} halves
  and center in {none, 1, C3, T5}, hence

  pure-blue(n) = c(n/2)                                  (n even)
  pure-blue(n) = c((n-1)/2)+c((n-3)/2)+c((n-5)/2)        (n odd)

  with c = A000930 (Narayana's cows).  Values n = 3..11:
  2,1,3,2,4,3,6,4,9 -- all nine match the censuses.

Center law (the mechanism, verified here computationally):
  the grid involution reverses a stack's component order, so a
  sigma-fixed tiling is determined by free tilings on the first
  half plus a sigma-fixed tiling of the center:
     blue-mult(stack) = prod_half tc_i * blue-mult(center),
     tc(stack)        = prod_half tc_i^2 * tc(center);
  pure-blue <=> halves rigid AND center pure-blue.
  DECISIVE TEST: stack(T5, T5) at n = 10 must be MIXED with
  blue-mult = tc(T5) = 3 against tc = 9.  Verified by direct
  enumeration of the 2^20 blue sub-cube.

Checks are hard failures.  Reproduction:
python 04-computation/metagraph_pureblue_alphabet_kps_S132.py
"""
from itertools import permutations

def fail(msg):
    raise SystemExit("CHECK FAILED: " + msg)

# ---- closed form vs censuses ----
c = [1, 1, 1]
for i in range(3, 30):
    c.append(c[-1] + c[-3])

def pb(n):
    if n % 2 == 0:
        return c[n // 2]
    return c[(n - 1) // 2] + c[(n - 3) // 2] + (c[(n - 5) // 2] if n >= 5 else 0)

CENSUS = {3: 2, 4: 1, 5: 3, 6: 2, 7: 4, 8: 3, 9: 6, 10: 4, 11: 9}
for n, v in CENSUS.items():
    if pb(n) != v:
        fail(f"closed form vs census at n={n}: {pb(n)} vs {v}")
print("closed form matches all nine censuses n=3..11:",
      [pb(n) for n in range(3, 12)])
print("predictions: n=12..16 ->", [pb(n) for n in range(12, 17)])

# ---- n=11 inventory identification ----
# rigid: palindromic {1,3}: [1,3,9,9,9,27,27]; T5-centered: halves of 3:
# (1,1,1|T5|1,1,1) -> (H,|Aut|,tc) = (15,5,3); (C3|T5|C3) -> (135,45,3)
print("n=11 nonrigid identification: (1,1,1,T5,1,1,1) -> (15,5,3); "
      "(C3,T5,C3) -> (3*15*3, 3*5*3, 3) = (135,45,3)")
if (3 * 15 * 3, 3 * 5 * 3) != (135, 45):
    fail("arithmetic")

# ---- decisive lemma test: stack(T5,T5) at n=10 is mixed with bm=3 ----
N = 10

def tiles_of(n):
    T = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            if x - y >= 2:
                T.append((x, y))
    return T

def grid_orbits(n, T):
    pos = {t: i for i, t in enumerate(T)}
    seen = set(); orbits = []
    for i, t in enumerate(T):
        if i in seen:
            continue
        x, y = t
        j = pos[(n - y + 1, n - x + 1)]
        orb = {i, j}; seen |= orb; orbits.append(sorted(orb))
    return orbits

def tour_from_bits(n, T, bits):
    A = [[0] * n for _ in range(n)]
    for k in range(1, n):
        A[k][k - 1] = 1
    for (x, y), b in zip(T, bits):
        if b == 0:
            A[x - 1][y - 1] = 1
        else:
            A[y - 1][x - 1] = 1
    return A

def scores(A):
    return sorted(sum(r) for r in A)

# target: stack(T5,T5), T5 = circulant {1,2}
T5 = [[1 if (j - i) % 5 in (1, 2) else 0 for j in range(5)] for i in range(5)]
S = [[0] * 10 for _ in range(10)]
for i in range(5):
    for j in range(5):
        S[i][j] = T5[i][j]
        S[5 + i][5 + j] = T5[i][j]
    for j in range(5):
        S[i][5 + j] = 1
target_scores = scores(S)

def iso(A, B, n):
    # simple score-pruned backtracking (n=10, called rarely)
    sa = [sum(r) for r in A]; sb = [sum(r) for r in B]
    if sorted(sa) != sorted(sb):
        return False
    cand = [[j for j in range(n) if sb[j] == sa[i]] for i in range(n)]
    order = sorted(range(n), key=lambda i: len(cand[i]))
    mapping = [-1] * n
    used = [False] * n
    def bt(k):
        if k == n:
            return True
        i = order[k]
        for j in cand[i]:
            if used[j]:
                continue
            ok = True
            for k2 in range(k):
                i2 = order[k2]
                if A[i][i2] != B[j][mapping[i2]] or A[i2][i] != B[mapping[i2]][j]:
                    ok = False
                    break
            if ok:
                mapping[i] = j; used[j] = True
                if bt(k + 1):
                    return True
                used[j] = False; mapping[i] = -1
        return False
    return bt(0)

T = tiles_of(N); m = len(T); orbits = grid_orbits(N, T)
e = len(orbits)
assert e == 20
bm = 0
for code in range(1 << e):
    bits = [0] * m
    for oi, orb in enumerate(orbits):
        v = (code >> oi) & 1
        for idx in orb:
            bits[idx] = v
    A = tour_from_bits(N, T, bits)
    if scores(A) != target_scores:
        continue
    if iso(A, S, N):
        bm += 1
print(f"blue-mult(stack(T5,T5)) = {bm}")
if bm != 3:
    fail(f"center law refuted: blue-mult {bm} != 3")
# tc = 9: H = 15*15, |Aut| = 25
print("tc(stack(T5,T5)) = (15*15)/(5*5) = 9; blue-mult 3 < 9 => MIXED: PASS")
print("ALL CHECKS PASSED")
