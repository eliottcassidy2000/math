#!/usr/bin/env python3
"""
monad-explorer-2026-07-07-S5 -- THE PERSPECTIVE FRAME and two new runner-tournament rules
(HYP-4937; owner directive: creative binary cutoffs on pair statistics; relativity of frames).

FRAME (owner's relativity directive, formalized).  With the observer as runner 0 (v=0),
the 14 runners see each other through the antisymmetric relative-velocity matrix
d_ij = v_i - v_j.  From i's frame every other runner j moves at d_ji; distances are
symmetric but their signed rates are antisymmetric.  Two exact structural facts:

  MIRROR LEMMA (local).  j,k sit at permanently mirrored positions in i's sky
  ({(v_j-v_i)t} = -{(v_k-v_i)t} for ALL t)  <=>  v_j + v_k = 2 v_i  (a 3-AP centered
  at i)  <=>  the balanced deficit vector (1,-2,1) on (v_j, v_i, v_k).  So the
  (1,-2,1) head-resonances that dominate the E[U]/mu deficit lattice (HYP-4907 d-perp
  probe) are EXACTLY the perspective-degeneracies: from the center's frame the two
  wings are one runner reflected.  The AP maximizes local mirrors (every interior
  element is a center); the GLOBAL mirror (reversal E -> max+min-E) is THM-639-A's
  reversal = tournament complementation.  Local mirrors refine the global one.

TOURNAMENT RULES (14 vertices incl. observer; every pair oriented):
  T1 CENTERING: c_ij = 1 if exists k with v_k = 2 v_i - v_j (i centers a mirror pair
     through j).  i -> j iff (c_ij, then A_ij fallback, then higher speed) wins.
     Arrows flow center -> wing.
  T2 EXCLUSIVE ATTENTION (pure difference dynamics): A_ij = meas{x in [0,1) :
     ||(v_i - v_j) x|| < 1/14  AND  ||(v_i - v_k) x|| >= 1/14 for all k != i,j}
     = the time-measure during which j is i's ONLY close companion.  Genuinely
     asymmetric (conditions on i's other neighbors, not j's).  i -> j iff
     A_ij > A_ji (i gives j more exclusive attention than j gives i; ties -> speed).
Invariants: score sequence, score variance, #3-cycles = C(14,3) - sum C(s_i,2),
SC-proxy (sorted scores == sorted complement scores), top/bottom vertices.
SEPARATION TEST: the winding tournament (mac-mini-S57) failed to separate tight from
loose.  Do T1/T2 invariants separate {AP, GW} (tight) from records/loose/hard-cores?

Tournament Analysis declaration:
  vertices: the 14 runners themselves (FIRST time in this project the LRC runners are
            literally the tournament vertices via a dynamical pair statistic);
  pairwise observable: exclusive-attention asymmetry A_ij vs A_ji / centering;
  switch/gauge: arrows toward the attended / the wing; tie path: speed order.
"""
from fractions import Fraction as F
from math import gcd
import numpy as np

BAND = 1.0 / 14

BANK = {
    "AP (tight)": list(range(1, 14)),
    "GW (tight)": list(range(1, 12)) + [13, 24],
    "record 2{1..11}+{11,13}": [2, 4, 6, 8, 10, 11, 12, 13, 14, 16, 18, 20, 22],
    "ds 2{1..12}+{13}": [2, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 22, 24],
    "S57 adversarial": [2, 6, 8, 10, 11, 12, 14, 16, 18, 20, 22, 26, 42],
    "sat block 14..26": list(range(14, 27)),
    "hardcore8 {9..13}+AP8@1009": [9, 10, 11, 12, 13] + [1009 - e for e in range(8)],
    "hardcore9 {10..13}+AP9@1009": [10, 11, 12, 13] + [1009 - e for e in range(9)],
    "random big": [61, 67, 74, 83, 89, 97, 104, 113, 122, 131, 140, 151, 163],
    "loose {2..14}": list(range(2, 15)),
}

def attention_matrix(v, res=30000):
    """A[i][j] = meas{x: ||(v_i-v_j)x||<1/14 and no other k has ||(v_i-v_k)x||<1/14}."""
    n = len(v)
    xs = (np.arange(res) + 0.5) / res
    D = np.array([[v[i] - v[j] for j in range(n)] for i in range(n)], dtype=np.float64)
    ph = np.multiply.outer(xs, D)                       # res x n x n
    close = np.abs(ph - np.rint(ph)) < BAND
    for i in range(n):
        close[:, i, i] = False
    deg = close.sum(axis=2)                             # res x n
    A = np.zeros((n, n))
    for i in range(n):
        only = (deg[:, i] == 1)
        for j in range(n):
            if i != j:
                A[i, j] = np.mean(close[:, i, j] & only)
    return A

def centering_matrix(v):
    n = len(v); S = set(v)
    c = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (2 * v[i] - v[j]) in S and 2 * v[i] - v[j] != v[i]:
                c[i, j] = 1
    return c

def build_tournament(v, rule, A=None, c=None):
    n = len(v)
    T = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if rule == "T1":
                key_i, key_j = c[i, j], c[j, i]
                if key_i == key_j and A is not None:
                    key_i, key_j = A[i, j], A[j, i]
            else:
                key_i, key_j = A[i, j], A[j, i]
            if key_i == key_j:
                win_i = v[i] < v[j]     # tie -> slower dominates (convention)
            else:
                win_i = key_i > key_j
            if win_i:
                T[i, j] = 1
            else:
                T[j, i] = 1
    return T

def invariants(T):
    n = T.shape[0]
    s = T.sum(axis=1)
    c3 = n * (n - 1) * (n - 2) // 6 - int(sum(x * (x - 1) // 2 for x in s))
    sc_proxy = sorted(s) == sorted(n - 1 - s)
    return dict(scores=sorted(int(x) for x in s), var=float(np.var(s)), c3=c3,
                sc=sc_proxy)

if __name__ == "__main__":
    print("=" * 96)
    print("PERSPECTIVE TOURNAMENTS: T1 centering / T2 exclusive-attention (14 vertices incl observer)")
    print("=" * 96)
    rows = []
    for nm, speeds in BANK.items():
        v = [0] + sorted(speeds)
        A = attention_matrix(v)
        c = centering_matrix(v)
        n_mirrors = int(c.sum())
        T1 = build_tournament(v, "T1", A=A, c=c)
        T2 = build_tournament(v, "T2", A=A)
        i1, i2 = invariants(T1), invariants(T2)
        # observer's role: its out-degree in T2 (whom it dominates in exclusive attention)
        obs_out = int(T2[0].sum())
        rows.append((nm, n_mirrors, i1, i2, obs_out))
        print(f"{nm:>28}: mirrors(c_ij)={n_mirrors:3d}  "
              f"T1: c3={i1['c3']:3d} var={i1['var']:.2f} SC~{i1['sc']}  "
              f"T2: c3={i2['c3']:3d} var={i2['var']:.2f} SC~{i2['sc']} obs_out={obs_out}")
    print()
    print("T2 score sequences (sorted):")
    for nm, nm_, i1, i2, oo in rows:
        print(f"  {nm:>28}: {i2['scores']}")
    print()
    print("SEPARATION READOUT: tight = {AP, GW}.  Do any invariants isolate them?")
    tights = {r[0] for r in rows if "tight" in r[0]}
    for inv, get in [("T1.c3", lambda r: r[2]['c3']), ("T2.c3", lambda r: r[3]['c3']),
                     ("mirrors", lambda r: r[1]), ("T2.var", lambda r: round(r[3]['var'], 2))]:
        tv = sorted(get(r) for r in rows if r[0] in tights)
        ov = sorted(get(r) for r in rows if r[0] not in tights)
        sep = tv[-1] < ov[0] or ov[-1] < tv[0]
        print(f"  {inv:>9}: tight={tv}  others range [{ov[0]},{ov[-1]}]  "
              f"{'SEPARATES' if sep else 'overlaps'}")
    print()
    print("MIRROR LEMMA numeric check: for 3 random 3-APs in the AP, verify permanent")
    print("mirroring {(vj-vi)x} = -{(vk-vi)x} at 5 random x (exact by construction):")
    import random
    random.seed(5)
    vAP = list(range(1, 14))
    for _ in range(3):
        i = random.choice(range(2, 12))
        d = random.choice(range(1, min(i - 1, 13 - i) + 1))
        vi, vj, vk = i, i - d, i + d
        ok = True
        for _ in range(5):
            x = random.random()
            a = ((vj - vi) * x) % 1.0
            b = ((vk - vi) * x) % 1.0
            ok &= abs((a + b) % 1.0) < 1e-9 or abs((a + b) % 1.0 - 1.0) < 1e-9
        print(f"  center {vi}, wings ({vj},{vk}): mirrored at all sampled x: {ok}")
    print("\nDONE.")
