#!/usr/bin/env python3
"""S93 lane seeds, FINISHER re-checks (written after the interruption; independent of the lane's scripts where possible).

F1  Theorem 4 (c): root of rho(s) = 1 for the towers of THM-5130, from the closed formula
      rho(s) = sum_{gamma in {5, 13}} [ sum_{m >= 0} (2 gamma 3^m - 1)^(-s) + (4 gamma)^(-s) / ((1 - 2^(-s))(1 - 3^(-s))) ].
F2  Theorem 2 (a), (b), (c) for X = TT_a, a <= 9, EVERY vertex set P, by direct enumeration of ordered set partitions
    (in a transitive tournament every vertex set carries exactly one path):
      W_X(P)  = sum over ordered pairs (e, f) of P of #{(S1, S2): nonempty, max S1 = e, min S2 = f}
      W3_X(P) = the same with a third nonempty block,   kappa = 2^a - 2 W - 2 W3
    against  W = sum over pairs of (2^(d-1) + [d = 1]),  kappa = 2^a (1 - sum phi(d)),  phi(d) = (1 + 3^(d-1)) / 2^d;
    kappa = 0 iff P is a pair at distance 1 or 2; and for r = +-1 and all 1 <= M <= 20000 (M >= 3 if r = -1) the
    transported offset r+ = (1 - W) r + kappa M is +-1 only for those pairs.
F3  b_small_slopes.out (its script was not saved): the seed / anti-seed slopes below 450 and the window counts,
    recomputed from b_census_minv.npz (table minv_10).
F4  d_two_types.out (its script was not saved in that form): number of unit pairs by type at n = 9, 10 from the npz files.
F5  one-apex: number of solutions of the twin equation of PROOFS section 6 (0-blocks) by length, from the closed
    forms only (no census): val(Y) = 1 + [y_1 = 1] + [y_n = 0] + J(y) + I(Y) + I(Z) for y = Y 00 Z.
"""
import itertools
from fractions import Fraction

import numpy as np


def rho(s):
    t = 0.0
    for g in (5, 13):
        t += sum((2 * g * 3 ** m - 1) ** (-s) for m in range(0, 400))
        t += (4 * g) ** (-s) / ((1 - 2 ** (-s)) * (1 - 3 ** (-s)))
    return t


lo, hi = 0.5, 1.0
for _ in range(80):
    mid = (lo + hi) / 2
    if rho(mid) > 1:
        lo = mid
    else:
        hi = mid
print(f"F1  root of rho(s) = 1: {lo:.6f}   rho(1) = {rho(1.0):.6f}   (PROOFS: 0.788507, 0.429511)")

# F2
ok_a = ok_b = ok_c = True
nP = 0
for a in range(1, 10):
    # Q[e][f], Q3[e][f] by enumeration of the assignments of the a vertices to blocks 1, 2, 3
    Q = [[0] * a for _ in range(a)]
    Q3 = [[0] * a for _ in range(a)]
    for asg in itertools.product(range(3), repeat=a):
        b1 = [v for v in range(a) if asg[v] == 0]
        b2 = [v for v in range(a) if asg[v] == 1]
        b3 = [v for v in range(a) if asg[v] == 2]
        if not b1 or not b2:
            continue
        e, f = b1[-1], b2[0]
        if b3:
            Q3[e][f] += 1
        else:
            Q[e][f] += 1
    for mask in range(1 << a):
        P = [v for v in range(a) if (mask >> v) & 1]
        nP += 1
        W = sum(Q[e][f] for e in P for f in P if e != f)
        W3 = sum(Q3[e][f] for e in P for f in P if e != f)
        kap = 2 ** a - 2 * W - 2 * W3
        ds = [j - i for i, j in itertools.combinations(P, 2)]
        Wf = sum(2 ** (d - 1) + (1 if d == 1 else 0) for d in ds)
        kf = Fraction(2 ** a) * (1 - sum(Fraction(1 + 3 ** (d - 1), 2 ** d) for d in ds))
        if W != Wf or kap != kf:
            ok_a = False
        neutral = (len(P) == 2 and ds[0] <= 2)
        if (kap == 0) != neutral or (kap > 0) != (len(P) <= 1):
            ok_b = False
        rX = 1 - W
        for r in (1, -1):
            M = np.arange(1 if r == 1 else 3, 20001, dtype=np.int64)
            rp = rX * r + kap * M
            unit = bool((np.abs(rp) == 1).any())
            if unit != neutral:
                ok_c = False
            if neutral and not bool((rp == -r).all()):
                ok_c = False
print(f"F2  TT_a, a <= 9, every P ({nP} gadgets): formulas for W and kappa: {ok_a}; kappa = 0 iff pair at distance 1, 2 and kappa > 0 iff |P| <= 1: {ok_b}; "
      f"unit transport only by those pairs (M <= 20000): {ok_c}")

# F3
mv = np.load("b_census_minv.npz")["minv_10"].astype(np.int64)
for sg, name in ((0, "seed"), (1, "anti-seed")):
    def is_slope(M):
        v = (M & -M).bit_length() - 1
        return mv[sg][(M >> v) >> 1] <= v
    sl = [M for M in range(1, 450) if is_slope(M)]
    wins = [(1, 8), (8, 16), (16, 32), (32, 64), (64, 128), (128, 256), (256, 450)]
    print(f"F3  {name} slopes below 450 (order <= 10, doubling closure): {len(sl)}; by window {[sum(1 for M in sl if lo <= M < hi) for lo, hi in wins]}; odd ones: {sum(1 for M in sl if M & 1)}")
    for k in (13, 14, 15):
        non = [M for M in range(1 << k, 1 << (k + 1)) if not is_slope(M)]
        print(f"      integers of [2^{k}, 2^{k + 1}) that are not {name} slopes: {len(non)}" + (f": {non}" if len(non) <= 25 else ""))
    full = [M for M in range(1, 1 << 16) if not is_slope(M)]
    # longest run of consecutive slopes below 2^16
    best = (0, 0, 0)
    start = None
    prev = 0
    for M in full + [1 << 16]:
        if M - prev - 1 > best[0]:
            best = (M - prev - 1, prev + 1, M - 1)
        prev = M
    print(f"      longest interval of consecutive {name} slopes below 2^16: [{best[1]}, {best[2]}]")
    oddk = {k: int((mv[sg][:(1 << k) // 2] == 0).sum()) for k in (13, 14, 15)}
    print(f"      odd M <= 2^k that are {name} slopes themselves: {oddk}")

# F4
for n, fn in ((9, "d_two_9.npz"), (10, "d_two_10_10.npz")):
    R = np.load(fn)[f"rec_{n}"]
    for sgn in (1, -1):
        cnt = {}
        for w, O, sg, M, H in R[R[:, 2] == sgn].tolist():
            u = (O >> n) & 1
            v = (O >> (n + 1)) & 1
            k = bin(O & ((1 << n) - 1)).count("1")
            key = ("u" if u else "") + ("v" if v else "") + f"+{k}"
            cnt[key] = cnt.get(key, 0) + 1
        print(f"F4  two-apex n = {n}, sign {sgn:+d}: {sum(cnt.values())} pairs; by type {dict(sorted(cnt.items()))}")


# F5: twin equation for 0-blocks of one-apex words, closed forms only
def I_of(bits):
    t = 0
    for i in range(len(bits)):
        if bits[i] == 1:
            for k in range(i + 1, len(bits)):
                if bits[k] == 0:
                    t += 2 ** (k - i - 1)
    return t


sols = {}
for LY in range(0, 9):
    for LZ in range(0, 9):
        for Yv in range(1 << LY):
            Y = [(Yv >> (LY - 1 - i)) & 1 for i in range(LY)]
            IY = I_of(Y)
            for Zv in range(1 << LZ):
                Z = [(Zv >> (LZ - 1 - i)) & 1 for i in range(LZ)]
                # twins sit in a maximal block of zeros of length >= 2: Y does not end with 0, Z does not start with 0
                if (Y and Y[-1] == 0) or (Z and Z[0] == 0):
                    continue
                y = Y + [0, 0] + Z
                n = len(y)
                J = sum(1 for k in range(n - 1) if y[k] == 0 and y[k + 1] == 1)
                C = (1 if y[0] == 1 else 0) + (1 if y[-1] == 0 else 0) + J + IY + I_of(Z)
                if Yv == 1 + C:
                    sols.setdefault((LY, LZ), []).append(("".join(map(str, Y)), "".join(map(str, Z))))
tot = sum(len(v) for v in sols.values())
print(f"F5  one-apex twin equation (maximal 0-block, |Y|, |Z| <= 8): {tot} solutions (Y, Z); by |Y| + |Z|: "
      f"{ {s: sum(len(v) for (a, b), v in sols.items() if a + b == s) for s in sorted(set(a + b for a, b in sols)) } }")
ex = [x for v in sols.values() for x in v][:14]
print(f"      examples (Y, Z): {ex}")
print(f"      solutions with Z empty or Z = 1^k (suffix of sinks): {sum(1 for v in sols.values() for (Y, Z) in v if set(Z) <= {'1'})}")
