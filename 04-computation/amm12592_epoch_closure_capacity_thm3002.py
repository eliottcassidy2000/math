#!/usr/bin/env python3
"""Referee for THM-3002: the epoch-closure normal form (pq)^R, its Bernstein
capacity criterion, and the verified dyadic closures at gamma = 1/2.

Frame: THM-2966 spine normal form for AMM 12592.  Block B = [m_lo, m_hi];
closure (E_B) = "sum over B's cells of delta * p^z q^o vanishes identically",
which by lane G5's Theorem G5-1 (sufficiency) assembles dyadic epochs into an
exactly fair extractor with deadline T(m) = m + 1 + d_m.

Checks:
  C1  (reduction, PROVED)  Closure of B forces  p^{m_lo-1} F = -q^{m_lo-1} G
      with F = sum_m p^{m-m_lo} Delta^0_m, G = sum_m q^{m-m_lo} Delta^1_m;
      hence q^{m_lo-1} | F, F = q^{m_lo-1} H and G = -p^{m_lo-1} H (same H).
      Verified on the audited [8,15] gamma=1/2 witness, where H = 1 EXACTLY.
  C2  (normal form)  H = 1 <=> the 0-side deficit sums to (pq)^{m_lo} and the
      1-side to -(pq)^{m_lo}: the block's whole imbalance is one middle-pair
      word 0^R 1^R against 1^R 0^R -- THM-2160's mechanism at block scale.
  C3  (mirror)  delta^1 = -delta^0 solves the 1-side whenever the 0-side
      problem (*) q^{R-1} = sum_i p^i Delta_i is solved: the sides decouple.
  C4  (capacity identity, PROVED)  max |[p^t] Delta| over the Lucas box of
      degree d equals binom(d,t) * 2^t, via
      sum_j binom(d,j) binom(d-j,t-j) = binom(d,t) 2^t.
  C5  (necessary criterion)  (*) forces, for every t,
      sum_{i<=t} binom(d_i, t-i) 2^{t-i}  >=  binom(R-1, t).
      Exact ledger: DEFICIENT for gamma < 1/2 at every scale; marginal and
      eventually DEFICIENT at gamma = 1/2 (t = 25 already at R = 64);
      uniformly ample (binding at t = 1, ratio ~1.2) for gamma >= 3/5.
  C6  (verified closures)  (*) solved exactly at gamma = 1/2 for
      R = 2, 4, 8, 16, i.e. every dyadic epoch through rows [16,31].
"""

import json
import os
import sys
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
WITNESS = os.path.join(REPO, "03-artifacts", "amm12592",
                       "laneG5-witness-epoch8-15-gamma-half.json")


def require(c, m):
    if not c:
        raise RuntimeError(m)


# ------------------------------------------------------------ polynomials --
def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def padd(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i, x in enumerate(a):
        r[i] += x
    for i, x in enumerate(b):
        r[i] += x
    return r


def pmul(a, b):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return trim(r)


def qpow(n):
    r = [1]
    for _ in range(n):
        r = pmul(r, [1, -1])
    return r


def basis_poly(d, k):
    """p^{d-k} (1-p)^k."""
    r = [0] * (d + 1)
    for s in range(k + 1):
        r[d - k + s] += comb(k, s) * (-1) ** s
    return r


def divq(a, n):
    """divide by (1-p)^n, asserting exactness"""
    a = list(a)
    for _ in range(n):
        b = [0] * (len(a) - 1)
        acc = 0
        for i in range(len(a) - 1):
            acc = a[i] + acc
            b[i] = acc
        require(a[-1] + acc == 0, "not divisible by q")
        a = b
    return trim(a)


# ------------------------------------------------ C1/C2: the [8,15] witness --
def c1_c2():
    W = json.load(open(WITNESS))
    meta, delta = W["meta"], W["delta"]
    dm = lambda m: m // 2
    LO, HI = 8, 15
    D0 = {m: [] for m in range(LO, HI + 1)}
    D1 = {m: [] for m in range(LO, HI + 1)}
    for (s, m, k, z, o), dl in zip(meta, delta):
        d = dm(m)
        require(abs(dl) <= comb(d, k) and (dl - comb(d, k)) % 2 == 0, "box/parity")
        require((z, o) == ((m + d - k, k + 1) if s == 0 else (k + 1, m + d - k)),
                "cell geometry drift")
        term = (pmul([0] * (d - k) + [1], qpow(k)) if s == 0
                else pmul([0] * k + [1], qpow(d - k)))
        tgt = D0 if s == 0 else D1
        tgt[m] = padd(tgt[m], [dl * c for c in term])
    F, G = [], []
    for m in range(LO, HI + 1):
        F = padd(F, pmul([0] * (m - LO) + [1], D0[m]))
        G = padd(G, pmul(qpow(m - LO), D1[m]))
    F, G = trim(F), trim(G)
    lhs = padd(pmul([0] * (LO - 1) + [1], F), pmul(qpow(LO - 1), G))
    require(trim(lhs) == [], "C1: p^{lo-1}F + q^{lo-1}G != 0")
    H = divq(F, LO - 1)
    require(H == [1], f"C2: H is not 1 (got {H})")
    require(trim(padd(G, pmul([0] * (LO - 1) + [1], H))) == [], "C1: G != -p^{lo-1}H")
    # C2 normal form: 0-side deficit total == (pq)^LO
    side0 = []
    for m in range(LO, HI + 1):
        side0 = padd(side0, pmul(pmul([0] * m + [1], qpow(1)), D0[m]))
    require(trim(side0) == trim(pmul([0] * LO + [1], qpow(LO))),
            "C2: 0-side deficit != (pq)^R")
    print("C1 reduction + C2 normal form on the audited [8,15] witness: "
          "H = 1, 0-side deficit = (pq)^8   OK")


# -------------------------------------------------- C4: capacity identity --
def c4():
    for d in range(0, 40):
        for t in range(0, d + 1):
            lhs = sum(comb(d, j) * comb(d - j, t - j) for j in range(t + 1))
            require(lhs == comb(d, t) * 2 ** t, f"C4 identity fails d={d},t={t}")
            # and it is attained: choose delta_k = +-binom(d,k) with the sign
            # matching (-1)^{t-d+k}, which is inside the box with right parity
            best = 0
            for k in range(d + 1):
                s = t - d + k
                if 0 <= s <= k:
                    best += comb(d, k) * comb(k, s)
            require(best == comb(d, t) * 2 ** t, "C4 attainment")
    print("C4 capacity identity max|[p^t]Delta| = binom(d,t) 2^t, d <= 39: OK")


# ------------------------------------------------- C5: the exact ledger ----
def c5():
    rows = []
    for (g1, g2) in [(1, 3), (2, 5), (1, 2), (3, 5), (2, 3), (3, 4)]:
        for R in [16, 64, 256]:
            d = [(g1 * (R + i)) // g2 for i in range(R)]
            worst, wt = None, None
            for t in range(1, R):
                cap = sum(comb(d[i], t - i) * 2 ** (t - i)
                          for i in range(min(t, R - 1) + 1) if t - i <= d[i])
                tgt = comb(R - 1, t)
                r = cap / tgt
                if worst is None or r < worst:
                    worst, wt = r, t
            rows.append((f"{g1}/{g2}", R, worst, wt))
            print(f"C5 gamma={g1}/{g2:<2} R={R:4d}: worst cap/target = {worst:12.6f}"
                  f" at t={wt:4d}  {'ample' if worst >= 1 else 'DEFICIENT'}")
    # assertions: the dichotomy
    for tag, R, w, t in rows:
        if tag in ("1/3", "2/5"):
            require(w < 1, f"C5: expected deficiency at gamma={tag}")
        if tag in ("3/5", "2/3", "3/4"):
            require(w >= 1, f"C5: expected ampleness at gamma={tag}")
    half = [(R, w) for tag, R, w, t in rows if tag == "1/2"]
    require(any(w < 1 for R, w in half if R >= 64),
            "C5: gamma=1/2 should be deficient by R=64")
    print("C5 dichotomy (deficient < 1/2 <= marginal < 3/5 <= ample): OK")


# ------------------------------------- C3/C6: solve (*) and rebuild closure --
def solve_star(R, g1=1, g2=2, D0v=0, node_cap=400000):
    d = [(g1 * (R + i)) // g2 + D0v for i in range(R)]
    nodes = [0]

    def step(sigma, di, want):
        if not sigma or abs(sigma[0]) != 1:
            return None
        res = list(sigma) + [0] * (di + 3)
        deltas = [0] * (di + 1)
        v = sigma[0]
        deltas[di] = v
        for t, c in enumerate(basis_poly(di, di)):
            res[t] -= v * c
        if di >= 1:
            cap = comb(di, di - 1)
            if want is None:
                v2 = max(-cap, min(cap, res[1]))
                if (v2 - cap) % 2:
                    v2 = v2 - 1 if v2 - 1 >= -cap else v2 + 1
            else:
                v2 = res[1] - want
                if abs(v2) > cap or (v2 - cap) % 2:
                    return None
            deltas[di - 1] = v2
            for t, c in enumerate(basis_poly(di, di - 1)):
                res[t] -= v2 * c
        for k in range(di - 2, -1, -1):
            cap = comb(di, k)
            v3 = max(-cap, min(cap, res[di - k]))
            if (v3 - cap) % 2:
                v3 = v3 - 1 if v3 - 1 >= -cap else v3 + 1
            deltas[k] = v3
            if v3:
                for t, c in enumerate(basis_poly(di, k)):
                    res[t] -= v3 * c
        if res[0] != 0:
            return None
        return deltas, trim(res[1:])

    def rec(i, sigma, acc):
        nodes[0] += 1
        if nodes[0] > node_cap:
            raise TimeoutError
        if i == R - 1:
            di = d[i]
            if not sigma or len(sigma) - 1 > di:
                return None
            res = list(sigma) + [0] * (di + 2)
            deltas = [0] * (di + 1)
            for k in range(di, -1, -1):
                cap = comb(di, k)
                want = res[di - k]
                if abs(want) > cap or (want - cap) % 2:
                    return None
                deltas[k] = want
                for t, c in enumerate(basis_poly(di, k)):
                    res[t] -= want * c
            return acc + [deltas] if not trim(res) else None
        for want in (1, -1, None):
            r = step(sigma, d[i], want)
            if r is None:
                continue
            deltas, nxt = r
            if not nxt or abs(nxt[0]) != 1:
                continue
            got = rec(i + 1, nxt, acc + [deltas])
            if got is not None:
                return got
        return None

    try:
        return rec(0, qpow(R - 1), []), d
    except TimeoutError:
        return None, d


def c3_c6():
    for R in [2, 4, 8, 16]:
        sol, d = solve_star(R)
        require(sol is not None, f"C6: (*) unsolved at R={R}")
        acc = []
        for i, deltas in enumerate(sol):
            for k, v in enumerate(deltas):
                require(abs(v) <= comb(d[i], k) and (v - comb(d[i], k)) % 2 == 0,
                        "C6 box/parity")
                if v:
                    acc = padd(acc, [v * c for c in
                                     pmul([0] * i + [1], basis_poly(d[i], k))])
        require(trim(acc) == qpow(R - 1), f"C6: (*) identity fails at R={R}")
        # C3: rebuild the FULL block closure from the mirror ansatz
        total = []
        for i, deltas in enumerate(sol):
            m = R + i
            for k, v in enumerate(deltas):
                if v:
                    total = padd(total, [v * c for c in
                                         pmul([0] * (m + d[i] - k) + [1],
                                              qpow(k + 1))])
                    total = padd(total, [-v * c for c in
                                         pmul([0] * (k + 1) + [1],
                                              qpow(m + d[i] - k))])
        require(trim(total) == [], f"C3: mirror closure fails at R={R}")
        print(f"C6 (*) solved and C3 mirror closure verified at R={R:3d} "
              f"(epoch [{R},{2*R-1}])")


# ------------------------------- C7: parity lemma and the mod-2 clock ------
def c7():
    """C7a: for deg <= d, the Lucas box parity condition on Bernstein-d
    coefficients is EQUIVALENT to Delta(p) == 1 (mod 2) coefficientwise.
    C7b: hence the residual recursion reduces mod 2 to the depth-free map
    sigma_i = (sigma_{i-1}+1)/p over F_2[p] from sigma_{-1} = (1+p)^{R-1};
    for every dyadic R this orbit survives all R-1 steps and ends at 1."""
    import itertools
    for d in range(0, 9):
        for coeffs in itertools.product(range(-2, 3), repeat=d + 1):
            # build Delta from Bernstein coefficients, test both sides
            uni = [0] * (d + 1)
            inbox = True
            for k, dk in enumerate(coeffs):
                if abs(dk) > comb(d, k):
                    inbox = False
                    break
                if dk:
                    for t, c in enumerate(basis_poly(d, k)):
                        uni[t] += dk * c
            if not inbox:
                continue
            lucas = all((dk - comb(d, k)) % 2 == 0 for k, dk in enumerate(coeffs))
            univ1 = (uni[0] - 1) % 2 == 0 and all(c % 2 == 0 for c in uni[1:])
            require(lucas == univ1, f"C7a fails d={d} coeffs={coeffs}")
    print("C7a parity lemma (Lucas box parity <=> Delta == 1 mod 2), d <= 8: OK")

    def T(f):
        if not f or f[0] != 1:
            return None
        g = list(f)
        g[0] ^= 1
        return g[1:] if len(g) > 1 else [0]

    for r in range(1, 12):
        R = 2 ** r
        sig = [comb(R - 1, j) & 1 for j in range(R)]
        for i in range(R - 1):
            sig = T(sig)
            require(sig is not None, f"C7b clock dies at R={R}, step {i}")
        while sig and sig[-1] == 0:
            sig.pop()
        require(sig == [1], f"C7b clock ends at {sig[:6]} != 1 for R={R}")
    print("C7b mod-2 parity clock closes for every dyadic R <= 2048 "
          "(depth-free: constrains every gamma identically): OK")


if __name__ == "__main__":
    c1_c2()
    c7()
    c4()
    c5()
    c3_c6()
    print("ALL THM-3002 REFEREE CHECKS PASSED")
