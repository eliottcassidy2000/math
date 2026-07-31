#!/usr/bin/env python3
"""Referee for THM-2976: binary-clock parity structure of the critical-run
deficit ledger (AMM 12592 / THM-2966 frame).

Objects. For a depth law d_M >= 0 and level M >= 1, the choice-independent
forced parity of the doubled homogenized deficit ledger is

    beta_M(x) = (1+x)^{A_M} + (1 + x^{M+1}) (1+x)^{d_M}   (mod 2),
    A_M = M + d_M + 1

(lane D closed form, re-derived in THM-2976 Lemma 0 from
2 D_M = 2 S_M - 1 + p^{M+1} + q^{M+1} with 2 S_M even-integral).

Claims refereed here, all in exact integer arithmetic via Lucas:

  T1 (checkpoint vanishing).  If M+1 = 2^r then beta_M == 0 identically,
      for EVERY d_M >= 0.  [Frobenius: (1+x)^{2^r+d} = (1+x^{2^r})(1+x)^d.]
  T2 (the clock).  If M+1 is NOT a power of two, the minimal positive
      forced-odd position is exactly o*(M) = 2^{v_2(M+1)}.
  T3 (corner timing).  o*(M) equals the isolated-corner position d_M + 1
      iff d_M = 2^{v_2(M+1)} - 1 (and M+1 is not a power of two).
  C4 (ladder).  For gamma = 1/(2^k - 1), D0 = 0, d_M = floor(gamma*M):
      every M = 2^r - 2^{r-k} - 1 (r > k) satisfies d_M = 2^{r-k} - 1 and
      v_2(M+1) = r-k, i.e. hits corner timing; k = 1 is the classical
      gamma = 1 dyadic clock.

Controls: hostile non-dyadic M+1 sweep (T1 does not overreach); exhaustive
(M, d) box; random large sample; beta recomputed two independent ways
(Lucas closed form vs direct F2 polynomial convolution).
"""

import random


def binom2(n: int, k: int) -> int:
    """binom(n,k) mod 2 by Lucas; 0 for k<0 or k>n."""
    if k < 0 or k > n:
        return 0
    while k:
        if (k & 1) and not (n & 1):
            return 0
        n >>= 1
        k >>= 1
    return 1


def beta_closed(M: int, d: int):
    """Forced parity vector via the closed form, positions 0..A."""
    A = M + d + 1
    out = []
    for o in range(A + 1):
        v = binom2(A, o) ^ binom2(d, o) ^ binom2(d, o - (M + 1))
        out.append(v)
    return out


def beta_conv(M: int, d: int):
    """Independent recomputation: F2 polynomial arithmetic by convolution."""
    A = M + d + 1

    def poly_pow_1px(n):  # (1+x)^n mod 2 as bit list
        return [binom2(n, i) for i in range(n + 1)]

    a = poly_pow_1px(A)
    b = poly_pow_1px(d)
    out = [0] * (A + 1)
    for i, c in enumerate(a):
        out[i] ^= c
    for i, c in enumerate(b):
        out[i] ^= c            # (1)*(1+x)^d
        if i + M + 1 <= A:
            out[i + M + 1] ^= c  # x^{M+1}*(1+x)^d
    return out


def v2(n: int) -> int:
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(msg)


def main() -> None:
    rng = random.Random(2976)

    # cross-validation of the two beta computations
    for _ in range(300):
        M = rng.randrange(1, 200)
        d = rng.randrange(0, 80)
        require(beta_closed(M, d) == beta_conv(M, d), f"beta mismatch {M},{d}")
    print("beta closed form == direct F2 convolution on 300 random (M,d): OK")

    # T1: exhaustive box + random large
    for r in range(1, 11):
        M = 2 ** r - 1
        for d in range(0, 130):
            require(all(v == 0 for v in beta_closed(M, d)), f"T1 fail {M},{d}")
    for _ in range(50):
        r = rng.randrange(1, 20)
        d = rng.randrange(0, 3000)
        M = 2 ** r - 1
        A = M + d + 1
        for _ in range(200):
            o = rng.randrange(0, A + 1)
            require(
                (binom2(A, o) ^ binom2(d, o) ^ binom2(d, o - (M + 1))) == 0,
                f"T1 rand fail {M},{d},{o}")
    print("T1 checkpoint vanishing (M+1=2^r, every d): OK")

    # T2: exact clock, exhaustive box
    for M in range(1, 260):
        if (M + 1) & M == 0:  # power of two: T1 case
            continue
        target = 2 ** v2(M + 1)
        for d in range(0, 90):
            b = beta_closed(M, d)
            pos = [o for o in range(1, len(b)) if b[o]]
            require(pos and min(pos) == target, f"T2 fail {M},{d}: {pos[:4]}")
    print("T2 clock o* = 2^(v2(M+1)) for all non-dyadic M+1 in box: OK")

    # T3: corner-timing classification on the same box
    for M in range(1, 260):
        if (M + 1) & M == 0:
            continue
        t = 2 ** v2(M + 1)
        for d in range(0, 90):
            corner = (t == d + 1)
            require(corner == (d == t - 1), "T3 tautology guard")
    # semantic check: corner position d+1 is forced-odd exactly then
    for M in range(1, 260):
        if (M + 1) & M == 0:
            continue
        for d in range(0, 90):
            b = beta_closed(M, d)
            pos = [o for o in range(1, len(b)) if b[o]]
            require((min(pos) == d + 1) == (d == 2 ** v2(M + 1) - 1),
                    f"T3 fail {M},{d}")
    print("T3 corner timing iff d = 2^(v2(M+1)) - 1: OK")

    # C4: ladder rates
    for k in range(1, 6):
        g_num, g_den = 1, 2 ** k - 1
        for r in range(k + 1, 22):
            M = 2 ** r - 2 ** (r - k) - 1
            d = (g_num * M) // g_den  # floor(gamma*M), D0 = 0
            require(d == 2 ** (r - k) - 1, f"C4 depth fail k={k},r={r}")
            require(v2(M + 1) == r - k, f"C4 clock fail k={k},r={r}")
            require(d == 2 ** v2(M + 1) - 1, f"C4 corner fail k={k},r={r}")
    print("C4 ladder gamma=1/(2^k-1) hits corner timing at M=2^r-2^(r-k)-1: OK")

    # C5 (corner-clocked rates, D0 = 0): infinitely many corner-timed levels
    # iff gamma = 1/J with J odd >= 3.  Proof in THM-2976; here we referee
    # (i) the positive family: for J in {3,5,7,9}, M = J*2^v - 1 is corner-
    #     timed for every v >= 1;
    # (ii) the forcing congruence: ANY corner-timed level M of ANY law
    #     floor(gamma*M) satisfies |gamma*(M+1) - 2^{v2(M+1)}| < 1 + gamma,
    #     so gamma*(odd part of M+1) -> 1 along any infinite hit sequence,
    #     forcing 1/gamma = that odd integer eventually;
    # (iii) hostile counts: non-unit-fraction and even-unit-fraction rates
    #     have only transient hits in a long window.
    for J in [3, 5, 7, 9]:
        for v in range(1, 18):
            M = J * 2 ** v - 1
            require(v2(M + 1) == v, f"C5i clock {J},{v}")
            require(M // J == 2 ** v - 1, f"C5i depth {J},{v}")
    print("C5(i) unit fractions 1/J, J odd: corner-timed at every M=J*2^v-1: OK")
    hits = {}
    for g_num, g_den in [(1, 2), (1, 4), (2, 5), (3, 8), (1, 3), (1, 5), (1, 1),
                         (2457, 6592)]:
        cnt = []
        for M in range(1, 200000):
            if (M + 1) & M == 0:
                continue
            v = v2(M + 1)
            if (g_num * M) // g_den == 2 ** v - 1:
                cnt.append(M)
                # C5(ii) forcing congruence, exact integer form:
                # |gamma*(M+1) - 2^v| < 1 + gamma  <=>
                # |g_num*(M+1) - 2^v*g_den| < g_den + g_num
                require(abs(g_num * (M + 1) - 2 ** v * g_den) < g_den + g_num,
                        f"C5ii fail {g_num}/{g_den} M={M}")
        hits[f"{g_num}/{g_den}"] = (len(cnt), cnt[-1] if cnt else None)
    print("C5 hits (count, last) M<200000:", hits)
    require(hits["1/3"][0] > 15 and hits["1/5"][0] > 12, "C5 positive density")
    # J even or gamma non-unit-fraction: only transient hits (finitely many,
    # PROVED via the forcing congruence; here: none past M=100 in the window)
    for key in ["1/2", "1/4", "2/5", "3/8", "2457/6592"]:
        require(hits[key][1] is None or hits[key][1] < 100,
                f"C5 transient fail {key}")

    print("ALL REFEREE CHECKS PASSED")


if __name__ == "__main__":
    main()
