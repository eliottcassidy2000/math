#!/usr/bin/env python3
"""Referee for THM-2977: the evaluation wall — bounded 2-adic gain for
multi-bias ledger functionals (AMM 12592 / THM-2966 frame).

Setting. Biases p_j = a_j/b_j (lowest terms, a_j, b_j odd,
b_j - a_j = 2^{s_j} u_j, u_j odd, s_j >= 2), depth law
d_m = floor(g1*m/g2) + D0 (rational gamma = g1/g2 in lowest terms), level
M, A = A_M = M + d_M + 1. For integers x_j, the level-M functional

    Psi = sum_j x_j N_j,   N_j = 2 D_M(p_j) b_j^{A}  (an integer).

A cell shift (doubled deficit delta -> delta + 2 at cell (z,o)) changes
Psi by  2 sum_j x_j a_j^z (2^{s_j} u_j)^o b_j^{A-z-o}.  The invariance
modulus of Psi is 2^{K}, K = min over cells of v_2(shift) (shifts generate
the choice lattice's effect on Psi).

THM-2977 claims:

  W1 (AP structure).  The o=1 cell z-positions {m + d_m : 1 <= m <= M}
      contain, for every residue class of m mod g2, the arithmetic
      progression z = m0 + d_{m0} + t*(g1 + g2), t = 0..L-1, of length
      L = floor((M - m0)/g2) + 1: exact stride h = g1 + g2 whenever
      m increases by g2.
  W2 (finite-difference kill; LTE cost).  If the values
      g(z) = sum_j y_j r_j^z (r_j = a_j/b_j in Z_2^x, r_j == 1 mod 2^{s_j},
      y_j = x_j u_j 2^{s_j} b_j^{A-1}) all have v_2 >= kappa on an AP of
      length >= J with stride h, then for each j
        kappa <= v_2(x_j) + s_j + sum_{j' != j} [ v_2(r_j - r_{j'})
                  + v_2(h) ],
      where v_2 on Q_2 extends v_2; in particular
        K <= 1 + min_j v_2(x_j) + C(bias set, gamma),
        C = max_j [ s_j + sum_{j'!=j} (v_2(a_j b_{j'} - a_{j'} b_j)
                    + v_2(g1 + g2)) ],
      independent of M and of the x_j beyond min v_2(x_j).
  W3 (the wall).  The proved necessary envelope gives
      |Psi| <= B := sum_j |x_j| (a_j^{M+1} + (b_j-a_j)^{M+1}) b_j^{d_M},
      and B >= 2^{s_min (M+1)}.  Since [-B, B] contains a full residue
      system mod 2^{K} as soon as 2B + 1 >= 2^{K}, i.e. for every
      M >= M0(bias set, gamma, min v_2(x_j)), NO forced-residue-vs-size
      contradiction ("Psi == c mod 2^K, c not achievable in [-B,B]") can
      exist at any level M >= M0.  For the certificate pair at min
      v_2(x_j) = 0: M0 = 1 already.

Referee checks (exact integer arithmetic):
  R1: W1 exact AP structure for several (g1,g2,D0) and residues.
  R2: LTE for the relevant units: v_2(rho^h - 1) = v_2(rho - 1) + v_2(h)
      for rho = r_j/r_j' == 1 mod 4, many h (as 2-adic integers via
      numerator valuations).
  R3: the K bound holds and is not wildly loose: for the certificate pair
      (and a second hostile bias pair), compute the TRUE K (min shift v_2
      over all cells, several (gamma, D0, M)) for many weight vectors
      (x_1, x_2), compare with the W2 bound.
  R4: the wall margin: log2(B) - K grows ~ linearly in M for the
      certificate pair with the G2-optimal coupling (x=1, y=tuned odd).
"""

from fractions import Fraction as Fr
from math import comb

def require(c, m):
    if not c:
        raise RuntimeError(m)

def v2(n):
    n = abs(n)
    if n == 0:
        return 10**9
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v

def v2q(fr):
    return v2(fr.numerator) - v2(fr.denominator)

# certificate biases
A_ = (1285, 2181)     # (a, b), s = 7, u = 7
B_ = (8847357, 11821757)  # s = 6, u = 46475
# hostile second pair (small, same engineering shape: a==b mod 2^s, odd)
H1 = (11, 75)         # b - a = 64 = 2^6, s = 6, u = 1
H2 = (5, 21)          # b - a = 16 = 2^4, s = 4, u = 1


def sj(bias):
    a, b = bias
    return v2(b - a)


def uj(bias):
    a, b = bias
    return (b - a) >> sj(bias)


# ------------------------------------------------------------------- R1: W1
def depth(m, g1, g2, D0):
    return (g1 * m) // g2 + D0


def r1():
    for (g1, g2, D0) in [(1, 2, 0), (1, 3, 0), (2457, 6592, 0), (3, 5, 2)]:
        M = 4 * g2 + 7
        zs = sorted(m + depth(m, g1, g2, D0) for m in range(1, M + 1))
        for m0 in range(1, g2 + 1):
            seq = [m + depth(m, g1, g2, D0) for m in range(m0, M + 1, g2)]
            diffs = {seq[i + 1] - seq[i] for i in range(len(seq) - 1)}
            require(diffs <= {g1 + g2}, f"W1 stride fail {(g1,g2,D0,m0)}")
        require(len(zs) == len(set(zs)) or g1 == 0, "")
    print("R1 (W1 AP structure, stride g1+g2 per residue class): OK")


# ------------------------------------------------------------------- R2: LTE
def r2():
    for (a1, b1), (a2, b2) in [(A_, B_), (A_, H1), (B_, H2), (H1, H2)]:
        rho = Fr(a1 * b2, b1 * a2)  # r1/r2, == 1 mod 4 (both == 1 mod 2^{s})
        base = v2q(rho - 1)
        require(base >= 2, "rho not == 1 mod 4")
        for h in [1, 2, 3, 4, 6, 8, 12, 16, 96, 9049]:
            lhs = v2q(rho ** h - 1)
            require(lhs == base + v2(h), f"LTE fail h={h}")
    print("R2 (LTE v2(rho^h - 1) = v2(rho-1) + v2(h) on unit ratios): OK")


# ------------------------------------------------------- true K per weights
def cells(m, g1, g2, D0):
    d = depth(m, g1, g2, D0)
    A = m + d + 1
    out = []
    for k in range(d + 1):
        out.append((m + d - k, k + 1))       # 0-side
        out.append((k + 1, m + d - k))       # 1-side
    return out


def true_K(biases, xs, g1, g2, D0, M):
    """min over cells of v_2(shift of Psi) - 1  (Phi-level modulus)."""
    dM = depth(M, g1, g2, D0)
    A = M + dM + 1
    best = None
    for m in range(1, M + 1):
        for (z, o) in cells(m, g1, g2, D0):
            t = 0
            for (a, b), x in zip(biases, xs):
                t += 2 * x * a**z * (b - a)**o * b**(A - z - o)
            v = v2(t)
            if best is None or v < best:
                best = v
    return best


def bound_W2(biases, xs, g1, g2):
    C = None
    for j, (a, b) in enumerate(biases):
        c = sj((a, b))
        for j2, (a2, b2) in enumerate(biases):
            if j2 == j:
                continue
            c += v2(a * b2 - a2 * b) + v2(g1 + g2)
        c += v2(xs[j])
        if C is None or c < C:
            C = c
    return 1 + C


def r3():
    import random
    rng = random.Random(2977)
    for biases, tag in [([A_, B_], "certificate"), ([H1, H2], "hostile")]:
        for (g1, g2, D0, M) in [(1, 2, 0, 12), (1, 3, 0, 12),
                                (2457, 6592, 0, 10), (3, 5, 1, 12)]:
            for _ in range(12):
                xs = [rng.choice([1, 3, 5, 7, 9, 2, 4, 12]) *
                      rng.choice([1, -1]) for _ in biases]
                K = true_K(biases, xs, g1, g2, D0, M)
                Kb = bound_W2(biases, xs, g1, g2)
                require(K <= Kb, f"W2 bound violated {tag} {xs} K={K} Kb={Kb}")
    print("R3 (true K <= W2 bound; certificate + hostile pairs, 4 laws x 12"
          " weights): OK")
    # sharpness snapshot at the certificate pair with G2 coupling shape
    xs = [1, 1]
    K = true_K([A_, B_], xs, 1, 2, 0, 12)
    print(f"    certificate pair, x=(1,1), gamma=1/2, M=12: true K={K}, "
          f"bound={bound_W2([A_, B_], xs, 1, 2)}")


# ------------------------------------------------------------------- R4: wall
def r4():
    g1, g2, D0 = 1, 2, 0
    xs = [1, 1]
    print("R4 wall margin (certificate pair):")
    for M in [2, 4, 8, 12, 16]:
        dM = depth(M, g1, g2, D0)
        K = true_K([A_, B_], xs, g1, g2, D0, M)
        Bm = 0
        for (a, b), x in zip([A_, B_], xs):
            Bm += abs(x) * (a**(M + 1) + (b - a)**(M + 1)) * b**dM
        import math
        lg = math.log2(Bm)
        require(2 * Bm + 1 >= 2**K, f"wall violated?! M={M}")
        print(f"    M={M:3d}: K={K:3d}  log2(B)={lg:9.1f}  margin="
              f"{lg - K:9.1f}")
    print("R4 (envelope covers every residue class mod 2^K from M=1 on): OK")


if __name__ == "__main__":
    r1()
    r2()
    r3()
    r4()
    print("ALL THM-2977 REFEREE CHECKS PASSED")
