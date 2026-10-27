"""Independent model checks (fable audit of opus S64), written before reading the
session's scripts.

A. sigma derivation (renewal-reward CLT, geometric valuation).
B. Heuristic exponent E(C) = max over (kappa, delta) with kappa ln2 - delta ln3 = 1,
   kappa + delta <= C of kappa H(delta/kappa) - delta ln 3, H in nats.
   Also a direct finite-X evaluation sum_{k,d: k+d <= C ln X, 2^k ~ 3^d X} C(k,d) 3^-d.
C. Proposition F constant C_F = (1 + (1 + log2 3)/2)/ln 2 and the entropy saving.
D. 3x-1 seeds r_j = (2^(2j+1)+1)/3: integrality, 3 r_j - 1 = 2^(2j+1),
   v_3(2^(2j+1)+1) = 1 + v_3(2j+1), permutation of Z/3^q for q <= 6.
E. Lagarias-Weiss constant 41.677 as the T-map limsup from the same word count.
"""
import math
from fractions import Fraction
from math import comb, log, log2, sqrt

ln2, ln3 = log(2.0), log(3.0)
out = []


def p(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    out.append(line)


def H(x):
    if x <= 0 or x >= 1:
        return 0.0
    return -x * log(x) - (1 - x) * log(1 - x)


# ---------------------------------------------------------------- A. sigma
p("# A. sigma derivation")
# a ~ Geometric on {1,2,...}, P(a) = 2^-a : E a = 2, Var a = 2  (check numerically)
Ea = sum(a * 2.0 ** -a for a in range(1, 200))
Ea2 = sum(a * a * 2.0 ** -a for a in range(1, 200))
p("E[a] = %.12f, Var[a] = %.12f" % (Ea, Ea2 - Ea * Ea))
Ec = 1 + Ea
Ed = Ea * ln2 - ln3
Cs = Ec / Ed
p("E[c] = 1 + E[a] = %.6f ; E[delta] = 2 ln2 - ln3 = ln(4/3) = %.8f ; C* = E[c]/E[delta] = 3/ln(4/3) = %.6f" % (Ec, Ed, Cs))
# Var(c - C* delta) = Var((1 - C* ln2) a + const) = (1 - C* ln2)^2 Var a
v = (1 - Cs * ln2) ** 2 * 2.0
p("Var(c - C* delta) = 2 (1 - C* ln2)^2 = %.6f ; / E[delta] = %.4f ; sigma = %.4f" % (v, v / Ed, sqrt(v / Ed)))
# numeric check of the variance by direct summation
vnum = sum(((1 + a) - Cs * (a * ln2 - ln3)) ** 2 * 2.0 ** -a for a in range(1, 200)) - (sum(((1 + a) - Cs * (a * ln2 - ln3)) * 2.0 ** -a for a in range(1, 200))) ** 2
p("direct numeric Var(c - C* delta) = %.6f (mean of c - C* delta = %.3e, should be 0)" % (vnum, sum(((1 + a) - Cs * (a * ln2 - ln3)) * 2.0 ** -a for a in range(1, 200))))
p("Mazur 523/50 = %.4f ; ratio to C* = %.6f ; i.e. %.3f%% above C*" % (523 / 50, (523 / 50) / Cs, 100 * ((523 / 50) / Cs - 1)))
p("renewal-reward CLT (e.g. Asmussen, Applied Probability and Queues, VI.3; Smith 1955):")
p("  cycles (delta_i, c_i) iid, N(t) = #cycles with cumulative delta <= t, R(t) = sum c_i:")
p("  (R(t) - r t)/sqrt(t) -> N(0, Var(c - r delta)/E delta), r = E c / E delta.  Here t = ln n.")

# ------------------------------------------------------------- E. L-W 41.677
p("")
p("# E. Lagarias-Weiss constant from the T-map word count (odd fraction p, gamma T-steps per ln n)")


def gamma_LW():
    # solve gamma (ln2 - H(p)) = 1 with p = (ln2 - 1/gamma)/ln3 (T-map: each step halves once; odd steps x3)
    lo, hi = 7.0, 60.0
    for _ in range(200):
        g = 0.5 * (lo + hi)
        pp = (ln2 - 1.0 / g) / ln3
        f = g * (ln2 - H(pp)) - 1.0
        if f < 0:
            lo = g
        else:
            hi = g
    return g, (ln2 - 1.0 / g) / ln3


g, pp = gamma_LW()
p("gamma_BP (T-steps) = %.6f, odd fraction p = %.6f ; ordinary-step version gamma (1+p) = %.4f" % (g, pp, g * (1 + pp)))
p("typical T-steps per ln n = 2/ln(4/3) = %.5f ; ordinary = 3/ln(4/3) = %.5f" % (2 / Ed, 3 / Ed))

# ------------------------------------------------------------- B. exponent E(C)
p("")
p("# B. heuristic exponent E(C)")
p("model: #{n < X : k halvings, d odd steps, 2^k = 3^d X} ~ C(k,d) 3^-d ; k = kappa ln X, d = delta ln X")
p("exponent(kappa, delta) = kappa H(delta/kappa) - delta ln3 with kappa ln2 - delta ln3 = 1")


def exponent_at_p(pp):
    kappa = 1.0 / (ln2 - pp * ln3)
    delta = pp * kappa
    return kappa * H(pp) - delta * ln3, kappa, delta


def E_of_C(C):
    # maximise over p in (0, ln2/ln3) subject to kappa (1 + p) <= C
    best = (-1e9, None, None, None)
    M = 200000
    for i in range(1, M):
        pp = i / M * (ln2 / ln3)
        e, kappa, delta = exponent_at_p(pp)
        if kappa + delta <= C + 1e-12 and e > best[0]:
            best = (e, pp, kappa, delta)
    return best


e, kappa, delta = exponent_at_p(0.5)
p("at p = delta/kappa = 1/2: E = %.6f, kappa = %.5f, delta = %.5f, tau/lnX = kappa + delta = %.5f" % (e, kappa, delta, kappa + delta))
# is p = 1/2 the unconstrained maximum?  scan
best = max((exponent_at_p(i / 10000 * (ln2 / ln3))[0], i / 10000 * (ln2 / ln3)) for i in range(1, 10000))
p("unconstrained max over p of the exponent: E = %.6f at p = %.5f" % best)
p("1/(ln2 - ln3/2) = %.5f" % (1 / (ln2 - ln3 / 2)))
# closed form on the boundary kappa(1+p) = C : p = (C ln2 - 1)/(1 + C ln3)
for C in (2.0, 2.5, 3.0, 3.307, 3.3074, 4.0, 5.0, 6.0, 8.0, 10.0, 10.428, 10.46, 13.0):
    e, pp, kappa, delta = E_of_C(C)
    pb = (C * ln2 - 1) / (1 + C * ln3)
    eb, kb, db = exponent_at_p(pb) if pb < 0.5 else exponent_at_p(0.5)
    p("C = %7.4f : E(C) grid = %.4f (p = %.4f, kappa = %.4f, delta = %.4f) | boundary closed form p_b = %.4f -> E = %.4f" % (C, e, pp, kappa, delta, pb, eb))

# finite-X direct sum (check the exponent empirically from the model, not from data)
p("")
p("finite-X direct model sum S(X, C) = sum_{d} C(k_d, d) 3^-d over k_d = round((ln X + d ln3)/ln2) with k_d + d <= C ln X")
for C in (3.307, 8.0):
    prev = None
    for k0 in (40, 60, 80, 100, 140):
        lnX = k0 * ln2
        S = 0.0
        for d in range(0, 2 * k0):
            k = round((lnX + d * ln3) / ln2)
            if k + d <= C * lnX:
                S += comb(k, d) * 3.0 ** (-d)
        slope = "" if prev is None else "  local exponent = %.4f" % ((log(S) - prev[1]) / (lnX - prev[0]))
        p("  C = %.3f, X = 2^%d : S = %.4e, ln S / ln X = %.4f%s" % (C, k0, S, log(S) / lnX, slope))
        prev = (lnX, log(S))

# ------------------------------------------------------------- C. Proposition F constant
p("")
p("# C. Proposition F constant and entropy saving")
CF = (1 + (1 + log2(3)) / 2) / ln2
p("C_F = (1 + (1 + log2 3)/2)/ln 2 = %.6f" % CF)
p("tau >= k0 + d0 + log2(m0) with m0 >= 3^d0  ->  tau >= k0 + d0 (1 + log2 3) ; with d0 = p k0, k0 ~ log2 X: tau/ln X >= (1 + p (1 + log2 3))/ln 2")
for pr in (0.5, 0.49, 0.48, 0.45, 0.4, 0.3):
    bound = (1 + pr * (1 + log2(3))) / ln2
    saving = 1 - (H(pr) / ln2)  # binary entropy in bits
    p("  p = %.2f : lower bound tau/ln X >= %.4f ; #{words of weight <= p k0} <= 2^(k0 H2(p)), H2 = %.5f, exponent saving c = %.5f" % (pr, bound, H(pr) / ln2, saving))
p("so for C < C_F choose p < 1/2 with (1 + p(1+log2 3))/ln2 > C: count of n in [X/2,X) with tau < C ln n is <= 2 * 2^(k0 H2(p)) = O(X^(1-c)).")
# entropy bound check: sum_{d <= p k} C(k,d) <= 2^(k H2(p)) for p <= 1/2
worst = 0.0
for k in range(1, 400):
    for d in range(0, k // 2 + 1):
        pr = d / k
        s = sum(comb(k, j) for j in range(0, d + 1))
        ratio = log2(s) - k * (H(pr) / ln2)
        worst = max(worst, ratio)
p("entropy bound log2(sum_{j<=pk} C(k,j)) - k H2(p) max over k < 400, p <= 1/2 : %.6f (must be <= 0)" % worst)

# ------------------------------------------------------------- D. 3x-1 seeds
p("")
p("# D. 3x-1 seeds r_j = (2^(2j+1) + 1)/3")


def v3(n):
    c = 0
    while n % 3 == 0:
        n //= 3
        c += 1
    return c


ok_int = all((2 ** (2 * j + 1) + 1) % 3 == 0 and ((2 ** (2 * j + 1) + 1) // 3) % 2 == 1 for j in range(0, 2000))
ok_rel = all(3 * ((2 ** (2 * j + 1) + 1) // 3) - 1 == 2 ** (2 * j + 1) for j in range(0, 2000))
ok_lte = all(v3(2 ** (2 * j + 1) + 1) == 1 + v3(2 * j + 1) for j in range(0, 2000))
p("r_j odd integer for j < 2000: %s ; 3 r_j - 1 = 2^(2j+1): %s ; v3(2^(2j+1)+1) = 1 + v3(2j+1): %s" % (ok_int, ok_rel, ok_lte))
# (2^(2j+1)+1)/3 mod 2: 2^(2j+1) = 2 mod 4... 2^(2j+1)+1 = 3 mod 4 ; divide by 3 (odd) -> odd. fine.
for q in range(1, 7):
    m = 3 ** q
    res = sorted(((2 ** (2 * j + 1) + 1) // 3) % m for j in range(0, m))
    p("q = %d : r_0..r_(3^q-1) permute Z/3^q: %s" % (q, res == list(range(m))))
# also the general statement r_j mod 3^q depends on j mod 3^q with period exactly 3^q?  (order of 4 mod 3^(q+1) is 3^q)
for q in range(1, 7):
    m = 3 ** q
    # r_j mod 3^q determined by 2^(2j+1) mod 3^(q+1); period of 4^j mod 3^(q+1) is ord(4) = 3^q
    o = 1
    x = 4 % (3 * m)
    while x != 1:
        x = (x * 4) % (3 * m)
        o += 1
    p("  ord_(3^%d)(4) = %d = 3^%d : %s" % (q + 1, o, q, o == m))
# (3y-1)/2^a <= 2y for y >= 1, a >= 1
p("(3y - 1)/2^a <= (3y-1)/2 < 3y/2 <= 2y for every odd y >= 1, a >= 1: True (and T(y) = (3y-1)/2 <= 2y)")
# S(y) <= 2y where S is the (3y-1)/2^a map: check numerically
ok = all(((3 * y - 1) >> ((3 * y - 1 & -(3 * y - 1)).bit_length() - 1)) <= 2 * y for y in range(1, 100000, 2))
p("numeric: S(y) = (3y-1)/2^(v2(3y-1)) <= 2y for odd y < 10^5: %s" % ok)

with open("audit_model_checks.out", "w") as f:
    f.write("\n".join(out) + "\n")
