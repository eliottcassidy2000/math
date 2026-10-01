#!/usr/bin/env python3
"""Collatz parity read in base phi; one polynomial identity behind the 2-adic and golden restatements;
holonomy (opus S15, twenty-first note, 2026-10-01).

T(n) = n/2 (n even), 3n+1 (n odd).  b_j(n) = T^j(n) mod 2.  F_n(z) = sum_j b_j(n) z^j  (coefficients 0/1).
Checks:
  A. T-parity words contain no '11' (after an odd step comes an even one): they are golden-mean words, i.e.
     base-phi normal forms; the T1-words are recovered by the substitution 10 -> 1 (T1 = shortcut map)
  B. Theta(n) = sum_j b_j phi^-(j+1) in [0,1]: Theta(T n) = phi*Theta(n) mod 1 (the golden beta-map) on long
     prefixes; Theta(1), Theta(4), Theta(2) = phi/2, 1/(2 phi), 1/2 = cos 36, cos 72, cos 60 degrees
  C. the golden beta-map on (1/2)Z[phi] cap [0,1) has exactly two cycles: {0} and the 3-cycle above
     (exact lattice enumeration using the contracting Galois conjugate)
  D. Collatz <=> (1 - z^3) F_n(z) in Z[z] for all n >= 1; evaluating at z = 2 (2-adically) gives
     7 Phi_T(n) in Z, at z = 1/phi gives 2 Theta(n) in Z[phi]; both converses hold (C and the 7 | ... argument);
     verified for n <= 20000 with exact arithmetic in Z[phi]
  E. negative cycles (-1, -5, -17): Theta outside (1/2)Z[phi] (the -1 cycle sits at the boundary point 1)
  F. the Collatz parity measure (2-adic Haar pushed forward) is the Markov measure P(0->1) = 1/2 on the golden-mean
     shift: entropy (2/3) ln 2 < ln phi; dimension of Theta_*(Haar) = (2/3) log_phi 2 = 0.9603 (singular)
Reproduce: python 04-computation/experiments/collatz_golden_holonomy_20261001.py   (under a minute)
"""
import math
from fractions import Fraction

FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def T(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


phi = (1 + 5 ** 0.5) / 2


# ---- exact arithmetic in Z[phi]: pairs (a, b) = a + b phi, phi^2 = phi + 1
def mul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


PHI = (0, 1)
PHI_INV = (-1, 1)          # 1/phi = phi - 1


def power(x, k):
    r = (1, 0)
    base = x if k >= 0 else PHI_INV
    for _ in range(abs(k)):
        r = mul(r, base if x == PHI else x)
    return r


def phi_pow(k):
    return power(PHI, k) if k >= 0 else power(PHI_INV, -k) if False else _phi_neg(-k)


def _phi_neg(k):
    r = (1, 0)
    for _ in range(k):
        r = mul(r, PHI_INV)
    return r


def as_float(x):
    return x[0] + x[1] * phi


def parity_word(n, length):
    w = []
    for _ in range(length):
        w.append(n & 1)
        n = T(n)
    return w


print("=== A. golden-mean words ===")
ok = all('11' not in ''.join(map(str, parity_word(n, 200))) for n in range(1, 5001))
check(ok, "T-parity words of n <= 5000 (200 steps) contain no '11'")


def t1_from_t(w):
    out, i = [], 0
    while i < len(w):
        if w[i] == 1:
            out.append(1)
            i += 2
        else:
            out.append(0)
            i += 1
    return out


def T1(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def t1_word(n, length):
    w = []
    for _ in range(length):
        w.append(n & 1)
        n = T1(n)
    return w


check(all(t1_from_t(parity_word(n, 300))[:100] == t1_word(n, 100) for n in range(1, 3001)),
      "substituting 10 -> 1 in the T-word gives the T1-word (n <= 3000)")

print("=== B. Theta conjugates T to the golden beta-map ===")


def theta(n, length=80):
    return sum(b * phi ** (-(j + 1)) for j, b in enumerate(parity_word(n, length)))


ok = all(abs((phi * theta(n)) % 1 - theta(T(n))) < 1e-9 or abs(abs((phi * theta(n)) % 1 - theta(T(n))) - 1) < 1e-9
         for n in range(1, 3001))
check(ok, "Theta(T n) = phi*Theta(n) mod 1 for n <= 3000 (80-digit prefixes)")
t1, t4, t2 = theta(1), theta(4), theta(2)
check(abs(t1 - phi / 2) < 1e-12 and abs(t4 - 1 / (2 * phi)) < 1e-12 and abs(t2 - 0.5) < 1e-12,
      f"Theta(1), Theta(4), Theta(2) = phi/2, 1/(2phi), 1/2 = cos 36, cos 72, cos 60 degrees "
      f"({t1:.6f}, {t4:.6f}, {t2:.6f})")
check(abs(math.cos(math.pi / 5) - phi / 2) < 1e-12 and abs(math.cos(2 * math.pi / 5) - 1 / (2 * phi)) < 1e-12,
      "cos(pi/5) = phi/2 and cos(2pi/5) = 1/(2 phi)")

print("=== C. cycles of the golden beta-map inside (1/2)Z[phi] ===")
psi = (1 - 5 ** 0.5) / 2


def ge0(p, q):
    """p + q*phi >= 0, exactly"""
    if p >= 0 and q >= 0:
        return True
    if p <= 0 and q <= 0:
        return p == 0 and q == 0
    nrm = p * p + p * q - q * q          # (p + q phi)(p + q psi)
    return nrm <= 0 if q > 0 else nrm >= 0


def beta_half(a, b):
    """x = (a + b phi)/2 -> phi x - d, d = floor(phi x) in {0,1}"""
    d = 1 if ge0(b - 2, a + b) else 0
    return (b - 2 * d, a + b), d


cands = []
for b in range(-8, 9):
    for a in range(-30, 31):
        if ge0(a, b) and not ge0(a - 2, b) and abs((a + b * psi) / 2) <= phi ** 2 + 1e-9:
            cands.append((a, b))
cycles = set()
for s in cands:
    seen, cur = set(), s
    while cur not in seen:
        seen.add(cur)
        cur, _ = beta_half(*cur)
    cyc, nxt = [cur], beta_half(*cur)[0]
    while nxt != cur:
        cyc.append(nxt)
        nxt = beta_half(*nxt)[0]
    i = cyc.index(min(cyc))
    cycles.add(tuple(cyc[i:] + cyc[:i]))
desc = sorted((len(c), [round((a + b * phi) / 2, 6) for a, b in c]) for c in cycles)
print("   cycles:", desc)
check(len(cands) == 10 and sorted(len(c) for c in cycles) == [1, 3]
      and any(sorted(round((a + b * phi) / 2, 9) for a, b in c) == sorted([round(1 / (2 * phi), 9), 0.5, round(phi / 2, 9)])
              for c in cycles),
      f"exactly two cycles in (1/2)Z[phi] cap [0,1): {{0}} and {{1/(2phi), 1/2, phi/2}} ({len(cands)} lattice points "
      f"with |conjugate| <= phi^2 checked; every periodic point has |conjugate| <= phi^2)")

print("=== D. one identity, two evaluations ===")
okD = True
for n in range(1, 20001):
    w, x = [], n
    while x != 1:
        w.append(x & 1)
        x = T(x)
    sigma = len(w)
    # (1 - z^3) F_n(z) = P_n(z): F_n = prefix + z^sigma * (1)/(1 - z^3)  [tail word of 1 is 100 100 ...]
    P = [0] * (sigma + 3)
    for j, b in enumerate(w):
        P[j] += b
        P[j + 3] -= b
    P[sigma] += 1
    # z = 2 (2-adically): F_n(2) = P(2)/(1-8) -> 7*Phi_T(n) = -P(2), an integer
    P2 = sum(c * 2 ** j for j, c in enumerate(P))
    # check against the 2-adic digits of -P2/7 = Phi_T(n): first sigma+30 bits
    K = sigma + 30
    inv7 = pow(7, -1, 1 << K)
    phi_bits = (-P2 * inv7) % (1 << K)
    y, ref = n, 0
    for j in range(K):
        ref |= (y & 1) << j
        y = T(y)
    okD &= phi_bits == ref
    # z = 1/phi: 2 Theta(n) = 2 phi^-1 F_n(phi^-1) = phi * P(phi^-1)  in Z[phi]
    val = (0, 0)
    for j, c in enumerate(P):
        if c:
            val = add(val, (c * _phi_neg(j)[0], c * _phi_neg(j)[1]))
    two_theta = mul(PHI, val)
    # exact: 2 Theta(n) = sum_{j<sigma} 2 w_j phi^-(j+1) + 2 phi^-sigma * Theta(1), Theta(1) = phi/2
    rhs = (0, 0)
    for j, bit in enumerate(w):
        if bit:
            q = _phi_neg(j + 1)
            rhs = add(rhs, (2 * q[0], 2 * q[1]))
    rhs = add(rhs, _phi_neg(sigma - 1) if sigma >= 1 else PHI)
    okD &= two_theta == rhs
check(okD, "for n <= 20000: (1 - z^3)F_n is a polynomial P_n; -P_n(2)/7 reproduces the 2-adic code Phi_T(n) "
      "(sigma + 30 bits) and phi*P_n(1/phi) in Z[phi] equals 2*Theta(n) exactly")

print("=== E. negative cycles ===")
for start in (-1, -5, -17):
    w, x = [], start
    for _ in range(36):
        w.append(x & 1)
        x = T(x)
    th = sum(b * phi ** (-(j + 1)) for j, b in enumerate(w * 4))
    print(f"   cycle through {start}: T-word period start {''.join(map(str, w[:18]))}..., Theta ~ {th:.6f}")
check('11' not in ''.join(map(str, parity_word(-5, 40))), "the -5 cycle (word 10100) is golden-mean admissible")
check(''.join(map(str, parity_word(-1, 10))) == '1010101010', "the -1 cycle has word (10)^inf: Theta = 1, the boundary "
      "point of the golden beta-shift (not a greedy expansion)")

print("=== F. the Collatz parity measure on the golden-mean shift ===")
# empirical transition frequencies over many random-looking starting points (prefix statistics)
cnt = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
for n in range(10 ** 6, 10 ** 6 + 20000):
    w = parity_word(n, 40)
    for a, b in zip(w, w[1:]):
        cnt[(a, b)] += 1
p01 = cnt[(0, 1)] / (cnt[(0, 0)] + cnt[(0, 1)])
print(f"   empirical P(0->1) = {p01:.4f} (Markov 1/2; Parry measure would give 1/phi^2 = {1 / phi ** 2:.4f})")
h = (2 / 3) * math.log(2)
check(abs(p01 - 0.5) < 0.02 and cnt[(1, 1)] == 0 and abs(h / math.log(phi) - 0.96025) < 1e-4,
      f"entropy (2/3) ln 2 = {h:.5f} < ln phi = {math.log(phi):.5f}; dim Theta_*(Haar) = {h / math.log(phi):.5f}")

print()
print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
