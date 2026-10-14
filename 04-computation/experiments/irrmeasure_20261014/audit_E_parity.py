"""Audit E (2026-10-14) of THM-4634 (parity-number dictionary): independent checks.
(1) the cycle formula x_0 = sum_i 3^(k-1-i) 2^(d_i) / (2^P - 3^k) on the named cycles {1,2}, {-1}, {-10,-7,-5},
    the 11-cycle through -17, the rational cycle 1/5 -> 4/5 -> 2/5, with d_i = positions of odd steps counted from x_0;
(2) the parity numbers rho of these cycle elements (rho(1)=2/3, rho(2)=1/3, rho(1/5) = ? -- THM-4634 (6) says 'word 100, rho = 1/7');
(3) lowest-terms denominator 3*2^j for n <= 100000 AND j = hitting time of {1,2} (sharper than THM-4634 (4));
(4) non-injectivity of rho: rho(-2) = rho(-1/3) = 1/2 (orbit reaching -1 vs 0);
(5) the Terras bijection x mod 2^k <-> first k parities, k <= 16, plus the isometry |Q(x)-Q(y)|_2 = |x-y|_2 on samples;
(6) the KL number 1 - H(log_3 2) and the density bound liminf (odd steps)/(steps) >= log_3 2 along long non-cyclic stretches;
(7) the Jarnik exponent formula mu = 2 + limsup log a_(n+1)/log q_n on a constructed example a_(n+1) = a_n^2 (mu = 3, not 4)."""
from fractions import Fraction as Fr
from math import log, log2

def T(x):  # on rationals with odd denominator: parity = parity of the numerator
    return x / 2 if x.numerator % 2 == 0 else (3 * x + 1) / 2

def orbit_word(x0, P):
    w = []; x = x0
    for _ in range(P):
        w.append(x.numerator % 2); x = T(x)
    return w, x

def cycle_formula(word):
    P, k = len(word), sum(word)
    pos = [i for i, b in enumerate(word) if b]
    S = sum(3 ** (k - 1 - i) * 2 ** pos[i] for i in range(k))
    return Fr(S, 2 ** P - 3 ** k)

def rho_periodic(word):  # rho of a purely periodic parity vector
    P = len(word); a = int("".join(map(str, word)), 2)
    return Fr(a, 2 ** P - 1)

print("(1)-(2) cycle formula and parity numbers")
for x0, P in ((Fr(1), 2), (Fr(-1), 1), (Fr(-5), 3), (Fr(-17), 11), (Fr(1, 5), 3), (Fr(4, 5), 3), (Fr(2), 2)):
    w, back = orbit_word(x0, P)
    cyc = []; y = x0
    for _ in range(P): cyc.append(y); y = T(y)
    print("  x0=%s: word %s, closes: %s, formula gives %s (match %s), cycle %s, rho(x0) = %s" % (
        x0, "".join(map(str, w)), back == x0, cycle_formula(w), cycle_formula(w) == x0, cyc, rho_periodic(w)))
print("  => rho(1/5) = 4/7 (word 100), rho(4/5) = 1/7 (word 001), rho(2/5) = 2/7 (word 010); THM-4634 (6) and parity_number_checks.py (d) mislabel 1/7 as the rho of word 100 / x0 = 1/5")
# the 11-cycle through -17 starting at its minimum -136? check the T-period and elements
cyc = []; y = Fr(-17)
for _ in range(11): cyc.append(int(y)); y = T(y)
print("  11-cycle elements from -17:", cyc, " period word from -17:", "".join(map(str, orbit_word(Fr(-17), 11)[0])), " k =", sum(orbit_word(Fr(-17), 11)[0]))

print("\n(3) lowest-terms denominators of rho(n), n <= 100000, and j vs hitting time L of {1,2}")
def parity_number(n):
    v = []; x = Fr(n); L = 0
    while x not in (1, 2):
        v.append(x.numerator % 2); x = T(x); L += 1
    r = sum(Fr(b, 2 ** (i + 1)) for i, b in enumerate(v)) + Fr(1, 2 ** L) * (Fr(2, 3) if x == 1 else Fr(1, 3))
    return r, L
bad = 0; jneL = 0
for n in range(1, 100001):
    r, L = parity_number(n)
    d = r.denominator; j = 0
    while d % 2 == 0: d //= 2; j += 1
    if d != 3: bad += 1
    if j != L: jneL += 1
print("  denominators not of the form 3*2^j: %d; cases with j != L: %d (expected 0: the last prefix digit is the parity of 4 or 1/3, so the numerator 3a+c is odd)" % (bad, jneL))
print("  rho(27) =", parity_number(27)[0], " L =", parity_number(27)[1])

print("\n(4) non-injectivity of rho")
def Q(x, k):
    v = []; y = x
    for _ in range(k): v.append(y.numerator % 2); y = T(y)
    return v
print("  Q(-2)[:8] =", Q(Fr(-2), 8), " Q(-1/3)[:8] =", Q(Fr(-1, 3), 8), " rho(-2) = 1/2 = rho(-1/3): 0111... vs 1000...; T(-1/3) =", T(Fr(-1, 3)), " T(-2) =", T(Fr(-2)), " T(-1) =", T(Fr(-1)))

print("\n(5) Terras bijection and isometry")
def Tz(x): return x // 2 if x % 2 == 0 else (3 * x + 1) // 2
import random
for k in range(1, 17):
    seen = set()
    for x in range(2 ** k):
        v = []; y = x
        for _ in range(k): v.append(y & 1); y = Tz(y)
        seen.add(tuple(v))
    assert len(seen) == 2 ** k, k
print("  bijection Z/2^k -> {0,1}^k verified for k <= 16")
def v2(n):
    if n == 0: return 10 ** 9
    c = 0
    while n % 2 == 0: n //= 2; c += 1
    return c
random.seed(1); ok = True
for _ in range(2000):
    x = random.randrange(0, 2 ** 40); y = random.randrange(0, 2 ** 40)
    vx = []; vy = []; a, b = x, y
    for _ in range(40): vx.append(a & 1); vy.append(b & 1); a = Tz(a); b = Tz(b)
    qx = sum(bit << i for i, bit in enumerate(vx)); qy = sum(bit << i for i, bit in enumerate(vy))
    if min(v2(x - y), 40) != min(v2(qx - qy), 40): ok = False
print("  |Q(x)-Q(y)|_2 = |x-y|_2 on 2000 random pairs (40 digits):", ok)

print("\n(6) information numbers")
p = log(2) / log(3)
H = -p * log2(p) - (1 - p) * log2(1 - p)
print("  log_3 2 = %.6f, H = %.6f bits, 1 - H = %.6f bits (THM-4634 says 'approx 0.049'; correct rounding 0.050)" % (p, H, 1 - H))
# the density bound: along the T-orbit of n, k ln3 - N ln2 + ln n + sum ln(1+1/(3 n_i)) = ln T^N(n) >= 0 if T^N(n) >= 1
n = 27; x = n; k = 0; N = 0; s = 0.0
while x != 1:
    if x % 2: k += 1; s += log(1 + 1 / (3 * x))
    x = Tz(x); N += 1
print("  n=27: N=%d T-steps, k=%d odd, k/N=%.4f; identity check: ln T^N = ln n + k ln3 - N ln2 + sum ln(1+1/(3n_i)) = %.3e (should be ln 1 = 0)" % (N, k, k / N, log(n) + k * log(3) - N * log(2) + s))

print("\n(7) Jarnik/Sondow exponent on a constructed continued fraction a_(n+1) = a_n^2 (essay claimed mu = 4)")
from mpmath import mp, mpf
a = [2]
for i in range(8): a.append(a[-1] ** 2)
p0, q0, p1, q1 = 1, 0, a[0], 1
ratios = []
for i in range(1, len(a)):
    p0, q0, p1, q1 = p1, q1, a[i] * p1 + p0, a[i] * q1 + q0
    if i + 1 < len(a): ratios.append(log(a[i + 1]) / log(q1))
print("  log a_(n+1)/log q_n for a_(n+1) = a_n^2:", ["%.4f" % r for r in ratios], " -> limsup 1, mu = 3 (not 4); general a_(n+1) = a_n^c gives mu = c + 1")
