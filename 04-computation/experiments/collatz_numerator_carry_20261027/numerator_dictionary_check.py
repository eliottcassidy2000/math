"""numerator_dictionary_check.py -- mac-mini-2026-10-27-potentials

Exact checks (fractions / integers only) of the Collatz numerator-coordinate dictionary for the native
compressed ROOT seed map of 05-knowledge/results/paid_root_cycle_budget_20261026/README.md (R1):

  seeds s: odd, s != 5 (mod 8);  a = 2 if s = 1 (mod 8) else 1;  (3s+1)/2^a = B^j(y), B(x) = 4x+1, y primitive;
  P(s) = y, native clock t = a + 2j.  Numerator coordinate X = 3s+1, Y = 3y+1.

Claims checked:
 (D1) 2^t Y = 3X + 2^a for every primitive seed s < 2*10^6 (the corrected R1 identity);
 (D2) the B-lift is multiplication by 4 in the numerator coordinate: 3B(x)+1 = 4(3x+1);
      the plain Syracuse map is X -> 3X/2^v + 1 (v = nu_2(3x+1)) for every odd x < 2*10^6;
 (D3) the per-edge contraction b = 3X/(2^t Y) and price lambda = (1-b)/(1+b) equal 2^a/(6X+2^a), and
      -log b = 2 artanh(lambda) = log(1 + 2^a/(3X)) (checked as exact rational identities (1+lam)/(1-lam) = 1/b);
 (D4) along any word of d native edges, 2^(T_d) X_d = 3^d X_0 + C_d with
      C_d = sum_i 3^(d-i) 2^(T_(i-1)+a_i) > 0, and prod b_i = 3^d X_0 / (2^(T_d) X_d) (exact, random words);
 (D5) cycle forcing identity: on a closed native word (X_d = X_0), prod (1+lam_i)/(1-lam_i) = 2^T/3^d, i.e.
      sum 2 artanh(lam_i) = T log 2 - d log 3; the trivial cycle s = 1 (X = 4, a = t = 2) has lam = 1/7, b = 3/4;
 (D6) single-edge fixed points: X* = 2^a/(2^t-3), s* = (2^a - 2^t + 3)/(3(2^t-3)); for j = 0 these are the
      anchored debts c_a = 1/(2^a-3) of THM-4600 (a = 2 gives s* = 1, a = 1 gives s* = -1);
 (D7) along the P-orbit of a primitive seed (sampled below 10^6; reaching 1 is implied by the known Collatz
      verification and recorded only to exercise the map) the forcing obeys the exact identity
      sum 2 artanh(lam_i) = log(1 + C_d/(3^d X_0)), hence sum lam_i < (1/2) log(1 + C_d/(3^d X_0)) (R2 eq. (7)).
"""
from fractions import Fraction
from math import log, atanh
import random


def nu2(x):
    return (x & -x).bit_length() - 1


def primitive(s):
    return s % 2 == 1 and s % 8 != 5


def native_edge(s):
    """Return (y, a, j, t) for a primitive odd seed s (s may be negative only in the fixed-point checks)."""
    a = 2 if s % 8 == 1 else 1
    z = (3 * s + 1)
    assert z % (1 << a) == 0
    z //= (1 << a)
    assert z % 2 == 1
    j = 0
    while z % 8 == 5:          # z = 4y+1 with y odd  <=>  z = 5 (mod 8)
        z = (z - 1) // 4
        j += 1
    return z, a, j, a + 2 * j


def main():
    N = 2 * 10 ** 6
    bad = 0
    for s in range(1, N, 2):
        if not primitive(s):
            continue
        y, a, j, t = native_edge(s)
        X, Y = 3 * s + 1, 3 * y + 1
        if (1 << t) * Y != 3 * X + (1 << a):
            bad += 1
        if not primitive(y):
            bad += 1
    print(f"(D1) 2^t Y = 3X + 2^a and P(s) primitive for all primitive s < {N}: failures {bad}")

    bad = 0
    for x in range(-999, 10 ** 4):
        if 3 * (4 * x + 1) + 1 != 4 * (3 * x + 1):
            bad += 1
    for x in range(1, N, 2):
        v = nu2(3 * x + 1)
        Tx = (3 * x + 1) >> v
        X = 3 * x + 1
        if Fraction(3 * X, 1 << v) + 1 != 3 * Tx + 1:
            bad += 1
    print(f"(D2) 3B(x)+1 = 4(3x+1); Syracuse X -> 3X/2^v + 1 for odd x < {N}: failures {bad}")

    bad = 0
    for s in range(1, 2 * 10 ** 5, 2):
        if not primitive(s):
            continue
        y, a, j, t = native_edge(s)
        X, Y = 3 * s + 1, 3 * y + 1
        b = Fraction(3 * X, (1 << t) * Y)
        lam = (1 - b) / (1 + b)
        if lam != Fraction(1 << a, 6 * X + (1 << a)):
            bad += 1
        if (1 + lam) / (1 - lam) != 1 / b or 1 / b != 1 + Fraction(1 << a, 3 * X):
            bad += 1
    print(f"(D3) lambda = 2^a/(6X+2^a), (1+lam)/(1-lam) = 1/b = 1 + 2^a/(3X): failures {bad}")

    rng = random.Random(20261027)
    bad = 0
    for trial in range(2000):
        s = rng.randrange(1, 10 ** 9, 2)
        while not primitive(s):
            s = rng.randrange(1, 10 ** 9, 2)
        X0 = 3 * s + 1
        d = rng.randrange(1, 40)
        T = 0
        C = 0
        prodb = Fraction(1)
        cur = s
        ok = True
        for i in range(1, d + 1):
            y, a, j, t = native_edge(cur)
            X, Y = 3 * cur + 1, 3 * y + 1
            C = 3 * C + (1 << (T + a))      # C_i = 3 C_(i-1) + 2^(T_(i-1)+a_i)
            T += t
            prodb *= Fraction(3 * X, (1 << t) * Y)
            cur = y
            if cur == 1 and i < d:
                pass
        Xd = 3 * cur + 1
        if (1 << T) * Xd != 3 ** d * X0 + C or C <= 0:
            ok = False
        if prodb != Fraction(3 ** d * X0, (1 << T) * Xd):
            ok = False
        bad += (not ok)
    print(f"(D4) 2^T X_d = 3^d X_0 + C_d, C_d > 0, prod b = 3^d X_0/(2^T X_d) on 2000 random words: failures {bad}")

    y, a, j, t = native_edge(1)
    X = 4
    lam = Fraction(1 << a, 6 * X + (1 << a))
    b = Fraction(3 * X, (1 << t) * (3 * y + 1))
    print(f"(D5) trivial cycle: P(1) = {y}, a = {a}, t = {t}, lambda = {lam}, b = {b}, (1+lam)/(1-lam) = {(1 + lam) / (1 - lam)} = 2^t/3 = {Fraction(1 << t, 3)};"
          f" 2 artanh(1/7) - log(4/3) = {2 * atanh(1 / 7) - log(4 / 3):.3e}")

    bad = 0
    for a in (1, 2):
        for j in range(0, 6):
            t = a + 2 * j
            Xs = Fraction(1 << a, (1 << t) - 3)
            ss = (Xs - 1) / 3
            if ss != Fraction((1 << a) - (1 << t) + 3, 3 * ((1 << t) - 3)):
                bad += 1
            # fixed point of the affine edge map X -> (3X + 2^a)/2^t
            if (3 * Xs + (1 << a)) / (1 << t) != Xs:
                bad += 1
            if j == 0 and ss != Fraction(1, (1 << a) - 3):
                bad += 1
    print(f"(D6) fixed points X* = 2^a/(2^t-3), s* = (2^a-2^t+3)/(3(2^t-3)), j = 0: s* = c_a = 1/(2^a-3): failures {bad};"
          f" c_1 = {Fraction(1, 2 - 3)}, c_2 = {Fraction(1, 4 - 3)}")

    # (D7) exact identity sum 2 artanh(lam_i) = log(1 + C_d/(3^d X_0)) at 60 digits on every 97th primitive seed < 10^6;
    # the inequality sum lam_i <= (1/2) log(...) then follows termwise from artanh(lam) > lam (lam > 0).
    import mpmath
    mpmath.mp.dps = 60
    bad = 0
    worst = mpmath.mpf(0)
    count = 0
    for s in range(3, 10 ** 6, 2 * 97):
        if not primitive(s):
            continue
        X0 = 3 * s + 1
        cur = s
        T = 0
        C = 0
        d = 0
        acc = mpmath.mpf(0)
        steps = 0
        while cur != 1:
            y, a, j, t = native_edge(cur)
            X = 3 * cur + 1
            lam = mpmath.mpf(1 << a) / (6 * X + (1 << a))
            acc += 2 * mpmath.atanh(lam)
            C = 3 * C + (1 << (T + a))
            T += t
            d += 1
            cur = y
            steps += 1
            if steps > 10000:
                bad += 1
                break
        err = abs(acc - mpmath.log(1 + mpmath.mpf(C) / (mpmath.mpf(3) ** d * X0)))
        worst = max(worst, err)
        count += 1
        if err > mpmath.mpf(10) ** -50:
            bad += 1
    print(f"(D7) every P-orbit sampled ({count} primitive seeds < 10^6) reaches 1; sum 2 artanh(lam) = log(1 + C_d/(3^d X_0)) at 60 digits: failures {bad}, worst error {mpmath.nstr(worst, 3)}")


if __name__ == "__main__":
    main()
