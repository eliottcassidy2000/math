"""Audit item 3: the exact law of Syrac(Z/3^J) by the recursion Syrac_J = (3 Syrac_(J-1) + 1) 2^-a (exact integers),
J = 1..7, and the largest primitive Fourier coefficient max_{3 !| xi} |E e(xi Y/3^J)| (mpmath, 30 digits) against
THM-4520 (5): 0.37792, 0.25224, 0.17700, 0.12927, 0.09611, 0.07587 at J = 2..7.  Also TV(Syrac, uniform on units)."""
from fractions import Fraction

import mpmath

mpmath.mp.dps = 30


def ordmod(b, m):
    k, x = 1, b % m
    while x != 1:
        x = x * b % m
        k += 1
    return k


def syrac_laws(Jmax):
    laws = {}
    # J = 1: Y = 2^-a mod 3 ; weights over D_1 = 2^2 - 1
    num = {0: 1}       # law of "Syrac_0" on Z/1 (a point), numerators over denominator 1
    den = 1
    for J in range(1, Jmax + 1):
        Mod = 3 ** J
        L = ordmod(2, Mod)
        h = pow(2, -1, Mod)
        D = (1 << L) - 1
        new = {}
        for y, w in num.items():
            u = (3 * y + 1) % Mod                     # 3y mod 3^J depends on y mod 3^(J-1) only
            z = u
            for a in range(1, L + 1):
                z = z * h % Mod
                new[z] = new.get(z, 0) + w * (1 << (L - a))
        num, den = new, den * D
        laws[J] = {y: Fraction(w, den) for y, w in num.items()}
        assert sum(laws[J].values()) == 1
    return laws


def main():
    thm4520 = {2: 0.37792, 3: 0.25224, 4: 0.17700, 5: 0.12927, 6: 0.09611, 7: 0.07587}
    laws = syrac_laws(7)
    tao9 = [Fraction(x, 63) for x in (0, 8, 16, 0, 11, 4, 0, 2, 22)]
    assert [laws[2].get(y, Fraction(0)) for y in range(9)] == tao9
    print('B2: Syrac(Z/9) = (0,8,16,0,11,4,0,2,22)/63 exactly (recursion)')
    # consistency: the law mod 3^(J-1) of Syrac_J equals Syrac_(J-1) (Tao's consistency)
    for J in range(2, 8):
        Mod = 3 ** (J - 1)
        red = {}
        for y, w in laws[J].items():
            red[y % Mod] = red.get(y % Mod, 0) + w
        assert red == laws[J - 1]
    print('    consistency Syrac(Z/3^J) mod 3^(J-1) = Syrac(Z/3^(J-1)) exact for J = 2..7')
    for J in range(1, 8):
        Mod = 3 ** J
        law = laws[J]
        units = [y for y in range(Mod) if y % 3]
        tv = sum(abs(law.get(y, Fraction(0)) - Fraction(1, len(units))) for y in units) / 2
        best, arg = mpmath.mpf(0), None
        if J >= 2:
            items = [(y, mpmath.mpf(w.numerator) / w.denominator) for y, w in law.items()]
            for xi in range(1, Mod):
                if xi % 3 == 0:
                    continue
                s = mpmath.fsum(w * mpmath.expjpi(2 * mpmath.mpf(xi * y) / Mod) for y, w in items)
                if abs(s) > best:
                    best, arg = abs(s), xi
        line = f'J = {J}: TV(Syrac, uniform on units) = {float(tv):.5f}'
        if J >= 2:
            line += (f'; max primitive |mu_hat| = {mpmath.nstr(best, 8)} at xi = {arg}'
                     f' (THM-4520 (5): {thm4520[J]}; agree to 5 decimals: {abs(float(best) - thm4520[J]) < 6e-6})')
        print(line)
    print('B2 DONE')


if __name__ == '__main__':
    main()
