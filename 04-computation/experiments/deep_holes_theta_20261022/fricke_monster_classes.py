"""The seven shell-1 Hauptmoduln are the Monster series 4A, 8A, 12H, 16C, 20F, 24I, 32A (session opus-2026-10-22-S53;
HYP-9368 (B), THM-4760).

fricke_cusp_certificate.py proves that Z_v - 24 is the normalized Hauptmodul of Gamma_0(e)+e for the D-type shell-1
points v. This script identifies it. For each class X, E_X = prod eta(delta tau)^{r_delta} is the eta quotient of the
OEIS entry for T_X (CITED; 4A and 8A derived below):
  - Ligozat: E_X is a modular function on Gamma_0(e) (weight 0, sum delta r = 0 and sum (e/delta) r = 0 mod 24,
    prod delta^r a square);
  - its orders at the cusps 1/d (d | e) are -1 at infinity and at 0 and >= 0 elsewhere;
  - r_delta = r_{e/delta}, so E_X o W_e = E_X exactly (the constant is prod delta^{-r/2} = 1).
Hence E_X - c_X is the normalized Hauptmodul of Gamma_0(e)+e, and equals Z_v - 24 by uniqueness. The script also checks
the q-expansions against the twining coefficients computed in fricke_cusp_certificate.out.
4A, 8A: the W_e-symmetric quotients with a simple pole at infinity, (eta(2)^2/(eta(1)eta(4)))^24 and
(eta(2)eta(4)/(eta(1)eta(8)))^8; their expansions q^-1 + 24 + 276q + 2048q^2 and
q^-1 + 8 + 36q + 128q^2 match OEIS A107080 (4A) and A007265 (8A).
Usage: python fricke_monster_classes.py"""
from fractions import Fraction as Fr
from math import gcd

FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


CLASSES = {
    # class: (e, {delta: r}, OEIS, twining coefficients [q^1, q^2] from the cusp certificate)
    '4A': (4, {1: -24, 2: 48, 4: -24}, 'A107080', [276, 2048]),
    '8A': (8, {1: -8, 2: 8, 4: 8, 8: -8}, 'A007265', [36, 128]),
    '12H': (12, {1: -4, 3: 4, 4: 4, 12: -4}, 'A058486', [14, 36]),
    '16C': (16, {1: -2, 2: -3, 4: 10, 8: -3, 16: -2}, 'A058516', [8, 16]),
    '20F': (20, {1: -2, 4: 2, 5: 2, 20: -2}, 'A058555', [5, 10]),
    '24I': (24, {1: -1, 2: -2, 3: -1, 4: 4, 6: 4, 8: -1, 12: -2, 24: -1}, 'A058579', [4, 6]),
    '32A': (32, {1: -2, 2: 3, 4: -1, 8: -1, 16: 3, 32: -2}, 'A058629', [2, 4]),
}
OEIS_TERMS = {   # first terms of the OEIS entries (offset -1): q^-1, q^0, q^1, ...
    '4A': [1, 0, 276, 2048, 11202, 49152],
    '8A': [1, 0, 36, 128, 386, 1024, 2488],
    '12H': [1, 0, 14, 36, 85, 180, 360],
    '16C': [1, 0, 8, 16, 34, 64, 112],
    '20F': [1, 0, 5, 10, 18, 30, 51],
    '24I': [1, 0, 4, 6, 11, 18, 28],
    '32A': [1, 0, 2, 4, 6, 8, 12],
}


def eta_quotient(r, L):
    """q-expansion (offset: q^{sum delta r / 24}) of prod eta(delta tau)^r, L terms."""
    s = [0] * L; s[0] = 1
    for d, ex in r.items():
        for n in range(1, L):
            if d * n >= L:
                break
            # multiply by (1 - q^{dn})^ex
            if ex > 0:
                for _ in range(ex):
                    for t in range(L - 1, d * n - 1, -1):
                        s[t] -= s[t - d * n]
            else:
                for _ in range(-ex):
                    for t in range(d * n, L):
                        s[t] += s[t - d * n]
    return s


def cusp_order(e, r, d):
    """Ligozat: order of prod eta(delta tau)^r at the cusp c/d (d | e), in the local parameter."""
    tot = Fr(0)
    for de, ex in r.items():
        tot += Fr(gcd(d, de) ** 2 * ex, de)
    return Fr(e, 24 * gcd(d, e // d) * d) * tot


for X, (e, r, oeis, twin) in CLASSES.items():
    assert all(e % d == 0 for d in r)
    w = sum(r.values())
    A = sum(d * x for d, x in r.items())
    B = sum((e // d) * x for d, x in r.items())
    P = Fr(1)
    for d, x in r.items():
        P *= Fr(d) ** x
    sq = P.numerator == int(round(P.numerator ** 0.5)) ** 2 and P.denominator == int(round(P.denominator ** 0.5)) ** 2
    lig = w == 0 and A % 24 == 0 and B % 24 == 0 and sq
    divs = [d for d in range(1, e + 1) if e % d == 0]
    orders = {d: cusp_order(e, r, d) for d in divs}
    poles_ok = orders[e] == -1 and orders[1] == -1 and all(orders[d] >= 0 for d in divs if d not in (1, e))
    symm = all(r.get(d, 0) == r.get(e // d, 0) for d in divs)
    ser = eta_quotient(r, 8)                       # starts at q^{A/24} = q^-1
    const = ser[1]
    T = [ser[0], 0] + ser[2:]
    match_oeis = T[:len(OEIS_TERMS[X])] == OEIS_TERMS[X][:len(T)]
    match_twin = T[2:4] == twin
    check(lig and poles_ok and symm and A == -24 and match_oeis and match_twin,
          f'{X:4s} e = {e:2d}: Ligozat {lig}; cusp orders {dict((d, str(o)) for d, o in orders.items())}; '
          f'W_e-symmetric {symm}; E - {const} = q^-1 + {T[2:6]} (OEIS {oeis} {match_oeis}, twining {twin} {match_twin})')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
