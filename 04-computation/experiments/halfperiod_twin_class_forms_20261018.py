#!/usr/bin/env python3
"""Which modular form governs each half-period twin class? (S43, after halfperiod_twin_classes_20261018.py)

For every sequence u in the twin classes found by halfperiod_twin_classes_20261018.py, test
    u((p-1)/2) == eps(p) * s(p) * a_f(p)   (mod p^2)    for all primes 5 <= p <= PMAX,
with f in a library of eta-quotients (the 21 M24 cycle-shape eta-products, plus the CM forms 27a, 32a, 64a and the
weight-3 level-12 form eta(2t)^3 eta(6t)^3), eps a quadratic character in {1, (-1/p), (2/p), (-2/p), (3/p), (-3/p),
(6/p), (-6/p)}, and s(p) = 2^(i n) 3^(j n) with n = (p-1)/2, -6 <= i <= 6, -3 <= j <= 3 (so 4^n = 2^(p-1) etc.).
Prints every surviving (f, eps, s) for each sequence; if none survives mod p^2, prints the mod-p survivors.
Usage: python halfperiod_twin_class_forms_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import os, sys, json
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m24_paley_golay_eta_20261018 import SHAPES
from halfperiod_twin_classes_20261018 import LIB


def eta_quotient(shape, N):
    """coefficients of prod_a prod_n (1 - q^(a n))^(e_a), up to q^N (exponents may be negative)."""
    b = [0] * (N + 1)
    b[0] = 1
    for a, e in shape.items():
        for _ in range(abs(e)):
            for n in range(1, N // a + 1):
                s = a * n
                if e > 0:
                    for m in range(N, s - 1, -1):
                        b[m] -= b[m - s]
                else:
                    for m in range(s, N + 1):
                        b[m] += b[m - s]
    return b


EXTRA = {'27a(w2,CM -3)': {3: 2, 9: 2}, '32a(w2,CM -4)': {4: 2, 8: 2}, '64a(w2,CM -4)': {4: -2, 8: 8, 16: -2},
         '12(w3,CM -3)': {2: 3, 6: 3}}


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def make_seq(name):
    fs = [LIB[c] for c in name.split(' ')]
    def u(n):
        tot = 0
        for k in range(n + 1):
            v = 1
            for f in fs:
                v *= f(n, k)
                if v == 0:
                    break
            tot += v
        return tot
    return u


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    primes = [p for p in range(5, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    forms = {}
    for cl, (shape, _) in SHAPES.items():
        sw = sum(a * e for a, e in shape.items())
        assert sw == 24
        forms[cl] = shape
    forms.update(EXTRA)
    coef = {}
    for name, shape in forms.items():
        assert sum(a * e for a, e in shape.items()) == 24
        b = eta_quotient(shape, PMAX)
        coef[name] = {p: b[p - 1] for p in primes}      # q^1 * prod: a(m) = b[m-1]
    chars = {'1': lambda p: 1, '(-1/p)': lambda p: legendre(-1, p), '(2/p)': lambda p: legendre(2, p),
             '(-2/p)': lambda p: legendre(-2, p), '(3/p)': lambda p: legendre(3, p), '(-3/p)': lambda p: legendre(-3, p),
             '(6/p)': lambda p: legendre(6, p), '(-6/p)': lambda p: legendre(-6, p)}
    here = os.path.dirname(os.path.abspath(__file__))
    classes = json.load(open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'halfperiod_twin_classes_20261018.json')))
    seqnames = []
    for c in classes:
        for nm in [c['base']] + [m['seq'] for m in c['members']]:
            if nm not in seqnames:
                seqnames.append(nm)
    report = {}
    for nm in seqnames:
        u = make_seq(nm)
        val = {p: u((p - 1) // 2) for p in primes}
        found = []
        modp = []
        for fname in forms:
            for cname, ch in chars.items():
                for i in range(-6, 7):
                    for j in range(-3, 4):
                        ok2 = True
                        ok1 = True
                        for p in primes:
                            n = (p - 1) // 2
                            pp = p * p
                            s = pow(2, i * n, pp) if i >= 0 else pow(pow(2, -i * n, pp), -1, pp)
                            s = s * (pow(3, j * n, pp) if j >= 0 else pow(pow(3, -j * n, pp), -1, pp)) % pp
                            rhs = ch(p) * s * coef[fname][p]
                            if (val[p] - rhs) % pp:
                                ok2 = False
                                if (val[p] - rhs) % p:
                                    ok1 = False
                                    break
                        if ok2:
                            found.append((fname, cname, i, j))
                        elif ok1 and i == 0 and j == 0:
                            modp.append((fname, cname))
        report[nm] = dict(modp2=found, modp_only=modp)
        print('u = sum %s' % nm)
        if found:
            for fname, cname, i, j in found:
                print('    == %s * 2^(%d n) 3^(%d n) * a_%s(p)  (mod p^2), all 5 <= p <= %d' % (cname, i, j, fname, PMAX))
        else:
            print('    no mod p^2 match; mod p only (eps, form):', modp[:12])
        sys.stdout.flush()
    json.dump(report, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'halfperiod_twin_class_forms_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
