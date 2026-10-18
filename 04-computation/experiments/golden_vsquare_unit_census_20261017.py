#!/usr/bin/env python3
"""Unit-difference census for a v^2 coordinate on clock_duality's shared golden source (S41).

clock_duality_oct15 (shared_golden_contact_20261017g.md) proposed adding v^2 coordinates to G = A + uC + vC, since
v^2 = -Zv - Z/conj(Z) introduces a uv coefficient outside the classified linear space. Before any large graph is built,
this script asks which differences

    delta = dA + u d1 + v d2 + v^2 d3,     dA in D(A) u {0},  d1, d2, d3 in D(C) u {0},  d3 != 0,

are unit vectors, D(X) being the difference set of X. If the only units are the pure v^2-axis units (dA = d1 = d2 = 0,
d3 a unit), then G' = A + uC + vC + v^2 C has no new mixed edges, every port pattern extends additively
(golden_additive_colourings_20261017.py), and the candidate cannot force at any size. The test is a multi-prime
modular sieve on the norm equation x conj(x) = 1 in K(w)(s), s^2 = D, followed by an exact characteristic-zero check of
every survivor with clock_duality's field arithmetic.
Usage: python golden_vsquare_unit_census_20261017.py
Session opus-2026-10-17-S41.
"""
import os, sys, json, time
from fractions import Fraction as Q
from itertools import product
import numpy as np
import sympy as sp
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from escaped_ports_20261017d import make, need


def places(f, count, start=1000):
    """Good primes p = 1 (mod 15) with sqrt(-11) and a root of the cubic; one embedding each."""
    out = []
    for p in sp.primerange(start, 10 ** 6):
        if p % 15 != 1:
            continue
        zs = [x for x in range(2, p) if pow(x, 15, p) == 1 and pow(x, 3, p) != 1 and pow(x, 5, p) != 1]
        js = sp.sqrt_mod(-11, p, all_roots=True)
        ws = [x for x in range(p) if (x ** 3 + 2 * x * x - 4 * x - 4) % p == 0]
        if not (zs and js and ws):
            continue
        z0, j0, w0 = zs[0], int(js[0]), ws[0]
        for poly, r in zip(f.polys, (z0, j0, w0)):
            need(sum(c * pow(r, k, p) for k, c in enumerate(poly)) % p == 0, 'defining relation at the place')
        emb = np.array([pow(z0, e[0], p) * pow(j0, e[1], p) * pow(w0, e[2], p) % p for e in f.basis], dtype=object)
        zc = pow(z0, 14, p)
        embc = np.array([pow(zc, e[0], p) * pow((-j0) % p, e[1], p) * pow(w0, e[2], p) % p for e in f.basis], dtype=object)
        out.append((p, emb, embc))
        if len(out) == count:
            break
    return out


def red(a, p, emb):
    s = 0
    for c, e in zip(a, emb):
        if c:
            need(c.denominator % p != 0, 'good prime for this coefficient')
            s += c.numerator * pow(c.denominator, -1, p) * e
    return s % p


@njit(cache=True)
def sieve(RA, R1, R2, R3, Dn, p):
    """Return index quadruples (iA, i1, i2, i3) whose summed residues satisfy the unit norm equations mod p."""
    out = []
    for i3 in range(R3.shape[0]):
        for i2 in range(R2.shape[0]):
            for i1 in range(R1.shape[0]):
                b0 = (R1[i1, 0] + R2[i2, 0] + R3[i3, 0]) % p
                b1 = (R1[i1, 1] + R2[i2, 1] + R3[i3, 1]) % p
                b2 = (R1[i1, 2] + R2[i2, 2] + R3[i3, 2]) % p
                b3 = (R1[i1, 3] + R2[i2, 3] + R3[i3, 3]) % p
                for iA in range(RA.shape[0]):
                    x0 = (b0 + RA[iA, 0]) % p; x1 = b1; x2 = (b2 + RA[iA, 2]) % p; x3 = b3
                    if (x0 * x2 + Dn * x1 % p * x3) % p == 1 and (x0 * x3 + x1 * x2) % p == 0:
                        out.append((iA, i1, i2, i3))
    return out


def main():
    t0 = time.time()
    f, z, j, alpha, ports, units, eta, conj, escape, unused = make('sqrt3')
    bar = lambda a: f.apply(a, conj)
    u = units[0]; d = f.add(f.one, alpha); Z = f.add(d, u)
    N = f.mul(Z, bar(Z)); D = f.sub(f.one, f.scale(f.inv(N), 4)); need(bar(D) == D, 'D real')
    rho = f.scale(f.add(f.scalar(5), j), Q(1, 6)); baseD = f.add(f.one, eta)
    M = [f.zero, f.one, eta, baseD, rho, f.mul(rho, eta), f.mul(rho, baseD)]
    A = list(ports); seen = set(A)
    for port in ports:
        for m in M:
            t = f.add(port, m)
            if t not in seen:
                seen.add(t); A.append(t)
    C = list(M)
    for p0 in ports[:5]:
        if p0 not in C:
            C.append(p0)
    need(len(A) == 54 and len(C) == 10, 'clock_duality source sizes (54, 10)')
    DA = sorted({f.sub(a, b) for a in A for b in A} - {f.zero})
    DC = sorted({f.sub(a, b) for a in C for b in C} - {f.zero})
    print('field ready; |D(A)| = %d, |D(C)| = %d [%.0f s]' % (len(DA), len(DC), time.time() - t0), flush=True)
    ZbarZ = f.mul(Z, f.inv(bar(Z)))
    v2a = f.sub(f.scale(f.mul(Z, Z), Q(1, 2)), ZbarZ); v2b = f.scale(f.mul(Z, Z), Q(-1, 2))
    # elements as pairs (a, b) meaning a + b s
    EA = [(x, f.zero) for x in [f.zero] + DA]
    E1 = [(f.mul(u, x), f.zero) for x in [f.zero] + DC]
    E2 = [(f.scale(f.mul(Z, x), Q(-1, 2)), f.scale(f.mul(Z, x), Q(1, 2))) for x in [f.zero] + DC]
    E3 = [(f.mul(v2a, x), f.mul(v2b, x)) for x in DC]
    # exact check of the v^2 identity
    v = (f.scale(Z, Q(-1, 2)), f.scale(Z, Q(1, 2)))
    vv = (f.add(f.mul(v[0], v[0]), f.mul(D, f.mul(v[1], v[1]))), f.scale(f.mul(v[0], v[1]), 2))
    need(vv == (v2a, v2b), 'v^2 = (Z^2/2 - Z/conj Z) - (Z^2/2) s')
    print('elements ready [%.0f s]' % (time.time() - t0), flush=True)
    pl = places(f, 3)
    survivors = None
    for p, emb, embc in pl:
        def R(E):
            return np.array([[red(a, p, emb), red(b, p, emb), red(a, p, embc), (-red(b, p, embc)) % p] for a, b in E], dtype=np.int64)
        RA, R1, R2, R3 = R(EA), R(E1), R(E2), R(E3)
        Dn = red(D, p, emb)
        if survivors is None:
            hits = sieve(RA, R1, R2, R3, Dn, p)
            survivors = set(hits)
        else:
            keep = set()
            for (iA, i1, i2, i3) in survivors:
                x = [(RA[iA, k] + R1[i1, k] + R2[i2, k] + R3[i3, k]) % p for k in range(4)]
                if (x[0] * x[2] + Dn * x[1] * x[3]) % p == 1 and (x[0] * x[3] + x[1] * x[2]) % p == 0:
                    keep.add((iA, i1, i2, i3))
            survivors = keep
        print('prime %d: %d survivors [%.0f s]' % (p, len(survivors), time.time() - t0), flush=True)
    units_found = []
    for (iA, i1, i2, i3) in sorted(survivors):
        a = f.add(f.add(EA[iA][0], E1[i1][0]), f.add(E2[i2][0], E3[i3][0]))
        b = f.add(f.add(EA[iA][1], E1[i1][1]), f.add(E2[i2][1], E3[i3][1]))
        na = f.sub(f.mul(a, bar(a)), f.mul(D, f.mul(b, bar(b)))); nb = f.sub(f.mul(b, bar(a)), f.mul(a, bar(b)))
        if na == f.one and nb == f.zero:
            units_found.append((iA, i1, i2, i3))
    pure = [t for t in units_found if t[0] == 0 and t[1] == 0 and t[2] == 0]
    mixed = [t for t in units_found if not (t[0] == 0 and t[1] == 0 and t[2] == 0)]
    unitDC = [i for i, x in enumerate(DC) if f.mul(x, bar(x)) == f.one]
    print('exact units with a nonzero v^2 part: %d (pure v^2-axis: %d; |unit differences of C| = %d; mixed: %d) [%.0f s]'
          % (len(units_found), len(pure), len(unitDC), len(mixed), time.time() - t0), flush=True)
    for t in mixed[:40]:
        print('  mixed unit: dA index %d, d1 index %d, d2 index %d, d3 index %d' % t)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(sizes=dict(DA=len(DA), DC=len(DC)), primes=[p for p, _, _ in pl], survivors=len(survivors),
                   units=len(units_found), pure_axis=len(pure), unit_differences_of_C=len(unitDC), mixed=mixed),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'golden_vsquare_unit_census_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
