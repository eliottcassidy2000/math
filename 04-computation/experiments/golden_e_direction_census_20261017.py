#!/usr/bin/env python3
"""Cross-edge census for clock_duality's degree-192 direction e = Zs/h on the golden source (S41).

clock_duality_oct15 (golden_square_contexts_20261017h.md, G2) adjoined h with h^2 = R = 4 - N to E = F(s), s^2 = D, and
introduced the unit direction e = Z s/h (e^2 = -Z/conj Z), so that sigma-mates become unit-distance pairs. Their stated
next target is a unit edge joining different port and e-coefficient fibres. This script asks, for the coefficient sets
of the golden source, which differences

    delta = dA + u d1 + v d2 + e d3,   dA in D(A) u {0},  d1, d2 in D(C) u {0},  d3 in D(C) u D(L),  d3 != 0,

are unit vectors in E' = F(s, h). Here L is their 13-point attachment set (the odd 13-cycle). Elements of E' are
written c0 + c1 s + c2 h + c3 sh over F, with conj: s -> -s, h -> h (h is real at the intended embedding). The norm
x conj(x) is computed componentwise. A multi-prime modular sieve is followed by an exact check of every survivor. If
only pure e-axis units survive (dA = d1 = d2 = 0), the e-direction adds no cross edges to the linear golden source.
Usage: python golden_e_direction_census_20261017.py
Session opus-2026-10-17-S41.
"""
import os, sys, json, time
from fractions import Fraction as Q
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from escaped_ports_20261017d import make, need
from golden_vsquare_unit_census_20261017 import places, red


@njit(cache=True)
def sieve(RA, R1, R2, R3, Dn, Rn, p):
    """Residue rows: [c0, c1, c3, conj c0, conj c1, conj c3] (c2 = 0 for every candidate)."""
    out = []
    DR = Dn * Rn % p
    for i3 in range(R3.shape[0]):
        c3 = R3[i3, 0]; c3b = R3[i3, 1]
        for i2 in range(R2.shape[0]):
            for i1 in range(R1.shape[0]):
                b0 = (R1[i1, 0] + R2[i2, 0]) % p; b0b = (R1[i1, 1] + R2[i2, 1]) % p
                c1 = R2[i2, 2]; c1b = R2[i2, 3]
                # h-component: c1 conj c3 + c3 conj c1 = 0 (independent of dA)
                if (c1 * c3b + c3 * c1b) % p != 0:
                    continue
                for iA in range(RA.shape[0]):
                    c0 = (b0 + RA[iA, 0]) % p; c0b = (b0b + RA[iA, 1]) % p
                    if (c3 * c0b - c0 * c3b) % p != 0:
                        continue
                    if (c1 * c0b - c0 * c1b) % p != 0:
                        continue
                    if (c0 * c0b - Dn * c1 % p * c1b - DR * c3 % p * c3b) % p != 1:
                        continue
                    out.append((iA, i1, i2, i3))
    return out


def main():
    t0 = time.time()
    f, z, j, alpha, ports, units, eta, conj, escape, unused = make('sqrt3')
    bar = lambda a: f.apply(a, conj)
    u = units[0]; d = f.add(f.one, alpha); Z = f.add(d, u)
    N = f.mul(Z, bar(Z)); D = f.sub(f.one, f.scale(f.inv(N), 4)); R = f.sub(f.scalar(4), N)
    need(bar(D) == D and bar(R) == R, 'D and R real')
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
    a5 = ports[5]                                    # (3 - j sqrt5)/8, |a5| = |a5 - 3/4| = 1
    Lp = [f.zero, f.one, f.sub(f.one, a5), f.scalar(Q(1, 4)), f.scalar(Q(5, 4)), f.sub(f.scalar(Q(5, 4)), a5), f.scalar(Q(1, 2))]
    L = list({x: 0 for x in Lp + [f.scale(x, -1) for x in Lp]}.keys())
    need(len(A) == 54 and len(C) == 10 and len(L) == 13, 'sizes (54, 10, 13)')
    DA = sorted({f.sub(a, b) for a in A for b in A} - {f.zero})
    DC = sorted({f.sub(a, b) for a in C for b in C} - {f.zero})
    D3 = sorted(({f.sub(a, b) for a in C for b in C} | {f.sub(a, b) for a in L for b in L}) - {f.zero})
    Rinv = f.inv(R); ZR = f.mul(Z, Rinv)
    print('field ready; |D(A)| = %d, |D(C)| = %d, |D(C) u D(L)| = %d [%.0f s]' % (len(DA), len(DC), len(D3), time.time() - t0), flush=True)
    # exact check that e = (0, 0, 0, Z/R) is a unit with e^2 = -Z/conj Z
    def mul4(x, y):
        a0, a1, a2, a3 = x; b0, b1, b2, b3 = y
        DR = f.mul(D, R)
        return (f.add(f.add(f.mul(a0, b0), f.mul(D, f.mul(a1, b1))), f.add(f.mul(R, f.mul(a2, b2)), f.mul(DR, f.mul(a3, b3)))),
                f.add(f.add(f.mul(a0, b1), f.mul(a1, b0)), f.mul(R, f.add(f.mul(a2, b3), f.mul(a3, b2)))),
                f.add(f.add(f.mul(a0, b2), f.mul(a2, b0)), f.mul(D, f.add(f.mul(a1, b3), f.mul(a3, b1)))),
                f.add(f.add(f.mul(a0, b3), f.mul(a3, b0)), f.add(f.mul(a1, b2), f.mul(a2, b1))))
    def conj4(x):
        return (bar(x[0]), f.scale(bar(x[1]), -1), bar(x[2]), f.scale(bar(x[3]), -1))
    e = (f.zero, f.zero, f.zero, ZR)
    need(mul4(e, conj4(e)) == (f.one, f.zero, f.zero, f.zero), 'e is a unit')
    need(mul4(e, e) == (f.scale(f.mul(Z, f.inv(bar(Z))), -1), f.zero, f.zero, f.zero), 'e^2 = -Z/conj Z')
    print('e verified [%.0f s]' % (time.time() - t0), flush=True)
    # building blocks: (c0, c1, c3) contributions
    BA = [(x, f.zero, f.zero) for x in [f.zero] + DA]
    B1 = [(f.mul(u, x), f.zero, f.zero) for x in [f.zero] + DC]
    B2 = [(f.scale(f.mul(Z, x), Q(-1, 2)), f.scale(f.mul(Z, x), Q(1, 2)), f.zero) for x in [f.zero] + DC]
    B3 = [(f.zero, f.zero, f.mul(ZR, x)) for x in D3]
    pl = places(f, 3)
    survivors = None
    for p, emb, embc in pl:
        def rows(B, which):
            out = []
            for c0, c1, c3 in B:
                if which == 'A1':
                    out.append([red(c0, p, emb), red(c0, p, embc)])
                elif which == '2':
                    out.append([red(c0, p, emb), red(c0, p, embc), red(c1, p, emb), red(c1, p, embc)])
                else:
                    out.append([red(c3, p, emb), red(c3, p, embc)])
            return np.array(out, dtype=np.int64)
        RA, R1, R2, R3 = rows(BA, 'A1'), rows(B1, 'A1'), rows(B2, '2'), rows(B3, '3')
        Dn, Rn = red(D, p, emb), red(R, p, emb)
        if survivors is None:
            survivors = set(sieve(RA, R1, R2, R3, Dn, Rn, p))
        else:
            keep = set()
            for (iA, i1, i2, i3) in survivors:
                c0 = (RA[iA, 0] + R1[i1, 0] + R2[i2, 0]) % p; c0b = (RA[iA, 1] + R1[i1, 1] + R2[i2, 1]) % p
                c1, c1b = R2[i2, 2], R2[i2, 3]; c3, c3b = R3[i3, 0], R3[i3, 1]
                ok = ((c1 * c3b + c3 * c1b) % p == 0 and (c3 * c0b - c0 * c3b) % p == 0 and (c1 * c0b - c0 * c1b) % p == 0
                      and (c0 * c0b - Dn * c1 * c1b - Dn * Rn * c3 * c3b) % p == 1)
                if ok:
                    keep.add((iA, i1, i2, i3))
            survivors = keep
        print('prime %d: %d survivors [%.0f s]' % (p, len(survivors), time.time() - t0), flush=True)
    found = []
    for (iA, i1, i2, i3) in sorted(survivors):
        c0 = f.add(f.add(BA[iA][0], B1[i1][0]), B2[i2][0]); c1 = B2[i2][1]; c3 = B3[i3][2]
        x = (c0, c1, f.zero, c3)
        if mul4(x, conj4(x)) == (f.one, f.zero, f.zero, f.zero):
            found.append((iA, i1, i2, i3))
    pure = [t for t in found if t[:3] == (0, 0, 0)]
    unitD3 = [i for i, x in enumerate(D3) if f.mul(x, bar(x)) == f.one]
    mixed = [t for t in found if t[:3] != (0, 0, 0)]
    print('exact units with a nonzero e part: %d (pure e-axis: %d; unit differences in D(C) u D(L): %d; mixed: %d) [%.0f s]'
          % (len(found), len(pure), len(unitD3), len(mixed), time.time() - t0), flush=True)
    for t in mixed[:40]:
        print('  mixed unit: dA %d, d1 %d, d2 %d, d3 %d' % t)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(sizes=dict(DA=len(DA), DC=len(DC), D3=len(D3)), primes=[p for p, _, _ in pl], units=len(found),
                   pure_axis=len(pure), unit_differences_D3=len(unitD3), mixed=mixed),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'golden_e_direction_census_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
