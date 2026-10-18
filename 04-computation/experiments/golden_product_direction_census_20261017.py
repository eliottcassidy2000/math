#!/usr/bin/env python3
"""Cross-edge census for product directions (u e, v e, e^2, u v) over clock_duality's golden source in F(s, h) (S41).

THM-4668 shows that the degree-192 direction e is Cartesian on K + uK + vK + eK. The remaining natural candidates for
cross edges are product coordinates. For a direction Y this script enumerates all differences

    delta = dA + u d1 + v d2 + Y d3,   dA in D(A) u {0},  d1, d2 in D(C) u {0},  d3 in D(C) u D(L),  d3 != 0,

writes them as c0 + c1 s + c2 h + c3 sh over F (s^2 = D, h^2 = R, conj: s -> -s, h -> h), sieves the four norm
equations x conj(x) = (1, 0, 0, 0) at three good primes, and checks every survivor exactly. Any survivor that is not a
pure Y-axis unit is a genuine cross edge between the Y-layer and the old coordinates.
Usage: python golden_product_direction_census_20261017.py ue|ve|ee|uv
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
    """Rows hold 8 residues: c0..c3 then conj c0..conj c3 (plain F-conjugates)."""
    out = []
    DR = Dn * Rn % p
    for i3 in range(R3.shape[0]):
        for i2 in range(R2.shape[0]):
            for i1 in range(R1.shape[0]):
                b = np.empty(8, dtype=np.int64)
                for k in range(8):
                    b[k] = (R1[i1, k] + R2[i2, k] + R3[i3, k]) % p
                for iA in range(RA.shape[0]):
                    c0 = (b[0] + RA[iA, 0]) % p; c1 = (b[1] + RA[iA, 1]) % p; c2 = (b[2] + RA[iA, 2]) % p; c3 = (b[3] + RA[iA, 3]) % p
                    e0 = (b[4] + RA[iA, 4]) % p; e1 = (b[5] + RA[iA, 5]) % p; e2 = (b[6] + RA[iA, 6]) % p; e3 = (b[7] + RA[iA, 7]) % p
                    sh = (-c0 * e3 + c3 * e0 + c1 * e2 - c2 * e1) % p
                    if sh != 0:
                        continue
                    hh = (c0 * e2 + c2 * e0 - Dn * ((c1 * e3 + c3 * e1) % p)) % p
                    if hh != 0:
                        continue
                    ss = (-c0 * e1 + c1 * e0 + Rn * ((-c2 * e3 + c3 * e2) % p)) % p
                    if ss != 0:
                        continue
                    one = (c0 * e0 - Dn * (c1 * e1 % p) + Rn * (c2 * e2 % p) - DR * (c3 * e3 % p)) % p
                    if one == 1:
                        out.append((iA, i1, i2, i3))
    return out


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else 've'
    t0 = time.time()
    f, z, j, alpha, ports, units, eta, conj, escape, unused = make('sqrt3')
    bar = lambda a: f.apply(a, conj)
    u = units[0]; d = f.add(f.one, alpha); Z = f.add(d, u)
    N = f.mul(Z, bar(Z)); D = f.sub(f.one, f.scale(f.inv(N), 4)); R = f.sub(f.scalar(4), N)
    DR = f.mul(D, R)
    def mul4(x, y):
        a0, a1, a2, a3 = x; b0, b1, b2, b3 = y
        return (f.add(f.add(f.mul(a0, b0), f.mul(D, f.mul(a1, b1))), f.add(f.mul(R, f.mul(a2, b2)), f.mul(DR, f.mul(a3, b3)))),
                f.add(f.add(f.mul(a0, b1), f.mul(a1, b0)), f.mul(R, f.add(f.mul(a2, b3), f.mul(a3, b2)))),
                f.add(f.add(f.mul(a0, b2), f.mul(a2, b0)), f.mul(D, f.add(f.mul(a1, b3), f.mul(a3, b1)))),
                f.add(f.add(f.mul(a0, b3), f.mul(a3, b0)), f.add(f.mul(a1, b2), f.mul(a2, b1))))
    def conj4(x):
        return (bar(x[0]), f.scale(bar(x[1]), -1), bar(x[2]), f.scale(bar(x[3]), -1))
    def sc(x, k):
        return tuple(f.mul(c, k) for c in x)
    one4 = (f.one, f.zero, f.zero, f.zero)
    U4 = (u, f.zero, f.zero, f.zero)
    V4 = (f.scale(Z, Q(-1, 2)), f.scale(Z, Q(1, 2)), f.zero, f.zero)
    E4 = (f.zero, f.zero, f.zero, f.mul(Z, f.inv(R)))
    need(mul4(E4, conj4(E4)) == one4 and mul4(V4, conj4(V4)) == one4, 'units v, e')
    Y = {'ue': mul4(U4, E4), 've': mul4(V4, E4), 'ee': mul4(E4, E4), 'uv': mul4(U4, V4)}[which]
    need(mul4(Y, conj4(Y)) == one4, 'Y is a unit')
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
    a5 = ports[5]
    Lp = [f.zero, f.one, f.sub(f.one, a5), f.scalar(Q(1, 4)), f.scalar(Q(5, 4)), f.sub(f.scalar(Q(5, 4)), a5), f.scalar(Q(1, 2))]
    L = list({x: 0 for x in Lp + [f.scale(x, -1) for x in Lp]}.keys())
    DA = sorted({f.sub(a, b) for a in A for b in A} - {f.zero})
    DC = sorted({f.sub(a, b) for a in C for b in C} - {f.zero})
    D3 = sorted(({f.sub(a, b) for a in C for b in C} | {f.sub(a, b) for a in L for b in L}) - {f.zero})
    BA = [sc(one4, x) for x in [f.zero] + DA]
    B1 = [sc(U4, x) for x in [f.zero] + DC]
    B2 = [sc(V4, x) for x in [f.zero] + DC]
    B3 = [sc(Y, x) for x in D3]
    print('direction %s ready: |D(A)| = %d, |D(C)| = %d, |D3| = %d [%.0f s]' % (which, len(DA), len(DC), len(D3), time.time() - t0), flush=True)
    pl = places(f, 3)
    survivors = None
    for p, emb, embc in pl:
        def rows(B):
            return np.array([[red(c, p, emb) for c in x] + [red(c, p, embc) for c in x] for x in B], dtype=np.int64)
        RA, R1, R2, R3 = rows(BA), rows(B1), rows(B2), rows(B3)
        Dn, Rn = red(D, p, emb), red(R, p, emb)
        if survivors is None:
            survivors = set(sieve(RA, R1, R2, R3, Dn, Rn, p))
        else:
            keep = set()
            for t in survivors:
                if sieve(RA[t[0]:t[0] + 1], R1[t[1]:t[1] + 1], R2[t[2]:t[2] + 1], R3[t[3]:t[3] + 1], Dn, Rn, p):
                    keep.add(t)
            survivors = keep
        print('prime %d: %d survivors [%.0f s]' % (p, len(survivors), time.time() - t0), flush=True)
    found = []
    for (iA, i1, i2, i3) in sorted(survivors):
        x = tuple(f.add(f.add(BA[iA][k], B1[i1][k]), f.add(B2[i2][k], B3[i3][k])) for k in range(4))
        if mul4(x, conj4(x)) == one4:
            found.append((iA, i1, i2, i3))
    pure = [t for t in found if t[:3] == (0, 0, 0)]
    mixed = [t for t in found if t[:3] != (0, 0, 0)]
    print('direction %s: exact units with a nonzero Y part: %d (pure Y-axis %d, mixed %d) [%.0f s]'
          % (which, len(found), len(pure), len(mixed), time.time() - t0), flush=True)
    for t in mixed[:40]:
        print('  mixed unit: dA %d, d1 %d, d2 %d, d3 %d' % t)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(direction=which, units=len(found), pure_axis=len(pure), mixed=mixed),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'golden_product_direction_census_%s_20261017.json' % which), 'w'), indent=1)


if __name__ == '__main__':
    main()
