#!/usr/bin/env python3
"""Referee check B4 (second S41 audit): modular sieve at a DIFFERENT prime set, with the referee's own place evaluation.

For a degree-one place phi of F (p = 1 mod 15, sqrt(-11) and a root w0 of the cubic in F_p), phi(conj y) is the place
phi_c: zeta15 -> z0^-1, j -> -j0, w -> w0.  Every element of E' = F(s, h) is c0 + c1 s + c2 h + c3 sh, with
  1 -> (1,0,0,0), u -> (u,0,0,0), v -> (-Z/2, Z/2, 0, 0), e -> (0,0,0,Z/R), ue -> (0,0,0,uZ/R),
  ve -> (0,0,Z^2 D/(2R), -Z^2/(2R)), uv -> (-uZ/2, uZ/2, 0, 0), e^2 -> (Z^2 D/R, 0, 0, 0) = -Z/conj Z.
x conj(x) = 1 is equivalent to the four F-equations
  c0 c0' - D c1 c1' + R c2 c2' - DR c3 c3' = 1,   c1 c0' - c0 c1' - R c2 c3' + R c3 c2' = 0,
  c0 c2' + c2 c0' - D (c1 c3' + c3 c1') = 0,        c3 c0' - c0 c3' + c1 c2' - c2 c1' = 0   (c' = conj c),
whose images under phi must hold mod p.  All sieve survivors are then identified exactly in K (audit2_field.py) and
any survivor that is neither a pure Y-axis unit nor in the claimed family is evaluated at 60 digits.
Usage: python b4_census_modular.py [directions]   (default: e ee uv)
"""
import sys, time, json
import numpy as np
import sympy as sp
from numba import njit
import audit2_field as F

T0 = time.time()


def find_places(count, start):
    out = []
    for p in sp.primerange(start, 10 ** 6):
        if p % 15 != 1:
            continue
        js = sp.sqrt_mod(-11 % p, p, all_roots=True)
        if not js:
            continue
        ws = [x for x in range(p) if (x ** 3 + 2 * x * x - 4 * x - 4) % p == 0]
        if not ws:
            continue
        zs = [x for x in range(2, p) if pow(x, 15, p) == 1 and pow(x, 3, p) != 1 and pow(x, 5, p) != 1]
        # take the LAST root of each list (the session took the first), to vary the place as well as the prime
        out.append((p, zs[-1], int(js[-1]), ws[-1]))
        if len(out) == count:
            return out


def fvals(p, z0, j0, w0):
    inv = lambda x: pow(x % p, -1, p)
    u = (w0 + 3 - j0) * inv(w0 + 3 + j0) % p
    ub = inv(u)
    a0 = pow(z0, 3, p)
    d = (1 + a0) % p; db = (1 + inv(a0)) % p
    Z = (d + u) % p; Zb = (db + ub) % p
    N = Z * Zb % p; D = (1 - 4 * inv(N)) % p; R = (4 - N) % p
    i2, iR = inv(2), inv(R)
    vec = {
        'one': [1, 0, 0, 0], 'u': [u, 0, 0, 0], 'v': [-Z * i2, Z * i2, 0, 0],
        'e': [0, 0, 0, Z * iR], 'ue': [0, 0, 0, u * Z * iR], 've': [0, 0, Z * Z * D * i2 * iR, -Z * Z * i2 * iR],
        'uv': [-u * Z * i2, u * Z * i2, 0, 0], 'ee': [Z * Z * D * iR, 0, 0, 0]}
    vecb = {
        'one': [1, 0, 0, 0], 'u': [ub, 0, 0, 0], 'v': [-Zb * i2, Zb * i2, 0, 0],
        'e': [0, 0, 0, Zb * iR], 'ue': [0, 0, 0, ub * Zb * iR], 've': [0, 0, Zb * Zb * D * i2 * iR, -Zb * Zb * i2 * iR],
        'uv': [-ub * Zb * i2, ub * Zb * i2, 0, 0], 'ee': [Zb * Zb * D * iR, 0, 0, 0]}
    vec = {k: [x % p for x in vv] for k, vv in vec.items()}
    vecb = {k: [x % p for x in vv] for k, vv in vecb.items()}
    # sanity: e^2 = -Z/conj Z at this place
    assert vec['ee'][0] == (-Z * inv(Zb)) % p
    return D, R, vec, vecb


@njit(cache=True)
def check(c, cb, Dn, Rn, p):
    c0, c1, c2, c3 = c[0], c[1], c[2], c[3]
    e0, e1, e2, e3 = cb[0], cb[1], cb[2], cb[3]
    if (c3 * e0 - c0 * e3 + c1 * e2 - c2 * e1) % p != 0:
        return False
    if (c0 * e2 + c2 * e0 - Dn * ((c1 * e3 + c3 * e1) % p)) % p != 0:
        return False
    if (c1 * e0 - c0 * e1 - Rn * ((c2 * e3) % p) + Rn * ((c3 * e2) % p)) % p != 0:
        return False
    if (c0 * e0 - Dn * ((c1 * e1) % p) + Rn * ((c2 * e2) % p) - ((Dn * Rn) % p) * ((c3 * e3) % p)) % p != 1:
        return False
    return True


@njit(cache=True)
def full_scan(rA, rAc, r1, r1c, r2, r2c, r3, r3c, U, Ub, V, Vb, Y, Yb, Dn, Rn, p):
    out = []
    c = np.zeros(4, dtype=np.int64); cb = np.zeros(4, dtype=np.int64)
    for i3 in range(r3.shape[0]):
        for i1 in range(r1.shape[0]):
            for i2 in range(r2.shape[0]):
                b = np.zeros(4, dtype=np.int64); bb = np.zeros(4, dtype=np.int64)
                for k in range(4):
                    b[k] = (r1[i1] * U[k] + r2[i2] * V[k] + r3[i3] * Y[k]) % p
                    bb[k] = (r1c[i1] * Ub[k] + r2c[i2] * Vb[k] + r3c[i3] * Yb[k]) % p
                for iA in range(rA.shape[0]):
                    c[0] = (b[0] + rA[iA]) % p; c[1] = b[1]; c[2] = b[2]; c[3] = b[3]
                    cb[0] = (bb[0] + rAc[iA]) % p; cb[1] = bb[1]; cb[2] = bb[2]; cb[3] = bb[3]
                    if check(c, cb, Dn, Rn, p):
                        out.append((iA, i1, i2, i3))
    return out


def main():
    dirs = sys.argv[1:] or ['e', 'ee', 'uv']
    S = F.golden_sets()
    d = S['d']; dbar = F.conj(d)
    DA = [F.K(0)] + F.diffs(S['A']); DC = [F.K(0)] + F.diffs(S['C'])
    D3 = F.uniq(F.diffs(S['C']) + F.diffs(S['L']))
    places = find_places(3, 5000)
    print('places (p, zeta15, j, w):', places, flush=True)
    res = {}
    for name in dirs:
        surv = None
        for (p, z0, j0, w0) in places:
            zc = pow(z0, -1, p); jc = (-j0) % p
            Dn, Rn, vec, vecb = fvals(p, z0, j0, w0)
            r = lambda L: np.array([F.to_mod(x, p, z0, j0) for x in L], dtype=np.int64)
            rc = lambda L: np.array([F.to_mod(x, p, zc, jc) for x in L], dtype=np.int64)
            rA, rAc, rC, rCc, r3, r3c = r(DA), rc(DA), r(DC), rc(DC), r(D3), rc(D3)
            A4 = lambda key: np.array(vec[key], dtype=np.int64)
            B4 = lambda key: np.array(vecb[key], dtype=np.int64)
            t = time.time()
            if surv is None:
                surv = set(full_scan(rA, rAc, rC, rCc, rC, rCc, r3, r3c, A4('u'), B4('u'), A4('v'), B4('v'), A4(name), B4(name), Dn, Rn, p))
            else:
                keep = set()
                for (iA, i1, i2, i3) in surv:
                    hit = full_scan(rA[iA:iA + 1], rAc[iA:iA + 1], rC[i1:i1 + 1], rCc[i1:i1 + 1], rC[i2:i2 + 1], rCc[i2:i2 + 1],
                                    r3[i3:i3 + 1], r3c[i3:i3 + 1], A4('u'), B4('u'), A4('v'), B4('v'), A4(name), B4(name), Dn, Rn, p)
                    if hit:
                        keep.add((iA, i1, i2, i3))
                surv = keep
            print('  %s: prime %d -> %d survivors [%.1f s]' % (name, p, len(surv), time.time() - t), flush=True)
        kinds = {'pure': 0, 'family': 0, 'other': []}
        for (iA, i1, i2, i3) in sorted(surv):
            dA, d1, d2, d3 = DA[iA], DC[i1], DC[i2], D3[i3]
            if F.is_zero(dA) and F.is_zero(d1) and F.is_zero(d2) and F.is_unit(d3):
                kinds['pure'] += 1
            elif name == 'uv' and F.is_zero(dA) and d1 == d2 and d3 == F.mul(dbar, d1) and F.is_unit(d1):
                kinds['family'] += 1
            elif name == 'ee' and dA == F.mul(d, d1) and F.is_zero(d2) and d3 == F.mul(dbar, d1) and F.is_unit(d1):
                kinds['family'] += 1
            else:
                kinds['other'].append((iA, i1, i2, i3))
        print('direction %s: survivors %d = pure Y-axis %d + claimed family %d + other %d %s' % (
            name, len(surv), kinds['pure'], kinds['family'], len(kinds['other']), kinds['other'][:10]), flush=True)
        res[name] = dict(survivors=len(surv), pure=kinds['pure'], family=kinds['family'], other=kinds['other'])
    json.dump(dict(places=places, results=res), open('b4_census_modular.json', 'w'), indent=1)
    print('total %.0f s' % (time.time() - T0))


if __name__ == '__main__':
    main()
