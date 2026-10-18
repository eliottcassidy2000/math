#!/usr/bin/env python3
"""Referee check B3 (second S41 audit): independent float census of the e and product-direction unit differences.

Method (independent of the session's modular sieve and of escaped_ports_20261017d.Field):
- A, C, L, D(A), D(C), D(C) u D(L) are built from their definitions with the referee's own exact K arithmetic
  (audit2_field.py); difference sets are deduplicated exactly.
- Every difference dA + u d1 + v d2 + Y d3 (dA in D(A) u {0}, d1, d2 in D(C) u {0}, d3 in D(C) u D(L), d3 != 0) is
  evaluated in complex double precision at the intended embedding (zeta15 = exp(2 pi i/15), j = +i sqrt 11, w real in
  (1, 2)); |x|^2 = 1 at this embedding is equivalent to x conj(x) = 1 in E' because conj is complex conjugation there.
- Candidates with ||x|^2 - 1| < TOL (TOL = 1e-8, far above the rounding error ~1e-13 for |x| <= 20) are identified
  EXACTLY in K: pure Y-axis (dA = d1 = d2 = 0, d3 conj(d3) = 1), or the claimed families
  uv: (dA, d1, d2, d3) = q (0, 1, 1, conj d);  e^2: q (d, 1, 0, conj d), with q conj(q) = 1 (exact in K).
  Those families are exact units by the identities u + v + conj(d) uv = uv conj(d + u + v) and
  d + u + conj(d) e^2 = Z conj(u)/conj(Z) (proved by hand in the report; checked numerically here to 1e-15).
- The smallest deviation among non-candidates is reported to show the separation.
Usage: python b3_census_numeric.py [directions...]   (default: e ue ve uv ee)
"""
import sys, time, json
import numpy as np
from numba import njit
import audit2_field as F

T0 = time.time()


@njit(cache=True)
def scan(XA, X12, X3, tol):
    out = []
    mind = 1e9
    tiers = np.zeros(5, dtype=np.int64)   # counts with dev < 1e-12, 1e-10, 1e-8, 1e-6, 1e-4
    for i3 in range(X3.shape[0]):
        y3 = X3[i3]
        for i12 in range(X12.shape[0]):
            y = X12[i12] + y3
            yr = y.real; yi = y.imag
            for iA in range(XA.shape[0]):
                zr = XA[iA].real + yr; zi = XA[iA].imag + yi
                dev = abs(zr * zr + zi * zi - 1.0)
                if dev < 1e-4:
                    if dev < 1e-12:
                        tiers[0] += 1
                    if dev < 1e-10:
                        tiers[1] += 1
                    if dev < 1e-8:
                        tiers[2] += 1
                    if dev < 1e-6:
                        tiers[3] += 1
                    tiers[4] += 1
                if dev < tol:
                    out.append((iA, i12, i3))
                elif dev < mind:
                    mind = dev
    return out, mind, tiers


def main():
    dirs = sys.argv[1:] or ['e', 'ue', 've', 'uv', 'ee']
    S = F.golden_sets()
    A, C, L, d = S['A'], S['C'], S['L'], S['d']
    assert len(A) == 54 and len(C) == 10 and len(L) == 13
    DA = F.diffs(A); DC = F.diffs(C)
    D3 = F.uniq(DC + F.diffs(L))
    print('sizes: |A| %d |C| %d |L| %d |D(A)| %d |D(C)| %d |D(C) u D(L)| %d; candidates per direction %d x %d x %d x %d = %d'
          % (len(A), len(C), len(L), len(DA), len(DC), len(D3), len(DA) + 1, len(DC) + 1, len(DC) + 1, len(D3),
             (len(DA) + 1) * (len(DC) + 1) ** 2 * len(D3)), flush=True)
    G = F.golden_numbers()
    u, v, e, Z = G['u'], G['v'], G['e'], G['Z']
    dn = F.to_complex(d)
    print('embedding: w = %.15f, N = %.15f, R = %.15f, D = %.15f, |u| %.15f |v| %.15f |e| %.15f |d+u+v| %.15f'
          % (G['w'], G['N'], G['R'], G['D'], abs(u), abs(v), abs(e), abs(dn + u + v)))
    print('identity checks: |u + v + conj(d) uv| - 1 = %.2e; |uv conj(d+u+v) - (u + v + conj(d) uv)| = %.2e; '
          '|d + u + conj(d) e^2 - Z conj(u)/conj(Z)| = %.2e; e^2 + Z/conj(Z) = %.2e'
          % (abs(u + v + dn.conjugate() * u * v) - 1, abs(u * v * (dn + u + v).conjugate() - (u + v + dn.conjugate() * u * v)),
             abs(dn + u + dn.conjugate() * e * e - Z * u.conjugate() / Z.conjugate()), abs(e * e + Z / Z.conjugate())))
    XA = np.array([0j] + [F.to_complex(x) for x in DA])
    XC = [0j] + [F.to_complex(x) for x in DC]
    X12 = np.array([u * a + v * b for a in XC for b in XC])
    X3c = np.array([F.to_complex(x) for x in D3])
    Ys = {'e': e, 'ue': u * e, 've': v * e, 'uv': u * v, 'ee': e * e}
    DAx = [F.K(0)] + DA; DCx = [F.K(0)] + DC
    dbar = F.conj(d)
    unitD3 = [i for i, x in enumerate(D3) if F.is_unit(x)]
    print('exact unit elements of D(C) u D(L): %d' % len(unitD3), flush=True)
    summary = {}
    for name in dirs:
        Y = Ys[name]
        X3 = Y * X3c
        t = time.time()
        cand, mind, tiers = scan(XA, X12, X3, 1e-8)
        kinds = {'pure': 0, 'family': 0, 'other': 0}
        others = []
        fam_vecs = []
        for iA, i12, i3 in cand:
            i1, i2 = divmod(i12, len(XC))
            dA, d1, d2, d3 = DAx[iA], DCx[i1], DCx[i2], D3[i3]
            if F.is_zero(dA) and F.is_zero(d1) and F.is_zero(d2) and F.is_unit(d3):
                kinds['pure'] += 1
                continue
            ok = False
            if name == 'uv':
                ok = F.is_zero(dA) and d1 == d2 and d3 == F.mul(dbar, d1) and F.is_unit(d1)
            elif name == 'ee':
                ok = dA == F.mul(d, d1) and F.is_zero(d2) and d3 == F.mul(dbar, d1) and F.is_unit(d1)
            if ok:
                kinds['family'] += 1
                fam_vecs.append((iA, i1, i2, i3))
            else:
                kinds['other'] += 1
                others.append((iA, i1, i2, i3))
        print('direction %-3s: %d candidates with ||x|^2-1| < 1e-8 -> pure Y-axis %d, claimed family %d, other %d; '
              'tiers(<1e-12,1e-10,1e-8,1e-6,1e-4) = %s; smallest other deviation %.3e [%.1f s]'
              % (name, len(cand), kinds['pure'], kinds['family'], kinds['other'], list(tiers), mind, time.time() - t), flush=True)
        for o in others[:20]:
            print('   OTHER candidate (needs exact check):', o)
        summary[name] = dict(candidates=len(cand), pure=kinds['pure'], family=kinds['family'], other=kinds['other'],
                             others=others, family_indices=fam_vecs, tiers=[int(x) for x in tiers], min_other=mind)
    json.dump(dict(sizes=[len(DA), len(DC), len(D3)], unitD3=len(unitD3), summary=summary),
              open('b3_census_numeric.json', 'w'), indent=1)
    print('total %.0f s' % (time.time() - T0))


if __name__ == '__main__':
    main()
