#!/usr/bin/env python3
"""Connectors between the translated coefficient copies of clock_duality's cubic recipe (mac-mini-2026-10-18).

golden_contact_circuit_20261017i.md section 4: the coefficient factors C + bM and C + cM (b, c with denominator 19) are
seven disjoint translated copies of the 10-point, 15-edge graph C; no unit edge joins different copies, so each copy is
coloured independently. The recorded next target is a coefficient context that connects these copies.

A k-step connector from x (in copy m) to y (in copy m') is a unit-distance path x = z_0, ..., z_k = y with every z_i in
the field F = Q(zeta15, sqrt(-11), w) (degree 48). For k = 2 the middle point is z = x + (t/2)(1 + s), t = y - x, with
s^2 = 1 - 4/N(t), N(t) = t * conj(t); it exists in F iff N(t)(N(t) - 4) is a square in F and |t| <= 2 in the geometric
embedding. Squareness is screened by quadratic characters at degree-1 places (one non-residue proves non-square) and
survivors are confirmed by an exact square root (sympy factorisation of x^2 - D over F is avoided: we confirm by finding
s in F with s^2 = D through the residue-CRT-free route of testing candidate s = (2z - t)/t for z on both unit circles is
not available, so survivors are reported as 'square at all tested places' and then certified separately).
k = 3 connectors: t = q1 + (two-step), q1 ranging over a unit set U0 (unit differences of C, and their products with
the 30 roots of unity in F), with the k = 2 test applied to t - q1.
Usage: PYTHONPATH=first_arrival_20261018/shim python3 first_arrival_20261018/copy_connectors.py [NPLACES]
"""
import os, sys, json, time, cmath, itertools
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
from golden_forcing_search_20261017 import Geo, recipe
from golden_vsquare_unit_census_20261017 import places, red

T0 = time.time()
def log(*a): print('[%5.0f s]' % (time.time() - T0), *a, flush=True)

def main():
    NPL = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    geo = Geo(); f = geo.f
    A, C = recipe(geo, 'base')
    d, db = geo.d, geo.bar(geo.d)
    w = f.gen(2); h = f.add(w, f.scalar(3))
    u = f.mul(f.sub(h, geo.j), f.inv(f.add(h, geo.j)))
    m2 = f.scale(f.add(f.scalar(20), f.scale(geo.j, 11)), Q(-1, 19))
    m1 = f.scale(f.sub(f.scalar(20), f.scale(geo.j, 11)), Q(1, 19))
    db2 = f.mul(db, db)
    cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2), f.one), f.mul(db2, m1)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
    cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1), db2))
    log('field degree', f.n, '; cubic coefficients a = conj(d), b, c built')
    # numeric embedding: find (zeta, j, w) images making every listed unit edge of C have modulus 1
    wr = [r for r in __import__('numpy').roots([1, 2, -4, -4]) if abs(r.imag) < 1e-12 and 1 < r.real < 2][0].real
    EC = geo.unit_edges(C)
    best = None
    for kz in range(1, 15):
        if kz % 3 == 0 or kz % 5 == 0: continue
        for js in (1, -1):
            zc = cmath.exp(2j * cmath.pi * kz / 15); jc = js * 1j * 11 ** 0.5
            def ev(x):
                return sum(float(c) * zc ** e[0] * jc ** e[1] * wr ** e[2] for c, e in zip(x, f.basis) if c)
            err = max(abs(abs(ev(f.sub(C[b2], C[a2]))) - 1) for a2, b2 in EC)
            if best is None or err < best[0]: best = (err, kz, js)
    err, kz, js = best
    zc = cmath.exp(2j * cmath.pi * kz / 15); jc = js * 1j * 11 ** 0.5
    def ev(x):
        return sum(float(c) * zc ** e[0] * jc ** e[1] * wr ** e[2] for c, e in zip(x, f.basis) if c)
    log('geometric embedding zeta -> e(%d/15), j -> %s i sqrt(11), w -> %.6f; unit-edge error %.2e' % (kz, '+' if js > 0 else '-', wr, err))
    pl = places(f, NPL)
    log('places:', [p for p, _, _ in pl])
    def nonsquare_somewhere(x):
        for p, emb, embc in pl:
            if any(c and c.denominator % p == 0 for c in x): continue      # not integral at this place: skip it
            r = red(x, p, emb)
            if r == 0: continue
            if pow(r, (p - 1) // 2, p) == p - 1: return True
        return False
    four = f.scalar(4)
    out = {}
    for name, co in (('b', cb), ('c', cc)):
        Nb = f.mul(co, geo.bar(co)); nb = abs(ev(co)) ** 2
        X = f.mul(Nb, f.sub(Nb, four))
        log('coefficient %s: |%s|^2 = %.6f; N(t)(N(t)-4) for t = %s*q is %s' % (name, name, nb, name,
            'NON-square (screened)' if nonsquare_somewhere(X) else 'square at all %d places' % len(pl)))
        # two-step connectors between points of different copies c + co*m and c' + co*m'
        cands = []; tested = 0
        for (i1, m), (i2, mm) in itertools.combinations(enumerate(geo.M), 2):
            for a1, x in enumerate(C):
                for a2, y in enumerate(C):
                    t = f.sub(f.add(y, f.mul(co, mm)), f.add(x, f.mul(co, m)))
                    nt = abs(ev(t)) ** 2
                    if nt > 4 + 1e-9 or nt < 1e-12: continue
                    tested += 1
                    Nt = f.mul(t, geo.bar(t))
                    if not nonsquare_somewhere(f.mul(Nt, f.sub(Nt, four))):
                        cands.append((i1, i2, a1, a2, nt))
        log('  two-step: %d point pairs across copies with |t| <= 2 tested; %d survive the place screen' % (tested, len(cands)))
        # three-step connectors through a first unit step q1 from a unit set
        U0 = sorted({f.sub(C[b2], C[a2]) for a2, b2 in EC} | {f.sub(C[a2], C[b2]) for a2, b2 in EC})
        if os.environ.get('ROOTS_OF_UNITY'):
            z15 = f.gen(0); roots = [f.one]
            for _ in range(14): roots.append(f.mul(roots[-1], z15))
            roots = roots + [f.scale(r, -1) for r in roots]          # the 30 roots of unity of F
            U0 = sorted({f.mul(r, q) for r in roots for q in U0 + [f.one]})
        if os.environ.get('U_DIRECTIONS'):
            # 19-adic unit directions: zeta^l * u^k, k = +-1, +-2 (u is a unit vector whose cubic over K has denominator 19)
            z15 = f.gen(0); roots = [f.one]
            for _ in range(14): roots.append(f.mul(roots[-1], z15))
            roots = roots + [f.scale(r, -1) for r in roots]
            ui = f.inv(u); pw = [u, f.mul(u, u), ui, f.mul(ui, ui)]
            U0 = sorted({f.mul(r, q) for r in roots for q in pw})
        if os.environ.get('K_DIRECTIONS'):
            # K-directions conj(gamma)/gamma with poles at the primes above 19 (incl. the h = 7 type), times roots of unity
            from kdirs import kdirections
            z15 = f.gen(0); roots = [f.one]
            for _ in range(14): roots.append(f.mul(roots[-1], z15))
            roots = roots + [f.scale(r, -1) for r in roots]
            kd = kdirections(geo); log('  K-directions:', [lab for lab, _ in kd])
            pw = [dl for _, dl in kd] + [f.inv(dl) for _, dl in kd]
            U0 = sorted({f.mul(r, q) for r in roots for q in pw})
        cands3 = []; tested3 = 0
        for (i1, m), (i2, mm) in itertools.combinations(enumerate(geo.M), 2):
            for a1, x in enumerate(C):
                for a2, y in enumerate(C):
                    t0 = f.sub(f.add(y, f.mul(co, mm)), f.add(x, f.mul(co, m)))
                    if abs(ev(t0)) > 3 + 1e-9: continue
                    for q in U0:
                        t = f.sub(t0, q); nt = abs(ev(t)) ** 2
                        if nt > 4 + 1e-9 or nt < 1e-12: continue
                        tested3 += 1
                        Nt = f.mul(t, geo.bar(t))
                        if not nonsquare_somewhere(f.mul(Nt, f.sub(Nt, four))):
                            cands3.append((i1, i2, a1, a2, U0.index(q), nt))
        log('  three-step via %d unit directions: %d tests; %d survive the place screen' % (len(U0), tested3, len(cands3)))
        out[name] = dict(norm=nb, two_step_tested=tested, two_step_survivors=cands, three_step_tested=tested3,
                         three_step_survivors=cands3[:200], three_step_survivor_count=len(cands3))
    suf = '_kdirs' if os.environ.get('K_DIRECTIONS') else ('_udirs' if os.environ.get('U_DIRECTIONS') else ('_roots' if os.environ.get('ROOTS_OF_UNITY') else ''))
    json.dump(out, open(os.path.join(HERE, 'copy_connectors%s.json' % suf), 'w'), indent=1, default=str)

if __name__ == '__main__':
    main()
