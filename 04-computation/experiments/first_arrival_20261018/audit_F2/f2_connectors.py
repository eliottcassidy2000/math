#!/usr/bin/env python3
"""Audit F2, check (5): independent re-screen of the copy-connector census (copy_connectors*.py / *.json).

Own field arithmetic (f2_field.py), own degree-1 places, own residue code.  Logic certified here:
  * two-step criterion: for t = y - x != 0 a middle point z in F with z conj z = 1 = (t - z) conj(t - z) exists iff
    X = N(t)(N(t) - 4) is a square in F and |t| <= 2 in one (hence every) complex embedding (F is CM);
  * a quadratic non-residue red_P(X) != 0 at a degree-1 prime P (p not dividing disc(O)) proves X is not a square.
Steps:
  A. reconstruct the session's direction lists (sorted exactly as the session sorted them) and check that the
     recorded |t|^2 of every survivor in copy_connectors{,_roots,_udirs}.json is reproduced;
  B. certify each recorded survivor as a non-square at my own places (all embeddings at each prime);
  C. recount ALL two-step and three-step tests (numeric prefilter in the session's embedding) and certify every
     test by a non-residue, so that the census claims are re-derived, not just the survivors.
"""
import sys, os, json, cmath, math, time, itertools
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2_field import E, ONE, J, W, EXP, N, ZPOW, build, roots_of_unity, conj, norm_cc

HERE = os.path.dirname(os.path.abspath(__file__))
SESSION = os.path.dirname(HERE)
T0 = time.time()


def log(*a):
    print('[%5.0f s]' % (time.time() - T0), *a, flush=True)


def primes_from(start):
    """Generator of primes p > start with p = 1 (mod 15) (trial division; p is small)."""
    p = start
    while True:
        p += 1
        if p % 15 != 1 or any(p % r == 0 for r in range(2, int(p ** 0.5) + 1)):
            continue
        yield p


def places_all(nprimes, start=1000):
    """All embeddings Z[zeta,j,w] -> F_p for primes p = 1 mod 15 where -11 is a square and the cubic has a root."""
    pl = []
    for p in primes_from(start):
        zs = [x for x in range(2, p) if pow(x, 15, p) == 1 and pow(x, 5, p) != 1 and pow(x, 3, p) != 1]
        js = [x for x in range(1, p) if (x * x + 11) % p == 0]
        ws = [x for x in range(p) if (x ** 3 + 2 * x * x - 4 * x - 4) % p == 0]
        if not (js and ws):
            continue
        embs = []
        for z0 in zs:
            for j0 in js:
                for w0 in ws:
                    embs.append((z0, j0, w0))
        pl.append((p, embs))
        if len(pl) == nprimes:
            break
    return pl


_TAB = {}


def _tab(p, emb):
    key = (p, emb)
    if key not in _TAB:
        z0, j0, w0 = emb
        _TAB[key] = [pow(z0, a, p) * pow(j0, b, p) * pow(w0, c, p) % p for (a, b, c) in (EXP[k] for k in range(N))]
    return _TAB[key]


def red(x, p, emb):
    assert x.d % p, 'denominator divisible by the place prime'
    tab = _tab(p, emb)
    s = 0
    for k, v in enumerate(x.n):
        if v:
            s += v * tab[k]
    return s * pow(x.d, -1, p) % p


def red_conj(x, p, emb):
    z0, j0, w0 = emb
    return red(x, p, (pow(z0, 14, p), (-j0) % p, w0))


def main():
    O = build()
    b, c, u, M, C = O['b'], O['c'], O['u'], O['M'], O['C']
    roots = roots_of_unity()
    # unit edges of C and the session's direction lists
    EC = [(i, k) for i in range(len(C)) for k in range(i + 1, len(C)) if norm_cc(C[k] - C[i]) == ONE]
    U0 = sorted({C[k] - C[i] for i, k in EC} | {C[i] - C[k] for i, k in EC}, key=lambda e: e.fracs())
    D90 = sorted({r * q for r in roots for q in U0 + [ONE]}, key=lambda e: e.fracs())
    ui = u.inv()
    D120 = sorted({r * q for r in roots for q in (u, u * u, ui, ui * ui)}, key=lambda e: e.fracs())
    log('C: %d points, %d unit edges; |U0| = %d, |D90| = %d, |D120| = %d' % (len(C), len(EC), len(U0), len(D90), len(D120)))
    log('all directions are unit vectors:', all(norm_cc(q) == ONE for q in U0 + D90 + D120))
    # numeric embedding used by the session: zeta -> e(2/15), j -> -i sqrt 11, w -> real root in (1, 2)
    wr = 1.5
    for _ in range(100):
        wr -= (wr ** 3 + 2 * wr ** 2 - 4 * wr - 4) / (3 * wr ** 2 + 4 * wr - 4)
    zc = cmath.exp(2j * cmath.pi * 2 / 15); jc = -1j * math.sqrt(11)
    pw = {(a, bb, cc_): zc ** a * jc ** bb * wr ** cc_ for (a, bb, cc_) in EXP.values()}

    def ev(x):
        return sum(v * pw[EXP[k]] for k, v in enumerate(x.n) if v) / x.d
    log('embedding check: |u| = %.15f, |rho| = %.15f, w = %.9f' % (abs(ev(u)), abs(ev(O['rho'])), wr))
    PL = places_all(40)
    log('places: %d primes from %d to %d, %d embeddings in all' % (len(PL), PL[0][0], PL[-1][0], sum(len(e) for _, e in PL)))

    def certificate(X, places=PL, want_all=False):
        """First (p, emb) where red(X) is a nonzero non-residue; None if none."""
        found = []
        for p, embs in places:
            for emb in embs:
                r = red(X, p, emb)
                if r and pow(r, (p - 1) // 2, p) == p - 1:
                    if not want_all:
                        return (p, emb)
                    found.append((p, emb)); break
        return found if want_all else None

    report = {}
    # ---------------- A + B: the recorded survivors
    for fname, D in (('copy_connectors.json', U0), ('copy_connectors_roots.json', D90), ('copy_connectors_udirs.json', D120)):
        Jn = json.load(open(os.path.join(SESSION, fname)))
        for name, co in (('b', b), ('c', c)):
            for (i1, i2, a1, a2, qi, nt) in Jn[name]['three_step_survivors']:
                x = C[a1] + co * M[i1]; y = C[a2] + co * M[i2]
                t = (y - x) - D[qi]
                Nt = norm_cc(t); X = Nt * (Nt - 4)
                my_nt = abs(ev(t)) ** 2
                cert = certificate(X, want_all=True)
                key = '%s %s copies (%d,%d) points (%d,%d) q#%d' % (fname, name, i1, i2, a1, a2, qi)
                ps = sorted({p for p, _ in cert})
                log('%s: |t|^2 recorded %.6f, recomputed %.6f ; non-residue at primes %s' % (key, nt, my_nt, ps[:8]))
                report[key] = dict(nt_recorded=nt, nt_recomputed=my_nt, nonresidue_primes=ps)
    # identity of the 22-direction survivors with the 90-direction ones
    q8, q20 = U0[8], U0[20]
    log('U0[8] == D90[39]: %s ; U0[20] == D90[88]: %s' % (q8 == D90[39], q20 == D90[88]))

    # ---------------- C: full recount
    for label, D in (('U0', U0), ('D90', D90), ('D120', D120)):
        evD = [ev(q) for q in D]
        for name, co in (('b', b), ('c', c)):
            pts = {(i, a): C[a] + co * M[i] for i in range(7) for a in range(len(C))}
            evp = {k: ev(v) for k, v in pts.items()}
            # residues of points and directions at the first 12 primes (all embeddings)
            P12 = PL[:12]
            resP = {}; resD = {}
            for p, embs in P12:
                for emb in embs:
                    resP[(p, emb)] = {k: (red(v, p, emb), red_conj(v, p, emb)) for k, v in pts.items()}
                    resD[(p, emb)] = [(red(q, p, emb), red_conj(q, p, emb)) for q in D]
            two_tested = two_unc = 0; three_tested = three_unc = 0; uncert = []
            near = []
            for (i1, i2) in itertools.combinations(range(7), 2):
                for a1 in range(len(C)):
                    for a2 in range(len(C)):
                        t0 = evp[(i2, a2)] - evp[(i1, a1)]
                        n0 = abs(t0) ** 2
                        if label == 'U0':
                            # two-step tests (independent of the direction list)
                            if abs(n0 - 4) < 1e-6 or n0 < 1e-6: near.append(('2', n0))
                            if n0 <= 4 + 1e-9 and n0 > 1e-12:
                                two_tested += 1
                                ok = False
                                for (p, emb) in resP:
                                    ry = resP[(p, emb)][(i2, a2)]; rx = resP[(p, emb)][(i1, a1)]
                                    nr = (ry[0] - rx[0]) * (ry[1] - rx[1]) % p
                                    X = nr * (nr - 4) % p
                                    if X and pow(X, (p - 1) // 2, p) == p - 1:
                                        ok = True; break
                                if not ok:
                                    two_unc += 1; uncert.append(('2', i1, i2, a1, a2))
                        if abs(t0) > 3 + 1e-9:
                            continue
                        for qi, eq in enumerate(evD):
                            n1 = abs(t0 - eq) ** 2
                            if abs(n1 - 4) < 1e-6 or n1 < 1e-6: near.append(('3', n1))
                            if n1 > 4 + 1e-9 or n1 < 1e-12:
                                continue
                            three_tested += 1
                            ok = False
                            for (p, emb) in resP:
                                ry = resP[(p, emb)][(i2, a2)]; rx = resP[(p, emb)][(i1, a1)]; rq = resD[(p, emb)][qi]
                                nr = (ry[0] - rx[0] - rq[0]) * (ry[1] - rx[1] - rq[1]) % p
                                X = nr * (nr - 4) % p
                                if X and pow(X, (p - 1) // 2, p) == p - 1:
                                    ok = True; break
                            if not ok:
                                three_unc += 1; uncert.append(('3', i1, i2, a1, a2, qi))
            # anything not certified at the first 12 primes: try all 40 primes exactly
            still = []
            for tup in uncert:
                if tup[0] == '2':
                    _, i1, i2, a1, a2 = tup; t = pts[(i2, a2)] - pts[(i1, a1)]
                else:
                    _, i1, i2, a1, a2, qi = tup; t = pts[(i2, a2)] - pts[(i1, a1)] - D[qi]
                Nt = norm_cc(t); X = Nt * (Nt - 4)
                if certificate(X) is None:
                    still.append(tup)
            log('%s %s: two-step tests %s (uncertified at 12 primes: %s) ; three-step tests %d (uncertified at 12 primes: %d) ; '
                'uncertified after 40 primes: %d ; near-boundary numeric cases: %d'
                % (label, name, two_tested if label == 'U0' else '-', two_unc if label == 'U0' else '-', three_tested, three_unc,
                   len(still), len(near)))
            report['%s %s' % (label, name)] = dict(two_step_tested=two_tested if label == 'U0' else None, three_step_tested=three_tested,
                                                   uncertified_after_12_primes=len(uncert), uncertified_after_40=len(still),
                                                   near_boundary=len(near))
    json.dump(report, open(os.path.join(HERE, 'f2_connectors.json'), 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
