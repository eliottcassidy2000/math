#!/usr/bin/env python3
"""Mixed-coordinate unit edges between Moser spindles in the sqrt3 escape F_s = Q(zeta15, sqrt-11, w) (S40).

clock_duality_oct15 (escape_renewal_idoneal_20261017d.md, section 5) found that the triple spindle sum M + uM + vM with
the carrier's cubic and mixed-golden units is exactly the Cartesian cube of the spindle: no unit edge changes two
coordinates. Here the rotation is chosen to create such edges. A unit vector e of F_s and spindle differences d1, d2 with
|e - d1| = |d2| give the unit-modulus U = (e - d1)/d2, and then d1 + U d2 = e is a unit edge of M + U M changing both
coordinates (a four-bar / codegree coincidence, THM-4663 (2)). The unit vectors e range over a finite ball:
+-zeta15^k times words of length <= 2 in rho, u, v and their inverses.

For each rotation U this script reports: whether U is a root of unity, the number of distinct points of the block
M + U M, its single-coordinate and mixed unit edges, and the block's chromatic number (SAT). Optionally (--ports) it
attaches the block at the eight ports (translated, and rotated by the powers of zeta5) and tests all 570 boundary
patterns by SAT, as in escaped_ports_20261017d.py: exact field arithmetic for every edge.
Usage: python mixed_spindle_rotations_20261017.py [--ports] [--maxU K]
Result (S40): 570 unit directions and 34 spindle differences give 34 rotations U, all of which lie in the base field
K = Q(zeta15, sqrt-11) (no w-component), as clock_duality's degree lemma (marked_contacts_root_current_20261017e.md (5))
forces: K(w)/K is cubic and a contact needs [K(U):K] <= 2. The blocks M + U M with strictly mixed edges (4 or 2 of them)
are four-chromatic, and attached at the eight ports (with the zeta5 rotations) they block none of the 570 boundary
patterns (output: 05-knowledge/results/mixed_spindle_rotations_20261017.out).
Session opus-2026-10-17-S40.
"""
import sys, os, tempfile, itertools, argparse, json
from fractions import Fraction as Q
from pysat.solvers import Solver
sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else __file__.rsplit('/', 1)[0])
from escaped_ports_20261017d import make
from golden_boundary_20261017c import PORT_EDGES, patterns


def chromatic(n, edges, kmax=7):
    for k in range(1, kmax + 1):
        with Solver(name='g3') as s:
            for a in range(n):
                s.add_clause([k * a + c + 1 for c in range(k)])
            for a, b in edges:
                for c in range(k):
                    s.add_clause([-(k * a + c + 1), -(k * b + c + 1)])
            if a_sym := (n > 0):
                s.add_clause([1])               # colour of vertex 0 fixed (symmetry)
            if s.solve():
                return k
    return None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--ports', action='store_true'); ap.add_argument('--maxU', type=int, default=0)
    ap.add_argument('--units', action='store_true', help='also rotate the attached blocks by the units u, v')
    ap.add_argument('--cache', default=os.path.join(tempfile.gettempdir(), 'mixed_spindle_rotations_cache.pkl'))
    args = ap.parse_args()
    f, z, j, alpha, ports, units, eta, conj, escape, reduce = make('sqrt3')
    norm = lambda x: f.mul(x, f.apply(x, conj))
    rho = f.scale(f.add(f.scalar(5), j), Q(1, 6))
    a1 = f.add(f.one, eta)
    M = [f.zero, f.one, eta, a1, rho, f.mul(rho, eta), f.mul(rho, a1)]
    D = {}
    for i in range(7):
        for k in range(7):
            if i != k:
                D.setdefault(f.sub(M[i], M[k]), (i, k))
    roots = []
    x = f.one
    for k in range(15):
        roots.append(x); roots.append(f.scale(x, -1)); x = f.mul(x, z)
    rootset = set(roots)
    u, v = units
    inv = lambda g: f.apply(g, conj)
    gens = [('rho', rho), ('u', u), ('v', v)]
    words = [((), f.one)]
    for name, g in gens:
        words += [(((name, 1),), g), (((name, -1),), inv(g))]
    for (n1, g1), (n2, g2) in itertools.combinations(gens, 2):
        for e1, h1 in ((1, g1), (-1, inv(g1))):
            for e2, h2 in ((1, g2), (-1, inv(g2))):
                words.append((((n1, e1), (n2, e2)), f.mul(h1, h2)))
    E = [((k, w), f.mul(r, wv)) for w, wv in words for k, r in enumerate(roots)]
    import pickle
    Nd = {d: norm(d) for d in D}
    Us = {}
    if os.path.exists(args.cache):
        Us = pickle.load(open(args.cache, 'rb'))
    else:
        for en, e in E:
            for d1, l1 in D.items():
                ned = norm(f.sub(e, d1))
                for d2, l2 in D.items():
                    if ned == Nd[d2]:
                        U = f.mul(f.sub(e, d1), f.inv(d2))
                        Us.setdefault(U, []).append((en, l1, l2))
        pickle.dump(Us, open(args.cache, 'wb'))
    print('unit directions %d, spindle differences %d, rotations U %d' % (len(E), len(D), len(Us)))
    out = []
    for idx, (U, prov) in enumerate(Us.items()):
        pts, seen, addr = [], {}, []
        for a in range(7):
            for b in range(7):
                p = f.add(M[a], f.mul(U, M[b]))
                if p not in seen:
                    seen[p] = len(pts); pts.append(p); addr.append({(a, b)})
                else:
                    addr[seen[p]].add((a, b))
        edges, mixed = [], 0
        for i in range(len(pts)):
            for k in range(i + 1, len(pts)):
                if norm(f.sub(pts[i], pts[k])) == f.one:
                    edges.append((i, k))
                    if all(sum(x != y for x, y in zip(s1, s2)) > 1 for s1 in addr[i] for s2 in addr[k]):
                        mixed += 1
        chi = chromatic(len(pts), edges)
        row = dict(index=idx, root_of_unity=U in rootset, points=len(pts), edges=len(edges), mixed_edges=mixed, chromatic=chi,
                   witness=str(prov[0]))
        out.append(row)
        print(row, flush=True)
    if args.ports:
        cand = [r for r in out if not r['root_of_unity']]
        cand.sort(key=lambda r: (-r['chromatic'], -r['mixed_edges']))
        Ulist = list(Us.keys())
        for r in cand[:max(args.maxU, 1)]:
            U = Ulist[r['index']]
            seed = sorted({f.add(M[a], f.mul(U, M[b])) for a in range(7) for b in range(7)})
            multipliers = [f.power(alpha, b) for b in range(5)]
            if args.units:
                multipliers = [f.mul(m, w) for m in multipliers for w in [f.one] + list(units)]
            points = list(ports); seenp = set(points)
            for p in ports:
                for m in multipliers:
                    for q in seed:
                        w = f.add(p, f.mul(m, q))
                        if w not in seenp:
                            seenp.add(w); points.append(w)
            import numpy as np
            res = np.array([reduce(w) for w in points], dtype=np.int64)
            P = 1319
            edges = []
            for a in range(len(points)):
                diff = res[a + 1:] - res[a]
                for off in np.flatnonzero((diff[:, 0] ** 2 + diff[:, 1] ** 2) % P == 1):
                    b = a + 1 + int(off)
                    if norm(f.sub(points[a], points[b])) == f.one:
                        edges.append((a, b))
            assert [e for e in edges if e[1] < 8] == sorted(tuple(sorted(e)) for e in PORT_EDGES)
            pats = patterns(8, PORT_EDGES)
            blocked = []
            with Solver(name='g3') as s:
                n = len(points)
                for a in range(n):
                    s.add_clause([6 * a + c + 1 for c in range(6)])
                for a, b in edges:
                    for c in range(6):
                        s.add_clause([-(6 * a + c + 1), -(6 * b + c + 1)])
                for pi, pat in enumerate(pats):
                    if not s.solve(assumptions=[6 * a + c + 1 for a, c in enumerate(pat)]):
                        blocked.append(pi)
            print('ports + block %d: vertices %d, edges %d, blocked boundary patterns %d of %d %s'
                  % (r['index'], len(points), len(edges), len(blocked), len(pats), blocked[:10]), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'mixed_spindle_rotations_20261017.json'), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
