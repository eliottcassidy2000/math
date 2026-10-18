"""S40 audit: are all contact rotations U = (e - d1)/d2 of mixed_spindle_rotations_20261017.py in K = Q(zeta15, sqrt-11)?

The session's script never tests field membership. This driver recomputes the rotation set from scratch (no cache),
using clock_duality's exact field arithmetic (escaped_ports_20261017d.make), and reports for every U the w-components
(basis exponents of w, the cubic generator of K(w)/K) and whether U is a root of unity, plus |U| = 1 and the provenance
identity d1 + U d2 = e.
Usage: python audit_rotations_field.py
"""
import sys, os, itertools, json, time
from fractions import Fraction as Q
EXP = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, EXP)
from escaped_ports_20261017d import make


def main():
    t0 = time.time()
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
    u, v = units
    inv = lambda g: f.apply(g, conj)
    gens = [rho, u, v]
    words = [f.one]
    for g in gens:
        words += [g, inv(g)]
    for g1, g2 in itertools.combinations(gens, 2):
        for h1 in (g1, inv(g1)):
            for h2 in (g2, inv(g2)):
                words.append(f.mul(h1, h2))
    E = [f.mul(r, w) for w in words for r in roots]
    Nd = {d: norm(d) for d in D}
    Us = {}
    for e in E:
        for d1 in D:
            ned = norm(f.sub(e, d1))
            for d2 in D:
                if ned == Nd[d2]:
                    U = f.mul(f.sub(e, d1), f.inv(d2))
                    Us.setdefault(U, (e, d1, d2))
    # w is generator index 2; basis exponent tuples (ez, ej, ew)
    wpos = [i for i, b in enumerate(f.basis) if b[2] != 0]
    rootset = set(roots)
    rows = []
    for U, (e, d1, d2) in Us.items():
        inK = all(U[i] == 0 for i in wpos)
        rows.append(dict(inK=inK, unit=(norm(U) == f.one), provenance_ok=(f.add(d1, f.mul(U, d2)) == e), root_of_unity=(U in rootset)))
    print('unit directions %d, spindle differences %d, distinct rotations %d (%.0f s)' % (len(E), len(D), len(Us), time.time() - t0))
    print('in K (no w-component): %d of %d; |U| = 1: %d; provenance identity: %d; roots of unity: %d'
          % (sum(r['inK'] for r in rows), len(rows), sum(r['unit'] for r in rows), sum(r['provenance_ok'] for r in rows),
             sum(r['root_of_unity'] for r in rows)))
    json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'audit_rotations_field.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
