"""P5(c) of PREREGISTRATION.md, executed separately (disclosed in the outcomes).

In confirm.py the ambient controls were guarded by `n <= 8`, but every P5(b) core has order >= 9
(k >= 3 plus at least three children of order >= 2), so P5(c) never ran. This script runs the
pre-registered prediction unchanged: for a core embedded as a module in an exterior of 2-3 extra
vertices, with both ports copying the core's exterior incidences, the half-contrast mod 2 restricted
to the core's kernel labels equals that of the bare core; hence the predicted functional values hold.
Writes p5c.json and prints a summary.
"""
import json, os, random
from confirm import (subst, TT, conv, random_tournament, rand_child, endpoint_parity, coeffs,
                     restrict_half, Lfun, HERE)

def main():
    rng = random.Random(4864)
    C = {}
    def bump(key, ok):
        C.setdefault(key, [0, 0]); C[key][0] += 1
        if not ok: C[key][1] += 1
    done = 0
    while done < 4:
        m = 4; Q = random_tournament(m, rng); j = rng.randrange(m); k = 3
        sizes = [k if i == j else 2 for i in range(m)]          # core order 9
        kids = [TT(k) if i == j else rand_child(2, rng) for i in range(m)]
        W, bl = subst(Q, kids); n = len(W)
        pj = endpoint_parity(Q)[j]
        lab = {}
        for i, b in enumerate(bl):
            for v in b: lab[v] = i
        lab[bl[j][0]] = 'u'
        for v in bl[j][1:]: lab[v] = 'beta'
        verts = [i for i in range(m) if i != j] + ['u', 'beta']
        def arcZ(a, b):
            qa = j if a in ('u', 'beta') else a; qb = j if b in ('u', 'beta') else b
            if qa == qb: return (a, b) == ('u', 'beta')
            return Q[qa][qb] == 1
        r = rng.choice((2, 3)); R = random_tournament(r + 1, rng)   # vertex 0 of R is the core
        for eps in (0, 1):
            kids2 = [((conv(TT(k)) if eps else TT(k)) if i == j else rand_child(2, rng)) for i in range(m)]
            Wp, _ = subst(conv(Q), kids2)
            var = list(range(n))
            CW, CWp = coeffs([(W, [0] * n, [0] * n, var), (Wp, [0] * n, [0] * n, var)])
            d, _ = restrict_half(CW, CWp, lab)
            TA, blA = subst(R, [W] + [[[0]]] * r); TB, _ = subst(R, [Wp] + [[[0]]] * r)
            N = len(TA); core = blA[0]
            a0 = [0] * N; c0 = [0] * N
            for e in range(1, r + 1):
                for v in blA[e]:
                    a0[v] = R[e][0]; c0[v] = R[e][0]     # v -> port iff v -> core (ports copy the core)
            CA, CB = coeffs([(TA, a0, c0, core), (TB, a0, c0, core)])
            dA, _ = restrict_half(CA, CB, lab)
            bump('P5c ambient restricted mod-4 contrast == bare core', dA == d)
            if eps == 1:
                bump('P5c ambient reversed L_u=1', Lfun(dA, arcZ, 'u', verts) == 1)
            else:
                lam = Lfun(dA, arcZ, 'u', verts) ^ Lfun(dA, arcZ, 'beta', verts) ^ dA.get((frozenset(['u', 'beta']), frozenset(['beta'])), 0)
                bump("P5c ambient unchanged Lambda'=1+p_j", lam == (1 + pj) % 2)
        done += 1
    json.dump(C, open(os.path.join(HERE, 'p5c.json'), 'w'), indent=1, sort_keys=True)
    for key in sorted(C): print('%-60s checks=%3d failures=%d' % (key, C[key][0], C[key][1]))
    print('TOTAL FAILURES', sum(v[1] for v in C.values()))

if __name__ == '__main__':
    main()
