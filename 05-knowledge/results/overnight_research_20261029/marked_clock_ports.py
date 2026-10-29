"""Exact decoder for labelled paired cokernels of matched C3 cross joins."""
from fractions import Fraction as F
from itertools import product
from math import lcm
from pathlib import Path
import json
from port_queries import invdet, mm, tr
from eisenstein_clock_joins import source, norm, A

GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def residue(x):
    return x % 1


def port_order(inverse, i):
    return lcm(*(row[i].denominator for row in inverse))


def main():
    vandermonde = [[F(t)**j for j in range(5)] for t in range(5)]
    vi, _ = invdet(vandermonde)
    T = [[A[i%2][j%2] if i//2 == j//2 else 0 for j in range(4)] for i in range(4)]
    ti, _ = invdet(T)
    dual_clock = tr(ti)
    identity = [[int(i == j) for j in range(4)] for i in range(4)]
    need(mm(mm(dual_clock, dual_clock), dual_clock) == identity, 'dual clock order3')
    tested = singular_probes = wrong_clock_witnesses = 0
    for d, e in product([x for x in range(-6, 7) if x], repeat=2):
        for a, b in product(range(-2, 3), repeat=2):
            K = source(d, e, a, b)
            ki, det = invdet(K)
            if not det:
                continue
            need(mm(tr(T), mm(K, T)) == K, 'primal integral clock')
            need(mm(tr(dual_clock), mm(ki, dual_clock)) == ki, 'dual cokernel pairing clock')
            lattice_action = mm(ki, mm(dual_clock, K))
            need(lattice_action == T, 'dual action preserves the relation lattice')
            wrong = mm(tr(T), mm(ki, T))
            wrong_clock_witnesses += any(residue(wrong[i][j]-ki[i][j])
                                         for i, j in product(range(4), repeat=2))
            coeff = []
            for axis in range(2):
                observations = []
                for w in range(5):
                    aa, bb = (a+w, b) if axis == 0 else (a, b+w)
                    _, detw = invdet(source(d, e, aa, bb))
                    need(detw == (d*e-norm(aa, bb))**2, 'group order or singular flag')
                    observations.append([detw])
                    singular_probes += not bool(detw)
                coeff.append([row[0] for row in mm(vi, observations)])
            l1, l2 = coeff[0][3]/2, coeff[1][3]/2
            aa, bb = (2*l1-l2)/3, (2*l2-l1)/3
            z = (l1*l1-coeff[0][2])/2 + norm(aa, bb)
            need((aa, bb, z) == (a, b, d*e), 'group orders recover response carrier')
            zero, _ = invdet(source(d, e, 0, 0))
            D, E = port_order(zero, 0), port_order(zero, 2)
            need((D, E) == (abs(d), abs(e)), 'labelled port orders after cancelling phase')
            L = 2*D+2*E+abs(z)+1
            distant, _ = invdet(source(d, e, L, 0))
            t = z-L*L
            need(abs(t) > 2*max(D, E), 'pairing residue has a unique signed lift')
            need(distant[0][1] == -F(e)/t and distant[2][3] == -F(d)/t,
                 'direct inverse verifies the two named pairing residues')
            found_d = [x for x in (-D, D) if residue(-F(x)/t) == residue(distant[2][3])]
            found_e = [x for x in (-E, E) if residue(-F(x)/t) == residue(distant[0][1])]
            need(found_d == [d] and found_e == [e], 'fully marked signed source decoder')
            need(source(-d, -e, -a, -b) == [[-x for x in row] for row in K],
                 'closed global sign operation')
            tested += 1
    need(wrong_clock_witnesses > 0, 'primal versus dual clock hostile is substantive')
    # Initial marked pairings cannot see global sign on a trivial group.
    left, _ = invdet(source(1, 1, 0, 0))
    right, _ = invdet(source(-1, -1, 0, 0))
    need([[residue(x) for x in row] for row in left] ==
         [[residue(x) for x in row] for row in right], 'initial trivial-group sign hostile')
    left, _ = invdet(source(1, 1, 5, 0))
    right, _ = invdet(source(-1, -1, 5, 0))
    need(residue(left[0][1]) != residue(right[0][1]), 'one shared future attachment separates signs')
    out = {'gates': GATES, 'nonsingular_sources_decoded': tested,
           'source_range': 'd,e=-6..6 nonzero; a,b=-2..2',
           'singular_branch_probe_values': singular_probes,
           'wrong_primal_clock_pairing_failures': wrong_clock_witnesses,
           'minimal_carrier': ['d','e','a','b'],
           'scope': 'four named generators, exact pairing, known induced dual clock; branching experiments'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
