#!/usr/bin/env python3
"""Independent exact checks of Session31 categorical and metagraph addenda.

No source-worker code is imported. No Git or external state is changed.
Class enumeration starts from fixed-path presentations; forests use parent
functions, Hamiltonian paths use permutations, and automorphisms use relabeling.
The asymptotic check uses small matching counts, not trigonometric sums.
"""
import hashlib
import itertools as it
import json
import math
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def categorical_checks():
    x, y, a = s.symbols('x y a')
    # Eliminate sqrt(2) from beta^4 = 3/sqrt(2) - 2.
    resultant = s.resultant(y*y-2, y*(x**4+2)-3, y)
    p = 2*x**8+8*x**4-1
    assert s.expand(resultant+p) == 0
    reciprocal = s.Poly(-s.expand(x**8*p.subs(x, 1/x)), x)
    assert reciprocal.as_expr() == x**8-8*x**4-2
    coeff = reciprocal.all_coeffs()
    assert coeff[0] == 1
    assert all(int(c) % 2 == 0 for c in coeff[1:])
    assert int(coeff[-1]) % 4 != 0
    # Eisenstein plus reciprocal substitution proves degree 8 without minpoly.
    norm_beta = s.Poly(p, x).TC()/s.Poly(p, x).LC()
    assert norm_beta == -s.Rational(1, 2)

    primes = list(s.primerange(2, 32))
    for prime in primes:
        pp = s.Poly(x**prime-prime, x)
        assert all(int(c) % prime == 0 for c in pp.all_coeffs()[1:])
        assert int(pp.TC()) % (prime*prime) != 0
        assert s.Poly(((x+1)**prime-prime).expand(), x).degree() == prime

    # A finite-order operator without parity compatibility is a hostile.
    G = s.Matrix([[a, 1], [1-a*a, -a]])
    parity = s.diag(1, -1)
    assert s.simplify(G*G) == s.eye(2)
    assert s.trace(parity*G) == 2*a
    assert G*parity != parity*G

    # The exact fourth-Chebyshev relation, as a polynomial identity.
    mu = s.symbols('mu')
    assert s.expand(-2*s.chebyshevt(4, mu/2)) == -mu**4+4*mu**2-2
    central = lambda g: s.cancel(1-6*(g-1)**2/g)
    branch_rows = []
    for m in [3, 4, 8, 12]:
        gd, gn = 2-s.Rational(4, m), s.Rational(4, m)
        branch_rows.append({'m': m, 'proposed_dilute_g': str(gd),
                            'proposed_dense_g': str(gn),
                            'first_c': str(central(gd)),
                            'second_c': str(central(gn)),
                            'labels_have_standard_g_ranges': bool(1 <= gd <= 2 and 0 < gn <= 1)})
    assert branch_rows[0]['labels_have_standard_g_ranges'] is False
    assert branch_rows[2]['first_c'] == '0'
    assert branch_rows[2]['second_c'] == '-2'
    return {'status': 'PASS with explicit scope repairs',
            'beta_elimination_resultant': str(resultant),
            'inverse_eisenstein_polynomial': str(reciprocal.as_expr()),
            'beta_degree': 8, 'beta_norm': str(norm_beta),
            'prime_degree_controls': primes,
            'parity_hostile': {'matrix': str(G), 'square': 'I', 'supertrace': '2*a'},
            'branch_controls': branch_rows}


def asymptotic_checks():
    z, t = s.symbols('z t', nonzero=True)
    n = 1/z
    v = z*z/(1+2*z/t)**2
    partial = 1
    # binom(n-1-l,l) = product(n-1-l-j,j=0..l-1)/l!.
    for ell in range(1, 5):
        partial += s.prod(n-1-ell-j for j in range(ell))/math.factorial(ell)*v**ell
    series = s.series(partial, z, 0, 5).removeO().expand()
    expected = (1+z-(s.Rational(3, 2)+4/t)*z*z
                +(4/t+12/t**2-s.Rational(10, 3))*z**3)
    assert s.expand(series-expected).coeff(z, 0) == 0
    assert all(s.expand(series-expected).coeff(z, j) == 0 for j in range(1, 4))
    fourth = s.factor(series.coeff(z, 4))
    assert s.expand(fourth-(s.Rational(85, 24)+26/t-4/t**2-32/t**3)) == 0

    # Rigorous rational finite-n remainder certificate, using v>=0 and nv<1.
    # The omitted ell>=4 sum is <= sum_{ell>=4}(nv)^ell/ell!
    # <= (nv)^4/(24*(1-nv/5)); ratio of successive terms is <=nv/5.
    rows = []
    for nn in (8, 16, 32, 64):
        for tt in (Q(1, 2), Q(1), Q(2), Q(7, 3)):
            aa = 1+Q(nn, 2)*tt
            vv = tt*tt/(4*aa*aa)
            terms = [math.comb(nn-1-ell, ell)*vv**ell
                     for ell in range((nn-1)//2+1)]
            actual = sum(terms, Q())
            tail = sum(terms[4:], Q())
            nv = nn*vv
            bound = nv**4/(24*(1-nv/5))
            assert 0 <= tail <= bound
            rows.append({'n': nn, 't': str(tt), 'ratio_exact': str(actual),
                         'ell_ge_4_tail': str(tail), 'certified_upper_bound': str(bound)})
    return {'status': 'PASS', 'independent_method': 'finite matching-binomial expansion and rational tail bound',
            'through_order_3': str(expected), 'order_4_coefficient': str(fourth),
            'finite_exact_remainder_checks': rows}


def forest_coefficients(out):
    """Count spanning rooted directed forests with arcs toward the root."""
    n = len(out)
    result = [0]*n
    for parents in it.product(*[(-1,)+tuple(sorted(row)) for row in out]):
        good = True
        for v in range(n):
            seen = set()
            while v != -1:
                if v in seen:
                    good = False
                    break
                seen.add(v)
                v = parents[v]
            if not good:
                break
        if good:
            result[sum(v != -1 for v in parents)] += 1
    return result


def tournament_checks(n):
    pairs = list(it.combinations(range(n), 2))
    index = {pair: j for j, pair in enumerate(pairs)}
    allbits = (1 << len(pairs))-1
    perms = list(it.permutations(range(n)))
    # Relabeling is represented by listing old vertices in the new order.
    transforms = []
    for perm in perms:
        transforms.append([(index[tuple(sorted((perm[i], perm[j])))], int(perm[i] > perm[j]))
                           for i, j in pairs])
    canonical = {}
    classes = {}

    def get_class(mask):
        if mask in canonical:
            return canonical[mask]
        orbit = set()
        for transform in transforms:
            image = 0
            for j, (oldbit, flip) in enumerate(transform):
                image |= (((mask >> oldbit) & 1) ^ flip) << j
            orbit.add(image)
        representative = min(orbit)
        canonical.update({m: representative for m in orbit})
        aut = math.factorial(n)//len(orbit)
        assert aut*len(orbit) == math.factorial(n)
        classes[representative] = {'automorphisms': aut, 'labelled_orbit_size': len(orbit)}
        return representative

    pathmask = sum(1 << index[(i, i+1)] for i in range(n-1))
    freebits = [index[pair] for pair in pairs if pair[1] != pair[0]+1]
    presentation_counts = Counter()
    for filling in range(1 << len(freebits)):
        mask = pathmask
        for j, bit in enumerate(freebits):
            mask |= ((filling >> j) & 1) << bit
        presentation_counts[get_class(mask)] += 1
    assert len(canonical) == 1 << len(pairs)

    def outgoing(mask):
        out = [set() for _ in range(n)]
        for j, (a, b) in enumerate(pairs):
            if mask >> j & 1:
                out[a].add(b)
            else:
                out[b].add(a)
        return out

    for representative, data in classes.items():
        out = outgoing(representative)
        hpaths = sum(all(b in out[a] for a, b in zip(perm, perm[1:])) for perm in perms)
        assert hpaths % data['automorphisms'] == 0
        tilings = hpaths//data['automorphisms']
        assert tilings == presentation_counts[representative]
        data.update({'representative_mask_bit1_lower_to_higher': representative,
                     'H': hpaths, 'h': tilings,
                     'forest_coefficients': forest_coefficients(out),
                     'signed_score_energy': sum((2*len(row)-(n-1))**2 for row in out),
                     'converse_representative': canonical[representative ^ allbits]})

    merged = []
    seen = set()
    mass = [Q() for _ in range(n)]
    energy_mass = 0
    total_weight = 0
    for representative in sorted(classes):
        if representative in seen:
            continue
        data = classes[representative]
        converse = data['converse_representative']
        other = classes[converse]
        seen.update((representative, converse))
        assert data['H'] == other['H']
        assert data['automorphisms'] == other['automorphisms']
        assert data['signed_score_energy'] == other['signed_score_energy']
        weight = data['h']*(1 if converse == representative else 2)
        fbar = [Q(a+b, 2) for a, b in zip(data['forest_coefficients'], other['forest_coefficients'])]
        total_weight += weight
        energy_mass += weight*data['signed_score_energy']
        mass = [v+weight*f for v, f in zip(mass, fbar)]
        merged.append({'representatives': sorted({representative, converse}), 'weight': weight,
                       'forest_average': list(map(str, fbar)),
                       'signed_score_energy': data['signed_score_energy']})

    d = (n-1)*(n-2)//2
    expected = [Q() for _ in range(n)]
    for ell in range((n-1)//2+1):
        for j in range(n-2*ell):
            expected[j+2*ell] += (1 << d)*Q(math.comb(n-1-ell, ell), 4**ell)*math.comb(n-1-2*ell, j)*Q(n, 2)**j
    assert mass == expected
    assert total_weight == 1 << d
    assert energy_mass == (1 << d)*(n*n-3*n+4)

    collision = None
    if n == 5:
        # Root's mask convention is the opposite: bit 1 means higher -> lower.
        records = []
        for rootmask in (40, 41):
            mask = allbits ^ rootmask
            out = outgoing(mask)
            H = sum(all(b in out[a] for a, b in zip(p, p[1:])) for p in perms)
            aut = classes[canonical[mask]]['automorphisms']
            fout = forest_coefficients(out)
            fin = forest_coefficients(outgoing(mask ^ allbits))
            records.append({'root_mask': rootmask, 'H': H, 'automorphisms': aut,
                            'out_forest': fout, 'in_forest': fin,
                            'arcs': [[a, b] for a in range(n) for b in sorted(out[a])]})
        assert records[0]['H'] == 15 and records[1]['H'] == 13
        assert records[0]['automorphisms'] == 3 and records[1]['automorphisms'] == 1
        assert all(r['out_forest'] == r['in_forest'] == [1, 10, 39, 70, 49] for r in records)
        collision = records

    return {'n': n, 'labelled_tournaments_covered': len(canonical),
            'unlabelled_classes': len(classes), 'merged_nodes': len(merged),
            'fixed_path_presentations': 1 << d, 'total_weight': total_weight,
            'weighted_forest_coefficients': list(map(str, mass)),
            'weighted_signed_score_energy': energy_mass,
            'classes': list(classes.values()), 'merged_inventory': merged,
            'full_forest_collision': collision}


def main():
    targets = ['session31-categories/CATEGORICAL-HONEYCOMB-MRB-AUDIT.md',
               'session31-categories/check_categorical_parameters.py',
               'session31-root/HAMILTONIAN_FOREST_SIZE_BIAS.md']
    def resolve_source(relative):
        first, rest = relative.split('/', 1)
        folder = first.removeprefix('session31-')
        portable = ROOT / (rest if folder == 'root' else folder + '/' + rest)
        return portable if portable.is_file() else ROOT / relative
    hashes = {p: hashlib.sha256(resolve_source(p).read_bytes()).hexdigest() for p in targets}
    report = {'status': 'PASS of addenda and corrected-scope categorical statements',
              'source_sha256': hashes, 'categorical': categorical_checks(),
              'asymptotic': asymptotic_checks(),
              'metagraph': [tournament_checks(n) for n in range(2, 7)]}
    dest = HERE/'audit_bridges.json'
    dest.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': report['status'], 'report': str(dest),
                      'counts': [{k: row[k] for k in ['n', 'unlabelled_classes', 'merged_nodes', 'fixed_path_presentations']}
                                 for row in report['metagraph']]}, indent=2))


if __name__ == '__main__':
    main()
