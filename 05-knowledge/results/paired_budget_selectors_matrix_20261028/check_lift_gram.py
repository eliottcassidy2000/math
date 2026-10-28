#!/usr/bin/env python3
"""Extend incoming Session89: full order-nine image, transports, signed Gram."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
from math import lcm, gcd
import sys
import json
import gzip
import hashlib

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/'session89_consecutive_tournament_pairs_20261028/generation'),
               str(HERE.parent/'root_ports_closed_carriers_20261028')]
import check_pairs as pairs
import check_lift as lift
import check_carriers as carrier

CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def transpose(A):
    return list(map(list, zip(*A)))


def mm(A, B):
    return [[sum(a*b for a, b in zip(row, col)) for col in zip(*B)] for row in A]


def det(A):
    A = [list(row) for row in A]; n = len(A)
    if not n:
        return 1
    old = 1; sign = 1
    for k in range(n-1):
        pivot = next((j for j in range(k, n) if A[j][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            A[k], A[pivot] = A[pivot], A[k]; sign = -sign
        a = A[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                value = A[i][j]*a-A[i][k]*A[k][j]
                need(value % old == 0, 'Bareiss exact division')
                A[i][j] = value//old
        old = a
    return sign*A[-1][-1]


def solve(A, b):
    n = len(A); a = [[Fraction(x) for x in row]+[Fraction(v)] for row, v in zip(A, b)]
    for i in range(n):
        k = next(k for k in range(i, n) if a[k][i])
        a[i], a[k] = a[k], a[i]
        z = a[i][i]; a[i] = [x/z for x in a[i]]
        for j in range(n):
            if j != i:
                z = a[j][i]; a[j] = [x-z*y for x, y in zip(a[j], a[i])]
    return [row[-1] for row in a]


def signed_core(n, tcode, scode, mask):
    T = carrier.old.matrix(n, tcode)
    L = [[int(i == j) for j in range(n)] for i in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            a = 2*((T[i] >> j) & 1)-1
            c = 2*(((tcode ^ scode) >> k) & 1)-1
            L[i][j], L[j][i] = a+c, -a+c
            k += 1
    r = [2*((mask >> i) & 1)-1 for i in range(n)]
    B = mm(L, transpose(L))
    G = [[B[i][j]+2*r[i]*r[j] for j in range(n)] for i in range(n)]
    return L, r, G


def algebra(n, key, W, polynomial=False):
    L, r, G = signed_core(n, *key)
    S = [[0 if i == j else 2*((W[i] >> j) & 1)-1 for j in range(2*n+1)] for i in range(2*n+1)]
    need(all(L[i][j] % 2 == int(i == j) for i in range(n) for j in range(n)), 'L = I modulo two')
    need(det(L) % 2 == 1, 'odd core determinant')
    y = solve(L, r); denominator = lcm(*(v.denominator for v in y))
    v = [denominator]+[int(denominator*x) for x in y]*2
    content = gcd(*v); v = [x//content for x in v]
    need(all(x % 2 for x in v), 'primitive integer kernel is odd at every named vertex')
    need(all(sum(a*b for a, b in zip(row, v)) == 0 for row in S), 'named integer kernel survives lift')
    need([sum(L[i][j]*v[1+j] for j in range(n))//v[0] for i in range(n)] == r,
         'kernel and core recover root orientation')
    base = mm(L, transpose(L))
    need([(G[i][0]-base[i][0])//(2*r[0]) for i in range(n)] == r, 'one pinned root sign decodes Gram')
    # Direct block calculation in the unnormalized +/- involution bases.
    for i in range(n):
        p, q = i+1, n+i+1
        need(S[0][p]-S[0][q] == 2*r[i], 'root antisymmetric block')
        for j in range(n):
            a, b = j+1, n+j+1
            need(S[p][a]+S[p][b]+S[q][a]+S[q][b] == 0, 'symmetric block vanishes')
            need(S[p][a]-S[p][b]-S[q][a]+S[q][b] == 0, 'antisymmetric block vanishes')
            need(S[p][a]-S[p][b]+S[q][a]-S[q][b] == -2*L[j][i], 'mixed block is minus L transpose')
    if polynomial:
        # Both polynomials have degree 2n+1; 2n+2 exact evaluation points.
        for x in range(1, 2*n+3):
            left = det([[x*int(i == j)-S[i][j] for j in range(2*n+1)] for i in range(2*n+1)])
            right = x*det([[x*x*int(i == j)+G[i][j] for j in range(n)] for i in range(n)])
            need(left == right, 'full characteristic-polynomial identity')
    return dict(core_determinant=det(L), primitive_kernel=v)


def strongly_connected(rows):
    n = len(rows)
    for v in range(n):
        seen = 1 << v
        while True:
            new = seen
            for i in range(n):
                if seen >> i & 1:
                    new |= rows[i]
            if new == seen:
                break
            seen = new
        if seen != (1 << n)-1:
            return False
    return True


def main():
    old, block = pairs.load_baseline(None)
    prep = pairs.prepare(4, old, block)
    bank, _ = pairs.core_first(prep, old)
    engine = carrier.Canon(); images = {}; profile = Counter(); core_counts = Counter()
    for index, (key, (_, group)) in enumerate(sorted(bank.items())):
        W, tau = lift.lift(4, *key)
        code, maps = engine.maps(W)
        need(code not in images, 'distinct aligned pairs remain distinct at order nine')
        images[code] = key; core_counts[key[0]] += 1
        inverse = pairs.inverse(maps[0])
        aut = {tuple(p[inverse[i]] for i in range(9)) for p in maps}
        centralizer = {g for g in aut if all(g[tau[i]] == tau[g[i]] for i in range(9))}
        lifted = {(0,)+tuple(1+p[i] for i in range(4))+tuple(5+p[i] for i in range(4)) for p in group}
        need(centralizer == lifted, 'all named centralizer transports retained')
        anti = sum(all(g[tau[g[tau[i]]]] == i for i in range(9)) for g in aut)
        need(anti*len(group) == len(aut), 'anti-involution conjugacy fibre')
        profile[(len(group), len(aut), anti)] += 1
        algebra(4, key, W, polynomial=index < 20)
    path = HERE.parent/'relation_tower_frobenius_ports_20261028/tournament_bank.json.gz'
    raw = path.read_bytes(); ordinary = json.loads(gzip.decompress(raw))
    need(ordinary['order'] == 9 and len(ordinary['rows']) == 191536, 'independent full ordinary bank')
    target = set(); eligible = 0
    for row in ordinary['rows']:
        row = tuple(row); degrees = sorted(x.bit_count() for x in row)
        if degrees != [8-d for d in reversed(degrees)]:
            continue
        eligible += 1
        code, _ = engine.maps(row)
        if code == engine.maps(lift.converse(row))[0]:
            target.add(code)
    need(set(images) == target, 'image equals full independent self-converse subset at nine')
    # Four-probe K4 kernel tested on actual strong lifts, not only raw bits.
    edges = [(i, j) for i in range(4) for j in range(i+1, 4)]
    def obs(g):
        return tuple(sum(g >> k & 1 for k, e in enumerate(edges) if v in e) % 2 for v in range(3))+(g.bit_count() % 2,)
    hostile = []
    for g in (0, 45):
        W, _ = lift.lift(4, 0, g, 0)
        need(strongly_connected(W), 'probe collision occurs in strong tournaments')
        hostile.append(dict(disagreement=g, probes=obs(g), canonical_code=engine.maps(W)[0],
                            **algebra(4, (0, g, 0), W, polynomial=True)))
    need(hostile[0]['probes'] == hostile[1]['probes'] and hostile[0]['canonical_code'] != hostile[1]['canonical_code'],
         'four probes lose a genuine strong-atom distinction')
    zero_fibre = {g: det(signed_core(4, 0, g, 0)[0]) for g in range(64) if obs(g) == (0, 0, 0, 0)}
    need(zero_fibre == {0: 1, 30: 1, 45: -31, 51: -31}, 'dyadic determinant restores exactly one kernel bit here')
    a, _ = lift.lift(1, 0, 0, 0); b, _ = lift.lift(1, 0, 0, 1)
    need(signed_core(1, 0, 0, 0)[2] == signed_core(1, 0, 0, 1)[2], 'Gram forgets global root sign')
    need(engine.maps(a)[0] != engine.maps(b)[0], 'forgotten root sign can change isomorphism class')
    controls = carrier.strong_controls()
    determinants = {}
    for n in range(1, 5):
        histogram = Counter()
        for tcode in range(1 << (n*(n-1)//2)):
            for scode in range(1 << (n*(n-1)//2)):
                L, _, _ = signed_core(n, tcode, scode, 0)
                d = det(L)
                need(d % 8 == 1, 'signed-tournament cycle-cover congruence')
                histogram[d] += 1
        determinants[n] = dict(sorted(histogram.items()))
    need(determinants[3] == {-7: 8, 1: 48, 9: 8}, 'triangle sign determinant split')
    # Norm-seven associate and norm-three ROOT phase in integer coordinates.
    import check_paired_budget as phase
    need(phase.mul((0, -1), (3, 2)) == (2, -1), 'correct conjugate associate of norm seven')
    need(phase.Q((1, -2)) == 7 and phase.Q((1, 2)) == 3, 'Eisenstein determinant factors')
    return dict(status='PROVED algebra; FINITE-EXACT full order-nine audit',
                aligned_pairs=len(bank), full_ordinary_bank=191536,
                necessary_score_filter_survivors=eligible, self_converse_image=len(target),
                exact_full_image_equality=True, per_K4_core=dict(sorted(core_counts.items())),
                centralizer_profiles=[dict(pair_aut=a, lift_aut=b, anti_involutions=c, classes=z)
                                      for (a, b, c), z in sorted(profile.items())],
                full_bank_sha256=hashlib.sha256(raw).hexdigest(),
                image_sha256=hashlib.sha256(repr(sorted(images.items())).encode()).hexdigest(),
                K4_probe_strong_hostiles=hostile, retained_81_27_controls=controls,
                K4_zero_probe_fibre_determinants=zero_fibre,
                labelled_signed_core_determinants=determinants,
                checks=CHECKS, source='incoming Session89 exact lift and aligned-pair generator; previous IR canonicalizer')


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
