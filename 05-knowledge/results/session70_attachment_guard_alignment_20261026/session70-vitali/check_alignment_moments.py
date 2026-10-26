#!/usr/bin/env python3
"""Independent native-tournament realization of the Boolean guard hierarchy.
No import from the inherited producer or its core. Explicit checks survive -O.
"""
import hashlib
import itertools
import json
import math
from fractions import Fraction

X = (2, 12, 41, 17, 39, 75, 31)
Y = (2, 20, 9, 51, 37, 71, 31)
SRC = (0, 1, 4, 2, 3, 5, 6)  # map X to Y
DST = (0, 2, 4, 1, 5, 3, 6)  # map tau(X) to tau(Y)
PROBE = (0, 1, 2, 3, 5, 6, 4)
COUNTS = {}


def need(ok, label):
    COUNTS[label] = COUNTS.get(label, 0) + 1
    if not ok:
        raise RuntimeError(label)


def edges(rows):
    return frozenset((i, j) for i, r in enumerate(rows)
                     for j in range(len(rows)) if (r >> j) & 1)


def check_tournament(e, n):
    return len(e) == n * (n - 1) // 2 and all(
        ((i, j) in e) != ((j, i) in e)
        for i in range(n) for j in range(i + 1, n))


def native(bitword, k):
    n = 7 * k + 2
    u, v = n - 2, n - 1
    e = set()
    for block in range(k):
        own = edges(Y if bitword & (1 << block) else X)
        e.update((7 * block + a, 7 * block + b) for a, b in own)
        e.update((a, b) for a in range(7 * block, 7 * (block + 1))
                 for b in range(7 * (block + 1), 7 * k))
    e.update((a, u) for a in range(7 * k))
    e.add((u, v))
    e.update((v, a) for a in range(7 * k))
    return frozenset(e)


def complement(e, n):
    hinges = {(i, i + 1) for i in range(n - 1)}
    return frozenset((a, b) if (a, b) in hinges else (b, a) for a, b in e)


def mapping(k, bitword, local):
    return tuple(7 * i + (local[j] if bitword & (1 << i) else j)
                 for i in range(k) for j in range(7)) + (7*k, 7*k+1)


def mapped(e, p):
    return frozenset((p[a], p[b]) for a, b in e)


def remarking(k, subset):
    return tuple(7 * i + j for i in range(k)
                 for j in (PROBE if subset & (1 << i) else range(7))) + (7*k, 7*k+1)


def legal(e, p):
    return all((a, b) in e for a, b in zip(p, p[1:]))


def strong(e, n):
    for start in (0,):
        for reverse in (False, True):
            reached = {start}
            while True:
                more = reached | {b for a, b in ([(y,x) for x,y in e] if reverse else e)
                                  if a in reached}
                if more == reached:
                    break
                reached = more
            if len(reached) != n:
                return False
    return True


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            top = row.bit_length() - 1
            if top in pivots:
                row ^= pivots[top]
            else:
                pivots[top] = row
                break
    return len(pivots)


def transform(weights, k):
    # Independent fast superset transform, not explicit native-guard multiplication.
    out = list(weights)
    for bit in range(k):
        for mask in range(1 << k):
            if not mask & (1 << bit):
                out[mask] += out[mask | (1 << bit)]
    return out


def invert(moments, k):
    out = list(moments)
    for bit in range(k):
        for mask in range(1 << k):
            if not mask & (1 << bit):
                out[mask] -= out[mask | (1 << bit)]
    return out


def run():
    need(check_tournament(edges(X), 7) and check_tournament(edges(Y), 7), 'seed tournaments')
    need(mapped(edges(X), SRC) == edges(Y), 'independent seed source isomorphism')
    need(mapped(complement(edges(X), 7), DST) == complement(edges(Y), 7),
         'independent seed target isomorphism')
    need(not legal(edges(X), PROBE) and legal(edges(Y), PROBE), 'seed guard distinguishes')
    seed_rigidity = {}
    for name, relation in [('X',edges(X)), ('tau_X',complement(edges(X),7))]:
        automorphisms = []
        for p in itertools.permutations(range(7)):
            if mapped(relation,p)==relation:
                automorphisms.append(p)
        need(automorphisms==[tuple(range(7))], 'exhaustive seed rigidity')
        seed_rigidity[name] = {'permutations_tested':math.factorial(7),'automorphisms':automorphisms}
    relative = tuple(DST.index(SRC[i]) for i in range(7))
    need(relative==(0,3,2,1,5,4,6), 'relative frame mismatch double transposition')
    need(mapped(edges(X),relative)!=edges(X), 'relative frame mismatch is not source automorphism')
    table = []
    for k in range(1, 8):
        n, m = 7*k+2, 1 << k
        base, base_target = native(0, k), complement(native(0,k), n)
        guards = [remarking(k, subset) for subset in range(m)]
        matrix = [[0] * m for _ in range(m)]
        states = set()
        for b in range(m):
            e = native(b, k)
            states.add(e)
            need(check_tournament(e,n), 'native tournament')
            need(legal(e, tuple(range(n))), 'standard marked path')
            need(strong(e,n), 'strong native source')
            degrees = [sum(a==i for a,z in e) for i in range(n)]
            need(degrees[-2:]==[1,7*k] and all(2<=d<=7*k-1 for d in degrees[:-2]),
                 'rigidity wrapper degree separation')
            need(mapped(base, mapping(k,b,SRC)) == e, 'paired rooted source fixed')
            need(mapped(base_target, mapping(k,b,DST)) == complement(e,n),
                 'paired rooted target fixed')
            need((mapping(k,b,SRC)[0],mapping(k,b,SRC)[-1]) == (0,n-1), 'source roots fixed')
            need((mapping(k,b,DST)[0],mapping(k,b,DST)[-1]) == (0,n-1), 'target roots fixed')
            # Literal colour: sigma=converse followed by positional reversal.
            sigma = frozenset((n-1-z,n-1-a) for a,z in e)
            need(sigma != e, 'all black')
            for subset,p in enumerate(guards):
                native_answer = int(legal(e,p))
                matrix[subset][b] = native_answer
                need(native_answer == int((subset & b) == subset), 'native guard equals Boolean monomial')
        need(len(states)==m, 'distinct native alignments')
        # Numerical subset order makes zeta upper triangular, determinant exactly 1.
        need(all(matrix[s][b] == 0 for s in range(m) for b in range(s)), 'integer zeta upper triangular')
        need(all(matrix[b][b] == 1 for b in range(m)), 'integer zeta diagonal one')
        packed = [sum(value << b for b,value in enumerate(row)) for row in matrix]
        need(binary_rank(packed)==m, 'independent full GF2 rank')
        ranks = []
        for t in range(k+1):
            expected = sum(math.comb(k,j) for j in range(t+1))
            rank = binary_rank([packed[s] for s in range(m) if s.bit_count()<=t])
            need(rank==expected, 'all truncation ranks')
            ranks.append({'t':t,'rank':rank,'missing_linear_rank':m-rank})
        parity = [[int(b.bit_count()%2==eps) for b in range(m)] for eps in (0,1)]
        moments = [[sum(matrix[s][b]*law[b] for b in range(m)) for s in range(m)] for law in parity]
        need(all(sum(law)==m//2 for law in parity), 'parity laws equal total mass')
        for s in range(m-1):
            need(moments[0][s]==moments[1][s]==1 << (k-s.bit_count()-1), 'all proper guard moments equal')
        need(moments[k%2][m-1]==1 and moments[1-k%2][m-1]==0, 'native top guard separates parity laws')
        # Every proper marginal, not merely all-one moments.
        for s in range(m-1):
            marginals = [{v:sum(law[b] for b in range(m) if b&s==v)
                          for v in range(m) if v&s==v} for law in parity]
            need(marginals[0]==marginals[1], 'all proper full marginals equal')
        weights = [b*b+3*b+1 for b in range(m)]
        direct = [sum(matrix[s][b]*weights[b] for b in range(m)) for s in range(m)]
        need(transform(weights,k)==direct, 'native matrix agrees with independent zeta transform')
        need(invert(direct,k)==weights, 'complete named moments invert')
        # Deliberately erase one native guard; hostile parity law must now collide.
        need(moments[0][:-1]==moments[1][:-1], 'top-deleted hostile collides')
        need(moments[0]!=moments[1], 'complete carrier repairs hostile')
        table.append({'k':k,'n':n,'alignments':m,'native_guard_tests':m*m,
                      'full_rank':m,'truncation_ranks':ranks,
                      'top_probability_gap':str(Fraction(1,m//2)),
                      'matrix_sha256':hashlib.sha256(bytes(v for row in matrix for v in row)).hexdigest()})
    # Boundary: k=0 has one state and only the constant guard; no parity pair.
    need(transform([7],0)==[7] and invert([7],0)==[7], 'zero-dimensional boundary')
    # Classical finite section obstruction, explicitly distinct from matrix rank.
    # C3 acts freely on its three Hamilton paths, so no fixed path; a reference does select.
    triangle = frozenset(((0,1),(1,2),(2,0)))
    hp = [p for p in itertools.permutations(range(3)) if legal(triangle,p)]
    need(len(hp)==3, 'triangle has three Hamilton paths')
    need(all(tuple((v+1)%3 for v in p)!=p for p in hp), 'no rotation-fixed physical path')
    need(legal(triangle,(0,1,2)), 'retained reference path is a lawful chart choice')
    return {'status':'PASS','scope':'all 2^k C4 states and all 2^k simultaneous block guards, k=1..7; theorem all k',
            'tables':table,'checks':COUNTS,'total_checks':sum(COUNTS.values()),
            'seed_rigidity':seed_rigidity,'relative_frame_mismatch':relative,
            'imports_inherited_producer':False}

if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
