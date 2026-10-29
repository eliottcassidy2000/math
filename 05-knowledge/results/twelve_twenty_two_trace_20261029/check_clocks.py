"""Exact power-clock, integral-observer and Frobenius/Paley certificates.

No producer imports. All rings, matrices, ports and exceptional primes are
specified here. This is a same-agent exact check, not an independent audit.
"""
from collections import Counter
from functools import reduce
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import json

from flint import fmpz_mat
from sympy import Matrix, eye, factorint, primerange

HERE = Path(__file__).resolve().parent
ONE, U = (1, 0, 0), (0, 1, 0)


def need(ok, why):
    if not ok:
        raise ValueError(why)


def add(a, b, modulus=None):
    v = tuple(x+y for x, y in zip(a, b))
    return tuple(x % modulus for x in v) if modulus else v


def scale(n, a, modulus=None):
    v = tuple(n*x for x in a)
    return tuple(x % modulus for x in v) if modulus else v


def mul(a, b, modulus=None):
    c = [0]*5
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    for i in (4, 3):
        # U^3 = -1 + U + 2 U^2.
        c[i-3] -= c[i]
        c[i-2] += c[i]
        c[i-1] += 2*c[i]
    v = tuple(c[:3])
    return tuple(x % modulus for x in v) if modulus else v


def power(a, n, modulus=None):
    out = ONE
    while n:
        if n & 1:
            out = mul(out, a, modulus)
        a = mul(a, a, modulus)
        n //= 2
    return out


def tr(a):
    return 3*a[0]+2*a[1]+6*a[2]


def sigma(a, modulus=None):
    image = (0, -2, 1)
    return add(add(scale(a[0], ONE), scale(a[1], image)),
               scale(a[2], mul(image, image)), modulus)


def observer(times):
    return [[tr(power(U, t+i)) for i in range(3)] for t in times]


def smith(rows):
    S = fmpz_mat([[int(x) for x in row] for row in rows]).snf()
    return [abs(int(S[i, i])) for i in range(min(S.nrows(), S.ncols()))]


def minor_gcd(rows, size):
    result = 0
    for rs in combinations(range(len(rows)), size):
        for cs in combinations(range(len(rows[0])), size):
            result = gcd(result, int(Matrix([[rows[i][j] for j in cs] for i in rs]).det()))
    return abs(result)


def lattice_index(rows):
    if len(rows) < 3:
        return 0
    d = smith(rows)
    return d[0]*d[1]*d[2]


def rank_mod(rows, p):
    A = [[x % p for x in row] for row in rows]
    rank = 0
    for j in range(len(A[0])):
        pivot = next((i for i in range(rank, len(A)) if A[i][j]), None)
        if pivot is None:
            continue
        A[rank], A[pivot] = A[pivot], A[rank]
        A[rank] = [x*pow(A[rank][j], -1, p) % p for x in A[rank]]
        for i in range(rank+1, len(A)):
            c = A[i][j]
            A[i] = [(x-c*y) % p for x, y in zip(A[i], A[rank])]
        rank += 1
        if rank == len(A):
            break
    return rank


def run():
    C = Matrix([[1, 0, 1], [1, 1, 1], [1, 1, 0]])
    ell, start = Matrix([[1,1,1]]), Matrix([1,0,0])
    count_observer = Matrix.vstack(ell,ell*C,ell*C*C)
    need(count_observer.det() == -1, 'actual safe-word readout is integrally observable')
    # Literal six-state grammar, with the inherited legality/singularity distinction.
    # States A,B0,B1,FA,FB,X; columns/letters a,b,c,d,e.
    transitions = ((0,2,5,5,5),(5,5,0,1,2),(5,5,0,1,4),
                   (3,4,5,5,5),(5,5,3,4,4),(5,5,5,5,5))
    grammar_counts = []
    for length in range(8):
        count = 0
        for word in product(range(5),repeat=length):
            state = 1
            for letter in word:
                state = transitions[state][letter]
            count += state in (0,1,2)
        need(count == int((ell*C**length*start)[0]), 'independent literal grammar enumeration')
        grammar_counts.append(count)
    D = C-eye(3)
    W = Matrix([[1, 0, -1], [1, -1, -1], [0, -1, -1]])
    S = Matrix([[0, -1, 2], [-1, -2, 4], [-1, -3, 5]])
    need(S.det() == 1 and W*S == S*D, 'inherited integral exterior-dual control')
    need(C == C**4-5*C**2+2*eye(3), 'squared counter retains the integer generator')
    need(D == 3*D**2-D**4-eye(3), 'squared dual retains its integer generator')
    need(power(U, 3) == (-1, 1, 2), 'declared monogenic ring')
    basis = (ONE, U, (0, 0, 1))
    for a in basis:
        need(sigma(sigma(sigma(a))) == a, 'order-three Galois action')
        for b in basis:
            need(sigma(mul(a, b)) == mul(sigma(a), sigma(b)), 'integral algebra action')
    gal = Matrix.hstack(*(Matrix(sigma(v)) for v in basis))
    need(gal.det() == 1 and tr(sigma(U)) == tr(U), 'unimodular trace-preserving phase')
    need(mul(U, mul(sigma(U), sigma(sigma(U)))) == (-1, 0, 0), 'unit norm minus one')
    need(mul(sigma(add(U, scale(-1, ONE))), U) == (-1, 0, 0), 'Galois dual is negative inverse')

    powers = []
    for k in (1, 2, 12, 22):
        L = [power(U, k*j) for j in range(3)]
        O = observer([0, k, 2*k])
        index = abs(int(Matrix(L).det()))
        ds = smith(O)
        native = Matrix.vstack(*(ell*C**(k*j) for j in range(3)))
        native_ds = smith(native.tolist())
        need([minor_gcd(O, i) for i in (1, 2, 3)]
             == [ds[0], ds[0]*ds[1], ds[0]*ds[1]*ds[2]], 'independent determinantal Smith check')
        need(abs(int(Matrix(O).det())) == 49*index, 'trace discriminant times power index')
        powers.append(dict(exponent=k, index=index, index_factors={str(p):e for p,e in factorint(index).items()},
                           smith=ds,counting_smith=native_ds))
    need([r['index'] for r in powers] == [1, 1, 2739659, 3282660883657], 'power-index controls')
    need(gcd(powers[2]['index'], powers[3]['index']) == 1, 'the two pure-clock lattices jointly saturate')

    times = [0, 12, 24, 22, 44]
    rows = [power(U, t) for t in times]
    bank_records = []
    for size in range(3, 6):
        for bank in combinations(times, size):
            bank_records.append(dict(times=bank, index=lattice_index([power(U, t) for t in bank])))
    need(min(len(r['times']) for r in bank_records if r['index'] == 1) == 4,
         'four readings necessary in the specified pure-clock bank')
    chosen = [0, 12, 22, 44]
    H, transform = fmpz_mat([power(U, t) for t in chosen]).hnf(transform=True)
    need(H == transform*fmpz_mat([power(U, t) for t in chosen]), 'integer reconstruction witness')
    need(H.table()[:3] == [[1,0,0],[0,1,0],[0,0,1]], 'chosen bank recovers the full generator lattice')
    need(smith(observer(chosen)) == [1, 7, 7], 'full primitive trace lattice restored')
    need(smith(Matrix.vstack(*(ell*C**t for t in chosen)).tolist()) == [1,1,1],
         'chosen bank reconstructs every actual counter state over every modulus')

    reachable = sorted({12*a+22*b for a in range(5) for b in range(3) if 12*a+22*b <= 48})
    horizons = [dict(time=t, index=lattice_index([power(U, s) for s in reachable if s <= t]))
                for t in reachable]
    need(next(r['time'] for r in horizons if r['index'] == 1) == 44, 'earliest complete integral horizon')
    triples = [dict(times=bank, determinant=int(Matrix([power(U, t) for t in bank]).det()))
               for bank in combinations(reachable, 3)]
    unit_triples = [r for r in triples if abs(r['determinant']) == 1]
    need(unit_triples == [dict(times=(44,46,48), determinant=-1)], 'first and unique triple through48')
    need(power(U, 12, 13) == ONE and rank_mod(observer([0,12,22,24,34,36]),13) == 2
         and rank_mod(observer(chosen),13) == 3, 'mechanism of the horizon44 obstruction')

    finite = []
    for modulus in range(1, 25):
        O, chosen_O = observer([0,1,2]), observer(chosen)
        mapping, primitive_images = {}, set()
        for source in product(range(modulus), repeat=3):
            y = tuple(sum(a*b for a,b in zip(row,source)) % modulus for row in O)
            z = tuple(sum(a*b for a,b in zip(row,source)) % modulus for row in chosen_O)
            need(z not in mapping or mapping[z] == y, 'common-source reconstruction at every tested modulus')
            mapping[z] = y
            primitive_images.add(y)
        expected = modulus**3 // gcd(modulus,7)**2
        need(len(mapping) == len(primitive_images) == expected, 'exact trace-carrier cardinality')
        finite.append(dict(modulus=modulus, states=expected))
    nil = (-3,1,0)
    need(power(nil,3,7) == (0,0,0) and power(nil,2,7) != (0,0,0), 'ramified cubic nilpotent')
    need(all(tr(mul(nil,power(U,j))) % 7 == 0 for j in range(12)), 'trace blind direction survives every primitive clock')
    N = C-3*eye(3)
    need(N**3 % 7 == Matrix.zeros(3) and C**7 % 7 == 3*eye(3),
         'incoming ramified-clock control')
    need([int((ell*N**j*start)[0]) for j in range(3)] == [1,0,-2],
         'native boundary sees the nilpotent second moment')
    for n in range(84):
        need(int((ell*C**n*start)[0]) % 7 == (pow(3,n,7)*(1-n*(n-1)//2)) % 7,
             'incoming all-degree binomial law, independent native boundary check')
    need(power(U,13,3) == (2,0,0) and power(U,26,3) == ONE, 'legal even steps omit odd clock phase')

    # The quadratic u=1+sqrt(2), represented in the basis (1,u).
    A = Matrix([[0,1],[1,2]])
    pell_primitive = [[int((A**(i+j)).trace()) for i in range(2)] for j in range(2)]
    pell_even = [[int((A**(i+2*j)).trace()) for i in range(2)] for j in range(2)]
    need(smith(pell_primitive) == [2,4] and smith(pell_even) == [2,8], 'Pell residual dyadic coordinate')
    need((A**12+eye(2)) % 11 == Matrix.zeros(2), 'Ramanujan22 twelfth clock collapses modulo11')
    z = [-1,1]
    need(all(sum(a*b for a,b in zip(row,z)) % 8 == 0 for row in pell_even)
         and sum(a*b for a,b in zip(pell_primitive[1],z)) % 8 == 4, 'explicit mod8 Pell phase witness')
    need(int((A**12).trace())*64+24 == 2508952, 'inherited Ramanujan22 integer trace')
    need((C**12-eye(3)) % 13 == Matrix.zeros(3), 'split prime13 power collapse')

    prime_controls = []
    dual = (-1,1,0)
    for p in primerange(2, 101):
        if p == 7:
            continue
        residue = p % 7
        phase = 0 if residue in (1,6) else 1 if residue in (2,5) else 2
        image = U
        for _ in range(phase):
            image = sigma(image,p)
        need(power(U,p,p) == tuple(a%p for a in image), 'Frobenius class from seventh-root exponent')
        target = power(U,2,p) if phase == 0 else dual if phase == 1 else sigma(sigma(dual,p),p)
        need(power(U,p+1,p) == tuple(a%p for a in target), 'norm/exterior dual clock law')
        need(power(U,2*p,p) == mul(image,image,p), 'doubled Frobenius clock')
        prime_controls.append(dict(prime=p, phase=phase))
    need(power(U,266,11) == ONE and all(power(U,266//q,11) != ONE for q in (2,7,19)), 'order266 in F1331')
    need(power(U,12,11) == sigma(sigma(dual,11),11), 'literal12-to-dual map at11')
    need(power(U,22,11) == sigma(sigma(power(U,2,11),11),11), 'literal22-to-square map at11')
    need(not any((x**3-2*x*x-x+1) % 11 == 0 for x in range(11)), 'declared cubic quotient is F1331')

    field8 = [power(U,j,2) for j in range(7)]
    need(len(set(field8)) == 7 and power(U,7,2) == ONE, 'counter generates F8 multiplicative points')
    need(power(U,3,2) == tuple(x%2 for x in dual), 'exterior shift changes the generator by a cube')
    squares = {1,2,4}
    arc = lambda i,j: (j-i) % 7 in squares
    need(all(arc(i,j) == arc(2*i%7,2*j%7) for i in range(7) for j in range(7) if i!=j), 'Frobenius preserves Paley orientation')
    need(all(arc(i,j) != arc(3*i%7,3*j%7) for i in range(7) for j in range(7) if i!=j), 'dual generator reverses Paley orientation')
    automorphisms = [perm for perm in permutations(range(7))
                     if all(arc(i,j)==arc(perm[i],perm[j]) for i in range(7) for j in range(i+1,7))]
    affine = {tuple((a*i+b)%7 for i in range(7)) for a in squares for b in range(7)}
    need(set(automorphisms) == affine, 'inherited Paley group literally generated by counter and Frobenius')
    hamilton = sum(all(arc(a,b) for a,b in zip(perm,perm[1:])) for perm in permutations(range(7)))
    need(hamilton == 189, 'inherited Paley Hamilton control')
    return dict(source_matrix=[list(map(int,row)) for row in C.tolist()],
                actual_counting_observer=[[int(a) for a in row] for row in count_observer.tolist()],grammar_counts=grammar_counts,
                galois_matrix=[list(map(int,row)) for row in gal.tolist()],
                power_observers=powers, banks=bank_records, reconstruction_times=chosen,
                reconstruction_transform=[[int(a) for a in row] for row in transform.table()],
                reachable_times=reachable, horizons=horizons, triple_certificates=triples,
                finite_carriers=finite, pell_primitive_smith=smith(pell_primitive), pell_even_smith=smith(pell_even),
                pell_mod8_witness=z, prime_controls=prime_controls, field11_order=266,
                field8_points=field8, paley_automorphisms=len(automorphisms), paley_hamilton_paths=hamilton)


if __name__ == '__main__':
    import sys
    result = json.loads(json.dumps(run()))
    path = HERE/'clocks.json'
    if '--write' in sys.argv:
        path.write_text(json.dumps(result, separators=(',',':'))+'\n', encoding='utf-8')
    else:
        need(result == json.loads(path.read_text(encoding='utf-8')), 'unchanged exact clock certificate')
    print('PASS cubic/Pell clock lattices; horizon44; triple48; prime11 dual; F8 orientation')
