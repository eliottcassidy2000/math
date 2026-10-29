"""Exact, source-named checks for tetrahedral gluing and unit squaring clocks.

This is a same-agent independent-algorithm check, not an independent referee.
All truth gates remain active under python -O.
"""
from collections import Counter
from itertools import combinations, permutations
from math import gcd, prod
from pathlib import Path
import argparse
import hashlib
import json
import sys
from sympy import Matrix, Poly, eye, factor, factorint, mobius, symbols, zeros
from sympy.ntheory.modular import crt

HERE = Path(__file__).resolve().parent
OLD = HERE.parent/'overnight_research_20261029'
sys.path.insert(0, str(OLD))
from three_run_glue import (FOURIER, three_runs, split_operator, diagonal,
                            local_layers, rref_mod, hostile, automorphisms)

GATES = 0
X = symbols('x')


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def mod(A, p):
    return A.applyfunc(lambda x: int(x) % p)


def rank(A, p):
    return len(rref_mod(A.tolist(), p)[1])


def polynomial(A, coeffs):
    out = zeros(A.rows)
    for c in reversed(coeffs):
        out = out*A+int(c)*eye(A.rows)
    return out


def coeffs(f):
    return list(reversed([int(c) for c in Poly(f, X).all_coeffs()]))


def right_inverse(H, p):
    piv = rref_mod(H.tolist(), p)[1]
    J = zeros(H.cols, H.rows)
    V = H[:, piv].inv_mod(p)
    for i, c in enumerate(piv):
        J[c, :] = V[i, :]
    need(mod(H*J, p) == eye(H.rows), 'quotient coordinate section')
    return J


def quotient_clock(Q, U, p):
    red, piv = rref_mod(Q.T.tolist(), p)
    free = [i for i in range(Q.rows) if i not in piv]
    H = zeros(len(free), Q.rows)
    for j, f in enumerate(free):
        H[j, f] = 1
        for i, c in enumerate(piv):
            H[j, c] = -red[i][f] % p
    J = right_inverse(H, p)
    C = mod(H*U*J, p)
    need(mod(H*Q, p) == zeros(H.rows, Q.cols), 'literal cokernel observer')
    need(mod(H*U-C*H, p) == zeros(H.rows, Q.rows), 'named clock on cokernel')
    return C


def t5():
    rows = [sum(1 << ((i+j) % 5) for j in (1, 2)) for i in range(5)]
    triples, ports, ml = three_runs(rows)
    P, N, W = map(Matrix, split_operator(ml, triples))
    M = Matrix(ml)
    F = Matrix.diag(*([Matrix(FOURIER)]*10))
    A = F*M*F.inv()
    need(all(x.q == 1 for x in A), 'integral sector operator')
    # Five oriented facets of the four-simplex, with their boundary signs.
    B = zeros(10, 5)
    for omitted in range(5):
        tetra = tuple(v for v in range(5) if v != omitted)
        for j in range(4):
            face = tetra[:j]+tetra[j+1:]
            B[triples.index(face), omitted] = (-1)**(omitted+j)
    K = B[:, :4]
    need(B*Matrix([1]*5) == zeros(10, 1), 'the five tetrahedral boundaries have one relation')
    need(P*B == zeros(10, 5) and B.T*P == zeros(5, 10), 'exact right and left tetrahedral kernels')
    need(P.rank() == 6 and K.rank() == 4, 'kernel exhausted')
    need(diagonal(K.tolist()) == [1]*4, 'tetrahedral kernel lattice is primitive')
    need(K.T*K == 5*eye(4)-Matrix.ones(4), 'kernel Gram matrix')
    need(diagonal((K.T*K).tolist()) == [1, 5, 5, 5], 'secondary five-primary clutch')
    Iplus = zeros(60, 10)
    for i in range(10):
        Iplus[6*i, i] = 1
    KF = Iplus*K
    need(A*KF == zeros(60, 4), 'kernel located in unsigned sector')
    native_kernel = F.inv()*KF*6
    need(all(native_kernel[6*i+j, c] == K[i, c] for i in range(10)
             for j in range(6) for c in range(4)), 'native kernel repeats coefficients on all run names')
    need(diagonal(native_kernel.tolist()) == [1]*4, 'native kernel is a primitive lattice')
    need(M.T*native_kernel == zeros(60, 4), 'native left and right kernels coincide')
    need(native_kernel.T*native_kernel == 6*K.T*K, 'six run markings scale the kernel Gram matrix')
    # Literal named vertex clock, transported through the chart.
    V = zeros(60)
    for j, port in enumerate(ports):
        V[ports.index(tuple((v+1) % 5 for v in port)), j] = 1
    U = F*V*F.inv()
    need(V*M == M*V and U*A == A*U, 'native order-five clock commutes with response')
    C = zeros(4)
    for j in range(3):
        C[j+1, j] = 1
    C[:, 3] = Matrix([-1]*4)
    need(U*KF == KF*C, 'tetrahedral clock keeps the facet phase')
    need(C**5 == eye(4) and coeffs(C.charpoly(X).as_expr()) == [1]*5,
         'augmentation quotient has cyclotomic clock Phi5')
    quotient_data = []
    for p, h in [(2, [[1, -1, 0, 0, 0, 0]]),
                 (3, [[-1, 0, 1, 1, 0, 0], [-1, 0, 0, 0, 1, 1],
                      [0, -1, 1, 0, -1, 0]])]:
        h = Matrix(h)
        j = right_inverse(h, p)
        H, J = Matrix.diag(*([h]*10)), Matrix.diag(*([j]*10))
        Q, R = mod(H*A*J, p), mod(H*U*J, p)
        Z = mod(H*KF, p)
        need(mod(H*A-Q*H, p) == zeros(H.rows, 60), 'entire gluing quotient intertwines')
        need(mod(H*U-R*H, p) == zeros(H.rows, 60), 'quotient clock intertwines')
        need(rank(Z, p) == 4 and mod(Q*Z, p) == zeros(Q.rows, 4),
             'kernel maps onto all four quotient kernel coordinates')
        need(rank(Q, p) == Q.rows-4 and rank(Q*Q, p) == Q.rows-4,
             'quotient zero eigenspace is semisimple: all powers covered')
        D = quotient_clock(Q, R, p)
        need(Poly(D.charpoly(X).as_expr(), X, modulus=p) == Poly(sum(X**i for i in range(5)), X, modulus=p),
             'cokernel retains the same order-five character')
        need(Poly(sum(X**i for i in range(5)), X, modulus=p).is_irreducible,
             'four-dimensional residue module is irreducible')
        quotient_data.append({'prime':p, 'matrix':Q.tolist(), 'clock':R.tolist(),
                              'kernel_embedding':Z.tolist(), 'cokernel_clock':D.tolist()})
    polynomials = [P.charpoly(X).as_expr()/X**3, N.charpoly(X).as_expr(), W.charpoly(X).as_expr()]
    # These polynomials are squarefree and pairwise coprime. For P, the
    # sole repeated characteristic root is zero and rank(P^2)=rank(P).
    need(P.rank() == (P*P).rank(), 'native zero eigenspace is semisimple')
    for f in polynomials:
        need(Poly(f, X).gcd(Poly(f, X).diff()).degree() == 0, 'sector minimum polynomial squarefree')
    for f, g in combinations(polynomials, 2):
        need(Poly(f, X).gcd(Poly(g, X)).degree() == 0, 'sector spectra disjoint')
    minimum = Poly(prod(polynomials), X)
    need(minimum.degree() == 37, 'exact universal matrix recurrence degree')
    need(polynomial(M, coeffs(minimum.as_expr())) == zeros(60), 'integer recurrence annihilates literal native matrix')
    unit = Poly(minimum.as_expr()/X, X)
    pseudo = abs(int(unit.nth(0)))*abs(int(W.det()))  # W occurs twice.
    need(pseudo == 52488000 == 2**6*3**8*5**3, 'complete nonzero characteristic product')
    need(pseudo == 324*int((native_kernel.T*native_kernel).det()), 'independent kernel-volume torsion identity')
    need(diagonal(ml) == [1]*52+[3]*2+[6]*2+[0]*4, 'inherited native Smith form')
    need(diagonal(P.tolist()) == [1]*6+[0]*4, 'unsigned response surjects its saturated image S')
    need(all(x % 25 == 0 for x in P**5), 'four-step five-primary lattice return is integral')
    need(diagonal((P**5/25).tolist()) == [1]*6+[0]*4,
         'P^4 S = 25 S as actual integer lattices')
    slope_certificates = []
    for name, B0, p, d in [('minus', N, 2, 4), ('standard', W, 2, 1), ('standard', W, 3, 4)]:
        f = Poly(B0.charpoly(X).as_expr(), X)
        zero_order = next(i for i in range(B0.rows+1) if int(f.nth(i)) % p)
        det, v = abs(int(B0.det())), 0
        while det % p == 0:
            det //= p
            v += 1
        need(zero_order == B0.rows-rank(B0, p) == v == d,
             'nonunit Fitting block is exactly p times an integral automorphism')
        slope_certificates.append({'sector':name, 'prime':p, 'dimension':d})
    five_bases = [[0]*6, [0]*3+[1]*3, [0]+[1]*4+[2], [1]*3+[2]*3]
    power_data = []
    for r in [*range(1, 10), 13, 17, 25]:
        Ap, Mp = A**r, M**r
        item = {'power':r, 'native':{}, 'unglued':{}}
        for p, expected in [(2, 6*r-4), (3, 8*r-4), (5, 3*r-3)]:
            precision = expected+6
            a = local_layers(Mp.tolist(), p, precision)
            b = local_layers(Ap.tolist(), p, precision)
            need(sum(a) == sum(b) == 56, 'all finite primary factors captured')
            need(sum(i*c for i, c in enumerate(a)) == expected, 'proved all-power torsion formula finite probe')
            need(sum(i*(v-u) for i, (u, v) in enumerate(zip(a, b))) == (4 if p in (2, 3) else 0),
                 'all-power gluing index finite probe')
            if p == 2:
                predicted = [0]*50+[r-1]*4+[r]*2
                unglued = [0]*50+[r]*6
            elif p == 3:
                predicted = [0]*48+[r-1]*4+[r]*4
                unglued = [0]*48+[r]*8
            else:
                q, s = divmod(r-1, 4)
                predicted = [0]*50+[e+2*q for e in five_bases[s]]
                unglued = predicted
            actual = [i for i, c in enumerate(a) for _ in range(c)]
            other = [i for i, c in enumerate(b) for _ in range(c)]
            need(actual == sorted(predicted) and other == sorted(unglued),
                 'complete proved all-power primary Smith profiles')
            for h in (1, r-1, r, r+1):
                if h < 1:
                    continue
                gap = sum(max(h-e, 0) for e in actual)-sum(max(h-e, 0) for e in other)
                need(gap == (4 if p in (2, 3) and h >= r else 0), 'exact precision of gluing visibility')
            item['native'][str(p)] = {str(i):c for i, c in enumerate(a) if c}
            item['unglued'][str(p)] = {str(i):c for i, c in enumerate(b) if c}
        power_data.append(item)
    return {'rows':rows, 'ports':ports, 'tetrahedral_boundaries':B.tolist(),
            'kernel_clock':C.tolist(), 'gluing_quotients':quotient_data,
            'sector_minimum_polynomials_low_first':[coeffs(f) for f in polynomials],
            'minimum_polynomial_low_first':coeffs(minimum.as_expr()),
            'slope_one_certificates':slope_certificates, 'five_primary_base_profiles':five_bases,
            'torsion_initial':324, 'torsion_multiplier':pseudo, 'powers':power_data,
            'matrix_sha256':hashlib.sha256(json.dumps(ml, separators=(',', ':')).encode()).hexdigest()}


def four_vertices():
    inherited = json.loads((OLD/'three_run_glue.json').read_text())['order4_complete_census']
    data = []
    for entry in inherited:
        rows = entry['rows']
        _, ports, ml = three_runs(rows)
        M, vacancy = Matrix(ml), zeros(24)
        for j, port in enumerate(ports):
            missing = next(v for v in range(4) if v not in port)
            for i, v in enumerate(port):
                target = list(port)
                target[i] = missing
                if rows[v] >> missing & 1:
                    vacancy[j, ports.index(tuple(target))] = 1
        need(M == vacancy, 'independent one-edge/vacancy model')
        need(sum(M) == 36, 'six edges times six marked lifts')
        need(all(sum(row) == 3 for row in (M+M.T).tolist()), 'star-transposition undirected lift is three-regular')
        need(diagonal(ml) == entry['full_smith'], 'the four inherited class controls')
        data.append({'rows':rows, 'rank':M.rank(), 'mass':int(sum(M))})
    # Distinct six-dimensional integral objects, from recovered canon.
    edges = list(combinations(range(4), 2))
    ps, ns = [], []
    for i, j in [(0, 5), (1, 4), (2, 3)]:
        p, n = zeros(6, 1), zeros(6, 1)
        p[i] = p[j] = 1
        n[i], n[j] = 1, -1
        ps.append(p); ns.append(n)
    unsigned = Matrix.hstack(Matrix.ones(6, 1), ps[0]-ps[2], ps[1]-ps[2], *ns)
    inc = zeros(4, 6)
    for j, (a, b) in enumerate(edges):
        inc[a, j], inc[b, j] = -1, 1
    cycles = []
    for a, b in combinations(range(1, 4), 2):
        c = zeros(6, 1)
        c[edges.index((0, a))] = c[edges.index((a, b))] = 1
        c[edges.index((0, b))] = -1
        cycles.append(c)
    signed = Matrix.hstack(inc[:3, :].T, *cycles)
    need(diagonal(unsigned.tolist()) == [1, 1, 1, 2, 2, 6], 'unsigned K4 clutch, THM-3045')
    need(diagonal(signed.tolist()) == [1, 1, 1, 1, 4, 4], 'orientation cut/cycle clutch, K(K4)')
    return {'representatives':data, 'unsigned_smith':diagonal(unsigned.tolist()),
            'signed_smith':diagonal(signed.tolist()), 'run_chart_smith':diagonal(FOURIER)}


def squaring():
    # Counts tail <= t and period dividing r, all units, odd modulus.
    def count(n, t, r):
        return prod(gcd((p-1)*p**(e-1), 2**t*(2**r-1)) for p, e in factorint(n).items())
    tests = 0
    for n in range(3, 152, 2):
        units = [x for x in range(n) if gcd(x, n) == 1]
        signatures = {}
        for x in units:
            seen, y = {}, x
            while y not in seen:
                seen[y] = len(seen)
                y = y*y % n
            signatures[x] = (seen[y], len(seen)-seen[y])
        for t in range(6):
            for r in range(1, 13):
                actual = sum(a <= t and r % b == 0 for a, b in signatures.values())
                need(actual == count(n, t, r), 'full CRT fixed-iterate count versus literal walks')
                exact = sum(int(mobius(r//d))*count(n, t, d) for d in range(1, r+1) if r % d == 0)
                need(exact == sum(a <= t and b == r for a, b in signatures.values()), 'Mobius exact-period recovery')
                tests += 1
        phis = [(p-1)*p**(e-1) for p, e in factorint(n).items()]
        alphas = [(v & -v).bit_length()-1 for v in phis]
        for t in range(6):
            need(count(n, t+1, 7) == count(n, t, 7)*2**sum(a > t for a in alphas), 'all-height branching recurrence')
    # Two local two-cycles require their relative phase; four states do not form one cycle.
    phases = [[int(crt([7, 13], residues)[0]) for residues in pair]
              for pair in [[(2, 3), (4, 9)], [(2, 9), (4, 3)]]]
    for a, b in phases:
        need(a*a % 91 == b and b*b % 91 == a, 'explicit relative-phase hostile')
    F5 = 2**32+1
    q = 6700417
    need(F5 == 641*q, 'Euler factorization positive control')
    local = next(pow(a, (q-1)//3, q) for a in range(2, 50) if pow(a, (q-1)//3, q) != 1)
    witness = int(crt([641, q], [1, local])[0])
    need(witness != 1 and pow(witness, 3, F5) == 1, 'explicit odd-order compositeness certificate')
    cycle = [witness, pow(witness, 2, F5)]
    need(cycle[0] != cycle[1] and cycle[1]**2 % F5 == cycle[0], 'two-cycle under squaring')
    for n in range(11):
        f = 2**(2**n)+1
        need(pow(2, 2**n, f) == f-1 and pow(2, 2**(n+1), f) == 1,
             'base-two source is blind to primality, prime and composite controls')
    return {'odd_moduli_inclusive':[3, 151], 'tail_period_probes':tests,
            'relative_phase_cycles_mod91':phases, 'F5_odd_order_witness':witness,
            'F5_squaring_cycle':cycle, 'status':'No Fermat completeness conclusion'}


def hostile_control():
    cyc = []
    old = json.loads((OLD.parent/'overnight_rank_return_20261029/two_run.json').read_text())
    ranks = []
    for flip, expected in [(False, 81), (True, 27)]:
        rows = hostile(flip)
        need(len(automorphisms(rows)) == expected, 'native hostile automorphism group')
        cyc.append([t for t in combinations(range(12), 3)
                    if all(sum((rows[v] >> w) & 1 for w in t) == 1 for v in t)])
        entry = old['objects'][-1 if flip else -2]
        need(rows == entry['rows'], 'hostile source labels match saved matrix')
        ranks.append(rank(Matrix(entry['response']), 3))
    need(cyc[0] == cyc[1] and len(cyc[0]) == 67, 'same cyclic-triple data do not identify the sources')
    need(ranks == [24, 22], 'port-sensitive hostile selector still separates')
    return {'automorphisms':[81, 27], 'common_labelled_cyclic_triples':67, 'two_run_mod3_ranks':ranks}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=HERE/'certificate.json')
    args = ap.parse_args()
    out = {'t5':t5(), 'four_vertices':four_vertices(), 'squaring':squaring(), 'hostile':hostile_control()}
    out['gates'] = GATES
    args.output.write_text(json.dumps(out, indent=2, default=int)+'\n', encoding='utf-8')
    print(json.dumps({'gates':GATES, 'output':str(args.output), 'sha256':hashlib.sha256(args.output.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
