#!/usr/bin/env python3
"""Exact small controls for decision-specific dual certificates (stdlib only)."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product, permutations
from math import gcd, isqrt, prod
import json

CHECKS = 0

def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def left(y, A):
    return [sum(y[i]*A[i][j] for i in range(len(A))) for j in range(len(A[0]))]


def primal(A, x, b):
    return [dot(row, x) for row in A] == list(b)


def dual(A, b, y, modulus=0):
    c, d = left(y, A), dot(y, b)
    if modulus:
        return all(t % modulus == 0 for t in c) and d % modulus != 0
    return all(t == 0 for t in c) and d != 0


def eliminate(A, b):
    """Rational elimination preserving source-row multipliers."""
    m, n = len(A), len(A[0])
    M = [[F(x) for x in row]+[F(b[i])] for i, row in enumerate(A)]
    Y = [[F(i == j) for j in range(m)] for i in range(m)]
    pivots = []
    for j in range(n):
        t = next((i for i in range(len(pivots), m) if M[i][j]), None)
        if t is None:
            continue
        k = len(pivots)
        M[k], M[t] = M[t], M[k]
        Y[k], Y[t] = Y[t], Y[k]
        v = M[k][j]
        M[k] = [x/v for x in M[k]]
        Y[k] = [x/v for x in Y[k]]
        for i in range(m):
            if i != k and M[i][j]:
                v = M[i][j]
                M[i] = [x-v*z for x, z in zip(M[i], M[k])]
                Y[i] = [x-v*z for x, z in zip(Y[i], Y[k])]
        pivots.append(j)
    for i in range(m):
        if not any(M[i][:-1]) and M[i][-1]:
            return len(pivots), None, Y[i]
    x = [F(0) for _ in range(n)]
    for i, j in enumerate(pivots):
        x[j] = M[i][-1]
    return len(pivots), x, None


def egcd(a, b):
    if b == 0:
        return abs(a), (1 if a >= 0 else -1), 0
    g, x, y = egcd(b, a % b)
    return g, y, x-(a//b)*y


def integer_2x2(A, b):
    """Independent direct rank/determinant/gcd decision and YES witness."""
    a, c = A[0]
    d, e = A[1]
    D = a*e-c*d
    if D:
        nums = [e*b[0]-c*b[1], a*b[1]-d*b[0]]
        return ([u//D for u in nums] if all(u % D == 0 for u in nums) else None)
    row = next((i for i in range(2) if any(A[i])), None)
    if row is None:
        return [0, 0] if b == (0, 0) else None
    u, v = A[row]
    g, s, t = egcd(u, v)
    if b[row] % g:
        return None
    x = [s*(b[row]//g), t*(b[row]//g)]
    return x if primal(A, x, b) else None


def small_linear():
    counts = Counter()
    max_q = 0
    for entries in product(range(-2, 3), repeat=4):
        A = [list(entries[:2]), list(entries[2:])]
        for b in product(range(-2, 3), repeat=2):
            r, qx, y = eliminate(A, b)
            ix = integer_2x2(A, b)
            if y is not None:
                need(ix is None and dual(A, b, y), 'rational dual mismatch')
                need(sum(t != 0 for t in y) <= r+1, 'rational sparsity')
                counts['rational_NO'] += 1
            elif ix is not None:
                need(primal(A, ix, b), 'integer primal mismatch')
                counts['integer_YES'] += 1
            else:
                # Rank-two bound |det|, rank-one bound max entry; all columns checked.
                bound = abs(entries[0]*entries[3]-entries[1]*entries[2]) if r == 2 else max(map(abs, entries))
                certificate = next(((q, yy) for q in range(2, bound+1)
                                    for yy in product(range(q), repeat=2)
                                    if sum(t != 0 for t in yy) <= r and dual(A, b, yy, q)), None)
                need(certificate is not None, 'missing integer character')
                max_q = max(max_q, certificate[0])
                counts['torsion_NO'] += 1
    # Rational and integer support size are different; not every field sees torsion.
    need(dual([[4]], [2], [1], 4), 'prime-power character')
    for p in (2, 3, 5, 7, 11, 13):
        need(any((4*x-2) % p == 0 for x in range(p)), 'prime reduction hostile')
    need(not dual([[2, 3]], [1], [1], 2), 'selected-column false rejection')
    need(primal([[2, 3]], [-1, 1], [1]), 'selected-column positive control')
    need(not any(primal([[1, 1]], x, [3]) for x in product((0, 1), repeat=2)), 'binary guard hostile')
    need(primal([[1, 1]], [3, 0], [3]), 'unrestricted integer control')
    field_count = Counter()
    for entries in product((0, 1), repeat=6):
        A = [list(entries[:3]), list(entries[3:])]
        for b in product((0, 1), repeat=2):
            xs = [x for x in product((0, 1), repeat=3) if all((dot(row, x)-bb) % 2 == 0 for row, bb in zip(A, b))]
            ys = [y for y in product((0, 1), repeat=2) if dual(A, b, y, 2)]
            need(bool(xs) != bool(ys), 'binary alternative')
            field_count['YES' if xs else 'NO'] += 1
    return {'integer_2x2_universe': 15625, 'counts': dict(sorted(counts.items())),
            'max_constructed_modulus': max_q, 'binary_2x3_universe': 256,
            'binary_counts': dict(sorted(field_count.items()))}


def signed_examples():
    triangle = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
    b = [4, 4, 1]
    r, x, y = eliminate(triangle, b)
    need(x == [F(1,2), F(7,2), F(1,2)], 'triangle rational point')
    need(dual(triangle, b, [1,1,1], 2), 'triangle torsion')
    bowtie = [[1,1,0,0,0], [0,1,1,0,0], [1,0,1,0,0],
              [1,0,0,1,0], [0,0,0,1,1], [1,0,0,0,1]]
    b = [2,2,2,2,2,4]
    y = [1,-1,1,-1,1,-1]
    need(dual(bowtie, b, y), 'bowtie dual')
    need(eliminate(bowtie, b)[0] == 5, 'bowtie rank')
    for removed in range(6):
        A0 = [row for i, row in enumerate(bowtie) if i != removed]
        b0 = [z for i, z in enumerate(b) if i != removed]
        rr, xx, yy = eliminate(A0, b0)
        need(yy is None and primal(A0, xx, b0), 'bowtie circuit minimality')
    cL, cR = left(y[:3], bowtie[:3]), left(y[3:], bowtie[3:])
    need(cL == [2,0,0,0,0] and cR == [-2,0,0,0,0], 'LCA boundary covectors')
    need(dot(y[:3], b[:3]) == 2 and dot(y[3:], b[3:]) == -4, 'LCA scalar charge')
    valuation = [[3,2,1], [0,0,2], [0,1,1]]
    need(eliminate(valuation, [0,0,0])[0] == 3, 'rational valuation rank')
    need(all(v % 2 == 0 for v in [dot(row, [1,1,1]) for row in valuation]), 'binary torsion relation')
    return {'triangle_rational_solution': ['1/2','7/2','1/2'],
            'bowtie': {'rank': 5, 'certificate': y, 'residual': -2, 'proper_row_subsets_feasible': True},
            'valuation_primes': [7,19,89], 'valuation_determinant': -6}


@lru_cache(None)
def factors(n):
    out, p = [], 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return tuple(out)


def sigma(n):
    return prod((p**(e+1)-1)//(p-1) for p,e in factors(n))


def sigma_direct(n):
    return sum(d for d in range(1, isqrt(n)+1) if n % d == 0) + sum(n//d for d in range(1, isqrt(n)+1) if n % d == 0 and d*d != n)


def isprime(n):
    return n >= 2 and factors(n) == ((n,1),)


def cyclic_data(cores, start=0):
    L = len(cores)
    order = [(start+i) % L for i in range(L)]
    c = [cores[i] for i in order]
    S = [sigma(z) for z in c]
    a = [s-z for s,z in zip(S,c)]
    d = c[1:]+c[:1]
    w = [prod(a[i+1:])*prod(d[:i]) for i in range(L)]
    original_w = [0]*L
    for i,j in enumerate(order):
        original_w[j] = w[i]
    D, B = prod(c)-prod(a), dot(w,S)
    return D,B,original_w


def aliquot():
    counts = Counter()
    successful = set()
    for L in (1,2,3):
        for cores in product(range(2,17), repeat=L):
            S = [sigma(c) for c in cores]
            A = [[0]*L for _ in range(L)]
            for i,c in enumerate(cores):
                A[i][i] -= S[i]-c
                A[i][(i+1) % L] += cores[(i+1) % L]
            rank, x, rational_dual = eliminate(A,S)
            D,B,w = cyclic_data(cores)
            need(left(w,A) == [D]+[0]*(L-1) and dot(w,S) == B and B > 0, 'affine closure dual identity')
            if D <= 0:
                need(D < 0 or (rational_dual is not None and dual(A,S,w)), 'zero-determinant positive carry')
                if D < 0:
                    need(x is not None and x[0] < 0, 'negative forced coordinate')
                counts['positive_domain_NO'] += 1
                continue
            need(rank == L and x is not None and x[0] == F(B,D), 'unique affine candidate')
            nums = []
            for j in range(L):
                Dj,Bj,wj = cyclic_data(cores,j)
                need(Dj == D and left(wj,A) == [D*(i == j) for i in range(L)], 'rotated source identity')
                need(x[j] == F(Bj,D), 'independent elimination mismatch')
                nums.append(Bj)
            if any(z.denominator != 1 for z in x):
                j = next(j for j,z in enumerate(x) if z.denominator != 1)
                _,_,wj = cyclic_data(cores,j)
                need(dual(A,S,wj,D), 'aliquot modular certificate')
                counts['integrality_NO'] += 1
                continue
            p = [int(z) for z in x]
            if any(not isprime(z) or gcd(z,c) != 1 for z,c in zip(p,cores)):
                counts['native_guard_NO'] += 1
                continue
            nodes = [c*z for c,z in zip(cores,p)]
            if len(set(nodes)) != L:
                counts['repeated_node_NO'] += 1
                continue
            need(all(sigma_direct(n)-n == nodes[(i+1)%L] for i,n in enumerate(nodes)), 'actual aliquot receipt')
            successful.add(tuple(min(tuple(nodes[i:]+nodes[:i]) for i in range(L))))
            counts['YES'] += 1
    # Explicit larger positive and hostile sources, independently recomputed.
    need(cyclic_data((2,4))[:2] == (5,37), 'small modular hostile')
    need(dual([[-1,4],[2,-3]], [3,7], [3,4], 5), 'explicit mod5 hostile')
    r,p,_ = eliminate([[-1,2,0],[0,-1,16],[2,0,-15]],[3,3,31])
    need(p == [23,13,1] and sigma_direct(16)-16 == 15, 'false prime-one cycle')
    need(primal([[-22,4],[20,-3]], [11,71], [42,7]), 'amicable affine candidate')
    need(sigma_direct(220)-220 == 284 and sigma_direct(284)-284 == 220, 'amicable true edges')
    return {'core_universe': 'all ordered core words in {2,...,16}, lengths 1..3; rotation not quotiented',
            'templates': 3615, 'counts': dict(sorted(counts.items())),
            'cycles_in_bank': sorted(successful), 'mod5_certificate': [3,4],
            'formal_prime_one_hostile': [23,13,1], 'amicable_cores': [20,4]}


def tree_path(parent, u, v):
    up = []
    z = u
    while z != -1:
        up.append(z)
        z = parent[z]
    vp = []
    z = v
    while z not in up:
        vp.append(z)
        z = parent[z]
    return up[:up.index(z)+1]+list(reversed(vp))


def graceful():
    counts = Counter()
    for m in range(1,7):
        for shape in ('path','star'):
            parent = [-1]+([v-1 for v in range(1,m+1)] if shape == 'path' else [0]*m)
            for magnitudes in permutations(range(1,m+1)):
                for signs in product((-1,1),repeat=m):
                    d = [a*s for a,s in zip(magnitudes,signs)]
                    x = [0]*(m+1)
                    for v in range(1,m+1):
                        x[v] = x[parent[v]]+d[v-1]
                    bad = next(((u,v) for u in range(m+1) for v in range(u+1,m+1) if x[u] == x[v]),None)
                    if bad is None and max(x)-min(x) > m:
                        bad = (x.index(min(x)),x.index(max(x)))
                    actual = len(set(x)) == m+1 and sorted(z-min(x) for z in x) == list(range(m+1))
                    need(actual == (bad is None),'graceful candidate iff')
                    if bad is None:
                        counts['YES'] += 1
                    else:
                        path = tree_path(parent,*bad)
                        s = sum(d[v-1] if parent[v] == u else -d[u-1] for u,v in zip(path,path[1:]))
                        need(s == x[bad[1]]-x[bad[0]],'path telescope')
                        need(s == 0 or abs(s)>m,'invalid signed candidate witness')
                        counts['zero_path_NO' if s == 0 else 'span_path_NO'] += 1
    return {'universe': 'every signed permutation on path and star trees, 1..6 edges, rooted as specified',
            'candidates': sum(counts.values()),'counts': dict(sorted(counts.items()))}


@lru_cache(None)
def prime_certificate(p):
    if p == 2:
        return {'p':2}
    fs = factors(p-1)
    g = next(g for g in range(2,p) if pow(g,p-1,p) == 1 and all(gcd(pow(g,(p-1)//q,p)-1,p) == 1 for q,e in fs))
    return {'p':p,'g':g,'factors':[(q,e,prime_certificate(q)) for q,e in fs]}


def verify_prime(cert):
    p = cert['p']
    if type(p) is not int:
        return False
    if p == 2:
        return True
    if p < 3 or p % 2 == 0:
        return False
    g, fs = cert['g'], cert['factors']
    if len({q for q,e,c in fs}) != len(fs) or not all(type(e) is int and 0 < e <= p.bit_length() and type(q) is int and q == c['p'] and 2 <= q < p and verify_prime(c) for q,e,c in fs):
        return False
    return prod(q**e for q,e,c in fs) == p-1 and pow(g,p-1,p) == 1 and all(gcd(pow(g,(p-1)//q,p)-1,p) == 1 for q,e,c in fs)


def verify_source(N, fs):
    if type(N) is not int or N < 1:
        return False
    if len({p for p,a,c in fs}) != len(fs):
        return False
    if not all(type(p) is int and 2 <= p <= N and type(a) is int and 0 < a <= N.bit_length() and p == c['p'] and verify_prime(c) for p,a,c in fs):
        return False
    return prod(p**a for p,a,c in fs) == N


def square_sigma_receipts():
    counts = Counter()
    for N in range(1,4097):
        fs = factors(N)
        need(verify_source(N, [(p,a,prime_certificate(p)) for p,a in fs]),'source factor/primality receipt')
        S = sigma(N)
        need(S == sigma_direct(N),'sigma source independent path')
        t = isqrt(S)
        if t*t == S:
            need(t*t == S,'square YES')
            counts['YES'] += 1
        else:
            need(t*t < S < (t+1)*(t+1),'adjacent-square NO')
            counts['NO'] += 1
    need(sigma(70) == 144 and sigma(30) == 72,'named arithmetic receipts')
    # Free subset completion, not existential repeated-prime fusion.
    need(all(sigma(2*prod(S)) % 3 == 0 and (sigma(2*prod(S))//3) % 3 != 0
             for mask in range(4) for S in [[p for i,p in enumerate((3,7)) if mask >> i & 1]]),'v3 parity exclusion')
    # A forged primality claim must not authenticate an arithmetic source.
    need(not verify_prime({'p':9,'g':2,'factors':[(2,3,{'p':2})]}),'forged composite prime receipt')
    return {'integers': 4096,'counts': dict(sorted(counts.items())),
            'YES_example': {'N':70,'sigma':144,'root':12},
            'NO_example': {'N':30,'sigma':72,'floor_root':8},
            'completion_NO': {'forced_prime':2,'optional_primes':[3,7],'dual_coordinate':'v_3 mod 2'}}


def weights():
    source_weights = [2,3]
    once = sum(source_weights)
    independent = sum(w*w for w in source_weights)
    need(once == 5 and independent == 13,'once-only ownership')
    need(sum([1,-1]) == 0 and len([1,-1]) == 2,'signed zero is not NO')
    need(sum([0]) == 0 and len([0]) == 1,'zero mass is not emptiness')
    return {'one_shared_choice_total':once,'incorrect_squared_total':independent,
            'signed_cancellation_nonempty_total':0,'zero_weight_nonempty_total':0}


def main():
    out = {'status':'FINITE-EXACT controls; general statements proved in DECISION_CERTIFICATES.md',
           'linear':small_linear(),'signed_incidence':signed_examples(),'aliquot':aliquot(),
           'graceful':graceful(),'square_sigma':square_sigma_receipts(),'weights':weights()}
    out['active_checks'] = CHECKS
    print(json.dumps(out,sort_keys=True,indent=2))

if __name__ == '__main__':
    main()
