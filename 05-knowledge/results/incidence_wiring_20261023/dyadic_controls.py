"""Test dyadic lift/rank rules against independent exact determinantal divisors.

The small-matrix oracle uses all minors and integer determinants, not the
binary lifting procedure. This is a same-agent independent algorithm check.
"""
from itertools import combinations
from pathlib import Path
import argparse, json, math, random
import check as core

HERE = Path(__file__).resolve().parent
need = core.need

def det(M):
    if not M: return 1
    A = [r[:] for r in M]; sign = 1; last = 1; n = len(A)
    for k in range(n-1):
        if not A[k][k]:
            p = next((i for i in range(k+1, n) if A[i][k]), None)
            if p is None: return 0
            A[k], A[p] = A[p], A[k]; sign = -sign
        pivot = A[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = A[i][j]*pivot-A[i][k]*A[k][j]
                need(numerator % last == 0, 'fraction-free determinant division')
                A[i][j] = numerator//last
            A[i][k] = 0
        last = pivot
    return sign*A[-1][-1]

def smith_valuations(A):
    previous = 1; vals = []
    for k in range(1, len(A)+1):
        divisor = 0
        for I in combinations(range(len(A)), k):
            for J in combinations(range(len(A)), k):
                divisor = math.gcd(divisor, det([[A[i][j] for j in J] for i in I]))
        need(divisor and divisor % previous == 0, 'determinantal divisor chain')
        s = divisor//previous; previous = divisor; v = 0
        while s % 2 == 0: s //= 2; v += 1
        vals.append(v)
    return vals

def dyadic_ranks(A):
    n = len(A)
    rows = [sum((a&1) << j for j, a in enumerate(row)) for row in A]
    V, _ = core.binary_kernel(rows, n)
    # Exhaustively audit the kernel dimension independently at these small sizes.
    brute = [v for v in range(1 << n) if all(not((v&r).bit_count()&1) for r in rows)]
    need(len(brute) == 1 << len(V), 'binary kernel misses solutions')
    vectors = [[(v >> i)&1 for i in range(n)] for v in V]
    def act(v): return [sum(a*b for a, b in zip(row, v)) for row in A]
    def form(vv, divisor):
        AV = list(map(act, vv)); M = []
        for v in vv:
            row = []
            for w in AV:
                z = sum(a*b for a, b in zip(v, w)); need(z % divisor == 0, 'form integrality')
                row.append((z//divisor)&1)
            M.append(row)
        return [sum(b << j for j, b in enumerate(row)) for row in M]
    Q2 = form(vectors, 2); R, rank2 = core.binary_kernel(Q2, len(V))
    lift4 = []
    for r in R:
        x = [sum(vectors[j][i] for j in core.bits(r))&1 for i in range(n)]
        Ax = act(x); y = core.binary_solve(rows, [(a//2)&1 for a in Ax], n)
        v = [a+2*((y >> i)&1) for i, a in enumerate(x)]
        need(all(a % 4 == 0 for a in act(v)), 'lift4 equation')
        lift4.append(v)
    Q4 = form(lift4, 4); S, rank4 = core.binary_kernel(Q4, len(R))
    lift8 = []
    for s in S:
        v = [sum(lift4[j][i] for j in core.bits(s)) for i in range(n)]
        Av = act(v)
        obstruction = [sum(a*b for a, b in zip(x, Av))//4 & 1 for x in vectors]
        coeff = core.binary_solve(Q2, obstruction, len(V))
        w = [sum(vectors[j][i] for j in core.bits(coeff))&1 for i in range(n)]
        v = [a+2*b for a, b in zip(v, w)]; Av = act(v)
        y = core.binary_solve(rows, [(a//4)&1 for a in Av], n)
        v = [a+4*((y >> i)&1) for i, a in enumerate(v)]
        need(all(a % 8 == 0 for a in act(v)), 'lift8 equation')
        lift8.append(v)
    Q8 = form(lift8, 8); rank8 = core.binary_rank(Q8)
    stages = {1: (vectors, Q2), 2: (lift4, Q4), 3: (lift8, Q8)}
    histogram = {1: rank2, 2: rank4, 3: rank8}
    for k in range(3, 13):
        old, q = stages[k]; surviving, _ = core.binary_kernel(q, len(old))
        if not surviving: break
        following = []
        for r in surviving:
            v = [sum(old[t][i] for t in core.bits(r)) for i in range(n)]; Av = act(v)
            for j in range(k-1, 0, -1):
                basis, Q = stages[j]
                rhs = [(sum(a*b for a, b in zip(w, Av))//(1 << k))&1 for w in basis]
                coeff = core.binary_solve(Q, rhs, len(basis))
                w = [sum(basis[t][i] for t in core.bits(coeff)) for i in range(n)]
                v = [a+(1 << (k-j))*b for a, b in zip(v, w)]; Av = act(v)
                need(all(a % (1 << k) == 0 for a in Av), 'general nested divisibility')
            y = core.binary_solve(rows, [(a//(1 << k))&1 for a in Av], n)
            v = [a+(1 << k)*((y >> i)&1) for i, a in enumerate(v)]
            need(all(a % (1 << (k+1)) == 0 for a in act(v)), 'general full lift')
            following.append(v)
        qnext = form(following, 1 << (k+1))
        stages[k+1] = following, qnext; histogram[k+1] = core.binary_rank(qnext)
    return [len(V), rank2, len(R), rank4, len(S), rank8, len(S)-rank8], {k: v for k, v in histogram.items() if v}

def main():
    p = argparse.ArgumentParser(); p.add_argument('--write', action='store_true'); args = p.parse_args()
    rng = random.Random(238413); cases = []
    for n in range(1, 6):
        for _ in range(12):
            B = [[int(i == j) for j in range(n)] for i in range(n)]
            for _ in range(3*n):
                a, b = rng.randrange(n), rng.randrange(n)
                if a != b: B[a] = [x+y for x, y in zip(B[a], B[b])]
            D = [2**rng.randrange(6) for _ in range(n)]
            A = [[sum(D[k]*B[k][i]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
            cases.append(A)
        accepted = 0
        while accepted < 12:
            A = [[0]*n for _ in range(n)]
            scale = 2**rng.randrange(3)
            for i in range(n):
                for j in range(i, n): A[i][j] = A[j][i] = scale*rng.randrange(-5, 6)
            if det(A): cases.append(A); accepted += 1
    records = []
    for A in cases:
        vals = smith_valuations(A)
        expected = [sum(v >= 1 for v in vals), sum(v == 1 for v in vals), sum(v >= 2 for v in vals),
                    sum(v == 2 for v in vals), sum(v >= 3 for v in vals), sum(v == 3 for v in vals), sum(v >= 4 for v in vals)]
        actual, histogram = dyadic_ranks(A)
        need(actual == expected, 'dyadic ranks disagree with independent minors')
        need(histogram == {k: vals.count(k) for k in set(vals) if k > 0}, 'complete nested lift disagrees with Smith factors')
        records.append({'matrix': A, 'smith_valuations': vals, 'dyadic_counts': actual, 'positive_valuation_histogram': histogram})
    out = {'scope': '120 nonsingular symmetric integer matrices, orders1..5, congruent diagonal and general signed controls',
           'method': 'all-minor determinantal divisors versus nested binary kernel lifts until no radical remains',
           'checks': core.CHECKS, 'cases': records}
    out = json.loads(json.dumps(out)); dest = HERE/'dyadic_controls.json'
    if args.write: dest.write_text(json.dumps(out, indent=2)+'\n')
    else: need(json.loads(dest.read_text()) == out, 'dyadic control receipt differs')
    print('exact dyadic controls:', len(records), 'checks:', out['checks'])

if __name__ == '__main__': main()
