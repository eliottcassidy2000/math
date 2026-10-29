"""Exact signed Paley pencils, including ramification and common content."""
import argparse
from collections import Counter
import json
from math import gcd, prod
from pathlib import Path

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise ValueError(why)


def diagonal(matrix):
    return [abs(int(x)) for x in smith_normal_form(Matrix(matrix), domain=ZZ).diagonal()]


def valuation(value, prime):
    need(value != 0, 'valuation requires a nonzero integer')
    result = 0
    while value % prime == 0:
        result += 1
        value //= prime
    return result


def local_layers(matrix, prime, precision):
    """Independent elimination over Z/(p^k); counts valuations 0,...,k-1.

    Unit row/column operations isolate every unit pivot. The residual block
    is divisible by p; divide it and recurse at one less digit. No Smith
    routine, rational rank routine, or predicted invariant is used here.
    """
    modulus = prime**precision
    a = [[x % modulus for x in row] for row in matrix]
    nr, nc, rank = len(a), len(a[0]), 0
    while rank < min(nr, nc):
        pivot = next(((i, j) for i in range(rank, nr) for j in range(rank, nc)
                      if a[i][j] % prime), None)
        if pivot is None:
            break
        i, j = pivot
        a[rank], a[i] = a[i], a[rank]
        for row in a:
            row[rank], row[j] = row[j], row[rank]
        inverse = pow(a[rank][rank], -1, modulus)
        a[rank] = [(x*inverse) % modulus for x in a[rank]]
        for i in range(rank+1, nr):
            coefficient = a[i][rank]
            a[i] = [(x-coefficient*y) % modulus for x, y in zip(a[i], a[rank])]
        # Column elimination changes only the pivot row: the pivot column
        # below it is already zero. The residual block is unchanged.
        rank += 1
    if precision == 1:
        return [rank]
    if rank == min(nr, nc):
        return [rank]+[0]*(precision-1)
    residual = [[a[i][j]//prime for j in range(rank, nc)] for i in range(rank, nr)]
    need(all(a[i][j] % prime == 0 for i in range(rank, nr) for j in range(rank, nc)),
         'residual block divisible by the prime')
    return [rank]+local_layers(residual, prime, precision-1)


def paley(q):
    squares = {x*x % q for x in range(1, q)}
    return [[int((j-i) % q in squares) for j in range(q)] for i in range(q)]


def predicted(q, a, b):
    r, content = (q-1)//2, gcd(a, b)
    if not content:
        return [0]*q
    x, y = a//content, b//content
    s = x+y
    if not s:
        return [content]*(r+1)+[content*q]*(r-1)+[0]
    norm = (s*s+q*(x-y)**2)//4
    scalar = abs(r*s)
    d = gcd(scalar, norm)//gcd(q, s)
    return [content]*r+[content*d]+[content*norm]*(r-1)+[content*scalar*norm//d]


def field27_control():
    # F3[z]/(z^3-z-1); the defining cubic has no root in F3.
    def digits(a):
        return [(a//3**i) % 3 for i in range(3)]
    def number(a):
        return sum((x % 3)*3**i for i, x in enumerate(a))
    def sub(a, b):
        return number([x-y for x, y in zip(digits(a), digits(b))])
    def mul(a, b):
        out = [0]*5
        for i, x in enumerate(digits(a)):
            for j, y in enumerate(digits(b)):
                out[i+j] += x*y
        for k in (4, 3):
            out[k-2] += out[k]
            out[k-3] += out[k]
        return number(out[:3])
    need(all((x**3-x-1) % 3 for x in range(3)), 'irreducible field cubic')
    need(all(any(mul(x, y) == 1 for y in range(1, 27)) for x in range(1, 27)),
         'every nonzero field element has an inverse')
    squares = {mul(x, x) for x in range(1, 27)}
    A = [[int(sub(j, i) in squares) for j in range(27)] for i in range(27)]
    need(all(A[i][j]+A[j][i] == int(i != j) for i in range(27) for j in range(27)),
         'actual Paley tournament over F27')
    records = []
    for a, b in ((1, 2), (1, 26), (1, -1)):
        E = [[a*A[i][j]+b*A[j][i] for j in range(27)] for i in range(27)]
        ds = diagonal(E)
        layers = local_layers(E, 3, 6)
        counts = Counter(valuation(d, 3) for d in ds if d)
        need(layers == [counts[k] for k in range(6)], 'F27 local check')
        need(layers[0] == 8, 'Lucas digit rank 2^3')
        need(ds != predicted(27, a, b), 'prime-only formula has a prime-power hostile')
        records.append(dict(a=a, b=b, smith=ds, three_layers=layers))
    return dict(field_polynomial='z^3-z-1 over F3', matrices=records)


def run():
    census, records = Counter(), []
    for q in (3, 7, 11, 19, 23, 31):
        A, r = paley(q), (q-1)//2
        Q = [[A[i][j]-A[j][i] for j in range(q)] for i in range(q)]
        need(all(sum(Q[i][k]*Q[k][j] for k in range(q)) == 1-q*int(i == j)
                 for i in range(q) for j in range(q)), 'quadratic character identity')
        need(local_layers(Q, q, 1) == [r+1], 'Vandermonde character rank')
        weights = {(a, b) for a in range(-4, 5) for b in range(-4, 5)}
        for t in range(1, 5):
            weights.update([(1, q**t-1), (2, q**t-2), (3, 3*(q**t-1))])
        for a, b in sorted(weights):
            E = [[a*A[i][j]+b*A[j][i] for j in range(q)] for i in range(q)]
            actual = diagonal(E)
            need(actual == predicted(q, a, b), f'Smith formula q={q}, a={a}, b={b}')
            nonzero = [d for d in actual if d]
            precision = 1+max((valuation(d, q) for d in nonzero), default=0)
            observed = local_layers(E, q, precision)
            counts = Counter(valuation(d, q) for d in nonzero)
            need(observed == [counts[k] for k in range(precision)], 'local independent elimination')
            g = gcd(a, b)
            kind = 'zero' if not g else 'singular' if a+b == 0 else 'ramified' if (a+b)//g % q == 0 else 'unramified'
            census[kind] += 1
            if a+b:
                s, delta = (a+b)//g, (a-b)//g
                norm = (s*s+q*delta*delta)//4
                need(prod(actual) == g**q*abs(r*s)*norm**r, 'determinant from spectral factors')
            records.append(dict(q=q, a=a, b=b, kind=kind, smith=actual, q_layers=observed))
    # The first failure of extending the old formula outside its q guard.
    need(predicted(3, 1, 2) == [1, 1, 9], 'minimal ramification witness')
    # Identical rational spectrum and determinant do not justify Z/3 + Z/3.
    return dict(universe='q in {3,7,11,19,23,31}; signed [-4,4]^2 plus three q^t probes for t=1..4',
                counts=dict(sorted(census.items())), cases=records,
                naive_ramified_diagonal=[1, 3, 3], actual_ramified_diagonal=[1, 1, 9],
                prime_power_hostile=field27_control())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE/'paley_pencil.json'
    if args.write:
        target.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    else:
        need(result == json.loads(target.read_text(encoding='utf-8')), 'receipt replay')
    print(json.dumps(dict(status='PASS', counts=result['counts'], cases=len(result['cases']))))


if __name__ == '__main__':
    main()
