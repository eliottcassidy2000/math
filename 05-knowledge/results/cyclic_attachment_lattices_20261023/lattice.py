"""Exact cyclic-chart attachment and integer-observation certificates.

The internal cyclic charts and the full named blocks are supplied inputs.
No graph-isomorphism oracle, floating point, or guessed port identification.
"""
from fractions import Fraction
from functools import lru_cache
from math import gcd, lcm, prod
from sympy import Matrix, Poly, Symbol, cyclotomic_poly
from sympy.matrices.normalforms import hermite_normal_form


def shifts(A, B):
    """All (a,b) with B[i+a,j+b] == A[i,j], on the two tori."""
    m, n = len(A), len(A[0])
    if len(B) != m or any(len(row) != n for row in A + B):
        raise ValueError('incompatible block dimensions')
    return [(a, b) for a in range(m) for b in range(n)
            if all(A[i][j] == B[(i+a) % m][(j+b) % n]
                   for i in range(m) for j in range(n))]


def encode_relation(m, n, relation):
    """Return a coset/lattice record, rejecting noncoset input relations."""
    R = set(relation)
    if not R:
        return None
    if any(not (0 <= a < m and 0 <= b < n) for a, b in R):
        raise ValueError('relation outside its named chart domains')
    ra, rb = min(R)
    H = {((a-ra) % m, (b-rb) % n) for a, b in R}
    h = min([a for a, b in H if b == 0 and a > 0] + [m])
    v = min([b for a, b in H if b > 0] + [n])
    s = 0 if v == n else min(a for a, b in H if b == v) % h
    if m % h or n % v or (s*(n//v)) % h:
        raise ValueError('invalid lattice periods')
    out = dict(m=m, n=n, offset=[ra, rb], h=h, s=s, v=v)
    decoded = {(a, b) for a in range(m) for b in range(n)
               if contains(out, a, b)}
    if decoded != R:
        raise ValueError('relation is not one subgroup coset')
    return out


def contains(record, a, b):
    if record is None:
        return False
    ra, rb = record['offset']
    h, s, v = (record[k] for k in ('h', 's', 'v'))
    return ((b-rb) % v == 0 and
            (v*(a-ra)-s*(b-rb)) % (h*v) == 0)


def block_record(A, B):
    return encode_relation(len(A), len(A[0]), shifts(A, B))


def rows_for_edge(record, i, j, size):
    """Two periodic integer congruences; no division or free port copies."""
    if record is None:
        return [([0]*size, 1, 2)]
    ra, rb = record['offset']
    h, s, v = (record[k] for k in ('h', 's', 'v'))
    row1 = [0]*size
    row1[j] = 1
    row2 = [0]*size
    row2[i], row2[j] = v, -s
    return [(row1, rb, v), (row2, v*ra-s*rb, h*v)]


def egcd(a, b):
    aa, bb = abs(a), abs(b)
    x0, x1, y0, y1 = 1, 0, 0, 1
    while bb:
        q, rr = divmod(aa, bb)
        aa, bb = bb, rr
        x0, x1 = x1, x0-q*x1
        y0, y1 = y1, y0-q*y1
    return aa, x0*(1 if a >= 0 else -1), y0*(1 if b >= 0 else -1)


def solve(moduli, equations):
    """All solutions as offset + an integer lattice, modulo fixed moduli.

    Elementary unimodular row reduction, with column-Hermite normalization.
    The determinant gives the exact number of solutions in the residue box.
    """
    n = len(moduli)
    if not n or any(m < 1 for m in moduli):
        raise ValueError('positive moduli required')
    x = Matrix.zeros(n, 1)
    B = Matrix.eye(n)
    for coeff, rhs, modulus in equations:
        if len(coeff) != n or modulus < 1:
            raise ValueError('invalid equation')
        if any(a*m % modulus for a, m in zip(coeff, moduli)):
            raise ValueError('equation does not descend to the stated tori')
        row = Matrix([coeff])
        w = [int(t) for t in row*B]
        target = (rhs-int((row*x)[0])) % modulus
        U = Matrix.eye(n)
        for j in range(1, n):
            u, v = w[0], w[j]
            if not u and not v:
                continue
            d, alpha, beta = egcd(u, v)
            c0, cj = U[:, 0], U[:, j]
            U[:, 0] = alpha*c0 + beta*cj
            U[:, j] = (-v//d)*c0 + (u//d)*cj
            w[0], w[j] = d, 0
        d = w[0]
        g = gcd(d, modulus)
        if target % g:
            return None
        step = modulus//g
        residue = 0 if step == 1 else (target//g)*pow(d//g, -1, step) % step
        B = B*U
        x += B[:, 0]*residue
        B[:, 0] *= step
        B = hermite_normal_form(B)
        x = Matrix([int(x[i]) % moduli[i] for i in range(n)])
    volume = int(B.det())
    total = prod(moduli)
    if volume <= 0 or total % volume:
        raise ArithmeticError('invalid period-lattice index')
    return dict(moduli=list(moduli), offset=list(map(int, x)),
                basis=[[int(B[i, j]) for j in range(n)] for i in range(n)],
                count=total//volume, index=volume)


def in_solution(certificate, point):
    if certificate is None:
        return False
    B = certificate['basis']
    n = len(point)
    z = [0]*n
    for i in reversed(range(n)):
        value = point[i]-certificate['offset'][i]-sum(B[i][j]*z[j] for j in range(i+1, n))
        if value % B[i][i]:
            return False
        z[i] = value//B[i][i]
    return True


@lru_cache(None)
def cyclotomic_coeffs(order):
    z = Symbol('z')
    return tuple(reversed([int(c) for c in Poly(cyclotomic_poly(order, z), z).all_coeffs()]))


def reduce_cyclotomic(poly, order):
    """Exact Q[z]/Phi_order(z) reduction, never a numerical root test."""
    f = cyclotomic_coeffs(order)
    out = [Fraction(poly.get(i, 0)) for i in range(max(poly, default=0)+1)]
    while len(out) >= len(f):
        c = out[-1]
        shift = len(out)-len(f)
        for j, a in enumerate(f):
            out[shift+j] -= c*a
        out.pop()
    while out and out[-1] == 0:
        out.pop()
    return tuple(out)


def charge_spectrum(certificate, weights):
    """Exact constancy test for sum weights[i]*least_residue(a_i).

    Restricted characters are keyed by their values on ALL lattice generators.
    Testing only individual generator charges would omit mixed carries.
    """
    if certificate is None:
        return dict(feasible=False)
    moduli, B, r = (certificate[k] for k in ('moduli', 'basis', 'offset'))
    if len(weights) != len(moduli):
        raise ValueError('wrong observation dimension')
    ds = [gcd(m, *row) for m, row in zip(moduli, B)]
    qs = [m//d for m, d in zip(moduli, ds)]
    order = lcm(*qs)
    terms = {}
    mean = Fraction(0)
    candidate_terms = 0
    for m, row, ri, wi, d, q in zip(moduli, B, r, weights, ds, qs):
        c, beta = ri % d, ri//d
        mean += wi*(c+Fraction(d*(q-1), 2))
        if not wi:
            continue
        for k in range(1, q):
            key = tuple(Fraction(k*b, m) % 1 for b in row)
            poly = terms.setdefault(key, {})
            candidate_terms += 1
            for t in range(q):
                exponent = ((order//q)*k*(beta-t)) % order
                poly[exponent] = poly.get(exponent, 0)+Fraction(wi*d*t, q)
    reduced = [(key, reduce_cyclotomic(poly, order)) for key, poly in terms.items()]
    nonzero = [(key, value) for key, value in reduced if value]
    return dict(feasible=True, constant=not nonzero, mean=str(mean),
                coordinate_orders=qs, cyclotomic_order=order,
                candidate_terms=candidate_terms, distinct_characters=len(terms),
                nonzero_characters=len(nonzero),
                first_nonzero=None if not nonzero else
                dict(character=list(map(str, nonzero[0][0])), coefficient=list(map(str, nonzero[0][1]))))


def necklace_relation(p, q):
    """Owned rotations comparing two matching blocks, not just necklaces."""
    n = len(p)
    if len(q) != n or sorted(p) != list(range(n)) or sorted(q) != list(range(n)):
        raise ValueError('two permutations of the same labels required')
    dp = [(p[(i+1) % n]-p[i]) % n for i in range(n)]
    dq = [(q[(i+1) % n]-q[i]) % n for i in range(n)]
    return [(a, (q[a]-p[0]) % n) for a in range(n)
            if all(dq[(i+a) % n] == dp[i] for i in range(n))]
