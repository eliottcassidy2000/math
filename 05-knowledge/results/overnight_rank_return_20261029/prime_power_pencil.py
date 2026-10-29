"""Prime-power Paley pencils: carry recurrence, Teichmuller sums, native matrices."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path

from paley_pencil import diagonal, local_layers, predicted as prime_predicted, valuation

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise ValueError(why)


def carry_trace(p, m):
    """tr W(z)^m, using the classical quadratic transfer denominator."""
    a = (p+1)//2
    previous, current = [2], [a, a]
    if m == 0:
        return previous
    for degree in range(2, m+1):
        out = [0]*(degree+1)
        for j, value in enumerate(current):
            out[j] += a*value
            out[j+1] += a*value
        for j, value in enumerate(previous):
            out[j+1] -= p*value
        previous, current = current, out
    return current


def digit_sum(n, p):
    out = 0
    while n:
        out += n % p
        n //= p
    return out


def carry_sources(p, m):
    q, r = p**m, (p**m-1)//2
    return {i: (digit_sum(i, p)+digit_sum(r, p)-digit_sum((i+r) % (q-1), p))//(p-1)
            for i in range(1, q-1) if i != r}


def primitive_layers(p, m, s):
    q = p**m
    f = carry_trace(p, m)
    if s == 0:
        f[-1] -= 2
        return Counter({e: count for e, count in enumerate(f) if count}), 1
    t = valuation(s, p)
    if t == 0:
        return Counter({0: q}), 0
    h = min(2*t, m)
    out = Counter({0: 2})
    out[t+h] += 1
    for c in range((m+1)//2):
        count = f[c]-2*int(c == 0)
        e = min(t, c)
        out[e] += count
        out[h-e] += count
    return +out, 0


def full_diagonal(p, m, a, b):
    q, r, g = p**m, (p**m-1)//2, gcd(a, b)
    if g == 0:
        return [0]*q
    x, y = a//g, b//g
    s, delta = x+y, x-y
    counts, free = primitive_layers(p, m, s)
    exponents = [e for e in sorted(counts) for _ in range(counts[e])]
    if free:
        return [g*p**e for e in exponents]+[0]
    norm = (s*s+q*delta*delta)//4
    scalar, d = abs(r*s), gcd(abs(r*s), norm)
    away = [1]*r+[d]+[norm]*(r-1)+[scalar*norm//d]
    away = [value//p**valuation(value, p) for value in away]
    need(len(exponents) == q, 'all invariant factors counted')
    result = [g*value*p**e for value, e in zip(away, exponents)]
    need(all(y % x == 0 for x, y in zip(result, result[1:])), 'Smith divisibility chain')
    return result


class Field:
    def __init__(self, p, coefficients):
        self.p, self.coefficients = p, coefficients
        self.m, self.q = len(coefficients)-1, p**(len(coefficients)-1)
        self.digits = [tuple((x//p**i) % p for i in range(self.m)) for x in range(self.q)]

    def encode(self, digits):
        return sum((a % self.p)*self.p**i for i, a in enumerate(digits))

    def subtract(self, x, y):
        return self.encode([a-b for a, b in zip(self.digits[x], self.digits[y])])

    def ring_multiply(self, x, y, modulus):
        out = [0]*(2*self.m-1)
        for i, a in enumerate(x):
            for j, b in enumerate(y):
                out[i+j] += a*b
        for degree in range(2*self.m-2, self.m-1, -1):
            coefficient = out[degree]
            for j in range(self.m):
                out[degree-self.m+j] -= coefficient*self.coefficients[j]
        return tuple(a % modulus for a in out[:self.m])

    def ring_power(self, x, exponent, modulus):
        out = (1,)+(0,)*(self.m-1)
        while exponent:
            if exponent & 1:
                out = self.ring_multiply(out, x, modulus)
            x = self.ring_multiply(x, x, modulus)
            exponent >>= 1
        return out

    def power(self, x, exponent):
        return self.encode(self.ring_power(self.digits[x], exponent, self.p))

    def matrix(self):
        need(all(self.power(x, self.q-1) == 1 for x in range(1, self.q)), 'all nonzero elements are units')
        squares = {self.power(x, 2) for x in range(1, self.q)}
        need(len(squares) == (self.q-1)//2, 'quadratic residue count')
        A = [[int(self.subtract(j, i) in squares) for j in range(self.q)] for i in range(self.q)]
        need(all(A[i][j]+A[j][i] == int(i != j) for i in range(self.q) for j in range(self.q)),
             'Paley tournament relation')
        frobenius = [self.power(x, self.p) for x in range(self.q)]
        need(len(set(frobenius)) == self.q, 'Frobenius permutation')
        need(all(A[frobenius[i]][frobenius[j]] == A[i][j]
                 for i in range(self.q) for j in range(self.q)), 'native Frobenius transport')
        return A

    def jacobi_check(self):
        p, m, q = self.p, self.m, self.q
        divisors = [d for d in range(2, q) if (q-1) % d == 0 and all(d % e for e in range(2, isqrt(d)+1))]
        generator = next(g for g in range(1, q) if all(self.power(g, (q-1)//ell) != 1 for ell in divisors))
        field_powers = [self.power(generator, k) for k in range(q-1)]
        need(len(set(field_powers)) == q-1, 'primitive field generator')
        logarithm = {x: k for k, x in enumerate(field_powers)}
        modulus = p**(m+1)
        lift = self.digits[generator]
        for _ in range(m+1):
            lift = self.ring_power(lift, q, modulus)
        one = (1,)+(0,)*(m-1)
        need(self.ring_power(lift, q-1, modulus) == one, 'Teichmuller root of unity')
        powers = [one]
        for _ in range(1, q-1):
            powers.append(self.ring_multiply(powers[-1], lift, modulus))
        sources, sums = carry_sources(p, m), {}
        for i, count in sources.items():
            acc = [0]*m
            for x in range(1, q):
                y = self.subtract(1, x)
                if y == 0:
                    continue
                sign = 1 if self.power(y, (q-1)//2) == 1 else -1
                value = powers[(-i*logarithm[x]) % (q-1)]
                for j in range(m):
                    acc[j] += sign*value[j]
            value = tuple(x % modulus for x in acc)
            v = min((valuation(x, p) for x in value if x), default=m+1)
            need(v == count, 'direct Jacobi valuation equals cyclic carry count')
            sums[i] = value
        r = (q-1)//2
        for i in sources:
            need(self.ring_multiply(sums[i], sums[(i+r) % (q-1)], modulus)
                 == ((-q) % modulus,)+(0,)*(m-1), 'opposite block product is minus q')
            need(sums[(p*i) % (q-1)] == sums[i], 'Frobenius permutes sources with equal Jacobi sum')
        unseen, orbits = set(sources), []
        while unseen:
            i, orbit = min(unseen), []
            while i not in orbit:
                orbit.append(i)
                unseen.remove(i)
                i = p*i % (q-1)
            need(i == orbit[0], 'closed source phase')
            orbits.append(dict(indices=orbit, valuation=sources[orbit[0]]))
        return dict(generator=generator, precision=m+1,
                    jacobi_values={str(i): value for i, value in sums.items()}, frobenius_orbits=orbits)


def field_check(field, large=False):
    p, m, q = field.p, field.m, field.q
    A = field.matrix()
    weights = {(a, b) for a in range(-4, 5) for b in range(-4, 5)} if not large else set()
    weights.update([(1, 0), (1, -1), (1, p-1), (1, p*p-1), (1, p**3-1), (p, p*(p-1))])
    if not large:
        for t in range(1, 6):
            weights.update([(1, p**t-1), (p, p*(p**t-1))])
    records = []
    for a, b in sorted(weights):
        E = [[a*A[i][j]+b*A[j][i] for j in range(q)] for i in range(q)]
        expected = full_diagonal(p, m, a, b)
        nonzero = [x for x in expected if x]
        max_exponent = max((valuation(x, p) for x in nonzero), default=0)
        observed = local_layers(E, p, max_exponent+1)
        counts = Counter(valuation(x, p) for x in nonzero)
        need(observed == [counts[e] for e in range(max_exponent+1)], f'native p-primary layers q={q}, a={a}, b={b}')
        if not large:
            need(diagonal(E) == expected, 'independent full integer Smith computation')
        records.append(dict(a=a, b=b, p_layers=observed, free_rank=q-len(nonzero)))
    return dict(p=p, degree=m, q=q, modulus_polynomial=field.coefficients,
                matrix_sha256=sha256(json.dumps(A, separators=(',', ':')).encode()).hexdigest(),
                actual_matrix_checks=records, jacobi=field.jacobi_check())


def run():
    recurrence_checks = []
    for p in (3, 7, 11):
        for m in range(1, 102):
            f = carry_trace(p, m)
            need(f == f[::-1] and min(f) >= 0 and sum(f) == p**m+1, 'trace polynomial identities')
            if p**m <= 20000:
                direct = Counter(carry_sources(p, m).values())
                adjusted = f[:]
                adjusted[0] -= 2
                adjusted[-1] -= 2
                need(adjusted == [direct[e] for e in range(m+1)], 'direct digit-sum census')
        for m in (1, 3, 5, 7, 9, 31, 101):
            counts, free = primitive_layers(p, m, 0)
            need(sum(counts.values())+free == p**m, 'singular size')
            need(sum(e*v for e, v in counts.items()) == m*(p**m-3)//2, 'skew determinant valuation')
            recurrence_checks.append(dict(p=p, degree=m, skew_layers=dict(sorted(counts.items()))))
    old = json.loads((HERE/'paley_pencil.json').read_text(encoding='utf-8'))
    for case in old['cases']:
        q, a, b = case['q'], case['a'], case['b']
        need(full_diagonal(q, 1, a, b) == case['smith'] == prime_predicted(q, a, b),
             'all556 prime cases are exact specializations')
    fields = [field_check(Field(3, [-1, -1, 0, 1])),
              field_check(Field(3, [-1, -1, 0, 0, 0, 1]), True),
              field_check(Field(7, [-2, 0, 0, 1]), True)]
    return dict(recurrence_checks=recurrence_checks, prime_specializations=len(old['cases']), fields=fields)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    encoded = json.dumps(result, indent=2)+'\n'
    target = HERE/'prime_power_pencil.json'
    if args.write:
        target.write_text(encoded, encoding='utf-8')
    else:
        need(json.loads(encoded) == json.loads(target.read_text(encoding='utf-8')), 'receipt replay')
    print(json.dumps(dict(status='PASS', prime_specializations=result['prime_specializations'],
                         fields=[dict(q=f['q'], native_pencils=len(f['actual_matrix_checks']),
                                      jacobi_sums=len(f['jacobi']['jacobi_values'])) for f in result['fields']])) )


if __name__ == '__main__':
    main()
