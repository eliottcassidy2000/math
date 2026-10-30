"""Exact finite guards for PROOFS P1--P6; no numerical convergence proof.

Run from any cwd; normal and -O outputs must agree. No third-party libraries.
All checks use explicit exceptions, including under python -O.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial, gcd, isqrt
from pathlib import Path
import argparse
import hashlib
import json


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


class Field:
    def __init__(self, q):
        self.q = q
        self.p = 3 if q == 27 else q
        if q != 27:
            need(q > 2 and all(q % d for d in range(2, isqrt(q)+1)), 'prime')
        self.add = [[self.addition(a, b) for b in range(q)] for a in range(q)]
        self.mul = [[self.multiplication(a, b) for b in range(q)] for a in range(q)]
        self.neg = [next(b for b in range(q) if self.add[a][b] == 0) for a in range(q)]
        self.sub = [[self.add[a][self.neg[b]] for b in range(q)] for a in range(q)]
        squares = {self.mul[a][a] for a in range(1, q)}
        need(len(squares) == (q-1)//2, 'square count')
        self.chi = [0] + [1 if a in squares else -1 for a in range(1, q)]
        for a in range(q):
            val = self.power(a, (q-1)//2)
            need(val == (0 if a == 0 else 1 if self.chi[a] == 1 else self.neg[1]), 'Euler versus squares')
        self.arc = [[int(self.chi[self.sub[b][a]] == 1) for b in range(q)] for a in range(q)]
        self.out = [sum(1 << b for b in range(q) if self.arc[a][b]) for a in range(q)]
        self.inc = [sum(1 << b for b in range(q) if self.arc[b][a]) for a in range(q)]
        for a, b in combinations(range(q), 2):
            need(self.arc[a][b]+self.arc[b][a] == 1, 'tournament orientation')

    def addition(self, a, b):
        if self.q != 27:
            return (a+b) % self.q
        return sum(((a//(3**i) + b//(3**i)) % 3)*3**i for i in range(3))

    def multiplication(self, a, b):
        if self.q != 27:
            return a*b % self.q
        aa = [(a//3**i) % 3 for i in range(3)]
        bb = [(b//3**i) % 3 for i in range(3)]
        c = [0]*5
        for i in range(3):
            for j in range(3):
                c[i+j] += aa[i]*bb[j]
        # F3[u]/(u^3-u+1); u^3=u+2. Cubic has no F3 root.
        for i in (4, 3):
            c[i-2] += c[i]
            c[i-3] += 2*c[i]
        return sum((c[i] % 3)*3**i for i in range(3))

    def power(self, a, e):
        v = 1
        while e:
            if e & 1:
                v = self.mul[v][a]
            a = self.mul[a][a]
            e >>= 1
        return v

    def poly(self, coeffs, x):
        y = 0
        for c in reversed(coeffs):
            y = self.add[self.mul[y][x]][c % self.p]
        return y


def root_counts(f, roots):
    counts = []
    for bits in range(1 << len(roots)):
        candidates = (1 << f.q)-1
        for i, a in enumerate(roots):
            candidates &= f.inc[a] if (bits >> i) & 1 else f.out[a]
        counts.append(candidates.bit_count())
    need(sum(counts) == f.q-len(roots), 'root partition')
    return counts


def root_bound(q, t, n):
    if t == 0:
        return n == q
    a = (t-2)*(1 << (t-1))+1
    b = t*(1 << (t-1))
    excess = max(0, abs((1 << t)*n-q)-b)
    return excess*excess <= a*a*q


def bank_bound(q, r, degree, z, n):
    if r == 0:
        return n == q-z
    a = ((1 << r)-1)*(degree-1)
    b = ((1 << r)-1)*z
    excess = max(0, abs((1 << r)*n-(q-z))-b)
    return excess*excess <= a*a*q


def root_census():
    rows = []
    for q in (3, 7, 11, 19, 23, 27, 31, 43):
        f = Field(q)
        tmax = min(q, 4 if q <= 23 else 3)
        records = []
        zeros = Counter()
        for t in range(tmax+1):
            for roots in combinations(range(q), t):
                counts = root_counts(f, roots)
                need(all(root_bound(q, t, n) for n in counts), 'linear Weil envelope')
                zeros[t] += counts.count(0)
                records.append((roots, counts))
        rows.append({'q': q, 'max_roots': tmax, 'configurations': len(records),
                     'sign_queries': sum(len(c) for _, c in records),
                     'empty_patterns_by_root_count': dict(zeros), 'sha256': digest(records)})
    return rows


def spans(t):
    todo = {frozenset({0})}
    done = set()
    while todo:
        s = todo.pop()
        done.add(s)
        for v in range(1 << t):
            ns = s | {x ^ v for x in s}
            if ns not in done:
                todo.add(ns)
    return sorted((tuple(sorted(s)) for s in done), key=lambda s: (len(s), s))


def basis(span):
    b, s = [], {0}
    for x in span:
        if x not in s:
            b.append(x)
            s |= {y ^ x for y in tuple(s)}
    return b


def polynomial_census():
    rows = []
    spaces = spans(3)
    need(len(spaces) == 16, 'all F2^3 subspaces')
    for q in (7, 11, 19, 27, 31, 43, 139):
        f = Field(q)
        # Distinct irreducibles X, X-1, X^2+1 over q=3 mod4.
        values = [[x, f.sub[x][1], f.add[f.mul[x][x]][1]] for x in range(q)]
        z = {x for x, v in enumerate(values) if 0 in v}
        need(z == {0, 1}, 'polynomial deleted set')
        records = []
        for space in spaces:
            b = basis(space)
            rank = len(b)
            active = 0
            for v in b:
                active |= v
            degree = sum(d for i, d in enumerate((1, 1, 2)) if (active >> i) & 1)
            counts = [0]*(1 << rank)
            for x, vals in enumerate(values):
                if x in z:
                    continue
                bits = sum((f.chi[v] == -1) << i for i, v in enumerate(vals))
                word = sum(((bits & v).bit_count() % 2) << i for i, v in enumerate(b))
                counts[word] += 1
                # Every product query in the row space is recovered from the basis.
                for bm in range(1 << rank):
                    v = 0
                    for i, row in enumerate(b):
                        if bm >> i & 1:
                            v ^= row
                    need((bits & v).bit_count() % 2 == (word & bm).bit_count() % 2, 'redundant query')
            need(all(bank_bound(q, rank, degree, len(z), n) for n in counts), 'polynomial envelope')
            # Exact all-Boolean relative contraction: average of selected deviations.
            max_dev = max(abs((1 << rank)*n-(q-len(z))) for n in counts)
            for menu in range(1, 1 << len(counts)):
                mass = sum(n for i, n in enumerate(counts) if menu >> i & 1)
                size = menu.bit_count()
                need(abs((1 << rank)*mass-size*(q-len(z))) <= size*max_dev, 'Boolean closure')
            records.append((space, degree, counts))
        rows.append({'q': q, 'spaces': len(records), 'sha256': digest(records)})
    return rows


def pattern_key(f, roots, vertices, color=None):
    colors = tuple(color[x] for x in vertices) if color else ()
    signs = tuple(f.arc[a][x] for x in vertices for a in roots)
    edges = tuple(f.arc[a][b] for a, b in combinations(vertices, 2))
    return colors, signs, edges


def sequential_census(q, roots, k, colored):
    f = Field(q)
    z = {0, 1} if colored else set()
    if colored:
        # rank2: chi(X), chi((X-1)(X^2+1)); their signed product is redundant.
        color = {}
        for x in range(q):
            if x not in z:
                a = f.chi[x]
                b = f.chi[f.mul[f.sub[x][1]][f.add[f.mul[x][x]][1]]]
                color[x] = (int(a == -1), int(b == -1))
    else:
        color = None
    allowed = [x for x in range(q) if x not in z and x not in roots]
    direct = Counter(pattern_key(f, roots, v, color) for v in permutations(allowed, k))
    recursive = Counter()
    frontier = [()]
    for j in range(k):
        newfront = []
        for vs in frontier:
            prior = tuple(roots)+vs
            # Independent mask intersections, preserving actual earlier endpoints.
            for bits in range(1 << len(prior)):
                mask = (1 << q)-1
                for x in z:
                    mask &= ~(1 << x)
                for i, a in enumerate(prior):
                    mask &= f.inc[a] if bits >> i & 1 else f.out[a]
                while mask:
                    low = mask & -mask
                    newfront.append(vs+(low.bit_length()-1,))
                    mask ^= low
        frontier = newfront
    for vs in frontier:
        recursive[pattern_key(f, roots, vs, color)] += 1
    need(direct == recursive, 'literal versus sequential joint entries')
    rank, degree = (2, 4) if colored else (0, 0)
    # Rational outer envelope for (4), using ceil(sqrt q).
    lo = hi = Fraction(1)
    sqceil = isqrt(q)+(isqrt(q)**2 != q)
    for j in range(k):
        t = len(roots)+j
        nt = 1 << (rank+t)
        mid = Fraction(q-len(z)-t, nt)
        err = Fraction(nt-1, nt)*((degree+t-1)*sqceil+len(z)+t) if nt > 1 else 0
        lo *= max(0, mid-err)
        hi *= mid+err
    need(all(lo <= n <= hi for n in direct.values()), 'joint envelope')
    possible = 1 << (rank*k+len(roots)*k+comb(k, 2))
    if len(direct) < possible:
        need(lo == 0, 'missing word versus positive lower bound')
    records = sorted((str(key), val) for key, val in direct.items())
    return {'q': q, 'roots': roots, 'k': k, 'colored': colored,
            'nonempty_joint_patterns': len(direct), 'formal_patterns': possible,
            'tuple_mass': sum(direct.values()), 'sha256': digest(records)}


def frobenius_control():
    f = Field(27)
    phi = [f.power(a, 3) for a in range(27)]
    need(sorted(phi) == list(range(27)), 'Frobenius bijection')
    need(all(phi[phi[phi[a]]] == a for a in range(27)), 'Frobenius order divides3')
    need(any(phi[a] != a for a in range(27)), 'nontrivial phase')
    for a, b in product(range(27), repeat=2):
        need(phi[f.add[a][b]] == f.add[phi[a]][phi[b]], 'Frobenius add')
        need(phi[f.mul[a][b]] == f.mul[phi[a]][phi[b]], 'Frobenius multiply')
        need(f.arc[a][b] == f.arc[phi[a]][phi[b]], 'joint arc phase')
    roots = (0, 3, 9)
    target = tuple(phi[a] for a in roots)
    need(root_counts(f, roots) == root_counts(f, target), 'transported root answers')
    # Moving x alone does not preserve named-root answers.
    hostile = next((a, x) for a, x in product(range(27), repeat=2)
                   if a != x and a != phi[x] and f.arc[a][x] != f.arc[a][phi[x]])
    return {'roots': roots, 'transported_roots': target, 'fixed_root_phase_hostile': hostile}


def apollonius_census():
    rows = []
    for q in (7, 11, 19, 31, 43, 139, 3511):
        counts = Counter()
        z = 0
        for t in range(q):
            a, b, c = (t*t-1) % q, (t*t-9) % q, (t*t-24) % q
            if t*a*b*c % q == 0:
                z += 1
                continue
            reduced = [-3*a, -15*a*c, -3*b*c, -15*b]
            defects = [-3*a*t*t, -15*a*c, -3*b*c, -15*b*t*t]
            chars = [pow(x % q, (q-1)//2, q) for x in reduced]
            need(chars == [pow(x % q, (q-1)//2, q) for x in defects], 'square-denominator removal')
            word = sum((v == q-1) << i for i, v in enumerate(chars))
            need(word.bit_count() % 2 == 0, 'Apollonius parity')
            counts[word] += 1
        words = [w for w in range(16) if w.bit_count() % 2 == 0]
        need(z <= 7 and sum(counts.values()) == q-z, 'Apollonius source mass')
        need(all(bank_bound(q, 3, 6, z, counts[w]) for w in words), 'Apollonius bound')
        rows.append({'q': q, 'excluded': z, 'even_word_counts': [counts[w] for w in words]})
    # At 23 the C polynomial equals A, so the polynomial rank drops to2.
    need(24 % 23 == 1, 'bad characteristic witness')
    return rows


def prime_bank_census():
    configs = [([2, 3, 6], 24), ([-1, 2, -2], 8),
               ([4, 2, 5, 3, 7, 13], 10920), ([2, 3, 5, 6, 10, 15, 30], 120)]
    records = []
    def jacobi(a, n):
        a %= n
        s = 1
        while a:
            while a % 2 == 0:
                a //= 2
                if n % 8 in (3, 5):
                    s = -s
            a, n = n, a
            if a % 4 == n % 4 == 3:
                s = -s
            a %= n
        return s if n == 1 else 0
    for bank, modulus in configs:
        cnt = Counter(tuple(jacobi(a, n) for a in bank)
                      for n in range(3, modulus, 4) if gcd(n, modulus) == 1)
        need(len(set(cnt.values())) == 1, 'exact CRT equal fibers')
        records.append({'bank': bank, 'Q': modulus, 'colors': len(cnt),
                        'classes_per_color': next(iter(cnt.values()))})
    # A square adds no bit but can remove an existing legal source.
    need(jacobi(3**2, 3) == 0 and jacobi(1, 3) == 1, 'square zero-prime guard')
    return records


def code_and_guard_controls():
    f = Field(139)
    roots = (0, 1, 2)
    # All seven nonzero columns of F2^3: simplex answers/Hamming relations.
    counts = Counter()
    for x in range(f.q):
        if x in roots:
            continue
        bits = sum((f.chi[f.sub[x][a]] == -1) << i for i, a in enumerate(roots))
        word = sum(((bits & col).bit_count() % 2) << (col-1) for col in range(1, 8))
        counts[word] += 1
    need(len(counts) == 8 and {w.bit_count() for w in counts} == {0, 4}, 'simplex realized')
    relation_count = 0
    for mask in range(128):
        xor = 0
        for col in range(1, 8):
            if mask >> (col-1) & 1:
                xor ^= col
        if xor == 0:
            relation_count += 1
            need(all((w & mask).bit_count() % 2 == 0 for w in counts), 'Hamming dual')
    need(relation_count == 16, 'Hamming dimension4')
    # Even all square multiples can delete every source without adding a bit.
    need(all(any(f.mul[f.sub[x][a]][f.sub[x][a]] == 0 for a in range(f.q))
             for x in range(f.q)), 'all-square bank has empty domain')
    # X and an alleged new root at0 are the SAME query, not independent.
    need(all(f.arc[0][x] == (f.chi[x] == 1) for x in range(1, f.q)), 'root support collision')
    return {'q': 139, 'simplex_word_counts': sorted(counts.items()),
            'Hamming_relation_count': relation_count,
            'support_collision_rank': 1, 'square_bank_rank': 0, 'square_bank_domain_size': 0}


def word_census():
    def step(n, sigma):
        v, d = 3*n+sigma, 0
        while v % 2 == 0:
            d += 1
            v //= 2
        return v, d
    records, rejected = [], Counter()
    for sigma in (1, -1):
        for length in range(1, 5):
            for word in product((1, 2, 3), repeat=length):
                total = sum(word)
                modulus = 1 << (total+1)
                valid = []
                for start in range(1, modulus, 2):
                    n, observed, values = start, [], [start]
                    for d in word:
                        n, actual = step(n, sigma)
                        observed.append(actual)
                        values.append(n)
                        if d != actual:
                            break
                    if tuple(observed) == word:
                        valid.append(values)
                need(len(valid) == 1, 'unique exact word cylinder')
                b = valid[0]
                prefix = 0
                a = [modulus]
                for j, d in enumerate(word, 1):
                    prefix += d
                    a.append(3**j * 2**(total+1-prefix))
                dets = [a[i]*b[j]-a[j]*b[i] for i, j in combinations(range(length+1), 2)]
                need(all(sigma*d > 0 for d in dets), 'nonzero same-sign determinants')
                charts = []
                for p in (7, 11, 19, 43, 139):
                    if any(d % p == 0 for d in dets):
                        rejected[p] += 1
                        continue
                    counts = [0]*(1 << (length+1))
                    excluded = 0
                    for t in range(p):
                        n = b[0]+a[0]*t
                        values = [n]
                        for d in word:
                            n, actual = step(n, sigma)
                            need(actual == d, 'retained source valuation')
                            values.append(n)
                        need(values == [aa*t+bb for aa, bb in zip(a, b)], 'affine intermediate values')
                        if any(n % p == 0 for n in values):
                            excluded += 1
                            continue
                        bits = sum((pow(n % p, (p-1)//2, p) == p-1) << i
                                   for i, n in enumerate(values))
                        counts[bits] += 1
                    need(excluded == length+1, 'distinct linear roots')
                    need(all(root_bound(p, length+1, n) for n in counts), 'word sign envelope')
                    charts.append((p, counts))
                records.append((sigma, word, a, b, dets, charts))
    return {'words': len(records), 'legal_prime_charts': sum(len(r[-1]) for r in records),
            'rejected_prime_charts': dict(rejected), 'sha256': digest(records),
            'scope': 'length1..4, each valuation1..3, both fixed signs; primes7,11,19,43,139'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'FINITE-EXACT; asymptotic laws proved separately',
              'root_census': root_census(), 'polynomial_census': polynomial_census(),
              'sequential': [sequential_census(7, (), 4, False),
                             sequential_census(11, (0, 1), 3, False),
                             sequential_census(11, (2,), 3, True),
                             sequential_census(27, (3, 9), 2, True)],
              'frobenius': frobenius_control(), 'apollonius': apollonius_census(),
              'prime_CRT': prime_bank_census(),
              'code_and_guard_controls': code_and_guard_controls(),
              'word_census': word_census(),
              'hostile_81_27_orbit_ratio': str(Fraction(factorial(12)//27, factorial(12)//81))}
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text, encoding='utf-8')
    print(text)


if __name__ == '__main__':
    main()
