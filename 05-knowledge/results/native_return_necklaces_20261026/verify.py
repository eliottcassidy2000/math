#!/usr/bin/env python3
"""Exact controls; no external packages and no removable asserts.

Run from any directory: python3 verify.py [--output receipt.json].
The finite universes below test the proofs in README.md, not universal Collatz.
"""
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
import argparse
import hashlib
import json


def require(condition, context):
    if not condition:
        raise RuntimeError(context)


@lru_cache(None)
def pairs(n):
    return tuple((i, j) for i in range(n) for j in range(i + 2, n))


@lru_cache(None)
def positions(n):
    return {pair: z for z, pair in enumerate(pairs(n))}


def edge(n, t, i, j):
    if i > j:
        return 1 - edge(n, t, j, i)
    return 1 if j == i + 1 else (t >> positions(n)[i, j]) & 1


def tau(n, t):
    return t ^ ((1 << len(pairs(n))) - 1)


def combine(a, x, b, y, cross):
    return sum((edge(a, x, i, j) if j < a else
                edge(b, y, i-a, j-a) if i >= a else cross) << z
               for z, (i, j) in enumerate(pairs(a+b)))


def reframe(n, t, p):
    require(all(edge(n, t, p[i], p[i+1]) for i in range(n-1)), 'illegal path')
    return sum(edge(n, t, p[i], p[j]) << z for z, (i, j) in enumerate(pairs(n)))


def rho(p):
    return sum(int(p[i] < p[j]) << z for z, (i, j) in enumerate(pairs(len(p))))


def separable(p):
    if len(p) <= 1:
        return True
    for k in range(1, len(p)):
        if max(p[:k]) < min(p[k:]) or min(p[:k]) > max(p[k:]):
            if separable(p[:k]) and separable(p[k:]):
                return True
    return False


def grammar_checks():
    g = {1: {0}}
    counts = {1: 1}
    for n in range(2, 11):
        g[n] = {combine(a, x, n-a, y, c) for a in range(1, n)
                for x in g[a] for y in g[n-a] for c in (0, 1)}
        counts[n] = len(g[n])
        require({tau(n, t) for t in g[n]} == g[n], ('tau closure', n))
        if n <= 7:
            independent = {rho(p) for p in permutations(range(n)) if separable(p)}
            require(independent == g[n], ('separable image', n))
    require(set(range(8)) - g[4] == {2, 5}, 'order-four exclusions')
    require(rho((2, 4, 1, 3)) == 2 and rho((3, 1, 4, 2)) == 5, 'forbidden witnesses')
    require(reframe(4, 1, (0, 2, 3, 1)) == 2, 'relocation escape')
    require(1 in g[4] and 2 not in g[4], 'escape membership')
    return {'counts': counts, 'permutation_cross_check_through_n': 7,
            'relocation_escape': {'input': 1, 'path': [0, 2, 3, 1], 'output': 2}}


MAPS = {'I': (1, 0), 'F': (1, 1), 'C0': (0, 0), 'C1': (0, 1)}


def compiler_checks():
    words = states = 0
    for length in range(8):
        for w in product(MAPS, repeat=length):
            words += 1
            for source in (0, 1):
                n, t, b, flips, clock = 3, source, source, 0, 0
                for letter in w:
                    old_n, old_t = n, t
                    a, c = MAPS[letter]
                    b = a*b ^ c
                    if letter == 'F':
                        t = tau(n, t)
                        flips ^= 1
                        clock += 1
                    elif letter in ('C0', 'C1'):
                        cross = int(letter == 'C1')
                        t = combine(n, t, 1, 0, cross)
                        if cross == 0:
                            primitive = tau(n+1, combine(n, tau(n, old_t), 1, 0, 1))
                            require(t == primitive, 'C0 native expansion')
                        n += 1
                        clock += 1 if cross else 3
                        for i, j in pairs(old_n):
                            require(edge(n, t, i, j) == edge(old_n, old_t, i, j), 'archive mutation')
                    require(edge(n, t, 0, n-1) == b, ('output', w, source))
                    require(edge(n, t, 0, 2) ^ flips == source, 'source recovery')
                    for j in range(2, n):
                        require(len({edge(n, t, i, j) for i in range(j-1)}) == 1, 'column carrier')
                    states += 1
                require(n == 3+w.count('C0')+w.count('C1'), 'vertex budget')
                require(clock == w.count('F')+w.count('C1')+3*w.count('C0'), 'primitive clock')
    return {'words_including_empty': words, 'sources_per_word': 2, 'prefix_states_checked': states}


def difference_checks():
    trials = 0
    for n in range(3, 6):
        for a in range(1 << len(pairs(n))):
            for b in range(a+1, 1 << len(pairs(n))):
                require(tau(n, a)^tau(n, b) == a^b, 'tau difference')
                for i in range(n):
                    for j in range(i+2, n):
                        p = list(range(n))
                        p.insert(j, p.pop(i))
                        if all(edge(n, t, p[k], p[k+1]) for t in (a, b) for k in range(n-1)):
                            aa, bb = reframe(n, a, p), reframe(n, b, p)
                            for u in range(n):
                                for v in range(u+1, n):
                                    require(edge(n, aa, u, v)^edge(n, bb, u, v) ==
                                            edge(n, a, p[u], p[v])^edge(n, b, p[u], p[v]), 'labelled difference')
                            trials += 1
    return {'all_pairs_through_n': 5, 'common_legal_relocations': trials}


def pmul(a, b, f):
    result, degree = 0, f.bit_length()-1
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        if a >> degree:
            a ^= f
        b >>= 1
    return result


def rank(columns):
    basis = {}
    for a in columns:
        while a:
            k = a.bit_length()-1
            if k not in basis:
                basis[k] = a
                break
            a ^= basis[k]
    return len(basis)


def polyrem(a, b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a


def irreducible(f):
    m = f.bit_length()-1
    return all(polyrem(f, d) != 0 for k in range(1, m//2+1)
               for d in range((1 << k)|1, 1 << (k+1), 2))


def field_checks():
    results = []
    for m in range(2, 10):
        f = next(f for f in range((1 << m)|1, 1 << (m+1), 2) if irreducible(f))
        for beta in range(1, 1 << m):
            basis = [beta]
            for _ in range(m-1):
                basis.append(pmul(basis[-1], basis[-1], f))
            if rank(basis) == m:
                break
        require(rank(basis) == m, 'normal basis missing')
        def encode(w):
            a = 0
            for i, b in enumerate(basis):
                if w >> i & 1:
                    a ^= b
            return a
        require(encode((1 << m)-1) == 1, 'trace normal generator')
        for w in range(1 << m):
            a = encode(w)
            rotation = ((w << 1) & ((1 << m)-1)) | (w >> (m-1))
            require(pmul(a, a, f) == encode(rotation), 'Frobenius rotation')
            trace, t = 0, a
            for _ in range(m):
                trace ^= t
                t = pmul(t, t, f)
            require(trace == (w.bit_count() & 1), 'trace parity')
            t, d = pmul(a, a, f), 1
            while t != a:
                t = pmul(t, t, f)
                d += 1
            word_period = next(k for k in range(1, m+1)
                               if all((w >> i & 1) == (w >> ((i+k) % m) & 1) for i in range(m)))
            require(d == word_period, 'degree-period')
        # A rank-one coordinate selector is not multiplication by beta.
        require(rank([pmul(beta, x, f) for x in basis]) == m, 'field multiplication rank')
        results.append({'m': m, 'irreducible_polynomial_bits': f, 'normal_generator': beta,
                        'words': 1 << m, 'selector_rank': 1, 'multiplication_rank': m})
    return results


def cyclic_checks():
    words = 0
    for n in range(1, 7):
        for w in product(MAPS.values(), repeat=n):
            a, c = 1, 0
            for m, b in w:
                a, c = m*a, m*c ^ b
            count = sum(all(((s >> ((i+1) % n)) & 1) == m*((s >> i) & 1)^b
                            for i, (m, b) in enumerate(w)) for s in range(1 << n))
            predicted = 1 if a == 0 else 2 if c == 0 else 0
            require(count == predicted, ('cyclic sections', w))
            columns = []
            for j in range(n):
                col = 0
                for i, (m, _) in enumerate(w):
                    col |= (int((i+1) % n == j) ^ (m*int(i == j))) << i
                columns.append(col)
            require(rank(columns) == n-int(a == 1), 'cyclic operator rank')
            words += 1
    def comp(w):
        return tuple(__import__('functools').reduce(lambda b, l: MAPS[l][0]*b ^ MAPS[l][1], w, s) for s in (0, 1))
    require(comp(('C0', 'F')) == (1, 1) and comp(('F', 'C0')) == (0, 0), 'cut hostile')
    return {'all_four_letter_words_length_1_to_6': words, 'cut_hostile': ['C0 F = C1', 'F C0 = C0']}


def reset_inverse_checks():
    cases = 0
    for n in range(1, 10):
        for mask in range((1 << n)-1):  # excludes the no-reset/all-one mask
            ms = [(mask >> i) & 1 for i in range(n)]
            zeros = [i for i, v in enumerate(ms) if not v]
            lengths = [((zeros[(j+1) % len(zeros)]-z-1) % n)+1 for j, z in enumerate(zeros)]
            # Independent real impulse recurrence around a reset, not matrix inversion.
            inverse = [[0]*n for _ in range(n)]
            cut = zeros[0]
            for source in range(n):
                value = 0
                for step in range(n):
                    i = (cut+step) % n
                    value = ms[i]*value+int(i == source)
                    inverse[(i+1) % n][source] = value
            # (S-D) inverse = I in Z, including wraparound and n=1.
            for i in range(n):
                for j in range(n):
                    require(inverse[(i+1) % n][j]-ms[i]*inverse[i][j] == int(i == j), 'integer inverse')
            row_sums = sorted(map(sum, inverse))
            require(row_sums == sorted(v for size in lengths for v in range(1, size+1)), 'reset block lengths')
            require(max(row_sums) == max(lengths), 'uniform memory norm')
            require(sum(x*x for row in inverse for x in row) == sum(s*(s+1)//2 for s in lengths), 'Frobenius norm')
            cases += 1
    # Direct determinant expansions independently check the two small spectral polynomials.
    from itertools import permutations as perms
    def determinant(a):
        n = len(a)
        return sum((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n)) *
                   __import__('math').prod(a[i][p[i]] for i in range(n)) for p in perms(range(n)))
    for ell in (2, 3):
        gram = [[ell-max(i, j) for j in range(ell)] for i in range(ell)]
        for x in range(-3, 9):
            value = determinant([[x*int(i == j)-gram[i][j] for j in range(ell)] for i in range(ell)])
            expected = x*x-3*x+1 if ell == 2 else x**3-6*x*x+5*x-1
            require(value == expected, 'golden/seventh-root polynomial')
    # rho=1/2; multiply every entry by a common power of two for integer checks.
    clock_tests = 0
    for size in range(1, 8):
        for increments in product(range(1, 4), repeat=size-1):
            times = [0]
            for dt in increments:
                times.append(times[-1]+dt)
            scale = 1 << times[-1]
            kernel = [[1 << (times[-1]-times[j]+times[i]) if i <= j else 0
                       for i in range(size)] for j in range(size)]
            require(max(map(sum, kernel)) < 2*scale, 'discount row bound')
            require(max(sum(row[i] for row in kernel) for i in range(size)) < 2*scale, 'discount column bound')
            clock_tests += 1
    return {'all_reset_masks_length_1_to_9': cases, 'norm_infinity': 'largest cyclic reset spacing',
            'Gram_polynomials': {'2': [1, -3, 1], '3': [-1, 5, -6, 1]},
            'integer_clock_discount_tests': clock_tests, 'discount': '1/2',
            'coefficient_order': 'ascending'}


def collatz_checks():
    words = points = 0
    for m in range(1, 5):
        for ks in product(range(1, 6), repeat=m):
            s = b = 0
            for k in ks:
                b = 3*b+(1 << s)
                s += k
            modulus = 1 << (s+1)
            c = (((1 << s)-b)*pow(3**m, -1, modulus)) % modulus
            require(c & 1, 'source parity')
            d = (3**m*c+b) >> s
            require(d & 1, 'target parity')
            for q in range(16):
                source = c+modulus*q
                x, clock = source, 0
                for k in ks:
                    x = 3*x+1
                    clock += 1
                    v = 0
                    while x % 2 == 0:
                        x //= 2
                        clock += 1
                        v += 1
                    require(v == k, ('valuation', ks, q))
                require(x == d+2*3**m*q and clock == s+m, 'endpoint or clock')
                require((x-1)//2 == (d-1)//2+3**m*q, 'normalized address')
                points += 1
            for r in range(1, 7):
                require(len({((d-1)//2+3**m*q) % (1 << r) for q in range(1 << r)}) == 1 << r,
                        'address bijection')
            words += 1
    return {'valuation_words_length_1_to_4_exponents_1_to_5': words, 'integer_sources': points,
            'address_bijections_moduli': [2**r for r in range(1, 7)]}


def rule30_checks():
    row = [1]+[0]*16
    machines = [(3, b) for b in row]
    for time in range(64):
        nxt = []
        for r in range(17):
            a, b = row[r-2] if r >= 2 else 0, row[r-1] if r >= 1 else 0
            mult, off = 1^b, a^b
            out = (mult*row[r]) ^ off
            n, t = machines[r]
            if mult:
                if off:
                    t = tau(n, t)
            else:
                t = combine(n, t, 1, 0, off)
                n += 1
            require(edge(n, t, 0, n-1) == out, ('Rule30 native lift', time, r))
            machines[r] = n, t
            nxt.append(out)
        row = nxt
    return {'fixed_left_front_depths': 17, 'integer_updates': 64, 'native_lifts': 1088}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    receipt = {'status': 'FINITE-EXACT PASS', 'grammar': grammar_checks(),
               'compiler': compiler_checks(), 'difference': difference_checks(),
               'normal_basis': field_checks(), 'cyclic_sections': cyclic_checks(),
               'reset_inverse': reset_inverse_checks(),
               'collatz': collatz_checks(), 'rule30': rule30_checks()}
    receipt['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = json.dumps(receipt, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(result)
    print(result, end='')


if __name__ == '__main__':
    main()
