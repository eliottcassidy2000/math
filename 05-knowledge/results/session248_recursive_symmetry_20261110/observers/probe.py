#!/usr/bin/env python3
"""Small exact probes of Session248 observer/code formulas.
Inherited physical response tables are inputs, not re-audited here.
"""
import itertools
import json
from collections import Counter, defaultdict
from fractions import Fraction


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def parity(x):
    return x.bit_count() & 1


def syndrome(b, m):
    out = 0
    for z in range(1 << m):
        if (b >> z) & 1:
            out ^= 1 | (z << 1)
    return out


def affine_mask(a, m):
    return sum(parity(a & (1 | (z << 1))) << z
               for z in range(1 << m))


def sector_means(values, labels):
    sums, sizes = defaultdict(int), Counter(labels)
    for val, label in zip(values, labels):
        sums[label] += val
    return {x: Fraction(sums[x], sizes[x]) for x in sizes}


def fraction_json(x):
    return x.numerator if x.denominator == 1 else str(x)


def main():
    records = []
    for m in (2, 3, 4):
        k = 1 << m
        labels = [syndrome(b, m) for b in range(1 << k)]
        sizes = Counter(labels)
        require(len(sizes) == 2 * k, 'surjectivity')
        require(set(sizes.values()) == {1 << (k-m-1)}, 'equal coset sizes')
        dual = [affine_mask(a, m) for a in range(1 << (m+1))]
        weights = Counter(w.bit_count() for w in dual)
        require(weights == {0:1, k//2:2*k-2, k:1}, 'affine dual spectrum')
        for s in range(1, 2*k):
            require(any(parity(a & s) and w.bit_count() == k//2
                        for a, w in enumerate(dual)), 'pair separator')
        records.append({'m':m, 'modules':k, 'sectors':len(sizes),
                        'sector_size':next(iter(sizes.values())),
                        'dual_weights':dict(sorted(weights.items())),
                        'sharp_root_threshold':k//4})

    # Verify all low-weight Walsh characters have constant coset means.
    m, k = 3, 8
    labels = [syndrome(b, m) for b in range(1 << k)]
    character_checks = 0
    for degree in range(4):
        for inds in itertools.combinations(range(k), degree):
            w = sum(1 << i for i in inds)
            vals = [(-1)**parity(b & w) for b in range(1 << k)]
            means = sector_means(vals, labels)
            require(set(means.values()) == ({1} if w == 0 else {0}),
                    'low-degree coset means')
            character_checks += 1

    # Sharp physical halfspace template, factored from credited rooted gadgets.
    f = [4161, 4297, 4933, 5085]  # bits (first,second): 00,01,10,11
    sharp = []
    for a in range(1, 2*k):
        w = affine_mask(a, m)
        if w.bit_count() != k//2:
            continue
        selected = [i for i in range(k) if (w >> i) & 1]
        vals = []
        for b in range(1 << k):
            h = 15**(k-len(selected))
            for i,j in zip(selected[::2], selected[1::2]):
                h *= f[2*((b >> i)&1) + ((b >> j)&1)]
            vals.append(h)
        means = sector_means(vals, labels)
        amplitude = 4**(k//4)*15**(k//2)
        baseline = sum(means.values(), Fraction())/(2*k)
        for s, mu in means.items():
            require(mu == baseline + amplitude*((-1)**parity(a&s)),
                    'sharp halfspace response')
        sharp.append({'dual_mask':w, 'mean_difference':2*amplitude})

    # Full-parity sharpness at every small module count, without graph census.
    parity_results = []
    for k in range(1, 9):
        vals, labels = [], []
        for b in range(1 << k):
            h = 1
            start = 0
            if k % 2:
                h *= [41,43][b & 1]
                start = 1
            for i in range(start, k, 2):
                h *= f[2*((b >> i)&1) + ((b >> (i+1))&1)]
            vals.append(h)
            labels.append(parity(b))
        means = sector_means(vals, labels)
        q = (k+1)//2
        predicted = (2*4**q if k%2 == 0 else -2*4**(q-1))
        require(means[0]-means[1] == predicted, 'parity threshold formula')
        parity_results.append({'modules':k, 'roots':q,
                               'even_minus_odd':fraction_json(means[0]-means[1])})

    # The address recursion in the proof, on the next-level m=3 bank.
    for b in range(1 << 8):
        b0, b1 = b & 15, b >> 4
        s0, s1 = syndrome(b0, 2), syndrome(b1, 2)
        expected = (s0 ^ s1) | ((s1 & 1) << 3)
        require(syndrome(b, 3) == expected, 'address recursion')

    print(json.dumps({'status':'PASS', 'scope':'new code transfer and response identities only',
                      'inherited_tables':{'pair':f, 'single':[41,43]},
                      'affine_banks':records, 'low_character_checks':character_checks,
                      'sharp_halfspace_templates':sharp,
                      'parity_thresholds':parity_results,
                      'recursive_address_cases':256}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
