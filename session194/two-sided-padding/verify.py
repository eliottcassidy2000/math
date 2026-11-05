#!/usr/bin/env python3
"""Exact standalone checks for source/sink response identities.

Subset path DP is adapted from Session 193 cut-digits/verify.py. No external
files, nonstandard packages, graph censuses or search completeness are needed.
"""
from itertools import permutations, product
from pathlib import Path
import hashlib
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length()-1
        mask -= low


def validate(rows):
    n = len(rows)
    full = (1 << n)-1
    for i, row in enumerate(rows):
        need(row >= 0 and row & ~full == 0 and not (row >> i & 1),
             'invalid outgoing row')
        for j in range(i):
            need((row >> j & 1)+(rows[j] >> i & 1) == 1,
                 'not a tournament')


def reverse(rows):
    full = (1 << len(rows))-1
    return [full ^ row ^ (1 << i) for i, row in enumerate(rows)]


def tables(rows):
    n = len(rows)
    count = 1 << n
    endings = [[0]*n for _ in range(count)]
    hp = [0]*count
    hp[0] = 1
    for i in range(n):
        endings[1 << i][i] = 1
    for mask in range(1, count):
        for i in bits(mask):
            value = endings[mask][i]
            hp[mask] += value
            if value:
                for j in bits(rows[i] & ~mask):
                    endings[mask | (1 << j)][j] += value
    return hp, endings


def hp(rows):
    return tables(rows)[0][-1]


def attach(rows, selected):
    n = len(rows)
    return [row | ((not (selected >> i & 1)) << n)
            for i, row in enumerate(rows)] + [selected]


def add_source(rows, selected, digit):
    n = len(rows)
    return list(rows)+[(1 << n)-1], selected | (digit << n)


def add_sink(rows, selected, digit):
    n = len(rows)
    return [row | (1 << n) for row in rows]+[0], selected | (digit << n)


def profile(rows, selected):
    validate(rows)
    H, endings = tables(rows)
    _, starts = tables(reverse(rows))
    full = (1 << len(rows))-1
    B = sum(H[a]*H[full ^ a] for a in range(full+1))
    R = sum(starts[a][i]*H[full ^ a]
            for a in range(1, full+1) for i in bits(a & selected))
    L = sum(endings[a][i]*H[full ^ a]
            for a in range(1, full+1) for i in bits(a & selected))
    s = sum(starts[-1][i] for i in bits(selected))
    e = sum(endings[-1][i] for i in bits(selected))
    F = hp(attach(rows, selected))
    return dict(h=H[-1], F=F, B=B, R=R, L=L, s=s, e=e, b=F-R+e)


def source_law(p, a):
    h, F, B, R, L, s, e = (p[k] for k in ('h', 'F', 'B', 'R', 'L', 's', 'e'))
    return dict(h=h, F=F+a*(B-L-s), B=2*B, R=R+a*B,
                L=2*L+a*h, s=a*h, e=e, b=p['b']-a*(L+s))


def sink_law(p, a):
    h, F, B, R, L, s, e = (p[k] for k in ('h', 'F', 'B', 'R', 'L', 's', 'e'))
    FF = F+(1-a)*(R+e)
    RR = 2*R+a*h
    ee = a*h
    # The source-dual laws for the unused boundary coordinates are also tested.
    LL = L+a*B
    return dict(h=h, F=FF, B=2*B, R=RR, L=LL, s=s, e=ee, b=FF-RR+ee)


def source_word(rows, selected, word):
    for digit in word:
        rows, selected = add_source(rows, selected, digit)
    return rows, selected


def sink_word(rows, selected, word):
    for digit in word:
        rows, selected = add_sink(rows, selected, digit)
    return rows, selected


def decrement(p, word):
    J = sum(a << j for j, a in enumerate(word))
    pairs = sum(word[i]*word[j]*((1 << (j-1-i))+(j == i+1))
                for i in range(len(word)) for j in range(i+1, len(word)))
    return p['L']*J+p['s']*(word[0] if word else 0)+p['h']*pairs


def strong(rows):
    full = (1 << len(rows))-1
    for rr in (rows, reverse(rows)):
        seen = 1
        while True:
            new = seen
            for i in bits(seen):
                new |= rr[i]
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


def literal(rows):
    return sum(all(rows[p[i]] >> p[i+1] & 1 for i in range(len(rows)-1))
               for p in permutations(range(len(rows))))


def is_power_two(value):
    return value > 0 and value & (value-1) == 0


def main():
    # Transitive, cyclic, weak, and holey-core controls are all explicit.
    controls = [([0], 1), ([2, 0], 1), ([2, 4, 1], 1),
                ([2, 12, 9, 1], 5), ([18, 8, 3, 5, 14], 6)]
    source_checks = sink_checks = literal_checks = word_checks = mixed_checks = 0
    for rows, selected in controls:
        p = profile(rows, selected)
        need(p['R'] >= p['h'] and p['L'] >= p['h'], 'positive boundary bounds')
        need(p['R'] < p['B'] and p['L'] <= p['B']-p['h'], 'total boundary bounds')
        for a in (0, 1):
            rr, ss = add_source(rows, selected, a)
            need(profile(rr, ss) == source_law(p, a), 'physical source law')
            source_checks += 1
            rr, ss = add_sink(rows, selected, a)
            need(profile(rr, ss) == sink_law(p, a), 'physical sink law')
            sink_checks += 1
        for length in range(5):
            for word in product((0, 1), repeat=length):
                J = sum(a << j for j, a in enumerate(word))
                pp = p
                for a in word:
                    pp = source_law(pp, a)
                need(pp['R'] == p['R']+p['B']*J, 'source response formula')
                D = decrement(p, word)
                need(pp['b'] == p['b']-D, 'exact decrement')
                need(p['L']*J <= D <= (p['L']+p['h'])*J+p['s'], 'decrement bounds')
                # Every source word of length <=2 receives full physical profile checks.
                if length <= 2:
                    rr, ss = source_word(rows, selected, word)
                    need(profile(rr, ss) == pp, 'source word physical interface')
                for t in range(4):
                    for sink_digits in product((0, 1), repeat=t):
                        fin = pp
                        for a in sink_digits:
                            fin = sink_law(fin, a)
                        A = sum(a << (t-1-j) for j, a in enumerate(sink_digits))
                        need(fin['R'] == (1 << t)*(p['R']+p['B']*J)+p['h']*A,
                             'normal-form response')
                        need(is_power_two(fin['R']) == (A == 0 and is_power_two(pp['R'])),
                             'exact dyadic return criterion')
                        word_checks += 1
        # Every sequence of 3 source/sink types and digits, including interleavings.
        for ops in product(((0, 0), (0, 1), (1, 0), (1, 1)), repeat=3):
            rr, ss = rows, selected
            pp = p
            for kind, a in ops:
                rr, ss = (add_source if kind == 0 else add_sink)(rr, ss, a)
                pp = (source_law if kind == 0 else sink_law)(pp, a)
            left = [a for kind, a in ops if kind == 0]
            right = [a for kind, a in ops if kind == 1]
            nn, ns = source_word(rows, selected, left)
            nn, ns = sink_word(nn, ns, right)
            need(profile(rr, ss) == profile(nn, ns) == pp, 'mixed-operation normal form')
            mixed_checks += 1
        if len(rows) <= 4:
            for w in ((), (1,), (1, 0), (0, 1)):
                rr, ss = source_word(rows, selected, w)
                child = attach(rr, ss)
                need(literal(child) == hp(child), 'literal Hamilton permutation control')
                literal_checks += 1

    seed = [12, 65, 26, 50, 35, 7, 61]
    selected = 72
    p = profile(seed, selected)
    need(p == dict(h=129, F=551, B=1536, R=512, L=335, s=47, e=32, b=71),
         'explicit seed full interface')
    need(strong(seed), 'seed strong')
    family = []
    for n in range(1, 5):
        word = tuple([1]+[a for _ in range(n-1) for a in (0, 1)])
        rr, ss = source_word(seed, selected, word)
        b = 24+86*n-421*((4**n-1)//3)
        pp = p
        for a in word:
            pp = source_law(pp, a)
        need(pp['R'] == 1 << (9+2*n) and pp['b'] == b, 'all-n formula finite controls')
        for t in (1, 2) if n <= 3 else (1,):
            out, mask = sink_word(rr, ss, (0,)*t)
            child = attach(out, mask)
            validate(child)
            k = 9+2*n+t
            value = hp(child)
            need(value == (1 << k)+b, 'physical family H')
            need(len(child) == k-2 and strong(child), 'physical family order and strongness')
            family.append(dict(n=n, k=k, order=len(child), H=value, offset=b,
                               rows=child, source_digits=list(word), zero_sinks=t))
    # A return that was impossible under sink-only padding actually occurs.
    need(family[0]['H'] == 3785 and family[0]['offset'] == -311, 'positive escape control')
    # Pure-power B creates no positive source return: h=5 core has B=32.
    hostile = profile([2, 12, 9, 1], 5)
    need(hostile['B'] == 32, 'power B hostile input')
    for J in range(1, 1000):
        need(not is_power_two(hostile['R']+hostile['B']*J), 'power B no-return control')
    # The simpler V13 returns are spaced by 12, not arbitrary exponents.
    v13 = profile([18, 8, 3, 5, 14], 6)
    need(v13['R'] == 32 and v13['B'] == 104, 'period-12 hostile input')
    for q in range(5, 100):
        need(((1 << q)-32) % 104 == 0 if (q-5) % 12 == 0 else ((1 << q)-32) % 104 != 0,
             'period-12 exact return test')
    # All-order formulas are also evaluated at large integer n without graph expansion.
    huge = []
    previous = None
    for n in range(1, 101):
        b = 24+86*n-421*((4**n-1)//3)
        if previous is not None:
            need(previous-b == 421*4**(n-1)-86, 'exact lacunary family gap')
        previous = b
        if n in (10, 100):
            huge.append(dict(n=n, response=512*4**n, offset=b))
    print(json.dumps(dict(
        status='PASS: finite exact checks; universal scope is proved in PROOF.md',
        source_law_physical_checks=source_checks,
        sink_law_physical_checks=sink_checks,
        two_sided_response_word_checks=word_checks,
        fully_physical_mixed_normal_form_checks=mixed_checks,
        literal_permutation_controls=literal_checks,
        seed_rows=seed, seed_mask=selected, seed_profile=p,
        strong_family_physical_checks=len(family), family=family,
        large_integer_formula_controls=huge,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    ), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
