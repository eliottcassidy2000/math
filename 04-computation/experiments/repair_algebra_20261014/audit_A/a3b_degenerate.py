#!/usr/bin/env python3
"""Audit A: which minimal vanishing {2,3}-sums occur inside primitive collisions at -1 (canonical pairs, exps <= E)?
For each pair, the padded difference is a vanishing sum of 2r terms +-3^a 2^b.  We list its minimal vanishing
subsums of smallest size (normalised by gcd and sign) and count, over all pairs, how often each normalised type is
the smallest minimal subsum.  Usage: python3 a3b_degenerate.py [E=16]"""
import sys, itertools, math
from collections import Counter, defaultdict

E = int(sys.argv[1]) if len(sys.argv) > 1 else 16


def K_of(s):
    K = -1
    for e in s:
        K = 3 * K + (1 << e)
    return K


reps = defaultdict(list)
for n in range(E + 1):
    for s in itertools.combinations(range(2, E + 1), n):
        reps[K_of(s)].append(s)


def norm(terms):
    g = 0
    for x in terms:
        g = math.gcd(g, abs(x))
    t = sorted((x // g for x in terms), key=lambda v: (-abs(v), -v))
    if t[0] < 0:
        t = [-v for v in t]
    return tuple(sorted(t, reverse=True))


def smallest_minimal_subsum(terms):
    n = len(terms)
    for size in range(3, n + 1):
        for idx in itertools.combinations(range(n), size):
            if sum(terms[i] for i in idx) == 0:
                return [terms[i] for i in idx]
    return None


types = Counter(); deg = 0; tot = 0; first_examples = {}
for k, rs in sorted(reps.items()):
    if k == -1:
        continue
    for s, t in itertools.combinations(rs, 2):
        if s[-1] == t[-1]:
            continue
        if len(s) > len(t):
            s, t = t, s
        sp = (1,) * (len(t) - len(s)) + s
        L = len(t)
        terms = []
        for j, (x, y) in enumerate(zip(sp, t)):
            if x != y:
                terms += [3 ** (L - 1 - j) * 2 ** x, -(3 ** (L - 1 - j)) * 2 ** y]
        sub = smallest_minimal_subsum(terms)
        tot += 1
        if len(sub) < len(terms):
            deg += 1
        ty = norm(sub)
        types[(len(sub) < len(terms), ty)] += 1
        first_examples.setdefault((len(sub) < len(terms), ty), (k, s, t))
print(f"E <= {E}: {tot} canonical primitive pairs, {deg} degenerate")
print("most frequent smallest minimal subsums (degenerate?, normalised terms): count, first example")
for (d, ty), c in types.most_common(14):
    print(f"   {d!s:5} {ty}: {c}   e.g. {first_examples[(d, ty)]}")
