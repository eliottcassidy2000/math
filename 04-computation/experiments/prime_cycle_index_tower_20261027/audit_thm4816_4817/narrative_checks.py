"""Histogram and narrative checks on the type tables (klein census n<=10, auditor tables 11-12, producer/auditor 13)."""
import os, sys
from collections import Counter
from math import gcd
from functools import reduce
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from audit_ledger import parse_types, KLEIN, PROD
src13 = 'my_nr13.txt' if os.path.exists(os.path.join(HERE, 'my_nr13.txt')) else os.path.join(PROD, 'nr13.txt')
tabs = {k: [r for r in parse_types(os.path.join(KLEIN, f'prime_types_n{k}.txt')) if r[2] > 1] for k in range(3, 11)}
tabs[11] = parse_types(os.path.join(HERE, 'my_nr11.txt'))
tabs[12] = parse_types(os.path.join(HERE, 'my_nr12.txt'))
tabs[13] = parse_types(os.path.join(HERE, src13) if not os.path.isabs(src13) else src13)
lcm = lambda a, b: a * b // gcd(a, b)
for k in range(3, 14):
    rows = tabs[k]
    hist = Counter(); sc = 0; orders = Counter(); z9 = 0; z32 = 0
    for (kk, cnt, G, aut, s, anti) in rows:
        hist[G] += cnt; sc += cnt * s
        for part, mult in aut:
            o = reduce(lcm, [l for l, m in part], 1)
            orders[o] += cnt * mult
        if G == 9:
            if any(reduce(lcm, [l for l, m in part], 1) == 9 for part, _ in aut): z9 += cnt
            else: z32 += cnt
        # anti-automorphism 2 mod 4 law and <= 1 fixed point; aut cycles odd
        for part, mult in anti:
            assert all(l == 1 or l % 4 == 2 for l, m in part) and dict(part).get(1, 0) <= 1, (k, part)
        for part, mult in aut:
            assert all(l % 2 == 1 for l, m in part)
            # cycle lemma for every element of prime order
            o = reduce(lcm, [l for l, m in part], 1)
            if o > 1 and all(o % q for q in range(2, o)):
                c = dict(part).get(o, 0); f = dict(part).get(1, 0)
                assert c >= 2 or f == 0, (k, part)
    el7 = sum(v for o, v in orders.items() if o % 7 == 0)
    el9 = {str(part): cnt for (kk, cnt, G, aut, s, anti) in rows for part, mult in aut if reduce(lcm, [l for l, m in part], 1) % 9 == 0}
    print(f"n={k}: non-rigid {sum(hist.values())} (SC {sc}); |Aut| hist {dict(sorted(hist.items()))}; |Aut|=9: Z_9 {z9}, Z_3^2 {z32}; "
          f"element orders (weighted by classes) {dict(sorted(orders.items()))}; elements of order divisible by 7: {el7}; cycle types of order-9 elements: {sorted(set(el9))}")
print("13 from:", src13)
