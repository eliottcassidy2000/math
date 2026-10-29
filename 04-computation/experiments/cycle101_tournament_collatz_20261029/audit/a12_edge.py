#!/usr/bin/env python3
"""THM-4849 (E) without the hypothesis p not | 2^v - 3: are there (w, v, p) with p | D(wv), all letters = v
mod ord_p(2), and NOT every marking cancelling p?  (Shows the hypothesis must stay in the title/summary.)"""
from math import gcd
from itertools import product
from sympy import factorint
from sympy.ntheory import n_order
def C_of(w):
    k = len(w); s = 0; t = 0
    for j in range(k):
        t += 3 ** (k - 1 - j) * 2 ** s; s += w[j]
    return t
def D_of(w): return 2 ** sum(w) - 3 ** len(w)
found = []; total = 0
for k in range(1, 5):
    for w in product(range(1, 9), repeat=k):
        if D_of(w) == 0: continue
        for v in range(1, 13):
            m = 2 ** v - 3
            Dwv = D_of(w + (v,))
            for p in factorint(abs(m)):
                if p in (2, 3) or Dwv % p: continue
                L = n_order(2, p)
                if all((x - v) % L == 0 for x in w):
                    total += 1
                    rots = [w[i:] + w[:i] for i in range(k)]
                    allc = all(gcd(C_of(r + (v,)), Dwv) % p == 0 for r in rots)
                    if not allc and len(found) < 5: found.append((w, v, p, allc))
print(f"cases with p | 2^v-3, p | D(wv), letters = v mod ord_p(2): {total}; first ones where not every marking cancels p: {found}")
