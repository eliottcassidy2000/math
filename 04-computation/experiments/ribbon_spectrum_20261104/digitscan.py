#!/usr/bin/env python3
"""digitscan.py -- search for digit blocks in the span-2 strong alphabet {1},{2},{1,2}:
blocks B (words of length <= 3) such that H(uB) - alpha*H(u) takes few values over all prefixes u
(|u| <= L), for alpha in {2, 3}.  Also records the affine relation statistics."""
import sys, itertools
from collections import Counter
from ribbon import transfer_system, letters, word_str

w = 2
order, states, mats, final = transfer_system(w)
F = frozenset
alph = [F([1]), F([2]), F([1, 2])]
finalset = set(final)

def apply(M, vd):
    nv = {}
    for (j, i), mult in M.items():
        if i in vd:
            nv[j] = nv.get(j, 0) + mult * vd[i]
    return nv

def Hof(vd):
    return sum(c for j, c in vd.items() if j in finalset)

L = int(sys.argv[1]) if len(sys.argv) > 1 else 9
# all prefix vectors (deduped) for |u| in [3, L]
prefixes = {}
cur = {((0, 1),): ()}
for n in range(1, L + 1):
    nxt = {}
    for vec, word in cur.items():
        vd = dict(vec)
        for a in alph:
            nv = apply(mats[a], vd)
            nvt = tuple(sorted(nv.items()))
            if nvt not in nxt:
                nxt[nvt] = word + (a,)
    cur = nxt
    if n >= 3:
        for vec, word in cur.items():
            prefixes[vec] = word
print(f"{len(prefixes)} distinct prefix vectors with 3 <= |u| <= {L}")
blocks = []
for k in (1, 2, 3):
    for B in itertools.product(alph, repeat=k):
        blocks.append(B)
results = []
for B in blocks:
    # product matrix as a function: apply letters in order
    resid = {2: Counter(), 3: Counter()}
    for vec, word in prefixes.items():
        vd = dict(vec)
        h = Hof(vd)
        v2 = vd
        for a in B:
            v2 = apply(mats[a], v2)
        h2 = Hof(v2)
        for alpha in (2, 3):
            resid[alpha][h2 - alpha * h] += 1
    for alpha in (2, 3):
        c = resid[alpha]
        results.append((len(c), alpha, B, c))
results.sort(key=lambda t: (t[0], len(t[2])))
for nres, alpha, B, c in results[:16]:
    common = c.most_common(6)
    print(f"block {word_str(list(B), w):>10} alpha={alpha}: {nres} distinct residuals; most common {common}")
