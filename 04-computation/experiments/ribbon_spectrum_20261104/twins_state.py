#!/usr/bin/env python3
"""twins_state.py -- are the profile state vectors (or at least all future H values) of u{1}{2} and
u{1}{1,2} equal for every prefix u?  Also: which other letter pairs are state-equivalent after a given
preceding letter?"""
import itertools
from ribbon import transfer_system, letters, count_with_matrices, word_str
w = 2
order, states, mats, final = transfer_system(w)
F = frozenset
L = letters(w)
finalset = set(final)
def apply(M, vd):
    nv = {}
    for (j, i), mult in M.items():
        if i in vd:
            nv[j] = nv.get(j, 0) + mult * vd[i]
    return nv
def Hof(vd): return sum(c for j, c in vd.items() if j in finalset)
prefixes = []
cur = [{0: 1}]
for n in range(1, 9):
    nxt = []
    for vd in cur:
        for a in L:
            nxt.append(apply(mats[a], vd))
    cur = nxt
    prefixes.extend(cur)
print(len(prefixes), "prefix vectors (all words, n<=8)")
pairs = [(a, b) for a in L for b in L if sorted(a) < sorted(b)]
for pre in L:
    for a, b in pairs:
        eqstate = all(apply(mats[a], apply(mats[pre], vd)) == apply(mats[b], apply(mats[pre], vd)) for vd in prefixes)
        # weaker: equal H for all suffixes of length <= 3
        eqH = True
        if not eqstate:
            for vd in prefixes[:3000]:
                va = apply(mats[a], apply(mats[pre], vd)); vb = apply(mats[b], apply(mats[pre], vd))
                for suf in itertools.chain.from_iterable(itertools.product(L, repeat=k) for k in range(0, 4)):
                    xa, xb = va, vb
                    for s in suf:
                        xa = apply(mats[s], xa); xb = apply(mats[s], xb)
                    if Hof(xa) != Hof(xb):
                        eqH = False; break
                if not eqH: break
        if eqstate or eqH:
            print(f"after letter {sorted(pre)}: letters {sorted(a)} and {sorted(b)} are", "STATE-equal" if eqstate else "H-equivalent for all suffixes (len<=3 tested)")
print("done")
