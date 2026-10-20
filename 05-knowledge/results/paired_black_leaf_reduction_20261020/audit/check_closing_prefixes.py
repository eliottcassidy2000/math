#!/usr/bin/env python3
"""Fresh literal closure and nonautomorphic-successor checks of prefix families."""
import itertools
import json
from pathlib import Path


def witness(n, S, depth=0):
    assert depth <= 2
    m = n // 2
    if 2 not in S and 3 in S:
        return [2, 0, 3, 1] + list(range(4, n))
    for k in range(4, m+1):
        if k in S and k-2 not in S:
            return [k-1] + list(range(1, k-1)) + [0] + list(range(k, n))
    if 2 in S:
        for k in range(3, m):
            if k not in S:
                if k % 2:
                    return [1,k,0] + list(range(2,k)) + list(range(k+1,n))
                return [1] + list(range(3,k+1)) + [0,2] + list(range(k+1,n))
    # With the local forbidden configurations absent, only A/B/C survive.
    positive = set(range(1,m+1)) & S
    assert positive in ({1}, set(range(1,m)), set(range(1,m+1)))
    step = n-2 if positive == {1} else 2
    assert step in S
    new_S = {j for j in range(1,n) if step*j % n in S}
    return [step*j % n for j in witness(n,new_S,depth+1)]


def verify(n, S, P):
    assert sorted(P) == list(range(n))
    assert all((v-u) % n in S for u,v in zip(P,P[1:] + P[:1]))
    r = dict(zip(P,P[::-1]))
    c = dict(zip(P,P[1:] + P[:1]))
    bad_reverse = [(i,j) for i in range(n) for j in range(n)
                   if (j-i) % n in S and (r[i]-r[j]) % n not in S]
    bad_successor = [(i,j) for i in range(n) for j in range(n)
                     if (j-i) % n in S and (c[j]-c[i]) % n not in S]
    assert bad_reverse and bad_successor


counts = {}
small7 = [
    ({1,4,5}, [0,1,2,6,4,5,3]),
    ({1,2,4}, [0,1,2,4,5,6,3]),
    ({1,3,5}, [0,1,2,5,3,6,4]),
    ({1,2,3}, [0,1,3,4,6,2,5]),
]
for S,P in small7:
    verify(7,S,P)
counts["order7"] = 4
for n in range(9,32,2):
    m = n//2
    count = 0
    for bits in itertools.product((0,1),repeat=m-1):
        S = {1} | {d if bit else n-d for d,bit in zip(range(2,m+1),bits)}
        verify(n,S,witness(n,S))
        count += 1
    counts[str(n)] = count
structured = 0
for m in range(4,101):
    n = 2*m+1
    for low in ({1}, set(range(1,m)), set(range(1,m+1))):
        S = low | {n-d for d in range(1,m+1) if d not in low}
        verify(n,S,witness(n,S))
        structured += 1
counts["residual_patterns_m4_to100"] = structured
print(json.dumps(counts,indent=2))
Path(__file__).with_suffix(".json").write_text(json.dumps(counts,indent=2)+"\n")
print("ALL CLOSING PREFIX AND NONAUTOMORPHIC SUCCESSOR CHECKS PASSED")
