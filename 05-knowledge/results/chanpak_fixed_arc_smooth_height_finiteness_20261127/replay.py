"""Finite controls for the all-budget skeleton-parity lemma.
Universal finiteness uses the proof and ESS, not these control cutoffs.
"""
from collections import Counter
from itertools import combinations, permutations

def patterns(m, b):
    for arcs in combinations(combinations(range(m), 2), b):
        if len(set(sum(arcs, ()))) == m:
            yield frozenset(arcs)

def coefficient_parities(m, arcs):
    support = set()
    for p in permutations(range(m)):
        if any(i > j and (j,i) not in arcs for i,j in zip(p,p[1:])):
            continue
        slots = list(zip((-1,)+p, p+(m,)))
        required = [i for i,e in enumerate(slots) if e in arcs]
        for omit in range(1 << len(required)):
            D = {r for j,r in enumerate(required) if (omit >> j) & 1}
            v = tuple(sum(i <= g < j for r,(i,j) in enumerate(slots) if r not in D)
                      for g in range(m-1))
            support.symmetric_difference_update({v})
    return support

def control(b):
    table = []
    for m in range(2,2*b+1):
        patterns_count = masks_count = 0
        for arcs in patterns(m,b):
            patterns_count += 1
            support = coefficient_parities(m,arcs)
            for mask in range(1 << (m-1)):
                # True means a gap is positive; adjacent reverse pairs are
                # required to have positive length only in a minimum order.
                if any(j == i+1 and not ((mask >> i) & 1) for i,j in arcs):
                    continue
                odd = set()
                for v in support:
                    w = tuple(v[g] for g in range(m-1) if (mask >> g) & 1)
                    if all(x%2 for x in w):
                        odd.symmetric_difference_update({w})
                assert odd == {(1,)*mask.bit_count()}, (m,arcs,mask,odd)
                masks_count += 1
        if patterns_count:
            table.append((m,patterns_count,masks_count))
    return table

def paley_budget():
    n=7
    adj=[sum(1<<j for j in range(n) if (j-i)%n in (1,2,4)) for i in range(n)]
    dp=[0]*(1<<n)
    for S in range(1,1<<n):
        dp[S]=min(dp[S^(1<<v)]+(adj[v]&S).bit_count()
                  for v in range(n) if (S>>v)&1)
    order=[4,1,5,2,6,3,0]
    witness=sum((j-i)%n in (1,2,4) for k,i in enumerate(order) for j in order[:k])
    assert dp[-1]==witness==7
    return witness

if __name__ == "__main__":
    print("Paley-7 minimum backward arcs:",paley_budget())
    for b in range(1,5):
        rows = control(b)
        print(b, rows, "total", sum(row[2] for row in rows))
