# Independent audit: information, unit-flow normalization and involution defect

2026-10-26. **PASS for the stated scope.** Read RESULT.md and check.py;
ran the producer in normal and optimized Python, obtaining byte-identical
output, including the final complete-defect-polynomial update. This review
does not import the producer. The parent applied the requested wording
clarification: signs of an injective potential represent the **transitive**
relation, not a cyclic tournament. No proof gap was found in the class-fibre,
entropy, normalization or defect arguments.

The exact distinction between three measures is essential: uniform path
chart, uniform labelled tournament and uniform class have different masses.
The free action on Hamiltonian paths proves fibre size H/Aut; the n32
capacity inequality concerns reconstruction from three statistics, not
the complexity of computing the aggregate class count. Neither this bound
nor the rank lane's finite-sensor obstruction proves impossibility for a
single fixed arithmetic source.

For the Hodge claim, differentiating the stated unit-flow objective gives
`n phi=d` in the zero-sum gauge, so the repaired energy x/n is correct.
The transitive baseline is necessary: acyclic orientation does not make its
unit-valued edge flow an exact gradient. The cited
[Jiang–Lim–Yao–Ye paper](https://arxiv.org/abs/0811.1067) indeed describes
the orthogonal gradient/divergence-free decomposition; the complete-graph
constants used here are derived directly. Its primary abstract was checked
on 2026-10-26. No new literature claim is needed for the elementary proof.

## Strict scope of the reflection bridge, with its complete defect enumerator

Let g be an involution with r transposed vertex pairs and f=n-2r fixed
vertices. Write

    s=(binom(n,2)-r-binom(f,2))/2.

There are r fixed unordered pairs whose endpoints are exchanged; these
necessarily reverse their arcs. There are binom(f,2) pointwise-fixed
unordered pairs; these never contribute. The remaining unordered pairs
form s orbits of size two. In each such orbit there are four orientation
assignments: two are compatible with g and contribute zero, and two
contribute exactly two. All these choices are independent. Therefore

    sum_(labelled tournaments T) t^(D_g(T))
      = 2^(r+binom(f,2)+s) t^r (1+t^2)^s.                (A)

This proves both the producer's minimum r and the full defect distribution.
For n=6,r=2 it gives s=6 and exactly

    number with defect 2+2j = 512 binom(6,j), 0<=j<=6.

In particular there are 512 minimizers, and no invariant tournament.
If A is the signed unit adjacency matrix, the normalization is

    ||A-P_g A P_g^T||_F^2 = 8 D_g(T),

because each changed unordered pair contributes two squared entries of
size four. Thus pair defect 2 corresponds to squared Frobenius defect16,
not2. The producer uses pair defect and does not conflate these values.

This applies to the realized reflection because its permutation has exactly
two transpositions and two fixed vertices. It retains the vertex set and
permutation, then asks a **different** predicate: invariance of a strict
complete binary relation. It does not derive a tournament from metric
data or make the free action on four frames into a tournament automorphism.
In particular it gives neither a planar colouring lower bound nor an
arithmetic rank.

## Independent exact executable gate

The following standalone block checks (A) on every involution cycle type
through n5 and the actual n6,r2 type, directly using adjacency matrices;
it also uses subset dynamic programming for the small Hamiltonian-count
collision. It has no assertions disabled by optimized Python.

```python
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb
import json

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def matrix(n, code):
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = 1 if code & (1 << k) else -1
        a[j][i] = -a[i][j]
    return a

checked = 0
case_count = 0
six = None
for n in range(7):
    for r in range(n // 2 + 1):
        if n == 6 and r != 2:
            continue
        g = list(range(n))
        for i in range(r):
            g[2*i], g[2*i+1] = 2*i+1, 2*i
        observed = Counter()
        for code in range(1 << comb(n, 2)):
            a = matrix(n, code)
            d = sum(a[i][j] != a[g[i]][g[j]]
                    for i, j in combinations(range(n), 2))
            norm = sum((a[i][j] - a[g[i]][g[j]])**2
                       for i in range(n) for j in range(n))
            need(norm == 8*d, "Frobenius normalization")
            observed[d] += 1
            checked += 1
        f = n - 2*r
        s = (comb(n, 2) - r - comb(f, 2)) // 2
        expected = {r + 2*j: 2**(r+comb(f,2)+s)*comb(s,j)
                    for j in range(s+1)}
        need(observed == expected, "full involution defect enumerator")
        case_count += 1
        if n == 6:
            six = dict(sorted(observed.items()))

def hp(a):
    n = len(a)
    dp = {(1 << i, i): 1 for i in range(n)}
    for mask in range(1, 1 << n):
        for end in range(n):
            ways = dp.get((mask, end), 0)
            if not ways:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and a[end][nxt] == 1:
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + ways
    return sum(dp.get(((1 << n) - 1, i), 0) for i in range(n))

a, b = matrix(5, 8), matrix(5, 10)
ps = list(permutations(range(5)))
need(not any(all(a[i][j] == b[p[i]][p[j]]
                 for i in range(5) for j in range(5)) for p in ps),
     "the scalar collision must be nonisomorphic")
for t in (a, b):
    need(hp(t) == 9, "independent subset-DP Hamiltonian count")
    need(sorted(sum(x == 1 for x in row) for row in t) == [1,1,2,3,3],
         "score collision")
    need(sum(all(t[i][j] == t[p[i]][p[j]]
                 for i in range(5) for j in range(5)) for p in ps) == 1,
         "automorphism collision")
    divergence = list(map(sum, t))
    phi = [F(x, 5) for x in divergence]
    residual = [[F(t[i][j]) - phi[i] + phi[j]
                 for j in range(5)] for i in range(5)]
    need(all(sum(row) == 0 for row in residual), "normal equations")
    need(sum((phi[i]-phi[j])**2 for i,j in combinations(range(5),2))
         == F(sum(x*x for x in divergence), 5), "gradient energy")

need(F(1,2)-F(-1,2) == 1, "true K2 gradient")
need(F(1,4)-F(-1,4) != 1, "half-normalization hostile")
print(json.dumps({"status":"PASS", "involution_cases":case_count,
                  "labelled_matrix_checks":checked,
                  "n6_r2_distribution":six,
                  "independent_HP_method":"subset dynamic programming"},
                 sort_keys=True))
```

The expected n6 distribution is `{2:512,4:3072,6:7680,8:10240,10:7680,
12:3072,14:512}`. The elementary orbit proof establishes the all-order
statement; exhaustive checks only cover the declared small universe.
Both normal and optimized execution of this independent block returned
`PASS`, with 13 involution types and 36,054 labelled matrix checks. The
producer's final normal and optimized outputs also agreed byte for byte.
