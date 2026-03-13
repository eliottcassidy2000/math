# THM-143: General Co-occurrence Linearity for Interval Tournaments

**Status:** PROVED
**Author:** kind-pasteur-2026-03-12-S59c
**Dependencies:** Vandermonde's identity, composition counting, inclusion-exclusion

## Statement

For the Interval tournament on Z_p with connection set S = {1, ..., m} where m = (p-1)/2
and any odd cycle length k in {3, 5, ..., p}, the vertex-set co-occurrence:

    co_occ_k(d) = #{k-vertex sets containing {0, d} that support a Hamiltonian cycle}

is exactly linear in d = min(v, p-v):

    co_occ_k(d) = a_k + b_k * d   for d in {1, ..., m}

with:

    Slope:     b_k = C(m-2, k-3)
    Intercept: a_k = C(2m-1, k-2) - C(m-1, k-2) - m * C(m-2, k-3)

Equivalently:

    co_occ_k(d) = C(2m-1, k-2) - C(m-1, k-2) - (m-d) * C(m-2, k-3)

## Special Cases

- k=3: co_occ_3(d) = d (slope 1, intercept 0). Recovers THM-141.
- k=5: co_occ_5(d) = C(m-2,2)*d + [C(2m-1,3) - C(m-1,3) - m*C(m-2,2)]
- k=p: co_occ_p(d) = 1 for all d (the unique p-vertex set is all of Z_p)
- Slope b_k = 0 when m < k-1, i.e., p < 2k-1. Co-occurrence is then constant.

## Proof

### Step 1: Gap Characterization (LEM-143a)

**Lemma:** A k-vertex set V on Z_p supports a Hamiltonian cycle in the Interval
tournament iff, when vertices are listed in natural cyclic order, ALL k consecutive
gaps are <= m.

**Proof:**
(⇒, "if"): If all gaps <= m, visiting vertices in natural cyclic order gives a
directed Hamiltonian cycle (each gap is in S = {1,...,m}). ✓

(⇐, "only if"): Suppose some natural gap g_i = v_{i+1} - v_i > m. The forward
cone of v_i is {v_i+1, ..., v_i+m} mod p. Since g_i > m, this arc is strictly
contained within the big gap (v_i, v_{i+1}), which contains no vertices of V.
Therefore v_i has NO outgoing arcs to any vertex in V, and no Hamiltonian cycle
can exist (every vertex needs a successor). ✓

Also verified empirically: all vertex sets with wrap-2 or wrap-3 directed cycles
also have wrap-1 cycles (tested at p=7,11,13 for k=3,5,7).

### Step 2: Composition Decomposition

For {0, d} in V, the remaining k-2 vertices divide into two groups:
- a vertices strictly between 0 and d (in the arc of length d)
- b = k-2-a vertices strictly between d and 0 (in the arc of length p-d)

The number of valid placements for a vertices in the arc of length d, creating
a+1 sub-gaps each in {1,...,m} summing to d:

    f(d, a+1) = #{compositions of d into a+1 parts in {1,...,m}}

Since d <= m, no part can exceed m (the maximum single part is d-a <= m), so:

    f(d, a+1) = C(d-1, a)    [unrestricted stars-and-bars]

For the arc of length p-d = 2m+1-d >= m+1, placing b vertices creates b+1 sub-gaps
summing to p-d. Using inclusion-exclusion (only one correction term since p-d-2m < 0):

    f(p-d, j) = C(2m-d, j-1) - j * C(m-d, j-1)

### Step 3: Vandermonde Summation

co_occ_k(d) = sum_{a=0}^{k-2} C(d-1, a) * [C(2m-d, k-2-a) - (k-1-a)*C(m-d, k-2-a)]
            = S1(d) - S2(d)

where:

S1(d) = sum_a C(d-1, a) * C(2m-d, k-2-a)
       = C((d-1)+(2m-d), k-2)    [Vandermonde]
       = C(2m-1, k-2)            [CONSTANT in d!]

S2(d) = sum_a (k-1-a) * C(d-1, a) * C(m-d, k-2-a)

Splitting (k-1-a) = (k-2-a) + 1:

S2 = sum_a (k-2-a)*C(d-1,a)*C(m-d,k-2-a) + sum_a C(d-1,a)*C(m-d,k-2-a)

Second sum: Vandermonde gives C(m-1, k-2).

First sum: absorption identity (k-2-a)*C(m-d, k-2-a) = (m-d)*C(m-d-1, k-3-a), then
Vandermonde gives (m-d) * C(m-2, k-3).

Therefore:

    S2(d) = (m-d)*C(m-2, k-3) + C(m-1, k-2)

And:

    co_occ_k(d) = C(2m-1, k-2) - (m-d)*C(m-2, k-3) - C(m-1, k-2)
                = [C(2m-1,k-2) - C(m-1,k-2) - m*C(m-2,k-3)] + d*C(m-2,k-3)

This is LINEAR in d with slope C(m-2, k-3).  QED.

## Verification

Exact match at all tested (p, k) pairs:

| p  | m | k | Slope C(m-2,k-3) | Intercept (formula) | Verified |
|----|---|---|-------------------|---------------------|----------|
| 7  | 3 | 3 | C(1,0)=1          | 0                   | YES      |
| 7  | 3 | 5 | C(1,2)=0          | 10                  | YES      |
| 11 | 5 | 3 | C(3,0)=1          | 0                   | YES      |
| 11 | 5 | 5 | C(3,2)=3          | 65                  | YES      |
| 11 | 5 | 7 | C(3,4)=0          | 126                 | YES      |
| 11 | 5 | 9 | C(3,6)=0          | 36                  | YES      |
| 13 | 6 | 3 | C(4,0)=1          | 0                   | YES      |
| 13 | 6 | 5 | C(4,2)=6          | 119                 | YES      |
| 13 | 6 | 7 | C(4,4)=1          | 455                 | YES      |
| 13 | 6 | 9 | C(4,6)=0          | 330                 | YES      |
| 17 | 8 | 3 | C(6,0)=1          | 0                   | YES      |
| 17 | 8 | 5 | C(6,2)=15         | 300                 | YES      |
| 17 | 8 | 7 | C(6,4)=15         | 2862                | YES      |
| 17 | 8 | 9 | C(6,6)=1          | 6426                | YES      |
| 17 | 8 | 11| C(6,8)=0          | 5005                | YES      |
| 17 | 8 | 13| C(6,10)=0         | 1365                | YES      |
| 19 | 9 | 3 | C(7,0)=1          | 0                   | YES      |
| 19 | 9 | 5 | C(7,2)=21         | 435                 | YES      |
| 19 | 9 | 7 | C(7,4)=35         | 5817                | YES      |
| 19 | 9 | 9 | C(7,6)=7          | 19377               | YES      |
| 19 | 9 | 11| C(7,8)=0          | 24310               | YES      |
| 19 | 9 | 13| C(7,10)=0         | 12376               | YES      |

**22 out of 22 test cases pass.**

## Correction

The earlier conjecture (HYP-524) proposed b_k = C(m-(k-1)/2, (k-1)/2). This agrees
with C(m-2, k-3) for k=3, k=5, and for all k when m <= 6. It DISAGREES at p=19 (m=9):
- k=7: correct C(7,4)=35 vs old C(6,3)=20
- k=9: correct C(7,6)=7 vs old C(5,4)=5

## Connection to Phase Transition

The slope b_k = C(m-2, k-3) is nonzero iff m >= k-1, i.e., p >= 2k-1.
For large p, b_k ~ m^{k-3}/(k-3)!, so longer cycles contribute LARGER co-occurrence
variance, amplifying the Interval tournament's disjointness advantage.

This is a key mechanism behind the Paley-to-Interval phase transition at p ~ 13:
the co-occurrence gradient drives disjointness excess that grows polynomially in p.

## Related

- THM-141 (k=3 case): co_occ_3(d) = d
- THM-142: Disjointness excess formula for 3-3 pairs
- HYP-524: CORRECTED — slope is C(m-2,k-3), not C(m-(k-1)/2,(k-1)/2)
- HYP-525: Paley co_occ constant (trivially follows from vertex-transitivity + edge-transitivity)
