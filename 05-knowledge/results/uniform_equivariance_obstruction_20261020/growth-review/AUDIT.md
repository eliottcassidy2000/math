# Independent audit: prime-power Hall-deficit growth

Date: 2026-10-20 UTC.

## Verdict and reviewed artifacts

**PASS.** The frozen producer artifact was read in full and its SHA-256
verified. The proof gives the stated
complete-bipartite lower bound at every integer n >= 5, and a quantitative
square-root lower bound along an explicit infinite subsequence. It does not
claim a square-root bound at every n, an exact matching deficit, or an
obstruction to arbitrary bijections between unlabelled type sets.

Frozen producer artifact:

- `session12-hall-growth/PRIME_POWER_HALL_GROWTH.md`
- SHA-256: `e600124b9147b138ae3a009365acb724bec4c8b09b15389b0553c8fdea7a2bdf`

Published dependency read directly:

- `05-knowledge/results/uniform_equivariance_obstruction_20261020/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md`
- SHA-256: `1ffbff76190de3c2603b9512d98aa5192c5454d5d97da358c101fd62ce562db2`
- Exact section: “Corollaries: a joint-prime obstruction and an explicit Hall family.”

The producer's computational code was not read or reused. This audit and all
its controls were written outside the repository checkout. Neither repository
files nor producer files were edited.

## 1. The exact inherited Hall inequality

Let A and B be respectively the source and target isomorphism-type sets, in
either of the two settings:

1. tournaments and Royle-sign-even graphs;
2. tournament switching classes and Royle-sign-even Eulerian graphs.

They have a common cardinality N. The stabilizer-containment graph has an edge
from source type a to target type X precisely when a source stabilizer embeds
in a conjugate of the target stabilizer, in their given actions on n vertices.
Let delta = N - nu, where nu is its maximum matching cardinality.

Write A_p and B_p for types whose stabilizers contain a full Sylow-p subgroup
of ambient S_n. The published proof establishes:

- |A_p| = |B_p| for each odd prime p;
- A_3 and A_5 are disjoint for n >= 5;
- stabilizer containment preserves the full ambient Sylow-p property.

For S = A_3 union A_5, its neighbor set is contained in B_3 union B_5. Hence

    |S| - |neighbors(S)|
      >= |A_3| + |A_5| - |B_3 union B_5|
       = |B_3 intersection B_5|.

Every matching leaves at least that many elements of S unmatched. Therefore

    delta >= |B_3 intersection B_5|.

This is a direct matching-count inequality; it requires no perfect-matching
assumption. Counting distinct target types in the intersection is sufficient.
The argument needs the published one-prime cardinality equality, not merely
equality of total prime valuations. That cardinality equality is proved in
the cited section from the inverse Walsh matrix and its p-unit determinant.
The published source-disjointness uses the no-double-transposition lemma and
the full-3-part/full-5-part incompatibility, both valid for n >= 5.

## 2. Complete bipartite targets, including full automorphism groups

Take 0 <= k <= floor(n/2), and let X_k = K_{k,n-k}.

For k = 0, X_k is the empty graph on n vertices. Its full automorphism group
is S_n, it is sign-even, and all degrees are zero. Thus it is Eulerian at
both even and odd n. This endpoint must not be discarded by a parity filter.

For 0 < k < n/2, the two parts have distinct degrees. Every automorphism
preserves the parts, and all permutations within either part are allowed.
Consequently the full group is S_k x S_(n-k), of order k!(n-k)!. Orient all
edges from the first part to the second. Every such automorphism reverses
zero edges, proving sign-evenness for the full automorphism group.

For 0 < k = n/2, the connected bipartite graph has its two parts determined
up to interchange. Its full group is (S_k x S_k) semidirect C_2, of order
2(k!)^2. A part interchange reverses all k^2 edges, and the within-part group
reverses none. Therefore this balanced graph is sign-even exactly when k is
even. The extra factor 2 does not affect the 3- or 5-adic arguments.

For positive k, the graph is Eulerian exactly when both k and n-k are even.
Its edge count k(n-k) is strictly increasing as k ranges from 0 to floor(n/2),
so different permitted indices give nonisomorphic target graphs.

In all three cases, for p = 3 or 5,

    v_p(|Aut(X_k)|) = v_p(n!) - v_p(binomial(n,k)).

Thus X_k contains a full ambient Sylow-3 and a full ambient Sylow-5 subgroup
if and only if gcd(binomial(n,k),15) = 1. Equality of these valuations is
equivalent to the full-ambient-Sylow condition by Sylow's theorem; it is not
merely the existence of an element of order 3 or 5.

## 3. General-n lower bound

For every n >= 5, define

    I_T(n) = {k : 0 <= k <= floor(n/2),
                    gcd(binomial(n,k),15) = 1,
                    and (2k != n or k is even)}.

Define

    I_S(n) = {k in I_T(n) : k = 0 or (k and n-k are both even)}.

Then the independent reasoning above proves

    delta_T(n) >= |I_T(n)|,
    delta_S(n) >= |I_S(n)|.

For odd n, I_S(n) is exactly {0}. For even n, it is the even-index subset of
I_T(n). These are lower bounds from a selected target family, not assertions
that all target types in B_3 intersection B_5 are complete bipartite.

## 4. Prime-power family and the short cancellation proof

Let a,b >= 1 be integers and put

    n + 1 = 3^a 5^b,
    M = min(3^a,5^b).

Both n and M-1 are even, n >= 14, and M >= 3. Since n+1 >= 3M,

    2(M-1) < n.

Thus all 0 <= k < M lie strictly below n/2, avoiding balanced-part exceptions.
For every 1 <= j <= k < M,

    v_3(j) < a,       v_5(j) < b.

If two positive integers have unequal p-adic valuations, their difference
has the smaller valuation. Therefore

    v_p(n+1-j) = v_p(j),       p = 3,5.

Apply this term by term to

    binomial(n,k) = product_(j=1)^k (n+1-j)/j.

The total valuations vanish. For k = 0 the same conclusion follows from the
empty product binomial(n,0) = 1. All M target graphs are sign-even and have
both full odd Sylow parts. Exactly (M+1)/2 indices are even, including zero,
and these give Eulerian graphs. Hence

    delta_T(3^a5^b-1) >= M,
    delta_S(3^a5^b-1) >= (M+1)/2.

Taking a=b=r gives n = 15^r-1, M=3^r, for every integer r >= 1. The first
case n=14 gives the promised lower bounds 3 and 2. They need not be sharp:
the general index formula actually gives 6 and 3 at n=14. Nothing depends on
extending this family to r=0, where n=0 lies outside the inherited Hall result.

## 5. Explicit square-root-growth subsequence

For each b >= 1, choose a least such that 3^a >= 5^b. Minimality gives

    5^b <= 3^a < 3*5^b.

Here M = 5^b and n+1 = 3^a5^b < 3M^2, so

    delta_T(n) >= M > sqrt((n+1)/3),
    delta_S(n) >= (M+1)/2 > (sqrt((n+1)/3)+1)/2.

As b increases, n tends to infinity. This proves a square-root lower bound
along an explicit infinite sequence, requiring no approximation theorem for
logarithms. It does not assert such a bound for every integer n.

## 6. The digit interpretation

The frozen note's additional base-p interpretation is correct. In
F_p[x], the factorization of (1+x)^n over the base-p digits of n has one
possible contribution to a coefficient precisely when each digit of k is
at most the corresponding digit of n. All digits selected from a factor
are below p, so the exponent representation has no carries and is unique.
The coefficient is a product of nonzero binomial coefficients with upper
entries below p. The binomial is a unit modulo p exactly under that digit
condition. This is equivalently the no-borrow condition for n-k, or the
no-carry condition for k+(n-k). The last a base-3 and b base-5 digits of
3^a5^b-1 are maximal. This supplies the stated second explanation.

The displayed b(n),e(n) examples are also correct, including (36,18) at
n=224 and (27,14) at n=674. They are explicitly labelled lower bounds,
not exact matching deficits.

## 7. Independent finite controls

Reproduce from the parent workspace with

    python3 -B session12-hall-growth-review/independent_controls.py

The control program uses exact integers, fresh code, and the following
explicitly bounded universes:

- It enumerates every vertex permutation for each of 25 individual graphs
  K_{k,n-k}, with 0 <= n <= 8 and 0 <= k <= floor(n/2). The actual full
  automorphism orders, edge-reversal signs, and degree parities agree with
  the formulas in Section 2. This is not a graph or tournament type census.
- It checks 2,592 complete-bipartite arithmetic cases, 5 <= n <= 100 and all
  0 <= k <= floor(n/2). Coprimality of the binomial coefficient agrees with
  full 3- and 5-adic automorphism valuations in every case.
- It additionally checks the two-base digit criterion against exact binomial
  coprimality in 5,151 cases, 0 <= n <= 100 and all 0 <= k <= n. It verifies
  the displayed b(n),e(n) values at n=224 and 674 by 451 further exact
  complete-bipartite arithmetic cases.
- It checks all constructed indices in the 16 prime-power cases 1 <= a,b <= 4,
  including exact binomial coprimality and numerator/denominator valuations.
- It checks six balanced prime-power pairs b=1,...,6, using sparse endpoint
  and midpoint samples for the p-adic identities and the exact scalar
  square-root inequality 3M^2 > n+1. It does not enumerate all k in these
  larger cases or compute any tournament census.

All checks pass. Hostile controls include K_{3,3} (sign-odd), K_{4,4}
(sign-even and Eulerian), and the empty graph at odd n. The first three
balanced-sequence cases are:

- (a,b,n,M) = (2,1,44,5), giving lower bounds 5 and 3;
- (3,2,674,25), giving 25 and 13;
- (5,3,30374,125), giving 125 and 63.

These calculations are finite controls of the proof, not evidence in place
of the all-order argument.
