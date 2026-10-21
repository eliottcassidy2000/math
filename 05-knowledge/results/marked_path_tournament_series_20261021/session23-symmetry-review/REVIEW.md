# Independent review: distinguished partial-path count series

Verdict: **PASS after one corrected explicit-example typo.** The identities,
normalized exact denominators, rooted central-orbit assertion, free-action
threshold, first correction, and raw-series non-D-finiteness arguments in the
final producer packet are mathematically sound in their stated domains.

This review is independent of the producer's implementation. The audit imports
no producer or repository counting module. It edits no repository files and
performs no Git operation. It makes no literature-priority claim.

## Reviewed artifact and scope

Producer packet: `session23-symmetry/README.md`, `check.py`, and its receipts.
Exact reviewed bytes are pinned in `review-manifest.json` alongside inherited
sources. The producer has corrected the explicit example to

    G_3(z) = (12z - 32z^2)/[(1-8z)(1-2z)].

Its earlier numerator `12z-48z^2` was inconsistent with the correct coefficient
`A_3(2)=88`; the final producer script now includes a regression for this.
The general partial-fraction formula was not affected.

The object is an isomorphism class of a tournament with an ordered directed
path on k distinct marked vertices, with `n >= k >= 1`. Isomorphisms preserve
every marked position. Converse reverses both all arcs and the marked order.
These domains and this action are essential. In particular, `k=0` is not an
instance of the displayed binomial prefactor, and pointwise two-root converse
is a different, empty fixed channel.

## Symbolic audit

Write `m=n-k`. After fixing the ordered marked labels, only `S_m` acts.
The path fixes k-1 internal arcs, leaving `binom(k-1,2)` internal arc bits;
every unmarked vertex has k independent adjacency bits to the marked labels.
This gives

    P(k+m,k) = 2^binom(k-1,2) F_m(1,k).

It also proves the labelled-scaffold version with prefactor
`2^[binom(k,2)-e]`: its shape is irrelevant only because all its vertices
remain individually distinguished and there is no quotient by scaffold
automorphisms.

For converse, put `k=2r+delta`, `m=2s+epsilon`. Mark reversal consists of r
transpositions and delta fixed points. The remaining permutation must have
cycle type `2lambda` plus epsilon fixed points, where lambda is an odd
partition of s; two fixed points are impossible. Its unordered-pair orbit
counts are, independently:

    marked-marked:       r(r+delta),
    unmarked-unmarked:   2E+(1+epsilon)K,
    marked-unmarked:     kK+r*epsilon.

The adjacent path pairs constitute r consistent orbit constraints. The
unmarked permutation centralizer is `2^K z_lambda`. These give exactly

    B(n,k)=0                                        if delta=epsilon=1,
    B(n,k)=2^[r(r-1+delta+epsilon)] F_s(2,k+epsilon)   otherwise.

The proof correctly includes anti-automorphisms that are not involutions.
The anti-invariant scaffold variant subtracts its actual number q of pinned
pair orbits rather than blindly subtracting the number of prescribed arcs.
Consistency of the pin assignments and the at-most-one-fixed-point condition
are both necessary.

For fixed m, the normalized ordinary sequence is a sum of distinct exponentials
`sum_K w_(m,K) (2^K)^k`. The occurring K are exactly `{0}` at m=0 and all
positive K at most m of the same parity as m otherwise. Existence follows
from the odd partition `(m-K+1,1,...,1)` with K parts. Each grouped weight is
strictly positive. Consequently no pole cancels, the reduced denominator is
exactly `product_K(1-2^K z)`, and the recurrence order is minimal. The converse
parity series repeats the argument with s in place of m and bases `4^K`.
The forbidden parity series has denominator 1, not that product. A positive
partial-fraction weight is not a positive analytic residue; the final proof
uses the noncancellation claim correctly.

The central-root theorem is objectwise, not merely equality of total counts.
For a self-converse tournament, average the number of fixed vertices over
the anti-automorphism coset of G=Aut(T). A G-orbit preserved by the converse
action contributes exactly one to this average; any other orbit contributes
zero. Each anti-automorphism has exactly `n mod 2` fixed vertices by its cycle
type. Thus there is exactly one preserved root orbit at odd n and none at
even n. This does not imply a unique vertex: in the directed 3-cycle all
three vertices occur as anti-involution centers and form one Aut orbit.

The weighted mass

    W(n,k) = sum_[T] h_k(T)/|Aut(T)|
           = 2^[binom(n,2)-k+1]/(n-k)!

is exactly the identity-permutation Burnside term. Every nonidentity
tournament automorphism moves at least three vertices. Hence equality
`P=W` holds for m<=2, whereas a positive `(3,1^(m-3))` term proves `P>W`
for every m>=3. This is a strict all-parameter statement, not just a first
counterexample. At m>=3 the ratio of the unique K=m-2 weight to the K=m
identity weight is

    c_m = m(m-1)(m-2)/3 * 2^(4-2m).

All remaining bases relative to the identity are at most 1/16, proving the
stated `1+c_m*4^-k+O_m(16^-k)` expansion. For m=3,4 there are no remaining
terms. This ratio is not an exact probability of nontrivial symmetry.

Finally the unnormalized coefficients have a positive Gaussian-scale lower
bound, while every P-recursive sequence is eventually bounded in magnitude
by `C^n (n!)^d` for suitable constants. This rules out D-finiteness of each
raw fixed-deficit ordinary series and each nonzero raw converse parity
series. The zero-radius observation alone would not suffice, and the proof
does not make that invalid inference. Normalized series remain rational.

## Independent finite-exact controls

`audit.py` contains two different direct models:

1. For every actual unmarked permutation, it builds affine equations on the
   unordered-edge orientation bits, including endpoint-order signs, the
   converse sign, and the path pins. A parity union-find solver counts the
   fixed tournaments. Ordinary and twisted Burnside averages agree with a
   separately implemented multiplicity-vector channel formula for all 45
   cells `1<=k<=n<=9`, involving 53,227 permutations and both signs.
2. It literally enumerates every constrained labelled tournament and obtains
   isomorphism orbits by breadth/depth search under adjacent unmarked
   transpositions, not by the producer's full-permutation orbit-list method.
   Direct reverse-complement membership decides converse-fixed orbits.
   All 28 cells `1<=k<=n<=7` agree, and every orbit is checked for pin
   preservation and orbit-stabilizer divisibility.

Both normal and `python -O` runs pass **499,356 explicit checks** and produce
identical receipts. Independent n=7 ordinary counts, in k=1,...,7 order, are

    3040, 8912, 22016, 43776, 65536, 65536, 32768.

The corresponding converse-fixed counts are

    88, 256, 256, 512, 512, 512, 512.

The producer was also rerun independently in normal and optimized modes.
Its 8,954 checks, including direct path DFS masses and all 76 ordinary
tournament classes through n=6, pass in both modes.

Reproduction, from this directory:

    python3 audit.py --receipt audit-normal.json
    python3 -O audit.py --receipt audit-optimized.json
    python3 ../session23-symmetry/check.py --receipt producer-normal.json
    python3 -O ../session23-symmetry/check.py --receipt producer-optimized.json

No assertions carry the tests, so optimization cannot disable them.

## Exact correction target and exclusions

The inherited source requiring repair is
`01-canon/theorems/THM-260-rooted-tournament-layer-decomposition.md`, not other
theorems sharing ID 260. In its General formula, replace the factor r! for
all nontrivial cycles by the product of factorials of the multiplicities of
each distinct nontrivial odd length. Equivalently its denominator is the
full centralizer

    z_lambda = m_1! * product_(l>=3 odd) l^(m_l) m_l!.

The correct deficit layer is `(n-m_1)2^E/z_lambda`. The displayed (5,3)
table entry at n=8, 128/15, is already correct; the erroneous general formula
would give 64/15. Do not double the table entry. Independent checks confirm
the complete printed rooted prefix through n=8 and total deficits
`0,0,2,4,12,40,152,784`.

The same old source's informal claim that layers form a graded ring is not
established by the Burnside formulas. The producer correctly does not inherit
it. This review approves the precise centralizer correction, not that
unspecified multiplication claim or every sentence in the legacy source.

No claim here classifies individual tournament automorphism groups, supplies
an unlabelled-scaffold quotient, proves a metagraph statement, or advances
the still-open LRC(14) theorem.
