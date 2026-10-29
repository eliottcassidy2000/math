# Number and operation atlas: recovered mechanisms, not matching labels

**PROVED elementary identities / CITED current repo results / OPEN transfers.**
This is a targeted recovery of old and concurrent work, not a claim to have
accounted for every occurrence of these integers in the repository.

## 1. Sources that actually change the interpretation

| Source | Mechanism retained | Boundary |
|---|---|---|
|May 3 [summand reflection](../../../07-reflections/summand-graph-fermat-zeckendorf.md), corrected Sept 17|Distinct parent pairs and the missing set {1,4,6}|The claimed Fibonacci reset was a truncation artifact|
|May 31 [operation shadows](../../../07-reflections/natural-operation-digraphs-and-product-sum-s365.md)|Two-input gates, not just one-input arrows|The full additive shadow loses the second summand; strict distinct-parent conventions also remove doubling edges|
|July 26 [THM-2422 operation fibres](../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md)|Swap-fixed diagonals; exact closure and synchronous frontier|Thin twin-center ancestry remains separately open|
|June 12 [power-balance towers](../../../07-reflections/triangular-power-balance-towers-and-additive-square-bridges.md), corrected in [HYP-2454](../../hypotheses/HYP-2454-triangular-power-balance-towers.md)|Exact ordinary/square rows and unsquared defect|The all-power unit-width bracket is refuted; do not re-prove or reuse it|
|August 1 [THM-3045 K4 lattice](../../../01-canon/theorems/THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch.md)|Unsigned edge lattice, three matching channels, index24|Unsigned edges are not the signed exterior-square representation|
|October 5 [atomic Collatz measure](../algorithmic_collatz_measure_20261005.md)|Every integer atom and its exact residue-cylinder weight|Haar almost-sure coverage does not supply atom decay|
|October 29 [Session103 grammar](../session103_words_necklaces_20261029/languages/ACTUAL_GRAMMAR_REFINEMENTS.md)|All-order unimodularity criterion, counts and native insertion ports|2065 is now proved; arbitrary tournament enumeration and arbitrary observer changes are not included|

## 2. The requested integers have several distinct roles

| Number | Recovered or explicit role | Map or reason |
|---:|---|---|
|10|First synchronous nonunique ancestry; order-ten grammar count182; size of P1(F9)|10=2+8=3+7 at the same birth stage; the other two roles have different objects|
|13|Summand target with three eventual parents; order with2065 observers; PG(2,3) point count|13=2+11=3+10=5+8; (3^3-1)/(3-1)=13|
|2065|Exact order-thirteen unimodular observer count|The cubic recurrence and separately proved grammar criterion; its factor7 is explained by WEDGE's all-order modular law|
|26|Live source26=11+15 at summand stage4; nonzero vectors of F3^3|26=27-1; antipodal quotient gives13 projective points|
|27|Summand frontier offset; ternary vector space size; a slow Collatz source|M_t-1=27*2^(t-4); the actual odd Collatz path has41 steps, independently checked|
|28|Summand stage-four frontier; P1(F27) point count|M_4=28; projective completion adds infinity to27 affine field points|

The projective map is exact: in F27, multiplication by -1 pairs its26
nonzero elements. After choosing an F3 basis these pairs are the13 points
of PG(2,3). Squaring identifies the pairs with the order13 norm-one subgroup
of F27*, because its kernel is {+1,-1}. Each projective line has four points,
and two distinct lines meet once. The checker constructs F27 as
F3[x]/(x^3-x-1), verifies a primitive element and all line incidences.

There is also an exact marked-set description:
`P1(F27) minus {0,infinity}` has26 points, and quotienting that complement
by z->-z yields13. This preserves a declared projective incidence only
after the F3-vector-space identification; it does not identify all the
28-point projective-line incidences with the13-point plane. P1(F9) has10
points by the same completion rule, but is a different field and geometry.

The earlier LRC hostile `26*{1,...,12} union {339}` has first good period27
(the witness is `2/27`; see
[THM-764, signed-pair deck](../../../01-canon/theorems/THM-764-covering-small-period-signed-pair-deck-and-q25-refutation.md)).
This is a period-existence obstruction, not a field-size consequence.

## 3. Recover the three chains before choosing another partition

Incoming Session104 recovered the matching historical source in
[its three-chain note](../session104_descent_exterior_20261029/summands/THREE_CHAINS_RECOVERED.md).
The September [square-filter result](../collatz_mod6_20260922_w6_summand_closure_square_filter.md),
S4–S9, proves the square-sum graph on 1..12 has paths

    10--6--3--1--8;    2--7--9;    11--5--4--12.

These are founded at {1,2,4}. Vertex 13 touches 3 and 12 and merges the
first and third components; vertex 14 touches 2 and 11 and completes
connectivity. Every n>=6 has an earlier neighbour because a square lies
strictly between n and 2n. Connectivity therefore persists after 14.
The strict-summand holes {1,4,6} are a different object. This recovered
source supersedes our initial ambiguity about which three-piece history
was intended; the operation-table partition below remains a separate model.

The deleted self-pair in this graph occurs exactly at vertex `2k^2`, since
twice that vertex is `(2k)^2`. This is exactly the starting term of the
even balance rows below. It is a repeating diagonal exclusion, not a
repeating leaf pattern: vertex 32 already has earlier neighbours 4 and 17.

For a precise operation-table model, adjoining label N+1 to the ordered
pair square [N]^2 adds three disjoint pieces:

    {(N+1,j):j<=N}, {(i,N+1):i<=N}, {(N+1,N+1)}.

Their sizes N,N,1 give `(N+1)^2-N^2=2N+1`. This is a separate explicit
operation-table growth rule, not the recovered three square-sum paths.
Swapping the two
parents exchanges the first two pieces and fixes the corner. After the
commutative quotient there are N new off-diagonal pairs and one diagonal.

For addition the swap-fixed equation is `2a=z`; for multiplication it is
`a^2=z`. THM-2422 already proves the resulting fibre counts. Thus doubles
and squares are genuinely the same fixed-point construction applied to two
different operations. The diagonal disappears in an alternating wedge;
it remains meaningful in the commutative operation fibre. These quotients
must not be identified without declaring which information is expendable.

## 4. The two balance towers, and what their equality forgets

Let the left side contain m+1 consecutive terms starting at a, and the
right side the following m terms. Their linear and square defects are

    D1(a,m)=a-m^2,
    D2(a,m)=(a-m(2m+1))(a+m).                         (10)

These are direct sum formulas and recover the earlier HYP-2454 identities.
The positive linear root a=m^2 partitions the shell
`[m^2,(m+1)^2-1]`; multiplying every term by two gives the requested
`2+4=6`, `8+10+12=14+16`, and subsequent rows.
The positive square root a=m(2m+1) gives

    3^2+4^2=5^2,
    10^2+11^2+12^2=13^2+14^2,
    21^2+22^2+23^2+24^2=25^2+26^2+27^2.

The square common sums are25,365,2030. In the square-balanced row the
ordinary unsquared difference is `m(m+1)`, not zero. This retained defect
is the useful coordinate; equality after a nonlinear observation does
not identify the original additive data. The higher-power bracket has
the explicit hostile `D_15(181,3)<0` at the formerly asserted upper wall;
the current repaired wider bracket and unresolved integer-root question
in HYP-2454 take precedence over its June reflection.

## 5. A literal ABC-adjacent map, with its normalization

For distinct coprime positive a,b, set c=a+b and
`h=4` if c is even, `h=1` otherwise. Then

    ((a-b)^2/h) + (4ab/h) = c^2/h                     (11)

is again a primitive positive additive triple. If c is even, both parents
are odd and dividing by four removes the full common factor. If c is odd,
the squared difference is odd and coprime to4ab. This proves primitivity
in both cases. For a=1,b=4 it gives9+16=25, hence the3,4,5 example.

Ellenberg's [Congruence ABC implies ABC](https://arxiv.org/html/math/9909098v1),
Theorem3 and Lemmas4–5, uses normalized even-power maps of this type to
force prescribed divisibility while controlling the height/radical defect.
The theorem is an equivalence between ABC and its congruence-restricted
versions, not a proof of ABC. Its proof retains prime valuations and the
power-of-two normalization; a picture of additive and multiplicative
shadows alone does not contain those bounds.

Equation (11) generally increases height. It is useful for a guarded
comparison or reduction, not automatically an integer descent. The same
distinction governs the Collatz work: a lawful re-encoding or a conserved
algebraic equation still needs an independently proved decreasing rank.

## 6. What the LRC and Collatz weight repairs have in common

The prior [positive LRC germ](../joint_port_germs_20261028/LRC.md), L1–L3,
already handles an isolated witness that ordinary length misses. If P,N
are its positive and negative active owner sets, threshold relaxation gives
it the positive first-order weight `1/max(P)+1/max(N)`. At threshold a/q
and speed height M, the certified relaxation is `epsilon=1/(4qM)` and a
nonempty good set gives relaxed length at least `1/(2qM^2)`.

This and the integer atoms in DESCENT preserve a common predicate:
an actual witness or exceptional source cannot disappear from the chosen
positive statistic. They use different constructions: a boundary derivative
for LRC, a countable atomic measure for Collatz. Neither detector establishes
its desired sign uniformly. The LRC theorem still needs positivity for every
speed row; Collatz still needs the independently proved decay that defeats
each fixed positive atom. Keeping LRC phase and Collatz original-source
residue is essential when later operations act on those labelled points.
