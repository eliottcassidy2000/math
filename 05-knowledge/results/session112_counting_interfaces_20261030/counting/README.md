# From divisor-prefix search to canonical tournament ranking

**PROVED:** the finite-set identities, weighted-observer boundary and reduction
below. **FINITE-EXACT:** all cut thresholds on the complete tournament cubes
of orders 3, 4 and 5, with independent adjacency-action controls; affine-rank
checks on selected cuts through order 4. **OPEN:** a compact, fast evaluator for
the resulting saturated cuts. No new A000568 value, factoring algorithm,
subexponential graph algorithm or independent Lean certification is claimed.

## Inheritance and the new question

The inspiration is Rudich and Rousseau,
[*Deterministic Factoring in Time N^{o(1)}*](https://hexagonmath.org/pdf/2610.00600v1),
Section 4, step 1: recover the least divisor by deciding successive binary
prefixes. Only that elementary search reduction is used here. The paper's
claimed counting engine and formal complexity theorem are not dependencies.

The closest repository mechanism is
[`anchor_content_flat_quotients_20261030/TOURNAMENTS.md`](../../anchor_content_flat_quotients_20261030/TOURNAMENTS.md),
at `b51ef5cd739375e6475ceafd42e36e86858c70ca`: forget construction bases by taking
an exact union, then apply Burnside to affine-flat intersections. Its hostile
is counting with the normalizer of a first bank after moves leave that bank.
Its unresolved sidecar is overlap between differently marked copies.

The question here is not another unlabelled total. Can that union-counting
interface **select and rank canonical representatives**, without canonicalizing
every labelled tournament first? The answer is an exact oracle reduction,
with a necessary change from prefix cylinders to saturated cumulative cuts.
The finite checks below certify the reduction, not an efficient oracle.

Concept board: prefix existence / orbit saturation / overlap / witness
multiplicity / canonical order / counting cost. The central distinction is
that positivity suffices to select a first witness, whereas ranking requires
unit mass for each distinct target object.

## 1. Why ordinary prefix recursion breaks on an orbit quotient

Let a finite group G act on an invariant finite subset X of `{0,1}^d`.
Use the numerical order on its binary strings, reading the most significant
bit first. A prefix u defines a cylinder C_u.

The quantity

    O(u) = number of G-orbits meeting X ∩ C_u

detects whether that prefix is possible, but generally

    O(u) ≠ O(u0) + O(u1).

**Small hostile.** Let G exchange the two coordinates and let X={01,10}.
There is one orbit. It meets both prefix children, so their counts are 1 and 1,
while the root count is 1. Greedy existence still returns 01; ordinary additive
rank/unrank incorrectly allocates two target objects.

Replacing G with a prefix stabilizer does not repair the quotient. Let C3
rotate the coordinates of X={001,010,100}. Its prefix-0 part is {001,010}.
The setwise stabilizer of that cylinder inside C3 is trivial, so it gives two
classes there, although those two points belong to one full G-orbit.
A stabilizer-chain construction needs exact transporter/overlap information;
the subgroup alone does not supply it.

## 2. Cumulative orbit-minimum theorem

Write `m(O)=min O`. For an integer t in `{-1,0,...,2^d-1}`, define

    A_t = {x ∈ X : x ≤ t},
    S_t = G A_t = union_{g∈G} g A_t,
    F(t) = |S_t/G|,       F(-1)=0.

Then

    F(t) = #{orbits O in X/G : m(O) ≤ t}.                         (1)

**Proof.** S_t is a union of complete orbits. An orbit occurs in this union
exactly when one of its elements lies at most t, exactly when its least
element does. Each such orbit contributes one. This also proves

    F(t)-F(t-1) ∈ {0,1},

with difference 1 precisely at canonical orbit minima.

For a binary prefix u of length l, let `[a_u,b_u]` be its full numeric interval.
The number of **canonical minima** with that prefix is exactly

    K(u) = F(b_u)-F(a_u-1).                                     (2)

These counts are additive under splitting u into u0 and u1. Unlike O(u),
they assign each orbit to the unique child containing its minimum.

Given zero-based rank r with `0 ≤ r < F(2^d-1)`, its canonical representative is

    min {t : F(t) > r}.                                         (3)

Binary search uses at most d calls to F. Its correctness follows from the
monotonicity and unit jumps in (1), or recursively from (2). Ranking an already
canonical string x uses `F(x-1)`. This does **not** give the rank of an arbitrary
labelled x without finding its orbit minimum or equivalent information.

The minimum of all nonempty X needs only ordinary prefix existence, as in the
factoring reduction. Saturated cumulative counts become essential when the
requested output is a specified rank among distinct orbit representatives.

## 3. Affine-flat implementation of the oracle

Suppose X is an invariant union of s affine flats in `F_2^d`, and every element
of G acts affinely. A cut `{x≤t}` is a disjoint union of at most d+1 coordinate
prefix flats: at each 1-bit of t, copy all earlier bits and choose 0 at that
bit; include the singleton t. Intersect these flats with the s input flats,
then apply all G elements. Their union is exactly S_t.

Merge intersections by exact inclusion-exclusion, obtaining an integer
indicator expansion

    1_{S_t} = sum_H c_H 1_H.

For g acting as `x ↦ L_g x+b_g`, its fixed set is the affine system
`(L_g-I)x=b_g`, possibly inconsistent. Thus

    F(t) = (1/|G|) sum_{g∈G} sum_H c_H
                      |H ∩ Fix(g)|,                            (4)

and each consistent intersection contributes `2^(d-rank)`.

**Proof.** The prefix partition is exact by the first differing bit. Affine
images and intersections are affine, so the constructed union is S_t.
Inclusion-exclusion is an equality of indicators, including every overlap.
Burnside applies because S_t is invariant; substitute the indicator expansion
into its fixed-point count. Equation (1) gives the claimed oracle semantics.

For tournaments on n labelled vertices, set `d=binom(n,2)` and record an edge
bit as 1 when its smaller-labelled endpoint points to its larger-labelled
endpoint. Relabeling is a coordinate permutation **plus endpoint-order sign
flips**. The proof retains these flips. Arbitrary converse, unsigned bank
identification, original-bank normalization and full vertex relabeling are
not interchangeable quotients. Named vertices use the intended pointwise
stabilizer and must remain part of X's type.

The construction takes at most `s(d+1)|G|` initial flats. Its intersection
arrangement can still be enormous; the straightforward formula has up to
`2^{s(d+1)|G|}` subset terms before merging, with potentially large integer
coefficients. Exact coalescing and containment pruning help but do not prove
a polynomial bound. Fixed-point affine ranks alone do not control the number
of intersections.

**Conditional complexity statement.** If an implementation of F has uniform
bit cost T(d), (3) costs `O(d T(d)+d^2)` with d-bit thresholds and ranks.
This is a reduction, not a bound on T. A hypothetical `2^{o(d)}` evaluator
would give `2^{o(d)}` selection, but for tournaments d is quadratic in n;
this is not a `2^{o(n)}` statement. No complexity claim follows merely from
replacing graph isomorphism calls by an exact sum.

## 4. Positive witness fibres versus exact counts

The paper's network construction has uniquely determined auxiliary values
once a factor pair is chosen. Unique extension is sufficient for exact
witness counts, but it is stronger than needed for **existence search**.

Let W be finite and ordered. Suppose the observer counts auxiliary records
above each witness w with integer multiplicity `μ(w)≥1`. For any prefix u,

    Q(u) = sum_{w∈W∩C_u} μ(w)

is zero iff that prefix has no witness. Consequently the bit-by-bit minimum
search returns `min W` for arbitrary positive multiplicities. If the total
multiplicity is bounded by B, any exact residue-based zero certificate for
integers in `[0,B]` can replace the full Q value. The bound must cover auxiliary
multiplicity, not merely the number of distinct witnesses.

Variable multiplicity invalidates distinct-witness rank/unrank. For W={01,10}
with multiplicities 2 and 1, minimum search returns 01 correctly, but weighted
rank 1 returns 01 again; distinct-witness rank 1 is 10. Known common
multiplicity can be divided out. Signed weights are worse: weights +1 and -1
can give zero for a nonempty set. Thus positivity is a mathematical resource,
not a cosmetic implementation detail.

The exact group version is useful for diagnosing a mistaken Burnside call.
For any rational point weight w(x),

    (1/|G|) sum_g sum_{x∈Fix(g)} w(x)
       = sum_{orbits O} (1/|O|) sum_{x∈O} w(x).                 (5)

**Proof.** Swap the finite sums. Each x contributes once per element of its
stabilizer, whose size is `|G|/|O|` by orbit-stabilizer.

Therefore:

- A saturated set indicator has mass 1 on each admitted orbit, giving (4).
- An unsaturated prefix indicator gives the fraction `|O∩C_u|/|O|` per orbit.
  For the C3 hostile it gives 2/3, not one orbit.
- Invariant positive integer multiplicity gives positive integer **weighted**
  orbit counts. It preserves zero tests, not distinct-orbit ranks.
- Noninvariant positive multiplicity can give fractions: weights 1,1,2 on the
  three points of the C3 orbit give 4/3.
- Merely summing indicators of overlapping construction flats counts records,
  not their union. This may preserve first-witness existence but not (1)–(3).

Equations (1), (4) and (5) explain precisely which information must be retained
or corrected: saturation removes partial-orbit mass, and union indicators
remove duplicate-construction mass. Neither operation is a free compression.

## 5. Exact checks and limits

Run `python3 check.py` and `python3 -O check.py`. They produce identical JSON.
The checks do not use Python assertions. Every relabeling of every tournament
at n=3,4,5 is compared with a separately constructed adjacency-matrix action.
For **all** 8, 64 and 1,024 numeric cuts, the script constructs saturated sets
literally and checks Burnside against independently enumerated orbit minima,
unit jumps and complete unranking. It recovers 2, 4 and 12 ordinary classes.

For n≤4, selected thresholds also use an independent affine-equation/RREF,
union-indicator and fixed-rank implementation of (4); their complete labelled
indicators are checked pointwise. The naive implementation merges identical reduced intersections but intentionally
lacks containment pruning and minimal-cover simplification: one four-vertex
full-cube cut leaves 728 signed terms although a
single full-cube flat describes the same set. This is a representation-cost
warning, not a lower bound or a claim that simplification is impossible.

The five-vertex canonical numeric representatives are
`0,2,4,5,8,9,10,11,12,40,41,76` in this specified edge-bit convention.
They are a finite control, not new enumeration data. The worthwhile open
continuation is a compact saturated-cut carrier or transporter theorem for
families already represented compactly. Without that, the factoring paper's
counting-to-search pattern supplies exact semantics but no imported speedup.
