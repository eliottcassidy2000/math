# Intrinsic Collatz swap order: a transitive relation with a nonlocal symmetry loss

**Status: PROVED** elementary statements S1–S5 below, with separate
**FINITE-EXACT** controls. No convergence claim, new-literature priority claim,
or equivalence between Collatz systems and tournament categories is made.
Source pin: `e871dd4513c7dcaa75cb29bb9847aafac35a795c`.

## Useful result first

There is a genuinely intrinsic pairwise relation on Collatz affine blocks:
compare the outputs obtained by exchanging their chronological order. It is
**always a total preorder**, never a cyclic tournament. Its scalar coordinate
is the negative reciprocal of the block's rational cycle point. This gives an
exact sorting theorem but also a sharp warning: equal comparison profiles
permit false nonlocal symmetries. Two commuting blocks can be interchanged
adjacently, yet not on opposite sides of a third block.

The denominator-101 examples have different swap coordinates `-101/19` and
`-101/31`. Their distinction is thus real before appending anything; reduced
denominator and aggregate clocks erased it. Appending `(1)` exposes it as
`35` versus `175`. A separate hostile example shows that even retaining the
whole unweighted swap relation does not repair arbitrary ordered assembly.

## Inheritance, concept board, and scope

Closest mechanisms:

* `05-knowledge/results/collatz_fermat_reframes_20261029/README.md`, C1–C3:
  marked cycle point, ordered carry composition, and the `101 -> 35/175` pair.
  That file is in the pinned Git tree even if omitted by sparse checkout.
* `05-knowledge/results/session100_selective_collatz_targets_20261029/atlas/AFFINE_ATLAS_FINGERPRINT.md`, §§1–4:
  pairwise order-difference blindness and the three-row multiplier ratio.
* `05-knowledge/results/session87_four_vertices_enumeration_20261028/EXACT_GENERATION.md`, §§6–7:
  identical cyclic triples and all subset-H values but automorphism orders
  `81/27`; corresponding binary attachment orbits `176/512`.
* `05-knowledge/results/session92_joint_order_carriers_20261029/runs/RECOVERED_RUN_WITNESS.md`:
  separate per-cursor transporters do not supply one common transporter;
  the even/odd interior run order affects the native partner Hamilton count.

The corrected near miss is to infer global substitution compatibility from
local observer equality. The least-used coordinate here is **scale inside a
commuting fixed-point family**. Portfolio: anchor = marked-cycle continuation;
niche = intrinsic swap relation; wildcard = sorting and fictitious symmetry.
The live board is carry / fixed point / swap rank / aggregate clocks / genuine
tournament transporters. The theorem tests all five rather than equating the
numerals in different lanes.

All comparisons are comparisons of rational affine charts. Distinct words need
not be legal trajectories from the same ordinary integer source. Actual source
cylinders and internal valuations remain separate guards. Every nonempty
positive valuation word does define a rational odd cycle, as inherited in C1.
“Same clocks” for the focal examples means `(length,total valuation)=(3,7)`,
not equality of the ordered valuation words `(1,1,5)` and `(1,3,3)`.

## S1. The canonical swap preorder and its exact tournament boundary

Let `A(x)=u x+v`, `B(x)=s x+t`, with `u,s,v,t>0`. Write `AB` for
chronological execution of A then B, hence `AB=B∘A`. Set

    r(A)=(u-1)/v.

Then

    AB(x)-BA(x)=v(s-1)-t(u-1)=v t [r(B)-r(A)].       (1)

The input cancels. Thus `A≺B` iff `AB>BA` is exactly `r(A)<r(B)`.
The vertices are labelled blocks, the pairwise observable is (1), and the
orientation gauge explicitly chooses the larger chronological output.
Ties are retained as ties, never arbitrarily oriented.

Consequences:

1. The relation is a total preorder. It has no strict directed triangle.
2. Quotienting by equal r gives a transitive tournament (possibly infinite).
3. `r(A)=r(B)` iff A and B commute exactly.
4. If `u!=1`, A has fixed point `x_A=v/(1-u)`, so `r(A)=-1/x_A`.
   The translation case `u=1` has r=0 and fixed point at infinity; it causes
   no exception to (1). Native nonempty odd-multiplier words have `u!=1`.

For Collatz multiplier b and word w of length k and height H,

    u=b^k/2^H, v=C_b(w)/2^H,
    r(w)=(b^k-2^H)/C_b(w).

No nontransitive tournament, in particular no C3 quotient or its `81/27`
module-symmetry mechanism, can be faithfully represented by this strict
relation while preserving its directed triples. This is a proved obstruction
to a proposed bridge, not the assertion that all possible arithmetic
observers are transitive.

## S2. Composition has an interval law; sorting is exact

The rank of chronological AB is

    r(AB)=[s v r(A)+t r(B)]/[s v+t].                  (2)

Both weights are positive; therefore it lies between the two ranks, strictly
between if they differ. This is the affine form of the older weighted-mediant
fixed-point phenomenon, not a new claim that fixed points average unweighted.

For any finite list of these blocks, the common slope product does not depend
on order. Interchanging adjacent A,B changes the final intercept by (1)
times the product of slopes of the suffix. This multiplier is positive;
any prefix cancels because AB and BA have the same slope. Consequently:

* increasing-r order maximizes the final output for every real input;
* decreasing-r order minimizes it;
* every maximizing/minimizing order is sorted up to tied adjacent blocks.

Proof: bubble-sort any inversion; its improvement is strict exactly when
the ranks differ. An unsorted sequence has an adjacent inversion. This gives
the stated equality boundary, not merely a sufficient schedule.

This optimizes unrestricted chart permutations. A native source-admissibility
constraint may remove the sorted schedule; no guarded Collatz optimization
claim follows without checking its cylinder.

## S3. All elementary pairwise swap data forget the multiplier

For an elementary branch `f_(b,a)(x)=(b x+1)/2^a`,

    r(f_(b,a))=b-2^a,
    f_(b,c)∘f_(b,a)-f_(b,a)∘f_(b,c)
        =(2^a-2^c)/2^(a+c).                          (3)

Hence even the **entire weighted elementary pairwise swap matrix**, on any
fixed valuation alphabet, is independent of odd multiplier b. Its strict
relation simply orders the valuation letters in reverse numerical order.
Retaining actual rank values would restore b; the pairwise differences have
forgotten their common additive rank origin.

Composition breaks this blindness: a common suffix `f_(b,e)` multiplies a
nonzero difference (3) by `b/2^e`. Alternatively the inherited same-clock
three-row atlas `112,121,211` gives ratio `b/2` without cross-length scaling.
The loss is therefore operation-specific, not a shortage of pair samples.
This refines the Session100 two-row lesson to an arbitrary elementary alphabet.

## S4. The focal denominator split is a loss of marked rank

At multiplier3:

| word | slope u | intercept v | rank r | cycle x | cycle after appending1 |
|---|---|---|---|---|---|
|115|27/128|19/128|-101/19|19/101|37/35|
|133|27/128|31/128|-101/31|31/101|221/175|

The appended carry is `3C+128`; the unreduced denominator is `256-81=175`.
The two gcds are `gcd(185,175)=5` and `gcd(221,175)=1`. This accounts for
35 versus175 exactly. Their real ordering already differs in the rank
coordinate; the arithmetic future additionally tests divisibility of the carry.

The marked cycle point together with the aggregate clocks is a complete
rational affine-block carrier: clocks determine `u`, and `v=(1-u)x`.
Equivalently clocks and nonzero r determine `v=(u-1)/r`. Thus one need not
store a redundant carry and marked point simultaneously. The reduced
**denominator alone** discards the numerator of the marked point and fails.
For native fixed-b positive valuation words, exact `(k,H,C)` also determines
the full ordered word: if k>1 then
`a0=v2(C-b^(k-1))`, and the remaining carry is
`(C-b^(k-1))/2^a0`; recurse with length k-1 and height H-a0.
The tail carry is odd, making this valuation exact; at length1 the
remaining letter is H. Thus exact clocks plus rank or marked point also
encode internal phases in this native word class, even if decoding them
requires work. This word-recovery claim does not apply to arbitrary
positive affine blocks. The parallel continuation lane proves and tests
this native injectivity statement independently.

## S5. Tied comparison symmetry is not nonlocal assembly symmetry

A sharper hostile tests a tempting repair: retain the entire intrinsic
comparison relation on blocks. Take

    A=f_(3,2)=(3/4,1/4),
    A²=(9/16,7/16),
    B=f_(3,1)=(3/2,1/2).

Here pairs denote `(slope,intercept)`. A and A² have the same fixed point1,
rank -1, and identical comparison profiles against every positive affine
block. They commute. The weak comparison graph permits exchanging their
labels while fixing B, whose rank is1.

But chronological products are

    A B A² = (81/128,119/128),   word (2,1,2,2),
    A² B A = (81/128,143/128),   word (2,2,1,2).

They have identical aggregate clocks `(4,7)` but different marked cycle
points `119/47` and `143/47`. **These two words are cyclic rotations, so
they represent the same unmarked rational cycle.** This is deliberately a
marked-composition obstruction, not a claim of different closed cycles.
Appending the same letter `(4)` makes the missing mark visible even to the
closed denominator observer: the resulting points are `97/361` and
`557/1805`, since the unreduced carries are485 and557 and the common
unreduced denominator is1805. Their gcds are5 and1. Thus cyclic re-marking
is a valid symmetry of a closed cycle but is not compatible with appending
a letter at an unspecified marked end. In general the intercept difference is

    V_(A B A²)-V_(A² B A)
       =u[(u-1)t-(s-1)v]=u v t[r(A)-r(B)].            (4)

Thus adjacent commutation is lawful, while transporting one tied block past
a third block is not a symmetry unless the third block also commutes.
The missing coordinate is scale `u`, not merely a strict/tied comparison bit.

Rank alone is not even a right-continuation congruence: A and A² both have
cycle1, yet after the same B their cycles are `-7` and `37/5`. The denominator
and positivity can both change. Adding aggregate clocks restores the affine
map in this native class, as S4 proves.

For a finite labelled bank with rank multiplicities m_i, the automorphism
group of the **unweighted weak comparison relation** is `product_i S_(m_i)`:
its finite ordered rank classes cannot move, and labels within a class may
be permuted. S5 shows why these observer symmetries must not automatically
be treated as symmetries of actual affine assembly. This is closely parallel
to the tournament warnings about preserving actual child types and common
transporters, but it is not an identification of the two automorphism groups.

## Bridge ledger and stopping boundary

| source → observer | preserved | lost | sufficient sidecar / decisive hostile |
|---|---|---|---|
| marked affine blocks → swap relation | every swap sign, adjacent commutation, optimal ordering | slope/intercept magnitudes and nonlocal assembly | actual `(u,v)`; S5 |
| fixed-clock cycle → reduced denominator | integrality after closed-cycle odd forcing divisibility | marked numerator and future appended gcd | marked rational point; 101→35/175 |
| elementary charts → all pairwise swap weights | complete elementary swap effects | multiplier | three-row atlas or common-suffix response; (3) |
| tournament child/ports → local profiles | specified local rooted/H observations | common global transport and ordinary child type | named common transporters; Session87/92 hostiles |

The shared theorem pattern is that an observer quotient need not be a
congruence for the native operation. Here it is proved by exact maps and
counterexamples in each source category, not by pretending the categories
are isomorphic. No claimed bridge transports Hamilton counts, automorphism
orders, Collatz convergence, or primality from one side to the other.

## Reproduction

Run from this folder:

    python3 check_swap_observer.py > receipt.json
    python3 -O check_swap_observer.py > /tmp/swap-optimized.json
    cmp receipt.json /tmp/swap-optimized.json

Both modes give identical receipts, 127,068 explicit checks. Universe:
odd multipliers3,5,7,9; every nonempty word of length≤3 over digits1..4 for
28,224 ordered block pairs; every triple of length1/2 words over digits1..3
for6,912 triple banks, all six permutations each. Controls include independent
closed carry sums versus recursive composition, rank interval identities,
commutation iff ties, both sorting extrema, and exact focal/S5 rationals.
No assertion is disabled by optimization. The unbounded theorems follow from
(1)–(4) and the proofs above, not this finite census.
