# What the five-part decomposition retains

**Status:** INHERITED PROVED mechanisms with source scopes below;
FINITE-EXACT independent controls; PROVED elementary distinctions.
The original observable defining the user's `1–5–4–5–1` partition has not
been identified in this bounded recovery. The explicit example in §4 is a
new reconstruction, not an attribution to previous work. No conclusion
about all tournament Hamilton counts follows from these cardinalities.

## 1. Inheritance and the closest recovered objects

The closest proved mechanism is
[THM-4733, divisor-plane fold](../../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md).
Its own status is **PROVED (pre-audit)**, with separate exact controls.
It relates a four-bit, sixteen-element cube to the twelve divisors of
`p²qr`; deletion of the two extreme divisors leaves ten. Its map preserves
a specified Klein-four action and gradings, and explicitly does **not**
preserve tournament isomorphism classes or Hamilton counts.

The canonical hostile for scalar induction is the inherited order-five
pair `1 -> 1 -> C3` and `1 -> C3 -> 1`. Both have `H=3`, automorphism
order3 and the same complete unmarked ordered path-cover profile, but
their one-vertex extension spectra differ. See
[G1–G2, the sharp order-five attachment obstruction](../../tournament_joint_parity_20261102/growth_obstruction.md).
The corrected near miss is identification of several `3–4–3` partitions:
prime/squarefree type, divisor rank, repeated-prime exponent and signed
prime difference have the same histogram but partition different objects.
The least-used relevant sidecar is the **address within a folded mirror
pair**, together with the two mirror-fixed tiling bits.

The local concept board is:

| Object | Retained data | Lost data under the proposed compression | Cheap decisive test |
|---|---|---|---|
| Four-bit cube | Sixteen states and complement | Meaning of five bin labels | Compare ordinary rank with the requested row |
| Proper divisors of `p²qr` | Exponents and prime names | One split bit in the fold | Lift both representatives and count paths |
| Blue order-five presentations | Mirror-fixed tournament matrices | Class multiplicity if only H is kept | Count all sixteen literal matrices |
| One-vertex attachment | Named cut and endpoint incidences | Joint incidence under a histogram | Reuse the two H3 parents above |
| Odd target recursion | Exact child value and reusable ports | Closure after each selected child | Keep the actual child rather than its cardinality bin |

The anchor is the odd tournament spectrum; the niche is the old divisor
fold; the wildcard is whether the requested symmetric row can itself
define a lossless tournament observer. The finite test below gives a
precise negative answer on one natural sixteen-state carrier, while
retaining the coordinate-level fold as a useful procedural generator.

## 2. Three recovered decompositions that must stay distinct

[Two decompositions of ten](../../trijection_reset_20261020/arithmetic.md)
uses proper nontrivial divisors `p^a q^b r^c`, with
`0<=a<=2`, `b,c in {0,1}`, after deleting `1` and `p²qr`.

1. Divisor rank `a+b+c=1,2,3` has counts **3,4,3**.
2. Signed height `a-b-c=-2,-1,0,1,2` has counts **1,3,2,3,1**.
3. Keeping only the sign of that height gives **4,2,4**.

Before folding, the sixteen binary words have ordinary Hamming-rank
counts **1,4,6,4,1**. The map

    (p1,p2,q,r) -> (a,b,c)=(p1+p2,q,r)

has eight singleton fibres and four double fibres. It turns sixteen
words into twelve divisor exponent vectors. Complementing all four
bits sends `(a,b,c)` to `(2-a,1-b,1-c)`; simultaneously exchanging the
two p bits and the q,r bits sends it to `(a,c,b)`. This is the exact
inherited symmetry connection, not an identification of the histograms.

In particular `1,5,4,5,1` is not the ordinary rank row of a four-bit
cube. Nor can it be the ordinary divisor-rank row of a number with
sixteen divisors: its rank-one coefficient would require five distinct
prime factors, hence at least `2^5=32` divisors. These are typing tests,
not a refutation of some other intended grading.

Related prior numerical work is routed in
[groups of five / PRIMEGAME / merged metagraph](../../groups_of_five_primegame_merged_metagraph_20261031.md).
It distinguishes the square-sum onsets `{15,16,17,23,25}`, five Fermat
primes, five Platonic solids, and five mixed order-five metagraph classes.
Its verdict is not a common five-object bijection: each proposed map
must retain a stated operation and observable. The divisor-fold theorem
above is its actual structural bridge to `p²qr`.

Two further recovered number families directly inform the induction.
[Boundary-carry returns](../../boundary_carry_returns_20261103/README.md)
proves all four two-step ternary descendants of
`h_m=8*3^m+6*2^m-5`, m>=1, retaining original vertices, clone addresses
and named pins. The outputs are `9h_m+8`, `9h_m+10`, `9h_m+14`,
`9h_m+16`. These are **two levels**, not arbitrary depth. The same packet
recovers the31/33 and32±1/13±1 comparisons, but the reusable coordinates
are the crossing matching, two-cover endpoints and separate2^m/3^m
coefficients. Its arbitrary-child substitution computes `3H+cH(S)`,
leaving the explicit extra obligation `c(H(S)-1)` for a constant offset.
That is a concrete reason why retaining only a selected path or number
is insufficient.

The [six-number atlas](../../six_roots_hidden_bits_20261103/NUMBER_ATLAS.md)
recovers `{1,23,25,63,65,67}` as the roots of the arithmetic ternary
tree after deleting7 and21. This is an actual shared provenance of
those six numbers; it does not establish that the associated tournament
interfaces recurse indefinitely. Its prior18,234-file lexical census
is inherited discovery infrastructure, not repeated here or represented
as an exhaustive semantic reading. Together these packets direct the
next step toward an invariant family of **actual return interfaces**,
rather than another matching of small cardinalities.

## 3. A literal witness of the fold's Hamilton-count loss

Use the fixed path `0 -> 1 -> 2 -> 3 -> 4`. A bit1 reverses its specified
nonspine arc relative to the transitive order. Let the moved coordinates
be

    (p1,p2,q,r)=((0,2),(2,4),(0,3),(1,4)),

and keep fixed coordinates `(0,4)` and `(1,3)` both zero. The two words

    (1,0,0,1)  and  (0,1,0,1)

have the same folded address `(1,0,1)`, and the same ordinary rank2.
Their tournament adjacency row masks are respectively

    (26,12,25,16,2), H=9;
    (30,12, 8,16,6), H=5.

Here bit j of row i indicates the arc `i -> j`; both are genuine loopless
tournaments containing the fixed path. Exhausting all four settings of
the fixed coordinates and all four choices of q,r gives sixteen double
fibre comparisons, eight of which change H. With both fixed bits1, the
same moved-word comparison gives H15 versus H9.

This strengthens the old note's general loss statement with tiny literal
certificates; it does not claim a new fold theorem. The mechanism is
that replacing one named reversed arc by its mirror changes which
other reversed arcs can coexist in a Hamilton order. The count
`p1+p2` forgets the incidence. To restore this particular map, one must
retain which p bit is1 on every `a=1` fibre and retain both fixed bits.
That restores the whole labelled six-bit tournament, so every subsequent
attachment is defined. Keeping only H after the lift would still lose
the next-cut response, by G1 above.

## 4. Why the histogram alone does not identify its grading

For a nonuniqueness control only, here is an arbitrary way to obtain the
requested row. It is **an example in this note**, not a recovered
historical definition or a proposed induction mechanism.
In the ordinary four-bit rank distribution, select one complementary
weight-two pair, say

    a=1100,  complement(a)=0011.

Move `a` to grade1 and its complement to grade3; leave every other word
at its Hamming weight. The five grade sizes become

    1,5,4,5,1.

Complement still acts by `grade -> 4-grade`. There are six choices of
the ordered middle word a. Thus the requested histogram is compatible
with a sixteen-state binary carrier and complement symmetry. Those
facts alone do not specify the operation that moved the pair, the
meaning of grade, or any Hamilton-count transport. The grading is not
the original rank and is not asserted to preserve cube adjacency.

There is an immediate obstruction to using *any* five-bin partition as
a lossless scalar Hamilton observer on the **blue order-five cube**.
The sixteen reflection-fixed marked-path tournaments have exact H
histogram

| H | 1 | 3 | 9 | 11 | 13 | 15 |
|---:|---:|---:|---:|---:|---:|---:|
| Presentations | 1 | 1 | 4 | 3 | 3 | 4 |

There are six distinct H values. If H were a function of five bin
labels, it could take at most five values. Therefore no five-bin
partition of those sixteen states, including one with sizes
`1,5,4,5,1`, determines H on its own. This is a proved pigeonhole
consequence of the exact census, not a statement about every possible
sixteen-element object. These are presentation counts, not sixteen
Hamilton paths and not sixteen unlabelled tournament classes.

The positive survivor is to keep the **bin plus its internal address**.
The bin can organize procedural generation; it does not replace the
named matrix/ports needed for a constructive return. Even a bin and H
together need not determine future attachment. The next useful question
is which internal addresses can be merged while preserving *all* allowed
return responses, rather than whether five bin sizes sum to sixteen.

## 5. Connection contract and verification

The recovered source is the named four-bit cube with its two mirror-fixed
bits; target is the divisor exponent box. The fold above preserves the
Klein action, rank and specified sign grading. It destroys the p-bit
address and can change H. The missing bit plus fixed bits restores the
source. The cheapest decisive test is the displayed H9/H5 pair.

For an induction step the source is instead a named tournament and an
attachment word; the target is its child tournament. The preserved
predicate must be the actual path count `3H+2` or `3H+4` **and** the
availability of the next allowed step. The old fold and numerical
decompositions establish neither predicate. This is precisely why the
order-five hostile matters before attempting to discard joint ports.

The independent [checker](check.py) uses no historical producer imports.
It compares literal permutation counting with endpoint subset DP on
all64 fixed-path order-five tournaments, extracts all16 blue states,
checks every double fold fibre at all four fixed-bit settings, and
checks every one of the six complement-symmetric regradings. There is
no random sampling and no inference from omitted larger orders.

    python3 -B 05-knowledge/results/decomposition16_odd_spectrum_20261110/history/check.py
    python3 -B -O 05-knowledge/results/decomposition16_odd_spectrum_20261110/history/check.py

The two stored outputs are byte-identical. They pin seven source hashes
and the checker hash. Exact source-pattern searches covered current
navigation, canon, result notes, reflections and archived inbox text;
five recent turns in the specifically relevant tiling chat were also
inspected. This was bounded recovery, not a claim that an unlocated
prior definition never existed. No mathematical result relies on that
negative search outcome.
