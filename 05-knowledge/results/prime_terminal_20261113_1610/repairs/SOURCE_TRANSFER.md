# A two-arc H79 repair feeds the next actual source

**PROVED:** the marked induced source transfer and its arbitrary-module
identity; the exact response boundary of the H79 repair. **FINITE-EXACT:**
the chosen-class minimum repair and the literal chain through H721.
Chan–Pak and universal source selection remain OPEN.

## Credit the inherited ray; certify the missing physical interface

The scalar transition 5·2^k−1→15·2^k−1 and endpoint frames are already in
[the incoming frame catalogue](../../prime_ternary_lift_20261113/frames.json.gz):
the `transitions` entry has `odd_amplitude=5`, `offset=-1`, `digit=2`,
`min_exponent=1`. Its source is (2,4,9,19,7), module1, and its target is the
incoming common source (2,20,9,35,45,7), module0. No new scalar ray is claimed.

Use the equivalent source and child

    A5=(2,12,1,21,7),        module0,
    B6=(2,12,33,53,7,19),    module0.

The source isomorphism to the catalogued five-core is (1,2,0,3,4). B6 contains
A5 literally on vertices0,…,4. Its sole newborn is vertex5, which beats
0,1,4 and loses to2,3. For an arbitrary nonempty module M at0 it therefore
beats every vertex of M, with all old arcs retained.

The complete response vectors are

    H(A5[M])=9K₁+5K₂ =5P₂(M)−H(M),
    H(B6[M])=29K₁+15K₂ =15P₂(M)−H(M).

The explicit isomorphism from B6 to the incoming common source S is

    B6 vertex (0,1,2,3,4,5) -> S vertex (0,1,2,4,5,3).

It fixes the entire variable slot0. Hence it lifts to every M by the identity
on its vertices. This is an actual source-to-source transition, not an
inference from equal scalar counts. In particular, for every k≥1,

    A5[TT_k] -> S[TT_k],        5·2^k−1 -> 15·2^k−1,

is an induced one-birth construction between strong tournaments. It can be
followed by either of the incoming certified children, yielding
45·2^k−1 or45·2^k+1. On an arbitrary M the digits are2H(M) and4H(M), as in the
inherited theorem; transitivity is what makes the digits equal to2 and4.

The marked capability at A5's old vertex1 is also retained in the receipt:

    G(A5[M],1)=11K₁+5K₂=H(A5[M])+2H(M).

Thus its +4 alternative is supplied by the generic two-birth lemma in
[TERNARY.md](TERNARY.md). The target B6 is used through the incoming interface,
not through an unsupported assertion that it still has this same G gap.

The apparent later chart40·2^m−1 is just A5[TT_(m+3)]: its source has four
transitive twins. Normalizing that module exposed this inherited source
family and the exact interface match.

## The smallest repair to the chosen H79 source class

The induced +4 child of the strong H25 core from [TERNARY.md](TERNARY.md) is

    T=(2,8,11,97,15,87,151,63),       H(T)=79.

Reverse exactly the arcs {0,3} and {5,6}. The resulting tournament is

    R=(10,8,11,96,15,23,183,63),      H(R)=79.

Both are strong. Reversing only the first arc gives H97; reversing only the
second gives H61. The operation is a count-preserving **batch** rewrite, not
a path of individual neutral flips.

R is A5[TT4] under the map from old names to its canonical inflated names

    (0,1,2,3,4,5,6,7) -> (2,3,1,4,0,5,6,7).

Its transitive module, in old names and order, is (4,2,0,1), and the marked
root is old vertex3. At that root (H,E,O)=(79,40,41), so G−H=2.

All8!=40,320 vertex identifications were checked. The minimum number of arc
reversals from T **to the isomorphism class of this chosen A5[TT4] target** is
two, and the minimizing identification is unique. This is not a minimum over
all H79 graphs or all possible rearmed interfaces.

## The exact all-module failure boundary

Keep the old variable port2 of the H25-child frame, and replace it by an
arbitrary nonempty M. Neither of the two repaired arcs touches that port, so
the same fixed rewrite is defined for every M. Its complete responses are

    before: 79K₁(M)+36K₂(M),
    after:  79K₁(M)+40K₂(M).

Consequently the exact count change is4K₂(M). It is zero if and only if M
is a singleton. For TT_m it is4(2^m−2), so the repair does not rearm the entire
12·2^m+1 source family. This positive base repair and its universal obstruction
are both needed: the fixed H79 success cannot be silently promoted to a
repeatable scalar return.

## A fully named composition through the next +4 child

In the same eight old names, append vertex8 with outgoing mask159, i.e. it
beats0,1,2,3,4,7 and loses to5,6. This is the A5→B6 induced birth above.
The actual output is

    (10,8,11,96,15,279,439,63,159),    H=239.

Its marked isomorphism to S[TT4] is

    (0,1,2,3,4,5,6,7,8) -> (2,3,1,4,0,5,7,8,6).

Apply the incoming +4 constructor W(10010110), inflated at its slot4 by TT4.
It keeps all these nine vertices, reverses only their arc {5,8}, and adds
three vertices9,10,11. With those names the output is

    (3594,3592,3595,1632,3599,3607,4023,3647,3775,3072,2048,8),
    H=721=7·103.

Literal tournament validity, strongness, all counts, both graph isomorphisms,
the two-arc minimum and the single old-arc change are independently checked.
The resulting graph-level chain is

    25 -> 79 --two-arc neutral repair--> 79 -> 239 -> 721,
    orders 6,8,8,9,12.

Only the 25→79 constructor and the exact H79 rearm are new entry work in this
packet; the scalar rays and final ternary family are inherited and explicitly
matched. The concurrent onward handoff is composed separately in
[INCOMING.md](../INCOMING.md). This note does not assert that every critical
target reaches that chart, or that the H79 repair works for arbitrary module
sizes.

Reproduce with

    python -B 05-knowledge/results/prime_terminal_20261113_1610/repairs/check.py
    python -O -B 05-knowledge/results/prime_terminal_20261113_1610/repairs/check.py

The normal and optimized receipts are retained separately and agree exactly.
