# Full slot equivalence still loses physical line wiring

**FINITE-EXACT:** every converse-merged tournament class at orders3,4,5;
the first failure is order5. The two full-Q collision pairs are inherited
from the incoming census and independently recounted here. The added
connection is to the literal complement-line dart matrix and its exact
transition probabilities. No larger census is independently audited here.

## Convention and incoming source

For a tournament T, let F_j(T) count **ordered** covers of its vertices by
j nonempty directed paths. Define

```text
Q_T(t)=sum_j F_j(T) binom(t,j).
```

Thus Q counts assignments to t labelled path slots, permitting empty slots;
Q_T(1)=H(T). It is invariant under converse. This is the convention of
[Session310, NEUTRAL_TAILS](../session310_positive_clocks_20261113/NEUTRAL_TAILS.md).
The incoming [THM-5153, exact slot pairings](../../../01-canon/theorems/THM-5153-exact-slot-pairings-of-merged-tournament-classes-and-self-converse-structure-to-order-nine.md)
and [S28 packet](../kps_s28_exact_pairing_sc_census_20261113/README.md) already
record eight Q-classes among the ten merged classes at order5, including
one SC–SC collision and one SC–NSC collision. Those facts are not claimed
as new. The incoming revision is
`5479c530c7907b5442dd885cbd6b47d76b41b8ba`; exact source-byte SHA256s are in
[the normal receipt](audit/slot_profiles.normal.json).

**INHERITED PROVED:** equality of the full Q polynomial allows substitution
inside every uniform exterior context for the scalar Hamilton-path
readout. Each exterior reduced itinerary with j visits has F_j fillings.
The theorem does not preserve selective internal ports or individual path
identities. The present result addresses precisely that excluded interface.

## Complete order-five profiles

Codes use the root census convention: bit1 means i→j at the corresponding
lexicographically ordered pair i<j; each merged code is the minimum under
relabeling and converse. SC means self-converse. Degree counts literal
darts, including doubled loop entries. Metadata and full incidence rows
are retained unchanged from [the root census](audit/small.normal.json).

| code | H | Aut order | SC | dart degree | Q_T(t) |
|---:|---:|---:|:---:|---:|---|
|0|1|1|yes|1|t^5|
|2|3|3|no|2|t^5+2t^3|
|4|5|1|no|10|t^5+4t^3|
|5|3|3|yes|1|t^5+2t^3|
|8|9|1|yes|9|t^5+6t^3+2t|
|10|9|1|yes|9|t^5+6t^3+2t|
|12|11|1|yes|11|t^5+8t^3+2t|
|40|15|3|yes|5|t^5+8t^3+6t|
|41|13|1|yes|13|t^5+8t^3+4t|
|76|15|5|yes|3|t^5+10t^3+4t|

The only nontrivial full-Q blocks are therefore `{2,5}` and `{8,10}`.
The first has one self-converse node and one node representing a converse
pair. The second has two self-converse nodes. In particular the two H15
nodes do **not** have equal full Q, despite their equal H.

For codes8 and10, both independent methods give the complete ordered
nonempty-cover vector

```text
(F_1,F_2,F_3,F_4,F_5)=(9,66,186,240,120).
```

They also have identical score multiset `(1,1,2,3,3)`, H=9, automorphism
order1, self-converse status, and dart degree9. Their physical incidence
to the transitive node0 differs:

```text
A_(8,0)=0,        A_(10,0)=1.
```

Consequently even the combined observer
`(Q, scores, H, |Aut|, SC status, dart degree)` does not determine a
class's complement-line response. This conclusion uses two exact
polynomials, not evaluations at selected slot counts.

## A full-Q quotient has the same two-step loss

Let A be the symmetric literal merged dart matrix and let
`P_(u,v)=A_(u,v)/degree(u)` be its random-walk matrix. The transitive node0
has its sole dart to node10. Thus its actual two-step return is

```text
(P^2)_(0,0)=1/9.
```

The full-Q block B={8,10} has total degree18 and exactly one dart to0.
The degree-weighted stationary quotient therefore has transitions
`0→B` with probability1 and `B→0` with probability1/18. Its predicted
return is1/18, losing exactly1/18. Equivalently, the equal-Q states8 and10
have different one-step transition probabilities to the singleton Q-block
of node0, respectively0 and1/9. The full-Q partition is not strongly
lumpable.

This does not prohibit an arbitrary history-dependent or enlarged-state
model. It disproves a memoryless transition rule determined only by full Q
that works for every starting class, and explicitly demonstrates the loss
in the degree-weighted stationary quotient.

The boundary is exact at small order. At orders3 and4, Q separates all
merged classes: there are respectively2 and3 Q-classes, with no collision.
At orders0,1,2 there is only one ordinary/merged class. Therefore order5
is the first order where full-Q equivalence can lose merged-class wiring;
the displayed pair proves that it already does. This boundary says nothing
about ordinary classes before converse merging, whose converse partners
necessarily have equal Q.

## Why uniform-context equivalence does not descend to wiring

The two maps preserve different information. Uniform substitution keeps
only how many fillings exist for each run count. In the converse-merged
gauge, a complement-line target can be computed by choosing an actual
Hamilton path and reversing its successive arcs: taking the full converse
then gives the native nonpath-complement target. Its output depends on
**which path and which internal arcs**, followed by ordinary-class and
converse identification. Equality of all F_j does not
transport those choices or their resulting class multiplicities.

Source: an actual tournament with its ordered path covers. Target1: the
full slot polynomial, preserving every uniform-exterior scalar count.
Target2: the physical complement-line multigraph, preserving actual dart
partners. The projection to Target1 loses the path-to-partner incidence
needed for Target2. The cheapest hostile is codes8 and10 above: all cover
counts agree, yet only one has a dart to the transitive class. The missing
sidecar is literal path/arc/partner incidence, not another unmarked slot
count. This extends the packet's score-plus-H obstruction to full Q while
respecting the incoming uniform-context theorem.

## Independent replay and scope

[audit/slot_profiles.py](audit/slot_profiles.py) imports no producer code.
For every representative in the root's complete n=3,4,5 merged census it
computes all F_j twice:

1. Enumerate every vertex permutation and every subset of its n−1 break
   positions; accept when every unbroken adjacency is directed forward.
   Concatenation and the breaks uniquely specify an ordered path cover.
2. Compute Hamilton-path counts on every vertex subset by endpoint dynamic
   programming, then convolve over ordered nonempty subset partitions.
   Each first subset contributes its path count and the remaining slots
   are counted recursively.

The methods agree on all fifteen merged representatives. Exact rational
conversion checks the entire binomial-to-monomial polynomial, and the
root H values provide a further first-coordinate control. The checker
preserves all root automorphism, SC, degree and neighbor metadata, groups
full polynomial coefficients, and verifies the exact return probabilities.
Normal and optimized receipts agree byte for byte:

```sh
python3 -I -B 05-knowledge/results/joint_metagraph_dyadic_20261113/audit/slot_profiles.py
python3 -I -B -O 05-knowledge/results/joint_metagraph_dyadic_20261113/audit/slot_profiles.py
```

S28's orders6–9 remain **INHERITED FINITE-EXACT**, not independently audited
by this check. This result supplies neither Chan–Pak completeness nor an
all-order classification of Q-classes or diagram symmetries.
