# Paired binary rewrites, actual source changes, and the remaining target equation

**11 November 2026. PROVED IN NOTES / FINITE-EXACT as specified below.
Chan–Pak remains OPEN.** This continuation produces actual paired
constructions and all-length identities. It does not complete the proof
that7 and21 are the only missing positive odd Hamiltonian-path counts.
The existing all-order exclusions of7 and21 are inherited from
`05-knowledge/results/exclusion_7_21_20261103/README.md`; this session
concentrates on existence and does not repeat that audit.

The main advance is a change in what can be manipulated. We now have a
common binary defect for **both** H13 targets, a digit operation that
cancels it for arbitrary prepared prefixes, and actual source switches
whose endpoint data remain explicit. We also have an exact comparison
under every future O/Z/J construction. These are usable operations rather
than a larger numerical coverage range.

## 1. Simultaneous rebasing through an explicit binary carry

Use the inherited function

    g(0)=g(1)=1,
    g(2q+1)=g(q),
    g(2q)=g(q)+q+(q mod2),

and the golden H13 pair

    L(B)=13g(B)+32B+6(B mod2)+18,
    L(B)+2.

Two actual C3-based constructions have counts

    M(B)=3g(B)+36B+30,
    M(B)+2.

Both differences, H13 minus C3, are the **same** function

    E(B)=10g(B)−4B+6(B mod2)−12.                   (1)

This gives a simultaneous construction problem with one retained integer
defect. Append the binary block1001, i.e. hexadecimal9. For

    B_m(q)=16^m q+9(16^m−1)/15,

the exact identity is

    E(B_m(q))=10g(q)−4q+4m−6.                    (2)

Given **any** prefix p≥0, choose

    q=4p+3,
    m=q−(5g(q)−3)/2 ≥2.

Equation(2) vanishes. Therefore both resulting H13 counts have literal
C3-source witnesses. The original binary prefix survives, followed by11
and a computed number of1001 blocks. The first prepared example is
B=921, with counts34293/34295.

The operation-bearing appearance of five is exact: the nonconstant
coefficient returns through the doubling orbit
`2/5→4/5→3/5→1/5→2/5`. Its four-step word is1001, and its constant
drift is4. This is a concrete binary-carry mechanism; it does not assert
that all the historical pentads are equivalent mathematical objects.

[Proof, literal C3 addresses, and the full response](square_graceful/ADDITIVE_CARRY.md).

**Scope:** this rules out simultaneous five-core failure on a language
containing every finite binary pattern as a prefix. Padding changes the
original target. It does not show that the original unpadded address,
or an arbitrary odd N, is covered. The aligned lower counts are0 modulo3
and the upper counts are2 modulo3; the residue1 lane requires a different
rewrite or another operation.

Eight additional exactly classified rewrite languages precede this
construction. In particular the H13 upper count at B=2^k−2, k≥2, has
H9 and H5 source realizations at the same target

    77·2^(k−1)−31.

The H13-to-H5 root-state difference is `(0,8,8,8)`, so equality survives
every subsequent Z. Both distinct H9 cores participate; they agree on
this one-sided response and separate again when a left address is added.
[Eight iff languages](square_graceful/GRAMMAR.md),
[ray, shell fibres, and the divisor fold](square_graceful/PROOF.md).

## 2. A common order for all later constructions

For a rooted tournament use `(H,E,U,d)`: total paths, paths ending at
the root, root-start ordered two-path covers, and root-deletion count.
The invertible coordinates

    q=(h,u,v,t)=(H,U,U−d,H+U−E)

turn the inherited physical operations into

    Oq=(h,2u,u+v,h+2u),
    Zq=(t,u+v,2v,t+2v),
    Jq=(t,h+u,u,h+2u).

**Theorem:** q(x)≥q(y) coordinatewise if and only if every common future
O/Z/J word has H(x)≥H(y). The necessity follows already from four
families of O/Z words; sufficiency follows from the three nonnegative
matrices. Thus arbitrary further root changes are included.

[Proof and exact equality languages](dynamics/PROOF.md).

This criterion changes the interpretation of source simplification.
Throughout the repeated-9 simultaneous rebase above, the H13 source is
cheaper under future continuations than its same-count C3 replacement:
the q difference is `(0,Δ,Δ−10,Δ−10)` with Δ≤−6. Reducing the deletion
count alone is therefore not a valid induction measure.

## 3. Actual changes of source and paired resynchronization

The two different birth cuts recovered from the H13 pair move its
difference from H into E:

    (31,7,45,13) / (33,7,45,13)
        -> (73,9,147,25) / (73,7,147,25).

One common Z reveals the difference again at211/213; one O erases it.
This operation arose from the old modular/Monster shared-quotient ideas:
the two cuts need compatible responses, not identical masks.

At211/213 there are now two actual source representations:

    old: (H,13,347,13),   new: (H,25,269,25).

Their difference in q is `(0,78,90,90)`. The new source is smaller under
every future word containing Z or J, even though its minimum deletion
count rises from13 to25. It also uses eight vertices instead of nine.

Different future addresses can restore the same target. Starting at211,
old address79 and new address102 both produce26649; the upper branches
both produce26651. Their U values are42843 and33157. This is a genuine
second paired synchronization, not an arithmetic sum of independently
available differences. A declared finite address box gives227 meetings,
152 of which improve U. Those finite counts are not a repeatability theorem.

[Native birth and all-address formulas](modular_cycles/PROOF.md),
[resynchronization and root-change transport](modular_cycles/RESTORATION.md).

The mixed-word search produced exact length14 loops whose13/25-weighted
linear drifts cancel. Repeated target agreement also requires matching
the exponential coefficient and the constant. None of the227 current
seed meetings satisfies both remaining equations for the tested1,568
loop pairs. The full three-equation iff criterion is retained in
[construction/PROOF.md](construction/PROOF.md), rather than replacing
seed compatibility with drift cancellation.

## 4. Every even difference is now a fixed-core physical switch

The original H13 construction has a pure-H difference
`3A−(A mod2)+2` at every pair of flank addresses. Its whole rooted
sidecar is explicit. One J turns this difference into a scalar multiple
of the rooted singleton response; repeated J therefore transports it
by the exact tribonacci recurrence1,1,3,5,9,17,... .

The asymmetric H25 birth is stronger. At B>0 its two sources have equal
E,U,d and named difference

    −3A−(A mod2)+2(B mod2).

Choosing A and B=1 or2 realizes **every positive even width**, with the
same deletion count25. For odd A and positive even B the birth reverses
the order of the H13 pair while preserving its width. This is an actual
operation outside common positive O/Z/J evolution.

[Full H13 response and tribonacci transport](dynamics/PORTABLE_GAPS.md),
[fixed-d25 all-even-width construction](construction/PROOF.md).

Width realization still comes with a particular baseline. It is not a
license to add or subtract that width from an arbitrary existing count.
The rooted-amalgamation calculation in
[ADDITIVITY.md](modular_cycles/ADDITIVITY.md) supplies a genuine weighted
addition operation and identifies the remaining compatibility conditions.
The H13 and pentagon switches give the actual pure-H rectangle
`2873,2899,2903,2929`, with independent increments26 and30 and one common
root sidecar. Finite chains give additive cubes with explicit baselines.
The same note extends the exact q order to arbitrary right-root amalgams,
and gives a sufficient seven-coordinate order for fixed contexts on both
sides. The latter retains a precise reason why fixed contexts alone cannot
repay the211 reset's count difference.

## 5. Incoming work supplied a different native repair

Session254's actual H37 arc reversal changes its path-colour polynomial by
`2t²(1−t²)`. Wrapping the six-vertex source with a retained global-source
root produces

    (37,0,376,37) -> (37,0,352,37).

It preserves H, the root, and all other arcs, and decreases all three
future-response coordinates by24. Appending the same arc patch as an
ordinal suffix to any rooted carrier X gives the corresponding U decrease
`24U_X`, at total count37H_X. This imports an incoming proved operation
into the present source order without re-auditing the older result or
silently importing its separate modular hypotheses.

[Actual rows, wrapper bijection, and arbitrary contexts](construction/PROOF.md).

## 6. What the broader historical search contributed

| Recovered work | Concrete transfer in this packet |
|---|---|
| Square-sum reflection orbits and paired base49 construction | Reflect shell fibres, retain both baselines, and search for one digit defect shared by two actual witnesses |
| Graceful gap-tail insertion and marked reconstruction | Keep the root response and the outstanding integer defect when transporting a target equality |
| `{p,p³,p²qr}` divisor folds | Four-bit energy descends to the p²qr exponent box; one lost orientation bit restores the actual count |
| Fifteen-prime / level-two supersingular work | Equal coarse endpoints can hide a selectable lift; use a fibre product of birth cuts over a common response |
| Monster/Baby Monster level refinements and common-C2 actions | Pair compatible local choices; the native73 pair hides the difference in E and a later step reveals it |
| Quadratic three-cycles, order-six lifts, and changing hexagonal forms | Seek one common positive frame for all source operations; this yields the exact q order |
| Collatz honeycomb expansion and shared halving tails | Retain the joint source data and separate clocks instead of identifying two branches solely by their endpoint |

The source paths, theorem IDs, actual maps, lost information, and tests
are recorded in the individual notes. These are methodological transfers
with resulting formulas and graph operations. No identity between these
different mathematical subjects is assumed.

## 7. The precise remaining proof obligation

The results solve several compatibility problems on infinite families,
but they do not imply an arbitrary-target descent. A completion still
needs a selection law which, given an arbitrary candidate odd N, either
constructs it directly or chooses a smaller admissible target together
with an actual operation returning **exactly N**.

The current carry cancels a common defect by changing the address and its
target. The current source reset preserves a target at its base but changes
most later targets. The current switch realizes any even difference at a
constrained baseline. These three operations now have exact response laws;
joining them at the same baseline is the unresolved implication.

The promising next object is a coupled target/defect transition with a
positive realization and a decreasing integer measure. The asymmetric
birth's order reversal is relevant because common positive contexts cannot
cancel certain directed source debts. A proof must specify the joint cut
and both resulting totals, not only their difference or their residues.

No new global coverage bound, conjecture completion, or canon promotion is
claimed. New identities have literal small-graph and exact arithmetic
controls; inherited large censuses were not rerun.
