# Selected-set module surgery: ternary multipliers, sharp mark positions, and relocated roots

**Session260, 2026-11-11. PROVED IN NOTES with exact finite coefficient certificates; bounded relocation atlas FINITE-EXACT. Chan–Pak remains OPEN.**

## Inheritance and the changed interface

Inputs: Session258's `seed_decision/CORE_CLASSIFICATION.md` and its definitions
of H,W,U; the selected-set tail identity; the latest
`spectrum_multiplicity_signed_tables_20261111` fixed-gadget and actual-return
results. We treat those results as inherited, rather than spend this session
rechecking their classification. A positive seed means

    D(B,O):=H(B)−W_B(O)=1,     U_B(O)>0.

Here W counts ordered two-nonempty-path covers with first end and second
start both selected; U permits an empty first path and requires the second
start selected. Such a seed supplies H=1+2^s U for every s≥1 by the inherited
literal tail constructor. If U=c2^j, the c-ray starts at exponent j+1.

Closest mechanism: the strong-core seed normal form. Hostile: fixed q=7
has no seed at any order, and the all-realizer +2 birth theorem prohibits
naive H→H+2 recursion above H=3. Corrected near miss: arbitrary integer
multipliers are not selected-set seeds, and a finite tail does not cover its
below-start head. Least-used sidecar here: **where the selected vertex sits
inside a transitive module**, plus permission to move the selected pair after
changing the underlying graph.

Concept board: strong seed core / ordered path-cover response / marked
position / ternary multiplier digit / multiple module replacements /
below-start heads. Anchor: produce actual new all-length constructors.
Niche: identify the exact mark-position obstruction. Wildcard: relocate both
marks after simultaneous replacements instead of preserving a root name.

All rows below are outgoing-neighbour bit masks, with zero-based vertices.
Define two strong five-vertex cores

    A=(20,1,2,7,14),      B=(24,1,3,6,14).

Replacing a vertex by a nonempty transitive module preserves strong
connectivity: every exterior connection is inherited, and directed paths
through the strong quotient connect all module vertices. A selected clone
is distinguished; unselected twins are not silently declared interchangeable
with it.

## 1. A finite coefficient certificate for unbounded module lengths

Let k exterior vertices be fixed, and replace an unselected vertex type by a
transitive chain of m vertices, all with identical exterior neighbourhood.
For a fixed ordered p-path-cover observable whose endpoint predicates refer
only to exterior vertices, compress each maximal chain run to one placeholder.
If there are r runs, r≤k+p: each path can have at most one more chain run than
its exterior vertices. A skeleton with r canonically ordered run places lifts
in exactly

    F_r(m)=r! S(m,r)

ways, where S is a Stirling number of the second kind. Indeed, partition the
m ordered chain vertices into r nonempty ordered parts and list each part in
its unique increasing order. Conversely restriction to the runs recovers
that ordered partition. Thus every such observable is

    Q(m)=sum_(r=0)^(k+p) a_r F_r(m),

with nonnegative integer skeleton coefficients. Differences such as D may
have signed coefficients. At m=0 use F_0(0)=1 and F_r(0)=0 for r>0.
The theorem allows an empty path when the observable does; p is an upper
bound on the number of paths, not a demand that all paths be nonempty.

There are two exact ways to obtain the coefficients:

1. Evaluate Q(m) for m=0,...,k+p and triangularly invert the integer matrix
   F_r(m), whose diagonal is m!.
2. Replace the module by r **mutually nonadjacent twins** instead. These are
   an auxiliary directed graph, not a tournament or a claimed construction.
   Every run is now one twin. Its observable equals r!a_r, because labels may
   be assigned to canonically ordered run places in r! ways.

`module_ladders.py` implements both and requires agreement. This is an exact
finite certificate for an unbounded identity, not extrapolation from an
arbitrarily truncated sequence. The degree bound and skeleton theorem are
what make the certificate conclusive.

For a first-selected chain, keep its first clone as a fixed exterior
vertex, and treat the remaining m−1 clones as the unselected module. This
changes the parameter and degree bound; it must not be confused with
unselected cloning. The related Hamilton-only substitution formula is an
inherited mechanism; using it jointly for H,W,U and moving the selection is
the new construction interface here.

## 2. Six explicit all-length strong seed families

The following are exact for every m≥1. In the first three rows replace the
specified unselected vertex by an unselected m-chain and retain O. In the
last two replace the selected vertex by an m-chain and select only its first
clone, keeping the other member of O.

| Core, O, replaced vertex | H | U |
|---|---|---|
| A, {1,2}, 3 | 3·2^m+3 | 7·2^m+6 |
| A, {1,3}, 2 | 3·2^m+3 | 4·3^m+9·2^m+2 |
| A, {3,4}, 2 | 3·2^m+3 | 8·3^m+12·2^m−4 |
| B, {1,2}, 1, first selected | 6·2^m−1 | 10·3^m−1 |
| B, {3,4}, 3, first selected | 6·2^m−1 | (50·3^(m−1))−3 |

Every row has D=1. The last two give the exact coupled recurrences

    (H,U) -> (2H+1,3U+2),
    (H,U) -> (2H+1,3U+6),

respectively. They are physical seed operations. The +2 or +6 is on the
**multiplier U**, not a forbidden +2 one-vertex Hamilton-path increase.

For transparency, the nonzero skeleton coefficient lists (append zeroes to
the degree bounds) are:

- First row, parameter m: H=(5,9,3), D=(1,1), U=(11,20,7).
- Second row, parameter m: H=(5,9,3), D=(−1,1), U=(15,32,21,4).
- Third row, parameter m: H=(5,9,3), D=(1,1), U=(16,44,36,8).
- Fourth row, parameter t=m−1: H=(11,23,12), D=(1,1), U=(29,89,90,30).
- Fifth row, parameter t=m−1: H=(11,23,12), D=(1,1), U=(47,147,150,50).

The deletion boundary in row2 has D=−1, whereas every positive chain length
has D=1. This is an instructive failure of extending a positive-length
identity to the absent-module boundary. The first three H formulas also do
not apply at deletion. The receipts retain those boundary values explicitly.

The proof is the skeleton certificate plus
F_1(m)=1, F_2(m)=2^m−2, F_3(m)=3^m−3·2^m+3 for m≥1.
All remaining coefficients vanish at the proved degree cutoff. Positivity
of U is immediate from the displayed formulas.

## 3. Complete marked-position theorem, not just a first-clone trick

In B keep original vertex2 selected. Replace original vertex1 by

    TT_j -> z -> TT_k,

where only z is selected in this module. All exterior arcs are inherited;
original vertex2 is the second selected vertex. Let 0^0=1 and 0^j=0 for j>0.
Then for every j,k≥0,

    H = 12·2^(j+k)−1,
    D = 0^j+2^j−1−6(3^j−2^(j+1)+1)2^k,
    U = (18·3^j+12·2^j)3^k−1.                 (1)

There is a second case with the same H,D: select original vertex4, replace
original vertex3 by TT_j→z→TT_k, and select z. Here

    U = (30·3^j+20·2^j)3^k−3.                 (2)

**Sharp consequence:** in either case, this is a unit-offset seed if and
only if j=0 or j=1. For j=0 the expression is1; for j=1 the coefficient of
2^k vanishes and 2^j−1=1. For j≥2 put f_j=3^j−2^(j+1)+1. We have
6f_2=12>2^2−2 and f_(j+1)=3f_j+2^(j+1)−2. Inductively
6f_j>2^j−2, hence D≤2^j−1−6f_j<1. No k can repair it.

Thus each marked carrier has two infinite ternary multiplier ladders:

    30·3^k−1,   78·3^k−1;
    50·3^k−3,  130·3^k−3.

These are genuine constructors but **not a free ternary digit alphabet**.
In a repeated cloning process that keeps the other selected original vertex,
the selected clone's prefix length never decreases unless the mark is
relocated. A first-to-second move can occur at most once while preserving
the seed predicate; a second-to-third move destroys it. The lost coordinate
in a purported scalar ternary recursion would be this prefix length.

### Exact symbolic proof of (1) and (2)

`marked_position.py` enumerates compressed path words on the five fixed
vertices and two repeatable symbols P,S. The only symbol-to-symbol edge is
P→S. Each fixed vertex appears at most once; a gap between fixed vertices
has at most two symbols. Therefore this is a finite symbolic enumeration,
not a truncation of arbitrary word length. It produces the joint skeleton
polynomial in F_r(j)F_s(k) for H,W,U. The two ordered path restrictions are
convolved over disjoint fixed-vertex subsets; the empty first path is added
only for U.

Expanding F_r(j)=sum_a (−1)^(r−a) binom(r,a)a^j gives the exact exponential
coefficient tables in `marked_position_receipt.json`. For H the nonzero
terms are −1 and12·2^j2^k. For D they are
0^j−1−6·2^k+2^j+12·2^j2^k−6·3^j2^k. For U they are exactly (1),(2).
This is a reproducible finite proof of the bivariate identity. An additional
25 literal chain replacements per case check the compiler's interface.

## 4. Root relocation after simultaneous module changes

Keeping the old pair fixed is restrictive. The separate experiment
`relocation_atlas.py` substitutes a transitive module for **every** vertex
of A or B, then finds every selected pair satisfying D=1. The finite universe
is explicit: every module length1,...,5, total order at most12, and selected
set size exactly2. These are labelled template presentations, not distinct
isomorphism classes and not all strong cores of those orders.

The exact output has627 seed presentations and212 distinct odd parts of U.
This is not a new claim that these scalar targets were previously unknown:
Session258's mixed constructor bank already covers many or all of these
individual slopes by other means. It is new evidence about **which pair
relocations realize the seed equations inside this particular module
carrier**.

Three useful literal examples (module lengths in original-vertex order;
selected indices refer to the expanded graph) are:

- c=7: A lengths(1,2,1,2,1), O={1,3}, H=33, U=112=7·2^4.
- c=13: A lengths(1,2,1,3,1), O={1,3}, H=61, U=208=13·2^4.
- c=21: A lengths(1,1,1,2,2), O={1,3}, H=37, U=168=21·2^3.

They therefore give the respective rays from exponents5,5,4. This does not
contradict the fixed-q=7 or q=21 exclusions: the exact multipliers differ.
The receipt also supplies literal tournaments for all below-start odd ray
values in these three examples, so their heads are not silently discarded.
The finite head search is a positive constructor search, not an absence
certificate or a universal bound on witness order.

The same box still misses slopes3,9,19,23,27 and others listed in the
receipt. This is a limitation of that finite box, not an all-order exclusion
for these two templates, much less an exclusion of the target counts.
In particular, the unrestricted position theorem above is only for one
specified selected pair/carrier; it does not rule out the relocated pairs.

### A relocated root becomes an infinite second binary ladder

The c=7 and13 examples lie on one all-length family, not two isolated hits.
Inflate original vertex1 of A to length2, selecting its first clone and
original vertex2. This gives the six-vertex marked core

    A2=(40,5,1,6,15,30),    O={1,3}.

Now replace its unselected vertex4 by an m-chain. For all m≥1,

    H=7·2^m+5,      D=1,      U=24·2^m+16.

The independent skeleton coefficient certificate is included as the sixth
row of `module_ladders_receipt.json`. Its nonzero lists are H=(11,19,7),
D=(1,1), U=(37,64,24), with deletion recorded separately. For m≥2,
U=16(3·2^(m−1)+1), so this constructs the odd-slope family

    c=3·2^r+1, r≥1,

always beginning at ray exponent5. In particular it gives c=7,13,25,49,97
and so on using a fixed seed interface. It satisfies the exact recurrence
(H,U)→(2H−5,2U−16); the negative affine constants are compatible with
positive underlying path-cover counts. They do not represent subtraction
of tournaments. **Incoming integration (48ffb162415d):** the now-current
`seed_recursion_births_20261111/PROOFS.md`, R3, proves every below-start
head for this entire family. Its explicit head kernels give the complete
grid H=1+(3·2^r+1)2^s for every r,s≥1; it also closes the separate
21·2^r+1 grid. We inherit that completed-head theorem and do not claim our
sixth family or scalar grid as a new coverage result. Our explicit
six-vertex coefficient certificate is an additional presentation of the
seed interface. The incoming R1/R2 classify fixed-mark stable inflations;
they expressly leave reselection outside their scope, matching our new
mark-position and selected-pair-relocation interface.

## 5. The next lemma needed for Chan–Pak

The experiments favor **changing the selected pair together with the
module graph**, rather than a fixed-root numeral machine. A concrete
remaining target is:

> Find a finite family of marked strong-core types and explicit module
> replacements plus selected-pair relocations, with a jointly maintained
> response vector, whose multiplier odd parts admit a terminating descent
> for every positive odd c. State and prove the admissible predecessor
> guards, rather than just a congruence covering by the scalar maps.

The six ladders and two-position classification supply actual edges and
an exact obstruction in such a graph. They do not yet supply this uniform
closure. A solution must also handle each resulting below-start ray head;
the three examples show how to discharge that finite obligation when an
actual start bound is known. General uniform seed existence, uniform
below-start completion, and Chan–Pak remain OPEN.

The complementary bridge packet supplies a prime/Carmichael congruence
constraint for integer Stirling expansions. It applies to our unselected
clone parameter without modification. In the first-selected rows the
variable is m−1, so that congruence compares to two total clones, not to the
original one-clone core. Indexing matters.

## Reproduction

From the repository root:

    python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/constructors/module_ladders.py
    python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/constructors/marked_position.py
    python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/constructors/relocation_atlas.py

No external packages or inherited computational imports are required.
Checks raise explicit errors and do not rely on Python assertions. The
first two are finite exact certificates for all-length identities. The
third is a finite parameter-box atlas and positive finite-head witnesses.
