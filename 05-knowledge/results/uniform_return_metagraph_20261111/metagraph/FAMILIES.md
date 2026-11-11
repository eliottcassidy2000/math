# Capped tournament metagraphs: finite rules and infinite source families

**PROVED IN NOTE, relative to the cited/inherited strong-size bound and five-core
classification. FINITE-EXACT dictionary and controls.** This classifies the
declared operations at every order under a fixed Hamilton-path cap. It is not a
classification uniform in the cap, or a proof of Chan–Pak.

## 1. Inheritance, carrier, and concept board

Closest mechanisms: the [fixed-H SCC species](../../spectrum_basins_20261020e/spectrum.md),
[actual native frame carrier](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md),
[native reset obstruction](../../adaptive_state_induction_20261103/METAGRAPH_RESETS.md),
and [complete five-core classification](../../chanpak_five_core_bridges_20261110/history/RECOVERY.md).
The existing [A000568 composition](../../session86_relations_attachment_carriers_20261028/A000568_COMPOSITION.md)
already gives unrestricted Burnside enumeration and aggregate prime-index
recovery; those are not new here. The older [82-sequence census](../../metagraph_aspect_sequences_20261111/README.md)
already distinguishes scalar target unions from compatible sources.

Concurrent incoming work in [infinite birth families](../../infinite_birth_families_20261111/README.md)
independently supplies the consecutive-SCC birth classification and compiles946
fixed-H birth families below125. That birth mechanism receives concurrent
credit here. The present operation package additionally retains arc reversals,
the finite capped native part, marked padding reachability, neutral components,
and the sharp excursion height. Its capped-birth series agrees with all50
tested rows of the incoming compiler; this is a cross-check, not a rerun of
the incoming primitive census.

Hostile: different H9 sources can have different operations despite equal H.
Corrected near miss: finite grammar saturation is not an all-height proof.
Underused sidecar: the ordered SCC word, including every singleton position,
together with a rooted local rule and its stabilizer.

| Concept | Retained structure | Decisive test |
|---|---|---|
| Bounded H | Actual strong SCC types and gaps | Can every move be localized? |
| Native complement | Actual Hamilton frame, loops, colour | Does padding preserve a move under the cap? |
| Birth/deletion | Newborn/root orbit and transport | Do orbit and named multiplicities agree? |
| Arc reversal | Actual unordered pair inside its SCC interval | Can its affected interval be unbounded? |
| Source return | Allowed intermediate H values | Are equal-H sources in the same neutral component? |
| Enumeration | Classes, rooted classes, frames separately | Do independent literal graphs match the series? |

An ordinary vertex is an isomorphism class of a tournament T, with raw path
count H(T). Operations are: one-vertex birth/deletion; reversal of one arc;
native K_P, which preserves a Hamilton path P and reverses every other arc;
and, when stated, converse. No ties are oriented into an artificial tournament.
For native lines we retain normalized frames, colour and literal multiplicity.
For birth edges we retain pointed children (U,v) up to simultaneous isomorphism.
These are different edge carriers and are never silently identified.

## 2. Every fixed cap has a finite alphabet and exact local rules

Fix B>=1. Let S_B be the finite set of nontrivial strong ordinary types C
with H(C)<=B. Write x for a singleton. Every tournament of H<=B is uniquely
an ordered word in S_B union {x}, with product of the letter H values <=B.

The alphabet is finite because Busch's inherited bound is

    b(3q)=3*5^(q-1), b(3q+1)=5^q, b(3q+2)=9*5^(q-1),
    H(C)>=b(|C|).

Let m_B be the largest possible nontrivial strong order under this bound
(take zero if B<3). There are at most floor(log_3 B) nontrivial letters.
Thus the language is a finite union of expressions

    x* C_1 x* ... C_r x* .                                  (1)

This normal form and its rational class series are inherited. The following
operation classification is the present extension.

**Birth/deletion.** Besides insertion/deletion of one x, all births are finite
rules

    SCCword(C-v) -> C,                                     (2)

where C is in S_B and v ranges over its vertex orbits. Deletion reverses these
rules. To prove completeness, adding a vertex cannot split an old SCC. Its new
SCC either is a singleton inserted at a cut, or contains the newborn and a
consecutive interval of old SCCs. Removing the newborn gives exactly (2), of
total old order at most m_B-1. Exterior components remain unchanged.

The two words in a rule are placed in the same exterior context, and both
whole-word products must satisfy the cap. Parent H never exceeds child H:
each parent Hamilton path can be extended by inserting the new vertex, and
erasing the inserted vertex recovers the original path, giving an injection.

**Arc reversal.** If its endpoints lie in one strong SCC, only that finite
letter changes. If they lie in different SCCs, their entire intervening SCC
interval becomes strong, except that two adjacent singletons merely exchange
their order. Indeed the reversed last-to-first arc, original forward arcs,
and paths within the endpoint SCCs connect every vertex of this interval in
both directions. For at least three interval vertices there is a forward
first-to-last path avoiding the reversed arc. Hence every cap-preserving
nontrivial cross-SCC reversal acts on an interval of total order <=m_B.
The exceptional rule is xx->xx, retaining its actual exchanged names.
Enumerating all such finite intervals and endpoint orbits is complete.

**Native complement.** Every native move from a nonstrong tournament of order
n>=3 lands in a strong tournament: its path's first vertex initially beats
its last, so K_P adds a closing arc to the retained Hamilton path. Thus at
least one endpoint of every native line is strong. If both H values are <=B,
then n<=m_B. The complete native part of the capped metagraph is finite.
This strengthens the inherited h<n neutral obstruction to the exact strong
size cutoff; the weak-to-strong lemma itself is inherited, not rediscovered.

Converse reverses the word and replaces each strong letter by its converse.
It need not act trivially on ordinary classes, even when each letter is
self-converse, because it reverses the gap positions.

## 3. Composition and retained names

Each birth/deletion/arc rule is a replacement of a bounded contiguous word
inside an unchanged left/right context. The cap is a finite multiplicative
state. Consequently a finite transducer recognizes the exact support relation
for each of these operations; finite native relations can be included. A
prescribed finite composition is obtained by composing
these transducers and projecting intermediate words. This supplies effective
tests for regular source families and fixed finite programs. Arbitrary-length
reachability is not inferred from closure under finite composition. Converse
is handled separately by reversal of the queried regular language and the
converse covariance of the local rules. A one-way transducer for arbitrary
word reversal is not being claimed.

There is also an exact finite quotient when free singleton padding is allowed.
Fix r distinct persistent names and a fixed induced tournament on those names.
Forbid deleting a name or reversing an arc between two names. Delete all
UNNAMED singleton SCCs to obtain a marked kernel. The nontrivial strong part
has at most floor(3 log_3 B) vertices by the inherited strong-product bound;
at most r additional marked singleton SCCs remain. Thus there are finitely
many marked kernels.

All states over one kernel are connected by insertion/deletion of unnamed
singleton SCCs, preserving H and the entire marked induced tournament. Every
other local step has a representative with at most m_B additional singleton
vertices: erase all unused exterior padding, retaining the bounded affected
interval. Native steps already have total order <=m_B. Hence a finite marked
kernel graph exactly determines reachability after these padding moves are
allowed. A kernel walk lifts: between consecutive edges, adjust padding inside
the current kernel fibre, then apply the stored actual rule. A concrete walk
projects by forgetting unnamed singleton gaps.

This quotient preserves source-return existence for the stated operations,
cap, persistent names and allowed H set. It loses gap lengths and step counts;
their actual words and transport maps restore execution. The construction can
be enormous, and no efficient universal-cap bound or 141-family theorem follows.
For a strong source with H=34*2^k+5, the cap itself grows with k.

## 4. The complete H<=9 language

Use the actual outgoing-row representatives

    A=C3=(2,4,1),             B=S4=(2,4,9,3),
    C=S5a=(2,4,9,19,7),       D=S5b=(2,4,25,3,11).

The inherited all-order classification gives strong letters x,A,B,C,D,
with H=1,3,5,9,9; their automorphism orders are 1,3,1,1,1. All are
self-converse. The complete ordered nontrivial kernels are precisely

    empty, A, B, C, D, AA.                                 (3)

Thus the five strong cores require a sixth composite family when classifying
all parents. The arbitrary gaps in AA are independent; they must not be
collapsed into one transitive block.

`classify_cap9.py` gives the COMPLETE local dictionary: 15 rooted strong-birth
rules, 32 directed arc-orbit rules (including xx->xx), and 14 literal native
lines of orders3..5. Counts here retain the stated marking and multiplicity.

The birth table, grouped only for display, is

| Rule | Root-orbit multiplicity | Named attachments per root-orbit |
|---|---:|---:|
| xx->A |1|1|
| A->B |2|3|
| xxx->B |2|1|
| B->C |2|1|
| Ax->C, xA->C |1 each|3 each|
| xxxx->C |1|1|
| B->D |3|1|
| xxxx->D |2|1|

The singleton insertion rule applies at every SCC gap. If U's root orbit has
stabilizer a_(U,v), and deleting it gives T, then its named attachment weight
is |Aut(T)|/a_(U,v), while its named deletion weight is |Aut(U)|/a_(U,v).
This follows by quotienting the isomorphisms T->U-v by the rooted stabilizer.
It distinguishes ordinary edge orbits from actual attachments.

The 14 native lines split by (order,colour) as

    (3,blue):1; (4,blue):2; (4,black):2;
    (5,blue):1; (5,black):8.

There are no capped native lines at any order >=6, including loops. For
example the C-D native lines at order5 cannot simply be padded and remain
under this cap. The physical C-D one-arc switch can be padded as an ordered
context. Native and local replacement composition are different operations.

## 5. Exactly five neutral source components

Consider only operations from section1 whose every intermediate count equals
the starting count, allowing arbitrary order, ordinary isomorphism, and converse.
Then the entire H<=9 metagraph has exactly FIVE connected components:

    {empty kernel}, {A}, {B}, {C,D}, {AA},                   (4)

where each brace includes all its singleton-gap families. For nonempty
tournaments the first component is the family TT_n, n>=1; allowing the empty
tournament adds it to the same component.

Proof: padding connects all gap vectors of each kernel. The finite dictionaries
show that C and D have a physical one-arc switch, while every other neutral
rule stays in its listed kernel. Native lines add only B's loop and the C-D
connection; converse reverses gaps and fixes these kernel classes. The complete
rule theorem excludes unseen larger-order bridges. This is a classification
of actual operation components, not a bijection with the five strong cores.

In particular the two-triangle H9 source is permanently separated from both
strong H9 source types if every intermediate H must remain9.

## 6. A sharp excursion height: 9 -> 17 -> 9

Relax neutrality by allowing H to rise, but require H>=9 at every step.
The minimum possible peak permitting an AA source to return to a C or D
source is exactly17, uniformly over all singleton gaps after free padding.

**Lower bound.** The three kernels AA, AB, BA form a closed component in the
band9<=H<=16. From AA, operating inside one triangle either preserves it,
drops H to3, or enlarges it to B, giving15. Merging the two triangles by a
birth gives a strong component of order>=7 and H>=25. Reversing an arc
between them gives17 when they are adjacent, and at least25 with an intervening
singleton. Creating any further nontrivial SCC multiplies H by at least3.
From AB or BA, an internal move either retains the kernel, returns to AA,
or drops below9. Any enlargement/merge has H>=25. Native moves from AB/BA
have strong order>=7 and H>=25; from adjacent AA they give17, and padded AA
has strong order>=7 and again H>=25. Converse only swaps AB and BA.

**Upper bound.** Start with the actual adjacent AA rows

    (58,60,57,16,32,8).

Reverse the arc between names0 and3, giving

    (50,60,57,17,32,8), H=17.

Delete name2; the remaining five-vertex tournament is C, with H=9.
The cross reversal creates exactly four new odd directed cycles; none can
be paired with a disjoint odd cycle. The inherited odd-cycle-collection
identity therefore gives H=9+2*4=17. Equivalently the four relevant cross
paths have left/right vertex lengths (1,2),(2,1),(2,3),(3,2).
All AA Hamilton paths are equivalent under Aut(C3)^2, so its native partner
also has H17. The literal checker verifies both readings and the deletion.

This route retains five original vertices but changes one arc between them;
it is not an immutable-five-vertex reset. If a count drop is allowed, delete
two vertices from one triangle and use xA->C or Ax->C: the route9->3->3->9
retains the remaining four old vertices and all their arcs. These are distinct
return budgets and distinct preservation predicates.

## 7. Exact refinements of A000568, at every order

Let a_n^(9) count ordinary classes with H<=9, with the empty tournament at n0.
The inherited SCC species specialized to (3) gives

    F(z)=1/(1-z)+(z^3+z^4+2z^5)/(1-z)^2+z^6/(1-z)^3.

For every n>=4,

    a_n^(9) = (n^2-n-4)/2.                                (5)

Rooted ordinary classes count vertex ORBITS. A kernel C_1...C_r contributes
v(C_1)+...+v(C_r) plus the number of singleton vertices, where v(C) is its
number of vertex orbits. Here v(A)=1,v(B)=4,v(C)=v(D)=5. Therefore, for n>=4,

    r_n^(9) = (n^3-5n^2+28n-72)/2.                        (6)

Every deletion stays under the cap, so r_(n+1)^(9) is exactly the number of
birth-edge orbits from order n to n+1, with the pointed-child convention.
It does not count distinct unmarked endpoint pairs or named masks.

The total normalized Hamilton-frame mass is

    f_n^(9)=sum_[T] H(T)/|Aut(T)|=(n^2+39n-156)/2, n>=4.  (7)

For n>=6, the automorphism strata are

    |Aut|=1: 3n-10;  |Aut|=3: n-2;  |Aut|=9: binom(n-4,2).

The self-converse class series is

    1/(1-z)+(z^3+z^4+2z^5)/(1-z^2)
             +z^6/((1-z)(1-z^2)).

At n>=6 its coefficient is n/2 for even n and (n+3)/2 for odd n.
Converse-merged counts are half the sum of this coefficient and (5).

| n | ordinary | vertex-rooted | normalized frames | converse-merged |
|---:|---:|---:|---:|---:|
|4|4|12|8|3|
|5|8|34|32|6|
|6|13|66|57|8|
|7|19|111|83|12|
|8|26|172|110|15|
|9|34|252|138|20|
|10|43|354|167|24|
|11|53|481|197|30|
|12|64|636|228|35|

There is a complete action inventory, not only an automorphism-order list.
At n>=6 the aggregate vertex-action cycle index is

    Z_n = (3n-10) p1^n
        +(n-2)(p1^n+2 p1^(n-3) p3)/3
        +binom(n-4,2)(p1^n+4 p1^(n-3) p3
                            +4 p1^(n-6) p3^2)/9.          (8)

The groups are respectively trivial, C3, and C3 times C3, with the last
two factors acting on the two intrinsic triangle components. In particular
the group of order9 is not an unidentified group of that order. Substituting
p_j=q counts q-coloured vertex orbits. Substituting p_j=1+t^j counts subset
orbits by subset size. These substitutions retain the actual permutation
action, so they give every finite-palette colouring refinement in this slice.

For r distinct ordered persistent names, let (n)_r be the falling factorial,
zero if n<r. Burnside gives the all-r formula, valid for n>=6,

    a_(n,r)^(9) = (3n-10)(n)_r
      +(n-2)((n)_r+2(n-3)_r)/3
      +binom(n-4,2)((n)_r+4(n-3)_r+4(n-6)_r)/9.             (9)

A permutation fixes a named tuple precisely when it fixes each named vertex.
The formula therefore counts named-root classes, not unordered root sets.
It recovers (5) and (6) at r0 and r1 and is eventually a polynomial of degree
r+2. The finite named-kernel state construction in section3 has a directly
enumerated sequence of state sizes, rather than an unspecified marking count.

Every arbitrary one-vertex attachment is a binary colouring of its parent.
Thus the number E_n of attachment ORBITS from all H<=9 ordinary parents,
with no restriction on the child's H, is Z_n(2,2,...), namely

    E_n = 2^(n-3)(n^2+19n-68), n>=6.                      (10)

For n6..10 this gives656,1824,4736,11776,28416. The capped-child count was
instead r_(n+1)^(9), for example111 at n6. Neither number is the count of
distinct ordinary unrooted children: forgetting the newborn can identify
different rooted cards. The aligned root-switch carrier inherited in
section1 is required for that final quotient.

These are structural refinements of A000568, not new unrestricted A000568
terms or a speed record. Increasing B requires new strong types and new local
rules. The unrestricted prime/module inventory is not replaced by these six
kernels. The marked-source obstruction, rather than total class count, is the
part directly relevant to Chan–Pak.

## 8. Reproduction and boundaries

Run `python -B classify_cap9.py --output cap9.json` from this directory.
The receipt retains every local row/vertex/arc witness. It checks all33,867
labelled tournaments of orders1..6 against the family description, compares
subset-DP and permutation path counts for every graph through order5, checks
1,220 literal deletions through order10, and evaluates all series through50.
The all-order cutoff uses the inherited audited classification; finite checks
do not supply that quantifier. `sequences.tsv` exports the ordinary, rooted,
frame, converse, arbitrary-attachment and two/three-named refinements through50.

`python -B audit_operations.py` independently uses a first-vertex recurrence
to check3,607 forward attachments and4,592 arc reversals, including502 and1,644
operations remaining under the cap. It verifies localization, both literal
9-return paths, and the neutral-kernel quotient. Its normal and optimized runs
agree. Direct graph automorphism enumeration and literal subset/named-tuple
orbit sets through order7 audit the cycle-index consequences. The incoming
fixed-H birth compiler comparison pins its source hash in the receipt.
`audit_defect.py` separately audits the root's fixed-defect grammar by
341 literal permutation/cut counts and independent strongness searches.

Connection contract: actual tournament -> ordered strong letters and singleton
gaps is lossless on ordinary classes. Forgetting gaps preserves H, kernel type,
and padding-adjusted reachability, but loses order and literal move incidence.
Rooted rules, gap words, stabilizers and named transports restore the executable
side. Projection farther to H loses the neutral-component distinction (4).
The sharp9/17 example is the cheapest test for a proposed source-return quotient.
