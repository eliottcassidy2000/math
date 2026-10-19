# The 105 optimal atlases: dual matchings, a tetrahedron, and the retained phase

**Status: PROVED elementary classification and operator identities;
FINITE-EXACT complete labelled dictionaries and inherited group-action replay.
Geometric selector realization remains OPEN.** No priority claim.

This is the follow-through from
[four charts and the mod-four codec](gray_selector_recursion_20261018e.md),
after integrating concurrent commit 5447ae726c:
[marked geometries and selectors](marked_geometries_selectors_20261018o.md).
The latter distinguishes four charts in one partition from four different
partitions in an SL2(3) orbit. That distinction exposes a larger exact object:
there are **105 optimal four-chart partitions**, not just the seven parallel
coset partitions. A marked primal/dual map connects them to the older 105
Hamming-replacement edges.

## 1. Inheritance and board

Anchor: optimal joint-boundary representations. Niche: the 105-matching
incidence graph. Wildcard: physical translation versus addition of a linear
function and their central phase.

- Closest proved mechanism: kappa(C)=dim(C)-nu(C) from the companion; the
  incoming SL2(3) action preserves the Hamming code but can change its atlas.
- Hostile: preserving a code does not preserve a chosen partition; preserving
  a partition does not fix each word or its phase.
- Corrected near miss: the preceding Gray carry can be omitted without
  changing either the code or its four cells, while corrupting its inverse.
- Sidecar: the physical address origin, the pairing between primal addresses
  and dual slopes, and the complete affine action on component phases.

Board: named codewords; maximum charts; slope quotient; physical matchings;
atlas-direction geometry; central phase. This packet supplies maps among
these objects. No unit-distance graph is inferred from their cardinalities.

Write a word of C=H8+I2 as

\[
f(X)=s+\ell\cdot X\quad(X\in\mathbf F_2^3),\qquad (t,t),
\]

where s,t are phase bits and ell is a three-bit slope. Keep the standard
address basis: the three linear functions have masks 170,204,240.
Let E=span(255,768). The quotient C/E is the eight-element **dual** slope
space; the four words in one quotient fibre vary only s,t.

## 2. All maximum charts, and all optimal partitions

**Proposition A1 (PROVED).** There are exactly 28 maximum component charts
in C and exactly 105 optimal four-chart partitions. Every maximum chart is
the inverse image of an unordered pair of distinct slopes in C/E.

The previous theorem bounds a chart by eight words. Its three disjoint
generators must have weights 2,4,4: without the unique weight-two word
768, three nonzero supports would require at least twelve coordinates.
The remaining two tetrads partition the first eight coordinates and are
complementary affine hyperplanes. Their sum is 255. Hence every maximum
direction space contains E and is E+<L_d> for a nonzero linear function L_d.

A coset of such a space projects to an affine line in F2^3, which is just a
two-element set. Conversely the preimage of every two-element set is a
maximum chart with the stated disjoint half-cube/pair generators. This gives
binomial(8,2)=28 charts.

Any exact four-chart cover of the 32-word code must use four disjoint
eight-word charts. Its projection is therefore a perfect matching on the
eight slopes. Conversely every perfect matching gives such a partition.
Thus the exact number is 7!!=105, with a lossless matching-to-atlas decoder.
There is no assumption that the four charts use the same direction space.

The seven parallel atlases pair every slope ell with ell+d for a fixed
nonzero d. They are exactly the coset partitions used in the first proof.
The other 98 optimal atlases mix directions.

**Proposition A2 (PROVED).** The complete decomposition by number of
different pair differences is

\[
\boxed{105=7+42+56.}
\]

The XOR of the four pair differences is zero, since all eight slopes occur
once. Their multiplicities can therefore only be 4; 2+2; or 1+1+1+1.
There are seven parallel matchings. For two distinct directions d,e, the
plane they span has two affine cosets; in one coset use d and in the other
use e, or reverse the choices. This gives 2*binomial(7,2)=42.
The remaining 56 matchings have four different directions. Those four
directions form an affine plane avoiding zero. Each of the seven such
direction sets occurs eight times: GL(3,2) is transitive on those planes,
and acts bijectively on the matchings. The complete dictionary also checks
this uniform count.

## 3. The old 105 edges become an explicit atlas dictionary

The old [ten-label matching/Hamming construction](ten_label_matching_hamming_20261018n.md),
Section 3, fixes the distinguished pair and takes perfect matchings M of
the eight remaining **physical coordinates**. The even unions of the four
pairs form a three-dimensional subcode E_M. Exactly two H8 codes extend it.
They form one edge in the 30-vertex degree-seven replacement graph, the
point-plane incidence graph of PG(3,2).

Choose the retained address origin and dot-product pairing to identify a
physical address X with the slope X dot (-). Then:

\[
\text{optimal atlas}
\longleftrightarrow
\text{matching of eight slopes}
\longleftrightarrow
\text{marked matching of eight physical addresses}
\longleftrightarrow
\text{old Hamming-replacement edge}.
\]

The checker constructs every arrow and its inverse. It independently forms
the three self-dual extensions of each E_M, recognizes the two with minimum
weight four, and reconstructs the 30 vertices and 105 distinct edges.
No arbitrary matching of two size-105 lists is used.

There is a sharper preserved predicate: **the seven parallel atlases map
exactly to the seven edges incident to the original H8 vertex**.
Indeed E_M is contained in H8 iff every affine function has the same
pair-XOR on all four pairs. Self-duality turns this into orthogonality
against every even union of pairs. For pair {X,Y}, its XOR is ell dot (X+Y);
equality for every ell is equivalent to all pair sums X+Y being equal.
That is precisely a parallel matching.

The chosen primal/dual identification is not equivariant for the full
physical affine group. A physical translation X -> X+b changes s to
s+ell dot b and leaves ell unchanged. It fixes each maximum chart as a
set, whereas it can move a matching of physical coordinates. For example,

~~~
{01,23,46,57} -- translation by 4 --> {02,13,45,67}.
~~~

The linear part acts contragrediently on slopes. Retain the basis, origin
and duality; do not replace them by an unqualified natural bijection.

## 4. Seven parallel directions split into four vertices and three matchings

The incoming SL2(3) action gives 24 actual permutations of the eight physical
ports. It preserves H8. Its induced action on slope directions has kernel
{I,-I} and image A4. The seven nonzero directions split into orbits of
sizes four and three. The four-element orbit is exactly the incoming four
parallel atlases, authenticated as partitions, not just counted.

For convenient quotient representatives modulo E, the complete space is

~~~
{0,15,51,60,85,90,102,105}.
~~~

The four-element orbit is O={15,60,90,105}; the remaining three form
W minus zero, where W={0,51,85,102}. Thus O=15+W. There is an exact
equivariant K4 dictionary: its four vertices are O, and each nonzero w in W
pairs vertices a,b precisely when a XOR b=w.

| Direction | Perfect matching of the four atlas vertices |
|---|---|
| 51 | 15--60 and 90--105 |
| 85 | 15--90 and 60--105 |
| 102 | 15--105 and 60--90 |

This makes the old four-vertex/three-matching structure reappear as an
actual action on **parallel atlas directions**. It is distinct from the
four T,+,-,S states selecting cells inside one atlas. The two uses of four
have different domains, although both admit small marked controllers.

The central element -I fixes every slope but still swaps opposite physical
ports: it sends the word 15 to 240=15 XOR255. Thus even the faithful A4
description of atlas motion forgets a genuine component-phase action.
Keeping only the tetrahedron is insufficient to reconstruct each word.

## 5. A second literal mod-four mechanism: two involutions and a phase

Let T_b be physical address translation, S_a addition of a linear function,
and Z complement of the first eight bits:

\[
T_b f(X)=f(X+b),\qquad S_a f(X)=f(X)+a\cdot X,\qquad Zf=f+1.
\]

All preserve C, with t unchanged. Both T_b and S_a are involutions, but

\[
\boxed{T_bS_a=Z^{a\cdot b}S_aT_b,\qquad
       (T_bS_a)^2=Z^{a\cdot b}.}
\]

Thus their product has order four when a dot b=1. The apparent commuting
operations on the slope quotient have a nontrivial phase upstairs.
At a=b=(1,0,0), the zero word goes to 85 under T_b S_a and to 170 under
S_a T_b; the two differ by 255.

**Proposition A3 (PROVED).** The group generated by these operators consists
of exactly 128 distinct maps

\[
g(a,b,c):(s,\ell)\mapsto(s+\ell\cdot b+c,\ell+a).
\]

Composition, with the right map applied first, is

\[
g(a,b,c)g(a',b',c')
 =g(a+a',b+b',c+c'+b\cdot a').
\]

The formula proves closure. The slope displacement determines a; evaluating
the constant displacement at ell=0 and the three basis slopes determines
c,b, proving faithfulness. The centre is exactly {1,Z}: commuting with all
slope additions forces b=0, and with all translations forces a=0.
There are 56 elements of order four and 71 nonidentity involutions.
Adjoining the independent distinguished-pair flip doubles the group to 256.

This is the finite binary Heisenberg group, with the ordered central
cocycle b dot a'. The incoming 27-point construction uses a ternary
Heisenberg group with an oriented-area coordinate. The connection is the
central-extension mechanism, not an identification of the two groups.
In characteristic two one cannot simply copy its alternating-area formula:
its commutator 2*area would vanish. The ordered cross term above retains
the missing phase.

This group is also different from the **abelian** quaternary source group
Z4^2 x Z2 in the companion. There, carries preserve addition in an explicit
Gray codec. Here, an ordered interaction of two involutions produces an
order-four operation. Both require a central binary record; their operations
and commutativity properties remain distinct.

## 6. Consequence and checks

The immediate gain is a complete search domain: 105 optimal partitions,
with 7 parallel, 42 two-direction and 56 four-direction alternatives.
Changing the partition can change the component graphs and which physical
attachments are convenient, while preserving exactly the same 32 boundary
words. Every candidate still faces the majority obstruction: no fixed
two-list graph can realize C, regardless of which of the 105 atlases is used.
The conditional three-colour-list gate is already constructed. Forcing its
required lists and gluing the listed gates at fixed physical ports remain
open, as does the unrestricted geometric selector.

The principle suggested by this recovery is now testable: **change the
decoder or partition while preserving the exact target relation; retain the
phase that distinguishes the operations used to change it.** The
primal/dual translation hostile is its counterindication to forgetting
the marking. No universal “least resistance” theorem across number theory
is asserted.

~~~
python -B -X utf8 04-computation/experiments/hamming_matching_atlases_20261018e.py
python -B -O -X utf8 04-computation/experiments/hamming_matching_atlases_20261018e.py --verify
~~~

The [JSON dictionary](hamming_matching_atlases_20261018e.json) retains all
28 charts, 105 matching/atlas/edge records, the 30 Hamming vertices, all
24 incoming port permutations, their quotient actions and the K4 dictionary.
It records a SHA-256 of the incoming geometry receipt. There are 531,206
active generation checks; optimized replay adds one equality check.
The all-size chart theorem is proved in the companion; this finite packet
classifies its length-ten equality case.

Independent agents enumerated all 105 partitions and seven parallel ones,
checked every SL2(3) action on the 4+3 dictionary, verified the physical
translation hostile and the seven incident-edge correspondence, and audited
the central group law. Full composition is replayed on all 128^2 pairs of
operators and all 32 actual source words. The masks 56,105,128 here have
specific counted meanings; equal integers elsewhere do not supply additional
maps.
