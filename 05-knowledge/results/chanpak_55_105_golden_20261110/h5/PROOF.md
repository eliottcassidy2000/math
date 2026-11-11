# A fixed-core attachment compiler and the sharp parent floor at 173

Status: PROVED IN THIS NOTE for the formulas and completeness of the inverse;
FINITE-EXACT for its evaluated target window 1..1001 and the explicit
certificates. Not independently canon-audited. Chan–Pak remains OPEN.

The complete one-vertex images of ALL tournaments of counts 1, 3, or 5
miss exactly **7, 21, 173, 335, 663** among the odds through 1001.
In particular every tournament of count 173 has every vertex-deletion
count at least 9. The bound is sharp: an explicit eight-vertex tournament
below has a deletion of count 9. This extends the previous sharp parent
floor 5 at target 55. Target 105 has a parent of count 3, but has no
parent of count 5 at any order.

## Inheritance and the scope of the compiler

Closest mechanism: the complete H=3 two-flank compiler in
`decomposition16_odd_spectrum_20261110/attachment/PROOF.md`. Its OCF path
expansion generalizes to an arbitrary fixed core. The seven-coordinate
ordinal law in `session250_sixteen_20261110/constructive/PROOF.md` already
retains the endpoint and two-cover sums used here; the present formula
compiles arbitrary binary flank words into two integers and their parity.
The signed cancellation is inherited rooted path reciprocity, proved in
`session216_overnight_constructions_20261107/carriers/reciprocity/UNIVERSAL_PAIR_RECIPROCITY.md`,
equations (2),(3). These are recovered mechanisms, not new universal laws.

Incoming Session251 also pointed back to the older
[Session197 one-core binary inverse](https://github.com/eliottcassidy2000/math/blob/8e097d715769dd942ab21807a5a5ece39b66bec2/session197/path-polynomials/ONE_CORE_BINARY_INVERSE.md),
read directly from the pinned Git object. Its equation(8) is exactly the
special case of (2) where O=V(K) minus one named vertex b. In that source's
notation (H,E,S,Q,L,C), our (h,q,e,s,u,v,t) becomes

    (H,C,L,C-H-Q,S,H-E,H-S+L),

and our flank integers A,B are its q,t. Substitution reproduces its full
formula, including the parity corrections. The general idea of a physical
two-integer compiler and a bounded target inverse is therefore inherited.
The present work supplies arbitrary core cuts, the complete H5-parent
image, the general reentry correction, and the sharp173 parent floor.
The older50-core bank and its462661 hostile belong to that specified
stationary language, not the complete low-parent image evaluated here.

The physical word conventions agree too: its prefix O/Z letters are our
left bits1/0, and its suffix Z indicators read least-significant first
are exactly our complemented, reversed right word. The retained root is
x, with core out-set all vertices except b. The checker verifies the
coefficient map at every named C3/Q4 port and the scalar equality for
all A,B in0..7. Thus this attribution is a literal restriction, not just
an analogy between recurrences.

Canonical hostile example: MISTAKE-756 in `01-canon/MISTAKES.md` disproves
collapsing a uniformly contacted module to its scalar H. A Hamilton path
can visit that module in two separated pieces. We retain its ordered
two-path-cover count. Least-used sidecar: the opposite endpoint sums in
the signed path reciprocity identity. Concept board: parent count;
named core cut; ordered path covers; binary flank words; module reentry;
target inverse. No external literature priority claim is made.

H counts ordered directed Hamilton paths, with H(empty)=1. An ordinal
join orients every cross arc from the earlier factor to the later one.
This note concerns one-vertex attachment, not path-complement edges of
the tiling metagraph. The attached vertex and core vertices are named.

## 1. General fixed-core formula

Fix a nonempty tournament K and a specified out-set O of a new vertex x
on K; put I=V(K) minus O. Start with

    TT_a -> K -> TT_b.

Encode x's incidences on the left flank by bits 1 iff x beats the old
vertex, in dominance order, and let their binary integer be A. Let B
be the binary integer obtained by complementing and reversing the
right-flank bits. Actual lengths a,b and all names remain sidecar data:
leading zero bits disappear from A,B but can be restored. Write
epsilon=A mod 2 and delta=B mod 2; empty flanks have integer zero.

Use the inherited one-apex count

    g(0)=g(1)=1,
    g(2r+1)=g(r),
    g(2r)=g(r)+r+(r mod 2)  for r>=1.                 (1)

This is the Hamilton count of a transitive source with the indicated
binary attachment word, with leading zeros immaterial. See THM-4875,
`01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`,
and the previous attachment proof for the derivation and conventions.

Define the following finite data on K:

* h=H(K), and q=Q_K(2), the number of ordered pairs of covering directed
  paths, allowing either path empty.
* e counts those ordered covers whose first path is nonempty and ends
  in I; the second path can be empty.
* s counts those ordered covers whose second path is nonempty and starts
  in O; the first path can be empty.
* u counts Hamilton paths of K starting in I; v counts those ending in O.
* t=H(K+x), with x's out-set exactly O.

**Theorem H5.1 (fixed-core compiler).** The resulting child has count

    F_(K,O)(A,B)
      = q A B + e A + s B + h(g(A)+g(B)-2)
        + t + u epsilon + v delta.                   (2)

Every term displayed after the first three is nonnegative: g>=1, t>=1,
and u,v>=0. This exact formula includes uniform core cuts O=empty/full.
It does not presume that K is a contiguous segment of a child path.

Proof. For a directed nonempty path p in K give it residual weight
H(K minus V(p)). Let W(z) sum these weights times z^|p| over all p;
E(z) restricts to paths ending in I; S(z) to paths starting in O;
D(z) imposes both restrictions. The one-variable z records path order.

The left-flank polynomial L(z) counts nonempty increasing paths starting
at an x-out-neighbor. The right polynomial R(z) counts increasing paths
ending at an x-in-neighbor. Directly summing the geometric choices gives

    L(1)=A, L(-1)=-epsilon, R(1)=B, R(-1)=-delta.

For example L(1)=sum_i bit_i 2^(a-1-i), whereas evaluation at -1 cancels
every first endpoint except the last vertex.

By THM-002, `01-canon/theorems/THM-002-ocf.md`, an odd directed cycle
containing x is x followed by an even-order path from an x-out-neighbor
to an x-in-neighbor in the deletion. The other odd cycles avoid that
path and contribute exactly its residual Hamilton weight. Thus the
increment in H is twice the total residual weight of such even paths.
For a polynomial P this twice-even contribution is P(1)+P(-1).

Paths avoiding K contribute

    h(g(A)-1)+h(g(B)-1)+h(AB+epsilon delta).

The first two terms lie wholly within one flank, and the third crosses
directly from left to right. Paths meeting K contribute the evaluations
at both signs of

    W L R + E L + S R + D.

The inherited squarefree identities P(-X)P(X)=1 and
S_i(-X)P(X)=-E_i(X), and the converse identity, imply

    W(1)=q-h, W(-1)=-h,
    E(1)=e, E(-1)=-u, S(1)=s, S(-1)=-v.

Here taking the squarefree top coefficient multiplies each chosen path
by H of its complement, exactly the residual weight above. Finally
D(1)+D(-1)=t-h by the same OCF expansion on K alone. Summing gives
(2); in particular the epsilon*delta coefficient cancels. QED.

**Corollary H5.2 (exact uniform reentry).** Suppose O is empty, and let
Q be the quotient obtained by replacing K by a single vertex with the
same uniform incidences. Then

    H(child)-h H(Q) = (q-2h) A(B+1).                 (3)

For O full the answer is (q-2h) B(A+1). The factor q-2h is exactly
the number P_2(K) of ordered covers by two nonempty paths.

Proof. For O empty, (e,s,u,v,t)=(q-h,0,h,0,h). Apply (2), and subtract
h times the singleton-core formula with (h,q)=(1,2). The full case is
its converse. QED. For the directed triangle P_2=6, recovering the
previous correction. For the strong four-vertex core P_2=22. Uniform
contact therefore does not even force divisibility by H(K).

## 2. A terminating target inverse for any fixed core

For a proper nonempty O, both e and s are positive: a singleton path
at any vertex of the required set contributes H(K-vtx)>0. Hence every
solution F_(K,O)(A,B)<=N satisfies

    0<=A<=(N-t)/e,
    0<=B<=(N-t-eA)/(qA+s),                           (4)

with integer floors. Evaluate (2) on this finite rectangle with varying
row lengths; retain exactly the target values. This is complete at every
order, not a vertex-order truncation.

For O empty, e=q-h>0, s=0, t=h. If A>=1 use

    1<=A<=(N-h)/(q-h),
    0<=B<=(N-h-(q-h)A)/(qA).                         (5)

If A=0 the count is h*g(B). Repeatedly strip trailing 1 bits by (1)
until B is zero or even B=2r. In the even case g(2r)>=r+1, so
r<=floor(N/h)-1. This gives a finite list of numerical witnesses for
all possible h*g(B)<=N. The full cut is the swapped construction.
Arbitrarily long all-1 words need this normalization; bounding their
unreduced integer directly would be an invalid inverse.

The inverse is for numerical existence. To retain every marked class
or original named-port configuration one must also retain widths and
the stripped-word history. For an actual construction keep the named
core, selected mask, both words, and a fresh name for x. These data
reconstruct every arc. They are not a closed replacement for the joint
17-coordinate observer or its F49/F343 and 81/27 coordinates.

## 3. Why this covers every H=5 parent

**Lemma H5.3.** Every tournament with H=5 is uniquely, up to its ordered
prefix/suffix lengths, TT_a -> Q4 -> TT_b, where Q4 is the unique strong
four-vertex tournament.

Proof. Hamilton counts multiply over the ordered strong components.
Every nontrivial strong component has H>=3, while a transitive singleton
has H=1. Primality of 5 leaves exactly one nontrivial strong component C,
of count 5. By OCF, count 5 permits exactly two directed odd cycles in
total; in particular there are at most two directed triangles.

Take a triangle D in C. If every exterior vertex contacted D uniformly,
partition them into L beating D and R losing to D. Strongness forces
an arc r->l with r in R,l in L, producing one new triangle l,d,r for
each of the three d in D, impossible. Thus some vertex contacts D in
both directions, giving a strong four-vertex subtournament Q. It has
two triangles and is isomorphic to Q4 (the one- and two-element cuts
on a directed triangle are isomorphic under rotation and relabelling).

Any additional vertex mixed toward Q creates another triangle: the
proper cut into its out- and in-neighbors has an arc from the former
to the latter because Q is strong. Therefore all additional vertices
contact Q uniformly. If any existed, strongness would again force an
R-to-L arc, producing further triangles. Consequently C=Q4. All other
strong components are singletons, proving the form. Conversely Q4 has
exactly two triangles and no longer odd cycle, so OCF gives H=5;
ordinal transitive padding preserves H. QED.

For explicit named calculations choose Q4 adjacency row masks

    (4,1,10,3),                                       (6)

on vertices 0,1,2,3, bit j in row i meaning i->j. It has h=5,q=32.
Its Hamilton start/end vectors are (1,1,1,2), (1,2,1,1), and its
two-cover start/end vectors are (5,6,7,9), (7,9,5,6).
The complete sixteen-cut coefficient table for (2) is:

| outmask O | e | s | u | v | t |
|---:|---:|---:|---:|---:|---:|
|0|27|0|5|0|5|
|1|20|5|4|1|9|
|2|18|6|4|2|9|
|3|11|11|3|3|9|
|4|22|7|4|1|11|
|5|15|12|3|2|11|
|6|13|13|3|3|13|
|7|6|18|2|4|9|
|8|21|9|3|1|13|
|9|14|14|2|2|15|
|10|12|15|2|3|11|
|11|5|20|1|4|9|
|12|16|16|2|2|15|
|13|9|21|1|3|13|
|14|7|22|1|4|11|
|15|0|27|0|5|5|

Every named cut is retained, even when some scalar coefficients coincide.
For example masks 9 and 12 have the same t=15 and u=v=2 but different
e=s (14 versus 16). Scalar child count and ordinary endpoint sums alone
therefore still fail to predict subsequent flank attachments.

## 4. The evaluated inverse and the first new obstruction

Apply the inverse to the singleton core, C3, and Q4. Their source
families are respectively all nonempty H=1, all H=3, and all H=5 tournaments.
The empty H=1 source adds only target1, which is already present.
The first two classifications are inherited; Lemma H5.3 gives the third.
For target odds through 1001 the three image sizes are 279,362,473;
their union misses exactly

    7, 21, 173, 335, 663.                             (7)

This is FINITE-EXACT target arithmetic using the all-order inverse.
It makes no claim about the first missed target beyond 1001 or a larger
finite parent bank. It refutes universality of the parent bank {1,3,5}.

The new 55 certificate is Q4, outmask 2, A=2,B=0, giving rows

    (62,124,80,4,104,76,9), H=55,
    deletion counts (9,31,9,27,15,23,5).              (8)

For target 105 the H5 inverse is empty at every order. A C3-parent
certificate, outmask1,A=2,B=3, gives

    (126,252,104,240,228,192,128,5), H=105,
    deletion counts (15,59,33,33,39,51,51,3).         (9)

Thus parent images are not nested by parent magnitude or divisibility:
the value 5 divides 105 but cannot be its one-vertex parent count.
Neither a product factorization nor a scalar spectrum coverage axiom
specifies an allowed attachment.

For 173, the empty H1/H3/H5 inverses, oddness, and the inherited exclusion
of 7 (THM-1370, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`)
show that every vertex deletion of every H=173 tournament has H>=9.
This lower bound is attained. Use the H=9 core (2,4,9,19,7), cutmask2,
A=B=1, producing

    (126,196,72,210,230,206,128,5), H=173,
    deletion counts (27,71,51,27,51,71,77,9).         (10)

The output8! Hamilton words and all eight deletion7! word counts are
recounted literally, independently of the compiler and subset DP.
This gives a sharp next step of the minimal-parent hierarchy:
55 first escapes {1,3}; 173 first escapes {1,3,5}, apart from7,21.
No infinite ascent of those floors is asserted.

As a reentry control, Q4 with A=B=1,O=empty has H=69. Its collapsed
quotient has H=5, so the false scalar module rule predicts 25. Formula
(3) accounts for the difference 44=22*1*(1+1) exactly. This also
explains why restricting uniform H5 cuts to multiples of5 would have
invalidated the target inverse.

## 5. Constructive consequence and next missing ingredient

Equation (2) is an exact instruction for generating physical children
from named core cuts and binary words. It offers a finite inverse for
each target and a reusable compiler for larger cores. The H9 core in
(10) is one such next instruction, not a new census assumption.
However a target may require a larger parent count, as 173 proves.
Increasing a finite bank or assuming target coverage does not show
closure under 3H+2/3H+4 with the required joint ports.

Source-to-target map: a fixed core and two ordered contact words map to
one attached tournament. Preserved predicate: literal Hamilton count,
including all reentries. Lost under numerical compression: zero-padding
names, full contact history, correlations with other marked ports.
Required sidecar: core labels, cut mask, flank words and widths, plus
joint attachment data for ports carried into a later step. Cheapest
decisive test: masks9/12 already separate equal t and u,v through e,s.
The squarefree reciprocity cancellation explains why two binary integers
are enough for one fixed-core attachment; it does not claim that finitely
many such integer states classify arbitrary descendants.

## Reproduction

From repository root:

    python3 05-knowledge/results/chanpak_55_105_golden_20261110/h5/check.py
    python3 -O 05-knowledge/results/chanpak_55_105_golden_20261110/h5/check.py

`check.json` and `check.optimized.json` must agree byte for byte. Checks
use exceptions, not Python assert. The explicit universe is every named
core cut for singleton/C3/Q4 and all flank widths with a+b<=3; 1274
subset-DP comparisons, of which314 also use literal permutations. The
64 labelled order4 tournaments independently give one H5 class. Every
one of the1114 selected inverse witnesses through1001 is graph-recounted.
The target bounds inspect9760,6318,4274 candidate parameter triples for
the respective cores; this is finite arithmetic, not graph enumeration
at a guessed order cutoff. The173 certificate has a separate literal
permutation check. Source file hashes are in `SOURCE_PINS.json`.
