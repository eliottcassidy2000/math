# What can be compressed, averaged, or transported?

2026-10-28. This is a targeted archive recovery and synthesis. It does not
claim a new universal theory of automata, Burnside rings, or conditional
expectation. The new concrete result is [M1--M5](THEORY.md). Proposed extensions
are OPEN. LRC(14) and the other major open targets receive no claimed solution.

## 1. The Borwein numbers identify a support mechanism

The recovered source is THM-4837
[`the-owners-borwein-ladder-first-nonvanishing-of-the-odd-sinc-cosine-transform-at-2-3-and-2pi-closed-form-first-failure.md`](../../../01-canon/theorems/THM-4837-the-owners-borwein-ladder-first-nonvanishing-of-the-odd-sinc-cosine-transform-at-2-3-and-2pi-closed-form-first-failure.md),
with its [corrected provenance account](../borwein_ladder_owner_numbers_vsc_twelve_20261028.md).
MISTAKE-692 matters: the thresholds and first-bite formulas are published
results, not discoveries of this session. Primary sources checked include
[Baillie--Borwein--Borwein 2008](https://www.tandfonline.com/doi/abs/10.1080/00029890.2008.11920606),
[Basel--Baillie 2016](https://arxiv.org/html/1510.03200), and
[Baillie's sinc-sum account](https://arxiv.org/abs/1105.3943).

Write O_m=sum_(j=0)^m 1/(2j+1). The first crossings of2,3,2pi occur at
m=7,56,40249, hence denominators15,113,80499. We independently recertify all
three crossings with directed decimal intervals and a rational alternating-
series bound for pi. No approximation to a tiny integral is needed.

For m>=1 and xi>=1, the inherited probabilistic formula is

    F_m(xi) = integral_R cos(xi x) product_(j=0)^m sinc(x/(2j+1)) dx
            = pi Pr(xi-1 <= U_1+...+U_m <= xi+1),

where U_j is uniform on [-1/(2j+1),1/(2j+1)]. Its density is positive in
the interior of its support. The window has positive overlap precisely when
O_m>xi. Equality gives zero mass. The m=0 endpoint exception is excluded.

**Elementary support principle, PROVED:** convolution of nonnegative
densities positive on interval interiors is positive on the interior of the
Minkowski sum of those intervals. Indeed every interior sum has an open set
of decompositions lying in both interiors, on which the convolution
integrand is positive. Induction proves the statement for any finite number
of factors. Hence a positivity query can discard the density's shape while
retaining support endpoints. Signed kernels or exact integral values require
more information; cancellation destroys this inference.

This gives a real, limited comparison to M1. Source objects are products of
sinc factors on one side and Hamilton presentations on the other. Their
reduced coordinates are respectively the real margin xi-O_m and the discrete
offset d plus polarity. Appending a factor subtracts its half-width from the
margin; singleton insertion changes d by one while testing polarity. The
preserved predicates are nonvanishing and blue, respectively. Lost data are
the convolution density and the interior relation. There is no proposed
isomorphism between the integrals and tournaments.

The useful research move is to prove why a response is a boundary-crossing
test before carrying its full history. M1 supplies that proof for blue;
M5 shows why the same shortcut fails for interior alignment queries. Long
agreement alone is insufficient:15,113,80499 are exact crossing locations,
and the padded seven-vertex pair supplies arbitrarily deep hidden alignment.

## 2. Oddness allows averaging; it need not allow choosing

There are several different oddness mechanisms in the archive.

1. Every tournament automorphism group has odd order: an involution would
   interchange a pair of vertices and reverse their arc. More generally every
   cycle of every automorphism has odd length.
2. The Hamilton-path number H(T) is odd by Redei's theorem. Aut(T) acts freely
   on ordered Hamilton paths, so the number M=H(T)/|Aut(T)| of literal path
   presentations of its class is odd as well.
3. An odd-order self-converse tournament has an anti-involution with one
   fixed vertex. Its anti-involutions form a single Aut(T)-orbit. This is
   the classical mechanism behind [Eplett's representation theorem](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/ED7BC2A1F6FF2D7F31547D771CE69B80/S000843950005743Xa.pdf/selfconverse_tournaments.pdf),
   recovered in the [consecutive-pair lift](../session89_consecutive_tournament_pairs_20261028/encoding/ODD_SELF_CONVERSE_LIFT.md).

For a nonempty finite fiber of odd size M, uniform averaging P=J/M is defined
over the rationals with odd denominators and descends modulo2. Since
J^2=MJ, its reduction is still an idempotent. This is an immediate instance
of the inherited [weighted-fiber theorem](../session88_symmetry_certificates_20261028/operators/TRANSFER.md),
not a new averaging theorem. Applied to item2, it explains why fresh path
reset has a valid mod2 projection at **every** vertex order. Applied to
item3, it averages anti-involutions without choosing one.

By contrast, the weighted root-switch operator satisfies Q^2=nQ: modulo2
it is idempotent at odd n and square-zero at even n. These statements concern
different fibers and different operators. The parity of the number of
vertices cannot substitute for the parity of a reset fiber.

An invariant averaging operator does not give an equivariant section.
A cyclic group acting regularly on three points has no invariant choice of
one point, despite the odd cardinality. The Paley7 example makes this concrete:
Aut has order21, there are seven anti-involutions, and each has centralizer
of order3. No anti-involution is fixed by the full group. Fixing the pole0
and anti-involution i->-i leaves i->2i, the order3 Frobenius phase. We replay
all these maps and H=189 exactly. A class bijection can therefore preserve
existence and counts while retaining only a centralizer of the full symmetry.

This also separates two older notions of oddness. THM-1560
[`the-halving-dictionary-and-the-mod-two-collapse.md`](../../../01-canon/theorems/THM-1560-the-halving-dictionary-and-the-mod-two-collapse.md)
shows that skew adjacency S reduces to J-I modulo2 for every tournament:
integer polynomials in S lose all orientation there. The conversion
A=(J-I+S)/2 cannot be used in characteristic2. Redei's odd path count is a
different phenomenon; it should not be inferred from that collapse.
THM-1380
[`free-Z2-but-no-odd-map.md`](../../../01-canon/theorems/THM-1380-free-Z2-but-no-odd-map.md)
is another useful warning: the free involution and the odd function in an
LRC analogy must concern the same action. Merely finding both separately
does not permit a Borsuk--Ulam argument.

## 3. Even graphs and trijections live at different levels of symmetry

Here an **even graph** means that every automorphism reverses an even number
of its edges under any orientation. It does not mean an Eulerian graph.
The equality of unlabelled tournament and even-graph counts is
[Royle--Praeger--Glasby--Freedman--Devillers' theorem](https://arxiv.org/abs/2204.01947).
The repo strengthens it in THM-4721
[`twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`](../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md):
signed even-subgraph counts give a Walsh basis over Z[1/2], and suitable
Sylow/prime-local refinements preserve more symmetry information.

This basis change is especially instructive beside the preceding paragraph.
Walsh inversion inverts2; odd-fiber averaging inverts odd numbers. They are
different arithmetic regimes. Neither licenses silently reducing the other
construction modulo2.

The inherited hostile to an unrestricted odd-group refinement is C15 acting
on10 vertices with orbit sizes3,5,1,1: there are32 tournament types but31
even-graph types admitting that action. See the exact
[C15 counterexample](../c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md).
Thus oddness alone does not make every symmetry-refined count agree.

For the user's graphs with exactly one fewer edge than vertices, use the
sharper theorem THM-4835
[`euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md`](../../../01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md).
Let D be an oriented carrier with n vertices and n-1 edges and A=Aut(D).
Its free completion bits always number C(n-1,2). An A-equivariant signed
pair-coordinate bijection to tournaments on an (n-1)-set exists exactly when

    chi(G^K)=1 for every subgroup K<=A,
    equivalently [V]=[E]+1 in the Burnside ring.         (inherited)

The identity subgroup tests only |V|-|E|=1. For a disconnected graph this
allows cycle rank=components-1, but says nothing about the other fixed-point
tests. Trees satisfy the full criterion, with a concrete map: a nonadjacent
vertex pair goes to the two end edges of its unique path. Functional
orientations give another explicit family, including a directed cycle plus
an isolated root. Arbitrary disconnected carriers need not qualify.

MISTAKE-689 prevents a second overreach: failure of this pair-coordinate
criterion does not automatically destroy weak equivariant bijections or
matching orbit counts. The archive has explicit separations of those tiers.
The cost order in M4 is the operational analogue: one marked port, one Boolean
guard and a selected reset are different interfaces, not interchangeable
measurements merely because each separates a chosen pair.

There is also a composition distinction. Strong components of a tournament
form an ordered sequence; components of an even graph form an unordered
multiset. THM-4721's EULER(connected-even)=INVERT(strong-tournament) identity
relates their counts without identifying these two composition operations.
The new blue offset is literally the difference of the two extreme singleton
strings after peeling; it uses the ordered-component structure. This connects
to THM-4836
[`coherent-carrier-sets-of-tournaments-and-the-pair-or-triple-law.md`](../../../01-canon/theorems/THM-4836-coherent-carrier-sets-of-tournaments-and-the-pair-or-triple-law.md),
whose path/cycle/star classification likewise detects extreme components.

## 4. A composition program with explicit boundaries

The cohesive proposal is to build a composition theory in layers:

| mathematical object retained | what it answers | obstruction to replacing it by less |
|---|---|---|
| a residual response such as B | every future blue verdict in M1 | change the selector language and the quotient changes |
| actual rooted relation types | structural attachment observations | endpoint decoder M2 makes these necessary |
| compatible maps or frames | joint transport and phase | separate pair orbits need not have a common witness |
| symmetry marks and weights | equivariant quotients and repeated forgetting | equal scalar counts need not lift to such maps |
| all ordered path runs and multiplicities | revisiting quotient-walk assembly | the inherited six-vertex two-endpoint counterexample |

These are typed interfaces, not five equivalent descriptions of every object.
The [LCA/frame theorem N1](../strong_attachment_tail_compiler_20261028/THEORY.md)
gives one exact interaction: on free frames, tree-edge transports determine
all joint alignments by path multiplication. For nonfree ports, multiplying
double cosets can invent incompatible intermediate witnesses. Frobenius
phases belong to the actual map layer. On six distinct normal-basis ports
only shifts0,2,4 can be tournament automorphisms; on eight only0 can. The
inherited F64/F256 checks remain in force. There is no claim that the new
blue polarity bit is the same bit as a Frobenius orientation character.

The 81/27 control is replayed:67 common labelled cyclic triples,462915 common
Hamilton paths, but176/512 binary attachment orbits and27/0 transporters
from0 to4. These tournaments already differ in q. They are a hostile test
of replacing relation types by cyclic triples or scalar counts, not evidence
that the blue coordinate alone distinguishes them.

The square-class/Hamming-code and Vitali lanes fit the same boundary only
after specifying their maps. Full binary syndromes preserve linear
dependencies; a blue observer returns just the truth of a conjunction of
reflection equations. Those are different outputs. The previous packet's
[syndrome witness](../singleton_pair_reconstruction_20261028/CONNECTIONS.md)
has different nonzero syndromes but equal false blue verdicts. In M1,
geometry limits which conjunctions can be realized, allowing compression;
named guards recover distinctions that this restricted geometry cannot.
M5 now gives an explicit family for the growing distance to a hidden alignment.
This is stronger evidence than a suggestive similarity of integer sequences.

## Incoming work: a quantitative alignment-memory mechanism

Before close-out, Session90's
[iterated triangle family](../session90_even_graphs_odd_descent_20261029/carriers/ITERATED_TRIANGLE_FAMILY.md)
and [transport theorem](../session90_even_graphs_odd_descent_20261029/carriers/TRANSPORT_AND_HOLONOMY.md)
arrived. They sharpen the proposed composition program and are credited to
that packet, not claimed as new here. We read the proofs and reran its exact
family checker through81 vertices; this replay is not an independent algorithm.

For W_k=C3[W_(k-1),W_(k-1),W_(k-1)] of order N=3^k, let C be the centralizer
of digitwise negation. Intrinsic top-block recovery proves the full wreath
group, and conjugating the two exchanged children forces their actions to
match. The resulting formulas are

    |Aut(W_k)| = 3^((N-1)/2) = N |C|^2,
    |C| = 3^((N-2k-1)/4).

Exactly |C| anti-involutions fix any one named vertex. At N=27 this is243;
specifying the actual reversal mark needs8 further bits. Asymptotically the
requirement is (N/4)log2(3)+O(log N) bits. At N=81 it already needs29 bits,
whereas the blue residual B has316 states and needs9 bits. These answer
different questions: an exact named reversal frame versus the restricted
future blue verdict. The small B does not compress every reversal frame.

This supplies a quantitative partner to M5. M5 hides alignment beyond every
fixed boundary width; the wreath family counts how many reversal alignments
remain even after a pole is supplied. It provides a concrete tournament
mechanism for the growing-memory theme, without identifying it with every
Vitali or arithmetic construction.

The same incoming packet supplies an exact way to avoid one search: for
supplied anti-involutions s,t, the odd permutation square root g(t,s)=sqrt(ts)
conjugates s to t. Compute it cyclewise. No search for a correcting
automorphism is necessary once the initial isomorphism and marks are known.
But g(u,t)g(t,s) can differ from g(u,s) by a nontrivial centralizer element;
the9-vertex wreath example realizes this error. Thus a local canonical
correction does not automatically give globally consistent transport.
A shared base frame makes the path products consistent, as in the LCA theorem.

The companion [coprime fixed-orbit proof](../session90_even_graphs_odd_descent_20261029/oddness/FIXED_ORBIT_DESCENT.md)
also clarifies the oddness comparison: for compatible p-group A and coprime
G actions, X^A/C_G(A) maps bijectively to (X/G)^A. An A-stable orbit, and the
transporter between two fixed points in it, both have size prime to p and
hence contain A-fixed elements. The classical descent argument applies to
the joint relational object; applying it separately to marginals discards
alignment. This is a precise sufficient condition for recovering a symmetry
choice, with the missing joint-data hypothesis visible.

## 5. Concrete open problems and next decisive tests

**OPEN A: least endpoint-decoding cost.** Four insertions suffice for arbitrary
presentations; two suffice for strong cores. Can one fixed word with at most
three insertions recover ordered endpoints for every presentation? Enumerate
the finitely many reduced words at that budget and find same-output/different-
endpoint witnesses. A finite failure rules out a word, while success on a
census alone does not prove it universal.

**OPEN B: moving versus fixed selectors.** M5 rules out every fixed boundary
width. Define a cursor that can move one position at a time while retaining
its vertex identity through complement. Determine its least extra state and
whether a single cursor with a local legality query reconstructs every
presentation. A cursor that is freshly resampled each step is a different
operation. The9-vertex buffered pair is the first required hostile probe.

**OPEN C: minimal carrier for general atom composition.** Replace singleton
contexts by a declared finite library of quotient-walk attachments. Determine
when a finite collection of homogeneous boundary defects still suffices and
when full joint ordered run data is forced. Start with the smallest prime
quotient that revisits a child; require the six-vertex counterexample and
81/27 control to survive. M3 settles only the singleton sublanguage.

**OPEN D: enumerating the behavioral quotient itself.** Count K_min classes
at8 using rooted source/partner types and B, and compare direct construction
with finite-horizon refinement. M3 gives a certified stopping horizon13,
but direct K construction avoids materializing X_13. Measure both storage
and canonicalization work. This is a new count channel, not A000568: at7
the respective counts are32272 and456. Existing A000568350-prefix and complete
order10 class banks remain inherited; this session does not extend them.

These targets turn the broad synthesis into falsifiable steps. A future
transfer to Collatz, arithmetic sieves, coloring or Jacobian problems must
name the represented object, allowed operations, preserved predicate and
missing sidecar. No universal scalar potential, canonical symmetry choice,
or flagship breakthrough follows from the correspondences established here.
