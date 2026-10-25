# Pair palettes, Fano completions, and the exact missing obstruction

**Status:** PROVED elementary transfer statements below; FINITE-EXACT controls.
This is a repository bridge application, not a historical-priority claim.
No six-colour unit-distance equality gadget, seven-chromatic unit graph, or
plane upper bound is supplied. Surface classifications are inherited from
the audited Session60 packet, not reproved here.

## 0. Inheritance

The closest mechanisms are the Oct 11 fixed-ternary parity compiler and
complete boundary relation, the genus-two cuff/seam codec, THM-4614's
six-palette doubling options, and Session60's thirty simplex deck kernels.
The hostile is a locally valid lift with incompatible global seams. The
corrected near miss is identifying a pairwise shadow, one character, or
a number of labels with a complete colouring relation. The least-used
sidecars are the *ordered colour pairing*, its endpoint frame, the full
joint boundary relation, and the actual unit geometry.

Pinned files are listed in `SOURCE_PINS.json`. In particular:

* `plane_cycles_pants_and_six_seven_20261011b.md` is the prior theme map.
* `plane_six_reframes_parity_20261011b.md` proves the fixed-label lift.
* `spindle_tournament_interfaces_20261011g.md` gives an induced geometric
  three-port counterexample and the exact six-equality/seven-witness equivalence.
* `THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md`
  retains 15 pair matchings and 90 ordered-pair options; it is a scoped
  extension theorem, not a whole-field three-colour theorem.
* Session60 `seven-lines/RESULT.md` gives 30 chosen quotient kernels and
  16/12/2 surface types; `recovered-bridges/THEOREM.md` gives the marked
  genus-two two-torsion and Richelot bridge.

The live concept board is **six-colour pair frames / exact boundary support /
Fano completion / deck character / graph cohomology / unit geometry**.

## 1. The precise cohomological lift for six colours

For a finite graph G and q:V(G)->{0,1,2}, retain the graph H_q with exactly
the edges whose endpoint q-labels agree. Regard a graph as a one-dimensional
cell complex. Its binary cochain complex is

    C^0(H_q;F2) --delta--> C^1(H_q;F2),
    (delta b)(uv)=b(u)+b(v).

There are no two-cells, so every one-cochain is a cocycle. Let 1_q be one
on every edge. Then

    c(v)=2q(v)+b(v) is proper  <=> delta b=1_q
                              <=> [1_q]=0 in H^1(H_q;F2).          (1)

Consequently G is six-colourable iff there EXISTS q satisfying (1).
The cohomology group itself changes with q. This is not a single class on
the original graph whose vanishing decides six-colourability. When nonempty,
the lift set is a torsor for H^0(H_q;F2), with 2^c lifts for c components.
A spanning forest computes either a lift or an actual odd-cycle obstruction.

More generally edge gains a lift iff every cycle has total gain zero.
The associated double cover has vertices (v,b) and edges
(u,b)--(v,b+a(uv)). A graph-theoretic section is exactly a potential b.
For connected H the cover is disconnected iff [a]=0, in which case it has
two components. Crucially the all-negative cover a=1 is ALWAYS bipartite,
partitioned by the second coordinate. The triangle's cover is C6, but the
triangle has no section. Thus observing bipartiteness of a sign cover is
vacuous for the original bit-lift problem; the marked section is indispensable.

This is a genuine common lifting formalism with deck/pants constructions.
It does not identify H_q with the branch complement, a K3 surface, or a
unit-distance graph extracted from a sign cube.

## 2. Fifteen pairings really are the shared object

Let P be six named palette colours and choose a perfect matching
M={{a1,b1},{a2,b2},{a3,b3}}. There are 15 such matchings. The same marked
pairing has three exact incarnations:

1. It partitions six colours into three two-colour channels. Ordering the
   pairs gives the 15*6=90 options in THM-4614's doubling lemma.
2. On the genus-two model

       W={x in F2^P: sum x=0}/<1>,
       B([x],[y])=x dot y,

   the three pair classes are the nonzero elements of the Lagrangian
   L_M=<[a1+b1],[a2+b2]>. Matching-to-Lagrangian is a bijection between the
   15 matchings and all 15 Lagrangians of this four-dimensional symplectic
   space. Session60's actual six-branch curve supplies its Jacobian meaning.
3. Add a seventh named point 0. A labelled Fano plane in which its three
   lines through 0 are {0,a_j,b_j} is a Fano completion of M.

**Theorem (missing completion bit).** Every M has exactly two Fano
completions. There is no S6-equivariant choice of one completion, even if
the three pairs are ordered.

**Proof.** Choose an endpoint bit on each pair. The four Fano lines not
through 0 are exactly the four transversal triples whose endpoint bits
have sum eta, for eta=0 or 1. These two line systems are Fano: choose
nonzero independent a,b,c in F2^3, use central a and pairs
{b,b+a},{c,c+a},{b+c,b+c+a}; the eta=0 system is additive incidence, and
swapping an endpoint gives eta=1. All off-central lines are transversals,
and the unique-line axiom forces this parity rule, proving completeness.
Swapping the two endpoints of one pair fixes M and its ordering but
exchanges eta. Therefore an equivariant selector cannot exist. QED.

The parity bit is not canonically 0 or 1: it is a torsor until endpoint
frames are chosen. Neither ordering the pairs nor identifying a matching
with a Lagrangian repairs this missing bit. This gives an exact connection
from the old six-palette interface to the Session60 30=2*15 kernels, while
explaining the data that fails to transfer automatically.

For completeness, the Lagrangian bijection is also elementary: every
nonzero W-class has a unique two-element representative. Distinct such
classes are orthogonal exactly when those pairs are disjoint. Their sum
is represented by the complementary third pair. Thus each Lagrangian's
three nonzero elements recover one and only one perfect matching.

## 3. A transfer theorem and a two-bit obstruction that it forgets

Fix an ORDERED matching and endpoint frames. Its independent pair reversals
form A=F2^3. They act on the two completions through

    p:A->F2,   p(x1,x2,x3)=x1+x2+x3.                            (2)

Let Gamma be any finite gluing graph and prescribe frame transport
omega in C^1(Gamma;A). This setup keeps pair labels fixed on every seam;
it does not silently allow permutations of the three channels.

**Theorem (completion-shadow lifting).** Compatible full pair frames exist
iff [omega]=0 in H^1(Gamma;A). Compatible Fano-completion choices exist iff
[p omega]=0 in H^1(Gamma;F2). The natural sequence is exact:

    0 -> H^1(Gamma;ker p) -> H^1(Gamma;A)
      -> H^1(Gamma;F2) -> 0,                                (3)

and the forgotten obstruction space has dimension 2*b1(Gamma).

**Proof.** A change of vertex frame adds the coboundary delta f. It changes
the completion choice by p(f), so the projected equation is exactly the
image of the original one. The constant-coefficient short exact sequence
0->ker p->A->F2->0 splits linearly, proving (3). Since ker p has dimension
two and graph H1 has dimension b1, the dimension follows. QED.

**Hostile.** On a triangle put gain (1,1,0) on every edge. Each projected
gain is zero, so completions glue, but the cycle gain is (1,1,0), so full
pair frames do not. The first two binary channels each have an odd-cycle
obstruction. As a colouring example these can be two disjoint triangles
with their two fixed ternary labels: that assignment fails, although the
graph itself is three-colourable. The existential q in (1) is not removed.

If pair permutations are allowed, A becomes a local coefficient system;
the elementary dimension 2*b1 formula is not asserted without that change
of setup. A single scalar obstruction cannot replace all three channels.

### Exact but marked connection to the six-root canonical character

The six-root projective deck is G6=F2^6/<1> and its unique canonical
two-form character is chi_all([s])=(-1)^(sum s_i). Choose one endpoint
r_j in each matched pair. Then

    i:A -> G6,  i(x)=[sum_j x_j e_(r_j)]

is injective (a word supported on three chosen endpoints cannot be the
all-six word), and chi_all(i(x))=(-1)^p(x). Thus (2) is exactly the
restriction of the canonical deck character under this specified map.

This map identifies binary group elements, NOT geometric operations:
swapping two palette labels is not flipping a root sign. It uses endpoint
choices, preserves their parity character, and discards the two-dimensional
even-flip subgroup when only that character is observed. It transports the
linear obstruction (3), not the proper-colouring predicate or unit metric.

## 4. The seven canonical characters: pairwise blind, cubic exact

Now take Session60 G7=F2^7/<1>. Because seven is odd, every class has a
canonical even representative

    s=v+(sum_i v_i)1 in E7.

The seven canonical characters are kappa_i([v])=(-1)^(sum_(j!=i) v_j)
=(-1)^s_i. A Fano quotient kernel K, represented in E7, is a simplex code

    S={ (a dot v_i)_(i=0)^6 : a in F2^3 },

where v_i bijectively enumerate the seven nonzero vectors. All 30 such
codes have eight words, with nonzero weights exactly four.

**Theorem (third-order recovery).** With uniform averaging over K,

    E_K kappa_i = 0,
    E_K kappa_i kappa_j = 0                    (i != j),
    E_K kappa_i kappa_j kappa_k
       = 1 if v_i+v_j+v_k=0, and 0 otherwise.              (4)

Consequently every one- and two-coordinate distribution of these seven
characters is IDENTICAL for all 30 kernels. All pair projections have four
states, each with two code lifts. The third-order data recovers the full
Fano incidence. For the actual hourglass crossings T=012,034,056,

    d = sum_T E_K product_(i in T) kappa_i.                (5)

Thus third-order canonical-character correlations recover the inherited
d=0,1,3 stratification (16,12,2 kernels), while every statistic determined
by first and second marginals is blind to that geometric distinction.

**Proof.** Average (-1)^(a dot w) over all a in F2^3. It is one when w=0,
zero otherwise. Distinct nonzero labels never sum to zero in orders one
or two; three sum to zero exactly on a Fano line. Equation (5) counts
precisely the dependent ACTUAL crossing triples. QED.

**Projective-character guard.** A three-subset is not itself a character
of G7. The product of its three canonical characters has support the
COMPLEMENTARY FOUR-subset, because (1+e_i)+(1+e_j)+(1+e_k)
=1+e_i+e_j+e_k. Formula (4) is therefore a legitimate weight-four
projective character, not an illegal odd-subset sign observable.

Nothing here calls the seven v_i vertex colours. Their role is inertia
labels of a chosen quotient. The geometric classification in Session60
also requires actual branch crossings, singularity resolutions and its
canonical divisor calculation; a pairwise/cubic statistic alone on an
unrealized seven-set is not a surface classification.

## 5. Boundary relations: where the transfer must stop

For fixed q and a consistent bit lift, the exact boundary relation is an
affine system of pair-parity equations, so named pair equations plus the
internal inconsistency flag and path witnesses suffice. This is the
inherited compiler, not a newly discovered compression. Once q or its
interior choices vary, the union of these relations need not itself be
affine or determined by pairwise projections.

The inherited induced unit-distance counterexample uses alpha=(1+i sqrt3)/2
and vertices (0,1,alpha,-alpha,2,2alpha). Its graph is an inner triangle
with one pendant port at each corner. In palette three its port relation
is exactly NOT-ALL-EQUAL: 24 assignments, with two lifts each. Every pair
projection is unrestricted, exactly as for three independent ports. The
script independently rebuilds all distances using norm(x+y alpha)=x²+xy+y²
and checks all 3^6 colourings. This hostile is at palette THREE, not six.

For a generic six-palette example, join K4 to three independent ports.
Its relation is exactly 'at most two boundary colours': 96 assignments,
and again every pair projection is unrestricted. Gluing a port triangle
gives K7 and no six-colouring. This abstract witness contains K4, so has
NO injective unit-distance realization in the Euclidean plane. It proves
failure of a general interface compression, not a Hadwiger--Nelson bound.

The code example in §4 and the boundary example share a rigorously stated
loss mechanism—projection of a joint relation—not an identification of
their constraints. Storing the full joint relation restores the boundary
predicate; retaining all three channels restores §3's vector obstruction.

Finally, a nonvacuous six-equality gadget with UNIT-separated omitted-edge
ports is already equivalent to the desired seven-chromatic unit witness:
adding that edge destroys six-colourability and at most one new colour
suffices; conversely, delete an edge from an edge-minimal non-six-colourable
unit graph. Relabelled spindle/diamond gadgets do not supply it. The Moser
field's inherited four-colouring separates all omitted unit pairs with
six colours, irrespective of construction size.

## 6. Map and lost-data ledger

| Source -> target | Exact map/preserved predicate | Destroyed information / needed sidecar |
|---|---|---|
| Fixed ternary-labelled graph -> H1 obstruction | b->delta b; proper six-colour lift iff [1_q]=0 | q is existential and changes H_q; retain graph and odd-cycle witness |
| Graph gain -> double cover | (v,b), lifted edges | Cover bipartiteness loses section; C3->C6 hostile |
| Six palette pairing -> genus-two Lagrangian | three even pair classes modulo complement | Palette graph/metric absent; retain actual six branch points for geometry |
| Pairing -> central Fano pencil | add seventh label, then choose parity completion | One completion bit absent; no equivariant selector |
| Three pair frame transports -> Fano completion transports | p(x)=sum x | H1(Gamma;ker p), two obstruction bits per cycle |
| Pair flips -> six-root deck | pick one root per pair, x->sum x_j e_rj | Marking required; permutations are not sign automorphisms |
| Fano kernel -> canonical-character pair marginals | uniform projection | All 30 collapse; cubic correlations restore incidence |
| Full boundary relation -> pair shadows | existential coordinate projections | Higher constraints; NAE or at-most-two hostile |
| Abstract graph/gluing -> unit graph | only after an exact injective realization | Code/torsor provides no coordinates or cross-unit incidence |

## 7. Reproduction and scope

Run `python3 check.py > normal.json` and `python3 -O check.py > optimized.json`,
then compare. The checker has no bare assertions. It exhausts all 76 graphs
of orders 0..4, their 761 binary edge-cochain instances, 5,422 fixed ternary
cases, all 30 Fano kernels (deduplicated from 5,040 label bijections), all
1,050 triple projections, all 15 matching/Lagrangian codecs, and the declared
geometric/abstract hostiles. It does not enumerate difficult six-colour
unit graphs or re-prove Session60 surface geometry. See output for the final
always-active gate count. No literature-wide novelty claim is made.

**Remaining precise target:** retain the full (not scalar) lift obstruction
and complete boundary relation while producing a new exact unit-distance
assembly that escapes inherited field/extension screens. The Fano and K3
connections specify which data survives; they do not provide that assembly.
