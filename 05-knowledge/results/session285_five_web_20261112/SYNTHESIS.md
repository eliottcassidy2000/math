# Session285: what generates five, and what survives when the scale changes

**Status: PROVED scoped results / inherited classical mechanisms / historical
recovery / explicitly labelled hypotheses.** This session does not prove
Chan–Pak or identify one object underlying all occurrences of five. It builds
several exact bridges that were unavailable from matching lists alone.

## 1. The supplied stress statement has an exact mechanism

Fix a regular local isothermal coordinate w=f(z), with h=(f')² nonzero. Set

    p=σxx+σyy,       q=σxx−σyy−2iσxy.

For equilibrated plane stresses without body forces and with harmonic trace,
the prescribed net is a principal-direction net exactly when q/h is real.
Compatibility then forces

    q=h[a|z|²+bz+bar(b)bar(z)+c],
    p=d−2 Re ∫h(z)[az+bar(b)] dz,

with a,c,d real and b complex. Thus **four real amplitude modes plus one
pressure anchor** explain the five parameters. The complexity of the fixed
analytic net h is not counted among them.

This recovers the classical Neményi mechanism, rather than claiming a new
five-parameter theorem. The [self-contained derivation](stress/PROOF.md)
also supplies exact five-state flat transport, a discrete four-mode analogue,
and a two-net decomposition of general compatible local stress, following
the classical Wegner mechanism. A fixed net need not represent every stress;
changing the net is a productive operation, not a failure of the theorem.

The amplitude is literally a 2×2 Hermitian form. Its real rank-one slice
R=|mz+n|² gives the Pythagorean null vector
(m²−n²,2mn,m²+n²). This is an actual map, with primitive/parity/chamber and
net/pressure information explicitly retained or lost. It is stronger than
saying that the triple(3,4,5) contains the same numeral.

## 2. Global topology turns five into two, three, four or five

On a multiply connected domain the pressure primitive must have zero real
periods. The [global theorem](bridges/PERIOD_FIBERS.md) gives

    dimension = 5 − rank_R[(a,Re b,Im b) ↦
                          (Re ∮γ h(z)(az+bar b) dz)_γ].

The rank is at most three, even with many holes. Rational isothermal
coordinates explicitly realize every dimension2,3,4,5. This is an exact
classification of stress gluing under the stated assumptions, not a claim
about globally single-valued displacement or arbitrary traction boundary data.

The full monodromy is an abelian rank-one shear of the pressure anchor;
its deviation algebra is square-zero. The number of holes alone is not the
sufficient statistic: the period operator is. This gives a precise version
of the recurring repository lesson that local compatibility needs a global
continuation sidecar.

## 3. Tournament attachments: real and physical fibers differ

Fix a deletion tournament D. Its root attachments have fixed deletion
coordinates and five variable fusion coordinates(H,U,V,d,e): four affine
responses and one quadratic path count. The [attachment theorem](tournaments/PROOF.md)
gives exact affine dimension rank(u,v,a,b)+1 for deletion order at least2.
The full dimension five first occurs at total order5.
This is a finite set spanning five dimensions, not five continuous physical
parameters and not an isomorphism with stresses.

For transitive D, equality of the four continuation responses gives an exact
real compatibility operator with stencil(2,−5,2). Its roots are2 and1/2,
not golden-ratio roots. Compatible real directions first acquire indefinite
count curvature at deletion order9, with exact witness−2. Nevertheless,
**every binary compatible fiber is a singleton at every order**, because
one response coordinate is a binary code of the cut.

There is a positive escape. In a vertex-transitive core every fixed-cardinality
root cut has identical four linear responses. Count variation becomes an
exact weighted induced-subset energy. Unequal pair-orbit energies create actual
compatible count switches; unordered-pair transitivity destroys all such
variation. The inherited41/43 gate is explained by C5's two pair energies.
Thus the useful middle regime is vertex symmetry without total pair symmetry.
Changing the cyclic core produces a literal ten-vertex full-central unit
switch with counts14,207 and14,209, deletion count3357, and identical
(U,V,d,e)=(18906,34455,2238,1119). An independent direct path counter confirms
both witnesses. The count3357 is coprime to the carousel7 deletion175; this
adds a concrete different interface, not a universal coprime supply theorem.
Full witness details and inherited/new distinctions are in the proof and receipt.

These results describe **where freedom lives**, not just whether the scalar
count changes. Neither continuous dimension, symmetry, multiple pair orbits
nor coefficient gcd ensures a gapless physical interval.

## 4. The five-state pentagon has a complete all-scale continuation

The old ready-tail pentagon is the k=2 case of a marked partial-permutation
operation graph. The [all-order grammar](holonomy/PROOF.md) has

    v_k=Σ_{m=0}^k k!/(k−m)!,       v_k=k v_(k−1)+1,
    v_0,v_1,... =1,2,5,16,65,326,1957,13700,... .

Fix a core and a signature for each potential tail name. A state is a
transitively ordered subset of names. Edges append, delete the last name, or
swap two adjacent names by one arc reversal. These graphs are the classical
stellohedra; the marked-tournament realization is explicit.

Two relation schemas suffice at **every scale**: contextual pentagons and
append/swap transport squares. Their two-dimensional complex collapses to
the canonical append tree, so the face boundaries form a minimal integral
cycle basis. The number of equations grows, but no new relation type is
needed. An explicit filling law bounds global additive error by local face
defects. Pentagon checks alone fail already at three labels.

This is a genuine compact recursive understanding of an infinite tournament
family. It is not a fixed-order unmarked atlas or a universal five-state
model. A five-vertex pentagon has five potential values, but the next graph
has sixteen; a further closure law is necessary for bounded-dimensional
memory. The literal hostile H(AB)=H(BA)=5 followed by H(ABC)=11 and
H(BAC)=13 shows why scalar counts cannot replace the state.

## 5. The historical web and its limits

The [mechanism atlas](history/MECHANISM_ATLAS.md) starts with five/fivefold,
sets and parameters of five, pentagons/quintics, Platonic/Fermat/Ogg themes,
and previous searches. It follows their neighbors into involutions,
curvature, integral modules, symplectic matchings, eta cancellation,
FRACTRAN, cellular automata, receiver interfaces and arithmetic height.
The pinned search reads34,864 of34,887 candidate prose paths and finds7,509
seed-hit paths plus10,069 expansion-only paths. These are lexical leads, not
counts of unique sessions read. About30 sources and four deterministic
lateral samples are examined in depth; the lateral rule excludes obvious
five-related filenames. Sources span retained March–November work; source prose is not
silently promoted to current theorem status.

The machine-readable [search summary](history/summary.json),
[compressed inventory](history/inventory.json.gz),
[compressed expansion](history/expansion.json.gz),
[commit-history sample](history/commit_history.json), and omitted-path list
state exactly what was searched. A pinned-tree lexical inventory is not a
claim to have read every hit or every deleted version. Shallow commit history
and unavailable/oversized blobs are explicitly separated.

One recovered exact bridge is the C5 golden integral frame: its constant
plus augmentation split has index5. The session extends that elementary
index statement to every n:

    Z^n / (Z(1,...,1) + {x:Σx_i=0}) ≅ Z/nZ.

The same n is therefore both a coordinate count and an integral obstruction
index, by an actual quotient map. Golden eigenvalues do not generalize to
every n. This illustrates how all numbers can be studied deeply without
assuming that every occurrence has one common cause.

## 6. What unifies the advances, and the next obligations

The proposed common language is **compatible-state fibers**:

    fixed background → continuation carrier → compatibility equations
                     → physical realization set → target observer.

The constructive question is whether an actual state lies over a prescribed
target, not merely whether the carrier has enough dimension. The stress
period law, tail relation spine and tournament binary-fiber obstruction are
three precise, different instances of this distinction.

The [research program](bridges/UNIFYING_PRINCIPLES.md) names four explicit
next routes: full-state relation completion; orbit-energy interval growth;
controlled changes of the background core; integral obstruction modules.
These are testable hypotheses and designs, not claimed reductions between
open problems. For Chan–Pak, a universal absolute-target section or a renewable interval-growth mechanism
sufficient for all allowed odd targets remains missing. **Chan–Pak remains OPEN.**

[Review](audit/REVIEW.md) · [Reproduction](REPRODUCE.md) · [Release](RELEASE.md)
