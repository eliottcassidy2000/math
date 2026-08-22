---
id: THM-3278
title: "Selector origin bit, weighted-core bipartition, and critical-character boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The two THM-3275 selector conflicts induce the same 5+7 availability
  partition on THM-3269's twelve-row core, and it is exactly the unique
  root-normalized bipartition of all 22 core edges.  The canonical primitive
  edge removes THM-3266's supplied pair-identity sidecar globally, leaving
  exactly one necessary and sufficient face-origin bit on both complete
  faces.  The six canonical phase-orbit edges turn the necessary face bit
  into a common orientation gauge covering every nonzero J12 residue.  The bit is an exact cut
  coordinate, not a critical-group character.  Its nonlinear rank-label
  difference spans index 35 in the integral augmentation lattice, and its
  primitive genuine-phase vector is cyclic for THM-3274's abstract
  norm-fibre transfer.  None of these abstract bridges is a physical walk.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-03
audit: >
  The assertion-independent exact companion pins the theorem, script and
  transcript triples for THM-3269, THM-3274 and THM-3275.  It reconstructs
  the 22-edge core and both conflict availability banks; exhausts every
  root-normalized F2 colouring; checks the primitive and six-edge selectors;
  verifies every orientation residue and all three character hostiles; and
  computes the exact resultants, circulant Smith forms, augmentation index
  and norm-fibre Krylov determinant.  It also checks the global canonical
  covering pair, the phase-backbone completion and every raw seam-coordinate
  projection through the sharp seven-count cut decoder.  Normal, optimized
  and stored replay plus LF hashes agree.  An independent hostile audit
  rederived both conflict cuts, all core crossings, the global pair and
  backbone completion, every sampler orientation and character hostile, the
  augmentation Smith calculation in an independent basis, the norm-fibre
  cyclic vector and both decoder censuses.  It found no repair.
depends_on:
  - THM-3266-all-common-two-face-row-pairs-require-one-origin-bit
  - THM-3269-scale-invariant-clutch-strength-and-canonical-weighted-bispanning-polarization
  - THM-3274-norm-fibre-constrained-phase-transfer-and-refinement-invoice
  - THM-3275-unrestricted-twenty-two-row-face-blind-selector-obstruction
  - THM-3277-weighted-critical-phase-geodesic-backbone-and-exchange-subatlas
related:
  - THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary
  - THM-3273-critical-group-c12-quotient-and-relative-c7-equivariance-boundary
script: 04-computation/fc3_selector_origin_bipartition_phase_bridge_thm3278.py
output: 05-knowledge/results/fc3_selector_origin_bipartition_phase_bridge_thm3278.out
script_sha256: 07cf5cc1056bdd978a3a1d1146fee8139bef9a84e3b83aad08a0522b4df63a0f
output_sha256: 5e67eec14698f5258d17c5ba1ed295de3b9ff355aa279adb6d523ce3df683c7a
hash_basis: LF-normalized bytes
---

# THM-3278 -- selector origin bit, weighted-core bipartition, and critical-character boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3275 proves that an untagged count vector admits no common response-row
selector on exactly two states.  THM-3269 independently selects a rooted,
integrally polarized twelve-row core and an eleven-coordinate phase sampler.
The two structures meet exactly, but the meeting is a cut coordinate rather
than the missing physical phase character.

## 1. The selector seam is the core bipartition

Let `V_0` be the twelve vertices of THM-3260's 22-edge core.  At each of
THM-3275's two conflict states

```text
(3,4,5),       (1,3,4,5),                              (1)
```

restrict the small- and full-face lawful-row sets to `V_0`.  Both states give
the same two sets:

```text
S={2,11,16,18,22},
F={3,7,10,13,17,19,21}.                                (2)
```

They partition `V_0`, and every one of the 22 core edges has one endpoint in
`S` and one in `F`.  Thus `(S,F)` is the graph's bipartition.  Because the
core is connected, fixing THM-3269's canonical root `17 in F` makes it
unique.

Equivalently, over `F_2` define the vertex zero-cochain

```text
c(v)=1_(v in S),       c(17)=0.                         (3)
```

Then

```text
dc(e)=1                                                     (4)
```

on every core edge, and `(3)` is the unique root-normalized solution of
`(4)`.  The origin bit is therefore an exact cut coordinate.  It is not a
cycle or holonomy class: the constant-one edge cochain is the coboundary
`dc`.

This internalizes but does not erase THM-3275's sidecar.  The same untagged
count vector still occurs on both faces, so the value small/full cannot be
recovered from that vector alone.

## 2. A canonical endpoint selector and orientation gauge

THM-3269 roots at 17 and proves that `(16,17)` is its unique incident tree
edge whose rooted difference generates the full critical group.  Equation
`(2)` places 16 in `S` and 17 in `F`; moreover, at both conflict states row 16
is lawful on the small face and row 17 is lawful on the full face.  Hence

```text
small face |--> row 16,       full face |--> row 17     (5)
```

is a canonical paired selector on the entire two-state dependency locus.

There is a stronger global conclusion.  THM-3266 proves that `(16,17)` is one
of the 24 row pairs covering both complete faces.  It formerly treated pair
identity as supplied.  THM-3269 now selects this pair intrinsically, so the
controller

```text
given face f and state n:
  choose row 16 if it is available on f at n; otherwise choose row 17       (5a)
```

is lawful on all `238+4318` nonreset decisions.  Its pair-identity sidecar has
cost zero.  One face-origin bit is sufficient by `(5a)` and necessary by
THM-3275 even when all 22 rows are allowed.  Thus the exact remaining static
sidecar cost is one bit.  The pair `(16,17)` itself has nineteen conflict
states, all with the single orientation small-16/full-17; the unrestricted
all-row dependency locus shrinks to the two states `(1)`.

THM-3277 further shows why phase geodesics alone do not select this bit.  Its
seven-edge phase-minimum backbone has five forest components and leaves
13, 16 and 17 isolated.  The exact tree completion is

```text
T_* minus E_back={(2,13),(2,17),(3,16),(16,17)}.        (5b)
```

Thus the origin carrier `(16,17)` is part of the minimal rooted gauge
completion omitted by phase-geodesic data.  THM-3273's delayed `(11,17)` then
adds the missing direct type-4 evaluation.  Phase-minimum forest, rooted
origin gauge and direct unimodular sampler are three distinct layers.

The relation extends to all six edges of THM-3269's unimodular sampler.  Each
crosses `(S,F)`.  Orient it from the small class to the full class and record
its genuine normalized `J_12` increment:

```text
type   edge       S -> F       increment
 1    (16,17)     16 -> 17        11
 2     (2,21)      2 -> 21        10
 3    (18,19)     18 -> 19         3
 4    (11,17)     11 -> 17         8
 5    (21,22)     22 -> 21         7
 6     (7,22)     22 ->  7         6.                  (6)
```

The reverse orientations contribute `1,2,9,4,5,6`; together the two
orientations cover every nonzero residue of `C_12`.  Thus the necessary face
bit is precisely a common sign/orientation gauge on the six phase-orbit
representatives.  Equation `(6)` is still a coefficient atlas, not a claim
that these graph edges are lawful physical transitions.

## 3. Why the bit is not a phase character

Root normalization removes affine constants.  A homomorphism from
`Jac(G_0)=C_74748` to `C_2` is either zero or parity in a primitive cyclic
coordinate.  But row 11 has

```text
c(11)=1,       exponent_Jac(11)=9088 == 0 mod 2.        (7)
```

So `(3)` does not factor through the full critical group.  The exact mismatch
set for critical parity is `{2,7,11}`.

The same witness survives both twelve-label reductions:

```text
j_12(11)=4,       ell(11)=4,       c(11)=1,             (8)
```

where `j_12` is the genuine group-linear phase and `ell` is THM-3269's
nonlinear rank label.  Their parity mismatch sets are respectively
`{2,7,11}` and `{2,3,10,11,13,18,21}`.  Hence neither a `J_12` character nor a
character of the rank-labelled `C_12` carries the origin bit.  The edge
orientation in `(6)` contains strictly different information from a vertex
phase character.

## 4. The nonlinear rank label gives an integral finite-index bridge

Although `ell` is not group-linear on vertices, it is a canonical bijection.
In

```text
R=Z[x]/(x^12-1),       N=1+x+...+x^11,                 (9)
```

the small-class indicator becomes

```text
f=x+x^2+x^4+x^5+x^10.                                 (10)
```

Exact polynomial and Smith calculations give

```text
gcd(f,x^12-1)=1,
Res(f,x^12-1)=175=5^2*7,
SNF(Circ(f))=(1^10,5,35).                              (11)
```

Apply THM-3260's cyclic-difference map:

```text
g=(x-1)f=-x+x^3-x^4+x^6-x^10+x^11.                    (12)
```

Then

```text
gcd(g,x^12-1)=x-1,
SNF(Circ(g))=(1^10,35,0).                              (13)
```

Therefore the cyclic orbit of `g` has rational rank eleven, is nonzero on
every charged Fourier mode, and spans an index-35 sublattice of
`Aug(Z[C_12])`.  Under the inverse cyclic-difference bridge it recovers the
cut potential `f` modulo constants.

This is an integral abstract bridge with a sharp `5*7` defect, not a physical
`C_7` carrier.  Its cyclic coordinate is the nonlinear rank label `ell`, not
the genuine vertex phase `j_12`.

## 5. A genuine-phase cyclic vector for the norm-fibre transfer

The canonical paired selector `(5)` has genuine phases

```text
j_12(16)=1,       j_12(17)=0.                           (14)
```

Thus its signed phase vector is `v=e_1-e_0`.  Let `Q_0` be THM-3274's
symmetric twelve-state quotient for a freshly varying increment in one norm
fibre.  Its characteristic polynomial is irreducible of degree twelve over
`Q`.  Consequently `Q^12`, as a `Q[Q_0]`-module, is simple: the annihilator of
any nonzero vector divides the irreducible minimal polynomial and hence must
equal it.  In particular,

```text
span_Q{v,Q_0 v,...,Q_0^11 v}=Q^12.                     (15)
```

The direct exact certificate is

```text
det[v,Q_0v,...,Q_0^11v]
 = -30396564468585450830924251521449984 != 0.           (16)
```

Thus one signed primitive-edge origin vector seeds the whole constrained
phase quotient.  This does not contradict Section 3: `Q_0` is an abstract
fresh-increment transfer, not a group character, and it need not preserve the
augmentation subspace.

There is also a sharp query boundary.  Pull the Boolean rank-label word

```text
beta(r)=1_(r in {1,2,4,5,10})                          (17)
```

back along the abstract `C_12` phase of THM-3274's seam source.  Ask which
coordinates of its raw target-phase histogram `H_N(x)` determine only
`beta(phi(x))`, rather than the whole source.  Exhaustion gives

```text
minimum raw count coordinates=7,       successful 7-sets=95,
first={0,1,2,3,4,5,11};

minimum with total degree appended=7,  successful 7-sets=122.              (18)
```

No set of at most six counts determines the cut bit, even with degree.  Thus
this particular Boolean sidecar is as hard in minimum coordinate width as
THM-3274's full 168-point decoder, although more seven-sets decode the bit.
Equation `(18)` is an abstract identification of two `C_12` coordinate sets;
it is not a physical identification of a response row with a Singer point.

## 6. Scope and failure boundary

The theorem concerns exactly the two promoted response faces.  Its all-row
bipartition statement concerns the two conflict vectors in `(1)` and the
12-row core; its one-bit controller conclusion uses THM-3266's complete
facewise pair cover.  Rows outside the core remain lawful alternatives and
are not classified by the bipartition.

The face-origin bit remains necessary.  Equations `(12)--(18)` neither infer
the face from an untagged count vector nor construct a physical response walk,
positive moment functional, owner current, or target-preserving transition.
They prove no row exclusion, no `FC(3)` or `SFC(3)`, no `LRC(14)` decrement,
and no new Gaussian Moment Conjecture consequence.

## 7. Exact verification

Run

```text
python 04-computation/fc3_selector_origin_bipartition_phase_bridge_thm3278.py
python -O 04-computation/fc3_selector_origin_bipartition_phase_bridge_thm3278.py
```

and compare LF-normalized bytes with

```text
05-knowledge/results/fc3_selector_origin_bipartition_phase_bridge_thm3278.out.
```

The companion uses exact integer, finite-field-transcript and polynomial
arithmetic only.  It has no assertion node, floating literal, randomness or
fitted recurrence.

QED.
