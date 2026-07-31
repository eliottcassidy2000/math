---
id: THM-2814
title: "Projective allocation-square holonomy and idempotent provenance no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Independent
  atom-preserving occupancy idempotents over a field have zero mixed
  face on every fourfold-co-supported raw atom, and arbitrary linear
  coarse-graining keeps their face in the joint-absent sector.  There are two
  distinct non-idempotent escapes.  On one physically normalized common line,
  a Segre square w(1,alpha,beta,alpha beta) has nonzero face exactly when both
  contrasts are nontrivial; its cross-ratio remains one.  With only four
  vertex lines modulo row/column rephasing, the complete invariant is instead
  kappa=v00v11/(v10v01), and kappa!=1 is the first intrinsic square holonomy
  on the nonzero field locus (equivalently, the all-units locus over a ring).
  THM-2593 supplies abundant coefficient contrasts but a vertex-coboundary
  gauge; THM-2779 supplies nontrivial holonomy abstractly; THM-2806/2807 and
  the THM-2771 local copy germ supply neither physical invoice.  Thus the
  remaining LRC gate is a typed physical realization, not another formal
  square statistic.
source: root/projective-allocation-holonomy-2026-07-28
depends_on:
  - THM-2593-charged-target-section-atlas-and-minimal-c91-holonomy-trivialization
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary
related:
  - THM-2556-reynolds-duty-curvature-and-fibre-covariance-mixed-cell
  - THM-2813-affine-lift-transvection-and-projective-horn-decoder
script: 04-computation/lrc14_projective_allocation_holonomy_thm2814.py
output: 05-knowledge/results/lrc14_projective_allocation_holonomy_thm2814.out
script_sha256: f5198812da7be1051037b2318e63d5ff5cc5d5173fba0e61c32fac49abf5ce7d
output_sha256: 8390fcea84c7df3247efd366b754ad27423988613a2253a806da93fbbfbae56a
secondary_script: 04-computation/lrc14_projective_allocation_holonomy_independent_audit_thm2814.py
secondary_output: 05-knowledge/results/lrc14_projective_allocation_holonomy_independent_audit_thm2814.out
secondary_script_sha256: 86585db9e7bc56ed221167e8414b7e626737c8a5652c28af0494195feec5fcad
secondary_output_sha256: 5abebb920af9c0a5b3ea97de4a5850d68baec825c8ce6501f9fffe67f6324921
hash_basis: LF-normalized bytes
---

# THM-2814 -- allocation contrast and projective square holonomy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2806 proves why the most literal fixed-sheet allocation misses
THM-2772's gate: the only raw atom supporting all four states is flat, while
the nonzero central face comes from atoms on which both carriers are absent.
This theorem classifies the first ways one could change that conclusion.

There are two answers because there are two different coefficient objects:

```text
Branch A: one physical coefficient line with fixed edge normalization;
Branch B: four vertex lines known only up to row/column rephasing.          (0)
```

In Branch A, ordinary non-idempotent source and target contrasts suffice and
the square remains rank one.  In Branch B, those contrasts are gauge choices;
the first intrinsic interaction is projective square holonomy.  Conflating
the branches would incorrectly reject THM-2772's Segre mechanism, whose
cross-ratio is necessarily one.

## 1. Idempotent common-atom and coarse-graining no-go

Let `K` be a field, let `x` be one full physical atom with bare amplitude
`w_x`, and let the independent source and target occupancy restrictions act
on that atom by idempotents `alpha_x,beta_x`.  Then

```text
(B,P,Q,H)_x=w_x(1,alpha_x,beta_x,alpha_x beta_x),

mu_x=B_x-P_x-Q_x+H_x
    =w_x(1-alpha_x)(1-beta_x).                           (1)
```

The only idempotents of a field are zero and one.  Fourfold co-support forces
`w_x,alpha_x,beta_x` nonzero, hence `alpha_x=beta_x=1` and `mu_x=0`.
Thus independent atom-preserving occupancy restrictions can never give a
nonzero mixed face on their common raw support.  THM-2806's wandering-delta
law is one instance of this general obstruction.

Nor can a linear pushforward repair the provenance.  On a finite atomic
space let `E,F` be Boolean masks, `f` any amplitude, and `Lambda` any linear
map.  Put

```text
B=Lambda(f),       P=Lambda(Ef),
Q=Lambda(Ff),      H=Lambda(EFf).                        (2)
```

Then

```text
B-P-Q+H=Lambda((1-E)(1-F)f).                             (3)
```

The source of `(3)` is the joint-absent sector.  Its restriction to `E=F=1`
is identically zero, regardless of whether `Lambda` is a positive average,
Reynolds sum, signed functional, orbit sum, or fibre pushforward.  The sharp
four-atom hostile is

```text
(B,P,Q,H)=(4,2,2,1),                 mu=1,               (4)
```

whose face is the unique both-absent atom.  All four coarse values being
nonzero does not make `(4)` a fourfold-co-supported raw face.

The companion exhausts every pair of four-bit masks and every weight vector
in `{-1,0,1}^4`: all `20,736` cases satisfy `(3)`.

## 2. Branch A: a physically normalized common line

Drop idempotence but retain one independently normalized physical coefficient
line.  If source and target toggles preserve the atom and act by scalars
`alpha,beta`, equation `(1)` remains both necessary and sufficient for an
independent rank-one square.  With all four entries nonzero,

```text
mu!=0  iff  alpha!=1 and beta!=1,                        (5)

kappa=(BH)/(PQ)=1.                                      (6)
```

Over `F_13`, there are exactly

```text
(13-2)^2=121                                             (7)
```

ordered nonzero contrast pairs.  This is exactly THM-2772's desired Segre
branch: its physical source and target amplitudes fix the edge normalization,
so nontrivial projective holonomy is neither necessary nor expected.  The
missing object is a lawful physical source of `alpha,beta!=1` on one raw
common atom.  THM-2806 has `alpha=beta=1` at its only such atom.

An external allocation sign illustrates why the normalization clause is
load-bearing:

```text
(1,-1,-1,1),                        mu=4,

Fourier_(C2^2)(1,1,1,1)=(4,0,0,0),
Fourier_(C2^2)(1,-1,-1,1)=(0,0,0,4).                    (8)
```

On a preidentified common line, a physically inherited sign would be a valid
contrast.  Without such an inheritance, `(8)` merely changes which external
allocation character is called mixed.

### 2.1 A sharp signed model and its positivity boundary

THM-2716 supplies the minimal integral coefficient model.  On `Z^2`, let

```text
J=[[0,-1],[1,0]],                    J^2=-I.              (9)
```

For `v=(1,0)`, using `J` on both toggles gives

```text
(v,Jv,Jv,J^2v),
(I-J)^2v=-2Jv=(0,-2)!=0.                                (10)
```

The readout `(x,y)->x+y` gives `(1,1,1,-1)` and mixed face `-2`.
After adjoining `i`, this is the rank-one square `(1,i,i,-1)`.
It is a genuine same-base-atom non-idempotent coefficient model, but
THM-2716 proves that `J^2=-I` preserves no nonzero pointed cone.  No theorem
attaches `J` to the positive THM-2791 allocation edges.

### 2.2 The THM-2593 coefficient atlas is nonvanishing but gauge-exact

In `R_7=F_13[z]/(Phi_7)`, THM-2593's actual coefficient transports are

```text
u_(q,a)=Y_(q+a)Y_q^(-1),       q in F_13, a!=0.          (11)
```

A fresh exact census gives

```text
156/156:       (1-u)^2 !=0;
152/156:       1-u is a unit;
24,336/24,336 ordered products (1-u)(1-v) are nonzero;
23,104=152^2 of those products are units.                (12)
```

The four nonunit contrasts are

```text
(q,a,q+a)=(9,5,1),(12,5,4),(1,8,9),(4,8,12).           (13)
```

The independent factorization

```text
Phi_7=(z^2+3z+1)(z^2+5z+1)(z^2+6z+1)                  (14)
```

explains the absence of annihilation.  Each exception vanishes only in the
last quadratic factor and survives in the other two; all other contrasts
survive in all three.

This is a large exact Branch-A coefficient control, not its physical
realization.  Equation `(11)` is a vertex coboundary, and THM-2593's gauge
`Y_q^(-1)` sends every edge transport to one.  Its endpoints are distinct
target sheets, and pairing arbitrary contrasts in `(12)` does not construct
one common atom or a two-carrier square.

## 3. Edge characters and the fixed-sheet gauge

The tempting rank-one local system on a marked `C_13 x C_13` orbit is

```text
P(a)=chi(a)B,       Q(b)=psi(b)B,
H(a,b)=chi(a)psi(b)B.                                   (15)
```

Algebraically its face is `B(1-chi(a))(1-psi(b))`.  It fails twice on the
THM-2791/2806 object.

First, physical presence has masks `delta_0` on both carrier axes.  The only
fourfold raw point is `(a,b)=(0,0)`, where every normalized character equals
one.  Filling the orbit by character-weighted translated sheets changes the
object; it is not a reweighting of the fixed sheet.

Second, THM-2806's representative gauge is

```text
(ell,a,b)->(ell+W,a+1,b-1).                             (16)
```

The four vertex lines in `(15)` acquire relative factors

```text
(1,chi(1),psi(-1),chi(1)psi(-1)).                       (17)
```

If all four values are required to descend on one common scalar line, the
nonzero bare vertex forces the common gauge factor to be one; the next two
vertices then force both characters trivial.  If instead the four graded
lines are retained, row/column parallel transport divides out `(17)` and
flattens the square.  Thus edge characters give either a chart-dependent
Branch-A number or a flat covariant square, not the missing physical face.

## 4. The THM-2806 associated-graded square is flat

The central fixed-sheet vector has

```text
v=w(13^2,13,13,1),
valuation(v/w)=(2,1,1,0).                               (18)
```

Its valuation curvature is `2-1-1+0=0`.  Dividing by the intrinsic powers
of thirteen gives four equal leading units; raising all vertices to grade
two gives `w(13^2,13^2,13^2,13^2)`.  Consequently the canonical Rees or
associated-graded comparison is flat.  The ungraded value

```text
(13-1)^2w=144w= w mod13                                 (19)
```

is the tensor product of two orbit-cardinality boundaries, not a
two-dimensional Bockstein curvature.  A filtered escape must add a
nonseparable sidecar rather than merely reduce `(18)` modulo thirteen.

## 5. Branch B: four vertex lines and projective holonomy

Let `K` be a field, and suppose the four nonzero coefficients live in vertex
lines whose
trivializations may be changed independently by a source-row and target-column
gauge:

```text
v_ij -> r_i s_j v_ij,                 r_i,s_j in K^*.    (20)
```

The cross-ratio

```text
kappa=(v_00 v_11)/(v_10 v_01)                           (21)
```

is invariant.  It is complete: take

```text
r_0=v_00^(-1),   r_1=v_10^(-1),
s_0=1,           s_1=v_00 v_01^(-1).                   (22)
```

Then `(20)` sends every square uniquely to

```text
(1,1,1,kappa).                                          (23)
```

The one-dimensional kernel `(r_i,s_j)=(t,t^(-1))` is the common scalar
redundancy.  Over `F_13`, the `12^4=20,736` nonzero squares split into twelve
cross-ratio fibres of size `12^3=1,728`; the companion also checks `248,832`
row/column rephasings.

The additive Mobius value, including its vanishing, is not invariant under
`(20)`.  Two sharp warnings over `F_5` are

```text
(1,2,2,4):                  kappa=1, mu=1,
row/column normal form:     (1,1,1,1), mu=0;

(1,2,2,3):                  kappa=2, mu=0,
row/column normal form:     (1,1,1,2), mu=1.             (24)
```

Equation `(24)` does not invalidate Branch A: there the common-line edge
normalization is physical and `(20)` is not a gauge freedom.  But when only
the four lines and the row/column-factorized gauge are supplied, `kappa!=1`
is the first intrinsic joint interaction.  The normalized chart `(23)` has
additive defect `kappa-1`; with an externally retained common scale `w`, this
is `(kappa-1)w`.  That defect is normalized/covariant, not a gauge-invariant
scalar: only `kappa` is intrinsic.  This branch departs from the rank-one
Segre square rather than improving it.

The scope is sharp.  Four unrelated vertex rephasings can change `kappa`.
Over a general ring, four merely nonzero entries need not be invertible, and
equal formal cross-products need not imply gauge equivalence.  The same proof
does extend verbatim from a field to the locus on which all four entries are
units.

Reversing the square orientation sends `kappa` to `kappa^(-1)`.  Hence any
tournament orientation extracted from the twelve nontrivial `C_13`
holonomies also needs an oriented frame or a chosen half-system; the unmarked
square supplies inverse pairs, not a canonical tournament.

## 6. The Heisenberg sidecar pays Branch B abstractly

THM-2779 acts on a thirteen-dimensional coefficient fibre by

```text
T e_r=e_(r+1),                  M e_r=zeta^(-r)e_r,
T M=zeta M T.                                             (25)
```

For the two chosen linear lifts, the oriented path multiplier `TM/MT` is
`zeta`.  After three lift trivializations its normalized scalar-holonomy
chart is

```text
(1,1,1,zeta),                    kappa=zeta!=1.           (26)
```

This is not a literal one-dimensional four-toggle square.  On the thirteen
basis lines, and after passage to `PGL`, the two shadows commute; the phase is
remembered only by the oriented comparison of the chosen linear lifts.  No
one-dimensional scalar pair realizes `(25)`, since scalars commute.  In
the endpoint plane, the phase of directions `s,t` is

```text
kappa(s,t)=zeta^det(s,t).                                 (27)
```

Exact enumeration over `F_13^2` gives `26,208` ordered nondegenerate frames
and `2,184` frames for each nonzero determinant/phase.  Swapping `s,t`
inverts the phase.  Thus THM-2779 contains exactly the joint interaction
missing from the edge-character proposal.

It remains abstract.  Its coefficient phase is killed by the thirteen-root
Boolean permutation action, while the faithful `169`-point endpoint carrier
has no proved same-ancestry action on the THM-2791 atom.  Equation `(26)` is a
target specification, not a physical LRC lift.

## 7. The commuting simplex and local cofiber copy germ

THM-2807's address edges are

```text
13,             169*4079=689351,             689364,
13+689351=689364.                                          (28)
```

They form an honest commuting translation triangle, so every character has
holonomy one.  For example, over `F_53` the character completion
`(1,10,15,44)` has `kappa=1` although its chart-dependent additive face is
`20`.  The residue `4079=10 mod13` is a retained full-address vertical edge,
not curvature.  Modulo `169`, the vertical endpoints coincide, and every one
of the thirteen `F_53`-valued quotient characters kills that edge.  The
thirteen affine lifts do not change this: they are address conjugacies, not
proved allocation-edge actions.

The two equal THM-2771 cofiber copies have full-address shifts

```text
-13,                    +1248=96*13.                    (29)
```

Their depth-one residues are `12,5`.  In
`F_13[epsilon]/(epsilon^13)`, their copy polynomial begins

```text
z^12+z^5=2+4epsilon+...,
9(z^12+z^5)=5+10epsilon+... .                            (30)
```

This is a unit germ, not an augmentation uniformizer.  More importantly, the
copies are disjoint full addresses: both fail native `E3`, the positive copy
also fails native `c2`, and their source carrier masks are empty.  Equal
content and ancestry labels do not make `(30)` a common-atom square.

## 8. Exact evidence and boundary

Run

```text
python 04-computation/lrc14_projective_allocation_holonomy_thm2814.py
python -O 04-computation/lrc14_projective_allocation_holonomy_thm2814.py
```

Both modes byte-match the stored output.  The dependency-free companion has
no Python `assert` nodes and independently reconstructs the THM-2593
`R_7/Phi_7` arithmetic, all row/column square classes, every symplectic frame,
the Heisenberg phase, the fixed-sheet Rees profile, the commuting address
triangle, and the local cofiber germ.  A second dependency-free implementation
also byte-matches its stored output under normal and optimized Python.  It
reconstructs the orbit theorem over `F_3,F_5,F_7,F_13`, tests `248,832`
linear provenance cases, and supplies the unrelated-vertex, nonunit-ring,
product-ring-idempotent, projective-shadow, and quotient-character hostiles.

Exactly proved are the field-idempotent common-atom and linear
provenance no-gos, Branch-A contrast classification, Branch-B cross-ratio
classification, and the stated exact applications.  Not proved are a
physically normalized non-idempotent toggle on one LRC atom, a same-ancestry
Heisenberg fibre, allocation-to-endpoint transport, root/Cech correction,
row exclusion, or LRC(14).

**QED.**
