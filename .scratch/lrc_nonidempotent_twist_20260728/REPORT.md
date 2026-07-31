# Non-idempotent same-atom allocation: exact classification and current boundary

**Status:** abstract lemmas proved below; finite applications are exact scratch
checks against the current THM-2806/2807 and THM-2771 cofiber constants.  This
is not a row exclusion or an LRC(14) proof.

## Inheritance and concept board

- Closest proved mechanism: THM-2806, whose raw fixed-sheet square is
  `w(1,delta_0,delta_0,delta_(0,0))` and whose only fourfold point is flat.
- Canonical hostile: THM-2772's post-Fourier one-corner vector
  `(0,0,0,h)`.
- Corrected near miss: THM-2806's central `D3=144w` is a bare-only marginal,
  not a common-atom face.
- Least-used sidecar: THM-2556's Reynolds multiplier, the canonical example
  of a non-idempotent amplitude produced by averaging an idempotent mask.
- Live objects: raw Boolean atoms; rank-one local systems; representative
  gauge; Rees leading units; Heisenberg cocycle; the THM-2807 address
  triangle; the two equal THM-2771 cofiber copies.

## 1. Universal atomwise classification

Let `R` be a commutative ring, let `x` be one full physical atom, and let
`w_x` be its bare amplitude.  If the source and target toggles preserve `x`,
they act on its coefficient line by scalars `alpha_x,beta_x`.  Independent
composition gives

```text
(B,P,Q,H)_x
 =w_x(1,alpha_x,beta_x,alpha_x beta_x),

mu_x=B_x-P_x-Q_x+H_x
     =w_x(1-alpha_x)(1-beta_x).                     (1)
```

This is both necessary and sufficient for a rank-one, atom-preserving,
independent square.

If the toggles are physical occupancy projectors, then
`alpha_x^2=alpha_x` and `beta_x^2=beta_x`.  Over a field their values are
zero or one.  Fourfold co-support forces both to be one, so `mu_x=0`.
This proves a general no-go strictly stronger than the wandering-delta
special case:

> **Idempotent common-atom no-go.**  No pair of independent, atom-preserving
> occupancy restrictions over a field can have both fourfold co-support and
> nonzero raw mixed face.

The exact escape in rank one is also clear: both `alpha_x,beta_x` must be
nonzero and different from one.  Over `F_13` there are exactly
`11^2=121` such ordered pairs.  This is an algebraic classification, not a
physical construction.

There are two logically different branches, which must remain separate.

**Branch A: physically normalized common line (the THM-2772 Segre branch).**
If the source and target amplitudes have an independently fixed physical
normalization on one common coefficient line, then (1) is already the exact
answer:

```text
alpha_x,beta_x in R^x\{1}
  => mu_x=w_x(1-alpha_x)(1-beta_x) !=0                (1a)
```

over a field, provided `w_x!=0`.  Its cross-ratio is still one:

```text
(B_x H_x)/(P_x Q_x)=1.                                (1b)
```

So a nontrivial projective interaction is neither necessary nor expected in
this branch.  This is precisely THM-2772's desired factorized square.  The
missing input is physical: THM-2806's common atom has
`alpha_x=beta_x=1`.

**Branch B: only four vertex lines modulo row/column rephasing.**  If no
physical edge normalization identifies the vertex lines, then `alpha_x` and
`beta_x` can be changed by those rephasings and the additive value in (1) is
not intrinsic.  The first surviving projective datum is the cross-ratio
treated in Section 6.  A nontrivial cross-ratio is an alternative joint
interaction, not a requirement for Branch A.

## 2. Positive averaging cannot change the provenance

Let `E,F` be arbitrary Boolean masks on a finite atomic space and let
`Lambda` be any linear coarse-graining (positive averaging, an orbit sum,
a fibre pushforward, or even a signed linear functional).  Then

```text
Lambda(1)-Lambda(E)-Lambda(F)+Lambda(EF)
 =Lambda((1-E)(1-F)).                                  (2)
```

The integrand on the right is supported exactly where both carriers are
absent.  Its restriction to the fourfold-common locus `E=F=1` is identically
zero.  Therefore no linear aggregation of idempotent masks can turn its
mixed face into a common-atom contribution.

The one-fibre conditional-expectation hostile is

```text
(B,P,Q,H)=(4,2,2,1),                mu=1.              (3)
```

All four coarse amplitudes are nonzero and the multipliers `1/2` are
non-idempotent, but the entire face in (3) comes from the unique
both-absent source atom.  This simultaneously classifies:

- orbit/Reynolds averages such as THM-2556;
- unions or positive multiplicity sums of translated indicator sheets;
- the THM-2806 central marginal.

Overlapping translates help only if the **raw** toggle ceases to be an
occupancy indicator or if distinct full atoms are first identified.  In the
latter case the result is a quotient face, not a same-atom face.

## 3. The proposed rank-one local-system template

On a marked `C_13 x C_13` orbit, the proposed amplitudes are

```text
P(a)=chi(a)B,   Q(b)=psi(b)B,   H(a,b)=chi(a)psi(b)B,

mu(a,b)=B(1-chi(a))(1-psi(b)).                         (4)
```

Algebraically, (4) works wherever both characters are nontrivial.  It does
not arise from the existing fixed THM-2791/2806 sheet:

1. the fixed sheet has present masks `delta_0` on both carrier axes;
2. hence its sole fourfold point is `(a,b)=(0,0)`;
3. every normalized character has `chi(0)=psi(0)=1`;
4. therefore (4) is zero at the only physically co-supported point.

Replacing `delta_0` by a character fills the entire orbit.  It is not a
reweighting of the fixed sheet; it is the phase-weighted direct sum of all
thirteen translated sheets.

There is a second, independent obstruction.  THM-2806's representative
gauge is

```text
(ell,a,b) -> (ell+W,a+1,b-1).                          (5)
```

Under (5), the four vertex lines in (4) acquire factors

```text
(1,chi(1),psi(-1),chi(1)psi(-1)).                      (6)
```

Since the bare vertex is nonzero, common scalar descent would force the
common factor to be one; the source and target vertices then force
`chi=psi=1`.  Thus a nontrivial edge-character square does not descend to
one coefficient line on which its ordinary Mobius sum is typed.

If instead the four allocation-graded lines are retained, their canonical
equivariant parallel transport divides out `chi(a)` and `psi(b)`.  The
transported square is `(B,B,B,B)` and its covariant mixed face is zero.
Equivalently:

```text
phase relative to the marked origin = chi(a-a_origin);
on the supported sheet a=a_origin, this phase is one.      (7)
```

This is the exact gauge/overlap obstruction.  An edge character can give a
nonzero chart-dependent number or a flat covariant square, but not a
gauge-covariant raw common-atom face.

## 4. Why an allocation sign is not yet the answer

Multiplying present-source and present-target vertices by `-1` gives

```text
(1,-1,-1,1),                         mu=4.              (8)
```

On a **preidentified, physically normalized common coefficient line**, the
constant sign is invariant under THM-2806's representative gauge and is
fourfold co-supported.  It is not invariant under independent allocation
row/column rephasing when that physical normalization is absent.  Moreover,
the allocation `C_2^2` Fourier audit shows exactly what happened:

```text
Fourier(1,1,1,1)=(4,0,0,0),
Fourier(1,-1,-1,1)=(0,0,0,4).                         (9)
```

The external alternating sign moves the flat square's trivial coefficient
into the mixed-character slot.  It does not detect a physical interaction;
it changes which allocation character is called mixed.  It becomes a
physical escape only if the sign/phase is inherited from an independently
defined carrier operation with fixed normalization.  No such sign occurs on
the positive THM-2791 sheet.

### 4.1 A genuine signed coefficient escape already exists abstractly

The repository's quarter-turn local-system reflection supplies the minimal
integral version of a fixed nontrivial phase.  On `L=Z^2`, let

```text
J=[[0,-1],[1,0]],                  J^2=-I.              (9a)
```

For `v=(1,0)`, use `J` on both allocation edges.  The same-atom
vector-valued square is

```text
(v,Jv,Jv,J^2v),
(I-J)^2v=-2Jv=(0,-2) !=0.                              (9b)
```

The scalar readout `(x,y)->x+y` gives the four nonzero values

```text
(1,1,1,-1),                         mu=-2.              (9c)
```

Equivalently, after adjoining `i`, the rank-one square
`(1,i,i,-1)` has mixed face `-2i`.  This is a correctly typed
non-idempotent same-base-atom coefficient construction.  It is intrinsically
signed: `J^2v=-v`, so no nonzero pointed cone is preserved.  Canon has not
attached this `J` to the THM-2791 allocation edges, and the conditional
quarter-turn holonomy elsewhere in the repo is not a transported edge action
on this sheet.  Thus it is a positive algebraic model and a sharp
nonnegativity boundary, not the physical lift.

### 4.2 THM-2593 gives a large exact coefficient positive control

THM-2593's thirteen units `Y_q` define the actual coefficient transports

```text
u_(q,a)=Y_(q+a)Y_q^(-1) in F_13[z]/(Phi_7).             (9d)
```

An independent exact census now gives:

```text
156/156 nontrivial edges have (1-u)^2 !=0;
152/156 have 1-u a unit;
all 156^2=24,336 ordered products (1-u)(1-v) are nonzero;
23,104/24,336 of those products are units.               (9e)
```

The four nonunit edge contrasts are

```text
(q,a,q+a)=(9,5,1),(12,5,4),(1,8,9),(4,8,12).           (9f)
```

The factorization

```text
Phi_7=(z^2+3z+1)(z^2+5z+1)(z^2+6z+1)
```

gives an independent CRT audit.  All `152` unit contrasts are nonzero in
all three quadratic components.  Each of the four exceptions vanishes only
in the `z^2+6z+1` component and survives in the other two.  Hence two
exceptional contrasts still cannot annihilate one another; this explains
structurally why every ordered product in (9e) is nonzero, and why the unit
products are exactly the `152^2=23,104` unit-by-unit pairs.

So the non-idempotent local-system mechanism is not merely hypothetical:
the current coefficient atlas supplies abundant exact mixed amplitudes.
But it also exhibits the obstruction perfectly:

```text
u_(q,a)=Y_(q+a)Y_q^(-1)
```

is a vertex coboundary.  THM-2593's gauge `Y_q^(-1)` sends every edge
transport to one.  Its endpoints are different target sheets, and the
theorem explicitly withholds a common physical carrier/root-target torsor
identification.  Hence (9e) is a nonzero comparison in a chosen global
coefficient trivialization, not a raw mixed face on one physical atom.

## 5. The THM-2806 Rees profile is canonically flat

The central square has

```text
v=(13^2,13,13,1)w,
valuation(v/w)=(2,1,1,0).                              (10)
```

Its valuation curvature is

```text
2-1-1+0=0.                                             (11)
```

After dividing each vertex by its intrinsic power of thirteen, all four
leading units are `wbar`; after raising all vertices to a common grade, all
four are `13^2 w`.  Hence every canonical common-grade comparison is flat.
The ungraded Mobius value

```text
(13-1)^2 w =144w congruent w mod13                    (12)
```

is the tensor product of two one-dimensional orbit-cardinality boundaries,
not a two-dimensional Rees/Bockstein curvature.  Any proposed "double
Bockstein" must distinguish (10) from an arbitrary separable valuation
square.  Canonical associated grading does not.

## 6. The first edge-rephasing-invariant alternative is a joint cocycle

Assume all four amplitudes are nonzero.  Write a square in the order

```text
v=(v00,v10,v01,v11) in (F_13^x)^4.
```

Row/column rephasing acts by

```text
(r0,r1;c0,c1).v
 =(r0 c0 v00, r1 c0 v10, r0 c1 v01, r1 c1 v11).       (13)
```

The simultaneous replacement
`(r0,r1;c0,c1)->(t r0,t r1;t^(-1)c0,t^(-1)c1)` is the
one-dimensional kernel, so the effective gauge group has size `12^3`.
The invariant joint datum is the cross-ratio

```text
kappa=(v00 v11)/(v10 v01).                              (14)
```

It is also complete.  Fixing the kernel by `r0=1`, the unique normalizing
factors are

```text
c0=v00^(-1),  c1=v01^(-1),  r1=v00 v10^(-1),           (15)
```

and they send every nonzero square to the unique representative

```text
(1,1,1,kappa).                                          (16)
```

Conversely the `12^3` anchored choices `(r1,c0,c1)` applied to (16) are
distinct and exhaust its cross-ratio fibre.  Thus the `12^4=20,736`
nonzero `F_13` squares split into exactly twelve gauge orbits, one for each
`kappa in F_13^x`, and every orbit has `12^3=1,728` elements.  The exact
enumeration reproduces all twelve counts.

This classification also forces a subtle correction: the ordinary additive
face is **not** a projective invariant.  The two squares

```text
(1,1,1,1),                    mu=0,
(1,2,3,6),                    mu=2 mod 13,              (17)
```

have the same `kappa=1` and lie in the same row/column orbit.  Indeed the
exact `kappa=1` orbit realizes every residue as `mu`: residue zero occurs
`276` times and each nonzero residue occurs `121` times.  Therefore the
formula

```text
mu=kappa-1                                                (18)
```

holds only in the anchored normalized chart (16); it does not descend from
four projective vertex lines.  This is precisely the difference between the
two branches:

- Branch A has an independently fixed physical common-line normalization,
  so the additive face of a `kappa=1` Segre square is meaningful and may be
  nonzero.
- Branch B has only projective vertex lines.  Its intrinsic target is
  `kappa`, not `mu`; a physical trivialization is still needed before an
  additive face can be read.

Independent rank-one toggles have `kappa=1`.  A genuinely mixed local
system has, in the anchored normalized chart,

```text
(v00,v10,v01,v11)=(w,w,w,kappa w),
mu=(kappa-1)w.                                         (19)
```

Thus, **when only the four vertex lines modulo independent row/column
rephasing are given**, `kappa!=1` is the cheapest intrinsic phase escape.
It is not itself a gauge-invariant additive face; it is square holonomy/a
two-carrier interaction.  It is an alternative to, not a prerequisite for,
the physically normalized `kappa=1` Segre mechanism in Branch A.

THM-2779 supplies the exact abstract model.  On the thirteen-dimensional
coefficient fibre,

```text
T e_r=e_(r+1),         M e_r=zeta^(-r)e_r,
T M=zeta M T.                                            (20)
```

The central phase `kappa=zeta` gives (19).  No one-dimensional
atom-preserving pair of scalar operators can satisfy (20) for
`zeta!=1`, because scalars commute.  Thus a genuine cocycle requires an
internal coefficient fibre (or THM-2779's 169-point faithful permutation
carrier).  Current canon has neither on one THM-2791 physical ancestry atom.

There is a complete two-dimensional frame census behind this model.  For
ordered vectors `u,v in F_13^2`, put

```text
kappa_H(u,v)=zeta^det(u,v) in mu_13 subset F_53^x.       (21)
```

Nontrivial holonomy is equivalent to `(u,v)` being an ordered basis.  There
are

```text
(13^2-1)(13^2-13)=168*156=26,208                        (22)
```

such frames: choose `u!=0`, then choose `v` outside its 13-point span.  For
each fixed nonzero determinant `d`, the equation `det(u,v)=d` has thirteen
solutions in `v` for every nonzero `u`, so each of the twelve nontrivial
phases `zeta^d` occurs exactly

```text
168*13=2,184                                             (23)
```

times.  Swapping the ordered frame negates the determinant and hence sends

```text
kappa_H(u,v) -> kappa_H(v,u)=kappa_H(u,v)^(-1).          (24)
```

The dependency-free enumeration checks all `26,208` swaps and all twelve
uniform phase fibres.  This makes the Heisenberg proposal exhaustive at the
`F_13^2` frame level: there is no exceptional nondegenerate direction, but
there is still no proved map from these abstract frames to two lawful
operations on one THM-2791 ancestry atom.

This sharpens the next target:

```text
source:   one full physical ancestry atom with a retained internal fibre;
map:      two lawful carrier operations with projective commutator kappa;
target:   intrinsic holonomy kappa, then an additive face only after a
          physical coefficient-line normalization;
preserve: full address, endpoint origin, clock, owner, and representative gauge;
loss:     none before forming the face;
test:     exhibit path-dependent H on the same atom and verify kappa!=1.
```

## 7. THM-2807 has no phase holonomy

The positive address triangle has edges

```text
13,                  169*4079,                  689364,
13+169*4079=689364.                                  (25)
```

In formal `(Z1,Z2)` exponents these are

```text
(1,0),                (0,4079),                 (1,4079). (26)
```

Every character of the honest translation action therefore has boundary
holonomy one.  The residue `4079=10 mod13` survives only as the vertical
edge forgotten modulo `169`; it is not curvature of the full triangle.
Treating the two quotient lifts as one atom while retaining phase `10`
simultaneously forgets and uses the same full-depth coordinate.

The thirteen affine lifts in THM-2807 do not change this conclusion: they
are address conjugacies and are not proved factor-covariant operations on
the selected physical cylinders.

## 8. The two THM-2771 cofiber copies

The equal right-cofiber pieces have full-address shifts from the THM-2806
common interval

```text
-13,                  +1248=96*13,                       (27)
```

so their depth-one residues are `-1=12` and `5` in `F_13`.  They are
literally disjoint intervals.  Their phase-weighted copy polynomial is

```text
z^12+z^5
 =2+4 epsilon+...             in F_13[epsilon]/(epsilon^13), (28)
```

and multiplying by the intrinsic copy Bockstein `9` begins

```text
5+10 epsilon+... .                                      (29)
```

The constant `5=9+9` is the observed total Bockstein.  Therefore the two-copy
factor is a unit germ, not an augmentation uniformizer or a hidden mixed
face.

More importantly, exact physical typing blocks a same-atom interpretation:

- both copies fail native `E3`; the `+1248` copy also fails native `c2`;
- both have empty source carrier-twist support, while the common interval
  has source mask `delta_0`;
- the endpoint masks are not a common translated copy;
- their intervals and full addresses are distinct, despite identical
  THM-2584 ancestry label sets and equal scalar coefficient.

Their phase sum is a lawful coefficient on a direct sum of two disjoint
atoms.  Identifying them with the common interval requires quotienting the
full address and the first native-factor ownership that distinguishes them.

## 9. Honest frontier

The explored candidates separate cleanly:

| candidate | algebraic face | why it is not the needed physical face |
|---|---:|---|
| idempotent occupancy | zero on common support | universal theorem |
| positive/Reynolds average | can be nonzero | provenance stays both-absent |
| allocation sign | `4w` | Fourier relabel of flat trivial mass |
| quarter-turn `J` | `-2Jv` | valid signed model; no pointed cone or physical edge action |
| THM-2593 unit atlas | all `24,336` ordered faces nonzero | vertex coboundary between distinct target sheets |
| carrier edge character | nonzero off origin | fixed sheet lives at origin; no common-line gauge descent |
| THM-2806 Rees unit | `w mod13` | canonical leading square is flat |
| THM-2807 vertical phase | residue `10` | honest triangle has holonomy one |
| THM-2771 two-copy phase | nonzero unit germ | disjoint addresses; native-factor/source-mask failure |
| Heisenberg joint cocycle | `(zeta-1)w` | correct abstract escape; physical same-ancestry fibre absent |

The cheapest additive escape is an independently normalized physical
constant phase on each allocation edge; the quarter-turn is its minimal
integral signed model.  The cheapest escape invariant under independent edge
rephasing is one same-ancestry internal fibre carrying a nontrivial joint
cocycle/projective commutator, matching THM-2779.  Neither structure is
currently attached to the THM-2791/2806 physical sheet.

Reproduction:

```bash
python3 .scratch/lrc_nonidempotent_twist_20260728/probe.py
python3 -O .scratch/lrc_nonidempotent_twist_20260728/probe.py
```
