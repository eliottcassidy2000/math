# Same-sheet common-carrier allocation audit

**Status:** `FINITE-EXACT` positive central post-carrier-DFT coefficient
face, together with `FINITE-EXACT` pre-Fourier common-atom and nonzero-step
allocation-to-endpoint no-go results.  Scratch only; no tracked theorem is
promoted here.

## Verdict

The selected rail-eight common carrier contains a literal Boolean ancestry
sheet, not merely a unit summand obtained by dividing an aggregate weight.
Keeping that sheet label through the endpoint construction gives a
fixed-origin central carrier-Fourier coefficient square whose mixed
coordinate `D3` is nonzero in both exact-order cyclotomic specializations.

This does **not** close the coefficient-side pre-Fourier common-atom gate in
THM-2772 Section 4.  On the primal pair of carrier-twist variables the four
raw states are

```text
(B,P,Q,H)=w*(1, delta_(a=0), delta_(b=0),
                delta_(a=0) delta_(b=0)).
```

Their sole fourfold co-supported point is `(a,b)=(0,0)`, where the vector is
`(w,w,w,w)` and its raw Möbius face is zero.  The nonzero central `D3=144w`
is instead the sum over the `12*12` complement, where only `B` survives.
Thus the exact result is a transformed marginal K4 co-support and,
simultaneously, a sharp fixed-sheet pre-Fourier common-atom no-go.

The literal present state keeps the outer selector `xi` fixed.  It meets
only carrier twist zero, so present is exactly `1/13` times bare at every
endpoint: the intrinsic endpoint steps are both zero and their determinant
face is zero.

The independently nondegenerate endpoint parallelogram comes instead from
moving `xi` through the carrier orbit.  Exhaustive full-fibre comparison
proves that orbit bank is not a translated/character-gauged bare bank.
Hence it cannot be identified with the physical allocation toggles.  The
map from the proved coefficient `D3` to a determinant/root-deck/Cech
correction remains open.

## 1. Aligned physical object

The computation uses

```text
A = source one-sided carrier,
B = pullback(target one-sided carrier),
C = A intersect B.
```

It never identifies either exclusive wing `A\B` or `B\A` with `C`.
Exhaustion of all `9*9*7=567` lawful cells reproduces `193` nonempty source
wings and `193` nonempty target wings and verifies that opposite wings are
disjoint in every cell.

The same exhaustion gives the nonempty-common-cell census by clock

```text
clock 0,1,2,3,4,5,6:  (0,81,56,56,0,0,0).
```

Thus this entire fixed rail-eight configuration bank has common support on
only three of seven clocks.  Even its central `144w` marginal cannot by
itself pay a seven-clock Cech invoice.

The selected cell is

```text
(s,t,clock)=(0,4,1),
I=[142004992589460,142005019034340),
weight(I)=27581135604.
```

Its target copy is exactly `I+431933040`.  The complete factor order is

```text
(E3, clock, q1, q2, c2, c3).
```

The source atom lies in all six native factors and all six pulled target
factors.  Its translated target lies in all six native factors and all six
pulled source factors.  Thus both simultaneous signatures are twelve
passes and both intrinsic first-failure labels are
`common_after_all_factors`.

The retained semantic states are

```text
source (residue,edge,carry,private-root)=(7,1,12,12),
target (residue,edge,carry,private-root)=(8,0,6,1),
delayed half kappa=1.
```

On unit weight, the delayed values are literally equal:

```text
source carry 12: (0,103478815440),
target carry 6:  (0,103478815440).
```

## 2. Literal Boolean sheet

The THM-2584 transfer profile unfolds into inverse-branch labels

```text
(a,b,e') in (Z/13^5) x (Z/13^2) x (Z/13^5).
```

At the source midpoint and its target translate, the complete labelled sets
agree, not only their cardinalities:

```text
|U_source|=|U_target|=966606,
|V_source|=|V_target|=28534,
source label set = target label set,
966606*28534=27581135604.
```

A concrete retained sheet is

```text
xi=(a,b,e')=(59162,26,56658).
```

All six inverse-branch conditions for this same triple give the joint
source-coordinate chamber

```text
[138281416853580,159555049051860).
```

The selected interval lies strictly inside it, with left/right slacks

```text
3723575735880, 17550030017520.
```

Thus the same labelled branch remains active at every point of `I` and
`I+SHIFT`; no seam or reselected marginal is involved.  The aggregate
weight forgets `(a,b,e')`.  Clock, semantic/carry state, simultaneous
factor signature, and endpoint masks are functions of the projected
coordinate and are constant on `I`; the companion retains the forgotten
triple separately.

This agrees with the independently audited THM-2791 partial-germ result:
the larger contributor chamber preserves all `27581135604` product sheets.
THM-2791 still withholds an endpoint-origin allocation, exactly as the
negative covariance result below requires.

## 3. Direct absent/present constructor

Let `1_xi` be the outer unit-sheet selector and `C_h` the source or target
carrier at twist `h`.  The exact interval containment check gives

```text
1_xi C_0=1_xi
```

on the source and, after translation, on the target.  Therefore the
carrier-absent atom-local raw bank is constructed directly by evaluating
the fixed `1_xi` interval with no carrier-twist variable.  The zero-twist
present raw bank is the same integral by idempotence, not by division of the
aggregate multiplicity.

To compare absent and present states in the harmonic-zero fibre, extend the
omitted variable constantly over the dual `F_13`.  Unnormalized inverse
Fourier then gives

```text
bare[address,k=0]=13*present_raw[address,h=0],
bare[address,k!=0]=0.
```

For the literal present state, retain the same outer selector and multiply
by the thirteen carrier twists.  Exact full-carrier intersections give

```text
1_xi C_0=1_xi,
1_xi C_h=0 for h=1,...,12
```

on both source and target.  Therefore present has only the twist-zero dual
entry and

```text
bare=13*present
```

pointwise on the entire endpoint plane.  Translating `xi` itself with
`C_h` is a different moving-orbit constructor and is kept separate below.

At the fixed pullback origin the complete address is

```text
(r,k,l,L,R)=(0,0,0,(0,0),(0,0)),
q=0, Delta=0.
```

Here `r=0` is the distinguished `e0` section.

## 4. Central fixed-origin carrier-Fourier square

At the fixed origin, the four factors are ordered as

```text
(P_bare(L),P_present(L),Q_bare(R),Q_present(R)).
```

For `p=352341050142921841`:

```text
factors =
(161934023005863791,
 175075409527953449,
 300649489018742460,
 50230041473974177)

products =
(346389504894336138,
 351883238969953710,
 351883238969953710,
 216790045382338969)

(D0,D1,D2,D3) =
(209922877787817004,
 129599459511997169,
 129599459511997169,
 211754122479689528).
```

For `p=956354278959359281`:

```text
factors =
(155560747474790541,
 11966211344214657,
 879469786637801433,
 361914377113479889)

products =
(17884554354397974,
 295638590014756546,
 295638590014756546,
 96307143767239679)

(D0,D1,D2,D3) =
(705468878151150745,
 877931689546517576,
 877931689546517576,
 479268797051483842).
```

Every displayed entry is nonzero.  In each field,

```text
D3=(P_bare-P_present)(Q_bare-Q_present) != 0,
D0*D3-D1*D2=0.
```

More intrinsically, if `b=P_present(L)` and `q=Q_present(R)`, then

```text
factors  =(13b,b,13q,q),
products =bq*(169,13,13,1),
Hadamard =bq*(196,168,168,144),
D3       =144bq != 0.
```

The mechanism is not special to 13.  For any finite twist set/group of
order `m`, if “absent” is constant and “present” is the delta function at
the identity on each side, then the central product profile and its
Hadamard transform are universally

```text
w*(m^2,m,m,1),
w*((m+1)^2,m^2-1,m^2-1,(m-1)^2).
```

The raw Möbius face is zero on both coordinate axes and equals `w` on
`(G\{e})^2`; central nonvanishing is exactly its `(m-1)^2`-point sum.  This
general lemma explains both the positive marginal and the common-atom
failure, rather than treating either as a field-specific accident.

This is the central coefficient obtained after independently summing the
source and target carrier-twist variables.  Before those sums, with
`w=bq`, the exact raw support census is

```text
state                    support count in F_13^2
B (both carriers absent)             169
P (source present)                    13
Q (target present)                    13
H (both present)                       1
all four co-supported                   1
raw Möbius nonzero                    144
```

At the unique all-four point `(0,0)`, the raw vector is `(w,w,w,w)` and
`B-P-Q+H=0`.  The raw Möbius function equals `w` precisely at the `144`
points with `a!=0` and `b!=0`, where its vector is `(w,0,0,0)`.  Hence

```text
central D3 = sum_(a,b) (B-P-Q+H) = 144w.
```

The nonzero `D3` is therefore a positive central transformed-marginal
certificate, not a positive pre-Fourier common-atom certificate.
The Pluecker identity is only rank-one decomposability; its zero value is
not the claimed nonvanishing.

The integral multiplier profile `(169,13,13,1)` has relative 13-adic
filtration degrees `(2,1,1,0)`.  After reduction modulo 13 it becomes the
one-corner profile `(0,0,0,1)`, while `144=12^2` is `1 mod 13`.  This gives
the familiar virtual one-corner hostile a canonical marked-sheet/Rees
filtration provenance.  A two-dimensional Bockstein or associated-graded
root-deck map is a concrete possible sidecar, but no such map is proved
here; forgetting the filtration again loses three vertices.

Both exact-order embeddings evaluate the same atom-indexed cyclotomic
integer constructor, so either nonzero residue already certifies
characteristic-zero nonvanishing.  The second field is an independent
specialization control.

## 5. Moving-orbit endpoint parallelogram and the failed identification

The distinct constructor that moves `xi` with all thirteen carrier twists
is nonzero on all `169` endpoint origins.  Sampling

```text
L=R=(0,0), source step=(0,1), target step=(12,0)
```

gives four nonzero moving-orbit products in both fields and nonzero sample
`D3`.  Their endpoint pairs have

```text
q=(0,0),(0,1),(1,0),(1,1),
Delta=(0,0,0,1),
Delta00-Delta10-Delta01+Delta11
 =det((0,1),(12,0))=1.
```

This is a genuine labelled-orbit endpoint parallelogram, but its four
physical supports are translated copies of `xi`, not the four allocation
flags on fixed `xi`.  The literal fixed-sheet allocation obeys

```text
present(x)=13^(-1) bare(x)
```

with intrinsic source and target endpoint steps both zero.  Its determinant
mixed face is therefore zero.

For the moving-orbit bank the companion searches, on every point of
`F_13^2`, every identity of the form

```text
present(x)=c*zeta^(u dot x)*bare(x+s)
```

for all `169` shifts `s`, all `169` characters `u`, and the forced scalar
`c`.  The answer is empty for source and target in both fields.  Hence no
allowed translation/character/scalar gauge turns either bare array into its
moving-orbit array.  In particular, `(0,1)` and `(12,0)` are not physical
allocation translations.

The independent constructor gives a cheaper obstruction on the inverse
endpoint-address side: the literal absent and present banks each have
support size `81`, while the transported-sheet replacement has support size
`121`, on both sides and in both fields.  Translation and character
modulation preserve support cardinality, so `81!=121` already forbids an
affine covariance before comparing values.  The full search above is the
stronger value-level control.

The translated bare/present quartet is numerically all-nonzero, but it is
only a rejected candidate because the full-fibre typing identity fails.
The moving-orbit quartet is a valid endpoint sample, not a fixed-sheet
bare/present allocation quartet.

Thus the exact faces and no-go boundaries are:

```text
central transformed coefficient face: fixed-xi D3 = 144w != 0;
primal all-four common-atom face:     raw Möbius = 0;
intrinsic geometric face:             zero-step determinant face = 0;
separate orbit label face:            determinant mixed face = 1;
proved boundary:                      orbit is not the allocation bank.
```

## 6. Gauge and arithmetic controls

For each exact field, all `2197` left and all `2197` right moving-orbit dual
cells obey

```text
(ell,a,b)~(ell+W,a+1,b-1).
```

The moving-orbit dual supports are `522/2197` on each side.  Holding the
carrier fixed under the representative change fails in `436` cells on each
side in each field, so the gauge test is non-vacuous.

For the literal mask, the same gauge transports

```text
(xi,h=0) -> (T_(-1/13)xi,h=1)     on the source,
(xi,h=0) -> (T_(-1/13)xi,h=12)    on the target.
```

The delta-mask sum remains one, so its `k=l=0` coefficient is invariant
provided the outer atom selector is transported as a sidecar.  The
companion records the result at the distinguished `e0` section; it does not
claim an allocation current after forgetting that selector.

The fixed-origin central square, translated fixed-sheet sample, and
moving-orbit square are each formed as a rank-one product before their
Hadamard coordinates are computed.  No four independently reselected
unlabelled marginals are used.  The explicit primal census nevertheless
shows that the central square collapses to a bare-only `12*12` complement;
this is why it cannot be typed as the required common-atom face.

## 7. Exact stopping boundary

The new positive result supplies one common Boolean ancestry label and one
nonzero central post-carrier-DFT coefficient face at one clock.  Its primal
all-four raw Möbius face vanishes, its allocation action is scalar rather
than a nonzero endpoint translation, and this fixed configuration bank has
common support on only three clocks.  It does not supply:

- a nonzero pre-Fourier common-atom mixed face;
- an allocation-to-endpoint translation;
- a map from coefficient `D3` to the determinant-sector face;
- seven coherently oriented clock faces;
- the THM-2591 root-deck/Cech functional or invoice;
- a global chronology/semantic arrival current;
- an all-`91` unit, row exclusion, or LRC(14).

The next sharp problem is to find a different allocation object whose
literal primal common atom has nonzero raw Möbius face and nonzero intrinsic
endpoint steps, or to construct a filtration-sensitive Bockstein/root-deck
functional that can extract information from the proved
`(169,13,13,1)` lift without mistaking its bare-only complement for a
common atom.  The moving-orbit shortcut and the direct seven-clock use of
this fixed rail-eight bank are now exactly excluded.

There is also a general Boolean obstruction behind the raw no-go.  For an
atom of weight `w` with indicator values `C_S,C_T in {0,1}`, its pointwise
allocation vector is

```text
(B,P,Q,H)=w*(1,C_S,C_T,C_S*C_T)
```

and therefore

```text
B-P-Q+H=w*(1-C_S)*(1-C_T).
```

Every common-present atom (`C_S=C_T=1`) necessarily has zero raw mixed
face; the mixed face of pure membership is supported exactly where both
carriers are absent.  A successful THM-2772 common-present construction
must therefore add a non-idempotent local amplitude—phase, signed weight,
cocycle, overlapping weighted translate, or a typed translation
groupoid—rather than merely refine Boolean carrier membership.

### Candidate connection to the graded address simplex

The incoming `THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary`
is still a reserved proof-complete candidate, not a proved dependency.  It
supplies exactly the complementary object: a positive full-depth address
two-simplex, a vertical `Z2^4079` lift with first residue `10`, and thirteen
affine lifts agreeing on one residue sheet, but no endpoint/allocation lift.
The present calculation supplies a canonically lifted allocation filtration
`(2,1,1,0)`, but no nonzero raw common-present face.

A precise next experiment is therefore to pull one retained ancestry label
over the three `THM-2807` cylinders and form a relative two-direction
Bockstein/Rees bicomplex:

```text
address direction:     pure edge / graded edge / vertical lift;
allocation direction:  filtration degrees (2,1,1,0);
sidecars retained:      ancestry, semantic owner, endpoint origin, affine-lift index.
```

The cheapest decisive tests are literal sheet survival at all three
vertices, factorwise transport on both edges, and whether the first
Bockstein differential depends on the affine-lift index.  If it does not,
the construction again collapses to the virtual one-corner hostile; if it
does, the lift index is the missing coordinate that could distinguish the
integral allocation square after quotienting.  No such compatibility has
yet been checked or claimed.

## Reproduction

```bash
python3 .scratch/lrc_atom_dependent_four_corner/atom_dependent_four_corner.py
python3 -O .scratch/lrc_atom_dependent_four_corner/atom_dependent_four_corner.py
```

The companion is assertion-free and pins all three imported top-level
dependencies.  Final replay certificate:

```text
script SHA-256:
  8f064ff8922a787211cb66ad7f7a9cdfe5a87ded84100d7dbd39e264e387482f

ordinary stdout:
  .scratch/lrc_atom_dependent_four_corner/
    atom_dependent_four_corner.normal.out

independent optimized stdout:
  .scratch/lrc_atom_dependent_four_corner_audit/
    witness.optimized.out

stdout SHA-256:
  a970776ed95128b5745c1fd370af768778b409d931a57b68006e04271e00813f

size:
  86 lines, 10862 bytes

normal versus optimized:
  byte-identical (cmp exit 0)

terminal line:
  ALL EXACT CHECKS PASSED
```

The ordinary run completed in approximately `205.4 s`.  The independently
launched optimized process was observed at approximately `329 s`; an exact
timer was not captured.  An AST audit finds zero Python `assert` nodes, so
optimized mode cannot silently remove a truth-bearing check.  The
independent hostile companion also passes in both modes and reproduces the
ancestry counts, field arithmetic, raw support census, sole fourfold point,
and 144-point bare-only boundary.
