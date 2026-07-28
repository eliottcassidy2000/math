---
id: THM-2751
title: "Root-zero clutch Mayer--Vietoris wing shear"
status: >
  RESERVED / REFUTED PROVISIONAL CANDIDATE.  The legacy natural-sheet
  constructor used below inserted E3 and the four target-sheet gates but
  omitted the source-one clock factor d_(1,ell).  It therefore cannot be
  compared with THM-2749's fully marked common section.  The internally exact
  12/2/7 calculation belongs to that clock-blind hostile carrier, not to the
  claimed physical Mayer--Vietoris decomposition.  No result may depend on
  this file.  The corrected clocked computation has coefficient-null left
  wing and a one-sided charged right wing; see MISTAKE-313.
source: root/root-zero-clutch-mayer-vietoris-2026-07-28
depends_on: []
related:
  - THM-830-b3-deletion-deck-mirror-current-calculus
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
  - THM-2750-arm-blind-clutch-no-go-and-minimal-marked-leakage
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
script: 04-computation/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py
output: 05-knowledge/results/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.out
script_sha256: b84d1d1a2942a3470cba86b7ab9d31c62a21451cd674977553cd008fa49db780
output_sha256: 92af039314386bd55d0d6e0854c13c1d03da86966d1fc913fc2751421b12520d
hash_basis: LF-normalized bytes
---

# THM-2751 -- the natural root-zero shear is a one-sided wing boundary

**RESERVED / REFUTED PROVISIONAL CANDIDATE.**

## Retraction boundary

The construction called `A` below is not the fully marked natural source
sheet.  Its helper `restrict_e3_and_sheet` never intersects the carrier with
the source-one clock comb `d_(1,ell)`, whereas THM-2749's common section uses
`source_present_section(...,source_clock=1,...)`.  The two delayed prefix banks
are identical; this missing factor is the first failed implication.

Consequently the weighted decompositions and `12/2/7` ratios below are exact
only for a clock-blind hostile carrier and cannot be combined with THM-2749.
They are retained as correction history.  The correctly clocked calculation
has source coefficient equal to the common coefficient, hence a physically
nonempty but coefficient-null left wing; the right wing carries the entire
one-sided shear.  See MISTAKE-313.  Nothing below this banner is a proved
dependency.

THM-2749 identifies an exact signed clutch after the natural root-zero sheet
is intersected with its translated pullback.  The same theorem records that
the unsplit natural source and target coefficients are unequal.  The point of
this candidate is to compare those calculations at the level of the actual
weighted interval functions.  The mismatch is neither numerical noise nor a
subtraction of unrelated tables: it is carried entirely by the two
one-sided complements of the common section.

## 1. The weighted Mayer--Vietoris cover

Use THM-2749's canonical row, translation and open seams

```text
(H,q1,...,q5,c1,c2,c3)
  =(1,14,27,40,53,66,13,13^3,2*13^5),

tau=7/13^6,
Sigma_-=(169,181)/182,
Sigma_+=(1,13)/182.                                    (1)
```

Fix target label `t=3`, rail `8`, the full `E3` source condition, lawful
sheet `U_(0,3)`, both relative-present factors, and the actual delayed
terminal fork `D^6 Q_(3,{1,2})`.  For each delayed clock `ell`, let

```text
A_ell = the natural source-coordinate E3/U_(0,3) weighted carrier,
B_ell = the natural forward target-coordinate E3/U_(0,3) weighted carrier.
                                                                  (2)
```

The source uses carry `12` and right root `12`; the target uses carry `6`
and lawful left root `1`.  Let `M^-_ell` be the source-coordinate common
carrier of THM-2749 and let

```text
M^+_ell=T_tau(M^-_ell)                                  (3)
```

be its direct forward-coordinate target copy.  Thus the symbol `M` below
means the one common physical section represented in the appropriate chart;
it does not identify arbitrary points in the complements.  Put

```text
L_ell=A_ell\M^-_ell,
R_ell=B_ell\M^+_ell.                                   (4)
```

The exact companion constructs all six profiles directly.  In every clock
their canonical weighted-piece counts are

```text
                     whole       common       wing
source A               1275         239       1036
target B               1275         239       1036,       (5)
```

and their weighted masses are

```text
mu_w(A)=929934280541992017600,
mu_w(M^-)=174321777293450297280,
mu_w(L)=755612503248541720320,

mu_w(B)=929934304688494607040,
mu_w(M^+)=174321777293450297280,
mu_w(R)=755612527395044309760.                           (6)
```

Before any coefficient evaluation, canonical interval comparison proves

```text
A_ell=M^-_ell disjoint-union L_ell,
B_ell=M^+_ell disjoint-union R_ell                      (7)
```

as weighted step functions, for all seven clocks.  In particular, `(7)` is
an additive Cech/Mayer--Vietoris decomposition, not a formal difference of
the two final vectors.

## 2. Exact additivity of the delayed functional

The natural `Q_(3,{1,2})` prefix and THM-2749's fully marked prefix agree as
exact interval prefixes for all seven clocks.  Let `S` denote the source
carry-`12` functional and `T` the target carry-`6` functional, using that
common prefix bank.  Since the delayed integral is additive on disjoint
weighted interval unions, `(7)` gives

```text
S(A)=S(M^-)+S(L),
T(B)=T(M^+)+T(R).                                      (8)
```

Direct exact evaluation of every summand gives

```text
S(A)=(0,a,a,a,a,a,a),
a=1812281403506324508838080,

T(B)=(0,b,b,b,b,b,b),
b=1826551436254490256030720,                            (9)

S(M^-)=T(M^+)=(0,C,C,C,C,C,C),
C=339633525654239542165440,                            (10)

S(L)=(0,l,l,l,l,l,l),
l=1472647877852084966672640,

T(R)=(0,r,r,r,r,r,r),
r=1486917910600250713865280.                           (11)
```

Equations `(8)`--`(11)` are checked twice: first as the actual weighted-piece
unions and then as the resulting integer coefficient identities.

## 3. The graded transition and the folded natural ratio

Divide by the inherited content `26`, normalize by the indicated nonzero
root in `F_13`, and reduce in `F_13[z]/(Phi_7)`.  The six classes are
constants, with multiplication determinants

```text
piece                     constant mod13       determinant
common source, root 12           9                  1
common target, root 1            4                  1
left wing, root 12               9                  1
right wing, root 1               5                 12
natural source, root 12          5                 12
natural target, root 1           9                  1.    (12)
```

Consequently

```text
common clutch gain       =4/9=12=-1,
formal wing ratio        =5/9=2,
folded natural ratio     =9/5=7             in F_13.       (13)
```

The exact coefficient grading is therefore

```text
(M,L)_source=(9,9),
(M,R)_target=(4,5)=(12*9,2*9).                         (13a)
```

This is an exact coefficient-diagonal **comparison** after naming separate
source and target summands, not an operator on one physical carrier.  Its
second entry is only a formal normalized wing ratio, because no physical
`L -> R` identification has been supplied.
No single scalar intertwines both additive summands in `(7)`: the common
section requires `12`, whereas the formally paired wing endpoints require
`2`.
After the non-equivariant fold `(x,y)->x+y`, the particular source vector
becomes `5` and the target becomes `9`, producing the effective ratio `7`.
Thus `7` belongs only to this unsplit folded vector.  It is not a third local
transition, and the fold does not conjugate `diag(12,2)` to a scalar map.

## 4. The exact wing boundary current

Subtracting the source functional from the target functional and using
`(8)` gives the exact identity

```text
T(B)-S(A)=T(R)-S(L)
          =(0,d,d,d,d,d,d),

d=14270032748165747192640,
v_13(d)=1,
(d/26) mod13=12.                                      (14)
```

Thus the natural single-sheet defect is supported entirely on the one-sided
wings.  The two-sided overlap itself has no raw shear: it is the signed
root-normalized mirror clutch of THM-2749.

This is the sharp positive conclusion and also the sharp categorical
boundary.  Equations `(4)` and `(7)` supply an additive common-section
cover, but neither THM-2744 nor THM-2749 gives a canonical physical map

```text
                         L_ell ---> R_ell               (15)
```

or a relative phase between the two wing functionals.  Therefore `(14)` is
a signed boundary current, not yet a physical cross-wing current.  A lawful
pairing, mapping-cone reference, or endpoint selector is still required.

## 5. Internal chart strata are not external C3 arms

There is a tempting but invalid relabelling.  If a future physical selector
assigned the three gains

```text
g=(12,2,2)                                             (16)
```

to three externally rotated equal arms, then over `F_13`, with primitive
cube root `omega=3`, its normalized Fourier profile would be

```text
g_hat=(1,12,12),
Pi_3 g=(1,1,1),
Q_3 g=(11,1,1).                                       (17)
```

Such a selector would charge both nontrivial `C3` modes.  Nothing in this
candidate constructs it.  The pieces `M,L,R` are internal source/target
chart strata.  If they are formally replicated in the same way on every
external arm, the only tautological external typing is

```text
                         I_arm tensor C_internal.       (18)
```

It commutes with external arm rotation, so directly

```text
Q_3 (I_arm tensor C_internal) Pi_3=0.                  (19)
```

This is an algebraic arm-blind boundary, not an existing physical operator:
neither the external carrier nor the internal `L -> R` wing map has been
constructed.  Formula `(19)` is proved independently by `Q_3 Pi_3=0` and is
the boundary proved abstractly by THM-2750; THM-2750 remains a related result,
not a dependency.  In particular, the natural ratio `7` is an
aggregate sheet scalar, not the external Fourier mean `1` in `(17)`.

There is one useful **conditional diagnostic**, not an LRC construction.  If
future work supplied one common external three-arm carrier, a lawful armwise
realization/pairing of the formal common/wing diagonal rule, and source totals

```text
S_i=M_i+L_i
```

were arm-invariant, then the internal diagonal rule `(13a)` would give

```text
T_i=12M_i+2L_i=2S_i+10M_i,
Q_3 T=10 Q_3 M.                                        (20)
```

Hence external charge would be exactly the variation of the common-section
share among the arms.  The present theorem supplies neither a common lattice
nor such a carrier, so `(20)` is a scoped target for a selector theorem, not
a consequence about the existing physical packet.

The same conditional statement has an exact THM-2348 rectangle form.  Regard
`(M_i,L_i)` as a formal `2 x 3` table.  Under the arm-invariant-total
hypothesis `M_i+L_i=S`, its rectangle difference for arms `i,j` is

```text
Delta_ij=M_i+L_j-M_j-L_i=2(M_i-M_j).                   (21)
```

Because `2` and `10` are units in `F_13`, all mixed rectangles vanish iff
the `M_i` are constant iff `Q_3T=0`; a nonzero rectangle is exactly the
arm-dependent common-share interaction detected by `(20)`.  The companion
exhausts all `13^4=28561` such formal tables.  This is a conditional ANOVA
corollary analogous to THM-2348, not evidence that the required physical
two-by-three table already exists.

## 6. Connection contract and scope

| item | exact content |
|---|---|
| source | natural source carrier `A`, carry `12`, right root `12` |
| target | natural forward target carrier `B`, carry `6`, left root `1` |
| map | common-section translation `M^- -> M^+` plus the additive decompositions `(7)` |
| preserved | rail `8`, `E3`, lawful `U_(0,3)`, relative present, terminal fork, prefix bank, open seam, and raw additivity |
| destroyed or absent | canonical `L -> R` point map, wing relative phase, external-arm label, owner endpoint |
| first failure | the common clutch gain is `12`, while the formal wing endpoint ratio is `2`, so one scalar cannot transport both |
| strongest survivor | exact signed wing current `(14)` with primitive content residue |
| cheapest next test | a lawful cross-wing pairing through the target/clock character bank, retaining a common endpoint reference |

No external `C3` selector, noncentral transporter, endpoint current,
owner/root provenance, row exclusion, scalar-ledger decrement, or proof of
LRC(14) follows.

## 7. Exact reproduction

Run

```bash
python 04-computation/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py
python -O 04-computation/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py
```

Both modes must byte-match

```text
05-knowledge/results/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.out.
```

The companion pins the THM-2749, root-zero overlap, and natural-sheet
reconstruction scripts by LF-normalized SHA-256.  It verifies the seven
prefix identities, reconstructs all `1275` whole and `239` common weighted
pieces per chart and clock, proves the two `1036`-piece complements are
literal disjoint differences, checks weighted mass and functional
additivity, recomputes the six unit profiles and determinants, verifies the
shear valuation, checks both sides of the external-`C3` boundary, and
exhausts the `28561` conditional ANOVA tables.  No
truth-bearing Python `assert` is used.

`QED` for the provisional theorem; the file remains reserved until
independent promotion.
