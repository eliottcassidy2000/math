---
id: THM-3068
title: "C3 escape inverse-pole ledger and reciprocal-cofactor shift"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A fixed identity
  sheet plus a punctured cubic Keller component gives an exact split local
  degree-four inverse pair with tame C3 escape.  Its polynomial eliminant has
  leading coefficient t, specialized degree one, and the exact full pole
  a_1=-1/t required by THM-2621.  Its supplied companion satisfies
  q=f_t b_u-f_u b_t=f_T modulo f; every branch Liouville form is exact, the
  trace form is exact, and every Laurent residue vanishes.  Nevertheless the
  integral reciprocal primitive has physical cofactor valuation -2.  The
  reconciliation is the exact quartic law D_Z=-P_z z_i^2 D_X and
  c_Z=-c_X/(P_z z_i^2): on the C3 orbit the missing reciprocal-coordinate
  factor has valuation -5 and moves w(q_X^-1)=3 to w(c_Z)=-2.  Hence the
  coefficient-pole and currently proved Laurent-exactness ledgers do not
  imply a mild reciprocal cofactor or obstruct the Keller congruence.  The
  hostile source is A2 disjoint-union (G_m x A1), not A2; its cubic target
  coordinate is Laurent and it has a nonconstant unit.  No polynomial Keller
  map or C3 component is constructed or excluded.
source: codex-jc-resolvent-bridge-2026-08-01
audit: >
  An independent audit rederived the punctured map's unit Jacobian, the split
  idempotent and branch companion, the inverse-spectral identity
  q_X/f_X'=J_inverse, all Liouville and trace identities, and the quartic
  reciprocal derivative sign and root-product factor.  It checked the fixed
  and C3 valuations directly and replayed normal and optimized companions
  byte-identically to the stored transcript with both declared LF hashes.
  The audit confirms that the source is disconnected and Laurent, so no
  polynomial A2 Keller realization is asserted.
depends_on:
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
  - THM-3057-tame-quartic-inertia-clutch-index-resonance
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
related:
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
script: 04-computation/c3_escape_reciprocal_cofactor_shift_thm3068.py
output: 05-knowledge/results/c3_escape_reciprocal_cofactor_shift_thm3068.out
script_sha256: 97072cec79f764752364557fb69fd4bbd7dd9b34e51592f3c66908736682dec5
output_sha256: 884dec30edb43d68149038bf857280d35eea626e4913896d893c2bdb951bfa26
hash_basis: LF-normalized bytes
---

# THM-3068 -- the reciprocal coordinate carries a five-step valuation shift

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The proposed implication and its type error

[THM-2621](THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger.md)
retains an affine primitive source coordinate `x`, its inverse quartic `f_X`,
the companion `y=b(x)`, and

```text
q_X=(f_X)_v b_u-(f_X)_u b_v,
J_phys=(f_X)'/q_X.                                      (1)
```

[THM-3064](THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary.md)
shows that an **integral** primitive on a tame `C3` orbit has derivative
generating the different and Keller cofactor generating the inverse
different.  Its ramified valuation is therefore `-2`.  The tempting route
is to use THM-2621's coefficient pole to prove

```text
w(q_X^(-1))>=-1                                         (2)
```

and read `(2)` as an exclusion of the inverse-different pole.  This theorem
proves that the reading is false.  The affine primitive `x` escapes, whereas
the integral primitive is its reciprocal `z=1/x`.  Their cofactors differ by
a load-bearing root-product Jacobian of valuation `-5`.

The hostile below also tests the stronger hope that THM-2621's coefficient
pole and polynomial-realization residue ledgers obstruct

```text
q_X=(q_X(a)/d)(f_X)'.                                   (3)
```

They do not: the hostile satisfies `(3)` with scalar one and satisfies every
presently proved branchwise and traced exactness condition.

## 2. A fixed sheet plus the punctured cubic map

Work in characteristic zero with target coordinates `(u,t)`.  On two source
components take

```text
fixed component A^2:       (x,y) |-> (u,t)=(x,y),

cubic component G_m x A^1:
                            (x,y) |->
                            (u,t)=(x^4 y/3,x^(-3)).      (4)
```

Both Jacobians, in the displayed source and target orders, are one.  The
cubic inverse branches over `t=s^3` are

```text
x=s^(-1),                    y=3u s^4.                  (5)
```

Thus the generic fibre has one fixed sheet and one transitive cubic orbit.
The cubic component is an honest algebraic morphism from `G_m x A^1`, where
`x^(-1)` is regular.  It is not a polynomial map from `A^2`: the target
coordinate `t=x^(-3)` is Laurent, and the source has the nonconstant unit
`x`.  The disjoint source is also not connected `A^2`.  These scope failures
are the point of the hostile.

Put

```text
R=Q(u)[[t]],                  K=Q(u)((t)),
f_X(T)=(T-u)(T^3-t^(-1)).                                  (6)
```

The polynomial eliminant with the same generic roots is

```text
R_X(T)=t f_X(T)=tT^4-utT^3-T+u.                         (7)
```

At `t=0`,

```text
R_X(T) mod t=-T+u,                  k_D=1.              (8)
```

Exactly one sheet survives and three escape.  The normalized monic
coefficients are

```text
a_3=-u,                 a_2=0,
a_1=-t^(-1),            a_0=u t^(-1).                  (9)
```

Hence THM-2621's load-bearing coefficient `a_(k_D)=a_1` has the exact full
base pole

```text
v_t(a_1)=-v_t([T^4]R_X)=-1.                            (10)
```

This is exactly the three-sheet escape pole ledger that occurs at a
`k_D=1` Jelonek component, not the finite-root collision used in THM-3064's
first split control.  The later scope audit explains why it is only a local
hostile to that ledger, not a realized polynomial Jelonek component.

## 3. Exact companion and inverse-spectral congruence

The fixed/cubic cross-resultant in the integral equation is

```text
t u^3-1 in R*,                                          (11)
```

so the completed split order contains the idempotent

```text
e(T)=(tT^3-1)/(tu^3-1),
e=(1,0) on fixed | cubic.                               (12)
```

Define the supplied companion

```text
b(T)=t e(T)+3u t^2 T^2(1-e(T)).                        (13)
```

On the two factors,

```text
b(u)=t,
b(s^(-1))=3u s^4,                                      (14)
```

exactly recovering `(4)--(5)`.  Coefficientwise differentiation at fixed
`T` gives

```text
q_X=(f_X)_t b_u-(f_X)_u b_t ==(f_X)'       modulo f_X. (15)
```

There are two proofs.  Direct reduction of `(13)` gives `(15)`.  More
structurally, THM-2621 identifies `q_X/(f_X)'` with the inverse Jacobian;
the fixed branch and `(4)` both have inverse Jacobian one, so the equality
holds on both coprime factors.

In the notation of THM-3064, let

```text
g_X=T^3-t^(-1),                  d_X=g_X(u).
```

The pointed decoder numerator is exactly

```text
d_X(T-u)g_X'-d_X(f_X)'=-d_X g_X,                       (16)
```

while

```text
Res(g_X,(T-u)g_X')=-27(tu^3-1)/t^3 !=0.                (17)
```

Here we have replaced `q_X mod f_X` by its unique reduced representative
`(f_X)'`, which has degree three.  Thus the shifted resultant vanishes with
every denominator gate live, and THM-3064 recovers `(15)` with `kappa=1`.
The coefficient-pole ledger does not obstruct the supplied-pair Keller
congruence.

## 4. The strongest existing Laurent exactness gates also survive

On the fixed branch,

```text
x dy-u dt=0.                                            (18)
```

On a cubic branch `(5)`, `t=s^3` gives

```text
dx wedge dy=du wedge dt,

x dy-u dt
 =3t du+3u dt
 =d(3ut)=d(xy).                                         (19)
```

Consequently:

1. THM-2621's branchwise Laurent residue of `x dy` is zero;
2. its branch potential is not merely rational but the polynomial `xy` in
   the punctured source coordinates;
3. summing one fixed branch and three cubic branches gives

```text
omega_F=Tr(x dy)-4u dt=d(9ut),                          (20)
```

so the trace--Liouville form is exact, regular, supported at no spurious
divisor, and has zero residue.

The Laurent packet is maximally sparse:

```text
x=s^(-1),                         y=3u s^4,
dx wedge dy=3s^2 du wedge ds=du wedge d(s^3).           (21)
```

Thus no hidden cancellation is being engineered.  The sole coefficients
`x_(-1)=1` and `y_4=3u` already provide the required tame `C3` Jacobian.
Equations `(18)--(21)` show that the coefficient pole, branchwise residue,
localized potential, and traced exactness ledgers can all coexist with
`(15)`.

## 5. Exact reciprocal derivative and cofactor law

Let the four nonzero affine-primitive roots be `x_i`, put

```text
z_i=x_i^(-1),                     P_z=product_i z_i,
D_X,i=product_(j!=i)(x_i-x_j),
D_Z,i=product_(j!=i)(z_i-z_j).                         (22)
```

Since

```text
x_i-x_j=-(z_i-z_j)/(z_i z_j),                          (23)
```

the three factors in the quartic derivative give

```text
D_X,i=-D_Z,i/(P_z z_i^2),
D_Z,i=-P_z z_i^2 D_X,i.                                (24)
```

Let `c_X,i,c_Z,i` denote the cofactors which convert the two derivative
stars into the **same original physical Jacobian**:

```text
J_i=c_X,iD_X,i=c_Z,iD_Z,i.                             (25)
```

Then

```text
c_Z,i=-c_X,i/(P_z z_i^2).                              (26)
```

For the THM-2621 affine primitive, `(1)` gives `c_X=q_X^(-1)`.  Formula
`(26)` proves that `q_X^(-1)` is not the reciprocal primitive cofactor.  The
destroyed coordinate is precisely the reciprocal root product and the
square of the marked root.

## 6. The exact five-step valuation shift

Normalize the cubic ramified valuation by

```text
w(s)=1,                         w(t)=3.                 (27)
```

The reciprocal quartic is the integral polynomial

```text
f_Z(Z)=(Z-u^(-1))(Z^3-t),       P_z=t/u.               (28)
```

On a cubic branch `z=s`, the exact valuation invoice is

| quantity | formula | `w` |
|---|---|---:|
| affine primitive | `x=s^(-1)` | `-1` |
| integral reciprocal | `z=s` | `1` |
| reciprocal root product | `P_z=t/u` | `3` |
| affine derivative | `D_X=3(s^(-1)-u)s^(-2)` | `-3` |
| inverse numerator | `q_X=D_X` | `-3` |
| affine cofactor | `c_X=q_X^(-1)` | `3` |
| integral derivative | `D_Z=3s^2(s-u^(-1))` | `2` |
| integral cofactor | `c_Z=D_Z^(-1)` | `-2` |

The factor in `(26)` has valuation

```text
w((P_z z^2)^(-1))=-3-2=-5,                             (29)
```

and indeed

```text
w(c_Z)=w(c_X)-5=3-5=-2.                                (30)
```

The value `2` for `D_Z` is the tame cubic different exponent; `c_Z` is its
inverse-different generator exactly as THM-3064 requires.  On the fixed
sheet, `w(D_X)=-3` but `w(D_Z)=0`, so `(24)` also restores the separated unit
derivative there.

Thus `(2)` is true in this hostile, with room to spare, while the desired
conclusion for the integral reciprocal cofactor is false.  Any future
valuation argument must transport the reciprocal factor `(P_z z_i^2)^(-1)`
rather than compare the two cofactors directly.

## 7. Exact stopping boundary

The hostile proves that the following local data do **not** exclude a tame
`C3` Keller principal part:

```text
specialized fibre degree k_D=1;
exact full inverse-coefficient pole a_1=-1/t;
fixed/cubic split idempotent and unit cross-resultant;
q_X=(q_X(u)/d_X)(f_X)';
branchwise symplectic form and zero Laurent residue;
polynomial branch potential and exact trace--Liouville form;
integral reciprocal different/inverse-different valuation (2,-2).       (31)
```

The last line of `(31)` is reached only after restoring `(26)`.  What the
hostile deliberately drops is global polynomial realization:

1. the source is `A^2 disjoint-union (G_m x A^1)`, not connected `A^2`;
2. the cubic target coordinate `t=x^(-3)` is Laurent, not polynomial on
   `A^2`; and
3. the cubic source has the nonconstant unit `x`, exactly the boundary in
   THM-2633's no-omitted-divisor mechanism.

Therefore the next viable `C3` exclusion must use a constraint which sees
polynomial target regularity or the constant-unit/connected affine source.
Another coefficient valuation, residue-zero equation, or untransported
bound on `q_X^(-1)` cannot suffice.

```text
PROVED HERE:       exact k_D=1 eliminant and full a_1 pole;
                   exact split companion and q_X=(f_X)' congruence;
                   branchwise and traced Liouville exactness;
                   quartic reciprocal derivative/cofactor law;
                   exact C3 valuation shift 3 -> -2;
                   punctured fixed-plus-C3 hostile to the local ledgers.

REFUTED HERE:      using w(q_X^(-1))>=-1 as a bound on c_Z;
                   excluding q_X proportional to (f_X)' from the current
                   coefficient-pole and Laurent-exactness ledgers.

NOT PROVED:        globalization of the hostile to a connected quartic;
                   polynomial-map realization on A^2;
                   existence or exclusion of a C3 Jelonek component,
                   A4, S4, G1, JC(2), or DC(2).                         (32)
```

## 8. Exact companion

Run

```text
python3 04-computation/c3_escape_reciprocal_cofactor_shift_thm3068.py
python3 -O 04-computation/c3_escape_reciprocal_cofactor_shift_thm3068.py
```

Both executions must LF-byte-match the stored transcript.  The companion
checks the eliminant and specialized degree, coefficient pole, split
idempotent and companion, coefficientwise inverse-spectral congruence,
shifted decoder denominator, punctured rational Jacobian, branchwise wedge
and Liouville potential, trace exactness, generic four-root reciprocal law,
fixed and cubic specializations, all declared valuations, and the tame cubic
discriminant.  Every truth-bearing check uses explicit runtime exceptions
rather than Python assertions.
