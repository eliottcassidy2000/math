# The affine-c base locus carries a nonsplit exceptional quadratic

**Status:** FINITE-EXACT + DUAL-DERIVATION + HOSTILE-REDUCTION PARTIAL SCOUT
on the fixed `C=c+x`, `d=k=1` slice in both THM-3212 accessory fields; its
pointed deck and true-gradient obstruction are canonized in
[THM-3309](../01-canon/theorems/THM-3309-exceptional-quadratic-deck-passport-and-gradient-unimodularity-obstruction.md),
while the universal trace/norm passport and independent replay are in
[THM-3312](../01-canon/theorems/THM-3312-exceptional-quadratic-trace-norm-and-cofactor-antidescent.md).
This note determines the exceptional quadratic and its first normal motion.
It does not construct a root section and proves no inverse or Keller mate.

**Subsequent continuation.**
[THM-3319](../01-canon/theorems/THM-3319-exceptional-quadratic-two-clutch-etale-persistence.md)
releases both external slopes `(d,k)` and proves that this fixed deck persists
as a connected finite-etale `C2` deck over a local algebraic etale germ.  It
does not produce a global component or rational section, and it does not yet
transport THM-3312's specific elimination-cofactor ratio across the germ.

THM-3312 universalizes the next native algebraic operation without choosing a branch.
For every `a+bt` in the quadratic field it records trace, norm, and conjugate-
difference square.  Applied to the critical `y`-root, first-normal velocity, and
elimination pair, it proves that the projective cofactor ratio itself generates
the quadratic extension: only its unordered trace/norm passport descends.
Although that elimination pair is unimodular, THM-3309 directly proves that
both true gradient components vanish on the critical fibre, so no Keller Bezout
row exists and the mate-integrability class is not entered.

The matching exact artifact is
[`jc_affine_c_exceptional_quadratic_blowup_scout_20260803.py`](../04-computation/jc_affine_c_exceptional_quadratic_blowup_scout_20260803.py),
with frozen
[`output`](../05-knowledge/results/jc_affine_c_exceptional_quadratic_blowup_scout_20260803.out).

## Inheritance pass and the typed question

The closest proved mechanism is
[THM-3306 -- affine-`c` square discriminant and transverse base locus](../01-canon/theorems/THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus.md).
In either accessory field `K_i`, it gives

```text
A_i=K_i[x]/(linear_a),             deg(linear_a)=36,
c_*=-b_0/b_1 in A_i,
linear_a=0, b_0+c_*b_1=0,
```

with a transverse coefficient ideal and a surviving quadratic
subresultant.  Its primary and independent companions prove that the common
`y`-gcd has degree exactly two.  They do not say whether those two roots are
defined over `A_i`.

The closest positive sidecar is the
[critical inverse-graph scout](jc-critical-inverse-cover-cofactor-jacobian-probe-agent-20260803.md):
away from this base locus, a unit linear subresultant selects the unique
critical root.  Its hostile example warns that a squarefree scalar resultant
does not force the quadratic pivot to remain a unit.

The least-used relevant mechanism is the first-normal-map discipline of
[THM-2985](../01-canon/theorems/THM-2985-multiparameter-normal-map-and-arc-factor-separation.md).
At a transverse base ideal, the next lawful object is its blow-up and strict
transform, not a formal division by one coefficient.

The live board was:

| Object | Native operation | Question / lost datum |
|---|---|---|
| degree-36 `linear_a` | blow up `(linear_a,b)` | which normal directions lie on the critical divisor? |
| surviving quadratic row | restrict to the exceptional divisor | split pair or quadratic deck exchange? |
| discriminant `delta` in `A_i` | norm and good-prime reduction | squareclass, not merely nonvanishing |
| first normal term | differentiate in a fixed blow-up chart | do the two directions move or remain stationary? |
| critical inverse graph | approach the base divisor | a root label is lost exactly at the quadratic extension |

The applicable method cards were **divide exceptional multiplicity before
judging a wall** and **audit and close sections under their next native
operation**.

## Two different pivots, kept different

There are two polynomials that earlier prose called `a`.  They are not the
same object.

1. `linear_a(x)` is the coefficient of `y` in the primitive linear
   subresultant after dividing the degree-44 boundary.  It has degree `36`, is
   irreducible, and defines `A_i`.
2. Here `P_2(x)` denotes the leading coefficient of the preceding quadratic
   row.  Before reduction modulo `linear_a`,

   ```text
   P_2=2V'+8V^2,
   ```

   and has degree `32`.

The companion uses the names `linear_a` and `P2` throughout and asserts the
two degrees separately.  No conclusion transfers merely because both are
PRS pivots.

## The exceptional-divisor equation

For

```text
L=y^2+y+CV,
R_1=2L(2y+1)+VA,
R_2=V^3+V^2y+L(-V'y+2V^2),
```

direct cubic cancellation gives the quadratic subresultant

```text
h=4R_2+V'R_1=P_2 y^2+Q_2 y+R_2^(0),                  (1)
```

where

```text
P_2=2V'+8V^2,
Q_2=2V'+12V^2,
R_2^(0)=V(4V^2+8CV^2+V'(2C+A)).                       (2)
```

To avoid overloading `R_2`, write the constant coefficient as `R_2^(0)` in
this note.  The script independently asks SymPy for the full generic
subresultant sequence and proves that its quadratic row is literally `(1)`,
not just a scalar multiple.

Let the next linear row be `(ell_1 y+ell_0)/4`.  A second direct derivation
proves the universal identity

```text
P_2 ell_0^2-Q_2 ell_0 ell_1+R_2^(0) ell_1^2
   =16P_2^2 Res_y(R_1,R_2).                            (3)
```

After owner/boundary localization, `(ell_1,ell_0)` differs from the fixed
primitive pair `(linear_a,b)` by one common unit.  Blow up in the chart

```text
u=linear_a,               b=u t.                       (4)
```

The strict transform of the linear row is `y+t=0`.  Dividing the double
base factor in `(3)` and setting `u=0` therefore gives the exact exceptional
equation

```text
F_0(t)=P_2 t^2-Q_2 t+R_2^(0)=0              in A_i[t]. (5)
```

Thus the surviving quadratic is not merely present at the base point: it is
the normal cone of the critical divisor.  Its roots are simultaneously the
two common `y` roots, with `y=-t`, and the two admissible blow-up directions.

The hostile charts are clean.  All three coefficients in `(5)` are units in
`A_i`.  Hence `t=0` is not a root.  In the other chart `s=1/t`, the equation
is

```text
P_2-Q_2s+R_2^(0)s^2=0,
```

and `s=0` is not a root because `P_2` is a unit.  There is no exceptional
direction at infinity.  Independently, the leading `y` coefficient of
`R_1` is the constant `4`, so the original cubics have no common projective
`y`-root at infinity.

## Exact nonsplitting and the degree-72 carrier

Put

```text
delta=Q_2^2-4P_2R_2^(0) in A_i.                        (6)
```

In both fields, `delta` is nonzero.  More strongly, its multiplication
characteristic polynomial over `K_i` is irreducible of degree `36`.  Thus
`delta` itself is a primitive generator of `A_i/K_i`.

The squareclass is decided by an exact norm obstruction.  If `delta=w^2` in
`A_i`, then `Norm_(A_i/K_i)(delta)=Norm(w)^2` in `K_i`.  Each norm is a unit
at the displayed simple degree-one prime of `K_i`, but reduces to a quadratic
nonresidue:

| accessory field | good reduction | norm residue | quotient factor degrees and characters |
|---|---|---:|---|
| `(4,1,1,1)` | `p=13`, `alpha=6` | `7`, non-square | `(1,-),(3,+),(15,-),(17,-)` |
| `(3,2,1,1)` | `p=23`, `alpha=16` | `22`, non-square | `(3,-),(13,+),(20,+)` |

The accessory polynomial has a simple root in each reduction, `linear_a`
remains squarefree of degree `36`, and `delta` is nonzero on every displayed
factor.  Euler exponentiation in each finite quotient gives the character
profile in the last column; its signs multiply to `-1`, independently
recovering the norm character.  These are exact characteristic-zero
nonsquare certificates, not statistical modular evidence.

Consequently, in both accessory fields,

```text
F_0 is irreducible and separable over A_i,
B_i=A_i[t]/(F_0)=A_i(sqrt(delta)) is a field of degree 2. (7)
```

The exceptional intersection is the connected closed point `Spec(B_i)` of
degree two over `Spec(A_i)`, hence degree `72` over `K_i`.  After fixing one
of the `36` geometric embeddings of `A_i` and base-changing that relative
fibre to an algebraic closure, it becomes two reduced directions exchanged by
the quadratic deck involution.  Over an algebraic closure of `K_i` there are
`72` geometric points in total, grouped as `36` such conjugate pairs.  There
is no `A_i`-rational exceptional direction and therefore no canonical
continuation of the off-divisor selected root over this base locus.

This is an algebraic `C_2` deck exchange on a fixed zero-dimensional slice.
Calling it global geometric or topological monodromy for a larger parameter
family would require a deformation and is not claimed.

## First normal motion in the fixed chart

The first normal term is well typed after fixing the same common
normalization used by THM-3306: `linear_a` is monic and `b` is divided by the
same leading unit.  Let

```text
b=b_0+cb_1,                 b_x=b_0'+c_*b_1'.
```

Because the base ideal is transverse, the path `(4)` has first derivatives

```text
dx/du=1/linear_a',
dc/du=(t-b_x/linear_a')/b_1.                            (8)
```

Here every `x` partial is taken with `c` fixed before substituting `c=c_*`.
Writing the strict transform as

```text
F(u,t)=F_0(t)+uF_1(t)+O(u^2),                           (9)
```

exact differentiation gives

```text
F_1(t)
 = (P_2,x/linear_a')t^2
 + (-Q_2,x/linear_a'+R_2,c^(0)/b_1)t
 + R_2,x^(0)/linear_a'
   -R_2,c^(0)b_x/(linear_a'b_1).                        (10)
```

All coefficients lie in `A_i`; `(10)` has degree two in `t`.  The exact
quotient computations give

```text
gcd(F_0,F_1)=1                                           (11)
```

in both fields.  Hence neither geometric exceptional branch is stationary
to first order.  If `t(u)=t_0+u dot(t)+...`, then

```text
dot(t)=-F_1(t_0)/F_0'(t_0).                             (12)
```

Since `(F_0')^2=delta` modulo `F_0`, the script reduces `(12)` to an explicit
linear form

```text
dot(t)=velocity_1 t+velocity_0 in B_i.                  (13)
```

Both coefficients are nonzero, and `velocity_1` is nonzero, so the two
conjugate velocities are distinct as well as nonstationary.

Because `F_0'(t_0)` is a unit in `B_i`, the formal implicit-function theorem
also gives a unique local power-series continuation `t(u)` through each
geometric exceptional direction after adjoining `sqrt(delta)`.  Equation
`(13)` is its first coefficient.  This is local continuation of the strict
critical divisor, not a global root label or a Keller inverse sheet.

The coefficient packet in `(10)` and the numerical scale in `(13)` depend
on the fixed normal coordinate `u`.  Multiplying the linear row by another
local unit reparametrizes `u`.  The invariant conclusions claimed here are
the exceptional quadratic, its nonsplit squareclass, absence of axis or
infinity roots, coprimality `(11)`, nonstationarity, and conjugate exchange;
the printed digests freeze the chosen normalization for audit rather than
promoting it to an intrinsic speed.

## Connection contract

```text
source:      transverse base ideal of the linear subresultant;
target:      degree-two closed point on the blow-up exceptional divisor;
map:         strict transform via (3), with t=b/linear_a and y=-t;
preserved:   the two common critical roots and their normal directions;
destroyed:   an individual root label over A_i and all Keller-mate data;
sidecar:     sqrt(delta), equivalently the quadratic field B_i;
test:        norm squareclass plus the first-normal gcd(F_0,F_1).
```

This sharpens the phrase “blow up the base ideal.”  The blow-up resolves the
projective coefficient map, but it does not split its critical root pair over
the original residue field.  The minimum lawful root-label sidecar is the
quadratic extension `B_i`, and its two labels move distinctly in the first
normal direction.

## Independent paths and reproduction

The companion does not import or execute either THM-3306 computation.  It
pins their hashes and reconstructs the two accessory responses from their
literal formulas.  Its two row paths are:

1. the full generic subresultant sequence; and
2. direct cancellation `4R_2+V'R_1`, followed by the explicit linear
   pseudo-remainder and identity `(3)`.

They agree exactly.  The nonsquare decision is checked both by the residue of
the characteristic-zero norm and by factorwise Euler characters in the
reduced quotient.  Normal and optimized runs agree byte for byte with the
stored transcript; the source has no Python assertions, floating literals,
or imported execution paths.

Run from the repository root:

```text
python3 04-computation/jc_affine_c_exceptional_quadratic_blowup_scout_20260803.py
python3 -O 04-computation/jc_affine_c_exceptional_quadratic_blowup_scout_20260803.py
```

The LF-normalized SHA-256 digests are

```text
script  6f050a583004172f812c3f7729427079d5df45c3a985c2e470b2a0d34ad8f337
output  0603517a0e97c1eb3d6b60051cfb02c2a8074ac3468028a36b97e7b8398f5670
```

## Scope and honest frontier

This result is sharply fixed-slice:

- only the two named THM-3212 accessory fields are covered;
- `C=c+x`, `d=k=1`, the THM-3289 localization, and its boundary quotient
  remain fixed;
- `delta` is the exceptional quadratic discriminant, not THM-3306's
  degree-36 parameter polynomial `D_i(c)`;
- the `C_2` exchange is algebraic on the fixed base locus, not a computed
  global braid or inertia group;
- separability gives local formal continuations after the quadratic base
  change, but the first-normal packet does not classify their higher
  coefficients or global gluing;
- the degree-119 residual repeated-`H` component remains unclassified;
- no primitive Keller cofactor, polynomial mate, inverse cover, `JC(2)`, or
  `DC(2)` follows.

The cheapest next deformation is to release one of `d` or `k`, retain the
quadratic squareclass over the resulting base divisor, and ask whether the
connected `C_2` cover ramifies, trivializes, or meets the owner walls.  A
root-label continuation should be attempted only after that cover, not by
choosing one root in a relative two-point geometric fibre by hand.
