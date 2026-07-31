---
id: THM-2759
title: "Split degree-six componentwise monicization closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every planar
  Keller pair in the complete chosen-sheet split polynomial exact-square-
  prefix reduced-degree-six terminal family is an automorphism.  The two
  weighted boundary points are excluded or slope-four monicized
  componentwise; q/h^3 then extends to zero on the complete source, forcing
  the vertical branch, whose full coefficient stratification is empty.  The
  argument includes constant U and reducible/nonreduced ambient response
  schemes.  It does not derive the chart for an arbitrary Keller pair, treat
  split degrees 10,14,18 or at least 26, or prove JC(2) or DC(2).
source: root/split-degree-six-monicization-2026-07-28
audit: >
  degree6-hostile-audit-2026-07-28 independently rebuilt the complete reduced
  bank by Laurent and multinomial expansion, audited every infinity valuation
  case including cancellation/equality/zero coordinates, checked weighted-
  index descent, polynomiality, vertical strata, constant U, and
  reducible/nonreduced boundaries, and returned SOUND.  Normal and optimized
  exact replays byte-match the stored output; the script has no assert nodes.
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2181-exact-square-prefix-compression-and-monic-depressed-quartic-closure
  - THM-2202-uniform-all-degree-quartic-pole-closure
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
related:
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2725-split-even-faber-nonzero-first-flux-unified-closure
  - THM-2745-highest-odd-faber-componentwise-exact-prefix-closure
  - THM-2752-all-even-zero-first-flux-response-regularization-closure
  - THM-2755-all-even-zero-flux-componentwise-global-regularity-closure
script: 04-computation/jc2_split_degree6_componentwise_monicization_thm2759.py
output: 05-knowledge/results/jc2_split_degree6_componentwise_monicization_thm2759.out
script_sha256: 7580c022febc86b79d0adda10b3976fb176a45fc6f8a1d76b8e9783c8542bc9d
output_sha256: 66d91a4448bac3490a22e240288519b349b4ec425f3a483aa80bc6447a619a09
hash_basis: LF-normalized bytes
---

# THM-2759 -- split degree six closes by componentwise slope-four monicization

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This closes exactly the chosen-sheet split polynomial exact-square-prefix
reduced-degree-`6` terminal family.  It does not derive that chart for an
arbitrary Keller pair, does not treat reduced degrees `10,14,18` or `>=26`,
and does not prove `JC(2)` or `DC(2)`.

## 1. Exact family and inherited dichotomy

Let `(P,Q)` be a polynomial Keller pair over `C`.  Enter the quartic,
twice-odd terminal chart and assume its canonical polynomial exact square
prefix is split.  THM-2202 gives

```text
P=H^2+L,
H=U^2 z^2+Bz+C,                 L=Az+E,              (1)
```

with `U,B,C,A,E in C[x]`, `U!=0`.  Choose the sheet `U` and put

```text
w=Uz+B/(2U),
d=C-B^2/(4U^2),
q=A/U,
s=AB/(2U^2)-E.                                      (2)
```

Then

```text
P=w^4+2d w^2+q w+(d^2-s)=(w^2+d)^2+(qw-s).          (3)
```

After the exact target-shear quotient of THM-2230, normalize the top reduced
coefficient.  The complete chosen-sheet degree-six Faber bank is

```text
Q=E_6+a_5E_5+a_3E_3+a_2E_2+a_1E_1,                 (4)
```

where all `a_j` are constants.  No parity column has been discarded.

THM-2723 gives constants `lambda,W,kappa`, with `kappa!=0`, such that

```text
Phi_Q=lambda,                  Psi_Q=W,
U R_Q'=kappa,                                            (5)
```

and exactly one of

```text
U in C*,
U=u_0(x-a)^m,           u_0 in C*, m>=2.             (6)
```

holds.  If `U` is constant, `(x,z)->(x,w)` is a polynomial source
automorphism and `(3)` is a monic depressed quartic with polynomial
coefficients.  THM-2181 makes `(P,Q)` an automorphism.  It remains to exclude
the second case of `(6)`.

## 2. The complete degree-six response bank

For `j=1,2,3,5,6`, the three THM-2129 observables are

```text
j   Phi_j                              Psi_j
1   2d                                 q
2   2q                                -2s
3   (3/2)(d^2-2s)                    -(3/2)dq
5   (5/8)(2d^3-4ds+q^2)             -(5/8)q(d^2+2s)
6  -3qs                              -(3/2)(dq^2-s^2)

j   R_j
1   (d^2-2s)/2
2  -dq
3   (4d^3-3q^2)/8
5   (5/32)(3d^4-4d^2s-4dq^2+4s^2)
6  -q(q^2-6ds)/4.                                      (7)
```

Adjoin `h` of weight one to `wt(d,q,s)=(2,3,4)`.  Equations `(5)` become

```text
F_7=Phi_6+a_5h Phi_5+a_3h^3 Phi_3+a_2h^4 Phi_2
             +a_1h^5 Phi_1-lambda h^7=0,

G_8=Psi_6+a_5h Psi_5+a_3h^3 Psi_3+a_2h^4 Psi_2
             +a_1h^5 Psi_1-W h^8=0,                  (8)

R_aff=R_9/h^9,
R_9=R_6+a_5h R_5+a_3h^3 R_3+a_2h^4 R_2+a_1h^5R_1.
```

The top forms are coprime.  Their common support at `h=0` is exactly

```text
P_d=[0:1:0:0],                 P_q=[0:0:1:0].        (9)
```

Indeed `qs=0` and `dq^2=s^2`.  At `P_q`,

```text
R_6(0,1,0)=-1/4.                                     (10)
```

Thus `R_aff` has a pole on every branch over `P_q`.

The coefficient map is nonconstant by `(5)`.  Its projective closure is an
integral curve inside `(8)`: the top gcd excludes a common hypersurface, and
the closure of a nonconstant integral image is a reduced irreducible curve
component.  After normalization, the rational source map extends to a finite
surjective map from `P^1_x`.

## 3. The `P_q` boundary contradicts polynomiality of `Q(x,0)`

Assume now that `U=u_0(x-a)^m`, `m>=2`.  THM-2723 gives

```text
R_Q=r_0+r_1(x-a)^(1-m),             r_1!=0.           (11)
```

Hence `R_Q` has exactly one pole, at the finite point `x=a`.  By `(10)`,
every source point above `P_q` would have to be `a`.  In particular all five
polynomial coefficients in `(1)` are regular there.

Work on a finite local `q=1` index cover over such a point and write
`v(h)=alpha>0`.  Put

```text
beta=h B/(2U).                                        (12)
```

The two exact prefix identities are

```text
beta^2+d=h^2C,                   q beta-s=h^4E.       (13)
```

At `P_q`, both `d` and `s` vanish and `q` is a unit.  The first identity and
regularity of `C` force `v(beta)>0`.  Therefore, in affine coordinates,

```text
v(q_aff)=v(q/h^3)=-3alpha,
v(beta_aff)=v(beta/h)>-alpha.                         (14)
```

Now evaluate the Faber polynomials at the original polynomial section
`z=0`.  Since `H(x,0)=C` and `L(x,0)=E`, exact expansion gives

```text
E_1(0)=beta_aff,
E_2(0)=C,
E_3(0)=(6C beta_aff-2beta_aff^3+3q_aff)/4,

E_5(0)=(15C^2 beta_aff-10C beta_aff^3+5Cq_aff
         +10E beta_aff+3beta_aff^5-5beta_aff^2q_aff)/8,

E_6(0)=(8C^3+12CE+3q_aff^2)/8.                       (15)
```

The last summand in `E_6(0)` has exact valuation `-6alpha`.  Every term in
`E_5(0)` has valuation strictly greater than `-5alpha`; every term in
`E_3(0)` has valuation at least `-3alpha` or strictly greater; `E_1(0)` is
strictly shallower than `-alpha`, and `E_2(0)` is regular.  Thus the
normalized top term contributes a unique pole to `Q(x,0)`.  This contradicts
`Q in C[x,z]`.  Consequently the physical image contains no `P_q` boundary.

## 4. A response-pole lemma at `P_d`

Consider a normalization branch over `P_d` on a finite `d=1` index cover.
Write

```text
v(h)=alpha>0,             x_q=v(q)>0,
x_s=v(s)>0,               b=min(x_q,x_s),             (16)
```

allowing either or both of `x_q,x_s` to equal `infinity`.  The following
implication uses only `(7)--(8)`:

```text
b<4alpha       implies       v(R_9)<9alpha.           (17)
```

In other words, every sub-slope-four branch makes `R_aff` polar.

Here is the complete valuation proof.  If `x_q<x_s`, the top of `G_8` has
order `2x_q`.  Its only possible competitors below slope four are the active
odd rows of orders

```text
alpha+x_q, 3alpha+x_q, 5alpha+x_q.                   (18)
```

Let `g in {1,3,5}` be the least active odd gap.  If `g alpha` is smaller or
larger than `x_q`, respectively the odd or top `G_8` term is unique.  At
equality, the corresponding `F_7` row has order `g alpha=x_q` and is unique,
again impossible.  The `a_2` and `W` rows are later because `x_q<4alpha`.
If no odd row is active, the top `G_8` term is itself unique.  Thus
`x_q<x_s` cannot occur.  The case `x_s<x_q` is analogous: if an odd row can
meet or precede the top `s^2` term, its gap satisfies
`g alpha< x_s`, so its unit `Phi` term is uniquely earliest in `F_7`.
With no active odd row, the top `s^2` term is unique.

Hence `x_q=x_s=b`.  If the leading coefficient of `q^2-s^2` is nonzero,
the same comparison between `2b` and `g alpha+b` gives a unique `G_8` or
unit `F_7` term; if no odd row is active, the top `G_8` term is unique.
Therefore

```text
s_0=epsilon q_0,                         epsilon in {1,-1}. (19)
```

The top of `F_7` is now the nonzero term `-3epsilon q_0^2` of order `2b`.
It must meet the least active member of

```text
(a_5,g=1), (a_3,g=3), (a_1,g=5), (lambda,g=7),
```

so `g alpha=2b`.  In the `lambda` case the leading response is simply
`(3/2)epsilon q_0^2`.  In the three odd cases, the constant
`(Phi_j,R_j)` pairs at `(d,q,s)=(1,0,0)` are

```text
j=5: (5/4,15/32),      j=3: (3/2,1/2),
j=1: (2,1/2).                                          (20)
```

Using the `F_7` cancellation, the residual response multipliers are

```text
3/2+3(15/32)/(5/4)=21/8,
3/2+3(1/2)/(3/2)=5/2,
3/2+3(1/2)/2=9/4.                                   (21)
```

All are nonzero.  Thus `v(R_9)=2b=g alpha<9alpha`, proving `(17)`.

## 5. Exact-prefix exclusion at the unique response pole

Let a source point over `P_d` satisfy `b<4alpha`.  By `(17)`, `R_Q` has a
pole there.  Equation `(11)` says that this point is the finite point `x=a`,
where `C,E` are regular.  Use `(12)--(13)` again.  Since `d` is a unit,
`beta` is a unit and its residue satisfies

```text
beta_0^2=-d_0.                                      (22)
```

The second identity in `(13)` and `b<4alpha` force

```text
v(q)=v(s)=b,                    s_0=beta_0 q_0.       (23)
```

Both top fluxes then have exact order `2b` and nonzero residues:

```text
Phi_6: -3 beta_0 q_0^2,
Psi_6: -(3/2)(d_0-beta_0^2)q_0^2=-3d_0q_0^2.         (24)
```

In `F_7`, the even lower row has order `4alpha+b>2b`; the three odd unit
rows and the scalar flux have distinct gaps `1,3,5,7`.  Therefore
regularity can hold only if the least active gap satisfies

```text
g alpha=2b.                                          (25)
```

But every corresponding row in `G_8` has order `g alpha+b>2b`; the even
row has order `4alpha+b>2b`, and `Wh^8` is later.  The nonzero top `Psi_6`
term in `(24)` is unique, a contradiction.

Combining this with `(17)` proves, at every source point over `P_d`,

```text
min(v(q),v(s))>=4v(h).                               (26)
```

Consequently

```text
v(q_aff)=v(q/h^3)>=v(h)>0.                           (27)
```

The coefficient strata in this argument are exhaustive:

```text
least constant Phi gap 1:  a_5 !=0;
least constant Phi gap 3:  a_5=0, a_3!=0;
least constant Phi gap 5:  a_5=a_3=0, a_1!=0;
least constant Phi gap 7:  a_5=a_3=a_1=0, lambda!=0;
no constant Phi gap:       a_5=a_3=a_1=lambda=0.      (27a)
```

Below slope four, the `a_2` row has order `4alpha+b` in both fluxes and
can never be first; `W h^8` is later as well.  The no-gap stratum leaves the
nonzero top `Phi_6` term unique at the finite response pole.  Thus no
coefficient specialization was omitted.

## 6. Global regularity forces the vertical branch

At a source point mapping to `h!=0`, `q_aff` is a regular affine coordinate.
Section 3 excludes `P_q`; equation `(27)` makes `q_aff` regular and zero at
every remaining boundary point.  A boundary point exists, since a
nonconstant complete curve cannot lie in the affine chart `h!=0`.

Thus `q_aff` is a global regular function on `P^1_x` and vanishes somewhere.
It follows that

```text
q_aff=0.                                             (28)
```

Equivalently `A=0`, and from `(2)`

```text
s=-E in C[x].                                        (29)
```

## 7. The vertical branch is empty

At `q=0`, the second flux in `(5)` is

```text
(3/2)s^2-2a_2s=W.                                    (30)
```

Hence `s` is algebraic over `C`, so `s in C`.  The first flux becomes

```text
(5a_5/4)d^3+(3a_3/2)d^2
 +(2a_1-(5a_5/2)s)d-3a_3s=lambda.                    (31)
```

Unless this polynomial in `d` is identically zero, `(31)` makes
`d in C`.  Then `R_Q` is constant, contrary to `(5)`.  Identity in `d`
occurs exactly when

```text
a_5=a_3=a_1=lambda=0.                                (32)
```

In that all-even exceptional case, every surviving response row at `q=0`
vanishes, so `R_Q=0`, again contradicting `(5)`.

Equivalently, the vertical coefficient strata are:

```text
a_5 !=0:                         cubic equation makes d constant;
a_5=0, a_3!=0:                  quadratic equation makes d constant;
a_5=a_3=0, a_1!=0:              linear equation makes d constant;
a_5=a_3=a_1=0, lambda!=0:       the first-flux equation is inconsistent;
a_5=a_3=a_1=lambda=0:           R_Q=0 identically.    (32a)
```

The parameters `a_2` and `W` are arbitrary throughout this split: `(30)`
still makes `s` constant, and `a_2` contributes no vertical response.

Therefore the nonconstant alternative in `(6)` is impossible.  The constant
alternative is an automorphism by THM-2181.  This closes the complete
chosen-sheet split polynomial exact-square-prefix reduced-degree-six
terminal family.

## 8. Weighted-index and normalization audit

The charts `q=1` at `P_q` and `d=1` at `P_d` may require finite cyclic index
covers of the weighted projective chart.  Every valuation in Sections 3--5
is taken after such a cover.  Finite ramification multiplies

```text
v(h), v(q), v(s), v(R_9), v(q/h^3)
```

by the same positive integer.  Hence the comparisons `<4v(h)`,
`v(R_9)<9v(h)`, and `v(q/h^3)>0` are invariant and descend to the coarse
normalization.  The functions `q_aff` and `R_aff` themselves are invariant
rational functions, not index-cover coordinates.

The use of source regularity is localized.  It is **not** assumed at an
arbitrary point over infinity.  First `(10)` or `(17)` proves that `R_Q` is
polar; only then does `(11)` identify the source point with the finite zero
of `U`, where `B,C,A,E` are regular.  At a source point where `R_Q` is
regular, the contrapositive of `(17)` already supplies slope four.

The normalized image need not be the whole complete intersection and the
complete intersection need not be irreducible.  The argument uses only the
reduced integral closure of the actual nonconstant image.  Coprime top forms
exclude a surface component; properness extends the source map; and finite
surjectivity ensures that every boundary point on that image has a source
preimage.  Nilpotents or unrelated components therefore do not enter.

## 9. Exact companion and false shortcuts

The exact companion

```text
04-computation/jc2_split_degree6_componentwise_monicization_thm2759.py
```

independently reconstructs `(7)` by Laurent recurrence and multinomial
extraction, checks `(9)--(10)`, all formulas `(15)`, the three nonzero
response multipliers `(21)`, and the vertical specialization `(30)--(32)`.

The proof uses two coordinates absent from the failed primitive `C_2` test:

1. polynomiality at the original section `z=0`, which kills `P_q`; and
2. the actual rational primitive, which localizes every sub-slope-four
   response pole at the unique finite zero of `U`, where the exact-prefix
   coefficient identities are regular.

The mechanism is special to degree six.  At higher degree the top boundary
and response may have more points and more cancellation strata; no uniform
all-degree monicization is claimed.

Several tempting shorter arguments are false or insufficient:

1. **The primitive `C_2` bank does not raise order.**  The exact-square vector
   `(128,32,1)` survives the first augmentation grade, and the fourth
   Hamiltonian row detects only a free additive constant.  The present proof
   does not reuse that refuted implication.
2. **One response pole does not force `U` constant.**  For example
   `U=(x-a)^2`, `R=-1/(x-a)` has `UR'=1`.  Here a response pole first locates
   the unique finite root; polynomiality and the exact prefix do the work.
3. **Pole capacity alone cannot close degree six.**  The boundary has only
   two points, within THM-2723's abstract `q`-pole capacity.  The decisive
   fact at `P_q` is the uncancellable `q_aff^2` term in the polynomial value
   `Q(x,0)`.
4. **`q=0` does not kill odd response rows by parity.**  `R_1,R_3,R_5` remain
   nonzero there.  The separate `Psi`-then-`Phi` constant-field argument in
   Section 7 is necessary.
5. **Constant `U` is not declared empty.**  It makes the chosen sheet a
   polynomial coordinate, after which THM-2181 proves automorphy.  The result
   excludes the nonconstant monicization obstruction, not polynomial
   automorphisms.
