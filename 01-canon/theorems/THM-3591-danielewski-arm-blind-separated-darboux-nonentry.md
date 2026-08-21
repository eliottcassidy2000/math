---
id: THM-3591
title: "Danielewski arm-jet interpolation and arm-blind separated Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On c^n e=Sigma(b), n>=2, every polynomial Darboux pair must interpolate
  the inverse derivative of Sigma in its first (c,e)-arm jet.  Hence
  no arm-blind pair exists; in particular, even after arbitrary
  Sigma-multiple corrections, two additively separated coordinates
  p(b)+u(c)+v(e) cannot form a Darboux pair.  The canonical calibrated
  two-piece attempt has bracket (Sigma J)'=1 mod Sigma and moves, rather
  than removes, the debt into the bulk ideal.  Explicit seven- and
  eight-piece candidates separate four distinct construction invoices:
  nonlinear central arms, inverse-derivative interpolation, the scalar
  bracket layer, and all remaining weight layers.  No Darboux pair and no
  counterexample to JC(2) is constructed.
source: root / planar-Jacobian construction-hostile session, 2026-08-21
audit: >
  An independent hostile audit rederived the intrinsic arm restriction,
  inverse-derivative sidecar, arbitrary Sigma-multiple and separated no-go,
  calibrated bulk debt, all seven/eight-piece weight invoices, and sharp
  degree, exponent, and characteristic boundaries.  Normal and optimized
  companions are byte-identical to the stored 455-gate output.
related:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3576-higher-exponent-belyi-keller-collision-tower
  - THM-3581-critical-value-multiarm-keller-compiler-and-A13-carrier
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
script: 04-computation/jc2_danielewski_arm_jet_interpolation_thm3591.py
output: 05-knowledge/results/jc2_danielewski_arm_jet_interpolation_thm3591.out
script_sha256: 95178037a86c52ecdf82234987acadfeb3ff3bdf478bbf205b40d541072fdf2f
output_sha256: a7ed1d7cdd54e7088c1d90462efb8082ce900dd069b7084fe39e90bf8dbfe6e1
hash_basis: raw LF bytes
---

# THM-3591 -- Danielewski arm-jet interpolation and arm-blind separated Darboux nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
gives a necessary arm-label sidecar for attempts to
construct a planar Jacobian counterexample through a smooth Danielewski
target.  It proves no general Darboux nonentry theorem.

Work over `C`.  Let `n>=2`, let `Sigma in C[b]` be squarefree of degree
`d>=2`, and put

```text
A=C[b,c,e]/(c^n e-Sigma(b)),                            (1)

{b,c}=c^n,          {c,e}=-Sigma'(b),
{b,e}=-n c^(n-1)e.                                      (2)
```

For every root `beta` of `Sigma`, the smooth affine arm is

```text
D_beta={b=beta,c=0}=Spec C[e].                           (3)
```

## 1. The exact arm-jet interpolation law

Suppose `P,Q in A` satisfy

```text
{P,Q}=kappa,                         kappa in C*.        (4)
```

Put

```text
Delta_PQ(b,e)
 =P_c(b,0,e)Q_e(b,0,e)-P_e(b,0,e)Q_c(b,0,e).            (5)
```

This is intrinsic: changing an ambient representative by a multiple of
`c^n e-Sigma` changes its `c`- and `e`-partials by terms vanishing on
`c=0` when `n>=2`.  Thus `(5)` is the determinant of the two differentials
in the tangent `(c,e)` frame along the reduced arm union.  On `(3)`, both
terms in `(2)` carrying `c^(n-1)` vanish because `n>=2`.  Therefore

```text
-Sigma'(beta) Delta_PQ(beta,e)=kappa                    (6)
```

for every root `beta` and every `e`.

Because `Sigma` is squarefree, there is a unique polynomial `J` of degree
less than `d` satisfying

```text
J(b)Sigma'(b)=1 mod Sigma(b).                            (7)
```

Equation `(6)` is equivalent to the global congruence

```text
Delta_PQ(b,e)=-kappa J(b) mod Sigma(b) in C[b,e].        (8)
```

Thus all positive `e`-coefficients of `Delta_PQ(b,e)` are divisible by
`Sigma`, while its constant coefficient must carry the nonconstant arm
label `J`.

The word **nonconstant** is automatic here.  If `J=lambda` were constant,
then every root would have `Sigma'(beta)=lambda^(-1)`.  But partial
fractions give

```text
sum_(Sigma(beta)=0) 1/Sigma'(beta)=0                    (9)
```

because the coefficient of `b^(-1)` in `1/Sigma(b)` is zero for `d>=2`.
If all derivatives were equal, the left side would be `d lambda`, nonzero
in characteristic zero.

For the normalized A13 target `Sigma=b(b^2+1)`, the sidecar is

```text
J=1+(3/2)b^2.                                            (10)
```

For the quadratic target `Sigma=b(b+4)`, it is

```text
J=(b+2)/8.                                               (11)
```

These are the cheapest exact witnesses that the different arms cannot be
treated by one unlabelled transverse jet.

## 2. Arm-blind and separated-coordinate nonentry

Call `(P,Q)` **arm-blind to first order** if the four functions

```text
P_c(beta,0,e), P_e(beta,0,e), Q_c(beta,0,e), Q_e(beta,0,e) (12)
```

are independent of the root `beta`.  Then `Delta_PQ(beta,e)` is one common
polynomial `H(e)`.  Equation `(6)` would force every `Sigma'(beta)` to be
the same nonzero constant, contradicting `(9)`.  Hence:

```text
no arm-blind polynomial pair has nonzero constant bracket. (13)
```

In particular, take arbitrary one-variable polynomials and arbitrary
corrections `R,S in A`:

```text
P=p(b)+u(c)+v(e)+Sigma(b)R,
Q=r(b)+s(c)+t(e)+Sigma(b)S.                             (14)
```

On every arm the `Sigma`-multiples disappear from the four transverse
derivatives, and

```text
Delta_PQ(beta,e)=u'(0)t'(e)-v'(e)s'(0),                 (15)
```

independent of `beta`.  Therefore no pair `(14)` satisfies `(4)`.  This is
an all-degree obstruction: raising the degrees of `u,v,s,t`, or appending
arbitrary terms visibly divisible by `Sigma`, cannot repair the missing
arm label.

## 3. The canonical calibrated near miss

The inverse derivative `(7)` also gives the cheapest positive repair of the
arm-blind failure.  Put

```text
P_0=c,                         Q_0=-J(b)e.               (16)
```

Its transverse determinant satisfies `(8)` on every arm.  Globally,

```text
{P_0,Q_0}=J Sigma'+Sigma J'=(Sigma J)'=1 mod Sigma.      (17)
```

For `d>=2`, `(Sigma J)'` is not constant: otherwise `Sigma J` would be
affine linear, impossible because `J` is nonzero and `deg Sigma>=2`.
Thus calibration does not solve the Darboux equation.  It moves the error
into the interior ideal:

```text
{P_0,Q_0}-1 is divisible by Sigma.                       (18)
```

This is the next constructive target.  A viable counterexample attempt must
retain the arm label `J`, preserve the two nonlinear conductor curves from
THM-3589, and cancel the bulk term `(18)` without collapsing back into an
arm-blind or stable-subalgebra ansatz.

### 3.1 A seven-piece candidate that pays every arm invoice

The arm calibration can be combined with the sharp nodal curves from
THM-3589.  On the normalized A13 target `Sigma=b(b^2+1)`, take `J` from
`(10)` and put

```text
P_*=c^3-Jc+e^2,
Q_*=c^2+e^3-e-(3/2)Jce.                                 (19)
```

The weight supports are

```text
supp(P_*)={-6,1,3},              supp(Q_*)={-9,-3,-2,2}. (20)
```

On the two central arms, `(19)` restricts to

```text
L_c: (c^3-c,c^2),                 L_e: (e^2,e^3-e).     (21)
```

Both curves are immersed and identify `+1` with `-1`; their tangent
directions at the arm crossing are independent.  More strongly, on every
simple arm `D_beta` one has

```text
Delta_(P_*,Q_*)=J,             so {P_*,Q_*}|D_beta=-1.  (22)
```

Thus `(19)` pays the nonlinear every-line debt, the transverse central-jet
debt, and the global inverse-derivative interpolation debt simultaneously.
It is nevertheless not a Darboux pair.  Its only scalar-weight channel is
the pair of weights `(1,-3)`, with coefficient functions `-J,-Sigma` on the
`c!=0` chart.  Its scalar coefficient is exactly

```text
W_(1,-3)(-J,-Sigma)
 =-J Sigma'-3J'Sigma=-1-(27/2)b Sigma(b).               (23)
```

No other layer can cancel a nonconstant scalar coefficient.  Thus `(23)` is
the scalar component of the bracket, while the scalar component of the
defect `{P_*,Q_*}+1` is `-(27/2)bSigma`.  The full quotient defect lies in
the arm-union ideal `(c,Sigma)`.  This is a first explicit seven-piece
near-counterexample.  Unlike an arm-blind ansatz, it fails only in the bulk.
Any attempted repair must cancel the scalar defect and the remaining
nonzero-weight layers without destroying `(21)--(22)`.

### 3.2 Paying the global scalar layer is not arm interpolation

There is a sharper seven-piece control which makes the scalar bracket
coefficient exactly constant.  Put

```text
q(b)=(225b^4+237b^2+8)/8,             v(b)=105b/16,

P_dagger=e^2+c+c^3,
Q_dagger=e^3+q(b)e+c^2+v(b)c^4.                        (24)
```

Its supports are again `3 x 4`:

```text
supp(P_dagger)={-6,1,3},
supp(Q_dagger)={-9,-3,2,4}.                             (25)
```

The central-arm restrictions are the sharp immersed nodal curves

```text
L_c: (c+c^3,c^2),                 L_e: (e^2,e+e^3).     (26)
```

On the `c!=0` chart the two scalar channels are `(1,-3)` and `(-6,4)`.
Their exact sum is

```text
W_(1,-3)(1,Sigma q)+W_(-6,4)(Sigma^2,v)=-1.            (27)
```

Thus `(24)` pays the nonlinear central-arm debt, the independent tangent
jet at their crossing, and the **global** scalar layer.  It does not pay
the labelled-arm interpolation bill.  Indeed

```text
Delta_(P_dagger,Q_dagger)=3e^2+q(b),
q(b)=J(b) mod Sigma,                                    (28)
```

so the bracket on `D_beta` is `-1-3Sigma'(beta)e^2`, not
constant.  It also fails immediately in a nonzero layer: the sum
`-6+2=-4` has a unique edge, whose coefficient is

```text
W_(-6,2)(Sigma^2,1)=4Sigma Sigma'!=0.                  (29)
```

No other edge can cancel it.  Thus scalar payment and central-arm geometry
do not imply the global labelled-arm condition.

### 3.3 One sidecar pays the arms and exposes the next layer

The missing arm interpolation in `(24)` has an exact one-piece repair:

```text
Q_ddagger=Q_dagger+(3/2)ce.                              (30)
```

This adds weight `-2`, so the supports have sizes `3 x 5`, eight pieces in
total.  On `c=0`, the new term contributes `Q_c=(3/2)e`, and therefore

```text
Delta_(P_dagger,Q_ddagger)
 =(3e^2+q)-2e(3e/2)=q=J mod Sigma.                      (31)
```

Consequently `{P_dagger,Q_ddagger}=-1` on every `D_beta`.  The new weight
has no scalar complement in `P_dagger`, so the exact scalar identity `(27)`
is unchanged.  The unique layer `(29)` is unchanged as well.  The repaired
eight-piece candidate therefore pays the two nodal central arms, their
independent tangent jet, the labelled inverse-derivative condition on every
arm, and the global scalar bracket coefficient, yet still fails at its first
uncancelled nonzero weight layer.  This is the sharpest explicit construction
control here: paying every boundary and scalar invoice is still not a Keller
pair.

## 4. Sharp boundaries and hostile examples

1. `deg Sigma>=2` is sharp.  For `Sigma=b`, the surface is `A2_(c,e)` and
   `(P,Q)=(c,-e)` has bracket one.
2. Characteristic zero is load-bearing.  In characteristic `p`,
   `Sigma=b^p-b` is squarefree with `Sigma'=-1` at every arm, and
   `(P,Q)=(c,e)` has bracket one.
3. The exponent `n>=2` is load-bearing for this proof.  At `n=1`, the term
   `-e(P_bQ_e-P_eQ_b)` survives on `c=0`; for example `(P,Q)=(b,e)` has
   arm bracket `-e`.  No exponent-one Darboux conclusion is asserted.
4. If `Sigma` has a multiple root, the surface and Poisson tensor are
   singular on that arm.  The smooth interpolation statement changes type.
5. Equation `(13)` excludes only arm-blind first jets.  It does not exclude
   pairs whose coefficients interpolate `J`, nor any general polynomial
   Darboux pair.

## 5. Preservation and loss ledger

```text
source       first transverse jets of (P,Q) on all arms D_beta
target       Delta_PQ in C[b,e]/(Sigma)
map          determinant of the (c,e)-jet
preserved    the constant symplectic bracket kappa
forced       -kappa*(Sigma')^(-1) modulo Sigma
lost         bulk terms divisible by Sigma and all higher normal jets
sidecar      J=(Sigma')^(-1) mod Sigma
cheap fail   one common arm-blind jet
positive     (c,-Je), exact on every arm but defective by Sigma in bulk
```

The inverse-derivative sidecar is necessary, not sufficient.  It is the
precise datum destroyed when a construction coalesces all arms before
solving the Darboux equation.

## 6. Exact verification contract

The companion uses exact rational arithmetic and explicit failure gates to
check:

- `(2)`, `(6)`, and the inverse-derivative congruence `(7)` for
  `2<=n<=7` and squarefree targets of degrees `2<=d<=7`;
- the quadratic and A13 sidecars `(10)--(11)`;
- separated and `Sigma`-corrected arm-blind controls;
- the calibrated near miss `(16)--(18)` and its nonconstant bulk debt;
- the seven-piece nodal-arm candidate `(19)--(23)`;
- the globally scalar-paid but nonscalar-layer-defective candidate
  `(24)--(29)` and its arm-interpolating eight-piece repair `(30)--(31)`;
- the degree-one, exponent-one, and positive-characteristic hostiles.

The proof is universal over `C`; the finite rows are controls, not an
extrapolation.

Reproduce with

```bash
python3 04-computation/jc2_danielewski_arm_jet_interpolation_thm3591.py
python3 -O 04-computation/jc2_danielewski_arm_jet_interpolation_thm3591.py
```
