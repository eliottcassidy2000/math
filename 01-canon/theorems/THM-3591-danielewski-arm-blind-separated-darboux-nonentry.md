---
id: THM-3591
title: "Danielewski arm-jet interpolation and arm-blind separated Darboux nonentry"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / PENDING INDEPENDENT AUDIT.
  On c^n e=Sigma(b), n>=2, every polynomial Darboux pair must interpolate
  the inverse derivative of Sigma in its first transverse arm jet.  Hence
  no arm-blind pair exists; in particular, even after arbitrary
  Sigma-multiple corrections, two additively separated coordinates
  p(b)+u(c)+v(e) cannot form a Darboux pair.  The canonical calibrated
  two-piece attempt has bracket (Sigma J)'=1 mod Sigma and moves, rather
  than removes, the debt into the bulk ideal.  No Darboux pair and no
  counterexample to JC(2) is constructed.
source: root / planar-Jacobian construction-hostile session, 2026-08-21
audit: >
  The universal proof and exact controls are present.  An independent
  hostile audit and immutable hashes remain pending.
related:
  - THM-3561-squarefree-danielewski-completion-and-poisson-descent-gate
  - THM-3576-higher-exponent-belyi-keller-collision-tower
  - THM-3581-critical-value-multiarm-keller-compiler-and-A13-carrier
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
script: 04-computation/jc2_danielewski_arm_jet_interpolation_thm3591.py
output: 05-knowledge/results/jc2_danielewski_arm_jet_interpolation_thm3591.out
script_sha256: PENDING
output_sha256: PENDING
hash_basis: raw LF bytes
---

# THM-3591 -- Danielewski arm-jet interpolation and arm-blind separated Darboux nonentry

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / PENDING INDEPENDENT
AUDIT.**  This theorem gives a necessary arm-label sidecar for attempts to
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

This is intrinsic: it is the determinant of the two differentials in the
normal `(c,e)` directions along the reduced arm union.  On `(3)`, both
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
- the degree-one, exponent-one, and positive-characteristic hostiles.

The proof is universal over `C`; the finite rows are controls, not an
extrapolation.

Reproduce with

```bash
python3 04-computation/jc2_danielewski_arm_jet_interpolation_thm3591.py
python3 -O 04-computation/jc2_danielewski_arm_jet_interpolation_thm3591.py
```
