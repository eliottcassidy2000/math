---
id: THM-3765
title: "Normalized three-consecutive-charge radial Keller classification"
status: >
  PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.  Over an
  algebraically closed characteristic-zero field, every smooth polynomial
  Q=X+chi(XT)+T psi(XT) has a rational constant-Jacobian mate exactly on one
  explicit degree-one generic-fibre boundary.  That boundary is
  Q=h+pT+X(1+gT/2)^2, with the sharp smoothness condition g=0 or p nonzero.
  A polynomial mate exists exactly when g=0, equivalently when Q is linear.
  All other smooth profiles have either a positive-genus holomorphic time
  form, a conic logarithmic residue, or the pure-profile residue obstruction.
  This is not proved canon until independent audit promotion and does not
  prove JC(2).
source: root + jc_quartic_c3_construct / 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The exact companion verifies the universal
  derivatives, torus resultant and absorbed-root criterion, axis boundary,
  generic quadratic model and Jacobian sign, pure-profile residue, complete
  degree-at-most-one coefficient classification, rational primitive, sharp
  smoothness boundary, and polynomial top-degree descent.  A direct
  coefficient-cube Groebner census and hostile Pell/absorbed profiles guard
  the formulas.  Normal and optimized runs byte-match the frozen transcript;
  independent hostile audit remains open.
depends_on:
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3754-affine-variable-euclidean-descent-classification
related:
  - THM-3741-radial-two-charge-keller-component-classification
  - THM-3757-pell-chebyshev-three-charge-hyperelliptic-obstruction-tower
  - THM-3759-central-product-profile-constant-flank-keller-classification
script: 04-computation/jc2_normalized_three_consecutive_charge_radial_classification_thm3765.py
output: 05-knowledge/results/jc2_normalized_three_consecutive_charge_radial_classification_thm3765.out
script_sha256: 76804229484cc0e2452081c766a69edb62775c87be89e4ff28729119de1a967b
output_sha256: c93f800c9f340252a76459d35f13afe5151c83693b27028afa306ddc8da03d62
semantic_sha256: a0400b6d12f55afc0019a2b30f1c4771f9d9b20cd805d713d495abfe7fa1c567
hash_basis: raw LF bytes
---

# THM-3765 -- the full normalized three-consecutive-charge radial ansatz

**PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.**  Give `X,T`
Euler weights `+1,-1`.  THM-3741 closes arbitrary profiles on the two
opposite charges.  Adding the intervening charge zero produces the first
genuine three-charge radial ansatz.  The new sector can absorb every affine
critical point, as THM-3757 demonstrates, but it cannot also make the generic
time form rationally exact except on one completely classified boundary.

Work over an algebraically closed field `k` of characteristic zero.  Let

```text
z=XT,                 chi,psi in k[z],
G=chi',               F=z psi,
Q(X,T)=X+chi(XT)+T psi(XT).                           (1)
```

Among the smooth polynomials `(1)`:

1. `Q` has a rational mate `P in k(X,T)` with `J(P,Q) in k*` if and only if

   ```text
   chi=h+gz,                 psi=p+(g^2/4)z,           (2)
   ```

   where `h,g,p in k` and the smoothness boundary is

   ```text
   g=0, or p!=0.                                      (3)
   ```

2. `Q` has a polynomial constant-Jacobian mate if and only if `g=0` in
   `(2)`.  These are exactly the nonzero linear forms

   ```text
   Q=X+h+pT.                                          (4)
   ```

3. If `g!=0,p!=0`, the nonlinear rational-exact boundary is

   ```text
   Q=h+pT+X(1+gT/2)^2.                                (5)
   ```

   It is smooth and has rational mates, but no polynomial mate.

Thus every smooth member of `(1)` has a polynomial Jacobian mate exactly
when it is linear.  This closes the normalized three-consecutive-charge
radial component design, not arbitrary three-charge supports or `JC(2)`.

## 1. Exact smoothness criterion

The universal derivatives are

```text
Q_X=1+TG+T^2 psi',                 Q_T=XG+F'.          (6)
```

The axis `T=0` is always safe because `Q_X=1`.  On `X=0`, simultaneous
vanishing is possible only when `psi(0)=0`; in that case the remaining
derivative is

```text
1+T G(0)+T^2 psi'(0).                                (7)
```

A nonconstant polynomial has a root over `k`.  Hence the `X=0` axis is safe
exactly when

```text
psi(0)!=0,
or psi(0)=psi'(0)=G(0)=0.                             (8)
```

On the torus eliminate `T=z/X`.  The critical equations become

```text
X^2+zGX+z^2psi'=0,                 GX+F'=0,            (9)
```

and their exact resultant is

```text
R(z)=F'(z)^2-F(z)G(z)^2.                              (10)
```

If `psi=0`, the torus is automatically safe: `Q_T=XG=0` forces `G=0`,
after which `Q_X=1`.  In this pure branch, `(8)` reduces to the exact
smoothness condition `G(0)=0`.

Suppose `psi!=0`.  Then smoothness forces `R!=0`.  At a nonzero root of
`R`, equations `(9)` have no common nonzero `X` exactly when

```text
psi(z)=psi'(z)=0.                                    (11)
```

Indeed, if `G!=0`, the unique common root is `X=-F'/G`; it is zero exactly
when `F'=0`, and `(10)` then forces `F=zpsi=0`.  Conversely
`psi=psi'=0` makes that common root zero.  If `G=0`, then `F'=0` and the
first equation is `X^2+z^2psi'=0`, which has no nonzero root exactly in the
same case.  Combining `(8),(10),(11)` gives a complete rootwise smoothness
test:

```text
psi!=0, Q smooth
iff (8), R!=0, and every nonzero root of R is a multiple root of psi. (12)
```

This is the precise absorbed-root sidecar missing from the scalar
resultant alone.

## 2. The generic fibre and automatic squarefreeness

Let `Lambda` be transcendental over `k` and work on `Q=Lambda`.  Define

```text
Y=X-Tpsi(z)=X-F(z)/X.                                 (13)
```

Since `Lambda-chi=X+F/X`, one has

```text
Y^2=Delta(z),              Delta=(Lambda-chi(z))^2-4F(z), (14)
X=[Lambda-chi+Y]/2,        T=z/X.                     (15)
```

For `psi!=0`, `(14),(15)` identify the generic function field; this is not
merely a quotient.  If `Q` is smooth, `Delta` is automatically squarefree
in `k(Lambda)[z]`.  For a common root `alpha` of `Delta,Delta'`, put
`W=Lambda-chi(alpha)`.  Then

```text
W^2=4F(alpha),                 WG(alpha)+2F'(alpha)=0. (16)
```

Squaring the second equation gives `R(alpha)=0`.  If `F(alpha)!=0`, the
torus criterion produces a critical point, so smoothness forces `F(alpha)=0`
and then `W=0`.  The nonzero polynomial `R` is fixed over `k`, hence
`alpha` is algebraic over `k`; algebraic closure gives `alpha in k`, and
`Lambda=chi(alpha)` is impossible.  Thus `Delta` is squarefree.

The same coordinates give the exact sign

```text
J(Q,z)=XQ_X-TQ_T=Y.                                  (17)
```

Consequently a rational mate `J(P,Q)=c in k*` would satisfy on the generic
fibre

```text
partial_z P=-c/Y,                    dP=-c dz/Y.       (18)
```

## 3. Every generic degree at least two is rationally obstructed

Put `d=deg_z Delta`.  If `d>=3`, then `dz/Y` is a nonzero holomorphic
differential on the smooth projective hyperelliptic model of `(14)`: finite
branch points are regular in the parameter `Y`, and the standard infinity
calculation has no pole.  A derivative of a rational function cannot be a
nonzero holomorphic differential.

If `d=2`, the smooth conic has two points at infinity over an algebraic
closure of `k(Lambda)`.  Writing `a` for the nonzero leading coefficient of
`Delta`, the residues of `dz/Y` there are

```text
+/-1/sqrt(a),                                          (19)
```

so the differential is again not rationally exact.  Therefore

```text
psi!=0, Q smooth, rational mate  =>  deg_z Delta<=1.  (20)
```

The pure branch `psi=0` has the rational generic coordinate
`X=Lambda-chi(z)`.  Equation `(17)` becomes

```text
dP=-c dz/[Lambda-chi(z)].                              (21)
```

If `chi` is nonconstant, the denominator is generically separable and every
root has a nonzero logarithmic residue.  Thus this branch has no rational
mate unless `chi` is constant, reproducing the diagonal residue obstruction
of THM-3551.

## 4. The complete rational-exact boundary

The coefficient of `Lambda` in `Delta` prevents cancellation across radial
degrees.  Condition `deg_z Delta<=1` first forces

```text
chi=h+gz.                                              (22)
```

Then `F(0)=0`, and cancellation of the quadratic coefficient gives exactly

```text
F=pz+(g^2/4)z^2,                  psi=p+(g^2/4)z.      (23)
```

Conversely `(22),(23)` make `Delta` linear or constant.  In source
coordinates they give `(5)`.  Its derivatives show the sharp boundary:

```text
Q_X=(1+gT/2)^2,
Q_T=p+gX(1+gT/2).                                    (24)
```

If `g=0`, `Q` is linear and smooth for every `p`.  If `g!=0`, a critical
point can only lie on `T=-2/g`, where `Q_T=p`; hence smoothness is equivalent
to `p!=0`.  This proves `(3)`.

For `g!=0,p!=0`, set

```text
D=g(Q-h)+2p,              Y=X-T[p+(g^2/4)XT].         (25)
```

Holding `Q` fixed, `Delta'=-2D` and `2YY_z=-2D`.  Thus

```text
P_0=cY/D,
J(P_0,Q)=c.                                            (26)
```

All rational mates differ from `(26)` by an element of `k(Q)`.  If `g=0`
and `p!=0`, one may use `P_0=c(X-pT)/(2p)`; if `g=p=0`, use `P_0=-cT`.

## 5. Polynomial exactness occurs only on the linear edge

The nonlinear rational-exact boundary `(5)` is affine in `X`.  Write

```text
A(T)=(1+gT/2)^2,                    f(T)=h+pT,
Q=f(T)+XA(T).                                          (27)
```

For a hypothetical polynomial mate

```text
P=sum_(j=0)^N p_j(T)X^j,              p_N!=0,          (28)
```

the top `X^N` coefficient of `J(P,Q)=c` is

```text
N A'p_N-Ap_N'=0,
(p_N/A^N)'=0.                                         (29)
```

Hence `p_N=lambda A^N`.  Subtracting `lambda Q^N` preserves the Jacobian
and strictly lowers the `X` degree.  Finite descent leaves `P=p(T)`, where

```text
J(p(T),Q)=-A(T)p'(T).                                 (30)
```

For `g!=0`, the nonconstant `A` cannot divide a nonzero scalar.  Thus no
nonlinear member of `(5)` has a polynomial mate.  This is the typed dual of
THM-3754's Euclidean descent.  When `g=0`, `(4)` is linear and the displayed
positive mates prove existence.  Equations `(20)--(30)` establish both
classifications.  **QED conditional on audit promotion.**

## 6. Exact controls and residual frontier

Reproduce with

```bash
python3 -B 04-computation/jc2_normalized_three_consecutive_charge_radial_classification_thm3765.py
python3 -B -O 04-computation/jc2_normalized_three_consecutive_charge_radial_classification_thm3765.py
```

The assertion-free companion checks the universal resultant, generic-fibre
equation and Jacobian sign, the rational-exact primitive and polynomial
descent, a complete low-degree coefficient cube against direct Groebner
smoothness, Pell-Chebyshev absorbed-root controls, and the first `R=cF^2`
absorbed profile.  These computations are hostile controls; the preceding
rootwise and differential proofs are all-degree.

THM-3757 is now a structured hostile subfamily rather than the source of the
classification: it proves that the high-degree branch is populated, smooth,
and genuinely three-charge.  THM-3759 treats the narrower constant-flank
slice.  The next live radial design must vary the `+1` flank as well, while a
nonradial design must leave the single product coordinate `z=XT` entirely.
