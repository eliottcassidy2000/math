---
id: THM-2865
title: "Gamma transverse null holotopy and uniform fourth exit"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The THM-2846 positive-cone
  factorial moment-three hostile extends through every Gamma functional
  L_alpha(s^m)=(alpha)_m with 9/10<=alpha<=11/10.  An explicit rational
  moving tube contains exactly one transverse null pair at every
  parameter, and its fourth moment has a uniform nonzero imaginary exit.
  Only alpha=1 is the one-complex-dimensional Gaussian radial law.
source: root/audit-2809-gamma-transverse-holotopy-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
related:
  - THM-2841-all-order-adjacent-difference-factorial-tensor-positivity
  - THM-2853-gamma-adjacent-tensor-cycle-weighted-positivity
  - THM-2855-shifted-positive-cone-transverse-asymptotic-family
script: 04-computation/gmc_gamma_transverse_null_holotopy_thm2865.py
output: 05-knowledge/results/gmc_gamma_transverse_null_holotopy_thm2865.out
script_sha256: ff32c787fb4fef58b77d06d84f3b032197744515ec03d6dcff8c1be6e70afc52
output_sha256: 6b943778460db66d48b50fe6d62f94feb6df62bbe08006df287ffcd4ca9ec0ec
hash_basis: LF-normalized bytes
---

# THM-2865 -- Gamma transverse null holotopy and uniform fourth exit

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The isolated factorial hostile of THM-2846 is not an accident of the
single radial law `L(s^m)=m!`.  It is one fibre of an exact one-parameter
family of transverse positive-cone nulls.  Here “holotopy” means this
certified continuation of a null section together with its fixed exit
direction; it does not assert a topological equivalence.

## 1. Gamma-adjacent tensors

For real `alpha>0`, define

```text
L_alpha(s^m)=(alpha)_m,
f_n=s^n/(alpha)_n,
d_n=f_(n+1)-f_n.                                      (1)
```

The Gamma integral

```text
L_alpha(F)
 =Gamma(alpha)^(-1) int_0^infinity
    F(s)s^(alpha-1)e^(-s) ds                          (2)
```

shows that `L_alpha(F^2)>0` for every nonzero real polynomial `F`.
Also `L_alpha(d_n)=0`.

For indices `n_1,...,n_k`, direct expansion gives the exact tensor

```text
L_alpha(prod_i d_(n_i))
 =
 [sum_(epsilon in {0,1}^k)
    (-1)^(k-|epsilon|)
    (alpha)_(sum_i(n_i+epsilon_i))
    prod_(epsilon_i=0)(alpha+n_i)]
 /
 [prod_i (alpha)_(n_i+1)].                            (3)
```

The companion derives `(3)` rather than importing a closed tensor table.
At `alpha=1` it independently checks every order-two, -three, and -four
entry with indices in `{0,1,2,3}` against the direct factorial expansion.
THM-2853 gives a complementary positive cycle model for these tensors;
the present theorem does not use that result as a dependency.

Put

```text
U=d_1+x d_3,                    V=d_2+y d_3.           (4)
```

Let

```text
g11=L_alpha(U^2),    g12=L_alpha(UV),    g22=L_alpha(V^2),

t111=L_alpha(U^3),  t112=L_alpha(U^2V),
t122=L_alpha(UV^2), t222=L_alpha(V^3).                 (5)
```

As in THM-2824, the quadratic moment of `U+zV` divides its cubic
moment exactly when

```text
I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22=0,
I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2=0.          (6)
```

Both rational functions in `(6)` have the same positive denominator

```text
D=alpha^4(alpha+1)^4(alpha+2)^4(alpha+3)^3.            (7)
```

Write

```text
F=D I1,                     H=D I2,
G=30H-(41-6alpha)F.                                      (8)
```

The row operation `(F,H)->(F,G)` has determinant `30`, so it neither
creates nor destroys a common zero.

## 2. An exact moving tube

Take

```text
9/10 <= alpha <= 11/10,

x_c(alpha)=26362/100000+(73/1000)(alpha-1),
y_c(alpha)=23420/1000000+(7/500)(alpha-1),             (9)

|x-x_c(alpha)| <= 1/10000,
|y-y_c(alpha)| <= 1/20000.                             (10)
```

Pull `(8)` back to the unit cube by writing

```text
alpha=9/10+t/5,
x=x_c(alpha)+(2u-1)/10000,
y=y_c(alpha)+(2v-1)/20000.                             (11)
```

Exact power-to-Bernstein conversion proves

```text
F(t,0,v)>0,                 F(t,1,v)<0,
G(t,u,0)>0,                 G(t,u,1)<0,                (12)

F_x<0,                      J=F_x G_y-F_y G_x>0        (13)
```

on the indicated faces or whole cube.  No floating-point inequalities
enter this certificate.  The pulled-back multidegrees are

```text
F       (15,4,3),             G       (16,4,4),
-F_x    (14,3,3),             J       (28,6,6).        (14)
```

Thus `(12)` checks respectively `64,64,85,85` signed face coefficients,
while `(13)` checks `240` and `1421` positive cell coefficients.
The companion separately evaluates the sparse affine pullback and its
Bernstein form at a rational interior point.

Fix `alpha`.  For every `y` in its moving interval, `F_x<0` and the two
vertical signs in `(12)` give exactly one root `x=x(y)`.  Along this
curve,

```text
d/dy G(x(y),y)=J/F_x<0.                                (15)
```

The two horizontal signs in `(12)` therefore give exactly one common
zero

```text
(x(alpha),y(alpha))                                   (16)
```

inside the moving rectangle.  Its Jacobian is nonzero.  The implicit
function theorem and uniqueness glue these points into one real-analytic
branch on the parameter interval (with analytic continuation across each
endpoint).

The whole tube is strictly positive: its smallest possible coordinates
already occur at the left parameter face and satisfy

```text
x >= 25622/100000,             y >= 21970/1000000.     (17)
```

Hence `(4)` consists of linearly independent nonzero positive adjacent
cones.  At `alpha=1`, the THM-2846 isolating rectangle lies inside
`(9)--(10)`, so uniqueness identifies `(16)` with its canonical
factorial hostile.

## 3. Nullity and the uniform fourth exit

For a point on `(16)`, put

```text
Q(z)=L_alpha((U+zV)^2)=q0+q1 z+q2 z^2,
C(z)=L_alpha((U+zV)^3).                                (18)
```

The positivity in `(2)` and linear independence in `(17)` make the
real Gram form `Q` positive definite.  Its two roots are nonreal.
Equations `(6)` say exactly that `Q` divides `C`.  For either root `z`,

```text
K=U+zV !=0,

L_alpha(K)=L_alpha(K^2)=L_alpha(K^3)=0.                (19)
```

The first equality follows from `L_alpha(d_n)=0`.

This null has a uniform exit at the next moment.  Write

```text
L_alpha((U+zV)^4)=k0+k1 z+k2 z^2+k3 z^3+k4 z^4,       (20)
```

where the binomial coefficients are included in the `k_i`.  Modulo
`Q(z)`, the coefficient `b` of `z` in the linear remainder satisfies

```text
q2^3 b =
 k1 q2^3
-k2 q1 q2^2
+k3(q1^2-q0 q2)q2
+k4(2q0 q1 q2-q1^3).                                 (21)
```

After substitution from `(3)`, the right side of `(21)` has denominator

```text
alpha^6(alpha+1)^6(alpha+2)^6(alpha+3)^4>0.            (22)
```

Its numerator pulls back to multidegree `(22,3,7)`, and all `736`
Bernstein coefficients are strictly positive.  Therefore `b>0`
throughout the moving tube.  At the upper-half-plane root of `Q`,

```text
Im L_alpha(K^4)=b Im(z)>0,                             (23)
```

and the lower root gives its nonzero conjugate.  In particular no member
of the branch also kills the fourth moment.

## 4. Meaning and boundary

The theorem proves a structurally stable obstruction to any argument
which tries to detect arbitrary pairs of positive adjacent cones using
only moments through order three.  It is stronger than a qualitative
implicit-function observation: the parameter interval, moving tube,
transversality, and fourth-moment exit are all exact and uniform.

It is not a counterexample to GMC2.  Only `alpha=1` in this interval is
the radial moment law of one standard complex Gaussian.  On that fibre,
THM-2846 converts `(19)` into a degree-seven Gaussian polynomial whose
moments one through six vanish and whose eighth moment is nonzero.
For `alpha!=1`, `(19)` is a positive Gamma-functional deformation, not
an additional two-dimensional Gaussian polynomial.

The result also does not classify all positive common zeros outside the
explicit moving tube.  The word “unique” always refers to the certified
zero inside `(9)--(10)`.

## 5. Exact verification

Run

```text
python 04-computation/gmc_gamma_transverse_null_holotopy_thm2865.py
python -O 04-computation/gmc_gamma_transverse_null_holotopy_thm2865.py
```

Normal, optimized, and stored-output replay must be byte-identical after
LF normalization.  LF-normalized SHA-256:

```text
script  ff32c787fb4fef58b77d06d84f3b032197744515ec03d6dcff8c1be6e70afc52
output  6b943778460db66d48b50fe6d62f94feb6df62bbe08006df287ffcd4ca9ec0ec
```

An independent hostile audit rederived the Gamma normalization, the two
division-free invariants, the affine pullback and separable
power-to-Bernstein formula, all face and cell counts, the Jacobian
orientation, analytic gluing, and the quartic-remainder identity.  It also
replayed normal and optimized execution against the stored transcript and
recovered the displayed LF hashes.  No owner-file edit was needed.

## 6. Connection contract

```text
source:
  THM-2846's isolated positive transverse factorial null;

target:
  a certified section over a nontrivial interval of Gamma radial laws;

map:
  exact Gamma-adjacent tensors, an alpha-dependent invertible row
  operation, and a rational moving-box Bernstein certificate;

preserved:
  positive adjacent-cone coordinates, binary quadratic/cubic
  divisibility, transversality, and nonzero fourth-moment exit;

destroyed / not transported:
  the literal one-complex-dimensional Gaussian interpretation away
  from alpha=1 and any all-moment nullity;

needed sidecar for GMC2:
  an all-order mechanism, not another finite-order continuation;

cheapest decisive extension test:
  enlarge the alpha interval until a face coefficient, Jacobian
  coefficient, quartic-exit coefficient, or positive coordinate reaches
  zero, then refine only the first failing face.
```

No Gaussian Moment Conjecture case is closed by this theorem.
