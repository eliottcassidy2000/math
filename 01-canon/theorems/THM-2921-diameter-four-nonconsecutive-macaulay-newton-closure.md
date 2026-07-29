---
id: THM-2921
title: "Diameter-four nonconsecutive Macaulay--Newton closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  n>=0, first-window SFC(4) holds on each of
  the three nonconsecutive four-slot supports in
  {n,n+1,n+2,n+3,n+4}.  Correct ordinary-monomial multinomial forms,
  exact common denominator clearing, and a fixed degree-seven Macaulay
  minor reduce the infinite depth parameter to three 197-term positive
  Gregory--Newton certificates.  Arbitrary three-slot detection makes
  the four-moment bound sharp with full support; at positive depth the
  corresponding five-monomial two-charge Gaussian bound is exactly eight.
source: codex-gmc-holotopy-extension-2026-07-29
depends_on:
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2849-four-slot-first-window-macaulay-box
  - THM-2908-consecutive-four-slot-projective-resultant-closure
  - THM-2917-all-three-slot-diameter-four-factorial-detection
script: 04-computation/gmc_diameter_four_nonconsecutive_macaulay_newton_thm2921.py
output: 05-knowledge/results/gmc_diameter_four_nonconsecutive_macaulay_newton_thm2921.out
script_sha256: 42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64
output_sha256: dadf97380759e55be9ac84d431806494b972d3b6d0fa3a12fc2fcb319d21fbba
hash_basis: LF-normalized bytes
---

# THM-2921 -- diameter-four nonconsecutive Macaulay--Newton closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let

```text
L:C[s] -> C,                         L(s^j)=j!.          (1)
```

Fix an integer `n>=0` and one of the three offset sets

```text
B_1=(0,1,2,4),       B_2=(0,1,3,4),       B_3=(0,2,3,4). (2)
```

For every nonzero polynomial

```text
H=sum_(b in B_i) c_b s^(n+b),                              (3)
```

at least one of

```text
L(H),                 L(H^2),                 L(H^3), L(H^4) (4)
```

is nonzero.  Equivalently, first-window SFC(4) holds on all three
nonconsecutive translated four-subsets of five consecutive exponents.

The bound is exact on every support in `(2)`: there is a polynomial `H`
with all four displayed coefficients nonzero such that its first three
moments vanish and its fourth does not.

## 2. Correctly typed moment forms

Put

```text
f_a=s^a/a!,                         L(f_a)=1.             (5)
```

Rescaling the four coefficient coordinates by nonzero factorials does
not affect projective nullity.  Since every vector in the mean-zero
hyperplane has a unique expression

```text
H=x(f_(n+b0)-f_(n+4))
 +y(f_(n+b1)-f_(n+4))
 +z(f_(n+b2)-f_(n+4)),                                 (6)
```

where `(b0,b1,b2,4)=B_i`, it is enough to study the homogeneous forms

```text
Q=L(H^2)/L(f_n^2),
C=L(H^3)/L(f_n^3),
F=L(H^4)/L(f_n^4)                                      (7)
```

in `R=Q[x,y,z]`.

The normalized symmetric tensor entry at offsets
`d_1,...,d_m` is

```text
T_m(d_1,...,d_m;n)
 :=L(prod_j f_(n+d_j))/L(f_n^m)
  =(mn+1)_(sum_j d_j)/prod_j (n+1)_(d_j).             (8)
```

If `alpha=(alpha_0,alpha_1,alpha_2)` has total degree `m`, repeat
direction `j` exactly `alpha_j` times.  The ordinary-monomial
coefficient of `x^alpha_0 y^alpha_1 z^alpha_2` in `(7)` is

```text
 m!/(alpha_0!alpha_1!alpha_2!)
 times the signed 2^m sum obtained by replacing any selected
 direction offset by the top offset 4.                         (9)
```

The multinomial factor in `(9)` is load-bearing.  A first scratch
constructor stored one symmetric tensor value per monomial but then
used ordinary polynomial multiplication, omitting the ordered copies.
Already the mixed quadratic coefficient has ratio `2:1` between the
correct and false constructors.  MISTAKE-323 records that invalid
minor and its repair.

## 3. Exact denominator clearing and the fixed Pluecker chart

For all three families in `(2)`, exact division in `Z[n]` proves that
the following common denominators clear every coefficient:

```text
D_2=(n+1)(n+2)(n+3),

D_3=(n+1)^2(n+2)^2(n+3)^2(n+4),

D_4=(n+1)^3(n+2)^3(n+3)^3(n+4)^2.                    (10)
```

The coefficient degrees of the scaled forms

```text
Q~=D_2 Q,                   C~=D_3 C,              F~=D_4 F
```

are at most

```text
3,                            7,                       11. (11)
```

All factors in `(10)` are nonzero at integer `n>=0`, so scaling does
not change the common projective zero set.

Let `R_d` denote the degree-`d` part of `R`.  As in THM-2849, use the
degree-seven Macaulay map

```text
Phi_n:R_5 direct_sum R_4 direct_sum R_3 -> R_7,
Phi_n(A,B,D)=A Q~+B C~+D F~.                            (12)
```

In the canonical monomial order, its stored transpose has

```text
21+15+10=46 rows,                       36 columns.     (13)
```

Select rows

```text
0,...,19;
21,...,29,35;
36,...,41.                                               (14)
```

Thus the selected square minor uses `20` quadratic, `10` cubic, and
`6` quartic rows.  Call its determinant `P_B(n)`.  Equations
`(11),(14)` give the rigorous degree invoice

```text
deg P_B <=20*3+10*7+6*11=196.                          (15)
```

This one row set works for every depth in every family in `(2)`.

## 4. Gregory--Newton nonvanishing

The companion evaluates each integer determinant exactly at the `197`
consecutive depths needed by `(15)` and takes exact forward differences.
The complete coefficient vectors are bound by these SHA-256 digests:

| offsets `B` | GN base | sign `P_B(0)` | scaled-form digest | GN-vector digest |
|---|---:|:---:|---|---|
| `(0,1,2,4)` | `1` | `-` | `8c133750387b4e2c285b1f1547e22b97dd5ba2851ddb0912568b55716f5bb936` | `222ef28eecedc4a59e33d13ff0a5bd9d45f218be0d1f3ad63681467cfe6848e2` |
| `(0,1,3,4)` | `1` | `-` | `6cfc17b9b99fe3eeffd549fedc0d2caa057f444ba873c410d13c1301f8d6d38a` | `beae481ee3937123a33f995469e856684988cc83b41d0c9fa08d7c12e2482817` |
| `(0,2,3,4)` | `0` | `+` | `2ba15eb5005f38bf4ffcca3124f810edca55f4d6238168cbdffdc824817ada88` | `7ee5fe137ad12cb1dc3874213b95ecce02d77ce8a44137ac7f79357b95ba33e0` |

Every one of the `197` forward differences in each row is strictly
positive.  For any polynomial `P` of degree at most `196`,

```text
P(n)=sum_(k=0)^196 Delta^k P(a) binom(n-a,k),
                                                        n>=a. (16)
```

Consequently the first two determinants are positive at every integer
`n>=1`; their separately computed nonzero negative value at `n=0`
closes depth zero.  The third determinant is positive for every integer
`n>=0` directly from base zero.  Therefore

```text
P_B(n)!=0                  for every B in (2), n>=0.   (17)
```

## 5. Projective consequence

Equation `(17)` gives `rank Phi_n=36`, so

```text
(Q~,C~,F~)_7=R_7.                                      (18)
```

If `[x:y:z] in C P^2` were a common zero, every form in the left side
of `(18)` would vanish there.  At least one coordinate is nonzero, but
its seventh power belongs to `R_7`, a contradiction.  Hence `(7)` has
no common nonzero projective zero.  Together with the mean elimination
`(6)`, this proves `(4)`.

## 6. Full-support sharpness and the Gaussian corollary

For a fixed four-slot envelope, THM-2173 gives a nonzero `H` such that

```text
L(H)=L(H^2)=L(H^3)=0.                                  (19)
```

Such an `H` cannot lie on a coordinate hyperplane: every coordinate
face is a three-slot support, and PROVED THM-2824 detects every arbitrary
three-slot polynomial by one of its first three moments.  Thus the
witness in `(19)` has all four coefficients nonzero.  The theorem just
proved forces

```text
L(H^4)!=0.                                             (20)
```

This proves exact full-support four-moment sharpness on all three
families.

If `n>=1`, write `H=s h(s)`, put `s=ZW` for a standard complex Gaussian
`Z` with `W=conj(Z)`, and choose `alpha!=0`.  For

```text
P=alpha W+Z h(ZW),                                     (21)
```

charge balance gives

```text
E[P^(2j+1)]=0,
E[P^(2j)]=binom(2j,j) alpha^j L(H^j).                 (22)
```

Therefore every genuinely two-charge lift `(21)` with nonzero radial
part is detected by one of the moments of orders `2,4,6,8`.  Applying
`(22)` to `(19),(20)` gives a five-monomial, full-radial-support
polynomial whose moments one through seven vanish and whose eighth is

```text
E[P^8]=70 alpha^4 L(H^4)!=0.                           (23)
```

Thus the exact uniform Gaussian detection depth inside each of these
positive-depth two-charge envelopes is eight.

## 7. Pluecker/holotopy interpretation and scope

The selected determinant is a Pluecker coordinate of the rank-`36`
Macaulay image.  In these three translated cells it never vanishes, so
one chart trivializes the whole discrete depth ray.  The invariant is
rank/projective emptiness; the chosen row address is only a chart.

This is the useful transfer from the endpoint cospan and aggregate-first
logic of THM-2452/2697.  On a later gap family, the vanishing of this
particular minor would mean only that the path left this chart.  The
correct next object would be a finite minor atlas, a no-common-zero
certificate between minors, or a positive Gram sum of minors—not a
claim of rank loss.

The proved scope here is exactly

```text
slots:       four;
offsets:     the three sets in (2);
window:      moments one through four;
depth:       every integer n>=0;
Gaussian:    n>=1, two charges, five monomials.         (24)
```

No arbitrary four-slot support, shifted moment window, SFC(5), full
Strong Factorial Conjecture, arbitrary-charge GMC(2), Gaussian nullcone
classification, or Jacobian-conjecture conclusion is asserted.  The
two consecutive four-subsets of a five-point window belong to THM-2908,
not to this proof.

## 8. Exact verification

The exact companion:

1. constructs every ordinary coefficient from `(8),(9)`;
2. proves every division in `(10)` inside `Z[n]`;
3. checks `(11)` and the degree invoice `(15)`;
4. computes all three fixed determinants with exact FLINT integer
   arithmetic and verifies all `591` positive Gregory--Newton
   coefficients;
5. checks the original unscaled depth-zero determinant against the
   scaled one; and
6. independently reconstructs the forms by direct four-variable
   multinomial expansion, then reproduces the selected determinant
   modulo `1000003` at seven hostile depths in each family.

The last route makes `21` exact cross-constructor minor checks and does
not use the rising-factorial form constructor.  The mixed-quadratic
`2:1` hostile freezes the MISTAKE-323 boundary.

Install the already pinned exact dependency with

```text
python -m pip install -r 04-computation/requirements-gmc-projective-resultant.txt
```

and run

```text
python 04-computation/gmc_diameter_four_nonconsecutive_macaulay_newton_thm2921.py
python -O 04-computation/gmc_diameter_four_nonconsecutive_macaulay_newton_thm2921.py
```

Normal and optimized executions byte-match the stored output with the
declared LF-normalized hashes.

An independent hostile audit reconstructed the ordinary coefficient forms
and the selected minor by a separate exact four-variable route at `12`
additional controls.  It rederived the denominator degrees, row allocation,
degree bound, all three Gregory--Newton base/sign arguments, Macaulay
surjectivity, projective emptiness, full-support sharpness, and the Gaussian
lift.  It also replayed normal, optimized, and stored output and reproduced
the declared LF-normalized hashes.  No defect remained.

**QED.**
