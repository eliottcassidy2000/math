---
id: THM-4028
title: "Sun 2-4-6-8 average-order and residue-class criticality"
status: >
  PROVED + VERIFIED-EXACT FINITE CONTROLS + INDEPENDENTLY HOSTILE-AUDITED.
  The cumulative representation count is asymptotic to
  V X^(25/24), V=24.31102486226095468...; for every fixed q and r the
  r-class constant is (sigma_q(r)/q)V. The mean exponent is only 1/24.
  No pointwise asymptotic, Poisson law, or eventual coverage is claimed.
source: root + independent anisotropic and residue-class audits, 2026-08-24
depends_on:
  - THM-4027-sun-two-four-six-eight-universal-modular-solubility
related:
  - THM-4026-sun-two-four-six-eight-binomial-counterexample
  - THM-2412-delta-exponential-and-central-newton-layer-split
script: 04-computation/sun_2468_average_order_thm4028.py
independent_audit_script: 04-computation/sun_2468_average_order_independent_hostile_audit_thm4028.py
independent_audit_output: 05-knowledge/results/sun_2468_average_order_independent_hostile_audit_thm4028.out
residue_script: 04-computation/sun_2468_residue_class_average_order_thm4028.py
residue_output: 05-knowledge/results/sun_2468_residue_class_average_order_thm4028.out
---

# THM-4028 -- barely-supercritical average order

**PROVED + VERIFIED-EXACT FINITE CONTROLS + INDEPENDENTLY HOSTILE-AUDITED.**
Let `a(n)` be the number of quadruples

```text
w>=2, x>=3, y>=5, z>=7                              (1)
```

such that

```text
C(w,2)+C(x,4)+C(y,6)+C(z,8)=n.                       (2)
```

Put

\[
S={1\over2}+{1\over4}+{1\over6}+{1\over8}={25\over24}. \tag{3}
\]

Then

\[
\sum_{n\le X}a(n)\sim V X^{25/24},                   \tag{4}
\]

where

\[
V={\prod_{k\in\{2,4,6,8\}}(k!)^{1/k}\Gamma(1+1/k)
\over\Gamma(1+25/24)}
=24.311024862260954682905938684\ldots .                \tag{5}
\]

Consequently

\[
{1\over X}\sum_{n\le X}a(n)\sim V X^{1/24}.          \tag{6}
\]

More strongly, for every fixed positive modulus `q` and residue `r mod q`,

\[
\sum_{\substack{n\le X\\n\equiv r\pmod q}}a(n)
\sim {\sigma_q(r)\over q}V X^{25/24},                 \tag{7}
\]

where `sigma_q(r)` is the exact complete-period local factor of THM-4027.

## 1. Why this is the right frontier quantity

Sun already highlighted the reciprocal sum `(3)` in the original
[MathOverflow question](https://mathoverflow.net/questions/323541/positive-integers-written-as-binomw2-binomx4-binomy6-binomz8-with).
It is only `1/24` above the volume threshold one. THM-4027 shows that the same
excess makes Cauchy--Davenport saturate every sufficiently large prime. The
connection is exact:

```text
source: four growth degrees (2,4,6,8)
target A: anisotropic lattice volume
target B: finite-field image-set sizes
preserved invariant: sum reciprocal degrees = 25/24
lost data: pointwise alignment, minor-arc cancellation, and exact fibres.   (8)
```

The closest repo coordinate is THM-2412's binomial/falling-factorial basis.
The canonical hostile is THM-4026 itself: growing average mass and universal
local solubility do not prevent an exact zero.

## 2. Anisotropic Riemann-sum proof

The left side of `(4)` counts exactly the lattice points in

\[
\mathcal R_X=\left\{(w,x,y,z)\text{ satisfying `(1)`}:
\sum_{k}\binom{t_k}{k}\le X\right\}.                  \tag{9}
\]

For fixed `k`, uniformly on `t=O(X^(1/k))`,

\[
{1\over X}\binom tk={1\over k!}\left({t\over X^{1/k}}\right)^k
+O(X^{-1/k}).                                           \tag{10}
\]

Scale `t_k=X^(1/k)u_k`. The grid cell has volume `X^-S`, and `(10)` makes the
scaled regions converge from inside and outside to the compact monotone body

\[
\mathcal R=\left\{u_k\ge0:\sum_k {u_k^k\over k!}\le1\right\}. \tag{11}
\]

Its algebraic boundary has Lebesgue measure zero. Upper and lower Riemann
sums therefore give

```text
#R_X / X^S -> vol(R).                                  (12)
```

The finite lower-index shifts in `(1)` remove only boundary faces and do not
change the limit.

To evaluate the volume, set `v_k=u_k^k/k!`. Then

\[
du_k=(k!)^{1/k}{1\over k}v_k^{1/k-1}dv_k,              \tag{13}
\]

and the Dirichlet simplex integral gives

\[
\operatorname{vol}(\mathcal R)
={\prod_k (k!)^{1/k}\Gamma(1+1/k)\over\Gamma(1+S)},    \tag{14}
\]

which is `(5)` and proves `(4)--(6)`.

## 3. The fixed-residue refinement

Fix `q`. For each role `k`, let `P_(q,k)` be the exact binomial period proved
in THM-4027. Split the `k`th input axis into its finitely many congruence
classes modulo `P_(q,k)`. Every fixed product class is an anisotropic
sublattice of index

```text
D_q=product_k P_(q,k),                                  (15)
```

and the same Riemann-sum sandwich gives it leading constant `V/D_q`.

Let `C_q(r)` be the number of product classes whose four binomial values sum
to `r mod q`. By definition,

```text
sigma_q(r)=q*C_q(r)/D_q.                                (16)
```

Summing `V/D_q` over those `C_q(r)` sublattices proves `(7)`. Since there are
`X/q+O(1)` target integers in the residue class, their mean representation
count is

\[
\sigma_q(r)V X^{1/24}(1+o(1)).                         \tag{17}
\]

Thus a local factor is rigorously an **average-order suppression factor**.
It is not a pointwise probability.

For THM-4026's target,

```text
N=20 mod 33,       sigma_33(N)=16/33,
N=53 mod 99,       sigma_99(N)=544/1089.               (18)
```

The class `20 mod 33` has only `16/33` of the uniform leading average mass,
which explains why the 2019 residue observation was a useful search prior.
It does not explain the exact zero by itself.

## 4. Exact finite controls

The primary companion counts `(9)` with integer arithmetic. Its convergence
table is

| `X` | exact cumulative count | count / `X^(25/24)` | fraction of `V` |
|---:|---:|---:|---:|
| `10^4` | `210751` | `14.358298484603` | `0.590608522921` |
| `10^5` | `2715445` | `16.807676380792` | `0.691360256345` |
| `10^6` | `33359460` | `18.759402944034` | `0.771641798333` |
| `10^7` | `395630618` | `20.212648692315` | `0.831419029302` |
| `10^8` | `4588576456` | `21.298285240348` | `0.876075170052` |
| `10^9` | `52415036766` | `22.103237731303` | `0.909185764752` |

The residue companion independently partitions the exact counts modulo 99
through `X=10^7`. At `X=10^7` it finds

```text
total=395630618,
class 53=2104568,      class 86=1951931,               (19)
```

with normalized ratios `0.8765` and `0.8916` toward the respective constants
in `(7)`. Normal and optimized runs retain all gates and agree byte-for-byte.

## 5. The hostile boundary

The formal derivative of the main term in `(4)` has constant

\[
S V=25.323984231521827794693686129\ldots .             \tag{20}
\]

This suggests an archimedean shell scale, but differentiating an average
asymptotic is not a pointwise theorem. Nothing here proves

```text
a(n) ~ S V n^(1/24) * an infinite local product,        (21)
```

a Poisson law, a zero-frequency estimate, or eventual representation.

The neighboring degree packet is a decisive hostile. Its reciprocal sum is

\[
{1\over2}+{1\over4}+{1\over6}+{1\over9}={37\over36}>1, \tag{22}
\]

yet exact enumeration gives

```text
1061619 != C(w,2)+C(x,4)+C(y,6)+C(z,9)                 (23)
```

even when zero triangular and higher atoms are allowed in the minimal-domain
convention. The companion checks all `28,037` feasible higher-term triples;
both neighbours are independently represented. Thus “supercritical” means
growing cumulative mass and eventual large-prime saturation, not a pointwise
additive basis theorem.

## 6. Scope and reproduction

Run

```text
python -B 04-computation/sun_2468_average_order_thm4028.py
python -B -O 04-computation/sun_2468_average_order_thm4028.py
python -B 04-computation/sun_2468_residue_class_average_order_thm4028.py
python -B -O 04-computation/sun_2468_residue_class_average_order_thm4028.py
python -B 04-computation/sun_2468_average_order_independent_hostile_audit_thm4028.py
python -B -O 04-computation/sun_2468_average_order_independent_hostile_audit_thm4028.py
```

The proof of `(4)` and `(7)` is the Riemann-sum argument above; the scripts
are exact finite controls, not extrapolations. THM-4028 says nothing about
the least counterexample, the density of zero fibres, or a Hardy--Littlewood
minor-arc estimate. **QED.**
