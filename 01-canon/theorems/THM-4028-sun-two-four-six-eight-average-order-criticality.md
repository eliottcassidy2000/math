---
id: THM-4028
title: "Sun 2-4-6-8 average-order and residue-class criticality"
status: >
  PROVED + VERIFIED-EXACT FINITE CONTROLS + INDEPENDENTLY PROOF-AUDITED.
  The cumulative representation count is V X^(25/24)+O(X^(11/12)),
  V=24.31102486226095468...; for every fixed q and r the r-class main
  constant is (sigma_q(r)/q)V with the same O_q(X^(11/12)) error. The mean
  exponent is only 1/24, and the shell law is proved only on intervals
  H=o(X), H/X^(7/8)->infinity. No pointwise asymptotic or coverage follows.
source: root + independent anisotropic, residue-class, and hostile-boundary audits, 2026-08-24
audit: >
  PASS. An independent proof audit checked the half-open unit-cube squeeze,
  shifted binomial faces, exponent 11/12, fixed-period cosets, constants V
  and J, and the mesoscopic threshold H/X^(7/8)->infinity. Constants may
  depend on fixed q; no uniformity for growing q is asserted.
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
quantitative_proof_audit: 05-knowledge/results/sun_2468_average_order_quantitative_proof_audit_thm4028.md
---

# THM-4028 -- barely-supercritical average order

**PROVED + VERIFIED-EXACT FINITE CONTROLS + INDEPENDENTLY PROOF-AUDITED.**
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
\sum_{n\le X}a(n)=V X^{25/24}+O(X^{11/12}),          \tag{4}
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
= {\sigma_q(r)\over q}V X^{25/24}+O_q(X^{11/12}),     \tag{7}
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

## 2. Quantitative anisotropic lattice proof

The useful statement is general. Fix positive degrees `d_1,...,d_s` and
fixed nonnegative lower bounds `L_i`, and put

\[
A_{\mathbf d,\mathbf L}(X)=
\#\left\{t_i\ge L_i:\sum_i\binom{t_i}{d_i}\le X\right\}. \tag{9}
\]

Let `sigma=sum_i 1/d_i`, `d_max=max_i d_i`, and

\[
V_{\mathbf d}=
{\prod_i(d_i!)^{1/d_i}\Gamma(1+1/d_i)\over
 \Gamma(1+\sigma)}.
\]

Then

\[
A_{\mathbf d,\mathbf L}(X)=V_{\mathbf d}X^\sigma
+O_{\mathbf d,\mathbf L}\!\left(X^{\sigma-1/d_{\max}}\right). \tag{10}
\]

For `t>=0`, exact integer binomials satisfy

\[
{(t-d+1)_+^d\over d!}\le\binom td\le {t^d\over d!}.   \tag{11}
\]

First count lattice points in the monotone pure-power body

\[
\mathcal R_X=\left\{u_i\ge0:\sum_i{u_i^{d_i}\over d_i!}\le X\right\}. \tag{12}
\]

Use half-open cubes `p+[0,1)^s`. If `u` lies in `R_X`, then `floor(u)` also
lies there. Conversely, every cube based at a lattice point of `R_X` lies in
`R_(X+C X^(1-1/d_max))`: increasing coordinate `i` by at most one changes
its pure-power term by `O(X^(1-1/d_i))`. The cubes are disjoint and have unit
volume, so scaling the body gives

\[
\#(\mathcal R_X\cap\mathbf Z^s)
=\operatorname{vol}(\mathcal R_1)X^\sigma
+O\!\left(X^{\sigma-1/d_{\max}}\right).               \tag{13}
\]

The upper inequality in `(11)` embeds every pure-power tuple into the
binomial sublevel set. In the reverse direction, map
`t_i` to `(t_i-d_i+1)_+`. Away from coordinate faces this is a fixed
translation and is injective; duplicate zero fibres, the fixed lower bounds,
and the finitely thick coordinate faces contribute at most

```text
O(sum_i X^(sigma-1/d_i))=O(X^(sigma-1/d_max)).
```

This proves the squeeze `(10)`. Finally, substituting
`v_i=u_i^{d_i}/d_i!` in `(12)` gives the Dirichlet volume

\[
\operatorname{vol}(\mathcal R_1)
={\prod_i(d_i!)^{1/d_i}\Gamma(1+1/d_i)\over
  \Gamma(1+\sigma)}=V_{\mathbf d}.                    \tag{14}
\]

Taking degrees `(2,4,6,8)` and the minimal bounds `(1)` makes
`sigma=25/24`, `d_max=8`, and `(10)` exactly `(4)`. Allowing every literal
index `>=2` changes only fixed zero-coordinate fibres and hence only the same
error term; the theorem keeps the OEIS minimal-domain multiplicity convention.

The general estimate also makes the reciprocal-degree threshold precise. If
`sigma<1`, the support below `X` has size `O(X^sigma)=o(X)` and cannot be
cofinite. If `sigma=1` and `V_d<1`, even the total tuple count is
at most `(1-epsilon)X` for some `epsilon>0` and every sufficiently large
`X`, so cofiniteness is again impossible. When `sigma>1`,
only average capacity grows: THM-4026 is the hostile showing that this
necessary capacity coordinate is not sufficient for universal coverage.

## 3. The fixed-residue refinement

Fix `q`. For each role `k`, let `P_(q,k)` be the exact binomial period proved
in THM-4027. Split the `k`th input axis into its finitely many congruence
classes modulo `P_(q,k)`. Every fixed product class is an anisotropic
sublattice of index

```text
D_q=product_k P_(q,k),                                  (15)
```

and the same lattice squeeze gives it leading constant `V/D_q`.
Quantitatively, tile each product coset by half-open boxes with side lengths
`P_(q,k)`. Dividing their union volume by `D_q` gives the coset count; the
translated orthant misses only fixed-width slabs along coordinate faces, and
inflating a box changes the body height by `O_q(X^(7/8))`. Moreover, away
from the zero-coordinate faces, the upper binomial map
`t_k -> t_k-k+1` sends a fixed input coset to a translated product coset.
Thus the face and zero-fibre losses remain `O_q(X^(11/12))`, and each coset is

```text
(V/D_q) X^(25/24)+O_q(X^(11/12)).
```

Let `C_q(r)` be the number of product classes whose four binomial values sum
to `r mod q`. By definition,

```text
sigma_q(r)=q*C_q(r)/D_q.                                (16)
```

Summing over the finitely many accepted sublattices proves `(7)`, including
its displayed error term. Equivalently one may use the safe periods
`q*k!`; exact minimality is not needed here. The implied constant is allowed
to depend on fixed `q`, and no uniformity for a growing modulus is claimed.
Since there are
`X/q+O(1)` target integers in the residue class, their mean representation
count is

\[
\sigma_q(r)V X^{1/24}(1+o(1)).                         \tag{17}
\]

THM-4027 gives `sigma_q(r)>0`, so the relative-error notation in `(17)` is
legitimate.

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

The primary companion counts the cumulative minimal-domain tuples with
integer arithmetic. Its convergence table is

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

## 5. Mesoscopic shell law and the hostile boundary

The derivative constant of the main term in `(4)` is

\[
J=S V=25.323984231521827794693686129\ldots .           \tag{20}
\]

The error term makes this a theorem on sufficiently long intervals. If
`H=o(X)` and `H/X^(7/8)->infinity`, then

\[
{1\over H}\sum_{X<n\le X+H}a(n)\sim JX^{1/24}.         \tag{21}
\]

Indeed, the main increment is `J H X^(1/24)(1+o(1))`, whereas the two
summatory errors total `O(X^(11/12))`; their ratio is
`O(X^(7/8)/H)`. For every fixed `q,r`, the same subtraction in `(7)` gives

\[
{1\over H}\sum_{\substack{X<n\le X+H\\n\equiv r\pmod q}}a(n)
\sim {\sigma_q(r)\over q}JX^{1/24}.                    \tag{22}
\]

Equivalently, after dividing by the approximately `H/q` eligible targets,
the residue-class mean is `sigma_q(r) J X^(1/24)`. This remains a
mesoscopic statement. At `H=1` the error dominates, so nothing here proves

```text
a(n) ~ J n^(1/24) * an infinite local product,          (23)
```

a Poisson law, a zero-frequency estimate, or eventual representation.

The neighboring degree packet is a decisive hostile. Its reciprocal sum is

\[
{1\over2}+{1\over4}+{1\over6}+{1\over9}={37\over36}>1, \tag{24}
\]

yet exact enumeration gives

```text
1061619 != C(w,2)+C(x,4)+C(y,6)+C(z,9)                 (25)
```

under the one-zero-representative convention for the higher atoms, even after
triangular zero is separately adjoined. The companion checks all `28,037`
feasible higher-term triples; both neighbours are independently represented.
Thus “supercritical” means
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

The proof of `(4)` and `(7)` is the quantitative lattice squeeze above; the
scripts are exact finite controls, not extrapolations. THM-4028 says nothing
about the least counterexample, the density of zero fibres, or a Hardy--Littlewood
minor-arc estimate. **QED.**
