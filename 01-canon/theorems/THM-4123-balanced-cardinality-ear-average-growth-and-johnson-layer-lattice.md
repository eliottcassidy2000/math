---
id: THM-4123
title: "Balanced-cardinality ear average, growth, and Johnson-layer lattice"
status: >
  PROVED ELEMENTARY FIXED-CARDINALITY DOUBLE COUNT + JOHNSON-LATTICE +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every cardinality-m cut-ear fibre
  has an exact binomial sum and mean. Central layers give a balanced-cut
  selected-bank recurrence that is coefficientwise stronger than THM-4115's
  reduced Cauchy recurrence for every n>=4. Its mean-only floor is
  incomparable with THM-4115's unreduced variance floor; THM-4127 later
  restores slice variance and proves central support dominance. Exchange
  increments generate the exact response lattice on each Johnson slice and
  sharpen mean rounding. This forces one large balanced child, not an
  interval or a balanced global maximizer.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-001-redei
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4115-uniform-ear-cut-walsh-variance-and-sharp-growth-refinement
  - THM-4118-ear-response-lattice-and-stateful-unit-component-intervals
related:
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - THM-4127-johnson-slice-hoeffding-variance-and-central-support-dominance
  - HYP-2879-strong-ear-atom-calculus
  - HYP-9029-strong-interval-tiling-law
script: 04-computation/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123.py
output: 05-knowledge/results/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123.out
independent_audit_script: 04-computation/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123_independent_audit.out
script_sha256: 39092103a67b33b92d3f434c585b2b3f4c9c3dae003fa671b899c20dd8b7790f
output_sha256: e9263501d1d8129a8460a0deacbdac5d12e528a7042dfb3e484dfba683cfe28e
semantic_sha256: 0f227c83e27b0f2beda7d432c46c2821318cf5cd1277c4c6bafd6e1de67963ec
independent_audit_script_sha256: 7221ca2a46291828333e487d67b112750bb52c2260683afc80aa9dd3595c174a
independent_audit_output_sha256: ee0d0bfa2aac74e39e90a623430df9834afa175805b4804540c92f6cf93926d3
independent_semantic_sha256: 0f227c83e27b0f2beda7d432c46c2821318cf5cd1277c4c6bafd6e1de67963ec
independent_extended_semantic_sha256: b0ee0ce581dccf3af6ba937473c45717b6132f7a396a2dee36be8b62f0efd178
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone subset Hamiltonian-path DP checks 33,864 literal ears,
  6,502 cardinality layers, and 124,468 old-vertex orderings over every 1,098
  labelled parent through order five. It verifies both slice formulas, every
  exchange increment, global-to-layer lattice divisibility, coset rounding,
  the recurrence comparison, and seven named boundary controls. Normal,
  optimized, and frozen outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room Start/End/Q boundary implementation imports no primary
  code and reproduces the shared semantic ledger, with literal child-ordering
  replays at C3 and codes 20 and 30571. It additionally scans all 22,320
  strong labelled order-six parents: 1,920 have no balanced global maximizer,
  and the balanced-layer lattice histogram is 22,080 at 2 and 240 at 6.
  Normal, optimized, and frozen outputs byte-match; no theorem failure occurs.
---

# THM-4123 -- balanced-cardinality ear averages and layer lattices

**PROVED ELEMENTARY FIXED-CARDINALITY DOUBLE COUNT + JOHNSON-LATTICE +
VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4111 averaged the full Boolean cut cube. THM-4115 restored its Walsh
variance. A different operation is unexpectedly stronger for recursive
growth: condition on the cut cardinality before averaging. The zero-sum
field disappears on every slice, while the probability of crossing an edge
is largest at a balanced cut. The same restriction exposes a new arithmetic
sidecar, the exact response lattice of the Johnson graph.

## 1. Exact fixed-cardinality sum and mean

Let `T` be a tournament on `V={1,...,n}`, `n>=2`, and put

```text
H=H(T),
F_T(S)=H(T+x_S),                 S subseteq V,             (1)
```

with THM-4097's convention `x->v` exactly when `v in S`. Let `F_1(T)` count
the old-vertex orderings with exactly one bad adjacency.

> **Theorem 1 (cardinality-slice identity).** For every `0<=m<=n`, with an
> out-of-range binomial coefficient interpreted as zero,
>
> ```text
> sum_(|S|=m) F_T(S)
>   =C(n,m)H+C(n-2,m-1)((n-1)H+F_1(T)).                   (2)
> ```
>
> Consequently, for uniform `S` of size `m`,
>
> ```text
> E_m F_T
>   =H+[2m(n-m)/(n(n-1))]W                               (3)
>   =[1+m(n-m)/n]H+[m(n-m)/(n(n-1))]F_1(T),              (4)
> ```
>
> where
>
> ```text
> W=((n-1)H+F_1(T))/2.                                   (5)
> ```

### Proof by exposed old-vertex words

Fix `|S|=m` and sum Hamiltonian paths over all such ears. If `x` is first,
deletion leaves a Hamiltonian path of `T`, and its first old vertex must lie
in `S`. This gives `C(n-1,m-1)H` paths. If `x` is last, the last old vertex
must lie outside `S`, giving `C(n-1,m)H`. Their sum is `C(n,m)H`.

If `x` is internal, deletion exposes an ordered adjacent pair `(a,b)`. The
old word is either a Hamiltonian path, with any of its `n-1` gaps exposed, or
an `F_1` word at its unique bad gap. In either case insertion requires

```text
a notin S,                 b in S.                         (6)
```

The remaining `m-1` members of `S` can be chosen in `C(n-2,m-1)` ways. The
number of exposed words is `(n-1)H+F_1`, proving `(2)`. Division by `C(n,m)`
and

```text
C(n-2,m-1)/C(n,m)=m(n-m)/(n(n-1))                        (7)
```

give `(4)`.

Equivalently, THM-4115 writes

```text
F_T(S)=H+sum_(i in S)h_i+cut_w(S),
sum_i h_i=0,                  sum_(i<j)w_ij=W.             (8)
```

On a uniform `m`-subset the field expectation is zero, and each unordered
edge crosses the cut with probability `2m(n-m)/(n(n-1))`. This independently
gives `(3)` and identifies `(5)`. **QED.**

## 2. The balanced bound and its exact comparison

Put

```text
q_n=floor(n^2/4),
B_n(T)=[1+q_n/n]H+[q_n/(n(n-1))]F_1(T).                   (9)
```

Because `W>0`, equation `(3)` is maximized exactly when
`m(n-m)=q_n`, namely on the central cardinality layer or layers.

> **Corollary 2 (balanced child).** Each central layer contains a cut `S`
> satisfying
>
> ```text
> F_T(S)>=oddceil(B_n(T)),                                 (10)
> ```
>
> where `oddceil` is the least odd integer at least its argument. Every such
> cut is nonconstant. If `T` is strong, the child is strong.

The mean itself is

```text
B_n(T)=
  (n+4)H/4+nF_1/[4(n-1)]                    n even,
  (n^2+4n-1)H/(4n)+(n+1)F_1/(4n)            n odd.         (11)
```

Compare this with THM-4115's reduced Cauchy floor

```text
R_n(T)=[(n+1)(n+2)/(4n)]H
       +[n(n-1)+2]/[4n(n-1)] F_1(T).                      (12)
```

Exact subtraction gives

```text
B_n-R_n=
  (n-2)/(4n) [H+F_1/(n-1)]                  n even,
  (n-3)/(4n) [H+F_1/(n-1)]                  n odd.         (13)
```

Thus the balanced bound is coefficientwise equal to `(12)` at `n=2,3` and
strictly stronger for every `n>=4`.

This comparison is deliberately with the **reduced** THM-4115 bound. Let

```text
E_n=H+W/2+mathcal_E/(2W),
mathcal_E=sum_i h_i^2+sum_(i<j)w_ij^2                     (14)
```

be its unreduced variance floor. Then

```text
B_n-E_n=
  [W^2/(n-1)-mathcal_E]/(2W)                 n even,
  [W^2/n-mathcal_E]/(2W)                     n odd.        (15)
```

Either sign occurs. A full-cut bank may therefore use
`oddceil(max(B_n,E_n))`; a balanced-only bank gets `(10)` without needing a
slice-variance theorem.

## 3. Selection-robust balanced-only growth

Start with a nonempty bank `mathcal B_n` of strong order-`n` tournaments. Expand only
the cuts in the central cardinality layer or layers of every retained parent,
then retain an arbitrary strong witness for each attained scalar value. Put

```text
M_n=max{H(T):T in mathcal B_n}.                            (16)
```

If `T_n` is a retained maximizing parent, Corollary 2 gives

```text
M_(n+1)>=oddceil(
  [1+q_n/n]M_n+[q_n/(n(n-1))]F_1(T_n)).                   (17)
```

This remains valid for the full-cut banks of THM-4111/4115 because they
contain the central layers. Using THM-4115's strong-parent floor

```text
eta_n=n-1 for even n,
eta_n=2   for odd n,                         n>=4,          (18)
```

gives the selection-independent refinement

```text
M_(n+1)>=oddceil([1+q_n/n]M_n+q_n eta_n/[n(n-1)]).        (19)
```

The displayed additive surplus is `n/4` in even order and `(n+1)/(2n)` in
odd order. At `n=3`, the strong parent is `C3`, where `F_1=0` and all six
central ears have value `5`; the factor `5/3` is exact. Ignoring helpful
surplus and odd rounding gives

```text
M_m>=M_n product_(r=n)^(m-1)[1+floor(r^2/4)/r].           (20)
```

This tends to infinity and strictly improves the reduced THM-4115 product
from every step `r>=4`, while expanding only balanced cuts.

## 4. The exact Johnson-layer response lattice

Write THM-4118's observable quadratic response as

```text
F_T(S)=H+sum_(i in S)A_i-sum_({i,j} subseteq S)K_ij,      (21)
A_i=F_T({i})-H,
K_ij=F_T({i})+F_T({j})-H-F_T({i,j}).                      (22)
```

Fix `1<=m<=n-1`. For distinct `i,j` and
`R subseteq V minus {i,j}`, `|R|=m-1`, define the exchange increment

```text
Delta_(i,j;R)
 =F_T(R union {j})-F_T(R union {i})
 =A_j-A_i-sum_(r in R)(K_jr-K_ir).                        (23)
```

Put

```text
d_(T,m)=gcd{Delta_(i,j;R)}.                               (24)
```

The gcd of an all-zero family is `0`.

> **Theorem 3 (Johnson-layer lattice).** For any fixed `m` and any
> `m`-subset `S_0`,
>
> ```text
> <F_T(S)-F_T(S_0):|S|=m>=d_(T,m) Z.                      (25)
> ```
>
> Thus a zero `d_(T,m)` means the layer is constant. Otherwise every layer
> response lies in the single coset
>
> ```text
> F_T(S_0)+d_(T,m)Z,                                      (26)
> ```
>
> and `d_(T,m)` is even.

Every generator in `(24)` is a layer difference. Conversely, the Johnson
graph `J(n,m)` is connected, and every path edge is exactly one exchange
`(23)`. Any two layer responses therefore differ by a sum of these
increments. This proves `(25)`. Hamiltonian-path counts are odd, proving the
evenness assertion. **QED.**

Let `ceil_(a+dZ)(y)` denote the least member of the coset `a+dZ` at least
`y`. Equations `(4)` and `(25)` sharpen averaging to

```text
max_(|S|=m)F_T(S)
 >=ceil_(F_T(S_0)+d_(T,m)Z)(E_m F_T).                     (27)
```

If the layer lattice is zero, its unique value equals its mean. The global
THM-4118 lattice

```text
d_T=gcd(A_i,K_ij)                                         (28)
```

divides every `d_(T,m)`, so a label-free weaker corollary is

```text
max_(|S|=m)F_T(S)
 >=H+d_T ceil((E_m F_T-H)/d_T).                           (29)
```

The cardinality coordinate can genuinely strengthen `(28)`: for code
`30571`, `d_T=6`, while the layer lattices from sizes zero through six are

```text
(0,18,6,6,6,18,0).                                       (30)
```

Its balanced mean is `738/5`; ordinary odd rounding gives `149`, whereas
the global and balanced-layer coset gives `153`. The actual balanced maximum
is `189`.

## 5. Sharp hostiles and finite boundary

The following rows separate the preserved statistics from what the slice
average destroys.

- Strong order-five codes `8` and `759` both have
  `(H,F_1,W)=(9,30,33)` and balanced mean `144/5`, but their balanced images
  are respectively

  ```text
  {17,19,23,25,27,29,33,37,41},
  {15,19,23,29,31,33,37,43}.                              (31)
  ```

  Hence the exact slice mean does not determine its image.

- At code `8`, the exact variance floor is `994/33>144/5`; at strong code
  `40`, it is `479/13<192/5`. This proves the mean-only incomparability
  asserted after `(15)`. THM-4127 resolves it after slice variance is retained.

- The regular `C5` class, first code `76`, has balanced mean `42` and image
  `{41,43}`. Thus the odd-ceiling `43` in `(10)` is attained, and no uniform
  extra `+2` can be added even for strong parents.

- A global maximizing cut need not be balanced. Strong order-six code `20`
  has cardinality-layer maxima

  ```text
  (19,91,133,131,123,91,19).                              (32)
  ```

  Its unique global maximum `133` is at weight-two mask `48`, while the
  balanced maximum is `131`. Among all `22,320` strong labelled order-six
  parents, exactly `1,920` have no balanced global maximizer.

The exhaustive order-two-through-five audit covers

```text
1,098 labelled parents,
33,864 literal ears,
6,502 cardinality layers,
124,468 old-vertex orderings.                              (33)
```

All `544` strong labelled order-five parents have balanced-layer lattice `2`.
At order six, the strong labelled balanced-lattice histogram is

```text
d=2: 22,080,                 d=6: 240.                    (34)
```

These are finite-exact boundaries, not all-order lattice classifications.

## 6. Scope and exact replay

The source of `(2)` is the full ear fibre, the operation is restriction to
the Johnson slice, and the preserved data are `H,W,F_1`, response parity,
and strongness of nonconstant ears. Averaging in this theorem destroys the
arrangement of the field and edge weights, slice variance, response
connectivity, and the location of a global maximum. THM-4127 later restores
the exact two-energy slice variance and central support dominance, but not an
interval or a balanced global maximizer. The layer lattice restores one
arithmetic coordinate but not interval connectivity; THM-4118's state-labelled
unit components remain the relevant sidecar for propagation.

Run

```text
python3 -B 04-computation/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123.py
python3 -B -O 04-computation/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123.py
python3 -B 04-computation/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123_independent_audit.py
python3 -B -O 04-computation/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_balanced_cardinality_ear_average_layer_lattice_thm4123_independent_audit.py
```

All six streams match their frozen outputs byte-for-byte. The theorem forces
one large response in each central layer and an exact response coset on every
layer. It proves no solid interval, overlap theorem, complete spectrum, or
balanced location for the global maximum. Those remain **OPEN**. **QED.**
