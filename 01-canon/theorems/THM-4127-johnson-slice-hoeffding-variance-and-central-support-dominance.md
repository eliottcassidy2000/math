---
id: THM-4127
title: "Johnson-slice Hoeffding variance and central support dominance"
status: >
  PROVED ELEMENTARY JOHNSON-HOEFFDING DECOMPOSITION + SHARP SLICE SUPPORT +
  CENTRAL DOMINANCE + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every
  cardinality-m tournament-ear fibre for n>=4 has an exact orthogonal
  degree-one and harmonic degree-two decomposition and variance. Its lower
  support yields a Johnson-slice maximum floor with exact layer-coset
  rounding. One middle layer in even order, or the better of both middle
  layers in odd order,
  always matches or improves THM-4115's full-cube support floor; the surplus
  is an explicit nonnegative disjoint-edge energy plus, in odd order, an
  absolute field/star correlation. This resolves THM-4123's former
  incomparability at the support-floor level, not at the raw-variance level.
  No slice interval or balanced global maximizer is proved.
source: codex-frontier-synthesis-creative-20260825h
depends_on:
  - THM-001-redei
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4115-uniform-ear-cut-walsh-variance-and-sharp-growth-refinement
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
related:
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - THM-4118-ear-response-lattice-and-stateful-unit-component-intervals
  - HYP-9029-strong-interval-tiling-law
script: 04-computation/tournament_johnson_slice_hoeffding_variance_thm4127.py
output: 05-knowledge/results/tournament_johnson_slice_hoeffding_variance_thm4127.out
independent_audit_script: 04-computation/tournament_johnson_slice_hoeffding_variance_thm4127_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_johnson_slice_hoeffding_variance_thm4127_independent_audit.out
script_sha256: 44be2322e9bf7a9687e7ba2294d026b3f60ebe65473d8342cf42b0a72ff9b0d2
output_sha256: 677670e4f549063ba8b1ecd0f1ffa1ce78fae7b6c2a7fe91507899885f85b907
semantic_sha256: 5a6ab550db55c636209cba0281d91dd9c2c6161cb5c532f14d78733754fdfd67
independent_audit_script_sha256: a9bf16d168486021161b68832f4f0543ece09275d7c5fe3d15bf9ed2211be0f7
independent_audit_output_sha256: 30bad1bf189eff8a66174250b205c5f973d5d85ca21cda11ede75ccedc6496cc
independent_semantic_sha256: 5a6ab550db55c636209cba0281d91dd9c2c6161cb5c532f14d78733754fdfd67
hash_basis: raw LF bytes
primary_audit: >
  PASS. A Start/End/exposed-gap implementation scans all 33,866 labelled
  tournaments through order six and exactly profiles all 1,098 parents
  through order five plus all 22,320 strong order-six parents. On every
  profile it checks the pointwise decomposition, orthogonal and raw variance
  formulas, support equality case, lattice rounding, central surplus, and
  total-variance identity. Literal child-DP replays and named boundary
  hostiles agree. Normal, optimized, and frozen outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room forward/reverse Hamilton-path implementation imports
  no primary code, reconstructs the cut field from ordered exposed gaps, and
  independently repeats the same census using raw slice incidences. It also
  checks all eight order-three parents, a one-sided odd-layer hostile, every
  preceding strong tie before the first order-seven variance split, and
  literal children at C3 and two order-six controls. Normal, optimized, and
  frozen outputs byte-match; no theorem failure occurs.
---

# THM-4127 -- Johnson-slice Hoeffding variance and central support dominance

**PROVED ELEMENTARY JOHNSON-HOEFFDING DECOMPOSITION + SHARP SLICE
SUPPORT + CENTRAL DOMINANCE + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4115 retained the complete Boolean-cube variance, while THM-4123 retained
the exact mean and lattice on each cardinality layer. Their floors were
previously incomparable because the conditioning step discarded variance.
The missing coordinate is an exact two-energy Hoeffding decomposition on the
Johnson graph. Restoring it makes a central slice dominate both earlier
floors without asserting that a global maximizing cut is balanced.

## 1. The cut field and its edge ANOVA

Let `T` be a tournament on `V={1,...,n}`. In THM-4115's notation, its ear
response is

```text
F(S)=H+sum_(i in S) h_i+cut_w(S),
sum_i h_i=0,                  w_ij>=0,
W=sum_(i<j)w_ij>0.                                    (1)
```

Here `H=H(T)`, and `W=((n-1)H+F_1(T))/2`. For `n>=3`, put

```text
d_i=sum_(j!=i)w_ij,
u_i=(d_i-2W/n)/(n-2),
wbar=2W/(n(n-1)),
z_ij=w_ij-wbar-u_i-u_j.                                (2)
```

For edge arrays, write

```text
||w||^2=sum_(i<j)w_ij^2,       ||z||^2=sum_(i<j)z_ij^2. (3)
```

> **Lemma 1 (edge ANOVA).** The coordinates in `(2)` satisfy
>
> ```text
> sum_i u_i=0,                 sum_(j!=i)z_ij=0 for every i,
> ||w||^2=2W^2/(n(n-1))+(n-2)||u||^2+||z||^2.           (4)
> ```

Indeed `sum_i d_i=2W`, which gives the first identity. For a fixed `i`, the
row sum of `z` is

```text
d_i-(n-1)wbar-(n-2)u_i=0.
```

The constant edge, star edge `u_i+u_j`, and zero-row-sum residual are
pairwise orthogonal. Also

```text
sum_(i<j)(u_i+u_j)^2=(n-2)||u||^2,
```

which proves the energy identity. **QED.**

## 2. Exact pointwise decomposition on a Johnson slice

Fix `0<=m<=n`, put `q=n-m`, and let `S` be uniform on the Johnson slice
`J(n,m)`. Define

```text
alpha=mq/(n(n-1)),
mu_m=H+2alpha W,
g_m=h+(n-2m)u.                                          (5)
```

> **Theorem 1 (pointwise Johnson-Hoeffding decomposition).** For every
> `m`-subset `S`,
>
> ```text
> F(S)=mu_m+sum_(i in S)g_m(i)
>             -2sum_({i,j} subseteq S)z_ij.              (6)
> ```

The constant edge `wbar` crosses an `m:q` cut exactly `mq` times. The star
part crosses with total

```text
q sum_(i in S)u_i+m sum_(i notin S)u_i
=(n-2m)sum_(i in S)u_i.                                  (7)
```

Finally, the zero row sums in `(4)` give

```text
sum_(i in S,j notin S)z_ij
=-2sum_({i,j} subseteq S)z_ij.                           (8)
```

Substituting `(7)--(8)` into `(1)` proves `(6)`. Both nonconstant summands
have slice mean zero, so `(5)` also recovers THM-4123's exact mean. **QED.**

## 3. Exact slice variance

For `n>=4`, define

```text
c_2=m(m-1)q(q-1)/(n(n-1)(n-2)(n-3)).                    (9)
```

> **Theorem 2 (two-energy variance).** Uniformly on `J(n,m)`,
>
> ```text
> Var_m F=alpha ||g_m||^2+4c_2||z||^2.                  (10)
> ```
>
> The degree-one and degree-two terms in `(6)` are orthogonal.

For a zero-sum vector `g`, direct sampling without replacement gives

```text
Var(sum_(i in S)g_i)=alpha||g||^2.                       (11)
```

Put `Z(S)=sum_({i,j} subseteq S)z_ij`. If `A` and `D` denote the sums of
`z_e z_f` over unordered incident and disjoint edge pairs, respectively,
the row sums and total sum of `z` give

```text
A=-||z||^2,                     D=||z||^2/2.             (12)
```

Writing `p_k=(m)_k/(n)_k`, expansion by equal, incident, and disjoint edge
pairs yields

```text
E[Z^2]=(p_2-2p_3+p_4)||z||^2=c_2||z||^2.                (13)
```

The covariance with `(11)` vanishes because it is a multiple of
`sum_j z_ij` in every row. The quadratic term in `(6)` is `-2Z`, so its
variance is **`4c_2||z||^2`**. This factor four is essential. **QED.**

There is also a useful raw form. Put

```text
H_2=sum_i h_i^2,       D_2=sum_i d_i^2,
Q_2=sum_(i<j)w_ij^2,  C_hd=sum_i h_i d_i,
delta=4c_2.                                                (14)
```

Substitution of `(2)--(4)` into `(10)` gives

```text
Var_m F=alpha H_2
 +2alpha(n-2m)C_hd/(n-2)
 +delta Q_2+(alpha-delta)D_2
 +(delta-4alpha^2)W^2.                                  (15)
```

Thus either the orthogonal packet or the raw aggregate packet determines the
complete slice variance.

## 4. Slice support and exact coset rounding

Let

```text
M_m=max_(|S|=m)F(S),                1<=m<=n-1.           (16)
```

THM-4115 proves `F(S)>=H` for every ear. Therefore

```text
(F(S)-H)(M_m-F(S))>=0.                                  (17)
```

Averaging `(17)` and using `mu_m-H=2alpha W>0` proves

```text
M_m>=J_m:=mu_m+Var_m F/(mu_m-H).                         (18)
```

For `n>=4`, Theorem 2 simplifies this to

```text
J_m=H+2alpha W+||g_m||^2/(2W)
 +2(m-1)(q-1)||z||^2/((n-2)(n-3)W).                     (19)
```

Equality in `(18)` holds exactly when every response on that slice belongs
to `{H,M_m}`; this includes the constant-slice case.

THM-4123 proves that a layer lies in one exact coset

```text
a_(T,m)+d_(T,m) Z.                                       (20)
```

Consequently `(18)` sharpens to

```text
M_m>=ceil_(a_(T,m)+d_(T,m)Z)(J_m).                       (21)
```

When `d_(T,m)=0`, the layer is constant and `J_m=a_(T,m)`. Since every
response is odd by THM-001, ordinary odd-ceiling rounding is always a weaker
fallback. Formula `(21)` is a maximum floor, not a claim that the layer image
is an interval.

## 5. Central slices dominate the full cube

Put `b=floor(n^2/4)` and denote THM-4123's common central mean by

```text
B_n=H+2bW/(n(n-1)).                                      (22)
```

For even `n>=4`, `m=n/2` makes `g_m=h`; the edge-star coordinate disappears:

```text
Var_(n/2)F=n||h||^2/(4(n-1))
            +n(n-2)||z||^2/(4(n-1)(n-3)),
J_bal=B_n+||h||^2/(2W)+(n-2)||z||^2/(2(n-3)W).           (23)
```

For odd `n>=5`, the lower and upper central layers have equal mean but
generally different variance. Their floors are

```text
J_down=B_n+||h+u||^2/(2W)+(n-1)||z||^2/(2(n-2)W),
J_up  =B_n+||h-u||^2/(2W)+(n-1)||z||^2/(2(n-2)W).        (24)
```

Compare these with THM-4115's exact full-cube support floor

```text
J_box=H+W/2+(||h||^2+||w||^2)/(2W).                     (25)
```

Define the nonnegative disjoint-edge energy

```text
D_4=sum_(e<f, e intersect f=empty) w_e w_f.              (26)
```

Expanding `W^2` and the degree squares gives

```text
2D_4=W^2+||w||^2-sum_i d_i^2.                           (27)
```

Straight substitution of `(4)` into `(23)--(25)` now yields

```text
J_bal-J_box=D_4/(W(n-3))>=0                    (n even), (28)

J_down-J_box=D_4/(W(n-2))+<h,u>/W,
J_up-J_box  =D_4/(W(n-2))-<h,u>/W              (n odd).  (29)
```

Hence

```text
max(J_down,J_up)-J_box
=D_4/(W(n-2))+|<h,u>|/W>=0.                              (30)
```

> **Corollary 3 (central support dominance).** One middle layer in even
> order, or the better of both middle layers in odd order, always gives a
> support floor at least as large as the full-cube floor. It also dominates
> its own mean floor. After applying `(21)`, it remains at least as strong as
> ordinary odd rounding of `(25)`.

The word **better** is necessary. For order-five code `11`,

```text
H=3, W=21, D_4=60, C_hd=-66,
J_down=421/21 < J_box=141/7 < J_up=155/7.                (31)
```

Thus a predetermined odd central layer has no unconditional dominance sign.
A selection-robust bank must expand both middle layers.

If `B_n` is a selected bank of strong parents and every central cut is
expanded before arbitrary value-representative selection, THM-4097 preserves
strongness and `(21)` gives the parent-specific recurrence

```text
M_(n+1)>=max_(m in {floor(n/2),ceil(n/2)})
            ceil_(a_(T,m)+d_(T,m)Z)(J_m(T))              (32)
```

for a maximizing retained parent `T`. This resolves THM-4123's former
comparison with `(25)` at the level actually used for growth.

## 6. Raw variance, boundaries, and information budget

Central support dominance is not raw-variance dominance. If `V_box` denotes
the full-cube variance, conditioning on `|S|` gives the exact identity

```text
V_box=sum_(m=0)^n binom(n,m)/2^n Var_m F
       +W^2/(2n(n-1)).                                   (33)
```

For `C_3`, both central variances vanish while `V_box=3/4`; nevertheless all
three support floors equal five. This is why the support denominator in
`(18)`, not a raw variance comparison, is the faithful consequence.

The boundaries are explicit.

- At `m=0,n`, the response is `H`, the variance is zero, and the quotient in
  `(18)` is undefined.
- At `m=1,n-1`, the harmonic degree-two residual is invisible.
- At `n=3`, zero row sums force `z=0`; only the degree-one energy remains.
- At `n=2,m=1`, `Var_m F=||h||^2/2` and `J_m=J_box`.
- At `n=1`, the response is the literal constant pair `(1,1)` and `W=0`.
- For a non-ear quadratic with `W=0`, the variance identities remain valid,
  but the support quotients are not licensed.

A sufficient invariant packet for the rational variance and support floor on
one slice is

```text
(H,W,||g_m||^2,||z||^2).                                 (34)
```

For every rational slice floor simultaneously it is

```text
(H,W,||h||^2,||u||^2,<h,u>,||z||^2).                    (35)
```

Exact lattice rounding additionally needs THM-4123's layer coset
`(a_(T,m),d_(T,m))`.

The smaller THM-4123 packet `(H,F_1,W)` is insufficient: order-five codes
`8` and `759` have the same balanced mean `144/5`, but variances `1469/25`
and `1529/25`. Equal odd central means also need not have equal variances;
the first strong labelled-code split occurs at order-seven code `34`:

```text
mu_3=mu_4=1209/7,
Var_3=234944/49,                 Var_4=1205856/245.       (36)
```

The sole preceding strong code ties, so “first” in `(36)` is an audited code
ordering statement.

## 7. Exact controls and census

The named exact controls are:

| parent | central mean | slice variance | `J_m` | odd / coset floor | `J_box` -> odd | central / global max |
|---|---:|---:|---:|---:|---:|---:|
| `C_3` | `5` | `0` | `5` | `5 / 5` | `5 -> 5` | `5 / 5` |
| code `8` | `144/5` | `1469/25` | `3145/99` | `33 / 33` | `994/33 -> 31` | `41 / 41` |
| code `40` | `192/5` | `261/25` | `505/13` | `39 / 39` | `479/13 -> 37` | `43 / 43` |
| regular `C_5`, code `76` | `42` | `1` | `1135/27` | `43 / 43` | `358/9 -> 41` | `43 / 43` |
| code `759` | `144/5` | `1529/25` | `287/9` | `33 / 33` | `333/11 -> 31` | `43 / 43` |
| code `20` | `422/5` | `17991/25` | `10399/109` | `97 / 97` | `9636/109 -> 89` | `131 / 133` |
| code `30571` | `738/5` | `4401/25` | `2837/19` | `151 / 153` | `2619/19 -> 139` | `189 / 189` |

Code `30571` shows the extra strength of its layer lattice `d=6`. Code `76`
attains the rounded floor `43`, so no universal extra `+2` survives. Code
`20` retains THM-4123's hostile: a global maximum need not be balanced.

The strong-parent census is:

| order | strong parents | central vs cube, rational | central vs cube, rounded | improves / ties old mean, odd | improves / ties old mean, coset | attains coset floor |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | `0 / 2` | `0 / 2` | `0 / 2` | `0 / 2` | 2 |
| 4 | 24 | `24 / 0` | `0 / 24` | `0 / 24` | `0 / 24` | 0 |
| 5 | 544 | `544 / 0` | `544 / 0` | `480 / 64` | `480 / 64` | 24 |
| 6 | 22,320 | `22,320 / 0` | `22,320 / 0` | `22,080 / 240` | `21,840 / 480` | 0 |

Here every pair is `strict / equal`. Both implementations also check
central/cube equality for all two order-two and all eight order-three labelled
tournaments, and strict rational improvement for all 64 order-four and all
1,024 order-five labelled tournaments.

The primary implementation obtains responses from Start/End/exposed-gap
tables. The clean-room audit instead builds independent forward and converse
Hamilton-path tables, reconstructs the cut field from raw ordered gaps, and
then uses fixed-slice incidences. Literal child Hamiltonian DPs agree at C3,
code `20`, and code `30571`. Run

```text
python3 -B 04-computation/tournament_johnson_slice_hoeffding_variance_thm4127.py
python3 -B -O 04-computation/tournament_johnson_slice_hoeffding_variance_thm4127.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_johnson_slice_hoeffding_variance_thm4127.py
python3 -B 04-computation/tournament_johnson_slice_hoeffding_variance_thm4127_independent_audit.py
python3 -B -O 04-computation/tournament_johnson_slice_hoeffding_variance_thm4127_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_johnson_slice_hoeffding_variance_thm4127_independent_audit.py
```

All six streams match their frozen outputs byte-for-byte. The theorem forces
a large central child and an exact layer-coset floor. It does not prove a
slice interval, a balanced global maximizer, overlap of recursively selected
images, or the global `H`-spectrum conjecture. **QED.**
