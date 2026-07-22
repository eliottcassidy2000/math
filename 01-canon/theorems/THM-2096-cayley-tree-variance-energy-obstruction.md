---
id: THM-2096
title: "Cayley-tree variance strengthens the centered triple-energy obstruction"
status: >
  PROVED. For the 21 nonnegative relative-Hunter edge weights on K_7, the
  uniform spanning-tree weight has variance
  (12 sum w_e^2-sum s_i^2)/49. Its maximum is at least its mean plus variance
  divided by the mean. Consequently THM-2091's necessary centered-energy
  threshold gains the explicit dispersion term (7/2)Var/mu. The refinement
  closes THM-2091's hostile negative-average control and raises the bounded
  bank's universal-threshold closures from 29 to 69. The complete maximum
  spanning tree can still be stronger; this is not LRC(14).
source: codex-2026-07-22-LRC-Cayley-tree-variance
depends_on:
  - THM-2081
  - THM-2091
related:
  - THM-1122
  - THM-1166
  - THM-2094
script: 04-computation/lrc14_cayley_tree_variance_codex_20260722.py
output: 05-knowledge/results/lrc14_cayley_tree_variance_codex_20260722.out
script_sha256: ceae5b5af684c60dbcdbc714d79755048e0ed6c130fb4fce984f82096044684d
output_sha256: b206120ac59964f916a7262ffab9b5a56280a8823ccc77e23f48c52fda62f5ec
hash_basis: repository blobs with LF line endings
---

# THM-2096 -- dispersion in the Cayley tree pays relative-Hunter margin

Let `w_e>=0`, `e in E(K_7)`, be arbitrary edge weights. Put

```text
W=sum_e w_e,
s_i=sum_(e incident to i) w_e,                         (1)
```

and let `T` be a uniformly random labelled spanning tree of `K_7`. Define

```text
X_T=sum_(e in T)w_e,
mu=E X_T,              V=Var(X_T),
tau=max_T X_T.                                           (2)
```

Then

```text
mu=(2/7)W,
V=[12 sum_e w_e^2-sum_i s_i^2]/49.                    (3)
```

If `mu>0`, then

```text
tau>=mu+V/mu.                                          (4)
```

## 1. Exact Cayley incidences

There are `7^5=16807` labelled trees. The standard forest-extension form of
Cayley's formula says that a fixed forest whose component sizes are
`n_1,...,n_c` lies in

```text
7^(c-2) product_j n_j                                  (5)
```

trees. Therefore a fixed edge occurs in

```text
2*7^4=4802                                             (6)
```

trees, while two distinct fixed edges occur together in

```text
3*7^3=1029       if they share a vertex,
4*7^3=1372       if they are disjoint.                 (7)
```

Equivalently,

```text
P(e in T)=2/7,
P(e,f in T)=3/49  for adjacent e,f,
P(e,f in T)=4/49  for disjoint e,f.                    (8)
```

The first identity in (3) follows immediately.

## 2. The variance identity

Let

```text
S_2=sum_e w_e^2,
A=sum_({e,f}:e,f adjacent)w_e w_f,
B=sum_({e,f}:e,f disjoint)w_e w_f.                     (9)
```

Equations (8) give

```text
E X_T^2=(14S_2+6A+8B)/49,
mu^2=(4S_2+8A+8B)/49.                                 (10)
```

Also

```text
sum_i s_i^2=2S_2+2A.                                  (11)
```

Subtracting the two lines of (10) and using (11) proves the second identity
in (3).

Since `0<=X_T<=tau`, pointwise

```text
X_T^2<=tau X_T.                                       (12)
```

Taking expectations and dividing by `mu>0` gives

```text
tau>=E X_T^2/E X_T=mu+V/mu,                           (13)
```

which is (4). QED.

## 3. Refined relative-Hunter gate

Apply the abstract result to THM-2091's restricted weights

```text
w_ij=measure(D_(q_i) intersect D_(q_j) intersect E_h^c).
```

Write

```text
Delta=2/7-sum_i measure(D_(q_i) intersect E_h),
A=mu-Delta.                                            (14)
```

THM-2081 and (4) give

```text
measure(G_Q minus E_h)
 >=tau-Delta
 >=A+V/mu.                                             (15)
```

Here the first inequality means that a positive right side supplies that much
Hunter safe mass; a negative lower bound is simply vacuous. Thus guard
containment, when `mu>0`, forces

```text
A<=-V/mu.                                              (16)
```

THM-2091 proved

```text
A>=2059/315315-(2/7)D,                                 (17)
```

where `D` is its centered guard/danger triple energy. Combining (16)--(17)
gives the sharper necessary condition

```text
D>=2059/90090+(7/2)V/mu.                               (18)
```

When `V=0`, this is exactly THM-2091. Positive dispersion measures how much
edge heterogeneity uniform averaging discarded. If `mu=0`, equation (18) is
not asserted and THM-2091's unrefined threshold remains available.

The modulo-seven charge caps in THM-2091 sharpen the first term in (17) in
exactly the same way: for `k=1,2,3,4`, replace `2059/90090` in (18) by the
corresponding row of THM-2091 (20c), then add `(7/2)V/mu`.

## 4. Hostile control and bounded bank

THM-2091's control

```text
h=1,            Q=(4,5,6,7,11,13,27)                  (19)
```

has exact data

```text
mu=57481/1203930,
Delta=181007/3783780,
V=919435861/5314640631300.                             (20)
```

Its mean margin is negative,

```text
A=-2467/26486460<0,                                    (21)
```

but the variance-refined margin is

```text
A+V/mu=488619049/138406200660>0.                       (22)
```

Thus (15) closes the precise row used to show that uniform averaging alone is
not equivalent to the maximum tree. The exact maximum-tree margin remains

```text
86717/1891890,                                         (23)
```

so one second moment does not reconstruct the optimum.

On the complete THM-2081 terminal bank through maximum speed `24`, all `1322`
scalar survivors have `V>0` and positive variance-refined margin. Applying
only the universal energy threshold closes `29` rows; applying (18) closes
`69`. This is exact bounded telemetry, not the all-height proof.

## 5. Higher forest moments and transfer scope

For every integer `r>=0` with `E X_T^r>0`, the same pointwise argument gives

```text
tau>=E X_T^(r+1)/E X_T^r.                              (24)
```

Because the tree set is finite, these ratios converge to `tau`. Their exact
expansion is indexed by unions of `r+1` forests and can be evaluated through
Cayley forest incidences. Equation (4) is the first nontrivial rung of a
systematic moment hierarchy between uniform averaging and the maximum
graphic-matroid basis.

This is the precise transfer from THM-1122/1166's multiplicity moments: retain
moments of the **tree-basis weight**, not only moments of pointwise danger
multiplicity. The carrier preserves edge location and therefore repairs the
pairwise-only LP loss documented in THM-1122.

The mechanism does not prove Gaussian Moment noncancellation. LRC's `w_e` and
`X_T` are nonnegative measures, which makes (12) possible. General GMC
channels and the unit-distance Bessel kernel are signed or complex; applying
absolute values there would erase the cancellation predicate. The legitimate
cross-domain lesson is narrower: once positivity exists, dispersion can turn
a failed average into an extremal certificate.

Candidate tournament vertices were runners, edges, spanning trees, and moment
obligations. Orienting edges by weight destroys the adjacent/disjoint pair
incidences in (8). The faithful object is the weighted graphic matroid with its
Cayley basis distribution. Tournament fingerprints add no invariant beyond
an arbitrary tie-break and are not used.

## Exact referee

The companion exhausts all `16807` Prüfer trees, verifies (6)--(13) on exact
weighted controls, replays all `4120` bounded terminal pairs and `1322` scalar
survivors, and checks (19)--(23). Runtime checks remain active under
optimization; normal and `python -O` transcripts byte-match the stored output
and end in `PASS`. QED.
