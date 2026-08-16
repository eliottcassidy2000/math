---
id: THM-3503
title: "Rule 30 odometer ultrametric regrading and orbit-closure dimensions"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The seed-orbit odometer
  conjugacy is shown to regrade
  every dyadic phase distance exactly by the next innovation depth.  This
  identifies the Hausdorff, packing, and box dimensions of the orbit closure
  with lower and upper innovation densities.  No Rule 30 prize is claimed.
source: root/rule30-sharp-unlocks/odometer-fractal-scale/2026-08-16
audit: >
  PASS (2026-08-16), independent proof, scope, and adversarial replay audit.
  The auditor rederived the all-Z_2 metric regrading, uniform ball masses,
  Hausdorff/packing/box dimension identities, innovation indexing, no-111
  bound, and wrap implication.  The repaired companion exhausts all 32,640
  unordered phase pairs in the width-24 quotient, uses integer-only scalar
  gates, and matches under ordinary and optimized execution.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3493-rule30-dyadic-wrap-atlas
related:
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3500-rule30-dyadic-section-cut-defect-and-cross-depth-valuation-carrier
script: 04-computation/rule30_odometer_ultrametric_dimension_thm3503.py
output: 05-knowledge/results/rule30_odometer_ultrametric_dimension_thm3503.out
script_sha256: e82e732c819ad2c0c04a5f206d6ff73cc05b06ce13dc440fb15a4fcfc259bdc5
output_sha256: 02fcc00d9ea2240f56235213f1ed5b9af1af5db8df8fa5ccf32542cb4ac19dbc
hash_basis: raw bytes
---

# THM-3503 -- Rule 30 odometer ultrametric regrading and orbit-closure dimensions

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3458 identifies the closure of the packed one-seed edge orbit with the
2-adic odometer, but treats that identification only topologically.  The
metric inherited from the ambient state space contains more information.  It
is not the ordinary phase metric: every phase scale is moved to the depth of
the corresponding Rule 30 innovation.

This gives an exact geometric meaning to the valuation sequence isolated in
THM-3493.  The sequence is the radial distortion function of the odometer
embedding, and its reciprocal slopes are precisely the fractal dimensions of
the orbit closure.  In particular, a positive Hausdorff-dimension bound would
force all sufficiently large depths into THM-3493's hard regime.  Proving such
a bound remains open.

## 1. Conventions and inheritance

Use the packed Rule 30 map and seed orbit

```text
Phi(x)=x xor ((2x) or (4x)),
R_t=Phi^t(1).                                           (1)
```

Let `P_w` be the least return time of the seed modulo `2^w`.  THM-3458 and
THM-3463 give

```text
P_(w+1)=2^(epsilon_w)P_w,
epsilon_w in {0,1},                                    (2)
```

and the periods are unbounded powers of two.  Enumerate the innovation
depths

```text
kappa_1<kappa_2<...,
epsilon_(kappa_r)=1.                                   (3)
```

For `m>=0`, retain THM-3493's valuation notation

```text
v_m=kappa_(m+1)=nu_2(R_(2^m)-1).                       (4)
```

Thus

```text
P_(v_m)=2^m,       P_(v_m+1)=2^(m+1).                 (5)
```

Put

```text
E_w=log_2 P_w=#{r:kappa_r<w}.                          (6)
```

The inheritance pass is:

1. the closest proved mechanism is THM-3458's odometer homeomorphism;
2. the canonical hostile is THM-3493's scalar period tower with very sparse
   innovations, which obeys all current period floors but has zero innovation
   density;
3. the corrected near miss is that topological conjugacy alone says nothing
   quantitative about the ambient 2-adic metric; and
4. the least-used sidecar is the complete innovation address `v_m`, rather
   than only the period value `2^m`.

The live objects are phase cylinders, state cylinders, innovation depths,
the pushforward Haar measure, covering numbers, and hard dyadic blocks.

## 2. Exact ultrametric regrading

Let

```text
X=closure{R_t:t>=0} subset Z_2.                        (7)
```

THM-3458 supplies the unique homeomorphism

```text
iota:Z_2 -> X,
iota(t+1)=Phi(iota(t)),       iota(0)=1.               (8)
```

Both copies of `Z_2` carry the standard valuation `nu_2` and metric
`|x-y|_2=2^(-nu_2(x-y))`.

### Theorem 2.1 (exact radial distortion)

For distinct `s,t in Z_2`, put `m=nu_2(t-s)`.  Then

```text
boxed:
nu_2(iota(t)-iota(s))=v_m.                             (9)
```

Equivalently,

```text
|iota(t)-iota(s)|_2=2^(-v_m)
whenever |t-s|_2=2^(-m).                              (10)
```

Thus `iota` is an isometry after regrading the phase radius `2^(-m)` to
`2^(-v_m)`.  This is stronger than a pointwise statement at the seed.

### Proof

Fix `m` and put `k=v_m`, `P=2^m`.  By (5), every state coordinate below `k`
has period dividing `P`.  The coordinate `b_k(u)=bit_k(R_u)` receives the
same lower-coordinate cocycle on the two time intervals starting at `u` and
`u+P`.  Hence

```text
D(u)=b_k(u+P)+b_k(u)                                   (11)
```

is independent of `u`.  At `u=0`, equations (2)--(5) give

```text
D(0)=epsilon_k=1.                                     (12)
```

So translation by `P` fixes every bit below `k` and complements bit `k`, at
every phase.  Translation by an odd multiple of `P` has the same first
effect.  If `nu_2(t-s)=m`, then `t-s` is such an odd multiple in `Z_2`.
Integer approximation and continuity of `iota` prove (9).  Equation (10) is
the definition of the 2-adic metric.  QED.

### Corollary 2.2 (uniform cylinder geometry)

Every phase coset modulo `2^m` maps to a subset of `X` of diameter exactly
`2^(-v_m)`, and its two child cosets are separated at that same state scale.
There is no exceptional phase with better or worse distortion.

This uniformity is the metric form of THM-3471's innovation-cube atlas.  It
does **not** make the nonlinear innovation coordinates a group homomorphism.

## 3. Exact covering and mass laws

At state precision `w`, the image of `X` modulo `2^w` is the embedded seed
cycle.  Therefore the least number `N_X(w)` of residue balls modulo `2^w`
covering `X` is

```text
boxed:
N_X(w)=P_w=2^(E_w).                                   (13)
```

Let `lambda` be Haar probability on the phase `Z_2` and put

```text
mu=iota_*lambda.                                      (14)
```

Every nonempty state ball of precision `w` has one phase residue modulo
`P_w`, hence

```text
boxed:
mu(B)=1/P_w=2^(-E_w)                                  (15)
```

for every precision-`w` ball `B` meeting `X`.

Equations (13)--(15) are exact at every scale, not asymptotic equidistribution
claims about the moving center observer.

## 4. Hausdorff, packing, and box dimensions

All dimensions below use the ambient 2-adic metric on `Z_2`.

### Theorem 4.1 (dimension dictionary)

The orbit closure has

```text
dim_H X = lower_box_dim X = liminf_(w->infinity) E_w/w,
dim_P X = upper_box_dim X = limsup_(w->infinity) E_w/w. (16)
```

In terms of the next-innovation valuations,

```text
boxed:
dim_H X = liminf_(m->infinity) m/v_m,
dim_P X = limsup_(m->infinity) m/v_m.                 (17)
```

### Proof

The box formulas follow immediately from (13), because a precision-`w`
2-adic ball has diameter `2^(-w)`.

Put `d_-=liminf E_w/w`.  Precision levels along a subsequence approaching
`d_-` give covers of total `s`-content

```text
P_w 2^(-sw)=2^(E_w-sw) -> 0                           (18)
```

for every `s>d_-`; hence `dim_H X<=d_-`.  Conversely, if `s<d_-`, equation
(15) gives, for all sufficiently small balls meeting `X`,

```text
mu(B)<=diam(B)^s.                                     (19)
```

The mass-distribution argument gives `dim_H X>=s`.  Letting `s` increase to
`d_-` proves the first equality in (16).

Every point has the same upper local dimension under `mu`, namely
`limsup E_w/w`.  The packing mass inequality gives the corresponding lower
bound for `dim_P X`, while `dim_P X<=upper_box_dim X` gives the upper bound.
Equivalently, one may apply the Baire definition of packing dimension: every
nonempty relative cylinder has the same tail covering exponent, so no
countable cover can lower the maximum upper-box exponent.  This proves the
second equality in (16).

Finally, `E_(v_m)=m`; between consecutive innovations `E_w` is constant.
The smallest ratios occur at the right ends of these plateaux and the largest
immediately after an innovation.  Thus

```text
liminf E_w/w=liminf m/v_m,
limsup E_w/w=limsup (m+1)/(v_m+1)=limsup m/v_m,        (20)
```

which proves (17).  QED.

## 5. Consequences and sharp stopping boundaries

### 5.1 The proved upper dimension bound

THM-3463 proves that the lift word `epsilon_1 epsilon_2 ...` contains no
`111`.  Consequently

```text
E_w <= (2/3)w+O(1),
boxed: dim_P X<=2/3.                                  (21)
```

The scalar word `(110)^infinity` shows that the constant `2/3` cannot be
improved from the no-`111` constraint alone.

### 5.2 Positive dimension would close all late wraps

By THM-3493, the dyadic block `[2^m,2^(m+1)-1]` contains a wrapped depth only
if

```text
v_m>=2^m.                                             (22)
```

If this happens for infinitely many `m`, then (17) forces

```text
dim_H X=0.                                            (23)
```

Equivalently,

```text
boxed:
dim_H X>0  ==>  only finitely many Rule-30 depths are wrapped. (24)
```

Indeed `dim_H X>0` is equivalent to the uniform linear bound `v_m=O(m)`,
which is much stronger than `v_m<2^m`.  The converse to (24) is false at the
level of scalar period towers: subexponential but superlinear innovation gaps
can have no late wraps and still give dimension zero.

The Mersenne innovation model

```text
kappa_r=2^r-1                                         (25)
```

obeys the dyadic period floor and no-`111` rule but has dimension zero.  It is
the canonical hostile to extracting positive dimension from the current
scalar constraints.

### 5.3 What this does not prove

The result changes the target but does not solve it.  It proves neither a
positive lower dimension nor any bound on the gaps `v_(m+1)-v_m`.  Even
positive dimension would only place almost all depths in the hard pointed-arc
regime; it would not determine the marked center bits there.  No temporal
balance, nonperiodicity, or random-access complexity conclusion follows.

In particular:

- Haar-typical phase geometry does not control the distinguished phase zero;
- fractal dimension of `X` does not turn the moving center into a fixed
  continuous observable;
- an upper entropy bound for a universal staircase carrier is not a lower
  dimension bound for this one orbit; and
- the finite innovation prefix is evidence, not an asymptotic dimension
  estimate.

## 6. Information ledger and next decisive tests

| map | preserves | destroys | required sidecar |
|---|---|---|---|
| phase `t` -> state `iota(t)` | complete orbit and addition dynamics | ordinary phase scale | innovation depth `v_m` |
| state -> precision-`w` residue | exact covering count `P_w` | later branching | period/innovation tower |
| Haar -> `mu` | uniform cylinder masses | marked phase owner | calibrated origin |
| `v_m` -> dimensions | lower/upper asymptotic innovation densities | finite heads and individual center bits | full gap sequence |
| phase atlas -> center observer | terminal-arc semantics if retained | observer depth under fixed truncation | moving address `k` |

The next sharp targets are therefore:

1. prove any uniform upper bound `v_m<=Cm`, equivalently `dim_H X>0`;
2. derive a renormalization inequality for the dyadic section-cut defect of
   THM-3500 that prevents arbitrarily long innovation gaps;
3. couple that defect to the calibrated marked terminal arc, since dimension
   alone remains point-blind; and
4. test whether a finite correlation state can control the OR-intersection
   current without silently replacing the distinguished orbit by Haar phase.

These are concrete stopping tests.  Another recurrence fit to the finite
`v_m` prefix is not one.

## 7. Exact companion

Run

```bash
python3 04-computation/rule30_odometer_ultrametric_dimension_thm3503.py
python3 -O 04-computation/rule30_odometer_ultrametric_dimension_thm3503.py
```

The companion:

- constructs every seed quotient through width `24`;
- checks the period-lift and innovation indexing exactly;
- exhausts all `32,640` unordered phase pairs in the width-24 seed quotient,
  including every available odd multiplier, for the metric law;
- verifies `2665` exact covering-state cells;
- freezes the innovation prefix `(1,3,4,6,7,9,15,16,24)`; and
- checks finite Mersenne-sparse and density-`2/3` scalar boundary models.

These finite checks are controls for the universal proofs, not extrapolations
of (16)--(17).  The independent audit rederived the metric and dimension
formulas, checked scope against all three Rule 30 prizes, and replayed the
ordinary and optimized companions byte-for-byte against the pinned output.
