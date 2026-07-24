---
id: THM-2143
title: "The Gibbs zero-count PGF is a complete certificate for strict loneliness; finite-temperature max certificates lose the tight boundary; and the signed log-Riesz pairing is unsound"
status: >
  PROVED (Gibbs zero-count equivalence, rational/artanh certificate, strict
  finite-temperature obstruction, all-n AP log-Riesz counterexample, and
  unrestricted-frequency unboundedness). VERIFIED-EXACT (small-n controls and
  eight deterministic h=3/41 defect-seven rows). This does not prove weak
  LRC(14): every measure/PGF statement here deliberately separates the
  positive-measure strict branch from isolated tight witnesses.
source: codex-2026-07-24-artanh-riesz-soundness
depends_on:
  - THM-406-covering-depth-master-object-factorial-moments-and-spectral-identity
  - THM-2047-phase-height-toric-arrangement-for-lrc
related:
  - HYP-9023
  - THM-515
script: 04-computation/lrc_artanh_gibbs_soundness_thm2143.py
output: 05-knowledge/results/lrc_artanh_gibbs_soundness_thm2143.out
---

# THM-2143 — the exact soundness boundary for artanh, Gibbs, and log-Riesz certificates

Let `T=R/Z`, let `V` be a nonempty finite set of positive integral speeds,
and fix a rational threshold `0<h<1/2`.  Up to a null set, it is immaterial
whether the danger combs

```text
D_v(h)={t: ||v t||<h}
```

use strict or weak endpoints.  Define the danger depth, its distribution, the
strict lonely measure, and the max-min gap by

```text
H(t)=sum_{v in V} 1_{D_v(h)}(t),
p_k=meas{t:H(t)=k},
p_0=meas{t:H(t)=0},
g_V(t)=min_{v in V} ||v t||,
Gap(V)=max_t g_V(t).
```

THM-406 identifies `(p_0,...,p_|V|)` as the spectral measure of the depth
observable.  THM-2047 explains the coordinate it loses: `p_0=0` does not
distinguish an isolated weak witness `Gap(V)=h` from failure `Gap(V)<h`.

## 1. The Gibbs zero-count map is sound and complete for the strict branch

For `0<z<1`, put

```text
Z_V(z)=integral_T z^{H(t)} dt=sum_{k=0}^{|V|} p_k z^k.       (1)
```

Then

```text
Z_V(z)-z
  =(1-z)p_0-sum_{k=2}^{|V|}(z-z^k)p_k.                       (2)
```

Consequently

```text
p_0>0
  iff there exists z in (0,1) with Z_V(z)>z
  iff there exists rational z in (0,1) with Z_V(z)>z.        (3)
```

**Proof.**  Equation (2) is obtained by separating the `k=0,1` terms in
(1).  If `p_0=0`, then `H>=1` almost everywhere and `z^H<=z` pointwise, so
`Z_V(z)<=z`.  Conversely, as `z` decreases to zero,
`Z_V(z)-z -> p_0`; hence `p_0>0` makes the difference positive on some
interval `(0,z_0)`.  A rational `z` lies in that interval. ∎

Thus a **single** temperature is only sufficient: the negative sum in (2) is
the exact multiple-overlap debt.  The whole zero-temperature family is
complete.  Equivalently,

```text
Z_V(z)=integral_T product_{v in V}
          (1-(1-z)1_{D_v(h)}(t)) dt,                         (4)
```

so this is a nonnegative indicator-Riesz product, or the full overlap PGF
with fugacity `z`.  Unlike `integral H log R`, its implication comes from the
pointwise full-cover inequality `z^H<=z`; no signed measure is introduced.
Expanding the product gives the exact Riesz/inclusion-exclusion form

```text
Z_V(z)=sum_{I subset V}(z-1)^|I| meas(intersection_{v in I}D_v(h)).
```

Thus overlap estimates, or the Fejer-regularized relation-lattice formula of
THM-2047 applied to each intersection, can estimate `Z_V`.  The price is
sharp: as `z` tends to zero this becomes the full alternating formula for
`p_0`, so the zero-temperature completeness does not bypass the all-orders
cancellation wall.

### What it preserves and destroys

- **Source:** the full depth distribution `(p_0,...,p_|V|)`.
- **Map:** evaluation of its PGF at `z`.
- **Preserved predicate:** `Z_V(z)>z` implies a positive-measure zero-depth
  set; allowing some rational `z` makes this an equivalence.
- **Destroyed information:** phase, components, endpoint owners, and every
  measure-zero height-`h` witness.
- **Needed sidecar for weak LRC:** the labelled phase-height/Euler carrier of
  THM-2047, or an equality/rigidity classification of the `p_0=0` locus.

## 2. The exact Cayley/artanh certificate

All endpoints of the depth cells have the form `(m+/-h)/v`.  Hence every
`p_k` is rational, and `Z_V(z)` is rational when `z` is rational.  Whenever
`Z=Z_V(z)>z`, set

```text
t=(Z-z)/(Z+z) in Q intersect (0,1).
```

Then exactly

```text
log(Z/z)=log((1+t)/(1-t))
        =2 artanh(t)
        =2 sum_{m>=0} t^(2m+1)/(2m+1)>0.                    (5)
```

Therefore the truncated odd-power lower bound and geometric-tail upper bound
used in HYP-9023 are a rigorous, float-free wrapper for quantitative Gibbs
free-energy comparisons.  For the mere sign in (3), however, the rational
comparison `Z>z` is already stronger and simpler: artanh certifies a
previously constructed gap; it does not create the LRC map.

The smallest hostile/positive pair is exact.  At `h=1/3`,

```text
V={1,2}:  (p_0,p_1,p_2)=(0,2/3,1/3),
           Z(z)=(2z+z^2)/3<z;                               (6)

V={1,3}:  (p_0,p_1,p_2)=(1/9,4/9,4/9),
           Z(z)=(1+2z)^2/9,
           Z(z)>z iff 0<z<1/4.                              (7)
```

At `z=1/10`, (7) gives `Z=4/25`, `Z/z=8/5`, and

```text
log(Z/z)=2 artanh(3/13).
```

This is an exact source-to-target instance of the requested mechanism.

### Verdict on the supplied two-log inequality

Write

```text
X=(2457/6592) log(8847357/2974400)-log(1285/896).
```

The supplied artanh sandwich proves `X>1/25` with the stated exact rational
margin.  But no map from an arbitrary speed row `V` to these two rational
arguments has been established.  Moreover `X` cannot itself be an
**unscaled single** rational Gibbs ratio `log(Z/z)`: the `3`-adic valuation of

```text
exp(X)=(8847357/2974400)^(2457/6592)/(1285/896)
```

is `2457/6592`, not an integer, whereas a nonzero rational has integral
prime valuations.  The obstruction disappears after multiplying by `6592`,
so a normalized free-energy difference is not ruled out; it still needs the
missing row-to-partition-function construction.

At the hostile LRC wider-gap threshold `h=3/41`, the external reciprocal
ratios

```text
z_A=896/1285,       z_B=2974400/8847357
```

are not plug-and-play Gibbs temperatures.  Exact endpoint sweeps on eight
fixed defect-seven rows give `Z(z_A)<=z_A` and `Z(z_B)<=z_B` in all eight
cases, while `z=1/10` certifies all eight.  This does not refute a derived
two-temperature inequality; it refutes using the two snippet arguments as
universal bare fugacities.

## 3. Sound max/free-energy maps, and their sharp equality obstruction

For any probability measure `mu` on `T` and any `beta>0`,

```text
A_mu(V)=integral g_V dmu <= Gap(V),                          (8)
F_{beta,mu}(V)=(1/beta) log integral exp(beta g_V) dmu
                 <= Gap(V).                                 (9)
```

Thus either `A_mu(V)>=h` or `F_{beta,mu}(V)>=h` is a sound
source-to-target certificate for a weak lonely point.

If `mu` has full support, both inequalities are strict.  Indeed `g_V(0)=0`
and `Gap(V)>0`, so continuity gives an open set on which
`g_V<Gap(V)`.  Full support gives that set positive `mu`-mass, making both
the average and the exponential average strictly smaller than their pointwise
upper bounds.

For the tight arithmetic progression

```text
V_n={1,...,n-1},       Gap(V_n)=1/n,                         (10)
```

where (10) follows from the witness `t=1/n` and the pigeonhole principle on
`0,t,...,(n-1)t`, every finite-temperature/full-support certificate in
(8)--(9) is strictly below `1/n`.  This is the sharp obstruction:
finite-temperature smoothing can handle a bulk separated from the boundary,
but cannot by itself certify the exact tight row.  A zero-temperature limit
with uniform quantitative control, or an exact rigidity sidecar, is required.

## 4. The signed log-Riesz pairing is not a loneliness certificate

The superficially attractive functional

```text
J_V(R)=integral_T H(t) log R(t) dt                           (11)
```

has no analogue of the pointwise implication in section 1.

For every `n>=3`, take the full-cover/tight row `V_n` at `h=1/n` and the
normalized one-factor Riesz density

```text
R(t)=1-(1/2)cos(2 pi t)>0,       integral R=1.               (12)
```

This factor uses the runner frequency `1`, not an extraneous mode.  Since
`log x<=x-1`,

```text
J_{V_n}(R)
 <= -(1/2) integral H(t)cos(2 pi t)dt
 =  -sin(2 pi/n)/(2 pi)<0.                                  (13)
```

For the equality, every `v>=2` danger comb is `1/v`-periodic and hence
orthogonal to `cos(2 pi t)`; the `v=1` comb contributes
`sin(2 pi/n)/pi`.  Nevertheless `p_0=0`.  Thus even a negative signed
log-Riesz pairing, with a legal runner-frequency factor, does **not** imply a
strict lonely point.

There is a stronger support warning.  If auxiliary Riesz frequencies are
unrestricted, (11) is unbounded below for every nonzero nonnegative
`H in L^1(T)`.  Fix `0<a<1` and put

```text
c_a=integral log(1+a cos(2 pi t))dt
   =log((1+sqrt(1-a^2))/2)<0.
```

By oscillatory averaging, frequencies `w_j` can be chosen recursively both
superincreasing and so large that

```text
integral H(t)log(1+a cos(2 pi w_j t))dt
   < (c_a/2) integral H.
```

Superincrease makes the frequencies dissociated, so

```text
R_N=product_{j=1}^N(1+a cos(2 pi w_j t))
```

is positive and has integral one.  Additivity of `log R_N` now gives
`J_V(R_N)->-infinity`.  A frequency/support sidecar is therefore mandatory
even for using (11) as a heuristic optimizer.

## 5. Final soundness ledger

| functional | implication | exact boundary | lost coordinate |
|---|---|---|---|
| `integral H R`, `R>=0` | the usual linear Riesz cover certificate | sound | concentration quality |
| `Z_V(z)=integral z^H` | `Z_V(z)>z => p_0>0`; complete after `exists z` | misses `p_0=0` weak witnesses | phase/owners/Euler class |
| `integral g_V dmu` | lower-bounds `Gap(V)` | strict for full-support `mu` | maximizer phase |
| `beta^-1 log integral exp(beta g_V)` | lower-bounds `Gap(V)` | strict at finite temperature | ground-state/equality data |
| `integral H log R` | no loneliness implication | false already on every tight AP | positivity of the test measure |

The artanh expansion is therefore genuinely useful, but in a precise role:
it certifies rational log-ratio inequalities after a sound carrier has
produced them.  The remaining LRC(14) dilemma is not transcendental
numeration.  It is to prove either strict `p_0>0` through a full-depth
partition estimate or nonemptiness of the labelled height-`1/14` boundary
through an exact rigidity sidecar.
