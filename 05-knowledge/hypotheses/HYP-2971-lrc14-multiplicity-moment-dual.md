---
id: HYP-2971
title: LRC14 multiplicity-moment dual
status: PROOF-INTERFACE / exact dual certificates for positive lonely mass; not a proof
source: codex-2026-06-24-S156
script: 04-computation/lrc14_multiplicity_moment_dual_codex_s156.py
result: 05-knowledge/results/lrc14_multiplicity_moment_dual_codex_s156.out
related:
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2967
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2960
  - HYP-2956
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2908
  - THM-523
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-2971: LRC14 Multiplicity-Moment Dual

This route complements HYP-2970's endpoint-credit winding-cycle dual,
HYP-2972's twist ladders, HYP-2973's danger-count moment duals, and HYP-2974's
Fourier-Toeplitz PSD scout by taking an even coarser Farkas shadow.  It
deliberately forgets apex charts, lift packets, endpoint owners, and
safe-interval locations.  For a 13-speed row `S`,
define the integer danger multiplicity

```text
X_S(t) = #{v in S : ||v t|| < 1/14}.
```

The strict lonely set has positive Haar measure exactly when

```text
Pr[X_S=0] > 0.
```

## Dual Lemma

Call a polynomial `P` admissible if

```text
P(0) < 0,
P(k) >= 0 for k=1,...,13.
```

If `E[P(X_S)] < 0`, then `Pr[X_S=0] > 0`, hence `M(S)>1/14`.

Equivalently, every strict counterexample must satisfy the infinite cone of
moment inequalities

```text
E[P(X_S)] >= 0
```

for every admissible `P`.  This is a Farkas-style dual certificate for lonely
mass.

## Computation

Script:

```text
04-computation/lrc14_multiplicity_moment_dual_codex_s156.py
```

Stored output:

```text
05-knowledge/results/lrc14_multiplicity_moment_dual_codex_s156.out
```

The script computes exact multiplicity histograms by sweeping all danger-arc
endpoints and then searches two bounded certificate families:

```text
B_m(x) = binom(x-1,m),  m odd,
```

with `B_m(0)=-1` and `B_m(k)>=0` for `k>=1`, plus root barriers with `P(0)=-1`
and selected integer roots, checked on all `k=1..13`.

Default audit:

```text
one-swap add<=160
two-swap add<=36
three-swap add<=24
q_min=14
root_degree<=7
```

Readout:

```text
rows = 17104
zero_safe_rows = 2
positive_safe_rows = 17102
certified_by_barrier_family = 17102
no_barrier_found_with_selected_family = 2
min_positive_safe_mass = 1/1260

certificate_degree_hist:
  binomial degree 13: 50
  root degree 4:      11
  root degree 5:      5116
  root degree 6:      10886
  root degree 7:      1039
```

The only zero-safe rows in the audit are AP and Goddyn-Wong `12->24`.

## Named Calibration

AP and GW have no admissible negative moment certificate because their strict
safe mass is exactly zero.

Positive rows do certify:

```text
12->36:       B13, E=-1/1260
10->20:       B13, E=-1/980
13->26:       B13, E=-1/182
12->84:       B13, E=-563/105105
12->168:      root 1,2,3,6,7,11,12, E=-209101/83243160
drop(2,6)->add(17,42):
               root 1,2,3,10,11,12,13, E=-148763/69369300
drop(4,6)->add(19,42):
               root 1,2,3,10,11,12,13, E=-75947/1318016700
```

This is a different proof signature from the endpoint-owner route: the
covering repairs that were apex-blocked still have global multiplicity
separators.

## Interpretation

A counterexample must be moment-indistinguishable from the AP/GW zero-safe
atoms before packet labels are even consulted.  In particular, it must defeat
all admissible low-degree or bounded-root separating polynomials.

The current finite evidence suggests a theorem target:

> Multiplicity moment rigidity for LRC14.  Outside the AP/GW boundary atoms,
> every primitive `qdiv>=14` source-core row has a negative admissible moment
> barrier, preferably from odd binomial barriers or bounded root barriers,
> unless it exposes a genuinely new NORK/F7 packet or constructs the
> HYP-2908/THM-572 state lift.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
integer multiplicity histogram
> dual moment barrier
> AP/GW zero-safe atom
> NORK pinch template
> qdiv witness
> apex aperture
> raw safe mass
> raw runner set
```

Pairwise observable:

```text
preserves no-lonely counterexample status, exact Haar distribution,
low-moment inequalities, AP/GW atom visibility, packet labels, and
resistance to scalarization.
```

Switch/gauge:

```text
majority over six retention coordinates; ties follow the displayed order.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
```

This quotient destroys endpoint location on purpose.  Its use is not to
replace HYP-2965/HYP-2966/HYP-2968/HYP-2970/HYP-2972/HYP-2973/HYP-2974, but to give a global
moment obstruction that any zero-open packet must satisfy before those labelled
routes begin.
