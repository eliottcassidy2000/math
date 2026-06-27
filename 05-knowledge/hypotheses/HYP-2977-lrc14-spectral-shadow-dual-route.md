---
id: HYP-2977
title: LRC14 spectral-shadow dual route
status: PROOF-ROUTE / global Fourier-dual certificate target; not a proof
source: codex-2026-06-24-spectral-shadow
script: 04-computation/lrc14_spectral_shadow_dual_codex_20260624.py
result: 05-knowledge/results/lrc14_spectral_shadow_dual_codex_20260624.out
related:
  - HYP-2976
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2971
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
  - HYP-2956
  - HYP-2949
  - HYP-2948
  - HYP-2889
  - HYP-2617
  - HYP-2614
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-2977: LRC14 Spectral-Shadow Dual Route

## Claim

The post-THM-571 LRC14 Moon core admits a second, global proof interface:
instead of locating a local endpoint pinch, treat the strict-safe open set

```text
U(S) = {t in R/Z : ||v t|| > 1/14 for every v in S}
```

as an indicator function `1_U` and search for a finite trigonometric dual
certificate:

```text
positive low-band Fejer / Beurling-Selberg shadow
  or controlled high-frequency relation-lattice tail
    -> Haar(U(S)) > 0.
```

This does not replace the labelled-packet route.  It supplies a different
closure theorem target for the same Moon core:

```text
qdiv>14, |S cap 14Z|<=6, not discharged by gK8/K33/fixed-margin packets
```

should either emit a positive Fourier minorant packet or force its spectral
mass into a labelled high-frequency relation tail already visible to
HYP-2968/HYP-2965/HYP-2614/HYP-2617.  After rebasing over the concurrent
HYP-2975/HYP-2976 work, read this as the strict-safe-set Fourier/Fejer member
of the dual stack: HYP-2970 endpoint-credit winding cycles, HYP-2971
multiplicity-moment barriers, HYP-2972 twist-ladder blockers, HYP-2973
danger-count moments, HYP-2974 Fourier-Toeplitz cover tests, HYP-2975
taut-bridge curvature, and HYP-2977 spectral shadows.  HYP-2976 is the
lineage/synthesis map around these proof interfaces rather than another
certificate.

## Computation

Script:

```text
04-computation/lrc14_spectral_shadow_dual_codex_20260624.py
```

Stored output:

```text
05-knowledge/results/lrc14_spectral_shadow_dual_codex_20260624.out
```

The script uses exact rational strict-safe components from the S146
Haar-Baire boundary scout, then computes numerical Fourier coefficients

```text
c_n = integral_U exp(-2*pi*i*n*t) dt
```

through bandlimit `H=224`.  It audits AP/GW, named low-frontier rows, covering
rows `6->98`, `12->84`, `12->168`, and the `18` tightest structured few-apex
rows from HYP-2968.

Default stored run:

```text
positive rows audited                         26
zero strict-mass rows                          AP, GW 12->24
smallest positive exact mass                   1/1260
smallest Fejer_14 midpoint value, positives    0.00604909
```

AP and Goddyn-Wong are the only zero-shadow atoms in the audit.  Every positive
row has a positive Fejer convolution at the midpoint of its largest exact
safe component.

The negative signal is just as important: this is not a cheap low-frequency
shortcut.  Many positive rows are not `90%` Parseval-captured by `H=224`.
For example:

```text
near 12->36        E<=224 = 0.177
petal 10->20       E<=224 = 0.226
covering 6->98     E<=224 = 0.518
covering 12->84    E<=224 = 0.414
```

The strongest spectral mass often sits in frequency packets tied to exact
resonance scales:

```text
near 12->36       dominant modes 41, 35, 76, 6, 117
petal 10->20      dominant modes 27, 54, 49, 76, 81
covering 12->84   dominant modes 6, 12, 18, 24, 30
covering 6->98    dominant modes 3, 6, 69, 66, 72
```

## Tournament Analysis

Vertices are Fourier bands, not runners:

```text
1-7, 8-14, 15-28, 29-56, 57-112, 113-224.
```

Pair observable:

```text
band A -> band B iff A has larger Parseval band energy than B on more
audited positive rows; ties follow the low-frequency Hamiltonian path.
```

Stored fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1}
c3=0
scc=(1,1,1,1,1,1)
hp=1
Hamiltonian path:
  113-224 > 57-112 > 29-56 > 15-28 > 1-7 > 8-14
```

The high-band dominance is a proof warning.  A spectral proof that only checks
the first few modes will miss the hard rows.  A viable spectral proof must
combine a one-sided Fejer/Beurling-Selberg packet with a relation-lattice tail
bound or a high-band packet classification.

## Assumption Challenge

Considered vertices:

```text
runners, gaps, fixed sections, section boundaries, wall-crossing events,
residues, cover arcs, lift packets, Fourier modes, and proof obligations.
```

Chosen quotient:

```text
Fourier bands and Fejer-dual packets.
```

Preserved:

```text
strict Haar mass, L2 mass of 1_U, spectral detectability,
and positive Fejer convolution witnesses.
```

Destroyed:

```text
endpoint-owner labels, C27 shell addresses, exact Farey maximizers,
and packet-family names.
```

Challenged assumption:

```text
the proof must be local at a boundary pinch.
```

HYP-2977 says the local and global views should meet as follows: local packet
labels explain which high-frequency resonances are allowed, while a global
dual minorant proves that any allowed nonzero shadow has positive safe mass.

## Proof Use

The theorem target is:

```text
Moon-core spectral dichotomy.

For primitive LRC14 rows in the post-THM-571 core, either
  (a) a bounded Fejer/Beurling-Selberg packet has a positive lower bound, or
  (b) the high-frequency energy lies in a relation-lattice packet that routes
      to HYP-2968 few-apex lifts, HYP-2965 boundary moments,
      HYP-2974 Toeplitz tests, HYP-2973 count moments,
      HYP-2971 multiplicity barriers, HYP-2970 endpoint-credit potentials,
      or HYP-2908/THM-572 state lift.
```

This connects the current labelled-packet proof stack with older support-6
Fourier-tail work: a high-band spectral proof should not sum absolute values
naively; it should use residue-addressed reciprocal cancellation and finite
relation-hyperplane ledgers.
