---
id: HYP-2974
title: LRC14 Fourier-Toeplitz dual scout
status: PROOF-INTERFACE / harmonic dual necessary condition; fills S157 stub; not a proof
source: codex-2026-06-24-S157 + codex-2026-06-24-S156
related:
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
  - HYP-2961
  - HYP-2956
  - HYP-2954
  - HYP-2953
  - HYP-2901
  - HYP-2908
  - THM-406
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_fourier_toeplitz_dual_scout_codex_s156.py
  - 05-knowledge/results/lrc14_fourier_toeplitz_dual_scout_codex_s156.out
---

# HYP-2974: LRC14 Fourier-Toeplitz Dual Scout

This hypothesis fills the S157 Fourier-Toeplitz PSD stub with the S156
computation.  It is the Fourier/Toeplitz member of the current
dual-certificate cluster, alongside HYP-2970 endpoint-credit winding cycles, HYP-2971
multiplicity-moment barriers, HYP-2972 twist ladders, and HYP-2973
danger-count moment duals.  Instead of classifying the residual first, look
directly at the covering formulation.

For a speed set `S`, define the danger multiplicity

```text
D_S(t) = sum_{v in S} 1_{||v t|| < 1/14}.
```

If `S` were a strict LRC14 counterexample, then the danger arcs cover the circle
and `F_S(t)=D_S(t)-1 >= 0` almost everywhere.  Hence every Toeplitz moment matrix

```text
T_K(S) = (hat F_S(i-j))_{0 <= i,j <= K}
```

must be positive semidefinite.  Any negative eigenvalue is a Fourier/Farkas
dual certificate that `S` has a safe open interval.  This is a necessary
condition for counterexamples, independent of endpoint-owner packet labels.

## Evidence

The default scout audits `52` rows: curated AP/GW, K33, petal, covering,
few-apex, and AP-mutation rows plus deterministic qdiv>=14 samples.

```text
rows audited                 52
boundary-only rows           2
positive-open rows           50
dual-certified rows          48
PSD-through-degree-90 rows    4
```

The two boundary-only rows, AP and Goddyn-Wong `12->24`, remain PSD through
degree `90`.  That is exactly what the dual condition should do: it does not
falsify known equality atoms.

The two positive-open rows not certified by degree `90` are also meaningful:
`K33 near 12->36` and `two-hole P10+GW`.  Those are already named exits in the
existing packet stack.  Thus the new harmonic lens appears to discharge many
covering/migration rows while naturally handing K33/petal exceptions back to
the state-lift and petal-rigidity machinery.

## Proof Target

A possible theorem shape is:

```text
After AP/GW boundary atoms and named K33/petal exits are removed,
every primitive qdiv>=14 LRC14 residual has a bounded-degree negative
Toeplitz moment certificate for F_S=D_S-1.
```

The bounded-degree phrase is currently empirical.  The proof would need a
structured harmonic argument, likely tracking the unit-apex residue band that
dominates the observed negative eigenvectors.

## Tournament Analysis

Vertices are harmonic proof lenses rather than runners, arcs, or packets.
The pairwise observable is which lens preserves the implication

```text
danger cover => F_S >= 0 => Toeplitz PSD
```

while retaining enough labelled information for known exits.  The default
fingerprint is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
c3=0
```

The challenged assumption is that the endgame must be packet-first.  HYP-2974
says a counterexample might first be ruled out by a trigonometric moment
obstruction, with packets reattached only for the few harmonic-invisible
exceptions.
