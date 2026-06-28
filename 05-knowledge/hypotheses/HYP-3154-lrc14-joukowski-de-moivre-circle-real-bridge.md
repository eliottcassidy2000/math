---
id: HYP-3154
title: The Joukowski/De Moivre bridge sends the LRC Lee-Yang circle to the tournament real-rooted class
status: VERIFIED geometry + structural reframe; not an LRC14 proof
source: mac-mini-2026-06-27-S73 plus codex-2026-06-28 integration
script: 04-computation/lrc_joukowski_resolvent_macmini_S73.py
result: 05-knowledge/results/lrc_joukowski_resolvent_macmini_S73.out
reflection: 07-reflections/the-joukowski-de-moivre-bridge-lrc-circle-is-the-tournament-real-rooted-class.md
related:
  - HYP-3150
  - HYP-3151
  - HYP-3152
  - HYP-3153
  - HYP-3161
  - HYP-3132
  - HYP-3147
  - HYP-3128
  - THM-577
  - OPEN-Q-108
---

# HYP-3154: Joukowski/De Moivre Circle-Real Bridge

## Claim

The map

```text
w = z + R^2/z
```

sends the Lee-Yang circle `|z|=R` for the LRC miss-count PGF to the real
axis.  Therefore the k=8 De Moivre/quartic packet should be read as a
circle-to-real-axis stability guardrail, not merely as a lucky solvable
quartic.

## Verified Payload

The S73 scout verifies the exact geometry:

- `z = R exp(i theta)` maps to `w = 2R cos(theta)`, hence real.
- The uniform apex-7 ideal `1 + z + ... + z^6` has roots at the nontrivial
  seventh roots of unity.
- Its Joukowski image is the de Moivre angle set
  `2cos(2*pi*j/7)` for `j=1,2,3`.
- In the stored tests, `consec {0..7}` has lower off-circle `max|Im w|` and
  higher coverage than the sampled non-consecutive 8-sets; the doubled row
  `{0,2,4,6,8,10,12,14}` has the same payload.

This aligns with HYP-3152: coverage is a radius payload `q0 = q6*R^6`, while
the dip is off-circle deviation.  It also refines HYP-3150/HYP-3153: the
bounded-core degree ceiling is a root-locus certificate, and the sidecar is the
whole PGF/root curve rather than a scalar coverage value.

## Compression Lesson

The bridge gives a sharper version of the controlled-forgetting rule.  A
compression may keep only the circle radius, or only the real resolvent roots,
only when the target proof predicate factors through that summary.  Otherwise
the missing coordinate must be named:

```text
root_radius_R
root_angle_spread
off_circle_imaginary_defect
cyclotomic_7fold_deviation
Joukowski_real_axis_status
ordered_odd_Worpitzky_sidecar
```

This is the same legality rule seen in the user's function quartet:
`a+b` and `a*b` survive unordered compression, while `a^b` and `b^a` require
an ordered sidecar.

## Tournament Analysis Hook

The natural tournament vertices are proof carriers, not runners or raw arcs:

```text
Lee-Yang circle packet
Joukowski real-axis packet
7-cyclotomic ideal packet
off-circle dip packet
Worpitzky odd-edge packet
biquadratic even-fold packet
Asano obstruction packet
real-rooted Omega packet
```

Pairwise observable: whether one carrier preserves the LRC14 coverage predicate
after quotienting better than another.

Switch/gauge: orient an edge toward the carrier that keeps more of
`root-locus + ordered odd sidecar + finite-address exit`.

Tie Hamiltonian path:

```text
7-cyclotomic ideal packet
-> Lee-Yang circle packet
-> Joukowski real-axis packet
-> biquadratic even-fold packet
-> Worpitzky odd-edge packet
-> off-circle dip packet
-> Asano obstruction packet
-> real-rooted Omega packet
```

Challenged assumption: the real-rooted tournament-side theorem is not itself an
LRC proof.  The quotient preserves the target only if the miss-PGF or its
single-runner factors can be routed through a stability-preserving
Grace-Walsh-Szego/Asano step; otherwise the off-circle defect remains named
debt.

## Remaining Obligation

The useful reframe is:

```text
dip >= 0
<=> the Joukowski image stays real-rooted enough
<=> the LRC miss-PGF remains in the relevant stability class
```

What remains is to prove that stability statement for the LRC miss-count
factors, or to reduce every exception to the existing k=8 bounded-core
degree<=4 packet plus the HYP-3147 Worpitzky odd-edge sidecar.
