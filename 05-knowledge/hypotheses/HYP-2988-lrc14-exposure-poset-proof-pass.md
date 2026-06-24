---
id: HYP-2988
title: LRC14 exposure-poset proof pass and no-hidden-kernel target
status: PROOF-INTERFACE / bounded exposure-channel audit and theorem target; not a proof
source: codex-2026-06-24-exposure-poset
artifacts:
  - 04-computation/lrc14_exposure_poset_creative_pass_codex_20260624.py
  - 05-knowledge/results/lrc14_exposure_poset_creative_pass_codex_20260624.out
  - 07-reflections/lrc14-exposure-poset-creative-pass-codex-20260624.md
related:
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2984
  - HYP-2983
  - HYP-2982
  - HYP-2981
  - HYP-2980
  - HYP-2979
  - HYP-2978
  - HYP-2977
  - HYP-2976
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2966
  - HYP-2965
  - HYP-2963
  - HYP-2956
  - HYP-2953
  - HYP-2908
  - THM-523
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-2988: LRC14 Exposure-Poset Proof Pass

This pass tries a different proof object from the previous holistic synthesis.
Instead of asking for one more invariant, it treats each known proof route as
an exposure channel and asks whether any audited source packet remains
unexposed after all channels are allowed to retain their labels.

The proposed theorem target is:

```text
Every primitive LRC14 source packet either has a q-witness, is one of the
AP/Goddyn-Wong taut boundary atoms, exposes a positive Haar bridge with a
familywise Fejer interval certificate, or carries a labelled K33/C27/petal
handoff whose state-lift theorem discharges it.
```

The strict-counterexample falsifier becomes:

```text
primitive qdiv>14 source packet
  + zero strict open mass
  + not AP/GW
  + no positive endpoint bridge
  + no Fejer/Toeplitz PSD violation
  + no K33/C27/petal state-lift label
  + no moment/twist/Ramanujan handoff.
```

## Computed Audit

The script
`04-computation/lrc14_exposure_poset_creative_pass_codex_20260624.py`
builds a bounded AP-neighborhood and hard-frontier bank, then computes exact
Haar safe components and explicit Fejer-vector PSD certificates.

Default run:

```text
rows                         12015
zero-safe rows               2: AP, GW 12->24
positive-safe rows           12013
positive rows with Fejer dual 12013
unexposed rows               0
rows carrying unknown kernel 0
```

The exposure channels are:

```text
Q_WITNESS
AP_GW_TAUT_BOUNDARY
OPEN_HAAR_BRIDGE
FEJER_PSD_DUAL
K33_STATE_LIFT
C27_PETAL_EXIT
LATE_COVERING_PRESSURE
HARD_FEJER_MARGIN
UNEXPOSED_SOURCE_KERNEL
```

Highest-severity q>=14 rows are exactly the expected hard packet faces:

```text
near/K33 12->36                 q=14, mu=1/1260, degree=159
P10+K33                         q=14, mu=4/2205, degree=124
drop(6)->63                     q=14, mu=1/294, degree=266
P10+GW                          q=14, mu=1/980, degree=280
drop(6)->86                     q=14, mu=1249/258258, degree=134
two drop(12,13)->add(14,29)     included as the S162 small-margin face
```

This does not prove LRC14, but it tightens the current proof obligation:
the hard rows are not mysterious; they are labelled positive-open packets
whose formal gap is a familywise interval-certificate assembly.

## Tournament Analysis

Vertices are certificate/exposure channels, not runners or arcs.

Pairwise observable:

```text
severity-weighted rows exposed by one channel and not the other,
with retention of q-threshold, boundary status, open interval,
harmonic dual, packet label, state label, and anti-scalarization.
```

The switch/gauge gives:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3_cycles=0
sccs are all singletons
```

The score order is:

```text
OPEN_HAAR_BRIDGE
FEJER_PSD_DUAL
Q_WITNESS
LATE_COVERING_PRESSURE
HARD_FEJER_MARGIN
C27_PETAL_EXIT
K33_STATE_LIFT
AP_GW_TAUT_BOUNDARY
UNEXPOSED_SOURCE_KERNEL
```

This should not be read as "open mass is more fundamental than AP/GW."
It is a finite exposure audit: on the sampled bank, open-Haar and Fejer
channels expose the broad positive frontier, while AP/GW and K33/petal labels
remain small but load-bearing boundary/state channels.

## Creative Proof Shape

The live proof should be written as a no-hidden-kernel lemma:

```text
If a labelled source packet has no q-witness and is not AP/GW, then it must
expose either positive open mass, a familywise Fejer certificate on that open
mass, or a named C27/K33 state-lift label.
```

This is deliberately stronger than a floating Fejer audit and weaker than
full LRC14.  It asks for the exact bridge between:

```text
positive endpoint bridge
  -> rational center / safe component
  -> divisor-curried Fejer atom bank
  -> interval-enclosed negative margin
  -> packet-family certificate
```

HYP-2981 already suggests the numerical burden is finite.  HYP-2988 says the
right family target is the exposure kernel, not the raw row list.

## Rebase Integration

After rebasing over the HYP-2982/HYP-2983/HYP-2984/HYP-2985/HYP-2986/HYP-2987 stack,
read this pass as the router sitting above those proof modules.  HYP-2982 supplies the analytic
weight guardrail: squarefree and inverse-unit quotients are useful only when
lost prime-power, endpoint, Ramanujan, Fejer, smoothing, or state labels are
restored.  HYP-2983 supplies the Kaczynski/exponential-sum module map: a
smoothed certificate must retain its boundary approach class.  The HYP-2984
kernel-homotopy lane says a deformation is safe only when it preserves a
labelled packet certificate or emits a boundary-defect atom, while the
HYP-3001 Farey scheduler says exact `M=p/q` and `e=14p-q` dispatch the
unit-excess packets into AP/GW, C27 petal/two-block, or K33/state-lift routes.
HYP-2985 adds the admissible-smoothing dispatcher: every smoothing policy must
name which labelled packet family it is allowed to certify.
The incoming HYP-2986 tope-wall/cocircuit and certificate-handoff lanes add
the endpoint-cell and zipper-arrow versions of the same question: which
labelled exit prevents a packet from becoming a no-tope/no-cocircuit
unexposed residual?
After the handoff atlas moved to HYP-2987, treat that HYP-2987 atlas as the
zipper-arrow sibling of this exposure audit.

In that language, `UNEXPOSED_SOURCE_KERNEL` is not just a missing row witness.
It is the simultaneous kernel of the analytic quotient, smoothing-boundary,
kernel-homotopy, Farey-dispatch, Haar-open, and Fejer/Toeplitz channels after
all required labels have been retained.  The bounded audit found no such
kernel.

## Caveats

- The audit is bounded and AP-neighborhood biased.
- Fejer signs are still floating in this script; HYP-2981 is the interval
  formalization route.
- Channel severity is a heuristic theorem-target finder, not a canonical
  invariant.
- The output intentionally destroys full endpoint-owner incidence.  That data
  must be restored before any formal proof.

## Next Steps

1. Turn the highest-severity rows into packet-family interval certificates.
2. Add endpoint-owner labels to each Fejer center, not only row/family labels.
3. Prove the largest safe component in each non-AP/GW q>=14 family supports a
   Fejer certificate before the first invisible speed, or else emits K33/C27
   state-lift data.
4. Treat any future `UNEXPOSED_SOURCE_KERNEL` row as the first real F7
   Johnson-harmonic/source-spectrum falsifier.
