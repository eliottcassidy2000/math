# LRC14 Exposure-Poset Proof Pass

I tried a different proof object for the LRC14 summit: not a new scalar, and
not another raw tournament class, but an exposure poset.

The strict counterexample is now phrased as an unexposed source packet:

```text
no q-witness
no AP/GW taut boundary atom
no positive Haar bridge
no Fejer/Toeplitz PSD failure
no K33/C27/petal state-lift label
no moment/twist/Ramanujan handoff
```

New artifact:

```text
04-computation/lrc14_exposure_poset_creative_pass_codex_20260624.py
05-knowledge/results/lrc14_exposure_poset_creative_pass_codex_20260624.out
05-knowledge/hypotheses/HYP-2988-lrc14-exposure-poset-proof-pass.md
```

Default bounded audit:

```text
rows                         12015
zero-safe rows               2: AP, GW 12->24
positive-safe rows           12013
positive rows with Fejer dual 12013
unexposed rows               0
```

The hard q>=14 front is exactly the expected one:

```text
P10+GW                    degree 280
drop(6)->63               degree 266
near/K33 12->36           degree 159
P10+K33                   degree 124
drop(6)->86               degree 134
two drop(12,13)->add(14,29) included as the small-margin S162 face
```

Tournament Analysis uses exposure channels as vertices:

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

The tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3_cycles=0
```

The proof target I would carry forward:

```text
Every primitive q>=14 non-AP/GW source packet either has positive open Haar
exposure with a familywise Fejer interval certificate, or carries a named
C27/K33 state-lift label.
```

This is not a proof of LRC14.  It is a sharper way to say what the next proof
must do.  HYP-2981 already suggests the interval arithmetic is finite-looking;
HYP-2988 says the exact object to certify is the no-hidden-exposure kernel,
not the raw list of rows.

Rebase note: this should now be read downstream of HYP-2982/HYP-2983/HYP-2984
and the HYP-2985 admissible-smoothing dispatcher, plus the incoming HYP-2986
tope-wall/cocircuit lane and HYP-2987 certificate-handoff atlas.
The no-hidden kernel is the simultaneous residual after analytic quotient
guardrails, Kaczynski smoothing labels, kernel-homotopy boundary-defect atoms,
Farey `e=14p-q` dispatch, Haar-open exposure, and Fejer/Toeplitz duals have all
had their labels restored.
