---
id: HYP-2618
title: OCF noise/Condorcet spectra as an LRC(14) sequence bridge
status: OPEN
source: codex-2026-06-19-S15
depends_on:
  - HYP-2617
  - HYP-2615
  - HYP-2614
  - THM-538
related:
  - HYP-2616
  - HYP-2608
  - OPEN-Q-108
---

# HYP-2618 - OCF Noise/Condorcet Spectra as an LRC(14) Sequence Bridge

## Stub Reservation

This packet reserves HYP-2618 / T866 / S15 for the user's prompt:

- decide whether the OCF `H(T)=I(Omega(T),2)` is a noise-stability functional
  at a specific parameter;
- reinterpret forbidden values `{7,21}` as forbidden Condorcet-cyclicity
  spectra when tournaments are read as majority relations;
- use the "large absolute mass but tiny signed mass" clue from HYP-2614 through
  HYP-2617 to extract integer/fractional sequence spines that may clarify the
  LRC(14) support-six tail.

Planned computation:

- `04-computation/ocf_noise_condorcet_lrc_sequence_bridge_codex_s15.py`
- `05-knowledge/results/ocf_noise_condorcet_lrc_sequence_bridge_codex_s15.out`

The expected guardrail is that OCF is certainly a hard-core partition function
at activity `2`; the open question for this session is whether a nontrivial
noise-stability normalization gives a useful new social-choice/LRC invariant,
or whether only a degenerate/biased-density reading survives.
