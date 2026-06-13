# Next Session Brief: Cartan Bridge + Tournament Probe + Adelic Space

## User's Request
Explore creative applications of the Cartan bridge, build and improve the tournament probe, and connect these to the adelic tournament space.

## Key Context

### The Cartan Bridge
The gl(4,R) decomposition splits any 4×4 matrix into:
- **Antisymmetric (6 dims)** = tournament structure (competition, who beats whom)
- **Symmetric (10 dims)** = metric structure (cooperation, self-knowledge)

This is the "Cartan decomposition" of the Lie algebra. The bridge connects:
- Tournament theory (our antisymmetric/competition side)
- Riemannian geometry (the symmetric/cooperation side)
- Napolitano's "dark modes" (the symmetric subspace carries all correctness info)

### The Tournament Probe
Three zero-parameter tools in `llm_tournament_tools_s93c.py`:
1. **TournamentConfidence**: conf = tanh(logit_gap). O(1) hallucination risk.
2. **AdaptiveDepthController**: skip layers when gap > threshold. 67% savings.
3. **IntransitivityDetector**: track top-k ranking stability. Cycle = confusion.

Plus the `TournamentWrapper` class that wraps any transformer.

### The Adelic Tournament Space
From `adelic_tournament_s91a.py` and `adelic_geometry_s91b.py`:
- A_T(n) = R × Π_{p|D} Z/p^e Z (truncated adelic space)
- D(n) = odd part of C(n,2) is the conductor
- The flip chain IS a Hecke operator on this space
- Eigenvalues decompose via CRT into local factors at each prime

## What To Explore

1. **Cartan decomposition of attention matrices**: For each attention head, decompose QK^T into antisymmetric (tournament) and symmetric (metric) parts. Does the metric part predict correctness?

2. **Tournament probe refinement**: Improve the IntransitivityDetector by tracking the FULL Cartan decomposition, not just the logit gap. Use the 10 dark modes for confidence.

3. **Adelic connection**: The local factors Z/p^e Z at each prime correspond to different "scales" of self-knowledge. At prime 3 (Eisenstein/ramified): structural consistency. At prime 7 (split/forbidden): contradiction detection. Map the 10 dark modes to specific primes.

4. **Practical implementation**: Build a PyTorch module that extracts the 6+10 decomposition from hidden states and outputs a calibrated confidence score.

## Key Files
- `04-computation/sixteen_dimensions_s93d.py` — the gl(4,R) decomposition
- `04-computation/llm_tournament_tools_s93c.py` — the tournament probe tools
- `04-computation/adelic_tournament_s91a.py` — adelic tournament space
- `04-computation/napolitano_assessment_s93b.py` — assessment of the "dark modes" paper
- `07-reflections/k-periodicity-synthesis.md` — the grand synthesis
