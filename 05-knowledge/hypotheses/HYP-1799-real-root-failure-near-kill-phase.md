---
id: HYP-1799
status: EXPLORATORY
source: opus-2026-05-30-S5
related:
  - THM-025
  - HYP-1785
  - HYP-1795
  - INV-189
---

# HYP-1799: Real-Root Failure Is Near-Kill Plus Phase Imbalance

## Statement

The first non-real-rooted Omega independence polynomials are controlled by a
mixed mechanism:

```text
small deletion residue with rank >= 2
plus an imbalanced root/phase channel.
```

The residue condition filters candidate failures.  The phase/root condition
explains why most small residues still do not fail.

## Evidence

THM-025 has

```text
I(Omega,x) = 1 + 94x + 10x^2 + x^3
max-loss deletion residue = (2,1,0)
residue rank = 2
near_kill_rank2_vertices = 1
```

The contrast table shows this differs from the two nearest misleading
families:

- `H=63` complete-core examples have exact kill residue rank `0`;
- Paley/interval examples have broad rank-2 residues, not near-kills.

So the suspected failure class is not "high deletion loss" alone and not
"rank 2" alone.  It is high deletion loss leaving a tiny rank-2 survivor, with
root coefficients then becoming too imbalanced for real-rootedness.

## Predictions

1. Non-real-rooted `n=9` examples should be enriched in
   `omega_near_kill_rank2_vertices > 0`.
2. Exact-kill complete-Omega examples should not by themselves predict
   non-real-rootedness.
3. Among near-kill rank-2 examples, Newton/root defects should correlate with
   the surviving pair/triple incidence pattern, not just `alpha_1`.

## Test Plan

1. Run an `n=9` real-root failure census or stratified sample.
2. Record residue features plus Newton inequality margins for `I(Omega,x)`.
3. Compare failure rates across exact-kill, near-kill rank-1, near-kill
   rank-2, and broad-residue classes.

## Sources

- `04-computation/residue_phase_incidence_contrast_s5.py`
- `05-knowledge/results/residue_phase_incidence_contrast_s5.out`
- `07-reflections/applied-residue-phase-incidence-programs.md`
