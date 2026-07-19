---
id: THM-1226
title: THE GCD-PERIOD PROJECTIVE-CHARGE OBSTRUCTION — exact vertex loads, a factorial strict-high counterfamily, and the disconnected-spectrum rescue
status: PROVED (exact arithmetic obstruction and protected-needle embedding; conditional disconnected-G_gt transfer using THM-1221's complete projective branches)
source: codex-2026-07-19-S82
depends_on: [THM-1166, THM-1221]
related: [HYP-7678, HYP-7870]
script: 04-computation/lrc14_gcd_period_projective_charge_obstruction_referee_codex_S82.py
output: 05-knowledge/results/lrc14_gcd_period_projective_charge_obstruction_referee_codex_S82.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCGCDPeriodProjectiveCharge.lean
script_sha256: PENDING
output_sha256: PENDING
formalization_sha256: PENDING
---

# THM-1226 — the gcd-period projective-charge obstruction

For an edge `s_i=g a`, `s_j=g b`, `(a,b)=1`, put

```text
eta_ij=rho_ij(1-rho_ij),
kappa_ij=eta_ij ab/(a+b).
```

Then the periodic positioning error factors exactly as

```text
eta_ij/g=kappa_ij(1/s_i+1/s_j).                         (1)
```

Consequently tree positioning error is an inverse-speed-weighted vertex-load
sum.  No absolute bound by `C sum_i 1/s_i` follows from the THM-1221 strict
spectrum: an explicit pairwise-coprime translated factorial family has
complete strict-high graph and makes every tree clear `15/154`, while every
tree error stays bounded below and its harmonic sum tends to zero.

This is an obstruction to the gcd-period discrepancy proof route, not to the
localized F7 statement.  The same family admits a covered core-safe interval
on which all seven danger combs are simultaneously active, exposing the phase
information discarded by the unsigned period error.

The exact finite projective branches of THM-1221 still give a conditional
positive result: when `G_gt={rho>1/63}` is disconnected, one may take

```text
C_disc=448916/194775
```

in the positioned-tree error bound.  The connected `G_gt` branch is exactly
the unbounded-projective-height frontier.

Detailed proof, executable certificate, formalization boundary, tournament
audit, and frozen hashes follow in this claimed namespace.
