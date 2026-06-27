---
id: HYP-2911
status: RESOLVED in the apex-majority branch by THM-571; retained as a guardrail for raw 14-shift sieves
source: codex-2026-06-22-S121
tags: [lrc14, apex-majority, halfstep, phase-collision, node3, open-q-108]
depends_on:
  - THM-568
  - THM-569
  - THM-570
  - THM-571
related:
  - HYP-2910
  - OPEN-Q-108
results:
  - 04-computation/lrc14_apex_majority_shift_guard_codex_s121.py
  - 05-knowledge/results/lrc14_apex_majority_shift_guard_codex_s121.out
---

# HYP-2911: half-step phase collision guardrail

THM-570 closes the raw fourteen-shift argument for

```text
S = 14Q union R,     |Q|>=7, |R|<=6,
```

unless at least two residual speeds in `R` are divisible by `7` but not by
`14`.

The reason is exact.  After the below-14 theorem gives a strict safe point `u`
for `Q`, every lift `t=(u+k)/14` keeps the `14Q` speeds safe.  Ordinary residual
speeds (`gcd(r,14) != 7`) forbid at most two shifts and at most one in each
parity.  A single half-step speed (`gcd(r,14)=7`) forbids at most one whole
parity class, so six residual speeds still leave at least two surviving shifts.

The obstruction starts exactly at two half-step speeds.  The exact guardrail is

```text
u = 2/49:
  r=7    forbids all even shifts,
  r=161  forbids all odd shifts.
```

So a proof cannot depend only on mod-14 residues or a raw fourteen-shift
pigeonhole.  THM-571 supplies the missing move for the actual apex-majority
branch: descend to the `7`-phase.  Once two half-step residuals exist, at least
nine speeds are multiples of `7`; after scaling by `7`, LRC<=13 gives a strict
safe phase and at most four non-`7` residual speeds remain.  On the seven lifts,
each such residual forbids at most one point, so a lift survives.

The remaining value of this hypothesis is as a guardrail.  Any proof that stays
at the fourteen-lift level must use at least one of:

1. the positive interval `J` around the scaled `Q` witness, not just one `u`;
2. arithmetic restrictions imposed by the actual `Q` block;
3. a finite phase-collision exclusion for two or more odd multiples of `7`;
4. an equidistribution or exact-period packet estimate on the half-step phases.

## Former concrete target, now discharged for apex-majority

Prove that for every `Q` with `7<=|Q|<=12` and every residual half-step packet
`H={7a_1,...,7a_h}` with `2<=h<=6`, the strict `Q`-safe interval cannot be
contained in the alternating parity-cover set

```text
{x : the half-step speeds forbid both parity classes among k=0,...,13}.
```

If this fails as stated, strengthen it by adding the ordinary residual speeds
and prove the full residual safe set still meets the `Q`-safe interval.

THM-571 takes the second route: it changes the phase variable from `14` to `7`
and proves the full residual safe set meets the lifted safe interval.

## Assumption challenge

The failed assumption was "residue class modulo 14 determines the obstruction."
The speed `161 = 7 + 14*11` has the same half-step residue as `7`, but its
offset includes an extra `11u`; at `u=2/49` that offset flips the forbidden
parity.  The right carrier is therefore a phase-labelled half-step packet,
not a residue-only tournament.
