---
id: THM-866
title: AXIS LEVEL NO-HOLES COMPLETENESS — the populated x-levels are exactly the full step-8 progression from the (near-)regular floor to the transitive ceiling (n³−n)/3, proved by the F3 exchange walk (tie-splitting flips)
status: CLAIMED mac-mini-2026-07-15-S109 (proof drafted; verification + writeup this session)
source: mac-mini-2026-07-15-S109; owner directive 2026-07-15 ("prove the no-holes completeness via the F3 exchange walk"); closes opus-S316 handoff ("the level-law no-holes completeness wants the F3 exchange-walk proof")
depends_on:
  - THM-855 F3 (per-flip drop law Δx = 4(d_l − d_w) + 8)
  - opus-S316 probe (i) (mod-8 congruence, two-line proofs; exact census n=4..8)
---

# THM-866 — axis level completeness (STUB, claim)

**Claim.** For every n, {x(T) = Σ_v (2s_v−(n−1))² : T tournament on n} = {x_floor, x_floor+8, ..., (n³−n)/3}, x_floor = 0 (n odd) / n (n even).

**Proof shape (drafted):** (1) any non-transitive T has two tied vertices (all-distinct scores ⟹ transitive); (2) flipping the arc between tied vertices has Δx = +8 exactly (THM-855 F3 with d_w = d_l); (3) so the F3 exchange walk from any T climbs in exact +8 steps and halts only at the transitive tournament — giving simultaneously the ceiling extremality, the mod-8 congruence, and (started from the regular/near-regular floor) the population of EVERY level. Floor minimality: x ≥ 0 (odd n), x ≥ n (even n, all d_v odd). Writeup + machine verification this session.
