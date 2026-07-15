---
id: THM-865
title: THE LOCKER PARITY LAW — H(D_n) ≡ 1 (mod 4) for the locker (divisibility) tournament, via the divisor-pairing involution on odd cycles
status: CLAIMED mac-mini-2026-07-15-S109 (proof in progress this session; reduction chain fixed, involution under construction)
source: mac-mini-2026-07-15-S109; owner directive 2026-07-15 ("prove the locker parity law via the divisor-pairing involution"); closes THM-853(III) evidence-log item
depends_on:
  - THM-002 (OCF)
  - THM-466 (2-adic digit tower; m=2 corollary H ≡ 1 + 2α₁ mod 4)
  - THM-853(III) (D_n definition + conjecture, kind-pasteur S128 c20)
---

# THM-865 — the locker parity law (STUB, claim)

**Claim.** For every n, the locker tournament D_n (vertices 1..n; for u < v: v→u iff u|v or v=u+1, else u→v) satisfies H(D_n) ≡ 1 (mod 4).

**Fixed reduction (proved already by THM-466 m=2):** H ≡ 1 + 2α₁ (mod 4), so the law ⟺ α₁(D_n) = total number of directed odd cycles of D_n is EVEN.

**Fixed stratification (this session):** D_{m} is the induced subtournament of D_n on {1..m}, so α₁(D_n) = Σ_{m≤n} t_m where t_m = #odd cycles through the top vertex m in D_m. Claim reduces to: t_m even for every m.

**What still needs evidence:** the fixed-point-free involution on odd cycles through the top vertex (divisor-pairing d ↔ m/d on the exit arcs m→d, d|m, plus the base-path exit m→m−1). Computation + proof this session.
