# THM-960 — The rung lock at the lonely-runner threshold (death-star-2026-07-17-S46)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCRungLock.lean`, standard
trio; the sharpness witness needs only `propext`). Source: HYP-7225. (THM-959
left for opus's block-structure renumber.)

## Statement

A runner fails the safe band at `p/q` iff some integer witness `w` has
`14·|v·p − w·q| < q`. Then:

1. `rung_lock`: two failing runners with EXACT integer speed ratio `v_j = m·v_i`,
   `1 ≤ m ≤ 13`, have exactly locked witnesses `w_j = m·w_i` — zero choices.
   The strict chain `14|w_j − m·w_i|·q < (m+1)·q ≤ 14·q` survives through
   `m = 13` because both band bounds are strict. **The lock threshold is the
   lonely-runner denominator itself: lock ⟺ m ≤ 13 = n** (for general LRC(n),
   band 1/(n+1) locks ratios ≤ n).
2. `rung_lock_sharp_at_14`: the break one step past 13, kernel-checked on
   `v = 1 → 14, q = 15, p = 1` (witnesses 0 vs 1).
3. `rung_lock_chain3`: multiplicative composition — on any comparable block
   with consecutive integer ratios ≤ 13, the bottom witness determines all.
4. `ray_fails_above`: an exact ray (`q ∣ v_i·p`) forces every integer multiple
   to fail — upward band propagation on towers.
5. `carrier_instant_eq_of_shared_witness`: the ALL-INTEGER Farey collapse — two
   moduli `q₁, q₂ ≤ D` sharing one witness for one runner with `D² ≤ 7·v_i`
   carry the same instant (`p₁q₂ = p₂q₁`), via the exact identity
   `v_i(p₁q₂ − p₂q₁) = (v_i p₁ − w q₁)q₂ − (v_i p₂ − w q₂)q₁`. No rationals.
   Composed with the lock, THM-956's component map is CONCRETE on exact-ratio
   strata: component = bottom failing witness.

## Honest scope (recon, `runglock_recon_deathstar_S46.out`)

Lock confirmed on 6221 double-failing pairs, zero violations; m=13 exhaustively
clean (q < 500); m=14 breaks found. The companion **divisor-collapse heuristic
is REFUTED on compressed block towers** (deep instants are common there; K up
to 1489 at D = 300; carriers do not divide the tower scale). The salvage that
survives: the **carrier-count pigeonhole** — worst 134 carriers vs 299 window
moduli over 60 block families, so a deep-free q exists in every recon window;
at such q the census needs only `liveCount > 0` (the nucleus, unchanged).
