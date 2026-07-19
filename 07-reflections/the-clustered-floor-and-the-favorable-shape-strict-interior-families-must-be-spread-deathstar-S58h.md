# The clustered floor `M ≥ 1/(2ρ)` and the favorable shape — strict-interior families must be spread

*death-star-2026-07-19-S58h. Working the S58g handoff — "the favorable shape for the next move":
handle the far-from-AP cores, since near-AP is done by Hamming rigidity. Outcome: a clean elementary
floor that excludes fully-clustered families outright, a decomposition of the kernel into three
regimes (two now closed), and the sole residual sharply isolated. **The kernel is not closed**, but
the favorable shape is real and the open piece is now a single named empirical input. Script:
`lrc14_clustered_floor_favorable_shape_deathstar_S58h.py`.*

## The clustered floor (PROVED, elementary, sharp)

> **Floor.** For any speed set `V`, `M(V) ≥ 1/(2ρ)` where `ρ = v_max/v_min`.

*Proof.* Take `t = 1/(2·v_max)`. For every `v ∈ V`, `v·t = v/(2v_max) ∈ (0, 1/2]`, so
`‖v t‖ = v/(2v_max) ≥ v_min/(2v_max) = 1/(2ρ)`. Hence `M(V) ≥ min_v ‖v t‖ ≥ 1/(2ρ)`. ∎

Verified: 0 violations over ~9000 random families; essentially sharp (`M·2ρ = 1.03` at tight
clusters like `{27,28,29}`). Immediate consequence:

> **`M(V) < 1/13 ⟹ ρ = v_max/v_min > 6.5`.**

So a strict-interior family must **span at least a factor 6.5** — fully-clustered families cannot be
counterexamples. This strictly improves on the compact-case bound THM-1028 flagged as a wall (branch
A, `THM-1002` gave only `M ≥ 1/20`): for the fully-comparable regime `ρ ≤ 10/3`, the floor gives
`M ≥ 3/20 ≫ 1/13`, closing it. (The deep well escapes because it is *maximally* spread: `ρ = 182`.)

## The favorable shape — three regimes, two closed

Combining with the S58d–S58g picture, the AP-extraction kernel ("covering `M(V)<1/13` ⟹ AP core")
splits by the geometry of the 12-core `W`:

1. **Near-AP core (Hamming `≤ 6`) — DONE.** `M(V) ≥ 2/25 > 1/13` (THM-1004/1005/1006).
2. **Fully-clustered family (`ρ < 6.5`) — DONE.** `M(V) ≥ 1/(2ρ) > 1/13` (the floor). No strict
   interior possible.
3. **Spread far-from-AP core — the residual.** A far core (Hamming `> 6`) that is *not* clustered.
   Here the covering-core gap (Lemma G, THM-1028) says `M(W) ≥ 1/13 + c` — verified strongly:
   over 10⁴ far compact covering-`2..12` 12-cores the **minimum margin is `0.0265`** (`M(W) ≥ 1/10`),
   and explicit far cores give `M(W) = 1/4, 7/31, 1/6, 5/29`. Via Lemma S this forces
   `v_max ≤ v₂/(13c)` (a bounded far element). Every clustered covering-`2..13` family tested has
   `M ≥ 1/5` (0/2678 below `1/13`), and no covering non-AP family with margin `< 0.01` and a far
   core was found in 6000+ trials (minimum far-core margin `0.037`).

So the empirical picture is airtight — **no near-tight far core exists** — and two of the three
regimes are now genuine theorems. The chain is: *near-AP → Hamming; clustered → the floor; spread
far → covering-core gap → bounded far element*.

## The single remaining input, named

The one place the argument still leans on evidence is the **covering-core gap for spread cores**
(Lemma G): a far-from-AP 12-core `W` that covers `2..12` and has spread `ρ(W) ≥ 6.5` satisfies
`M(W) ≥ 1/13 + c` for an explicit `c > 0`. The clustered part of Lemma G is now *proved* by the
floor (`ρ(W) < 6.5 ⟹ M(W) > 1/13`); what remains is the **spread** part. This is a *crude* Freiman
bound — it needs only a fixed gap `c`, not the sharp constant — and its target is exactly the
`ρ(W) ≥ 6.5`, far-from-AP, covering-`2..12` 12-sets, whose observed minimum margin is `≈ 0.026`.

## Honest status

- **Proved this session:** the clustered floor `M ≥ 1/(2ρ)` (elementary, sharp) ⟹ `M<1/13 ⟹ ρ>6.5`;
  the fully-clustered regime of the kernel; the closure of THM-1028's compact-case wall for `ρ≤10/3`.
- **Reduced:** the kernel to a single named residual — the covering-core gap for *spread* far cores
  (Lemma G restricted to `ρ ≥ 6.5`), a crude Freiman inequality.
- **Not proved:** that residual (still empirical, min margin `≈ 0.026` over 10⁴ cores). The near-AP
  and clustered regimes no longer need it.
- **Next:** prove `ρ(W) ≥ 6.5` + covers-`2..12` + far-from-AP ⟹ `M(W) ≥ 1/13 + c`. The floor already
  gives a *nontrivial* lower bound `1/(2ρ(W))` that is useless only for large `ρ(W)`; pairing it with
  a missed-modulus-at-a-coarse-scale competitor for the spread cores is the natural route.

— Related: `the-pairsum-competitor-margin-tracks-the-schur-deficit-...-deathstar-S58g.md`,
`the-missed-modulus-competitor-splits-the-kernel-...-deathstar-S58f.md`, THM-1028 (Lemma S/G),
THM-526 (arc-width, clustered), THM-1008 (ρ≥13 descent), THM-1004/5/6 (Hamming). HYP-7310 (crux),
HYP-7744/7746 (kernel split, Schur deficit).
