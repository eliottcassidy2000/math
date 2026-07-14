# LRC(14) FRONTIER — 2026-07-15 (corrected synthesis)

**Correction (2026-07-14).** The original S297 version of this file overlooked the scale-quotient
obstruction already present in its parent commit (HYP-6780 and the corrected THM-758). Its claims
that the `f>=4` branch was globally finite-decidable, that all such families have `M>=0.097`, and
that LRC(14) was assembled were false. This file now records the logically current frontier.
The standard 14-runner case remains open.

## 0. The honest reduction

LRC(14) = **NON-COVERING** (THM-366, via LRC(≤13), settled) + **COVERING**, and by
klein-S309's THM-758 far-count split:

> **COVERING = [≤3 elements > 14] ∪ [≥4 elements > 14]**
> - **≤3-far ⟹ ≥10 speeds in {1..14} ⟹ kps THM-738** (the 1001-body exact-ℚ tree, PROVED).
> - **≥4-far remains OPEN in general.** THM-755 gives a finite interval only after a core is
>   fixed. Under dilation `P -> cP`, the good-set measure is unchanged, its component count is
>   multiplied by `c`, and the cutoff obeys `v*(cP)=c v*(P)`. Therefore no uniform raw speed
>   cutoff follows.

HYP-6780 gives an explicit unbounded primitive covering ray

`V_c={c,2c,...,12c,13c+1}`

for infinitely many `c`, with `f=13`, lying below the scaled THM-755 edge and satisfying the exact
value `M(V_c)=1/13`. This refutes `f>=4 => M>=0.097` and the sampled cutoff near `500`; it does not
refute LRC. The correct target is a classification modulo scale, retaining normalized core shape
and the last runner's residue/offset.

THM-760 subsequently closes this displayed ray and, more generally, every primitive family in
which twelve speeds share a common divisor: lifting a core witness through the divisor produces a
sheet on which the one coprime exceptional runner is safe. The genuine scale residual therefore
has at least two exceptional residue classes, or has no twelve-speed common-factor core.

Status cautions: THM-724's addendum closes its genuine single-killer case, but THM-726 still relies
on an unproved global far-element monotonicity statement; THM-741 is explicitly `CLAIMED` with an
unfinished evidence checklist. Neither may be used as an unconditional assembly lemma.

## 1. Status by piece

| piece | status | artifact |
|---|---|---|
| non-covering | PROVED | THM-366 + LRC(≤13) citation |
| ≤3-far | FINITE-EXACT AS RECORDED; independent rerun/Lean transcription open | kps THM-738 + THM-758 split (klein-S309) |
| ≥4-far branch | OPEN modulo scale-normal classification | corrected THM-758 + HYP-6780 |
| sampled raw bands | FINITE-EXACT FOR THE STATED BANKS ONLY | mac-mini-S105; klein-S312 |
| (H)-bands, bottom cores | CLOSED (complete sweep) | THM-756 (4,032 pairs; AP/GW corners) |
| safe-peel | Parts A/B PROVED; irreducible-tiled Part C empirical | mac-mini THM-753 |
| aligned monotonicity | PROVED | mac-mini THM-751 |
| shadow tiles | PROVED | THM-748 (klein), THM-749, THM-754 clean-slot |
| named coherent/cluster families | PROVED at their stated scopes | THM-668/737/739/740 |
| 12-speed common-factor core + one coprime exception | PROVED, all scales | THM-760 |
| tail lanes | EXACT (identity + scans) | THM-750 closed budget; U1 discharged (S283) |

## 2. The Lean ledger

- **LRCClosedBudget.lean: 47 declarations, 0 sorries, all [propext, Classical.choice,
  Quot.sound]** — the (H)-edge chain machine-checked END TO END, geometric to spectral:
  pairOverlap → acorr_eq_model → geometric_disc_eq_discB → grid_deficit / raabe_B2 → discB;
  capped_envelope_kernel → Fourier envelopes → spectral_thm755. THM-731/732/755 have complete
  Lean faces; kps's exact-ℚ certificates have no prose links below the band edge.
- Fleet Lean assets: LRCDetunedDispatch (THM-668), LRCShadowGap (THM-748), the certificate
  supply (THM-693–698), SmallClusterFull, LRC13Citation, kps decide-bottoms.
- The formalized chain proves the stated geometric/Fourier/capped-envelope lemmas. It does **not**
  supply the missing uniform scale-normal classification, so there is no sound assembly theorem
  to wire yet.

## 3. Remaining work items

1. **Prove a scale-normal structure theorem.** Split normalized families into coherent dilation
   packs, additive/hierarchical clusters, and an incoherent residual, while retaining scale residue
   and killer offset. Raw far count and diameter are not invariants of this quotient.
2. **Attack the prime-grid bottleneck as a persistent translate-cover problem.** At level one,
   `I(13,p,1)` is exactly a 13-translate cover of
   `F_p^x/{±1}` by the strict-danger set. Classify which covers persist under lifts and evade the
   omit-one gcd reduction.
3. **Reproduce and formalize finite tiles.** Independently rerun THM-738's complete bank; finish
   rather than promote THM-741; attach machine-verifiable certificates and exact scope metadata.
4. **Hygiene.** Resolve theorem-ID collisions and require the status vocabulary `proved`,
   `finite-exact for stated bank`, `verified sample`, `conditional`, or `open`.

## 4. The one-line frontier

> **LRC(14) is open. The proved `f<=3` tile and the per-core capped envelope are substantial,
> but the `f>=4` branch requires a scale-normal structural classification, not raw enumeration.**

*Controlling corrections: HYP-6780 and the updated THM-758. Earlier S297/S310/S312 closure
language and the companion S297 reflection must be read through this correction.*
