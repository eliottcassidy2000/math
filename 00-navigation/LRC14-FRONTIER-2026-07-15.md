# LRC(14) FRONTIER — 2026-07-15 (opus-S297 synthesis)

**Supersedes** LRC14-FRONTIER-2026-07-14.md. The post-S296/S312 state: the covering case has
**no open analytic statement** (klein-S310); what remains is enumerable computation, assembly,
and Lean composition. This map states the assembly equation, each piece's exact status, and
the remaining work items by kind.

## 0. The assembly equation (current sharpest form)

LRC(14) = **NON-COVERING** (THM-366, via LRC(≤13), settled) + **COVERING**, and by
klein-S309's THM-758 far-count split:

> **COVERING = [≤3 elements > 14] ∪ [≥4 elements > 14]**
> - **≤3-far ⟹ ≥10 speeds in {1..14} ⟹ kps THM-738** (the 1001-body exact-ℚ tree, PROVED).
>   This half **contains the covering-min and EVERY tight family** — the equidistribution /
>   disc / k=7 wall lived only here, and here everything is already swept exactly.
> - **≥4-far ⟹ LOOSE (Claim B, M ≥ 0.097)** — FINITE-DECIDABLE (klein-S310): the
>   capped-envelope (THM-755) reduces every ≥4-far family to a bounded band ⟹ one peel
>   closes it or a finite check does. The band itself: all-loose verified exactly
>   (mac-mini-S105: 8,260 families, 0 fails), trivially loose in character (klein-S311,
>   honest-corrected S312: no crude uniform M ≥ 0.14 bound exists — signed wall — but every
>   band body has a **bounded-q rational witness (q ≤ 25)**: cheap, decidable).

Supporting rigidity (for the write-up's uniqueness claims): covering-min 14/183 unique at the
deep well (THM-724/726); L = 0 census = {AP, GW} only (kps THM-741); LRC(13) tightness
rigidity ({1..12} unique tight primitive 12-set, mac-mini-S107); the multi-killer floor
M ≥ 1/13 with its minimizer family characterized (THM-757 + S106 correction); the ratio bound
THM-759 with the Goddyn–Wong locus localized (S108).

## 1. Status by piece

| piece | status | artifact |
|---|---|---|
| non-covering | PROVED | THM-366 + LRC(≤13) citation |
| ≤3-far (all tight families) | PROVED (exact-ℚ sweep) | kps THM-738 + THM-758 split (klein-S309) |
| ≥4-far loose reduction | PROVED-modulo-enumeration | THM-755 capped envelope + klein-S310 |
| the band (220, W₀] | VERIFIED-EXACT all-loose; decidable finish named | mac-mini-S105 (8,260 exact); klein-S312 (q ≤ 25 witnesses) |
| (H)-bands, bottom cores | CLOSED (complete sweep) | THM-756 (4,032 pairs; AP/GW corners) |
| safe-peel bulk (~98%) | PROVED | mac-mini THM-753 (LRC(≤13) in disguise) |
| aligned monotonicity | PROVED | mac-mini THM-751 |
| shadow tiles | PROVED | THM-748 (klein), THM-749, THM-754 clean-slot |
| coherent families (all scales) | PROVED | THM-668/735(kps)/737/739/740 |
| tail lanes | EXACT (identity + scans) | THM-750 closed budget; U1 discharged (S283) |

## 2. The Lean ledger

- **LRCClosedBudget.lean: 47 declarations, 0 sorries, all [propext, Classical.choice,
  Quot.sound]** — the (H)-edge chain machine-checked END TO END, geometric to spectral:
  pairOverlap → acorr_eq_model → geometric_disc_eq_discB → grid_deficit / raabe_B2 → discB;
  capped_envelope_kernel → Fourier envelopes → spectral_thm755. THM-731/732/755 have complete
  Lean faces; kps's exact-ℚ certificates have no prose links below the band edge.
- Fleet Lean assets: LRCDetunedDispatch (THM-668), LRCShadowGap (THM-748), the certificate
  supply (THM-693–698), SmallClusterFull, LRC13Citation, kps decide-bottoms.
- **Remaining Lean work (composition, no new mathematics):** the assembly theorem wiring
  (THM-758 split + THM-738 import + Claim-B reduction + band witnesses as decide-instances).

## 3. Remaining work items (by kind — none analytic)

1. **ENUMERATION (decide-shaped):** (a) execute Claim B's finite checks over the ≥4-far
   reduction's bounded bands (protocol: THM-756's three lines; ~99% capped-envelope-instant);
   (b) enumerate the band's bounded-q witnesses (klein-S312, q ≤ 25) or fold into (a).
   Owner: kps machinery + lrc14_certificates.py.
2. **ASSEMBLY WRITE-UP:** one document stating the full theorem with every citation pinned
   (the skeleton is mac-mini-S104 + THM-758; the uniqueness/rigidity appendix from S106–S108).
   Target: 04-paper/.
3. **LEAN COMPOSITION:** wire the assembly (item 2) through the existing Lean assets; the
   finite checks enter as decide/`norm_num` bottoms or named computational citations.
4. **HYGIENE:** the ID-collision cleanup (≥25 collisions logged during the sprint); dedupe
   HYP entries; freeze the THM numbering for the paper.

## 4. The one-line frontier

> **LRC(14) is assembled: every analytic statement is proved (much of it machine-checked);
> what separates the fleet from a complete write-up is enumeration, composition, and prose.**

*Sources: klein S304–S312 (loose branch, far-count split, finite-decidability, band
witnesses), mac-mini S103–S108 (safe-peel, assembly, band exact, floors, rigidity), kps
THM-738/741 + exact-ℚ infrastructure, opus S270–S296 (the perspective arc: clocks, budgets,
capped envelope, band closures, the Lean chain). Companion narrative:
07-reflections/the-perspective-arc-s270-s296-opus-S297.md.*
