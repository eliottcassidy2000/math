# THM-979 — The two-circle deep certificate, Part I (death-star-2026-07-17-S52)

**Status:** Part I PROVED (Lean, kernel-pure — `TournamentH7/LRCTwoCircle.lean`,
standard trio ×8); the full iff is recon-EXACT and Part II (⟹) is fully
structured. Source: HYP-7265 (the "wagner-circle certificate" hint decoded).

## The certificate (recon: 1185 moduli, every multiplier, ZERO mismatches)

On the canonical family v = (1,…,13):
**deep (bandCount ≥ 6) ⟺ the multiplier lies on one of TWO resonance circles**
- (I) the integer circle: `84p < q ∨ 84(q−p) < q`
- (II) the half circle: `84|2p − q| < q` (the mirror-invariant circle).

Exact counts follow: `#I = 2⌊(q−1)/84⌋`; `#II` closed-form by parity
(recon-verified); overlap only at tiny q. CoverageCapped and the census race
on the canonical family become closed-form.

## Part I (Lean ×8, this session)

`divisor_descent` (failures descend to the reduced ray — the witness-range
workhorse); `circleI_low_fails`/`circleI_high_fails`/`circleII_fails` (the
explicit witnessed failures); `bandCount_ge_of_card` + `not_inBand_of_witness`
(assembly tools); **`circleI_deep` + `circleII_deep`** (the full ⟸: circles ⟹
six failures ⟹ deep).

## Part II (⟹, structured; named next)

Case analysis on the smallest failing speed k₀: k₀=1 nests via the hub lock
(largest of six failing speeds ≥ 6 gives the 84-width); k₀=2 forces the six
evens (parity locks kill odd ≤ 11; the 13-branch ray estimate
`13|2p−q| > 6q/7` kills 13-coexistence; 12 ∈ F delivers circle II); k₀ ∈
{3,4,5,6,8} die by parity/divisor lock contradictions (each caps |F| ≤ 5);
k₀ = 7 dies by the residue gem: the branch bound forces `l·w₇ ≡ ±1 (mod 7)`
for every failing l — at most 2 residue classes among {8,…,13}, so |F| ≤ 3;
k₀ ≥ 9 is trivial (≤ 5 speeds remain). Every step is an instance of the
THM-967/972 lock machinery.

## Why "circles"

The deep set is a union of Bohr-type circles indexed by resonance denominators
m with ≥ 6 multiples ≤ 13 — only m = 1, 2 qualify. Circle II is invariant
under the mirror p ↦ q−p with the half-twist (the Möbius/Wagner echo).
