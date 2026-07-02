# THM-607: The certificate-box theorem (the exact V-region of fixed ladder data, and the region/arithmetic factorization)

**Status:** PROVED (one-line consequences of THM-606 + the slack computation; exact-rational
verification: 400 random sharp-region tuples, 0 failures — `lrc14_certificate_box_theorem_opus_20260702_S37.py` + `.out`)
**Author:** opus-2026-07-02-S37 (HYP-3906)
**Studies:** the V-independence of THM-606 certificates, composed with the Module-0
(`RatIntervals`, kps HYP-3961) view; the multi-cluster half of the shape universe.

---

## The factorization principle

A ladder certificate datum is `(window [lo,hi]; offsets offs_ℓ; targets c_ℓ)` — **regions on the
τ-line**. The fast speeds `V_ℓ` enter the ladder ONLY through the ruler condition (H-len) and the
drift budget — **two rational inequalities per level**. So the certified set factors:

    REGIONS carry the offset data   (V-free; interval calculus: clip/translate/scale/union —
                                     exactly Module 0's decidable layer);
    ARITHMETIC carries the V data   (two inequalities per level; pure `decide`).

This is *why* one datum certifies infinitely many configurations, and why the Lean consumption
shape is a `LadderBox` structure (schedule + certs + slacks) whose family theorem quantifies over
the integer points of a rational box — the playbook's certificate-carrying discipline, verbatim.

## (A) Intrinsic slack

For a certificate `c` of offsets `offs` on `[lo, hi]`, the **intrinsic slack** is
`μ̄ := min_{o ∈ offs} [ dist of the arc (o·lo − c, o·hi − c) to its cell boundary ] − h` —
a rational computed on the region side alone. `arcSafe(h + μ, o, c)` holds iff `μ ≤ μ̄`.
S36 instance: `μ̄ = 793/28000, 43/1750, 1/14000` (each sharper than the schedule-μ used there;
level 3's tiny slack `1/14000` is why the S36 margin was razor-thin: `143/2000 − 1/14`).

## (B) The fixed-schedule box

For fixed schedule `δ_1 > … > δ_d = 0`, THM-606 certifies **exactly**

    V_ℓ ∈ ( 1/(δ_{ℓ−1} − δ_ℓ) , μ̄_ℓ/δ_ℓ ]  for ℓ < d,   V_d ∈ ( 1/δ_{d−1} , ∞ ).

*Proof:* (H-len) is the lower endpoint; the drift budget `V_ℓ δ_ℓ ≤ μ̄_ℓ` is the upper; (A) says
these are all that is needed. ∎  S36 instance: `V₁ ∈ {41..50} × V₂ ∈ {1819..2010} × V₃ > 81818`
— 10 × 192 × ∞ integer tuples from one datum.

## (C) The sharp per-tuple region

Choosing the schedule per tuple (`δ_ℓ := μ̄_ℓ/V_ℓ`), the datum certifies **every** tuple with

    V_1 > (1 + μ̄_1)/δ_0    and    V_ℓ/V_{ℓ−1} > (1 + μ̄_ℓ)/μ̄_{ℓ−1}   (ℓ = 2..d).

*Proof:* saturate the budget; (H-len) becomes the ratio inequality. ∎  These are the fixed
datum's SHARP multiplicative thresholds — S36 instance: `V₁ > 41.13`, `V₂/V₁ > 36.18`,
`V₃/V₂ > 40.70`; verified end-to-end on 400 random tuples (0 failures, worst margin `143/2000`).

## (D) Mesh tiling (decidable form of (C))

(C)'s schedules are tuple-dependent (not finite data). The geometric mesh `δ_ℓ^{(j)} = μ̄_ℓ·2^{−j}`
restores finiteness: each box has multiplicative width ≥ 2 per bounded axis, so the mesh covers
every tuple with ratios `≥ 2(1+μ̄_ℓ)/μ̄_{ℓ−1}` (S36: ratios ≥ 73, 82). Covering `V ≤ 10⁹` on the
`d−1` bounded axes needs `≤ (log₂ 10⁹)^{d−1} ≈ 900` boxes at `d = 3` — **one certificate datum
plus a 900-row ledger of pure rational inequalities certifies every sharp-ratio tuple to 10⁹**,
with no new region computation per row. The gap between (C)'s sharp ratios and (D)'s mesh ratios
(factor `r = 2`) is the price of finite data; any `r > 1` works, trading rows for threshold.

## Consequences for the shape universe (F5 refined)

The multi-cluster leg of the universe splits as: cluster ratios above the mesh threshold —
covered by finitely many `LadderBox` rows per shape; ratios below — the merge/split finite cases
(F0/F4 of THM-606). Both halves are decide-shaped. The region/arithmetic factorization means the
swarm's Module-0 interval engine and the box ledgers never block each other: semantic region
lemmas (`mem_inter` invariant-free) discharge the certificate half; kernel `decide` on two
inequalities per level discharges the V half.

-> THM-606, HYP-3905, HYP-3961 (Module 0), HYP-3864 (playbook), HYP-3960 (star-census), OPEN-Q-108.
