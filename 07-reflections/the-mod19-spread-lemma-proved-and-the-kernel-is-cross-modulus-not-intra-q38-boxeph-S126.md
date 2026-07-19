# The mod-19 antipodal-spread lemma (PROVED), and why the "second interior aligned gap" kernel is cross-modulus, not intra-q=38

*boxeph-2026-07-19-S126. Owner: prove the isoperimetric spread bound on the 1/19 alphabet; prove the
kernel (forbid a second interior aligned gap). Result: (1) the spread bound is a genuine PROVED theorem —
`M(C) < 2/19 ⟹ the residues mod 19 antipodally cover ℤ/19` (Lemma 1), a translation-sensitive necessary
condition on the 1/19 alphabet; (2) the kernel, examined honestly, does NOT close 3/38 as an intra-q=38
statement: the q=38 conditions (band + covering + mod-19 spread) are FEASIBLE, and the witnesses that beat
the q=38 hole are CROSS-MODULUS (the S124 medium needles at q=24, 29, …). So the "no second interior
aligned gap" kernel is provably about the origin gap at q=38 only, and the residual obstruction is
irreducibly cross-modulus — the analytic (unbounded-modulus) core. Verified S126.*

## Lemma 1 — the isoperimetric spread bound on the 1/19 alphabet (PROVED)

> **Lemma 1.** Let `C` be a 12-set of distinct positive integers with `M(C) < 2/19`. Then either `19 ∣ v`
> for some `v ∈ C`, or the residues `{v mod 19 : v ∈ C}` cover **all 9 antipodal unit-pairs**
> `{±1, ±2, …, ±9}` of `ℤ/19`.

*Proof.* For `n = 1, …, 18` consider the witness `t = n/19`. For each speed `v`,
`‖v·n/19‖ = ‖(vn mod 19)/19‖` is a multiple of `1/19`, hence so is `min_v ‖v·n/19‖`. Since
`M(C) = max_t min_v ‖vt‖ ≥ min_v ‖v·n/19‖` and `M(C) < 2/19`, this minimum is `< 2/19`, and being a
multiple of `1/19` it is `≤ 1/19`. So for every `n` there is a speed `v` with `‖v·n/19‖ ≤ 1/19`, i.e.
`vn ≡ 0, 1,` or `−1 (mod 19)`. Now suppose `19 ∤ v` for all `v`. Since `n ∈ {1,…,18}` gives `19 ∤ n`, the
case `vn ≡ 0` is impossible, so `vn ≡ ±1`, i.e. `v ≡ ±n^{-1} (mod 19)`. As `n` ranges over the units
`(ℤ/19)^*`, so does `n^{-1}`, so for every unit `u` some speed satisfies `v ≡ ±u (mod 19)`. Hence
`{±(v mod 19)} ⊇ (ℤ/19)^*` — all 9 antipodal pairs. ∎

**This is the "isoperimetric spread bound":** a family that is tight enough to enter the gap regime
(`M < 2/19`, which covers `3/38`, `2/25`, and the whole gap) must have its residues **maximally spread**
across `ℤ/19` — they antipodally cover the entire multiplicative group; they cannot be concentrated. It is
**translation-sensitive** (it is a statement about residues relative to the origin `0 (mod 19)`), so it is
on the correct side of opus's triage (translation-invariant invariants cannot see the LRC floor).

*Verified* on every genuine small-`M` family: `{1,…,12}` (`M=1/13`), `{1,…,11,24}` (`2/25`), `{1,…,11,36}`
(`3/37`), `{1,…,11,48}` (`4/49`), `2·{1,…,12}` (`1/13`), `{1,…,10,11,13}` (`1/12`) — all antipodally
spread mod 19, `0` violations.

**Lean (kernel-pure, added S127).** `LRCMod19Spread.lean` (`namespace LonelyRunner`), the direct mod-19
analogue of `LRCMod13Blocking`, all `[propext, (Classical.choice,) Quot.sound]`, no `sorry`:
- `mod19_middle_far` — integer core: `r∈[2,17] ⟹ 2 ≤ |19k+r|`.
- `sieve19_single` / `sieve19_middle_witness` — if `(v·b) mod 19 ∈ [2,17]` then `t=b/19` gives
  `‖v·(b/19)‖ ≥ 2/19`; family form gives `M ≥ 2/19`.
- `no_middle_band_of_close` — contrapositive: a `<2/19`-close runner empties the middle band.
- `antipodal_spread` — the lemma: `¬(19∣c_i)` for all `i`, plus a `<2/19`-close runner at every scale `b`
  (which holds when `M<2/19`), forces, at every `b` with `19∤b`, some runner with residue `±1` mod 19 —
  the per-scale form of the antipodal covering (via `b ↦ b⁻¹` on the units).

## The kernel — "forbid a second interior aligned gap" — is cross-modulus, not intra-q=38

The natural kernel would close `3/38` at the `q=38` maximizer if a second `q=38`-aligned deep gap could be
forbidden. Examining it honestly:

**Structural fact (proved).** At the `s=38` maximizer `t* = m/38`, any pair `(w_1, w_2)` with `w_1+w_2 = 38`
lands at **antipodal positions about the origin**: `(w_1+w_2)t* = 38·(m/38) = m ∈ ℤ`, so `w_1 t* ≡ −w_2 t*`.
Hence **every sum-38 pair straddles the same origin gap** — there is no *second* origin gap that a sum-38
(determinant-carrying) pair could open. The origin gap is the unique `q=38`-aligned deep hole.

**But this does not close 3/38.** A family fails `M=3/38` not by a second `q=38` gap but by a **deeper hole
at a different modulus**. Concretely, the band-filled covering families with the `(3,35)` pair have their
`3/38` hole at `t=1/38`, yet their true maximizer sits at a **cross-modulus** point:
`{3,5,7,8,9,10,11,12,13,15,21,35}` peaks at `t=1/24` (`M=1/8`), `{…,17,21,24,35}` at `t=2/29` (`M=5/29`).
The witnesses that beat `3/38` have denominators `24, 29, 21, 23` — **not 38**. This is exactly the S124
medium-needle covering.

**And the `q=38` conditions are feasible.** The intra-`q=38` necessary conditions — band `⊂[3,35]` (mod 38),
covering, and the mod-19 spread (Lemma 1) — are jointly **satisfiable**: of `990` band-filled covering
families, `271` also satisfy the mod-19 spread. None reaches `M=3/38`, but not because the `q=38`
constraints clash — they don't. The obstruction lives entirely in the **other** moduli.

> **Conclusion (the honest kernel).** The `q=38`-intra-modulus analysis is *complete but insufficient*:
> Lemma 1 gives a real necessary condition (antipodal spread on the 1/19 alphabet), and the origin gap is
> the unique `q=38`-aligned deep hole (no second intra-`q=38` gap). But `3/38` is defeated **cross-modulus**
> — by a deeper hole at some `q' ∉ {38-divisors}`. There is no single-modulus kernel that forbids it; the
> residual is the S124 adaptive needle-covering, whose modulus is unbounded (the escape tail). That is
> precisely why `3/38` sits in the analytic core.

## What this settles, and what it does not

- **Settled (proved):** the mod-19 antipodal-spread bound (Lemma 1) — a genuine necessary condition for any
  gap-regime family, translation-sensitive; and the uniqueness of the `q=38`-aligned origin gap (no second
  sum-38 gap).
- **Settled (negative, honest):** a `q=38`-intra-modulus kernel **cannot** close `3/38` — the `q=38`
  conditions are feasible; the defeat is cross-modulus. This sharpens the S124/S125 picture: the obstruction
  is not a single missing constraint at `38` but the *joint* (adaptive, unbounded-modulus) needle-covering.
- **Not proved:** `3/38` unachievability. It requires the cross-modulus argument (all `q'` at once, up the
  escape tail) — the analytic core, unchanged. Lemma 1 is one more proved necessary condition on the way in.

Cross-links:
[[the-rank-identity-at-q38-cover-debt-must-spread-and-3-38-wastes-it-boxeph-S125]],
[[the-kakeya-needle-obstruction-to-3-over-38-medium-modulus-needles-cover-the-band-boxeph-S124]],
[[the-determinant-stratified-gap-numerator-two-is-excluded-and-3-38-is-the-depth-minimal-target-boxeph-S123]],
opus THM-1185/1220 (translation-invariance triage), macmini-S27 (mod-19), kps-S12 (gap empty [1,26]),
`lrc14_mod19_spread_kernel_boxeph_S126.py`.
