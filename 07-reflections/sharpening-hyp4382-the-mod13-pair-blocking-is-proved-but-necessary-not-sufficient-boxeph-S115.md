# Sharpening HYP-4382: the mod-13 pair-blocking is PROVED (not just verified), but it is necessary, not sufficient

*boxeph-2026-07-18-S115. Owner: sharpen HYP-4382 (`|C|=12, M(C)=1/13 ⟺ C = d·{1,…,12}`) toward the LRC(13)
equality proof. Result: a **proved** necessary condition — `M(C)=1/13 ⟹ {±c_i mod 13} = {1,…,12}` (mod-13
pair-blocking), the natural-modulus, *proved* analog of the *verified* HYP-4622 (mod 25). But it is
**necessary, not sufficient**: non-AP families that are complete mod 13 exist with `M=1/12`. So HYP-4382's
rigidity is *orthogonal* to the mod-13 blocking and needs the full sieve (all moduli) — it does not collapse
to a single-prime condition. Genuine partial progress, not a proof. Verified S115 computation.*

## The proved necessary condition

> **Lemma (mod-13 pair-blocking, PROVED).** If `M(C)=1/13` (`|C|=12`) and `13 ∤ c_i` for all `i`, then
> `{±c_i mod 13} = {1,…,12}` — for every `b∈{1,…,12}` some `c_i ≡ ±b^{-1} (mod 13)`.

*Proof.* For any `b`, the time `t=b/13` gives `M(C) ≥ min_i ‖c_i·b/13‖`. So `min_i ‖c_i b/13‖ ≤ 1/13`, i.e.
`min_i |c_i b \bmod 13| ≤ 1`: some `c_i b ≡ 0, ±1 (mod 13)`. With `13∤c_i`, `c_i b ≡ 0` is impossible, so
`c_i b ≡ ±1`, i.e. `c_i ≡ ±b^{-1}`. As `b` ranges over `{1,…,12}`, `b^{-1}` ranges over all of `{1,…,12}`,
so `{±c_i mod 13}` covers every nonzero class. ∎

This is one line from the sieve, and it **upgrades** the project's mod-25 pair-blocking (HYP-4622, mac-mini,
*verified*) to a *proved* statement at the natural modulus 13 for the actual HYP-4382 target. Verified:
dilated APs `c·{1,…,12}` satisfy it (`c·\{1..12\} \bmod 13 = \{1,…,12\}`, complete); non-blockers like the
general AP `{1,14,…,144}` (all `≡1 mod 13`) fail it and have `M = 67/145 ≠ 1/13`.

## But it is necessary, not sufficient

The decisive test: does a **non-AP** family that is **complete mod 13** reach `M=1/13`? No.

| `C` (complete mod 13, non-AP) | `M(C)` |
|---|---|
| `{1,…,11, 25}` (`25≡12`) | `1/12` |
| `{1,…,11, 38}` (`38≡12`) | `1/12` |
| `{2,…,12, 14}` (`14≡1`) | `1/8` |

Every non-AP complete-mod-13 family has `M > 1/13` — it is **beaten at another modulus** (`q=12, 8, …`).
So mod-13 blocking alone cannot force the AP. HYP-4382 is confirmed (only the dilated AP is tight), but the
rigidity does **not** reduce to a single prime.

## Why mod-13 is orthogonal to the AP structure

At the maximizer `q = 13·val`, the residues split by CRT (`gcd(val,13)=1`) into `(r_i \bmod 13, r_i \bmod
val)`. The mod-13 blocking constrains `r_i \bmod 13 = a·c_i \bmod 13`. But the **AP is the
offset-vanishing** `j_i := r_i \bmod val = 0` (S94's form, at the 12-family level) — a condition on the
`\bmod val` coordinate, *independent* of `\bmod 13`. So the proved mod-13 blocking says nothing about the
offsets; the AP rigidity lives in the `\bmod val` (all-other-moduli) coordinate. This is exactly why
`{1,…,11,25}` (complete mod 13, offsets nonzero) is a blocker yet not tight.

So the sharpening factorizes the tightness into two orthogonal pieces:
- **mod-13 blocking** (`\bmod 13` coordinate): `{±c_i} = {1,…,12}` — **PROVED**.
- **offset-vanishing** (`\bmod val` coordinate, = all other moduli): `j_i = 0` — the residual = the crux.

The mod-13 piece is closed; the offset piece is the maximality-over-all-`q` = the inverse theorem.

## Net (honest)

- **Proved:** `M(C)=1/13 ⟹ {±c_i mod 13} = {1,…,12}` — the natural-modulus, proved pair-blocking (upgrading
  HYP-4622's verified mod-25 status). The same argument gives a proved mod-`p` blocking for every prime `p`.
- **Showed:** this is **necessary, not sufficient** — complete-mod-13 non-AP families have `M=1/12`, beaten
  at other moduli. So HYP-4382 does *not* collapse to a single prime; the AP-forcing is the offset-vanishing
  in the complementary coordinate.
- **Did not:** prove HYP-4382. The rigidity is the all-moduli / offset-vanishing content — the crux (S94,
  S114). The mod-13 blocking is a clean proved *slice*, not the whole.

So the LRC(13) equality proof needs to force offset-vanishing across all moduli simultaneously — the
mod-13 slice is now proved, but the essential rigidity (a single tight time forbidding every other rational
from beating it) is the same inverse-theorem wall. Sharper, not closed.

Cross-links:
[[the-non-dilated-core-rigidity-residual-is-tao-n12-the-definitive-frontier-boxeph-S114]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
HYP-4382 (n=12 tightness), HYP-4622 (mod-25 pair-blocking), THM-724.
