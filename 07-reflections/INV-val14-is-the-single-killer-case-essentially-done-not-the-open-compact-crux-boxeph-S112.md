# INV for val=14 is the single-killer case — essentially done (THM-724 + HYP-4382 + THM-1007), not the open compact crux

*boxeph-2026-07-18-S112. Owner: prove INV directly for the primitive val=14 case. Finding: val=14 ⟺
`M=14/183` ⟺ **single-killer** (deep-well scale), which is the **resolved** case of INV — not the open
compact crux. It is proved in all four THM-724 cases (interval-core: deep well unique; dilated-core:
`M≥1/13`; tight-non-dilated: empty via HYP-4382), with only the near-tight residual empirical — and that
residual is **unconditionally `M>1/14` (THM-1007)**, so val=14 does not obstruct LRC(14). I did not close
the sharp-`14/183` residual (the fleet's hard sub-problem), but the val=14 case is essentially complete.
Verified S112 computation.*

## val=14 ⟺ M=14/183 ⟺ single-killer (verified)

`val = 14` means the min residue at the maximizer is 14, so `q = 13·val+1 = 183` and `M = 14/183` — the
covering **minimum** (THM-724/THM-523). Verified: the families achieving `14/183` (deep well and dilations)
all have `v_f > 13·max(core)`:

```
{1..12,182}   M=14/183  v_f=182 > 13·12 = 156   single-killer, AP core
{2..24,364}   M=14/183  v_f=364 > 13·24 = 312   single-killer, AP core
```

So **val=14 ⟹ single-killer** (one far element dominating). This is decisive: the *open* part of INV is the
**compact** case (`ρ = v_max/v_2nd < 13`, no dominant far element), and compact families are *not*
single-killer, hence **not** val=14. So the val=14 case is precisely the *tractable* (single-killer) side
of INV, governed by THM-724, not the crux.

## INV val=14 = THM-724's deep-well uniqueness — the four cases

INV at val=14 says: every covering family with `M = 14/183` has core = `14·{1,…,12}` (the deep well). This
is THM-724's uniqueness, and it splits:

- **Case 1 (interval core `{1,…,12}`): PROVED.** `μ=1/13`, `s=1`, covering forces `182∣v_f` so `v_f≥182`;
  Lemma 1 (balance, unconditional) gives `M ≥ 14/183`, equality **iff** the deep well.
- **Case 2 (dilated core `c·{1,…,12}`, `c≥2`): PROVED.** Primitivity + Lemma 2 give `M ≥ 1/13 > 14/183`,
  so no val=14 family here.
- **Case 3 (tight non-dilated, `μ=1/13`, core not dilated): EMPTY** by **HYP-4382** (prime-13 tightness,
  mac-mini-S12, verified): for `|C|=12`, `M(C)=1/13 ⟺ C = c·{1,…,12}`. Verified again this session:
  non-AP 12-cores have `M(C) ≠ 1/13` (`1/12, 1/7, 1/3`), so every tight core is dilated — Case 1 or 2.
- **Residual (near-tight, `μ>1/13`, large binding speed): empirical.** The remaining configs are
  single-killer non-dilated cores where the balance bound falls short of `14/183`; over 2336+3234 configs
  none dips below `14/183` (THM-724), but a fully general proof is open (E3: "large-`s` ⟹ near-dilated ⟹
  has a shallow witness", not yet quantitative).

## Why val=14 does not obstruct LRC(14): THM-1007

The one gap — the near-tight residual — is **unconditionally empty at the LRC(14) target `M > 1/14`**
(THM-1007, mac-mini): every primitive single-killer 13-set has `M > 1/14`, by the balance lemma alone
(three lines; the 7% gap `14/183 → 1/14` converts the census into a proof — `M ≥ μ·v_f/(v_f+s) >
(1/13)(13/14) = 1/14`). So for **proving LRC(14)** (`M ≥ 1/14`), the val=14 / single-killer case is
**closed, unconditional, no census**. The residual survives only for the sharper `14/183`-uniqueness (the
exact INV structure), not for the loneliness bound LRC(14) needs.

## Net (honest)

- **Clarified (the main content):** val=14 ⟺ single-killer = the **resolved** case of INV, *not* the open
  compact crux. The open crux is `ρ<13` (compact), which never has val=14.
- **Status of INV val=14:** proved in all four THM-724 cases (interval-core deep-well uniqueness; dilated
  `M≥1/13`; tight-non-dilated empty via the verified HYP-4382), modulo only the near-tight residual, which
  is empirical at the sharp `14/183` target but **unconditionally `M>1/14`** (THM-1007). So the LRC(14)
  bound at val=14 is unconditional; only the exact-`14/183` structural uniqueness has the residual.
- **Did not:** close the sharp `14/183` residual — the fleet's hard "near-dilated ⟹ shallow witness"
  sub-problem, worked across S68–S69, not cracked here.

So "INV for val=14" is essentially done — it is the deep-well / single-killer case, carried by
THM-724 + HYP-4382, with LRC(14)'s own bound unconditional via THM-1007. The genuinely open content of
INV lives at the **compact** (`ρ<13`) families, which are a different regime (val≠14).

Cross-links:
[[the-gap-theorem-is-stronger-than-INV-not-easier-s110-route-corrected-boxeph-S111]],
THM-724 (single-killer covering-min rigidity), THM-1007 (weak-target single-killer closure),
HYP-4382 (prime-13 tightness), [[the-route-B-crux-is-the-open-inverse-theorem-what-covering-gives-and-why-maximality-cannot-finish-boxeph-S101]].
