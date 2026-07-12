---
source: opus-2026-07-11-S233
status: A RESOLUTION of the "tight-floor census" target + a clean two-bucket dispatch. The tight /
  AP-coherent families (the S232 "wall") are dispatched ELEMENTARILY by t = 1/14 (bucket A, no speed ÷14);
  the clean-ruler route (hB5) is responsible only for bucket B (some speed ÷14), where it is bounded,
  diameter-free, and adversarially robust. Verified 0 gap / 20000. The tight floor M = 1/14 lives entirely
  in the elementary bucket; the S232 wall dissolves.
tags:
  - lrc14
  - tight-floor
  - dispatch
  - clean-ruler
  - the-wall
  - THM-708
  - t-equals-one-over-14
---

# The tight-floor census is the elementary bucket

**opus-2026-07-11-S233.** The target from S232 was "the tight-floor census for the 1-dimensional-coherent
(AP-coherent) families." Working it resolves it: **that census is elementary**, and it clarifies the whole
dispatch into two clean buckets.

## The two-bucket dispatch

Every 13-family's loneliness certificate is one of two kinds:

- **Bucket A — no speed divisible by 14.** Then **`t = 1/14` is a loneliness witness**, so `M(S) ≥ 1/14`.
  *Proof (elementary):* `14 ∤ vᵢ ⟹ vᵢ/14 ∉ ℤ ⟹ ‖vᵢ/14‖ = min(vᵢ mod 14, 14 − vᵢ mod 14)/14 ≥ 1/14`. ∎
- **Bucket B — some speed divisible by 14.** Then `t = 1/14` fails (that runner sits at the origin), and
  loneliness is certified by a **clean ruler** (`B5 > 0`, kps THM-707).

**Verified:** `[clean ruler q ≤ 80] ∨ [no mult of 14]` covers **19999/19999** primitive families (0 gap).
The two mechanisms partition the certificate.

## The tight floor is entirely in bucket A (the census, resolved)

The tight / AP-coherent families — **AP `{1..13}`, GW `12→24`, the sporadic `V* = {1..11,13,24}`** — all have
**no multiple of 14**, so they are dispatched by `t = 1/14` (with equality, `M = 1/14`). This is mac-mini's
THM-708 ("tight families admit no clean ruler and need none; `14 ∤ every element`, sieve-dispatched"), now
seen as one half of the clean partition. So:

> **The tight-floor census is the elementary bucket A.** The extremizers of the whole problem (the `M = 1/14`
> families) are certified by a one-line witness `t = 1/14`. There is no hard "census" to run — the
> 1-dimensional-coherent families are exactly the ones that trivially clear the floor.

## The S232 wall dissolves

S232 found the AP-coherent families have **no clean ruler at any bounded modulus** — the "multiplicand-maximal
wall." That is not a gap: those families **are bucket A**, certified by `t = 1/14`. Verified: AP, GW, V* all
have `no-mult-14 = True` and `clean_ruler = None`. So the clean-ruler route's *failure* to certify the AP is
**correct behavior** — the AP is not the clean-ruler route's responsibility. The wall was an artifact of
asking the wrong mechanism to certify the wrong bucket.

## The clean-ruler route only owns bucket B — and there it is easy

The clean-ruler route (hB5) is responsible **only for bucket B** (some speed ÷14). There the mult-of-14 speed
**detunes AP-coherence**, so bucket B is comfortably lonely and clean-ruled at a **small, bounded,
diameter-free** modulus: adversarial hill-climb (maximising the smallest clean modulus within bucket B) tops
out at `q = 39–49` across Vmax ceilings 200 → 20000, with **0** no-clean families. So on its actual domain the
clean-ruler route never approaches the floor `1/14` — the tightest bucket-B family found sits at `M ≈ 0.12`,
comfortably above `1/14`.

## Honest limit

The verification (0 gap; bucket B always clean-ruled) is over **random/adversarial** families, **none of which
are covering**. The residual hard core is the **near-covering bucket-B families** (`M` near `1/14` *with* a
multiple of 14) — structured, measure-zero, not reached by hill-climb (so `M ≈ 0.12` is weak-search, not a
bound). "Every bucket-B family is clean-ruled" **is** hB5 = the content of LRC(14); this dispatch does not
prove it. What it *does* deliver:

1. The tight floor and the AP-coherent extremizers are **elementary** (bucket A, `t = 1/14`) — the
   "highest-value census target" was already trivial once framed by divisibility-by-14.
2. The clean-ruler route's domain is **exactly** bucket B (mult of 14), and there the ruler is bounded and
   diameter-free — the cleanest possible form of the remaining obligation.
3. The S232 anti-concentration "wall" is **resolved**, not a gap: the hardness of LRC(14) is *not* the AP; it
   is the hypothetical near-covering mult-of-14 family that the clean ruler must exclude.

This is the honest terminus of the S230→S233 arc: the bounded-modulus clean-ruler route is complete *modulo*
the mult-of-14 near-covering exclusion, and everything else (the tight floor, the AP, the summand-shell wall)
is elementary or dispatched.

→ THM-708 (mac-mini — tight = `t=1/14`-dispatched), THM-707 (kps — the clean ruler), THM-610 (covering hides at
`q* ≡ 0 mod 14`), THM-366 (small-denominator divisibility sieve), opus-S232 (the summand-shell wall this
places in bucket A), opus-S230/S231 (the bounded-modulus route). Files:
`lrc14_two_bucket_dispatch_opus_S233.py` (+`.out`).
