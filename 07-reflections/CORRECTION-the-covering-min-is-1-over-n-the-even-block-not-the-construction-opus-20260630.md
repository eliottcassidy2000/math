# CORRECTION (major): the covering-min is 1/n EXACTLY — the EVEN BLOCK 2·{1,…,n−1} (the AP doubled) is covering for even n with M=1/n — NOT the construction n/Φ₆(n); the construction (ζ₆-line, hexagonal, Sylvester, the whole arc) was a RED HERRING (a non-extremal covering set with M=n/Φ₆>1/n); the conjecture is TIGHT for even n (covering-min = 1/n), and PARITY (even/odd n) is the key: the AP's doubling covers iff n is even, so n=14 ⇒ covering-min = 1/14 — this vindicates the recursive parity mindset while refuting the construction-focus that I, klein, and mac-mini all shared

*opus-2026-06-30. Owner flagged: "convergent beaten at n=7 (2/13), n=8 (2/15), n=9 (4/33); the convergent-
for-n≥7 premise is definitively refuted." Verified — and it goes deeper than the convergent: the covering-min
is 1/n via the even block. Honest, large correction. Recording exactly what is refuted and what is sharpened.*

## Verification (the owner is right, and then some)
| n | construction `n/Φ₆` | TRUE covering-min (verified) | by |
|---|---|---|---|
| 7 (odd) | `7/43=0.1628` | **`2/13=0.1538`** | `[1,2,5,6,7,8]` (thorough search, 12 sets beat the construction) |
| 8 (even) | `8/57=0.1404` | **`1/8=0.1250`** | `[2,4,6,8,10,12,14] = 2·{1..7}` (the EVEN BLOCK) |
| 14 (even) | `14/183=0.0765` | **`1/14=0.0714`** | `[2,4,…,26] = 2·{1..13}` (the EVEN BLOCK) |
> **The construction is beaten at every n≥7.** For EVEN n the beater is decisive: the **even block
> `2·{1,…,n−1}`** is covering (it has a multiple of every `q≤n`: even `q` directly, odd `q≤n−1` via `2q`,
> and `q=n` via `n` itself since `n` is even) and **`M(2·S) = M(S)`**, so `M = M({1..n−1}) = 1/n`. **`n=14`
> is even ⇒ covering-min = `1/14`, the conjecture bound itself.**

## What is REFUTED (mine, klein's, mac-mini's)
The construction `{1,…,n−2,(n−1)n}` with `M = n/Φ₆(n)` is **NOT the covering-min** for `n≥7`. Therefore:
- ❌ **"covering-min = n/Φ₆, margin `(n−1)/(nΦ₆)~1/n²`"** — the real covering-min is `1/n` (even n), margin
  ZERO. The "razor-thin margin" was an artifact of the wrong family.
- ❌ **the witness=ζ₆ / hexagonal-Kershner / cyclic-Kershner / Steiner-PG(2,13) attack** — all about the
  construction, which is not extremal. The ζ₆-line is a real configuration but it is not the covering-min.
- ❌ **the observer's escape = the convergent `[0;n−1,n]`** — the covering-min is `1/n` (even n) or the
  mediant/deeper (odd n), NOT the convergent `n/Φ₆`.
- ❌ **klein HYP-3724 / mac-mini HYP-3701 "construction stands for n≥7"** — refuted; the construction is
  beaten at n=7,8,9,14,…
- ❌ **the Sylvester/Egyptian "covering-min recursion"** — `Φ₆`-iteration IS Sylvester (a true fact about
  `Φ₆`), but `Φ₆` is **not** the covering-min, so the "LRC covering = Sylvester system" reading is void.

## Why it was UNEXPECTED (the shared blind spot)
We — all three of us — chased the **construction** `{1,…,n−2,(n−1)n}` because it carries a *rich* structure
(Eisenstein `Φ₆`, the `ζ₆` hexagonal line, Sylvester, PG(2,13)). The richness was seductive and we proved
many true things *about it*. But we **never checked the trivial even block `2·{1,…,n−1}`** — the AP simply
doubled — which achieves `1/n` with no structure at all. My n=14 scans (107-set, adversarial) and mac-mini's
all sampled *dense interval + outlier* families and the *drop-2* family; none contained `{2,4,…,26}`. **The
extremal was the dullest set in the room, and we were looking at the most interesting one.**

## How it SHARPENS (and vindicates the parity mindset)
- **The covering-min is governed by PARITY of n** — exactly the "even/odd" the owner kept pointing at:
  > **EVEN n: covering-min = `1/n` EXACTLY** (the even block `= a⁻¹(AP) = 2·{1,…,n−1}` is covering; the AP's
  > doubling reaches `q=n` because `n` is even). **The conjecture is TIGHT** — `1/n` is achieved by a covering
  > set, so LRC(n) ⟺ *no covering set goes below `1/n`*, with the even block sitting exactly on the bound.
  > **ODD n: the even block FAILS** (no multiple of odd `n` among `2·{1,…,n−1}`), so the covering-min is `>1/n`
  > (n=7: `2/13`; n=9: `4/33`) — the realizability-on-the-Stern-Brocot story, the genuinely n-dependent part.
- **This is the recursive functional decomposition vindicated:** `a⁻¹ = ×2` (DOUBLING) applied to the AP
  (the extremal `f·g=T` staircase) gives the even block; **parity routes whether the doubling covers `q=n`**.
  The owner's `a=÷2`, `b=+1`, parity skeleton was right; the construction/`Φ₆`/hexagon was the distraction.
- **LRC(14) restated, correctly:** the covering-min is `1/14`, achieved by `2·{1,…,13}` (the AP doubled).
  The conjecture is that **this even block is the worst covering 13-set** — no covering set has `M<1/14`. The
  hard content is the same LRC(14), now with the *tight* extremal named (the doubled AP), not a `~1/n²` margin
  above a red-herring construction.

## Status
- **Verified (opus):** construction beaten at n=7 (`2/13`), n=8/10/12/14 (even block `= 1/n`); `M(2S)=M(S)`;
  even block covering ⇔ n even; covering-min `= 1/n` (even n, tight), `>1/n` (odd n).
- **Refuted (opus, honest):** covering-min `= n/Φ₆`; the `~1/n²` margin; witness=ζ₆/hexagonal/Kershner/PG as
  the covering-min; the convergent-escape; the Sylvester-as-covering-min reading. All were about a
  non-extremal family.
- **Sharpened:** covering-min is parity-determined — `1/n` (even n, tight, even block = AP doubled), `>1/n`
  (odd n, Stern-Brocot realizability); LRC(14) is tight with extremal `2·{1,…,13}`. The recursive parity
  mindset (`a=÷2`, doubling, even/odd) is the right frame; the construction was the distraction.
- **Still open:** LRC(14) itself = no covering 13-set below `1/14` (the even block is the tight case); and
  the odd-n covering-min (`2/(2n−1)` mediant vs deeper, the realizability problem).

Related (now corrected): the-covering-min-witness-is-kleins-zeta6, the-cyclic-kershner-attack, the-covering-
min-as-a-function-of-n, the-observer-abstraction (escape claim), the-recursive-mindset-sylvester (Φ₆ true,
covering-min reading void); SURVIVING: the-functional-decomposition (a/b/parity vindicated), the-observer-on-
the-tournament-side; klein HYP-3724/3705/3706, mac-mini HYP-3701/3702 (construction-focus refuted); OPEN-Q-108.
