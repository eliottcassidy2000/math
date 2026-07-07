# The Fan–Sun gcd template IS the fleet's covering system — "divisibility = clearing" unifies the two frames

*klein-2026-07-06-S147 (HYP-4621). Owner: work the residual; find past work / analogous connections
that unlock the route. opus-S116 already identified the external home of our crux (Fan–Sun,
arXiv:2306.10417); I pulled the actual template and found that **its gcd/divisibility case-split is
literally the fleet's small-modulus covering layer** — so the two parallel proof frames the fleet has
been running (Fan–Sun spectrum-order vs kps's covering vs mac-mini's bands) are ONE proof, with a
clean two-layer structure and a bounded modulus.*

## The residual, and the two frames it was being attacked in

The open node (klein-S144): **compressed non-AP 12-family ⟹ cleared at some `q ≤ 39`** — the covering
completeness of `(C)`. Two frames were running in parallel:
- **Fan–Sun / spectrum-order** (opus-S116): `(C)` = the `n=12` first-gap Lonely Runner Spectrum case;
  the amended spectrum is `ML = s/(ns+k)`, `k ≤ n`; the window `(1/13,2/25)` needs order `k ≥ 2`.
  Fan–Sun's `n=4` gap-emptiness proof is a **gcd/divisibility case analysis**.
- **kps covering** (S43–46): every non-AP compressed family clears at some `q ≤ 39` via an explicit
  `rational_point_margin` cert; the small-`q` layer `LRCSmallModFloor` (no multiple of `q ⟹ M ≥ 1/q`).

## The unlock: they are the same proof

**Order = the covering denominator.** For `n = 12`, write `M = r/Q` in lowest terms; then the
Fan–Sun form `M = s/(12s+k)` has `s = r` and **order `k = Q − 12r`** (`lrc14_fansun_gcd_is_covering`):
`1/13` is `k=1,s=1`; `2/25` is `k=1,s=2` — the two window *endpoints* are Kravitz rungs `k=1`; a value
*strictly* inside needs `k ≥ 2` (opus-S116 `k<s<2k`). So the window is exactly "order `≥ 2`."

**The gcd case-split IS the small-`q` covering layer.** Fan–Sun split on divisibility; kps's
`LRCSmallModFloor` says a family missing modulus `q` (no speed divisible by `q`) has `M ≥ 1/q`. These
are the same step, and the near-AP families make it vivid: the AP `{1,…,12}` covers every `q ≤ 12`
(each `q` divides `q`), with `q ∈ {7,8,9,10,11,12}` each having a **unique** multiple. A "13-lift"
`v_i ↦ i + 13k_i` (the near-AP moat, `≡ AP mod 13`) that moves the unique carrier of some
`q ∈ {7..12}` **breaks that divisibility ⟹ misses `q` ⟹ `M ≥ 1/q > 2/25`, cleared instantly.**

> **Verified (39,987 non-AP 13-lifts): every family that breaks a small-modulus divisibility clears
> *at that broken modulus* — 33,444/33,444 = 100%.** The Fan–Sun gcd branch = the covering's `q ≤ 12`
> layer, on the nose.

## The clean two-layer structure of (C)

The 100% fact partitions the crux into exactly two layers:

1. **gcd / small-`q` layer (`q ≤ 12`)** — a compressed non-AP family that **misses** some `q ∈ {2..12}`
   is cleared by `LRCSmallModFloor` (`M ≥ 1/q`). This is Fan–Sun's divisibility case, already GREEN
   in Lean. It disposes of every family whose lift breaks a divisibility.

2. **the near-AP moat (`13 ≤ q ≤ 32`)** — the survivors are exactly the **divisibility-preserving**
   families (cover *all* of `{2,…,12}`, `≡ AP mod 13` at the unique-carrier moduli). These the
   small-`q` layer cannot see (they look like the AP mod every small `q`). **Verified: 24,772
   divisibility-preserving non-AP lifts all clear at `q ∈ [13, 32]`** (histogram peaked at `17`;
   max `32`). This is kps-S46's "13-lift residual cleared by a fixed covering" — refined: the range
   is `[13,32]`, not the circulated `{11..23}` (same over-optimism I flagged in S144).

3. **AP exception** — `{1,…,12}` covers all `q ≤ 12` *and* is `≡` itself mod 13, so it is cleared by
   nothing: `M = 1/13`, the sole survivor (kps-S1 HYP-4096 / klein-S140).

So `(C)` = **[break a divisibility → cleared at `q ≤ 12`] ⊕ [preserve all divisibilities → near-AP
moat, cleared at `q ∈ [13,32]`] ⊕ [AP]**. The order-`k` formula finitizes the moat (opus-S116
O-korder): moat values have bounded order, hence lie among finitely many forms `s/(12s+k)`.

## Why this unlocks the route

- It **collapses three parallel attacks into one** (Fan–Sun gcd = kps small-`q` layer = mac-mini's
  "pins" in THM-619; the residual moat = mac-mini's "bands" = opus's decorrelation = kps's 13-lifts).
  The fleet can stop maintaining separate frames and formalize the *single* two-layer covering.
- It gives the **exact clearing-modulus budget**: `q ≤ 12` (gcd layer, GREEN) and `q ∈ [13,32]` (moat).
  The covering set is `{2,…,32}` (not `{6..39}` or `{11..23}`) — sharper and Lean-sized.
- It imports Fan–Sun's **proof template** for the moat: their `n=4` argument bounds the order `k` and
  handles the finitely-many generalized-AP exceptions with a sub-AP cap — exactly opus-S115/S116's
  subfamily-cap. Carrying their order-bound to `n=12` is the surviving mathematical obligation, and it
  is a *finite* classification, not an analytic rigidity.

## Honest residual

The moat (divisibility-preserving lifts clear at `q ∈ [13,32]`) is verified, not proved uniformly.
The proof = (a) Fan–Sun order bound `k ≤ K₀(12)` (finitely many moat forms) + (b) each form's
generalized-AP exception capped by a sub-AP rung `≥ 2/25` (opus subfamily cap) + (c) the peel for
non-compressed (klein-S144). The gcd/`q ≤ 12` layer is done (`LRCSmallModFloor`). No new Lean this
session — this is the synthesis that says *which* known template closes the moat and how it wires to
the fleet's covering.

## Links

- Scripts: `04-computation/lrc14_fansun_gcd_is_covering_klein_S147.py` (+ `.out`),
  `lrc14_residual_orderk_klein_S147.py`/`.out`. HYP-4621.
- Builds on / credits: **opus-S116** (identified (G) = Fan–Sun first-gap; O-korder/O-gcd obligations),
  **Fan–Sun arXiv:2306.10417** (amended spectrum + gcd template), kps-S43/44/46 (covering,
  `LRCSmallModFloor`, 13-lift residual), opus-S124/126/127 (mod-25/13 factoring, decorrelation),
  mac-mini-S33/S48/49 (THM-633/634, THM-619/620 pins+bands), klein-S140 (rigidity)/S144 (compressed
  uniformity). Open: the moat order-bound + sub-AP cap (Fan–Sun template at `n=12`) + wiring.
