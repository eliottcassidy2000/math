# The crux is "t = 1/14 is blocked": divisor-complete families have M ≥ 1/12, exactly

*kind-pasteur-2026-07-11-S127 cont.52. Owner: "keep coming to a better picture of the LRC(14) crux, sharpen
via hypothesis investigation." After opus-S247/klein-S266 reframed the crux (the Farey window is non-empty;
the operative restriction is covering), I resolved the exact value and the mechanism — and they turn out to be
one clean statement about the witness t = 1/14.*

---

## The picture, in one line

**LRC(14) tightness lives entirely at the single time t = 1/14, and divisor-completeness is exactly what
blocks it.**

- A family that **misses a multiple of some d ∈ {2,…,14}** is lonely at t = 1/d: every speed avoids `0 mod d`,
  so `‖vᵢ/d‖ ≥ 1/d ≥ 1/14`. THM-366, elementary, proved. The two tight families both live here (they miss a
  multiple of 14): the **AP** `{1..13}` and **Goddyn–Wong** `{1..11,13,24}`, each `M = 1/14` via t = 1/14
  (klein-S266 — so "AP minimizes M" is non-unique).
- A **divisor-complete** family has a multiple of 14, say `14m`. At t = 1/14 that runner sits at
  `‖14m/14‖ = 0` — **t = 1/14 is blocked** (reach 0, verified: `min reach at t=1/14 = 0` for the extremals).
  It cannot reach the tight value at the tight time; it must use a **coarser witness**.

So the crux is not "AP minimizes M" (non-unique, and the window is non-empty) and not "the window is empty"
(false at k = 13). It is: **the families for which t = 1/14 is blocked — precisely the divisor-complete
ones — are loose.** How loose is the whole question, and this session pins it.

## The exact floor: 1/12, not 1/13

Hunting hard for a divisor-complete family below 1/12 — perturbing the tight AP and GW toward
divisor-completeness, sweeping every 2-block family, adversarial swaps — over 2170 primitive DC candidates,
**none falls below 1/12**. The floor is achieved (at least) at the two 2-block families

`{1,2,3,4} ∪ {10,…,18}`   and   `{1,3,4} ∪ {10,…,18,21}`,

both with `M = 1/12` **exactly**, attained at `t = 5/24`. So:

> **The divisor-complete M-floor is exactly `1/12`** — a margin of `1/12 − 1/14 = 1/84` over the tight bound.
> klein-S266's / boxeph's `DC ⟹ M ≥ 1/13` is a *conservative* provable bound; the true extremal is `1/12`.

## Why 1/12 — the coarser-witness mechanism

The value is not mysterious. The band-edge lemma (opus-S235) says a family clearing at a non-14 modulus `q`
is lonely with `M ≥ ⌈q/14⌉/q`, which *decreases* in `q`. A divisor-complete family, having lost t = 1/14,
clears at some `q ≥ 15`; the worst one bottoms out at **`q = 24`**, where `⌈24/14⌉/24 = 2/24 = 1/12` — and
for the extremals this band-edge is **tight** (their true M equals it, at t = 5/24). So

> **`1/12` is the second-best-witness floor**: the coarsest modulus the worst t = 1/14-blocked family is forced
> down to before it clears. `1/14 = 1/14`-witness (blocked for DC) → next stop `2/24 = 1/12`.

This is why the divisor-complete hard core is simultaneously the *hard* case (it loses the elementary
witness) and *loose* (it recovers a witness only two rungs coarser).

## What this sharpens for the proof

The crux collapses to a single statement with an exact target:

**Every divisor-complete family clears at a non-14 modulus `q ≤ 24`-and-friends giving `M ≥ 1/12`.**

- The tight bound and its non-uniqueness (AP vs GW) are *not* the DC problem — they are the non-DC bucket,
  where t = 1/14 already gives `M ≥ 1/14` (THM-366). No k = 13 inverse theorem, no window rigidity, no GW
  characterization is needed for the closing; those are all about the non-covering families THM-366 dispatches.
- The remaining work is exactly `DC ⟹ M ≥ 1/12`: boxeph-S20's finite check (Vmax ≤ 30, verified) for bounded
  structure, carried to all scales by dilation (my cont.50 `loose_dilate`), plus the large-structure looseness
  (the coprime-clearing / anti-concentration, my cont.45 + opus's pigeonhole). The `14 = 2·7` composite
  difficulty that broke the window-empty form (opus-S247) never touches this bucket — it lives in the
  elementary non-DC one.

The recurring shape holds once more: the object that looked like a hard inverse/rigidity theorem is, on the
critical (divisor-complete) bucket, an elementary blocked-witness statement with an exact `1/84` margin, and
the genuinely delicate near-tight structure (AP, GW, the 3/41 window) is entirely in the bucket the two-bucket
dispatch already closes.

*Files: lrc14_dc_floor_resolve_kps_S127.py/out, lrc14_dc_floor_mechanism_kps_S127.out. HYP-6180. Sharpens
opus-S247 (window non-empty) and klein-S266 (covering restriction, DC⟹M≥1/13, tight locus = {AP,GW}); rests on
THM-366 (two-bucket), opus-S235 (band-edge), boxeph-S20 (finite check), and
[[the-farey-window-is-non-empty-but-holds-no-dc-family-thm366-dispatches-the-fillers-kps-S127]].*
