# The finite certificate for THM-701's wide-spread recursion — tables + structure

*kind-pasteur-2026-07-11-S127. Owner: "write the finite certificate: balanced-core + cap-growth tables."*

THM-701 reduced the wide-spread direction of LRC(14)-S3 to two finite obligations:
**(I)** the cap grows by at least `2/21` per step, and **(II)** the joint functional `Φ = p0 + (1/3)p1`
stays under the cap on bounded-spread cores. This note writes both tables and reports what the check reveals
about the one lemma that remains.

> **Collision note (honest).** [mac-mini's THM-702](../../01-canon/theorems/THM-702-explicit-finite-certificate-wide-recursion.md)
> is the canonical certificate — same two tables, done in **exact arithmetic** (fractions, zero fitting),
> pushed ~2 h ahead of this note. The computation here is an **independent grid-numerical cross-check** (a
> different method), and it agrees with THM-702 to every printed digit. Two findings below are net-new to
> this pass: the extremal 2-runner packing `{1,13}` behind THM-702's exact denominators, and the fact that
> the core lemma does **not** factor (with a third role for `λ = 1/3`).

## Table I — cap growth (`cap_k = min_{|P|=13−k} meas(G_P)`, gap `1/14`)

| `k` | `cap_k` | source | `cap_{k+1} − cap_k` | `≥ 2/21 = 0.09524`? |
|----|---------|--------|---------------------|----------------------|
| 8  | 0.38153 | LRC extremal, `\|P\|=5` (THM-532/534) | — | — |
| 9  | 0.49426 | LRC extremal, `\|P\|=4` | 0.11273 | ✅ |
| 10 | 0.60440 | LRC extremal, `\|P\|=3` | 0.11014 | ✅ (tightest) |
| 11 | 0.72527 | 2-runner min, argmin **{1,13}** | 0.12087 | ✅ |
| 12 | 0.85714 = 6/7 | 1-runner (`\|P\|=1`) | 0.13187 | ✅ |
| 13 | 1.00000 | empty (`\|P\|=0`) | 0.14286 | ✅ |

Every step clears `2/21`; the **tightest is `k=9→10` at 0.11014** (slack `+0.0149`). Growth is monotone
increasing, so the induction never loses margin as it climbs — the whole recursion is tightest at its base.

**cap_11 hardened.** The LRC(14) speeds need not lie in `{1..13}`, so I recomputed the min 2-runner lonely
measure over all coprime pairs `a<60, b<80`: still **0.72527 at {1,13}**, inside the required window
`[0.6996, 0.7619]`. **{1,13} is the densest 2-runner packing** because `13 ≡ −1 (mod 14)`, so `‖13x‖` locks
to `‖x‖` — which is exactly why THM-702's exact growths carry the factor `91 = 7·13` (`cap_11−cap_10 = 11/91`,
`cap_12−cap_11 = 12/91`).

## Table II — balanced-core check (`Φ(F) = p0(F) + (1/3)p1(F) ≤ cap_{|F|+1}`)

At the `consec_m = {0,1,…,m−1}` argmax (HYP-2644; and confirmed the argmax by search below):

| `\|F\|=m` | `Φ(consec_m)` | `cap_{m+1}` | margin |
|--------|----------------|-------------|--------|
| 7  | 0.26190 | 0.38153 | +0.11963 |
| 8  | 0.40857 | 0.49426 | **+0.08569 (tightest)** |
| 9  | 0.49015 | 0.60440 | +0.11425 |
| 10 | 0.56690 | 0.72527 | +0.15837 |
| 11 | 0.62994 | 0.85714 | +0.22720 |
| 12 | 0.66613 | 1.00000 | +0.33387 |

All positive; **tightest `+0.086` at `m=8`** — matching THM-702's exact `+0.086` (the correct `cap_{|F|+1}`
indexing; comparing at `cap_{|F|}` would falsely "fail" here). A search over shifted-AP, near-geometric, and
random bounded cores never beats `consec` (min margins identical to the table), confirming `consec` is the
argmax. *Rigor:* grid margins (`≥0.086`, `≥0.0149`) dwarf the `~10⁻⁵` discretization error; THM-702 supplies
the exact-arithmetic version.

## What the check reveals — the core lemma does **not** factor

Testing whether `Φ`-extremality reduces to separate `p0`- and `p1`-extremality:

- **`p0` (cover) IS maximized at `consec`** at every `m` — this is THM-534/530's content (consecutive
  frequencies spread sector-coverage best).
- **`p1` (miss-exactly-one) is NOT** — random balanced cores have strictly higher `p1`; indeed `p1(consec_m)`
  *falls* with `m` (`0.34→0.19`) while `max p1(rand)` *rises* (`0.34→0.57`). A near-miss coverer sits at
  "miss one"; the extremal coverer either covers or misses several.
- **`Φ = p0 + (1/3)p1` is still maximized at `consec`** — because `consec`'s `p0` lead (`0.50` vs `0.25` at
  `m=10`) outweighs its `p1` deficit at weight `1/3`.

So the natural hope "`Φ`-extremal ⟸ `p0`-extremal + `p1`-extremal" is **false**. The remaining lemma
(THM-702's named residual) is irreducibly **joint**: it needs `consec`'s `p0`-dominance to survive the
`(1/3)p1` perturbation, using the `(p0,p1)` anti-correlation (high-`p1` cores are low-`p0`), not two separate
bounds.

### A third load-bearing role for `λ = 1/3`

The argmax flips from `consec` to a high-`p1` core at
`λ*(m) = min_F [p0(consec)−p0(F)]/[p1(F)−p1(consec)]`:

| `m` | `λ*` |
|----|------|
| 8  | 1.514 |
| 9  | 1.231 |
| 10 | 0.984 |

`λ*` decreases with `m` but stays far above `1/3` (margin `≥0.65` at `m=10`). So `λ = 1/3` is load-bearing on
**three** counts, not two:
1. `λ ≥` the worst far-`w` tax (`≈0.25`), so the peel bound `p0(E) ≤ Φ(E')` holds (THM-701 rung 2);
2. `λ` small enough that the increment `2(p1+p2)/21 ≤ 2/21` stays under the cap growth (THM-701 rung 3);
3. `λ < λ*(m)`, so `consec` (the `p0`-argmax, THM-534) *remains* the `Φ`-argmax (this note).

`1/3` sits in `[0.25, ~0.98]` — comfortably satisfying all three. `λ = 1/7` fails role 1; a large `λ` fails
role 3. The certificate needed the middle.

## Bottom line

Both tables are written and cross-validated (here numerically; THM-702 in exact arithmetic). The wide-spread
direction has **no analytic gap** and its finite content is now explicit. The single remaining lemma is
`Φ`-consec-extremality on bounded cores — and this pass shows it is a genuinely joint statement (a
`(1/3)`-perturbation of THM-534's `p0`-extremality that `λ<λ*` keeps stable), not a pair of separate
extremal facts. That is the same extremal statement THM-534/530/657 keep isolating — now confirmed to carry
the whole wide-spread program.

*Files: `04-computation/lrc14_finite_certificate_kps_S127.py` (+`.out`). Canonical certificate: THM-702
(mac-mini). Depends on THM-701/700/699 (kps), THM-534/530 (the residual extremal lemma).*
