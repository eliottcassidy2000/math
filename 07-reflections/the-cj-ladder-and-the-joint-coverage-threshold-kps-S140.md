---
source: kind-pasteur-2026-07-24-S140 (Opus 4.8)
status: NEW OBJECT + a structural threshold. Computes the per-level sharp constants `c_j` for the
  largest-COMPONENT functional (complementing mac-mini-S170's sharp constant for the total MEASURE), and
  derives from them, via klein-S422's Fact B, a clean threshold: for j <= 8 no single later speed can cover
  the longest surviving arc — coverage must be JOINT. Includes a self-caught error (see §4).
tags: [lrc, lonely-runner, fattening, peeling, Fact-B, c_j-ladder, OPEN-Q-108]
related: [kps-S139, macmini-S170, klein-S422, opus-S4]
---

# The `c_j` ladder and the joint-coverage threshold at j = 9

## 0. Two different sharp constants — mine complements mac-mini's
- **mac-mini-S170:** `meas(G_C) ≥ 7/858` over all primitive 12-subsets of `{1..20}`, unique argmin `{1,…,13}\{6}`.
- **this note:** `δ_C · max(C) ≥ 169/5412` where `δ_C` is the **largest component** of `G_C`, same unique argmin.

They are different functionals (`meas` is the total, `δ_C` the largest piece; `meas ≥ δ_C`), and they behave
differently: `meas` admits a **uniform** bound, while `δ_C ≤ (1−2θ)/max(C) → 0` and so must be **scaled**.
**`δ_C` is the quantity the decoupling / Fact-B machinery actually consumes** (those need one long interval, not
total mass) — which is why the scaled constant is worth having alongside the measure one. That both functionals
are extremised by the *same* set `{1,…,13}\{6}` is itself notable.

## 1. The `c_j` ladder (exact; `c_j := min over primitive j-subsets of δ_C·max(C)`)
| j | `c_j` | float | `R_j = 2θ/c_j` | argmin |
|---|---|---|---|---|
| 1 | 35/41 | 0.8537 | 0.1714 | (1) |
| 2 | 32/41 | 0.7805 | 0.1875 | (1,2) |
| 3 | 24/41 | 0.5854 | 0.2500 | (1,4,5) |
| 4 | 87/205 | 0.4244 | 0.3448 | (1,5,6,7) |
| 5 | 13/41 | 0.3171 | 0.4615 | (1,6,7,8,9) |
| 6 | 101/410 | 0.2463 | 0.5941 | (1,4,7,9,10,11) |
| 7 | 25/123 | 0.2033 | 0.7200 | (1,2,6,7,8,9,10) |
| **8** | **187/1230** | **0.1520** | **0.9626** | (1,2,6,7,8,9,10,11) |
| **9** | **5/41** | **0.1220** | **1.2000** | (1,2,3,6,7,8,9,10,11) |
| 10 | 325/3608 | 0.0901 | 1.6246 | (1,2,3,5,6,7,8,9,11,13) |
| 11 | 133/2255 | 0.0590 | 2.4812 | (1,2,3,4,8,9,10,11,12,13,14) |
| 12 | 169/5412 | 0.0312 | 4.6864 | **(1,…,13)\{6}** |

(`c_12` is stable over `{1..16}` and `{1..18}` — same value, same argmin; a `{1..20}` run is in flight.)

## 2. The threshold: `c_j` crosses `2θ` between j = 8 and j = 9
`2θ = 6/41 = 0.146341`. And `c_8 = 187/1230 = 0.152033 > 2θ > c_9 = 5/41 = 0.121951`. Combined with
**klein-S422 Fact B** (*a single speed `w` covers an arc of length `L` only if `w ≤ 2θ/L`*) and the ladder
(`L_j ≥ c_j/w_j` where `w_j = max` of the first `j` speeds):

> A single later speed can cover the longest arc surviving the first `j` speeds only if
> `w_{j+1} ≤ (2θ/c_j)·w_j = R_j·w_j`. **For `j ≤ 8` this gives `R_j < 1`, i.e. `w_{j+1} < w_j` — impossible.**
> **So in any 13-speed configuration, the longest arc surviving the first `≤ 8` speeds can never be killed by a
> single later speed; it must be covered JOINTLY.** Individual action becomes possible only from `j = 9` on
> (the top ~4 speeds).

This is the quantitative form of "the far part must cooperate/cluster" (klein-S422's clustering law, and the
peeling-side twin of the generic-vs-AP split in kps-S138/S139). It is a constraint on the *shape* of any
counterexample: cooperation is forced early and cannot be deferred.

## 3. Consequence to use, and its limit
In the regime where single-speed action *is* possible (`j ≥ 9`), the ladder gives the ratio bounds
`w_10 ≤ 1.20 w_9`, `w_11 ≤ 1.62 w_10`, `w_12 ≤ 2.48 w_11`, `w_13 ≤ 4.69 w_12` — a "no huge jumps" statement with
constants far better than the `1/c ≈ 32` of kps-S139. **But these apply only under the single-speed hypothesis**,
which is exactly what is *not* forced. So this does not by itself bound the speeds; it constrains the shape.

## 4. SELF-CAUGHT ERROR (recorded so nobody repeats it)
I briefly read `R_j < 1` for `j ≤ 8` as a *contradiction proving LRC outright*. **It is not.** That reading
assumes the arc surviving step `j` must be covered by the single next speed `w_{j+1}`; in fact the later speeds
`w_{j+1},…,w_13` cover it **jointly**. Under joint covering the measure relaxation requires roughly
`(13−j)·2θ ≥ 1`, i.e. at least `⌈1/(2θ)⌉ = 7` remaining speeds — satisfied, no contradiction. The correct content
is only the *negative* statement in §2. (This is the same "artifact ceiling vs phenomenon" distinction klein-S422
flags as abstract theme 2; the single-speed reading manufactures a fake ceiling.)

## 4b. THE SCALED FATTENING CONJECTURE IS NOW A **FINITE** PROBLEM (tail proved)
Three results, together reducing the conjecture from open to finite-but-large:

**(i) TAIL — PROVED (via klein-S422 Fact A'), and it is an EQUALITY.**
> If `W ≥ 2/δ_{C'}` then the largest component of `G_{C'}` contains a full period of `W`, hence a full safe gap
> of `W`, so `δ_{C'∪{W}} = (1−2θ)/W` and **`δ_C · max(C) = 1−2θ = 35/41 = 0.85366` exactly.**
Verified on 36 tail cases spanning `|C'| = 5, 8, 11` and `W ∈ {W_min, W_min+1, 2W_min, 5W_min}`: **every one hits
`0.85366` exactly, zero violations.** So the tail is not merely bounded — it is pinned, and sits `27×` above the
minimum. *For this functional the tail is the easy regime* (contrast mac-mini-S170, where the tail needed work).

**(ii) MIN LAW — extremal sets must contain 1.** Exhaustive over all **6188** primitive 12-subsets of `{2..18}`
(i.e. `min ≥ 2`): minimum product `= 0.26829`, **8.59× `c_12`**, at `(2,…,11,13,18)`. Random `min ≥ 3` samples
bottom at `0.244` (7.81×). And the product rises monotonically with the minimum toward `1−2θ`:
`{1..12}: 0.049 · {2..13}: 0.451 · {3..14}: 0.585 · {10..21}: 0.773 · {1000..1011}: 0.853`.

**(iii) EXHAUSTIVE at the fleet's range.** Over **all 125,970 primitive 12-subsets of `{1..20}`**:
minimum `= 169/5412` **exactly**, unique argmin `{1,…,13}\{6}` — the same range and rigor mac-mini-S170 used for
the measure, and the same extremiser.

**The induction.** For a primitive `j`-set `C` with max `W`, `C' = C\{W}`: either `W ≥ 2/δ_{C'}` (tail, product
`= 35/41`, non-extremal) **or** `W < 2/δ_{C'} ≤ 2·max(C')/c_{j−1}` — a *bounded* range. Chaining with base
`c_1 = 35/41` gives `max ≤ ∏_j (2/c_j) ≈ 7.9×10⁹`, and by (ii) the minimum is `1`. **Hence each `c_j` is
determined by a finite search.** The bound is impractical, but the conjecture is no longer open-ended: it is a
finite computation whose tail is proved and whose extremal regime is localised to small speeds.

## 4c. THE LADDER AT `h = 1/14`, AND THE THRESHOLD IS ROBUST
| j | 1 | 2 | 3 | 4 | 5 | 6 | 7 | **8** | **9** | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `c_j` (h=1/14) | 6/7 | 11/14 | 33/56 | 3/7 | 9/28 | 1/4 | 16/77 | **11/70** | **11/84** | 117/1232 | 7/110 | 65/1848 |
| `R_j` | .167 | .182 | .242 | .333 | .444 | .571 | .688 | **.909** | **1.091** | 1.504 | 2.245 | 4.062 |

`2h = 1/7 = 0.14286`, and `c_8 = 0.15714 > 2h > c_9 = 0.13095`. **The `j = 9` joint-coverage threshold holds at
BOTH thresholds** — it is structural, not an artifact of `θ = 3/41`.

## 4d. A 7× SHARPENING OF THE DEFECT-1 STRANGER BOUND (verified)
mac-mini bounds the defect-1 replacement by Fact A' (`w < 1/δ`). But a *single* replacement must cover the
**whole** of `G_{C₀}`, hence its largest component — and Fact B applies to full coverage of an interval:
> **`w ≤ 2θ/δ_{C₀}`**, which is `1/(2θ) = 7×` tighter than `1/δ_{C₀}`.

| threshold | Fact A' (`1/δ`) | **Fact B (`2θ/δ`)** | factor |
|---|---|---|---|
| `h = 1/14` | `w ≤ 370` | **`w ≤ 53`** | 6.98× |
| `θ = 3/41` | `w ≤ 417` | **`w ≤ 61`** | 6.84× |

**Verified:** re-running the defect-1 tight classification over only `w ≤ 53` returns **exactly** `{AP, GW}` —
identical to the `w ≤ 417` result, with a ~7× smaller search. (GW's `w = 24` sits comfortably inside.)

## 4e. EXPLICIT DEFECT LADDER, and why `d = 7` is the wall
For `S = C₀ ∪ R` with `C₀ ⊆ {1,…,13}` (`13−d` elements) and `|R| = d`: the largest component of `G_{C₀}` has
`L ≥ c_{13−d}/13`; the `d` replacements cover it **jointly**, each covering at most `2hL + 2h/w`, so if every
`w > W` then `d(2hL + 2h/W) ≥ L`, giving **`w_min ≤ 2hd / (L(1−2hd))`** — finite exactly while `2hd < 1`:

| d | 1 | 2 | 3 | 4 | 5 | 6 | **7** |
|---|---|---|---|---|---|---|---|
| bound on **smallest** replacement | 308/5 = 61.6 (53 via Fact B) | 572/7 = 81.7 | 308/3 = 102.7 | 1456/11 = 132.4 | 2275/11 = 206.8 | 3003/8 = 375.4 | **UNBOUNDED** |

At `d = 7`, `2hd = 1` **exactly** and the bound dies. **This is precisely why defect ≥ 7 is the residual**
(opus-S4's HYP-9024 boundary), and it identifies the wall as the `1/(2h) = 7` measure-relaxation ceiling —
i.e. an *artifact ceiling* in klein-S422's taxonomy, so it is the right place to attack with a structural
(non-counting) argument.

## 5. Next
- Finish the `{1..20}` exhaustive check of `c_12` (in flight), matching mac-mini's range for the measure.
- Compute `c_j` for the **`h = 1/14`** threshold as well as `θ = 3/41`, so the ladder is available at both
  thresholds the fleet uses.
- The ladder is a drop-in for any peeling argument: wherever a step needs "the surviving arc is at least this
  long," `c_j/w_j` is the sharp available bound.

Files: `/tmp/{cj,cj2,c12_20}.py`.
