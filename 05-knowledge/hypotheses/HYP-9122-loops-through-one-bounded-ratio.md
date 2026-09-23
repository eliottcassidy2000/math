---
id: HYP-9122
title: "Loops through 1 of every length with ratio below 3 (robust 1-escape)"
status: >
  OPEN HYPOTHESIS with FINITE-EXACT support to s=40. For every s>=1 there is
  a forward E-cycle through 1 with s multiplications and K halvings,
  2^K<3^(s+1). Together with the PROVED 1-escape lemma this handles every m
  with v_3(m-1)>=2 in Q2. For s>=2 every loop has ratio >= 13/9, and the
  observed minimal ratios lie in [1.517, 2.96], thin at s=11 and s=23. The
  robust form allows other exit words and is what the E-SCC induction
  actually needs.
source: collatz-procgen-20260922 (mac-mini)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260922_choice_ladder.md
---

# HYP-9122 -- bounded-ratio loops through 1

Refutation form: an s for which every loop has 2^K >= 3^(s+1), together
with no alternative exit of bounded ratio.

## Structure of the minimal loops (FINITE-EXACT, s<=22)

The minimal loops are built from a few primitive cycles:

* the trivial loop at 1 (`k=2`, ratio `4/3`);
* the loop at 5, `5->13->4->5` (`k=3,0,2`, ratio `32/27`);
* the loop at 85, `85->28->37->49->16->85` (`k=0,2,2,0,4`, ratio
  `256/243`), which is the convergent clock `8/5` of `log_2 3`.

For example, the `s=9` and `s=14` loops are
`1 -(8)-> 85 [-> 85-loop]^t -> 28 -> 37 -> 49 -> 16 -> 5 -> 13 -> 4 -> 1`.
The resulting family has length `9+5t+3u+v` and ratio
`1.665 (256/243)^t (32/27)^u (4/3)^v`. It grows like `1.0535^t` and passes
`3` near length 100.

A bounded family therefore needs cycles whose excess vanishes. Candidates
are positive E-cycles at the upper convergent clocks of `log_2 3`
(`8/5, 65/41, 485/306, ...`), entered from 1 and exited back to 1 at
bounded total cost. What remains is a carry-covering statement: at each
such clock some legal E-word has carry divisible by `2^K-3^L`, with an
integer fixed point. The heuristic count is about `2^(0.9L)` hits. It is
not proved.

## Exits A/B

At every depth `3<=k<=41` of the 1-neighbourhood, at least one exit
works: (A) a loop of length `k-1` with ratio `<3`, or (B) a loop of length
`k-2` with ratio `<9/4` followed by the mod-27 descent `(4x-4)/9`.

## Sharpened by the loops lane (collatz_procgen_20260922_loops_and_escapes.md)

* **PROVED.** Every loop through 1 with `s>=2` has `2^(K+1)>3^(s+1)`. So
  HYP-9122 is **equivalent** to "`minK(s)=K0(s)=floor((s+1) log_2 3)` for all
  `s>=2`", and the 1-escape at depth `k` is the only descent within `k`
  moves.
* **PROVED (record reduction, Theorem 4.2).** HYP-9122 follows from one
  base loop `B_q` and one hub cycle `C_q` per record (upper
  best-approximation) denominator `q` of `log_2 3`: `3, 5, 17, 29, 41, 94,
  147, 200, 253, 306, 971, 1636, 2301, ...`. These are spliced in Ostrowski
  order.
* **FINITE-EXACT.** `minK(s)=K0(s)` for all `s<=6000`. Explicit verified
  loops for all `s<=2500` come from 13 base loops and 11 cycles through the
  hubs `1, 4, 13, 40, 121, 1093`.
* **PROVED obstruction.** Loop height is at least about `s/(3 eps(s+1) ln 2)`,
  so no finite family works, and record base loops cannot come from
  splicing. The hard base loops climb, wander at height about `s/eps`, and
  land on a trunk value `R_j=(4^j-1)/3`, for example `341` or `5461`.
* **OPEN.** Records `q>5626`.

