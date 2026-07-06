# The CRT-fold of the Farey-cell voids: composite q reduce to small-modulus multi-color problems

**kind-pasteur-2026-07-05-S16 (HYP-4197).** A creative tool for mac-mini's THM-622
(the second-value gap `(1/13,2/25)` = a Farey-cell void, "(G) ⟺ no 12-set attains
`c/q` with `c ≥ 3`, `q ∈ (12.5c, 13c)`"). It extends my S11 parity-split from `q = 2p`
to general composite `q`, folds the small-`c` cases to small moduli, and characterizes
the max-attained side via pinning at the prime factors.

## The CRT-fold (extends S11)

A **clearance-`c` witness at `q`** is a unit dilation `A` with all runners
`dist(v·A, q) ≥ c`. The danger band `{|r|_q < c}` is a set of `2c−1` residues. For
`q = a·b` (coprime), by CRT each residue is `(r mod a, r mod b)`, so the danger band
factors into a set of **joint color-classes** `(r mod a, r mod b)`, and the witness is a
choice of dilation avoiding those classes — a **multi-color avoidance** on `Z/a × Z/b`.

- **`q = 2p` (parity-split, S11):** even runners avoid `E_p`, odd runners avoid `O_p`.
  For `c = 3`, `E_p = {0, ±2}`, `O_p = {±1}` mod `p`. **Verified EXACT** for `c=3, q=38`
  (`p=19`): brute-force `≡` parity-split, 0/50,000 mismatches.
- **`q = 3·17 = 51` (`c=4`):** the danger band `{0,±1,±2,±3}` mod 51 factors into 7 classes
  `(0,0),(0,3),(0,14),(1,1),(1,15),(2,2),(2,16)` on `Z/3 × Z/17` — a 3-coloring by `mod 3`
  with `mod-17` avoidance. **Verified.**

## The small-`c` cells all fold; the prime cores start at q = 89

The Farey-cell interior points `c/q` (`c ≥ 3`, `q ∈ (12.5c, 13c)`), by increasing `c`:

| c | q (all interior) | factorization | reducible? |
|---|---|---|---|
| 3 | 38 | 2·19 | ✔ parity-split |
| 4 | 51 | 3·17 | ✔ 3-color |
| 5 | 63, 64 | 9·7, 2⁶ | ✔ |
| 6 | 77 | 7·11 | ✔ |
| 7 | 88, **89**, 90 | 8·11, **prime**, 2·45 | 89 IRREDUCIBLE |

So **every Farey-cell void with `c ≤ 6` (the shallowest, first to attack) folds to a
small modulus**; the first prime core is `q = 89` (`c = 7`). The mediant `3/38` (mac-mini's
first target) is the `c=3` case — it folds to a **mod-19** problem.

## The max-attained side, folded: `M = 3/38 ⟹ pinned at ±1 mod 19`

The fold simplifies the *witness* side, but (mac-mini) the contradiction is on the
*max-attained* side. That side folds too: `M(B) = 3/38 = 1.5/19` forces the covering-min
*restricted to denominator 19* to be `≤ 1/19`, i.e. **`B` is pinned at `±1 mod 19`**
(every unit dilation `a` has some runner `≡ 0, ±1 mod 19`). Similarly pinned mod 2. So a
`3/38`-attainer is:

> pinned at `±1 mod 19` **and** carries a clearance-3 (parity-split) witness mod 38 **and**
> covers `2..12` — a joint condition on `Z/2 × Z/19` plus the integer covering.

By the S11 **parity dividend**, mod-19 pinning does *not* by itself block the mod-38
witness (satisfiable) — matching mac-mini's "residue side satisfiable." The contradiction
must use the *other* denominator-38 competitors (`t = 3a/38`, the `a/38 ± δ` neighbours)
together with covering. The fold reduces all of these to conditions on `Z/2 × Z/19`, i.e.
a **finite mod-38 problem** — the same THM-619-style band problem, one cell down, now with
the parity structure explicit.

## Farey map near the second value

The second value `σ₂ = 2/25`. Its Farey neighbours: `1/13` (below — the gap floor, with
`3/38` the mediant *inside* the cell), and `3/37`, `1/12` (above). **Structured covering
`38`-merge sets attain `3/37` (near-AP, tightest above `2/25`) or `2/25` exactly (the
sporadic `{1,2,3,5,7,8,9,10,11,12,17,19}`), never the cell interior `3/38`** — the
attainment quantizes onto the cell *boundary*, as THM-622 predicts.

## Status

- **Verified/tool (Lean-able, like the S11 atom):** the CRT-fold of composite Farey-cell
  voids to small-modulus multi-color avoidance; the `c≤6` reducibility table; the
  `M=3/38 ⟹ pinned ±1 mod 19` characterization.
- **Not closed:** the `c=3, q=38` void itself (the covering + max-attained contradiction
  on `Z/2 × Z/19` + competitors — mac-mini's covering-pigeonhole lane, now folded).
- **Lead:** the prime cores (`q = 89, 101, 103, …`, `c ≥ 7`) are the irreducible residue
  problems; the composite ones (`c ≤ 6`) should be dispatched by the fold + a bounded
  small-modulus check.

Connects S11 (parity-split), mac-mini THM-622 (Farey-cell reduction), THM-621/HYP-4177
(the ladder attains only boundaries), kps HYP-4167 (Jain–Kravitz relative-spectrum voids).
Script: exploration in the S16 session log.
