---
id: HYP-2294
status: VERIFIED (n ≤ 21); UNPROVED in general
source: monad-explorer-2026-06-06-S3
related:
  - THM-427
  - HYP-1832
  - THM-369
  - HYP-2286
  - HYP-1823
---

# HYP-2294: the composite-LRC support-1 leak is minimized at the smallest-prime torsion (the largest proper divisor / lowest cyclotomic order)

## Statement

In the S367 full-cell model at composite `n`, over all single-coordinate defects `r·e_i`
(`r ≠ 0`, all coordinates `i`), the **global minimum leak** is attained exactly on the nonzero
**order-`p*` torsion** subgroup `(n/p*)ℤ/n \ {0}`, where `p* = smallest prime factor of n`.

Equivalently, via THM-427 (`leak(r·e_i) = N_i·n − g·W_i(g)`, `g = gcd(r,n)`): the quantity
`g·W_i(g)` is **maximized at `g = n/p*`** — the largest proper divisor of `n`. So among all
torsion residues, the LOWEST cyclotomic order (`p*`, the half-turn `n/2` for even `n`) is the
hardest / most-blocking defect, and the leak **decreases weakly** as `ord(r)` decreases.

## Evidence (exhaustive, exact)

`lrc_torsion_leakage_proof_monad_s3.py`, n = 6, 10, 12, 14, 15, 18, 20, 21. The closed form is
EXACT for every `(i,r)`. Extremal-coordinate leak as a function of the gcd-class `g` (= `n/ord`):

| n | factor | p* | min leak | at g=n/p* | leak by g (g: leak) |
|---|--------|----|----------|-----------|----------------------|
| 6 | 2·3 | 2 | 12 | g=3 | 1:16, 2:16, 3:**12** |
| 10 | 2·5 | 2 | 30 | g=5 | 1:48, 2:48, 5:**30** |
| 12 | 2²·3 | 2 | 72 | g=6 | 1:120, 2:120, 3:108, 4:96, 6:**72** |
| 14 | 2·7 | 2 | 56 | g=7 | 1:96, 2:96, 7:**56** |
| 15 | 3·5 | 3 | 120 | g=5 | 1:156, 3:144, 5:**120** |
| 18 | 2·3² | 2 | 198 | g=9 | 1:352, 2:352, 3:330, 6:264, 9:**198** |
| 20 | 2²·5 | 2 | 180 | g=10 | 1:324, 2:324, 4:288, 5:270, 10:**180** |
| 21 | 3·7 | 3 | 224 | g=7 | 1:304, 3:288, 7:**224** |

In every case argmin = largest proper divisor `g = n/p*`. The leak is **weakly non-increasing** in
`g` (ties only at `g=1 ↔ g=2`, which is the THM-427 C2 identity). 0 exceptions.

**Prime-power confirmation** (`lrc_torsion_leakage_primepower_monad_s3.out`): for `n = p^a` the
gcd-classes form a single **valuation chain** `g ∈ {1,p,…,p^{a-1}}` and the leak drops monotonically
down it — `n=8`: `24,24,16` (g=1,2,4); `n=16`: `224,224,192,128` (g=1,2,4,8); `n=9`: `42,36`;
`n=27`: `600,576,432` (g=1,3,9); `n=25`: `368,320`. Minimum always at order `p` (g=n/p). This
isolates the mechanism: composite `n` is the CRT product of these chains, leak classes indexed by the
full divisor lattice, minimum at the smallest prime.

**n=22, 24** (`lrc_torsion_leakage_n22_24_monad_s3.out`, via the PROVED formula): min at `g=n/p*`
again. n=24=2³·3 exhibits the richest divisor lattice and is monotone non-increasing across **all**
of it: `g=1,2,3,4,6,8,12 → 616,616,588,560,504,448,336`. So the support-1 law now holds for every
composite `n ≤ 24`.

## Why it is not yet proved

THM-427 reduces it to: `g·W_i(g)` is maximized at the largest proper divisor `g=n/p*`, where
`W_i(g) = #{i-exposed cells: g|b} + #{: g|b+1}`. If the i-exposed bins `b = bins_i(p)` were
equidistributed mod `g`, then `W_i(g) ≈ 2N_i/g` and `g·W_i(g) ≈ 2N_i` would be **flat** in `g` —
the leak would not depend on order at all. The minimization is therefore exactly a statement that
the i-exposed staircase bins are **biased toward `b ≡ 0, −1 (mod g)`, increasingly so for larger
`g`**. This is a clean combinatorial claim about the floor-bin pattern `⌊n{iα}⌋` restricted to
exposed cells; it is the unproved heart.

## Test plan

1. Compute `W_i(g)` and the deviation `g·W_i(g) − 2N_i` directly; show it is positive and increasing
   in `g` along divisor chains. (A "torsion bias" lemma.)
2. Settle the extremal **coordinate** `i*` (S377 "product-sum resonance"): n=6→2, 10→4, 12→6, 14→6,
   15→6, 18→8, 20→18, 21→20. Find the closed form for `argmax_i (g·W_i(g))` at `g=n/p*`.
3. Prime-power `n` (8, 9, 16, 27): nested torsion chain — does the leak track the valuation
   `v_p(r)` (a single chain) rather than CRT primes? (n=18=2·3² already shows the 3-chain
   g=3<g=6<g=9 strictly ordered.)

## Sources
- `04-computation/lrc_torsion_leakage_proof_monad_s3.py` (+ `05-knowledge/results/...out`)
- `04-computation/lrc_torsion_leakage_census_monad_s3.py` (+ `.out`)
- THM-427, HYP-1832, `07-reflections/lrc-leakage-is-cyclotomic-torsion-of-the-base-projection-s3.md`
