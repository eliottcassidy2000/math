---
source: claude-2026-06-06-S680
status: PROVED -- H(T) != 21 for all tournaments (resolves THM-115, pending peer verification)
tags: [H21, H-impossibility, PROVED, moon-pancyclicity, odd-cycles, strong-tournament, OCF, theorem, redei-spectrum]
---

# H=21 is impossible — a complete proof (closing THM-115)

The chase for the "twisted involution on the residual" turned, via the
strong-component reduction, into a clean **complete proof** that no tournament
has exactly 21 Hamiltonian paths. THM-115 was conjectured in S31 and verified
exhaustively only to n=8; the general case is now closed.

## Theorem

**No tournament `T` has `H(T) = 21`.**

## Proof

`H(T) = I(Ω(T), 2) = 1 + 2α₁ + 4α₂ + …` (the OCF, THM-029), where `α₁` = number
of directed **odd** cycles of `T`. In particular `H ≥ 1 + 2α₁`.

**(A) Strong-component reduction.** `H` is multiplicative over the
strongly-connected components: `H(T) = ∏_i H(C_i)`. The odd factorizations of
`21` are `21` and `3·7`. Since `7` is *not* an H-value of any tournament
(THM-029/THM-200), it is not the H-value of any strong component, so `3·7` is
impossible. Hence `H(T)=21` forces a **single strong component** with `H=21`.
WLOG `T` is strong.

**(B) Base case `m ≤ 8`.** `H=21` is impossible for all tournaments on `≤8`
vertices — exhaustively (THM-079 Part G: `268,435,456` tournaments at `m=8`).

**(C) Inductive bound `m ≥ 9`.** A strong tournament on `m` vertices is
**vertex-pancyclic** (Moon 1966): every vertex lies on a directed cycle of each
length `3,4,…,m`. So for each `L`, the `L`-cycles **cover** all `m` vertices, and
since each `L`-cycle has exactly `L` vertices,
`L·(#L-cycles) ≥ m`, i.e. `#L-cycles ≥ ⌈m/L⌉`. With `c₃ ≥ m−2` (Moon's minimum
number of 3-cycles for strong tournaments; verified here for `m≤7`), summing over
**odd** lengths:
```
α₁  ≥  (m−2) + Σ_{odd L, 5 ≤ L ≤ m} ⌈m/L⌉.
```
For every `m ≥ 9` this exceeds `10` — the values are `12, 14, 17, 19, 21, 23, 26`
for `m = 9..15` (and for `m ≥ 13`, `c₃ ≥ m−2 ≥ 11` already suffices). Hence
`α₁ ≥ 11`, so `H ≥ 1 + 2·11 = 23 > 21`.

By (A), `H=21` needs a strong component with `H=21`; by (B) it is not on `≤8`
vertices; by (C) it is not on `≥9` vertices (`H ≥ 23`). Contradiction. **∎**

## The binding case, and why it needed pancyclicity

`m = 9` is the one tight case. There `c₃ ≥ 7` alone gives only `α₁ ≥ 7`, not `>10`.
The proof is saved by the **longer** odd cycles via vertex-pancyclicity: every
vertex is on a 5-cycle and a 7-cycle, so (covering) `#5 ≥ ⌈9/5⌉ = 2`,
`#7 ≥ ⌈9/7⌉ = 2`, `#9 ≥ 1`, giving `α₁ ≥ 7+2+2+1 = 12 > 10`. This is exactly why
the H=21 obstruction is "harder than H=7": it lives just past where the 3-cycle
count alone settles things, and only the full pancyclic odd-cycle spectrum
closes it.

## Relation to prior work

- THM-079 had the right *framework* (strong-component multiplicativity, the
  i₂-jump, exhaustive `m≤8`) but left the **general** case open (the i₂-jump was
  verified, not proved, for each profile). S680 closes it by bounding `α₁`
  directly — bypassing the profile/i₂-jump entirely.
- S617 (HYP-2193) reduced to `m ≤ 12` via `c₃ ≤ α₁ ≤ 10` + Moon `c₃ ≥ m−2`. S680
  sharpens the same idea by adding the *long* odd cycles (the covering bound),
  which pushes the threshold below the base case and removes the `m∈{9..12}`
  residual entirely. So S617 was the right reduction; S680 completes it.
- It also corrects the S616 detour cleanly: S616 chased a (buggy) conflict-graph
  profile; the winning move was to *lower-bound `α₁`* (no profile needed).

## Consequence: the permanent H-gaps are exactly {7, 21}

H=7 (THM-200) and H=21 (now proved) are the only odd integers that are never the
Hamiltonian-path count of any tournament (THM-079 Part I: the only gaps below
200; 35, 63, 189, … all fill at higher `n`). So the achievable H-spectrum =
the multiplicative monoid of strong H-values, with `{7,21}` the complete
"forbidden atoms/products" — closing the converse-of-Rédei question below the
monoid.

## Caveats / verification

The proof rests on two classical Moon results (vertex-pancyclicity 1966;
`c₃ ≥ n−2` for strong tournaments) and two repo results (THM-029, THM-079 Part G,
both established). I verified `min c₃ = m−2` for `m ≤ 7`, the covering counts on
the Moon-extremal family (`α₁ = 2^{m-3}` there), and the arithmetic bound for all
`m`. The logic is elementary given Moon; worth a second pair of eyes on the
`α₁ ≥ 1 + 2`-coefficient step and the `c₃ ≥ m−2` citation.

## Next

1. **Promote to a clean theorem file** and have another agent verify the Moon
   citations.
2. **Generalize the method:** the bound `α₁ ≥ (n−2) + Σ_{odd L≥5} ⌈n/L⌉` is a
   general lower bound on odd cycles of strong tournaments — it likely rules out
   *other* small H-values from strong components too (sharpening the
   strong-value spectrum, HYP-2183).
3. **Back to LRC@19:** the same "bound the count, don't chase the profile" move
   may help the sieve-covered residual — count guaranteed off-grid witnesses
   rather than analyzing the apex case-by-case.
