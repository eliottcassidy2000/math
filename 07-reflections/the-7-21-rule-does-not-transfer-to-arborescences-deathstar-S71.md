# The {7,21} rule does NOT transfer to arborescences — and exactly why

**death-star-2026-07-20-S71** (HYP-8545). Owner: how is the Ham-path `{7,21}` forbiddenness (S70)
reflected for the **arborescence count** `A(T)=Σ_r a_r` (out-arborescences, Matrix-Tree, poly-time),
what rules are shared, and how do they differ exactly? Answer: **it does not transfer** — `7` is a
perfectly good arborescence count. The two invariants differ on every structural property that makes
`{7,21}` special, and the arborescence "forbidden set" is a different, infinite, growth-driven object.

## The headline: 7 is achievable, {7,21} is a Ham-path phenomenon only

Arborescence spectra (exhaustive `n≤6`, `a_r=det((D_in−A)` with row/col `r` deleted`)`):
`n=3:{2,3}`, `n=4:{6,7,9,10}`, `n=5:{24,26,28,30,36,40,42,43,48,49,55}`, `n=6:{120,…,333}`.

**`7` appears at `n=4`.** And `21`? It is absent — but only because it falls in the *gap between bands*
(`n=4` maxes at 10, `n=5` starts at 24), a trivial arithmetic absence, not a structural prohibition.
So the Ham-path fact "`H≠7`, `H≠21`" has **no arborescence analog**.

## The three properties that make {7,21} special — and how arborescences fail each

| property | Ham paths `H=I(Ω,2)` | arborescences `A=Σ_r a_r` |
|---|---|---|
| **parity** | **always odd** (Rédei) ⟹ every even auto-forbidden | **both parities** — transitive gives `(n−1)!` (even, `n≥3`), 3-cycle gives 3 (odd). No parity law. |
| **scale** | **n-stable**: small values recur ∀n ⟹ one global spectrum `= odd\{7,21}` | **grows**: `A≥(n−1)!`, min = transitive; each `n` is a band `[(n−1)!,·]` = 1,2,6,24,120,… |
| **algebra** | **multiplicative monoid**: `H(T⊕S)=H(T)H(S)` (ordinal sum) | **not multiplicative**: `A(T₂⊕T₂)=6≠1·1`. Composes *factorially* (tree parent-choices), not by concatenation. |
| **forbidden set** | **finite & structural**: exactly `{7,21}` (7 = unique bad prime; 21 = its composite casualty) | **infinite & growth-driven**: `{4,5,8,11–23,25,27,29,31–35,37–39,41,44–47,50–54,…}` — the gaps between/within `(n−1)!`-bands, permanent because higher bands start higher (`n=7` min `=720`). |
| **complexity** | `#P`-hard | poly-time (Matrix-Tree determinant) — the repo's THM-1580 point |

Every entry is a genuine difference. The **oddness is the crux**: Rédei's parity (via the path-reversal
involution) forces `H` odd, which auto-forbids all evens and concentrates the whole question on two odd
gaps. Arborescences count spanning *trees*, have no such involution, take both parities, and their
counting is a *determinant* — smooth, poly-time, and monotone-in-`n` (min `(n−1)!`), so their spectrum
stratifies into sparse bands instead of a dense stable line.

## What is shared (the real analogies)

1. **The transitive tournament is the unique minimizer of both** — `H=1`, `A=(n−1)!`. (For `H` this is
   Moon/Rédei; for `A` it is the product-of-in-choices `∏(k−1)=(n−1)!`.) Both invariants measure
   "how far from transitive," from opposite mathematical directions.
2. **Both have forbidden values / spectral gaps** — but of opposite character: `H` has a *finite*
   structural obstruction (`{7,21}`), `A` has *infinitely many* arithmetic gaps from factorial growth.
3. **Both count spanning substructures** (a Hamiltonian path is a spanning *out-path*; an arborescence
   is a spanning *out-tree*). The path is the degenerate tree, so `H`'s multiplicativity is the "linear"
   shadow of `A`'s factorial branching.

## The one-sentence answer
`{7,21}` is a **parity+monoid** phenomenon of the odd, multiplicative, `#P`-hard Hamiltonian-path count;
the poly-time arborescence count is even-permitting, factorially-growing, and non-multiplicative, so it
has no `{7,21}` — instead its forbidden set is the infinite family of gaps left by the `(n−1)!` bands
(smallest: `4,5,8,11,…`). The two invariants share only the transitive minimizer and the bare fact of
having gaps.

## Honest status
Exhaustive `n≤6` (arborescences) + the `n≤8` Ham data from S70. The arborescence "permanently forbidden
`<56`" list is rigorous (bands are monotone: `n=7` min `=720`, so nothing `≥56` from a later band drops
below). Not claiming the arborescence spectrum has a *clean* description — the point is precisely that,
unlike Ham's `odd\{7,21}`, it does not.

## Credit
death-star-S70 (the Ham `{7,21}` spectrum this contrasts), klein-S355 / THM-1580 (arborescence vs
Ham-path separation, the poly-vs-#P framing), canon THM-029/079 (`H≠7,21`); Tutte/Matrix-Tree, Rédei.

## Cross-links
S70 (`the-hamiltonian-path-spectrum-is-odds-minus-7-and-21`), THM-1580 (HYP-8415), definitions.md (`I(Ω,2)`),
`04-computation/arborescence_spectrum_deathstar_S71.py`, HYP-8545.
