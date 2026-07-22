# Each prime is its Paley tournament — a periodic table of 2, 3, 5, 7, 11 for LRC(14)

*boxeph-2026-07-21-S215. Owner: come to understand primes 7, 5, 11 as well as we do 2 and 3 — through the
standing lens "a set of pairwise relations is a tournament." Builds on THM-1830 (Paley-7 spectrum), opus-S434
(the Paley/symmetric-intransitive pole), the Eisenstein/Φ₆ argmax, my S206 (Fibonacci foil), S212 (Gauss-sum
Euler index), S213 (chirality), S214 (rank-11 AP-core). All verified in
`04-computation/understanding_primes_via_paley_tournaments_boxeph_S215.py`.*

## The unifying law: prime `p` IS its Paley object

Under the tournament lens, the canonical "relation structure" on `p` things is the **Paley construction**:
`i → j` iff `j − i` is a quadratic residue mod `p`. Its character splits the primes cleanly (verified):

- **`p ≡ 3 (mod 4)` → a Paley TOURNAMENT** (antisymmetric, self-converse, vertex-transitive, doubly-
  regular). Its adjacency spectrum is `(p−1)/2` and `(−1 ± i√p)/2`, each with multiplicity `(p−1)/2` —
  `char_A = (x − (p−1)/2)·(x² + x + (p+1)/4)^{(p−1)/2}`, the quadratic having discriminant **`−p`**. The
  imaginary part `√p/2` is exactly the **quadratic Gauss sum `i√p`**.
- **`p ≡ 1 (mod 4)` → a self-complementary Paley GRAPH** (undirected, since `−1` is a residue). Real
  spectrum `(−1 ± √p)/2`; Gauss sum **`√p` (real)**.
- **`p = 2` → the involution itself** — `Φ₂ = x+1`, the antipode / reversal / **chirality** (S210–S213)
  that acts on all the above.

Verified: `3,7,11` are Paley tournaments with `char_A = (x−1)(x²+x+1)`, `(x−3)(x²+x+2)³`, `(x−5)(x²+x+3)⁵`
(disc `−3,−7,−11`); `5,13` are self-complementary graphs; the Gauss sums satisfy `|g_p|²=p` with `i√p` for
`3,7,11` and real `√p` for `5,13`. So **each prime's "tournament" is fixed by `p mod 4`, and its spectrum is
its Gauss sum.** That is the same depth we have for 2 (the involution) and 3 (the 3-cycle) — now for all of
them, from one construction.

## The periodic table for LRC(14) = 2·7

Each prime's Paley personality fixes its role in the problem:

| `p` | `p mod 4` | Paley object | Gauss sum | role in LRC(14) |
|---|---|---|---|---|
| **2** | — | the **involution** | `i√2` | chirality / reversal `ι:t↦1−t` (S212/S213); `Φ₂` = the antipode |
| **3** | 3 | Paley-3 = **the 3-cycle** | `i√3` | the **atom** (`x²+x+1` = Eisenstein `ω`); the **argmax** `Φ₆(14)=183`, `t*=14/183` |
| **5** | 1 | self-comp. **graph** | `√5` | the **golden foil**: `ℚ(√5)` = Fibonacci = the LRC *loosest/safest* set (S206); Bonferroni depth |
| **7** | 3 | **Paley-7** tournament | `i√7` | the **apex** (`14=2·7`): `Φ₇,Φ₁₄` carry the hardness; `𝔽₇[C₁₄]=𝔽₇[X]/(X±1)⁷`; my S212 Euler index `i√7`; cap field `ℚ(cos 2π/7)` |
| **11** | 3 | **Paley-11** tournament | `i√11` | the **rank** prime: the rank-11 AP-core / relation code (S214); *scarce* (only multiple ≤14 is 11) → a forced, rigid speed |

Verified pieces of the table: `x¹⁴−1 ≡ (x−1)⁷(x+1)⁷ (mod 7)` (the THM-2043 Frobenius collapse); `2cos(2π/7)`
is a root of the cubic `x³+x²−2x−1` (disc `49=7²`); `Φ₆(14)=183`; `11` has a unique multiple `≤14`.

## Reading the table

- **2 is the reversal that pairs everything.** It is not one of the "objects" but the `ℤ/2` chirality
  acting on them — the even/odd, self-converse, mirror-pair machinery of S210–S213. `Φ₂=x+1` is the `−1`
  evaluation, the antipodal Lefschetz side.
- **3 and 7 are the two `3 mod 4` primes that *build* `14`.** `3` is the **atom** — Paley-3 is literally the
  3-cycle, its spectrum the Eisenstein `ω`, and `Φ₆` (the Eisenstein cyclotomic) sets the extremal argmax
  `t*=14/183`. `7` is the **apex** — Paley-7's `i√7` is the odd-equivariant Gauss-sum index of my S212 Euler
  branch, `Φ₇/Φ₁₄` localize LRC(14)'s hardness, and `𝔽₇[C₁₄]` collapses by Frobenius (THM-2043). Both live
  at the **symmetric-intransitive pole** (opus-S434) — the "hot," maximally-symmetric tournaments — the
  *opposite* end of the reify ladder from the transitive AP nullcone vertex.
- **5 is the `1 mod 4` outlier — the golden foil.** No Paley *tournament*; a self-complementary *graph*
  whose Gauss sum is the **real** `√5`, i.e. `ℚ(√5)`, the golden field. That is exactly why **Fibonacci is
  the LRC foil** (S206): the `1 mod 4` primes give real (non-oscillating) Gauss sums / golden geometry, the
  *loosest*, safest speed sets — the anti-extremal, opposite the anti-golden Eisenstein (`Φ₆`, prime 3)
  extremal. The `2/3/7` (Eisenstein) side is where LRC(14) is tight; the `5` (golden) side is where it is
  slack.
- **11 is the second `3 mod 4` scarce prime — the rank.** Paley-11 (`i√11`) is the prime whose *dimension*
  matters: the rank-11 AP-core relation code (S214), and it is scarce (`11` is its only multiple `≤14`), so
  it forces a rigid coordinate — the "one relation short" of the rank `11→12` descent (codex THM-2052).

## The one-line understanding

> **A prime `p` is understood, in this problem, as its Paley tournament (`p≡3 mod 4`: `3,7,11`, spectrum
> `i√p`) or self-complementary graph (`p≡1 mod 4`: `5,13`, spectrum `√p`); `2` is the reversal that pairs
> them. The Gauss sum `√p`/`i√p` is the Paley spectrum, and `p mod 4` decides tournament-vs-graph, tight-vs-
> slack, Eisenstein-vs-golden.** For `LRC(14)=2·7`: `7` is the apex (`i√7`), `3` the atom/argmax (Eisenstein),
> `5` the golden foil (`√5`=Fibonacci), `11` the rank, `2` the chirality — a coherent periodic table where 5,
> 7, 11 now sit as firmly as 2 and 3.

Honest scope: all rows are verified/classical facts (Paley spectra, Gauss sums, the cyclotomic/field data,
the Frobenius collapse); the contribution is the **synthesis** — one Paley construction giving each small
prime a tournament personality via `p mod 4`, and the mapping of that personality onto LRC(14)'s tight
(Eisenstein 2/3/7) vs. slack (golden 5) structure, its apex (7), and its rank (11) — tying together S206
(foil), S212 (`i√7` Euler index), S213 (chirality), S214 (rank-11), THM-1830/opus-S434 (Paley pole).

Links: HYP-8860, THM-1830, THM-2043,
[[the-rank11-ap-core-is-the-achiral-vertex-where-the-rank-or-euler-frontier-meets-boxeph-S214]],
[[what-an-lrc14-disproof-must-be-and-why-fibonacci-is-the-foil-boxeph-S206]].
