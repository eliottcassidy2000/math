# Three avenues through the recursive observer lens, and they lock into one GF(2) structure: the EVEN-GRAPH dual is the observer's CYCLE space (dim C(n−1,2) = the staircase tile count m), the SANDPILE/chip-firing is the observer's CUT space with the SINK = the observer (critical group of K_n = (ℤ/n)^{n−2}, order n^{n−2} = #spanning trees, Cayley), and the ZETA/Euler product is the observer's ANALYTIC trace (residue 1 at s=1 = the observer's irreducible 1, ζ(−1)=−1/12 = the triangle regularized); cut⊕cycle = sandpile⊕even-graph = base-path⊕wiggly is the project's GF(2) split, the observer (sink = empty even graph = residue 1) is the shared baseline, and halve (a=÷2) + add-one (b=+1) routed by parity builds all three

*opus-2026-06-30. Owner: take the recursive mindset into the even-graph dual and the zeta/Euler product, and
find a 3rd avenue to explore in tandem. The 3rd that compelled me: chip-firing/sandpile — it has a built-in
SINK (the marked observer) and it is the GF(2)-DUAL of the even-graph cycle space. The three lock together.*

## (1) Even-graph E_n = the observer's CYCLE space
The even subgraphs of `K_n` (every vertex even-degree) ARE the **cycle space** `GF(2)^{C(n−1,2)}`:
| n | edges `C(n,2)` | cut dim `n−1` | **CYCLE dim `C(n−1,2)` = staircase m** | A002854 iso classes |
|---|---|---|---|---|
| 3..7 | 3,6,10,15,21 | 2,3,4,5,6 | **1, 3, 6, 10, 15** | 2, 3, 7, 16, 54 |
> **The cycle-space dimension is exactly `C(n−1,2)` = the staircase tile count `m`** — the triangle again
> (`f·g=T`). The observer is the **empty even graph** (the cycle-space `0`); `b` = XOR-in-a-cycle builds the
> space. (Aside: A002854 starts `2,3,7` — sharing the Sylvester/apex start — then diverges `16,54`; flag as
> coincidence pending a reason.)

## (2) Sandpile/chip-firing = the observer's CUT space; the SINK is the observer
The chip-firing model on `K_n` has a marked **sink** — *literally the observer* (the marked vertex). Its
sandpile/critical group is on the **cut space** (dual to the cycle space):
| n | reduced-Laplacian det = #spanning trees | `n^{n−2}` (Cayley) | critical group |
|---|---|---|---|
| 2..7 | 1, 3, 16, 125, 1296, 16807 | same | `(ℤ/n)^{n−2}` |
> **The sink = the observer; the critical group of `K_n` is `(ℤ/n)^{n−2}`, order `n^{n−2}` = the number of
> spanning trees (Cayley).** Toppling is the chip-firing *descent* (the discrete `a`); the recurrent
> **identity configuration = the observer's baseline**. The cut space has dim `n−1` (the score hierarchy =
> the base-path arcs); it is the GF(2)-orthogonal complement of the even-graph cycle space.

## (3) Zeta/Euler = the observer's ANALYTIC trace
`ζ(s) = ∏_p (1−p^{−s})^{−1}` — the **Euler product is the descent over primes** (`a` = the `p=2` factor, the
2-adic). Two observer signatures:
> **The residue of `ζ` at `s=1` is `1`** — the **observer's irreducible `1`** (the pole/baseline, the same
> `1` as Rédei's parity, the Farey hair, the descent fixed point, the Sylvester Egyptian sum). And
> **`ζ(−1) = −1/12`** = the regularized `1+2+3+…` = the **"total triangle"** (Bernoulli `−B₂/2`; triangular
> GF `= x/(1−x)³`). The graph version is the **Ihara zeta** `∏_{prime cycles}(1−u^{ℓ})^{−1}`, whose Bass
> formula ties the cycle space (the prime cycles = the even-graph side) to the determinant (the trees = the
> sandpile/cut side) — so the zeta is literally the spectrum joining avenues (1) and (2).

## The three lock: cut ⊕ cycle = sandpile ⊕ even-graph, zeta = the spectrum
| | observer's `1` | the descent `a=÷2` | the `+1` `b` | what's the triangle |
|---|---|---|---|---|
| **even-graph (CYCLE)** | empty even graph (cycle `0`) | 2-adic cycle reduction | XOR-in a cycle | cycle dim `C(n−1,2)` = `m` |
| **sandpile (CUT)** | recurrent identity (the sink) | chip-firing toppling | mark the sink | cut dim `n−1` (scores) |
| **zeta (ANALYTIC)** | residue `1` at `s=1` | Euler factor at `p=2` | the `+1` in the pole | `ζ(−1)=−1/12` (reg-total) |
> **`cut ⊕ cycle` is the project's GF(2) split** — `base-path ⊕ wiggly` (CLAUDE.md) `=` `sandpile ⊕
> even-graph`. The **even-graph is the cycle half** (the wiggly, the staircase `C(n−1,2)`), the **sandpile is
> the cut half** (the base-path, the scores `n−1`, the sink = observer), and the **zeta is their analytic
> spectrum** (the Euler product, residue `1` = the observer, `ζ(−1)=−1/12` = the triangle). All three are
> built by **halve (`a`) + add-one (`b`) routed by parity (GF(2))**, and all three share the **observer's
> `1`** as baseline — the sink, the empty even graph, and the pole's residue are the same marked point.

## What this buys (the mindset extended)
- **The GF(2) cut/cycle split has names on both halves now:** sandpile (cut, the sink-observer, `n^{n−2}`
  trees) and even-graph (cycle, the staircase `m`), with the **zeta (Ihara/Bass) as the bridge** — prime
  cycles (even-graph) and trees (sandpile) are the two terms of one zeta. The recursive lens unifies the
  combinatorial halves and the analytic whole.
- **The observer's `1` is now quadruply realized:** Rédei parity · LRC Farey hair / Sylvester Egyptian sum ·
  sandpile recurrent identity (the sink) · `ζ` residue at `s=1`. One irreducible baseline, four faces.
- **The triangle recurs as the cycle dimension** `C(n−1,2)=m` (even-graph) and as `ζ(−1)=−1/12` (the
  regularized total) — the staircase is both the GF(2) cycle rank and the zeta-regularized `Σn`.

## Status
- **Computed/verified (opus):** even-graph cycle dim `= C(n−1,2)` = staircase `m`; sandpile critical group
  `(ℤ/n)^{n−2}`, order `n^{n−2}` = trees (Cayley, verified n=2..7); `ζ` residue `1` at `s=1`, `ζ(−1)=−1/12`
  (`−B₂/2`); the GF(2) cut⊕cycle = sandpile⊕even-graph = base-path⊕wiggly identification.
- **The synthesis:** three avenues, one structure — cycle (even-graph) + cut (sandpile, sink=observer) =
  GF(2) halves, zeta = their spectrum (Ihara/Bass bridges trees↔cycles); the observer's `1` (sink = empty =
  residue) is the shared baseline; `a` (halve) + `b` (add-one) + parity builds all three.
- **Suggestive (to pin):** A002854's `2,3,7` start vs Sylvester; the precise Ihara-zeta of `G_n`/`K_n` as
  the project's "metazeta"; the observer's `1` as a single object across all four realizations.

Related: the-recursive-mindset (Sylvester/Egyptian), the-functional-decomposition (a/b/triangle), the-
observer-abstraction + observer-on-the-tournament-side; even-graphs-as-first-class, even-graphs-through-the-
metagraph (the cycle space), the-zeta-function-and-the-ocf-read-complementary-halves; OPEN-Q-039/108.
