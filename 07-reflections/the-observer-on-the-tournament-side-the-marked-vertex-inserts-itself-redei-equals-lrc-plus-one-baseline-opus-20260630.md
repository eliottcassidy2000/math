# The observer on the tournament side: the OCF's marked vertex IS the observer — it inserts itself into each Ham path with count inshat = 1 + 2·#(its directed 3-cycles), the "+1 baseline + even odd-cycle correction" that is EXACTLY the LRC observer's Farey escape n/(n(n−1)+1); deleting the observer changes H by an even 2·Σμ(C) over the odd cycles THROUGH it (parity preserved ⇒ Rédei), and the transitive class (H=1, no odd cycles) is the pure observer baseline = the metagraph's identity/origin; so Rédei (H odd) and LRC (escape > 1/n) are ONE principle — the observer is irreducibly present (+1), with an even correction (2× its odd-cycle/blocking interaction)

*opus-2026-06-30. Owner: apply the observer lens to the tournament/metagraph side (the marked vertex and
its odd cycles). It lands perfectly — the OCF's vertex-deletion IS the observer, and its "+1" is the SAME
atom as the LRC observer's Farey hair. The two halves of the project are one observer principle.*

## The OCF's marked vertex IS the observer (computed)
The project's vertex-deletion / OCF setup (definitions.md, Claim A) is literally the observer abstraction.
The marked vertex `v` is the OBSERVER; it **inserts itself** into each Hamiltonian path `P'` of `T−v`:
- **the observer's VIEW of a path** = the **signature** `s` (binary: `s[j]=1` iff `v→P'[j]`, else `0`);
- **the observer's self-insertion count** `inshat(v,P') = 1 + 2·#{Type-II} = ODD` (Lemma 5.3), where
  **Type-II positions = the observer's directed 3-cycles** `(v,a,b)` with `(a,b)` consecutive in `P'`;
- **`H(T) = Σ_{P'} inshat(v,P')`** (verified: rotational `T₅`, `H=15 = 3+3+5+1+3` over the 5 paths of
  `T−v`), and **`H(T) − H(T−v) = 2·Σ_{C∋v} μ(C)` is EVEN** (verified `=0` transitive, `=10` rotational).
> So the observer reads each path through its signature, and **inserts itself an odd number of ways: `1`
> (a baseline) plus `2×` its own 3-cycles.** Deleting the observer perturbs `H` only by an EVEN amount (its
> odd cycles, weighted `2`), so **parity is preserved — inducting down to one vertex (`H=1`) gives Rédei.**

## The "+1 baseline" is the SAME atom as the LRC observer's Farey hair
Both observers carry an **irreducible `+1` baseline plus an even correction `2×(interaction)`**:
| | the observer's invariant | the `+1` (irreducible) | the `2×` correction |
|---|---|---|---|
| **LRC** | escape `M = n/Φ₆ = n/(n(n−1)+1)` | the `+1` in `Φ₆=n(n−1)+1` (Farey hair / antipodal killer) | the blockings `n(n−1)` (the runners returning) |
| **tournament** | self-insertion `inshat = 1 + 2·#Type-II` | the `+1` (the irreducible self-insertion) | `2× #(its 3-cycles)` (the odd cycles through it) |
> **`Rédei` and `LRC` are ONE principle: the observer is irreducibly present.** Its local invariant can
> never be trivial — there is always the `+1`. In the LRC that `+1` makes the **escape POSITIVE** (`M>1/n`,
> the observer always clears its blockings by a hair). In the tournament that `+1` makes the **count ODD**
> (`H` odd, the observer always inserts once more than an even number). Same atom, two faces: a positive
> escape and an odd count are both "the observer's `1` survives." The even correction is, on both sides,
> `2×` the observer's interaction with the **odd-cyclic structure** — the LRC's harmonic blockings (the
> runners' returns) and the tournament's directed 3-cycles (the odd cycles through `v`). The factor `2` is
> the antipodal pairing on both sides (`{±1}` blocking pairs; the in/out signature flip of a 3-cycle).

## The metagraph is the observer's view from the transitive baseline
- **The transitive class has `H=1` and NO odd cycles** — it is the **pure observer baseline**, the metagraph's
  identity/origin class (the `0` of the H-gradient).
- **The H-gradient is the odd-cycle accumulation:** every other class is `H = 1 + 2·(its odd-cycle
  collection)` (the OCF, `H=I(Ω,2)=1+Σ 2^j α_j`). So **the merged metagraph `G_n/Z₂` IS the observer's view
  of tournament space from the transitive origin** — distances/edges measure how many odd cycles separate a
  class from the baseline. (This is why the principal line runs from the transitive `H=1` outward, and why
  the H→1 corner binds at the 3-cycle = the LRC doublet's mirror, mac-mini HYP-3585: the minimal odd cycle
  is the observer's first correction on both sides.)

## What the observer lens unifies (both halves, one frame)
| object | the observer | its local view = the invariant |
|---|---|---|
| **LRC** | the origin (deleted speed-0 runner) | escape `1 + 2×blockings`-form, positive by a hair |
| **tournament** | a marked vertex `v` | self-insertion `1 + 2×(its 3-cycles)`, odd |
| **metagraph** | the transitive class `H=1` | the H-gradient `= 1 + 2×(odd cycles)`, the origin outward |
| **Dirac comb** | the missing tooth | the gap (the `+1` that is the observer itself) |
> **One principle across the project: the observer's invariant is `1` (irreducibly present) + `2×` its
> interaction with the odd-cyclic structure.** Rédei = the `1` is odd; LRC = the `1` is a positive escape;
> the metagraph = the `1` is the transitive origin. The whole project is the study of this single `+1` and
> its even odd-cycle correction, read from the marked point.

## Status
- **Computed/verified (opus):** the OCF vertex-deletion = the observer (marked `v`); `inshat=1+2·#Type-II`
  (odd) with Type-II = the observer's 3-cycles; `H(T)=Σ inshat`; `H(T)−H(T−v)` even (Rédei via parity); the
  transitive observer = baseline `H=1`.
- **The unification:** LRC escape (`+1` Farey hair) and tournament insertion (`+1` self-insertion) are the
  SAME irreducible-`+1`-plus-even-odd-cycle-correction atom; Rédei (odd) and LRC (positive escape) are one
  observer principle. The metagraph = the observer's view from the transitive origin.
- **Open (unchanged):** the LRC realizability node (observer's escape = the convergent for `n≥7`) and the
  full OCF→LRC bridge; but they now share one lens — the marked point's irreducible `1`.

Related: the-observer-abstraction (the LRC observer), the-cyclic-kershner-attack (the grid/empty-tooth);
definitions.md (OCF/Claim A), CONJ-001 (OCF proved), mac-mini HYP-3585 (3-cycle ↔ doublet); the-isomorphism-
class-graph, everything-is-the-triangle; OPEN-Q-039/108.
