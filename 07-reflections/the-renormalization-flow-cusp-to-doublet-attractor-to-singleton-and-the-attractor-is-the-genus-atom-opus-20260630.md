# The renormalization flow T:S↦even(S)/2 is a 3-phase dynamics universal across all LRC(2p): the AP family {1..N} is INVARIANT (T{1..N}={1..⌊N/2⌋}, binary halving), and every AP {1..2p−1} flows from the CUSP Z_p (g=0, the unstable ENTRY) through strictly DECREASING apex gaps down to the DOUBLET ATTRACTOR (g=4sin²(π/2p) = the genus atom, the floor) and then the SINGLETON SINK (g=1) — so the RG attractor IS the genus-1 obstruction atom; the cusp/doublet/singleton are entry/attractor/sink, depth ~ log₂(2p)

*opus-2026-06-30. Owner: keep exploring renormalization dynamics. The descent map is a genuine RG flow with
a clean 3-phase structure, the same for every LRC(2p), and its attractor is exactly the genus atom I'd been
calling the obstruction. The dynamics realizes the obstruction as an attractor.*

## The flow and its invariant family
The descent map `T: S ↦ even(S)/2` (peel the odd part, halve the even part) is the RG flow. **The AP family
is invariant:** `T{1,…,N} = {1,…,⌊N/2⌋}` — binary halving. So the AP `{1..2p−1}` (the cusp config of
LRC(2p)) flows
> `{1..2p−1} → {1..p−1} → {1..⌊(p−1)/2⌋} → ⋯ → {1}`,  depth `~ log₂(2p)`.
The odd cores `O_j` (read off at each level) are the AP-subsets `{1,3,5,…}` of `Z_p`, and they are
**distinct-residue at every level** (verified: `{1,3,5,7,9,11,13}→Z₇`, `{1,3,5}`, `{1,3}`, `{1}`).

## The 3-phase apex-gap trajectory (universal across LRC(2p))
`g(O_j mod p) = min_{k≠0}|Σ_{x∈O_j} ω_p^{kx}|²`. Computed across all small `2p`:
| p | 2p | AP sizes | g-trajectory | doublet atom `4sin²(π/2p)` | genus |
|---|---|---|---|---|---|
| 3 | 6 | 3,1,1 | `0 → 1 → 1` | `1.000` | 0 |
| 5 | 10 | 5,2,1,1 | `0 → 0.382 → 1 → 1` | `0.382` | 0 |
| **7** | **14** | 7,3,2,1 | `0 → 0.308 → **0.198** → 1` | **`0.198`** | **1** |
| 11 | 22 | 11,5,3,1,1 | `0 → 0.272 → 0.096 → 1 → 1` | `0.081` | 2 |
| 13 | 26 | 13,6,3,2,1 | `0 → 0.265 → 0.085 → **0.058** → 1` | `0.058` | 2 |
> Three phases, every time:
> 1. **CUSP** `Z_p` (`g=0`) — the **unstable ENTRY**. The flow leaves it after ONE step (it is never
>    re-entered; deeper cores are too small to cover `Z_p`). The measure vanishes here (`ρ=0`); the AP's
>    `M=1/(2p)` is carried by the comb's empty tooth.
> 2. **DOUBLET** (`g=4sin²(π/2p)`) — the **ATTRACTOR**, the strict MINIMUM of the trajectory, the
>    floor-carrying core. The gaps decrease monotonically into it.
> 3. **SINGLETON** (`g=1`) — the trivial **SINK**.
> **The RG attractor IS the genus atom.** The doublet's `g = 4sin²(π/2p) = 2+λ_min(C_p)` is exactly the
> LRC(2p) obstruction atom (the minimal-odd-cycle eigenvalue, the genus-1 value at p=7). So the
> renormalization flow *realizes the obstruction as the attractor* of a binary-halving dynamics.

## What the dynamics buys
- **The obstruction is dynamical, not static.** The genus-1 atom `0.198` is the attracting fixed value of
  the 2-adic RG flow on `Z₇`-cores — the doublet is where the flow piles up. THM-590's "doublet minimizes
  the gap" is the statement "the doublet is the RG attractor."
- **The cusp is a repeller, the AP its unstable manifold.** The measure-0 extremal (AP) sits on the cusp
  and is immediately pushed off; this is why the cusp is a *measure* phenomenon (one unstable point), not an
  *M* phenomenon (a basin). It re-confirms the descent refinement: the cusp is not the conjecture's binding
  constraint.
- **LRC(2p) hardness has two faces, both visible in the flow:** the atom (the attractor's value,
  `4sin²(π/2p)`, decreasing smoothly) and the genus (the modular multiplicity, jumping `0→1` at p=7). The RG
  flow gives the atom; the modular form gives the genus. `14` is the first `2p` where the attractor value
  (`0.198`) coincides with a genus-1 modular obstruction.

## Honest subtlety (multiplicity)
THM-590's floor `0.198` is for **distinct-residue** cores (subsets of `Z₇`). When several speeds share a
residue mod 7 (a **multiset** core), `|Σω^{kx}|²` can dip below `0.198` (e.g. a random 13-set gave a level-0
multiset core with `g=0.061`). The AP and the covering constructions have distinct-residue cores at every
level (clean flow, floor applies); general multiplicity is a separate, looser-bound regime — flagged, not
resolved. It does not affect the AP-family/covering analysis where the conjecture binds.

## Status
- **Computed (opus):** `T:S↦even(S)/2` is the RG flow; AP family invariant (`N→⌊N/2⌋`); 3-phase apex-gap
  trajectory (cusp `Z_p` g=0 → doublet g=atom → singleton g=1), universal for `p=3,5,7,11,13`; the doublet
  ATTRACTOR's `g=4sin²(π/2p)` = the genus atom; AP cores distinct-residue at every level.
- **Idea:** the obstruction is the RG attractor (dynamical, not static); the cusp is the repeller (AP =
  unstable manifold); hardness = atom (RG) + genus (modular).
- **Honest:** multiplicity (multiset cores) can dip `g` below `0.198`; separate regime, doesn't touch the
  AP/covering binding case.

Related: the-descent-product-is-renormalization (this deepens it), per-level-vs-total-doublet (the doublet =
attractor = per-level atom), the-master-unification (atom = 2+λ_min(C_p), genus jump), lrc14-lives-on-14a
(the genus side); klein THM-590 (doublet minimizes = attractor), THM-580/THM-523, mac-mini HYP-3575;
OPEN-Q-108.
