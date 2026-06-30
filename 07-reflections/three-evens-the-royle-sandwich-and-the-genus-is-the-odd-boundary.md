# Three "evens," the Royle sandwich, and the genus is the odd boundary: Eulerian = bulk, blue = the obstruction

*klein-2026-06-29-S12. Asked to understand tournament/even-graph equinumerosity and its particularities, to find past Eulerian-section work on the merged metagraph, and to work the genus-dimensional column creatively. The three converge: the genus is the odd boundary, and codex-S675b already named it.*

## "Even graph" is three different things (the equinumerosity particularity)

There is no single tournament/even-graph equinumerosity — because "even graph" is polysemous (opus-S260),
and the three meanings give three Burnside sums on the edges of `K_n`, equal only at `n=3`:

| count | meaning | Burnside form | n=3..7 |
|---|---|---|---|
| **A000088** | ALL graphs | `(1/n!) Σ_{all σ} 2^{E_orb(σ)}` | 4,11,34,156,1044 |
| **A000568** | tournaments = **Royle-even** | `(1/n!) Σ_{|σ| odd} 2^{E_orb(σ)}` (= signed cycle index `P_n(1)`, THM-587) | 2,4,12,56,456 |
| **A002854** | degree-even = **Eulerian** = cycle space | `(1/n!) Σ_{all σ} 2^{E_orb − rank_{GF2}(deg map)}` | 2,3,7,16,54 |

Verified (`even_graph_equinumerosity_three_counts_klein.py`). The **particularities**:

1. **The tournament equinumerosity is with ROYLE-even, not Eulerian.** Tournaments `= A000568 =` Royle-even
   graphs (Devillers–Freedman–Glasby–Praeger–Royle 2022) `=` the **odd-order restriction** of the all-graph
   Burnside `=` the signed cycle index `P_n(1)` (THM-587; `|σ| odd ⟺ all cycle lengths odd ⟺ the
   orientation-preserving / all-odd-cycle Burnside). Degree-even/Eulerian (A002854) is a *different* count.
2. **The Royle/odd-order sandwich:** `A002854 ≤ A000568 ≤ A000088` — **Eulerian ≤ tournaments ≤ all graphs**.
   Tournaments are the *odd-order slice* of the full edge-space; Eulerian is the *cycle subspace quotient*.
   They are two different reductions of `GF(2)^{C(n,2)} mod S_n`, not the same one.
3. **The GF(2) parity obstruction is real.** The Eulerian count needs the `GF(2)` rank of the degree map
   (`E=Cut⊕Cycle iff n odd`, opus-S260); the naive `E_orb−V_orb+1` is wrong at even `n` (gives 4,19 instead
   of 3,16 at `n=4,6`). The cut/cycle entanglement at even `n` is the same parity subtlety the whole project
   keeps meeting.
4. **All three `= 2` at `n=3`** — the classic polysemy trap (they look equinumerous at the smallest case;
   run the persistence test, they diverge at `n=4`).

So the clean statement: **tournaments are the signed/odd-order object; Eulerian graphs are the cycle space;
the even-graph "shadow" of a tournament (`T_cycle = (I+L)T mod 2`, odd `n`) is its cycle-space/bulk part.**

## The Eulerian sections of the merged metagraph (codex-S675b, found)

The past work the question pointed at is `merged-line-parity-even-odd-s675b` (codex-S675b). In the corrected
complement-line lift inside the fixed-base tiling cube `Q_m`:

> **black is always an even graph in the strong Eulerian sense (boundary-zero over `F_2` = in the cycle-space
> kernel); blue is the odd-degree boundary.** Verified through `n=7`. Mechanism: complement pairs tilings
> `{x, C(x)}`; grid-reflection-fixed tilings = the blue support; over an SC class `odd−odd` (even), over an
> NS merged pair `odd+odd` (even) ⟹ **black Eulerian**; the blue endpoint count is odd ⟹ **blue =
> odd-boundary coset**. And the punchline S675b already wrote: *"find a line lift where the black side is an
> Eulerian certificate and the blue odd boundary marks the exact floor atoms."*

"Odd is not bad-even; **odd is even plus a named boundary**." The black Eulerian sits in the cycle-space
kernel; the blue is an affine coset of the same cycle space, with boundary one on the SC support.

## The genus is the odd boundary (the creative column)

Now the genus-dimensional global column from the last sessions snaps onto this. The modular split is
**Eisenstein (bulk) ⊕ cusp form (obstruction)**; `dim(cusp) = genus` (HYP-3587). Match it to S675b and the
equinumerosity:

| BULK / local / black / even | OBSTRUCTION / global / blue / odd |
|---|---|
| Eisenstein series | cusp form `f_14` |
| Eulerian = cycle-space kernel (black) | odd-boundary coset (blue) |
| boundary-zero over `F_2` | "even **plus a named boundary**" |
| the even-graph shadow of the tournament | the part the shadow misses |
| `dim = ` (cycle space) | **`dim = genus`** |

So the **genus is the dimension of the odd boundary** — the number of "named boundary atoms" (blue) beyond
the Eulerian/even/cycle-space bulk (black). The cusp form `f_14` IS the blue odd-boundary; S675b's *"blue
odd boundary marks the exact floor atoms"* is, in the modular chart, *"the cusp form is the floor
obstruction."* They are the same statement, found independently (codex 2026 on the metagraph; klein S10/S11
on the modular curve).

For **LRC(14)**: `genus = 1` ⟹ **one** named-boundary atom = **one** blue floor atom = the **doublet** (the
binding `2`-residue core, THM-578, `4cos^2(3pi/7)`) = the single cusp form `f_14`. For `genus 0` (LRC(6)):
zero blue atoms — the boundary is empty, the black Eulerian certificate alone closes it (no cusp form).
For `genus ≥ 2`: several blue atoms, irreducibly complex (HYP-3587). The genus literally counts the blue
floor atoms.

## The unification

Three threads, one statement:
- **equinumerosity**: tournament `=` signed/odd-order object `=` (even-graph cycle-space bulk) `+` (the
  odd/oriented structure); Eulerian `=` the cycle-space bulk alone. The gap `A000568 − A002854 =
  0,1,5,40,402` is the obstruction-bearing surplus.
- **Eulerian sections** (S675b): black `=` cycle-space/even bulk; blue `=` odd-boundary `=` the floor atoms;
  "odd `=` even `+` named boundary."
- **genus column**: bulk `=` Eisenstein/Eulerian/black; obstruction `=` cusp form/blue/odd-boundary;
  `dim(obstruction) = genus =` #blue floor atoms; LRC(14) `= 1` (the doublet).

The master dichotomy gains its homological row: **the local/global split is the even/odd split is the
cycle-space-kernel/boundary-coset split is the Eulerian-black/blue split is the Eisenstein/cusp-form split,
and the genus is the dimension of the odd-boundary obstruction.** Everything computable (the even/Eulerian
cycle-space bulk, the metagraph, `CV(H)`) is black; the one missing thing is the blue boundary — the genus
counts it, and for `N=14` it is a single doublet.

See HYP (this session), [[merged-line-parity-even-odd-s675b]] (codex, the black/blue Eulerian split),
[[even-graphs-through-the-metagraph]] (the three "evens", opus-S260),
[[the-permanent-is-the-unsigned-twin-of-the-spectrum]] (signed/unsigned),
[[the-genus-is-the-local-global-gap-and-the-one-master-dichotomy]] (HYP-3587), THM-587 (signed cycle index),
THM-578 (the doublet), HYP-3586 (cusps=Klein, hardness=genus).
