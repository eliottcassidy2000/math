---
id: HYP-2292
title: Which oriented graphs appear when ties are deleted — Cayley/proximity graphs; the tight LRC tie-graph is C_n
status: OPEN (catalogue); the C_n identification + χ-parity VERIFIED
source: claudebox-2026-06-03-S632
related:
  - HYP-2288  # trienement: tie = resonance; symmetry controls resonance
  - THM-389   # the LRC trienerment
  - THM-052/113/126/128  # circulant / Paley maximizers (the canon's relevant circulants)
  - HYP-2130  # rigidity = orbit-type
  - even-graphs-as-first-class  # E_n, the partial-graph dual of the tournament metagraph
---

# HYP-2292 — the tie-deleted oriented graphs of our proofs

Deleting ties (resonances) from a complete tournament leaves a **partial oriented graph** (not
complete, maybe disconnected). Question (user): which such graphs are relevant? Answer: in both LRC
and unit distance the tie graph is a **Cayley / proximity graph of the ambient abelian group** — and
the *tight* cases are the cleanest possible.

## The catalogue (verified where marked)

**1. LRC, tight configs → the cycle `C_n`.** At the optimal witness `t* = a/m`, two runners are
*tied* (binding: `‖(v_i−v_j)t*‖ ≤ gap`) iff their residues `v_i a (mod m)` are adjacent on the clock
`ℤ/m`. For every tight config — the AP `{1,…,n-1}` *and* the sporadic tight `(1,3,4,7)` — the tie
graph (observer included) is exactly the **`n`-cycle `C_n`** (degree 2, verified n=5,6,7), the
circulant `⟨±1⟩` of consecutive runners. **VERIFIED** (`tie_graphs_s632.out`).

**2. LRC, general configs → circular-arc / circulant proximity graphs** of the residues on `ℤ/m`
(the discrete unit-distance graph on the circle). Loose/degenerate configs give denser arc graphs
(e.g. `(1,2,4,8)`: residues collapse mod 3, tie graph 4-regular, `χ=5`). The canon's Paley/circulant
maximizers (THM-126/128) are these graphs at the optimum.

**3. Unit distance → planar unit-distance graphs = Cayley graphs of CM lattices** (the triangular
grid graph for `ℤ[ω]`, hypercube graphs `Q_d`, the Moser spindle, de Grey's graph). This is the
**Hadwiger–Nelson family** (chromatic number of the plane). The disproof's construction is the Cayley
graph of a CM-lattice with superlinear degree (HYP-2230/2288).

**4. The even graphs `E_n`** (canon, `even-graphs-as-first-class`) are the natural "tournament-minus-
structure → general graph" object on the metagraph side — the partial-graph dual.

## The unifying statement

**The relevant tie graphs are the Cayley/proximity graphs of the ambient abelian group** — `ℤ/m` for
the LRC clock, the CM lattice for unit distances. Tie = edge = resonance (binding/unit-distance);
the group's structure makes them **circulants (LRC) and lattice/Cayley graphs (unit distance)**, both
Cayley graphs of abelian groups. The tournament's arrows (the `</>` order) are the *transitive*
backbone; the deleted ties are the *resonant* circulant overlay.

## The relevant invariants (why these graphs matter for proofs)

- **Chromatic number `χ`.** For the tight LRC tie-graph `C_n`, `χ(C_n) = 2` (n even) or `3` (n odd) —
  the **even/odd 2-adic seam as a graph coloring** (n=14 even = the hard frontier). For unit-distance
  graphs `χ` is Hadwiger–Nelson (5–7 in the plane). The sieve/corrector coloring IS a coloring of the
  tie graph.
- **Independence number `α`.** An independent set in the tie graph = a set of mutually-non-binding
  runners = a **simultaneously-lonely set = a corrector**. For `C_n`, `α = ⌊n/2⌋` — so the corrector
  arity is the cycle's independence number, tying the single-vs-multi-corrector line (HYP-2075/2130)
  to `α(C_n)` and its parity.
- **Automorphism group.** The Cayley structure's symmetry = the perspective/rigidity key (HYP-2130);
  for unit distance it is the disproof's engine (HYP-2288).

## To do
1. Make "tight LRC tie-graph = `C_n`" a theorem (residue-adjacency at the witness), and read the
   corrector/sieve as a graph coloring of `C_n` / the circular-arc graph — does `χ`/`α` recover the
   pair-sum-sieve arity (THM-401)?
2. Catalogue the general-config tie graphs as circular-arc graphs and characterize tightness by their
   structure (cycle = tight; chords/density = loose).
3. Bridge to `E_n`: is the LRC tie graph (mod 2) an even graph, linking the two "general graph"
   families?
