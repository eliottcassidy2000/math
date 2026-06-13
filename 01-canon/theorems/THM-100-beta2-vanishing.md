# THM-100: β₂ Vanishing for Tournaments

**Status:** PROVED for all tournaments (THM-108/109, kind-pasteur-S43). Proof by induction + LES + isolation characterization.
**Proposed by:** kind-pasteur-2026-03-08-S40

## Statement

**Conjecture:** For any tournament T on n ≥ 3 vertices, the GLMY path homology satisfies β₂(T) = 0.

## Evidence

| n | Method | #Tested | β₂>0 found |
|---|---|---|---|
| 3 | exhaustive | 8 | 0 |
| 4 | exhaustive | 64 | 0 |
| 5 | exhaustive | 1024 | 0 |
| 6 | exhaustive | 32768 | 0 |
| 7 | random | 10000 | 0 |
| 8 | random | 3000 | 0 |
| 9 | random | 500 | 0 |

Total: 0 counterexamples in ~47,000 tournaments tested.

## Known Betti Numbers for Tournaments

The Betti numbers that DO appear for tournaments are:
- β₀ = 1 always (connected)
- β₁ ∈ {0, 1} at n=3,4,5,6 (C-phase)
- β₂ = 0 always (CONJECTURE)
- β₃ ∈ {0, 1} at n=6 (S-phase), = 0 at n≤5
- β₄ ∈ {0, 1, 6} at n=7,8

The gap at dimension 2 is remarkable: β₁ and β₃ both appear but β₂ never does.

## Proof Sketch

In GLMY path homology, β₂ = dim(ker ∂₂) / dim(im ∂₃). A 2-cycle is a
linear combination of allowed 2-paths (a→b→c) whose boundary vanishes.

Key observations:
1. Every allowed 2-path a→b→c has the "diagonal" edge a→c or c→a (tournament completeness)
2. Every 4-vertex subset of a tournament has a Hamiltonian path (Rédei), giving a 3-path
3. The 3-path space is rich enough to fill any 2-cycle

More precisely: the GLMY ∂₃ map from allowed 3-paths to allowed 2-paths
is surjective onto ker(∂₂) because tournament completeness ensures that
every 2-cycle can be expressed as a boundary of 3-chains.

## Significance

1. **Explains dimension gap:** The hereditary topology chain
   β₁ → β₃ → β₄ skips dimension 2 because β₂ never exists.

2. **β₃ appears "de novo":** S-phase (β₃=1) at n=6 cannot come from
   hereditary β₂ at n=5 (since β₂=0 always). The 3-dimensional topology
   is a genuinely new phenomenon at n=6.

3. **Structural constraint:** Any tournament-theoretic result depending
   on β₂ is vacuously true.

## Connection to Complete Digraphs

The complete directed graph K_n (both directions on every pair) has
trivial path homology (contractible). The tournament T is obtained by
keeping exactly one direction per pair. The β₂=0 result suggests that
even after this "halving," enough connectivity survives to kill dim-2 holes.

## Scripts

- `04-computation/beta2_vanishing.py`
- `04-computation/maximizer_hereditary_topology.py`

## Open Questions

1. Prove β₂(T) = 0 for all tournaments algebraically
2. Does β₂ = 0 generalize to "almost tournaments" (digraphs with ≥50% edge density)?
3. What is the minimum edge density required for β₂ to appear?
