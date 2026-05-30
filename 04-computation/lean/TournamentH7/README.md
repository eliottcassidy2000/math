# TournamentH7 — Lean 4 formalization of THM-343 (H(T) ≠ 7)

A Lean 4 + Mathlib formalization of the theorem

> **THM-343.** For every tournament `T` on any number of vertices `n`, the
> Hamiltonian path count `H(T)` is never equal to 7.

This is the project's flagship result that the value 7 is permanently
forbidden from the H-spectrum. The original informal proof is in
`/home/ubuntu/math/01-canon/theorems/THM-343-H7-impossible.md`.

## Structure

```
TournamentH7/
├── lean-toolchain            leanprover/lean4:v4.30.0
├── lakefile.toml             mathlib v4.30.0
├── TournamentH7.lean         root: imports everything
└── TournamentH7/
    ├── Basic.lean            Tournament structure
    ├── Cycles.lean           DirectedCycle T k
    ├── SCC.lean              Reachability, IsSCC, IsHamiltonianPath, H
    ├── OCF.lean              7 axioms (OCF + Moon-Moser + Moon-Camion + …)
    ├── H7.lean               The proof:  Tournament.H_ne_seven
    ├── GoodCuts.lean         good-cut bucket constraints and exact spectrum
    ├── StaircaseConnectivity.lean
    │                          concrete tilings -> tournaments; top bucket iff SC
    ├── BucketBalance.lean    abstract finite bucket half-line and unordered balance
    └── Verify.lean           #print axioms audit
```

## How to build

```bash
cd 04-computation/lean/TournamentH7
lake exe cache get      # one-time: download Mathlib oleans (~6 GB)
lake build TournamentH7
```

`lake build` will print the `#print axioms` output from `Verify.lean`,
showing exactly which axioms `H_ne_seven` depends on.

## Axiom audit

`#print axioms Tournament.H_ne_seven` yields:

```
propext, Classical.choice, Quot.sound,              -- Lean foundational (always present)
ocf,                                                -- Grinberg-Stanley 2024
moonMoser, moonCamion_oddSize,                      -- Moon-Moser 1962, Camion 1959
oddCyclesIn_size3, oddCyclesIn_size4,               -- small-SCC enumeration
omegaTriangleLocalises,                             -- folklore: cycles ↪ SCC partition
alpha_subset_bound                                  -- independence polynomial subset bound
```

No `sorry` is used anywhere. The core H≠7 proof uses the cited axioms in
`OCF.lean`; later project modules introduce additional documented axioms
for H≠21, finite H≠63, tiling score formulas, and counting constants. See
`ARCHITECTURE.md` and `SUBMISSION.md` for the current full audit.

`GoodCuts.lean` also contains the formal axiom-free core of THM-336: the
all-down tiling is exactly bucket 0, any upward tile forces at least two good
cuts, no tiling has exactly one good cut, the strengthened dichotomy
`goodCutCount = 0 ∨ 2 ≤ goodCutCount`, grid reflection preserves the
bucket index, and for `n >= 3` the exact bucket spectrum is
`{0} ∪ {2,...,n-1}`. `StaircaseConnectivity.lean` constructs the
tournament induced by a concrete tiling and proves the top bucket
`goodCutCount = n - 1` is exactly strong connectivity of that induced
base-path tournament. `BucketBalance.lean` isolates the finite
internal/escaping half-line count used by quotient-bucket arguments:
`|selfHalf| + |crossHalf| = |fiber| * |moves|`. It also formalizes the
partner-map layer for unordered balance: internal half-lines are closed under
involutive moves, fixed-point-free moves give no self-partners, finite
fixed-point-free involutions have even cardinality, and the unordered balance
follows for fixed-point-free involutive move systems. These are the Lean cores
of THM-348 and THM-350; the remaining bridge toward full THM-346 is the
Boolean-mask specialization.

## Proof sketch

By `ocf`: H(T) = 1 + 2α₁ + 4α₂ + 8α₃ + 16α₄. Setting H = 7 gives
α₁ + 2α₂ + 4α₃ + 8α₄ = 3. Non-negative solutions are (3,0,0,0) and
(1,1,0,0). The case (1,1,0,0) requires a vertex-disjoint pair of odd
cycles with only one odd cycle total, contradicting `alpha_subset_bound`.

So α₁ = 3 and α₂ = 0: three odd cycles, all pairwise vertex-meeting.
By `omegaTriangleLocalises` they lie in a single SCC `S` of size `s ≥ 3`,
with `oddCyclesIn T S = 3`.

Case-split on `s`:

| s     | bound                         | source                  | contradicts oddCyclesIn = 3 |
|-------|-------------------------------|-------------------------|-----------------------------|
| 3     | `oddCyclesIn T S = 1`         | `oddCyclesIn_size3`     | 1 ≠ 3                       |
| 4     | `oddCyclesIn T S = 2`         | `oddCyclesIn_size4`     | 2 ≠ 3                       |
| 5     | `oddCyclesIn T S ≥ s − 1 = 4` | `moonCamion_oddSize`    | 4 > 3                       |
| ≥ 6   | `oddCyclesIn T S ≥ s − 2 ≥ 4` | `moonMoser`             | 4 > 3                       |
| ≥7 odd| stronger but unneeded         | `moonCamion_oddSize`    |                             |

Q.E.D.

## What this proves vs. what it assumes

The Lean proof verifies that **conditional on** the seven mathematical
axioms in `OCF.lean`, H(T) = 7 is impossible. Each axiom is either:

* a published external result (OCF, Moon-Moser, Moon-Camion), or
* a folklore/computational fact (size-3/4 SC enumeration, SCC partition
  of cycles, independence-poly subset bound).

The session-S5 informal proof (`THM-343-H7-impossible.md`) treats all of
the above as known; the Lean file makes this dependency explicit and
machine-checks the case-split logic that combines them.

## Future work (de-axiomatisation roadmap)

1. **Discharge `oddCyclesIn_size3` / `oddCyclesIn_size4` by `decide`**:
   for small SC tournaments these are decidable propositions; with a
   computable definition of `oddCyclesIn` these axioms can become lemmas
   provable by `native_decide`.
2. **Prove `alpha_subset_bound` directly** from a computable definition
   of `alphaCount` as `Finset.card` of independent sets.
3. **Prove `omegaTriangleLocalises`** from a Lean development of the SCC
   partition of digraphs. Mathlib has `SimpleGraph.ConnectedComponent`
   but no direct SCC theory for digraphs yet.
4. **The two hard axioms** — `ocf` (deep combinatorial identity) and
   `moonMoser` / `moonCamion_oddSize` (classical structural theorems) —
   require substantial Mathlib development.
