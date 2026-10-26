# Fixed twists, conductors over non-UFD bases, the reverse valuation potential, Antihydra, and LCA merges

2026-10-25, klein (worktree `codex/session-tournament-twist-unify`, second half). Exploration session
on the owner's directive; statuses are per lane. **PROVED** statements have complete proofs in the
lane notes; **FINITE-EXACT** receipts were produced normally and under `-O` and agree byte for byte;
**CITED** and **DERIVED-UNDER-MODEL** items are labelled. No independent agent audit was performed
(the owner asked for exploration rather than audit). Canon entries: THM-4772 (conductors) and
THM-4773 (fixed-twist lifts); the Antihydra lane takes no ID (opus-S56 reservation).

| Lane | Result | Status |
|---|---|---|
| [Conductors over normal bases](conductors/DEDEKIND_CONDUCTORS.md) | Normalization of twisted multi-root algebras is `⊕_z c_z^{-1} u_z` with `c_z` the square-root ideal of `b^z`; mixed reduced roots appear exactly at primes with two odd valuations; per-character (port) ideal classes; conductor valuation matrix `v_p(d_z) = max(0, max_{z'}(v_p(c_{z'}) - v_p(b^{z∧z'})))`; the full normalized-character compatibility bank is retained; a shared-prime shortcut was corrected | PROVED; FINITE-EXACT on `(15,21,55)` over `Z[1/2]` and `Z[√-21][1/2]`; prime-2 hostile |
| [Fixed-twist boundary lifts](twists/FIXED_TWIST_LIFTS.md) | Pointwise squareclass criterion; open central line: three pair-port conditions plus a conic; complete central divisor: zero, eight or infinitely many points (corrected Oct26); Pythagorean twist gives `(0:3:3:4:4:5:5)`; outer lines are genus-five (finite) | PROVED (central over infinite fields; over Q via Hilbert symbols), CITED (genus five, Faltings); 300 random twists agree with brute force |
| [Reverse valuation potential](counting/REVERSE_IMPLEMENTATION.md) | Implemented on the production ordinary A000568 graph (13,445 states at `N=100`); `B = floor(N/2)`; normalized digits exact; fixed-precision pruning: `K=16` keeps 51%/25% of states/edges at `N=100`, `K=1` keeps 550 states | IMPLEMENTED + FINITE-EXACT; pruning PROVED |
| [Antihydra](antihydra/ANTIHYDRA.md) | `floor(3a/2) = -T(-a)` with `T` the THM-2228 Mahler carry map; odd first-return map `3^{v-1} S^-(a)`; valuation sequences form a full shift with cylinders = 2-adic balls; halting iff `sum (2v_i - 3)` reaches `-6` (rational slope `3/2` vs Collatz's `log_2 3`); model halting probability `phi^{-6}`; 2×10^6 exact steps | PROVED (A1–A4, A6), DERIVED-UNDER-MODEL (A5) |
| [LCA and Petersen](lca/LCA_AND_PETERSEN_MERGES.md) | Valuation partial sums are trie LCA depths; the Gale tournament of a linked `K6` is the clone hostile `C3[TT2,TT2,TT2]` with depth triple `(1,1,3)`; ten 2-subsets of a 5-set need five hierarchies (`chi(L(K5)) = 5`) | PROVED (small), inheritance surveyed |

Not pursued: the six order-12/14 ledger failures of
[LOW_ORDERS](../port_residual_integral_20261025/LOW_ORDERS.md). They require replacement word
templates inside the marked-path/K-move framework, which this session did not re-implement; they
remain template failures, not connectivity obstructions, exactly as that packet states.

## Inheritance and board

Closest proved mechanisms: BOUNDARY B1/CONTINUATION D2 (conductor and one-step test),
TORIFICATION T1 / THM-4765 (root surfaces and twist fibres), REVERSE_POTENTIAL (the construction),
THM-2228 (Mahler carry conjugacy), Session56 (valuation-LCA). Canonical hostiles: the repeated
radicand (`b_1 = b_2`), the prime 2 (`(u+v)/2` integral over `Z`), the naive initial subtraction of
`B`, and the i.i.d. model's silence about the actual orbit. Corrected near misses: twisting breaks
normality (so the test bank must target `A'`), and `B ≠ v_2(a(N))` at even `N`. Least-used sidecars:
per-character fractional ideals; the reduced-cost zero core of the counting graph; the trie reading of
valuation sequences. See [CONCEPT_BOARD.md](CONCEPT_BOARD.md).

## Reproduction

```
cd conductors && python3 check_dedekind_roots.py; python3 check_conductor_matrix.py     (SymPy)
cd twists && python3 check_central_conic.py; python3 check_boundary_lifts.py             (stdlib)
cd counting && python3 reverse_potential.py; python3 reverse_prune.py                    (stdlib; ~15 s)
cd antihydra && python3 antihydra_walk.py 300000; python3 antihydra_sim.py 2000000        (stdlib; the 2e6 run takes ~4 min)
cd lca && python3 check_lca_merges.py                                                     (stdlib)
```

Each produced `*-normal` and `*-optimized` receipts that must agree byte for byte; `SHA256SUMS` pins
the files.
