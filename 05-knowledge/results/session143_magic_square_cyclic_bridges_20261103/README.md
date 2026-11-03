# Session143: magic-square and cyclic-structure bridges

**Current checkpoint: 3. The magic square of nine distinct integer squares remains OPEN.**
Session began **2026-11-02 at 23:55 UTC** (user-local November2); source checks
and release continued **2026-11-03 UTC**. The suffix records the release date.

## Main results

1. **Complete real cubic classification.** For `k=-4,...,4`, the values of
   `Ak^3+Bk^2+Ck+D` can be rearranged into a 3-by-3 magic square exactly when
   `B=0` and either `A=0` or `C=-15A`. This excludes every genuinely quadratic
   assignment and every squared nonconstant root arithmetic progression.
   The surviving nonlinear family is an ordinary magic array, not a solution
   with nine integer-square entries.
2. **Exact cyclic/algebraic maps.** The given nine-cycle needs ternary carry,
   not an affine plane or cube relabeling. Signed octonion coordinates give
   a chosen subalgebra lattice, a common multiquadratic character action,
   and an explicit order-four lift of the grid half-turn preserving all
   eight magic-line equations. None proves existence. The user's intended
   **lambda cube remains unresolved**; type-theory comparisons are conditional.
3. **Repaired n=5 spectral bridge.** The census is 12 classes, 8 self-converse,
   10 modulo converse, but **nine adjacency spectra and one regular class**.
   `5:4` acts on four orientations, not as one tournament's automorphism
   group. The actual operator `rho(g)+rho(g)^2` has Fourier multiplicities;
   cited Baby Monster data specialize them without identifying characters
   with eigenvalues or constructing a missing representation/intertwiner.

4. **Explicit puzzle coordinates.** Two legal closed Conway-puzzle words on
   `PG(2,3)` generate a chosen `5:4` subgroup. Its five-point orbit gives the
   regular tournament and its four orientations. This is a genuine finite
   construction, not a projective collineation or Baby Monster intertwiner.

5. **Sharp integral refinement.** Ordinary polynomial calculus has
   `[Z[P]:Z[P+P^2]]=11`, with generation modulo `m` exactly when11 is a
   unit. This is **not** loss from the full adjacency matrix: retaining
   transpose/all-ones data or a Hadamard operation recovers `P` integrally.

6. **Abstract Baby Monster corollary.** The explicit puzzle `5:4` embeds
   abstractly into the cited normalizer, selecting class5A. The4371-module
   restriction has870 four-dimensional and891 total linear constituents;
   whole five-dimensional permutation-module multiplicities do not follow.
   No commuting outer involution or explicit Baby Monster intertwiner is assumed.

## Read

- [Final arithmetic and abstract-group addenda](CHECKPOINT3.md) and [independent audit](independent_integral_order/REVIEW.md)
- [Checkpoint 2 overview](CHECKPOINT2.md) and [synthesis](SYNTHESIS_CHECKPOINT2.md)
- [Magic/cycle theorem](magic_cycle/THEOREM.md)
- [Cube/octonion proofs](lambda_octonions/PROOF.md)
- [Tournament census and four-orientation action](tournaments/RESEARCH.md)
- [C5 operator and cited Baby Monster specialization](sporadic/report.md)
- [Puzzle-coordinate construction](puzzle_coordinates/report.md) and [independent audit](independent_puzzle/REVIEW.md)
- [Independent full audit](independent/REVIEW.md)
- [Verification and limits](VERIFICATION.md)
- [Historical checkpoint 1 overview](CHECKPOINT1.md)

From any working directory, run `python3 /path/to/packet/replay_checkpoint3.py`.
It uses the Python standard library, fresh relocation, exact output comparison,
normal and optimized runs, content pins and hostile corruption controls.
`replay_checkpoint1.py` and `replay_checkpoint2.py` remain valid on their
immutable scientific source sets.
GAP/CTblLib reconstruction was **not run**. No novelty priority, unspecified
Hecke map, general sporadic correspondence or unrelated conjecture closure is
claimed. Historical source snapshots remain provenance, not blanket endorsement.
