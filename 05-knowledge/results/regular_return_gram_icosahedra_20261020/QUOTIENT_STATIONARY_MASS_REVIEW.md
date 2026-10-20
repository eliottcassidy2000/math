# Independent review of the stationary-mass corollary

**Verdict: PASS; READY as a separate addendum. No required fixes.**

## Structural checks

1. The inverse diagonal matrices are legitimate. All source row degrees are positive raw Hamiltonian-path counts, and the column universe contains only realized intermediates, so all column degrees are positive.
2. The two-stage walk is correctly normalized. Its conductance matrix is `B D_c^(-1) B^T`, symmetric without any regularity assumption. Detailed balance gives stationary measure H on every nonempty fixed labelled score fibre.
3. The regular fibre, unlike a general fixed score vector under arbitrary relabellings, is invariant under the full symmetric group. Both incidence and column degrees transform equivariantly. Orbit-row sums are therefore representative-independent, proving strong rather than merely stationary-distribution lumpability.
4. Summing the reversible measure over an ordinary orbit gives `(n!/|Aut(T)|)H(T)`. The n! factor is common at fixed order. This proves the exact class mass h, not H, without assuming an asymmetric representative.
5. Converse acts on paths by word reversal and commutes with the nonpath complement. It preserves regular sources, realized columns, H, c, and transition probabilities. A merged orbit contains one ordinary class if self-converse and two otherwise. The factor kappa is therefore correct.
6. The relation to fixed-path tiling fibres uses the established free path-orbit bijection, not a new formula for H. `THM-781-hamiltonian-path-inverse-metagraph-fibre.md` was read directly to verify the ordinary and converse-merged normalization conventions.
7. No uniqueness, irreducibility, global connectivity, or deterministic-complement lumpability claim is made. The distinction from the row-normalized path-pair matrix W and from the unrestricted two-step metagraph walk is explicit.

## Independent finite evidence

`check_quotient_stationary.py` independently implements the walk using directed arc sets and depth-first literal path generation. It imports the separate audit's arc-set utilities, not producer code.

- All labelled regular states at orders 1, 3, and 5 pass detailed balance: respectively 1, 4, and 576 matrix entries.
- The full relabelling orbits of the three inherited regular order-seven representatives are rebuilt as disjoint sets of sizes 1680, 720, and 240.
- Every incidence column in the three representative rows is generated from an actual path. Its degree is computed by enumerating all matching-endpoint continuations, each with its unique regular target.
- Exact quotient rows are computed, then independently checked after each of two generators of S_7 and after converse. This confirms covariance and the ordinary lumpability implementation.
- The stationary ratio 57:25:9 and the displayed quotient kernel pass exact rational detailed balance and stationarity.
- Normal and optimized Python runs produce byte-identical receipts. Explicit exceptions retain every check under optimization.

The completeness of the three regular order-seven class representatives is inherited from the audited `spectrum_basins_20261020e/regular.md` census; this checker does not claim a new all-order class census. The all-order proof does not depend on this finite example.

## Preservation

The V1 ready manifest, the four original frozen proof/code files, and all original package files remain byte-identical. This corollary and its checker are additions only. `READY_MANIFEST_V2.json` lists the combined package; `PUBLICATION_SHA256SUMS_V2` supplies package-relative verification while the V1 hash files remain preserved.
