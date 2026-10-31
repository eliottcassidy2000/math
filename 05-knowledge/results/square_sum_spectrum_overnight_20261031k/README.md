# A finite source-specific coverage kernel on the105 boundary atlases

**Status: FINITE-EXACT full-matrix certificate, with PROVED all-time consequences for this finite graph only.** This is a successor to [packet j](../square_sum_spectrum_overnight_20261031j/README.md). Universal tournament coverage, cyclic arithmetic seam realization and interval absorption remain OPEN.

The vertices are all105 perfect matchings on eight **marked incidence slots**. Two vertices are adjacent when their union is one alternating eight-cycle. Packet j identifies each neighbourhood with the48 possible single-cycle boundary types relative to its marked port matching. The older matching/Hamming dictionary turns these into marked atlas objects, with its basis/origin retained. It does not produce a unit-distance gadget.

Let A be the adjacency matrix. Exact integer arithmetic gives:

    every row sum is48;
    every entry of A^2 is20,24 or48;
    spectrum(A) = 48^1, 4^56, (-2)^14, (-6)^14, (-8)^20.

Each row of A^2 has60 entries20,44 entries24 and one entry48. Thus **every ordered source-target pair has at least20 two-step routes**. There is no zero-mass exceptional state.

For K=A/48 and uniform stationary law u, this implies

    K^2(s,t) >= 20/48^2 = 5/576 = (175/192) u(t)

for every actual source s and target t. Consequently K^2=(175/192)U+(17/192)R for a stochastic matrix R, where every row of U is u. Iterating gives a total-variation contraction by at most(17/192)^m after2m steps. Symmetry and the full spectrum give the sharper Euclidean estimate

    ||mu K^t-u||_2 <= 6^(-t) ||mu-u||_2

for every initial probability row mu and integer t>=0. These are all-time finite-state consequences, not extrapolation to growing state spaces. Source-specific positivity follows directly from the two-step entries, independently of an average or stationary limit.

## Exact audit

`atlas_kernel.py` enumerates all105 matchings and builds the full105-by105 integer matrix. It checks every two-step entry and every row histogram. It then verifies

    (A-48I)(A-4I)(A+2I)(A+6I)(A+8I)=0

as an equality of full integer matrices. Since A is real symmetric, its eigenvalues lie among these distinct roots. Exact traces of powers0..4 evaluate the Lagrange spectral-projector traces and give the displayed multiplicities, including a one-dimensional stationary eigenspace. This avoids inferring the global spectrum merely from an equitable quotient or floating-point eigenvalues. A complete five-orbit quotient discovered the candidate roots; the full-matrix audit is the proof certificate.

Normal and optimized Python outputs are identical. Reproduce from the repo root:

```bash
python3 05-knowledge/results/square_sum_spectrum_overnight_20261031k/atlas_kernel.py
python3 -O 05-knowledge/results/square_sum_spectrum_overnight_20261031k/atlas_kernel.py
python3 05-knowledge/results/square_sum_spectrum_overnight_20261031k/verify_packet.py
```

## Research consequence and stopping boundary

This gives an exact finite model of the desired research pattern: retain each source, make every target reachable with a pointwise lower bound, and only then use an averaged norm estimate. The new bound is usable for computations whose legal operations really are these matching transitions. It cannot pay Collatz reset debt, force geometric palettes or synthesize tournament counts until an authenticated map preserves those operations and predicates. A tree of native square-sum carriers has such a paid arithmetic construction; an arbitrary circuit of these abstract atlas transitions does not yet have one.

**Anchor:** interval ownership with exact boundary states. **Niche:** arithmetic realization of a small cyclic matching controller. **Wildcard:** test whether a tournament source transition admits this finite boundary kernel, with its source labels and clocks preserved. The cheap hostile is a proposed transition that changes the actual port pairing while preserving only its cardinality.
