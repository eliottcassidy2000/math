# Synthesis: preserve the expensive coordinates

**Anchor:** explicit-NP unrestricted B2 lower bounds. **Niche:** exact Ramsey/tournament
transfer algorithms. **Wildcard:** affine symplectic geometry and inner product.
LRC(14) remains open; this user-directed session does not assert progress on its remaining case.

Five live concepts were compared: input-coordinate cost, colored affine flats, twin
projection, local amortized potentials, and certificate-preserving symmetry reduction.
The earlier odd-eighth/Cayley-Dickson work supplies motivation for retaining signs,
coordinates and multiplication data, not a circuit lower-bound theorem by analogy.

## Exact maps and their losses

- Input GL maps preserve full colored affine-flat geometry and Walsh multisets of
  point defects, but need not preserve B2 size. The missing information is the
  cost of realizing the linear input change. Translation, permutation and absorbed
  pin negations are safe here; arbitrary linear maps are not free wiring.
- Projection of twin-containing doubled chains preserves base transitivity. The
  lost primed interval is restored by its two endpoints, giving k(k+1)/2 lifts.
  Twin-free subsets need their absent/plain/prime state. Scalar maximum size alone
  cannot reconstruct the full polynomial.
- The provisional row table maps potential weights to coefficients. Its exact
  projection is proved, but coverage of actual circuit states, safe-core invariants,
  termination and uniform disperser construction are additional obligations.
- A rectangle in the fixed x/y input partition is not a clique. Inner-product
  rectangle/rank bounds lose the ability of unrestricted gates to mix both halves.

## Concurrency integrated, not overwritten

Commit 25626f02207500a2abb897049456d82bb63f1ece introduced the separate
[Ramsey saturated-port packet](../ramsey_saturated_ports_20261030/README.md).
Its stated faithful rooted restriction and joint attachment inventories concern
hypothetical (3,3,5;56) colorings. We read its README and structural proofs, preserve
its changes, and credit that work rather than present it as this session's result.
Its cited classification inputs and full certificates are not independently audited
here. It does not claim a new Ramsey bound.

The useful common rule is precise: quotient only with enough marked incidence data
for reconstruction. Rooted port compatibility and the doubled-chain interval are
concrete examples. Our affine counterexample shows that preserving a combinatorial
invariant still does not preserve circuit implementation cost.

## Next decisive experiment, not an asserted breakthrough

The complete potential frontier rules out improvement by retuning the same three
weights. A meaningful next target is a strengthened structural deletion lemma or a
finite-state/global potential proving that the worst cleanup and one-child gate
rows cannot occur at arbitrary joint frequency. A proposed certificate must retain:

1. The restricted domain, paired/unpaired coordinates and affine safe-core dimension.
2. The actual DAG, shared children, repeated pins and terminal-output exceptions.
3. A transition soundness proof and a reconstruction of every legal elimination case.
4. A fixed uniform P/NP family at every actual input length, with the resulting error term.

Before optimizing a transition graph, prove that its states are sufficient. Otherwise
an attractive cycle mean is merely another unsound quotient. No such strengthened
transition theorem has been established in this session. The completed results give
exact positive controls and obstructions against which to test it.
