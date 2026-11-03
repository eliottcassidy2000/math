# What is classified, what is constructed, and what remains open

Session144, 2026-11-03. All three lanes recover earlier work before extending it. Mathematical proofs, exact finite checks, numerical checks, external claims, and remaining questions are separated below.

## Jacobian: four different classification questions

1. **Generic degree values in ambient dimension n>=3:** inherited THM-3438 classifies the spectrum as `{1,3,4,5,...}` and constructs symmetric-monodromy primitive maps in every degree at least three. This is a spectrum of values and existence, not a list of maps up to coordinate equivalence. See the [inheritance and scope atlas](jacobian_followup/classification_scope_and_repairs.md).
2. **Fibers of one explicit map:** the [independent fixed-map audit](jacobian/independent/INDEPENDENT_AUDIT.md) completely classifies complex fibers using
   `L=27a²g²−18abg+16a+b³g−b²` and `Q=27ag²−9bg+8`.
   There are three points if `L!=0`, one if `L=0,Q!=0`, and none if `L=Q=0`. The zero set of Q outside L is a failure of the x projection to distinguish points, not ramification of the Keller map. Restoring the discarded wall equation repairs the old conic elimination.
3. **A universal equivariant iteration law:** the [one-positive-weight theorem](jacobian_followup/universal_equivariant_boundary_corollary.md) gives smooth, pairwise disjoint nonempty orbit-zero levels for every nonconstant scale in its stated category. After m iterates there are at least m irreducible components. The particular [degree-three tower](jacobian/equivariant_composition_tower.md) additionally has generic degree `3^m` and exact scale degree `(3*5^(m-1)-1)/2`. The universal theorem also describes all fibers of the positive-weight coordinate: its zero fiber is a disjoint affine hyperplane plus `G_m times V(l)`, while each nonzero fiber is the quotient complement. That component is neither a stable coordinate nor a locally nilpotent dual direction when l is nonconstant; no all-coordinate claim is made.
4. **Full map classes:** a complete classification under arbitrary polynomial source/target automorphisms, under equivariant equivalence, or after stabilization is not established. These are distinct relations. Necessary data, complete degree spectra, and a family with a proven operation law do not solve their realization and equivalence problems.

Two hostile examples sharpen the question. The [equivariant shear](jacobian_followup/equivariant_gauge_hostile.md)
`T_k=(x,y+x^(2k-1)z^k,z)` turns the primitive map G's affine scale into `2−3v−t−3t^k`, of arbitrary degree k, without changing its equivalence class. Thus even primitive nonlinear-scale **existence** is not an open problem. This supersedes the raw-degree formulation of the primitive question in the frozen first checkpoint. The meaningful question concerns primitive classes inequivalent to the known seeds. The composition tower's component-count obstruction is invariant under equivariant gauge, but those iterates are composite.

The [quartic boundary audit](jacobian_followup/independent/QUARTIC_BOUNDARY_AUDIT.md) supplies the other necessary correction: a known degree-four map has only one plane sheet plus two fold sheets over a generic point of its invariant target hyperplane. A fourth sheet escapes. Hence `fold degree = total degree − plane degree` needs a no-loss-at-infinity hypothesis. Smooth source boundary and full target sheet count are different predicates.

## ETH: inspect the theorem actually imported

The [ETH audit](eth/ETH_AUDIT.md) reads the pinned manuscript, formal target, encoding, solution endpoint and public replication reports. Eight source files were independently matched byte for byte; [the receipt](eth/INDEPENDENT_SOURCE_BYTES.json) records the exact hashes and upstream commit. An initially extra trailing newline in local materialization was corrected and is disclosed in the source ledger rather than silently erased.

The formal scope uses the weaker exponent `n(loglog n)^2/sqrt(log n)`, while the manuscript states sharper general and bounded-density bounds. The formal machine is chosen for each fixed width before the input formula and size. Unary literal fields have polynomial total encoding length. None of these source inspections is a new full Lean build or a reproof of the imported SAT algorithm; that large result remains CITED/EXTERNAL here.

The new [bounded-integer Keller collision interface](eth/ETH_AUDIT.md#7-a-bounded-integer-keller-collision-interface) is exact for its specified coefficient, degree and integer-witness bounds. Its determinant grid and gate-by-gate encoding preserve the intended predicate and count labelled ordered collision triples. They do not count polynomial-equivalence classes or cover unbounded complex maps. The solver's variable count is `N=O(W+G)`, including verification gates, not only W free witness bits. A merely polynomial reduction need not turn the cited bound into subexponential time in W. There is no automatic global classification or `P=NP` conclusion.

The exact small matrix experiment supplies a genuine information-loss test: deleting a rank-two correction changes a permanent from 1 to 81. Keeping one too few interpolation jets permits a nonzero invisible polynomial that changes the desired value by 72. These are hostile controls for unjustified recovery, not counterexamples to the imported paper.

## MRB: a complete filter criterion and a sharp nonuniqueness boundary

Let `f(n)=n^(1/n)−1`. The [MRB proof](mrb/RESEARCH.md) classifies convergence for every fixed periodic coefficient sequence: the sum `sum a_n f(n)` converges exactly when the period mean is zero. All odd-modulus filtered alternating masks converge; even moduli need equal even/odd weight totals. In particular an even unit mask destroys cancellation and has an explicit logarithm-squared divergence.

Cutoff and exponentially damped finite parts are specified operations. For mean mu, the Abel expansion has the linear-log term `−mu gamma log(1/t)`, and its finite part differs from the cutoff constant by `mu(gamma²/2+pi²/12)`. Dropping that term is not a legitimate regularization.

The stronger new obstruction is constructive: for every finite prefix, an entire smooth-Laplace perturbation can preserve that prefix, the entire branch jump, and all algebraic-order asymptotics uniformly on the right half-plane while changing the alternating constant by any chosen real amount. This is nonuniqueness in the enlarged analytic class, not a claim that the modified function still satisfies the defining equation `f(z)=exp(Log z/z)−1`. The [independent audit](mrb_independent/MRB_INDEPENDENT_AUDIT.md) and [precision audit](mrb/audit/README.md) verify this distinction.

The fixed MRB constant is polynomial-time computable in its requested output precision. The new evaluator uses explicit rational enclosures for log and exp and exact dyadic weights, not free transcendental oracles. Normal and optimized executions certify p=1,8,16,32; no p=64 certificate is claimed. This numerical-computability theorem does not recognize a modular identity or settle irrationality, transcendence, or SAT.

## The actual shared method

The three lanes illustrate one discipline, not a causal reduction between their subjects. Specify the object, quotient or encoding, preserved predicate, discarded information, reconstruction data, and cost. For Jacobian fibers the missing data include the wall equation and escaping sheets; for ETH encodings they include correction matrices, jets and gate counts; for MRB they include the global defining profile and the exact regularization operation. In each case a concrete hostile example identifies the first failed implication.

The useful stopping boundary is explicit: degree spectra, fixed-map fibers, equivariant boundary laws, periodic-filter convergence, and bounded encoding interfaces have exact scoped results. Full map-equivalence classification, a new certification of the complete ETH algorithm, and the arithmetic/modular nature of MRB remain outside what this session establishes.
