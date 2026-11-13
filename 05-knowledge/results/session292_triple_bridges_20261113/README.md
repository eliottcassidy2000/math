# Session 292: guarded fibres, sparse certificates and constructive boundaries

2026-11-13. **PROVED / FINITE-EXACT / CITED / OPEN**, distinguished below.
Chan–Pak completeness remains OPEN. No closure of the 473 repair obligations,
RH proof, or new A000568 enumeration is asserted.

## Main advance: classify what the boundary cannot distinguish

With named interior vertices, fixed port incidences and port–port arcs, the
complete squarefree path-response matrix is equal for two tournaments **if
and only if their augmented labelled cyclic-triple inventories agree**.
Only cubic interior observations and quadratic port observations are needed
for this relative classification. This is not reconstruction from an unmarked
Hamilton count, nor a polynomial-time algorithm for computing it.

The [proof](triples/PROOF.md) leads to an exact
[colored modular-fibre grammar](triples/COLORED_FIBERS.md). Color vertices by
their port-incidence vectors. The number of distinct labelled tournaments
with that complete boundary is

    2^(monochromatic prime nodes) × product (pure same-color linear run size)!.

One deliberately constructed port makes this fibre a singleton. At most n−1
selected quadratic observations suffice after the interior triple inventory
and port incidences are given. Random k ports instead leave a nontrivial fibre
with probability at most (n−1)2^(−k); transitive tournaments give the sharp
logarithmic scale. These are physical tournament extensions, not count-
preserving operations. They expose how much endpoint information a proof
must deliberately retain or allow to change.

Keeping only selected quadratic guards gives a more flexible exact grammar:
prime quotient reversals survive unless locked by a retained comparison;
linear quotient orders become the linear extensions of a comparison poset.
Thus prospective compensating moves can be generated as actual tournaments
instead of formally subtracting path counts. Exact target gains and renewable
endpoint readiness remain the constructive obligations.

## Other proved bridges

- [Sierpiński and endpoint algebra](triples/PROOF.md): destination ideals remove
  Schur denominators; two-pivot elimination has five positive terms. Divided
  powers give a sharp factorial-valuation truncation modulo prime powers,
  attained by tournament cycle denominators. Full-observer seam proposals
  fail by rigidity, not merely through unsuccessful finite searches.
- [First-25-prime wheel](arithmetic/PROOF.md): the inherited ray 99·2^k+23 has
  least sieve period **61,209,720**, with **32,353,278** surviving residues.
  A conditional CRT product explains the shared-clock dependence. Every
  surviving cylinder contains subcylinders forcing arbitrarily many distinct
  additional prime factors. Survival is not primality or target landing.
- [Square-sum reflection words](square_sum/PROOF.md): exact prefix domains,
  collisions and first-birth generating functions. Consecutive square colors
  t²,(t+1)²,(t+2)², t≥4, first realize the repeated-three-color six-cycle at
  floor((t²+6t+5)/2)+1 vertices. The first instance is at 23. This concerns
  the specified word family, not all cycles or Hamilton paths.
- [Analytic inspiration and finite observation](manuscript/PROOF.md): any
  eight distinct lengths distinguish the fixed eight-mode observer; seven
  need not. Escaping-zero examples and an inherited endpoint hostile identify
  what finite data cannot establish. The linked manuscript supplies ideas,
  not an assumed RH theorem. The stated totient inequality is traced to
  Nicolas's published criterion in [sources](manuscript/SOURCES.md).

## Research loops, provenance and checks

[Connection web and next proof obligations](SYNTHESIS.md) records the exact
maps and their limits, rather than identifying structures because both have
three or five elements. [History](history/CONNECTIONS.md) separates inherited
mechanisms from this session's advances.

Independent reports: [rigidity](audit/rigidity_reviewer.md),
[colored fibres](audit/colored_fibers_review.md),
[geometry/arithmetic](audit/manuscript_reviewer.md),
[endpoint/triple proofs](audit/square_reviewer.md), and
[manuscript lane](audit/triples_reviewer.md). Scripts and JSON receipts live
beside the claims. Controls include 33,864 tournaments through order six for
the relative classifier, 248,832 three-color cases at order five, all 1,099
engineered-port tests through order five, the complete 61-million-residue
wheel, and all-word/birth checks. Normal and optimized Python agree.

The geometric reviewer repaired a missing minimum-cycle-length condition
before publication; see [boundary record](square_sum/BOUNDARY_REPAIR.md).
The final manifest identifies the exact delivered files. Verification is a
check on the proofs and declared finite scope, not evidence of global coverage.
