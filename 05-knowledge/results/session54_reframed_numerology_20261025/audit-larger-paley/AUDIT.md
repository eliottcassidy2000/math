# Independent larger Paley certificate audit

Frozen producer manifest: 5a2aa000c7c1feff629b3f15587d7a46a6b354e0aa533ec7ef37a999f0bd0a60.

Verdict: PASS for P19, P23 and P31, subject only to the explicitly proved sufficient criterion reproduced with this packet.

## Independent checks and exact scope

The auditor imports no producer code. It constructs the tournaments and entire square-affine groups independently, checks primality of each order and 101 by trial division, checks group identity/composition/inverses, checks every affine permutation on every arc, and checks transitivity on arcs. Every supplied Hamilton cycle and forward 0-to-1 Hamilton path is validated. Every selected cycle belongs to the first seed's affine orbit. The connected bipartite support, indegrees, outdegrees, and both common-neighbor counts are checked directly.

Column elimination, separately implemented rather than imported or copied from the producer's row-elimination routine, gives determinant residues 59, 44 and 38 modulo 101 for minor sizes 134, 208 and 404. Nonzero residues imply nonzero integer determinants and rational lower bounds. Every difference is explicitly checked to have zero sums at all tails and heads. Connected support makes the potential space dimension 2n-1 and supplies the matching rank upper bound. Hence the rational ranks are exactly 134, 208 and 404. This is not a unimodularity or integral-lattice saturation conclusion.

The finite verification examines 135, 209 and 405 listed cycles respectively and the groups of 171, 253 and 465 elements. Its dominant determinant arithmetic is cubic in the listed rank, with quadratic matrix storage. It does not enumerate Hamilton cycles, Hamilton paths, or arbitrary deletion subsets at these orders. It also constructs only the three predicted Hall/star families to verify their cardinalities and distinctness, not to infer completeness or test blocking. The optional producer P7/P11 exhaustive modular probe is bundled as provenance but its complete-cycle census is outside this independent audit's claims.

## Why the theorem yields complete counts

The copied cycle-potential proof and prior independent audit establish the general criterion, not just the P11 instance. Reviewing the proof again: arc-orbit averaging applies to each Hamilton cycle separately, so a minimum blocker intersects every cycle in exactly two arcs. The rank equality therefore puts its incidence vector in the rational row-column potential space. Connectivity and integer arc values make the potentials integral by propagation along a support tree, without needing a unimodular cycle minor.

Normalized integer thresholds split the binary blocker into disjoint compatible cuts with strictly positive integer deficiencies totaling two. Zero deficiency at an active level contradicts support connectivity. Thus there is either one deficiency-two cut or two nested deficiency-one cuts. The common-neighbor Cauchy-Schwarz inequality has the stated denominator and both polynomial comparisons are correct. Their strict positivity, also checked at every relevant integer t here, excludes all interior deficiency-one and deficiency-two possibilities for r=9,11,15.

The one-level alternative leaves exactly the two-source and two-sink families. The two-level alternative leaves stars or mixed source/sink words; two rows or two columns would overlap. Arc transitivity transports each supplied forward endpoint path to exclude every non-star mixed word. Every remaining Hall deletion or star genuinely blocks a spanning path. Thus completeness follows analytically, not from examining the predicted family alone. Their sizes are 171+171+19=361, 253+253+23=529 and 465+465+31=961, with blocking numbers 18, 22 and 30.

## Discovery, portability, and hostiles

The producer's seeded randomized DFS witness discovery is heuristic. The frozen witness verification is fully deterministic and independent of discovery success rates or timing. A first successful seed orbit at each named order does not imply that every seed orbit works or that other Paley orders satisfy the criterion.

Normal and optimized independent outputs are identical. Replay verifies the exact closed file set and all SHA-256 digests, then runs from an arbitrary working directory. Its self-test rejects invalid cycles, reversed endpoint paths, wrong rank targets, duplicate minor columns, false first-seed metadata, and a composite modulus; it also checks singular and column-swap determinant controls. Addition, mutation, and removal of files are rejected, and a relocated copy replays without external workspace paths. All files are text. No producer or repository files were edited. There is no all-Paley-order extension, no P23 Hamilton-path count, and no literature-priority claim.
