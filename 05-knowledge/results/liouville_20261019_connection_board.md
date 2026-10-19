# Connection board: from the square-sum seed to Pythagorean trees, tournaments, and elliptic curves

**Updated 2026-10-19.** This is a research map of the session's proved statements, cited inputs, exact experiments, and remaining questions. No major named conjecture or literature-priority claim is asserted. All linked packets give their precise quantifiers and reproducible controls.

## 1. The organizing example

The labels `(482,3362,6242)` have common difference `2880` and pair sums `62²,82²,98²`. Their half-roots satisfy `31²+49²=2·41²`, yielding the primitive Pythagorean triple `(9,40,41)`. In the repository's Euclidean-ratio convention this is the third left descendant of `(3,4,5)` in the Berggren ternary tree.

This is a complete classification route, rather than an isolated numerical analogy. Every positive AP square-sum triangle is uniquely

`2g²(C²−4AB, C², C²+4AB)`

for a primitive Pythagorean triple with `A>B` and `C²>4AB`. The classification proves that `6242` is the smallest possible largest label. See [AP/PPT classification](ap_square_pythagorean_20261019.md).

## 2. Three structures on the same arithmetic objects

### Ancestry and positivity

Berggren's three branches generate the primitive triples; positivity selects two narrow angular cones. Inside the admitted subgraph only an outer child continues, giving disjoint unary rays. Nevertheless membership cannot be recognized by a finite automaton on the full branch word: the irrational cone walls have degree four, producing infinitely many distinguishable residual cuts. See [nonregular branch language](berggren_ap_square_nonregular_20261019.md).

### Prime factorization and tournament paths

For a hypotenuse with `r` distinct prime divisors, all congruent to 1 modulo 4, primitive triples are parametrized by choosing one Gaussian prime above each rational prime, modulo simultaneous conjugation. These choices are exactly the complement-paired nonexceptional Hamilton-path addresses of the endpoint tournament `E_(r+2)`.

Thus the `E_3,E_4,E_5` fibers have sizes `1,2,4`. Retaining the prime exponents and angular phases makes this a bijection and a weighted identity; retaining only the scalar path count loses the positivity information. See [tournament/PPT fibers](endpoint_pythagorean_fibers_20261019.md).

### Elliptic continuation and graph gluing

The gap squareclasses are exactly congruent-number classes. The associated elliptic curve supplies infinitely many rational continuations and, after a common square scaling of any finite family, integer graph constructions. This is a third coordinate system, not a claim that tree ancestry, factorization, and elliptic addition are identical.

See [congruent gaps](congruent_square_gap_classes_20261019.md), [PPT port bridges](square_sum_ppt_bridges_20261019.md), and [common neighbors](ap_triangle_common_neighbors_20261019.md).

## 3. What this yields beyond the seed

1. Every admissible primitive petal has infinitely many compatible opposite-port neighbors with spreads arbitrarily close on either side. This constructs arbitrary finite AP-paired fans and selected-label square-sum cycles. Rank compression gives a separate graceful Walecki labeling; it does not preserve square sums.
2. Every AP square-sum triangle, after a scale chosen for the finite target, extends to `K_3 join I_k` for every finite `k`. In particular it supplies `K_4` and `K_{3,3}` square-sum subgraphs.
3. Three rational AP anchors have infinitely many positive rational common neighbors exactly when their gap has congruent-number squareclass. Four AP anchors have none, by Fermat's four-squares theorem.
4. The original seed has **no distinct integer common neighbor**. Its sole integer solution of the three square equations is the middle anchor itself. Rescaling is essential.
5. No uniformly bounded number of translates in the specified elliptic contact search succeeds for all seeds. At the singular limit the translation becomes the Berggren middle map `u -> 2+1/u`; its Pell iterates stay outside the positivity cones. This proves a concrete connection between the tree and elliptic obstruction, while identifying precisely where the connection holds. See [contact-index theorem](square_sum_contact_step_bound_20261019.md) and [algebra supplement](square_sum_contact_step_algebra_20261019.md).

The bounded search succeeds for all 528 admitted primitive triples of hypotenuse at most 10000 by absolute translate index 11 when both positive-spread charts are accepted. That finite fact is not a universal bound.

6. The complementary quantitative result bounds a successful signed contact index by `O(1+log(1/f))`, hence `O(log C)` in the primitive input hypotenuse. It follows from a real-circle pigeonhole argument and a uniformly long one-sided admissible arc. A fully rational stopping budget avoids uncertified numerical integration. Directional/window-constrained bounds are also explicit. See [logarithmic recurrence bound](square_sum_contact_recurrence_bound_20261019.md). This controls the number of elliptic additions, not their bit cost or the output hypotenuse.

7. For the prescribed plus-port search, the preceding upper order is sharp near either zero-spread cusp: the first successful signed index is `Theta(log(1/f))`, including odd indices and both accepted charts. The proof uses two bounded admissibility arcs separated by half a diverging real period. Along either outer Berggren spine of depth `j`, this becomes `Theta(log j)`, while the input hypotenuse is `Theta(j²)`. See [matching cusp lower bound](square_sum_contact_cusp_lower_bound_20261019.md).

## 4. Odd zeta and graph geometry: exact maps

The signed primitive-triangle Dirichlet series is

`P(s) = 1/2 [zeta(s) beta(s) / ((1+2^(-s)) zeta(2s)) − 1]`, for `Re(s)>1`.

At odd positive arguments 5,7,9 this is a rational multiple of `zeta(s)/pi^s` minus one half. Positive AP triangles require angular Hecke terms; simply taking one third of this weighted series is false. See [PPT odd-zeta series](ap_square_odd_zeta_20261019.md).

This norm-height series is distinct from the Eisenstein scattering coefficient at indices 3,4,5, which involves zeta at 5,7,9. The [E atlas](e_atlas_20261019.md) keeps those interpretations, root systems, and tournament notation explicit.

The Petersen boundary intersection lattice has a separate integral incidence description via `K_5`, with `A_4=E_4` perpendicular to the anticanonical class. Stable-tree boundary products and tensor signatures are given in [odd-zeta boundary graphs](odd_zeta_boundary_graphs_20261019.md). These geometric identities do not prove zeta irrationality.

## 5. Quadratic dynamics and graceful exchanges

The rational period-three orbit at `c=−29/16` maps affinely to the seed. Real order-three graph realizations have exact orbit-capacity constraints; rational realizations need additional conditions. At `c=−7/4` the parabolic three-cycle is not rational. See [period-three bridge](period3_graph_square_bridge_20261019.md).

The session also formulated two precise open exchange conjectures:

- **GP3:** all gracefully labeled paths on fixed labels `0,...,n−1`, modulo reversal, are connected by moves replacing at most three edges, with endpoints allowed to change. Exhaustive connectivity holds through `n=14`.
- **GT1:** all gracefully labeled trees on fixed labels are connected by replacing one edge with another of the same difference, while allowing tree shape to change. Exhaustive connectivity holds through `n=11`.

See [graceful exchange packet](graceful_square_exchange_20261019.md). Neither conjecture by itself establishes that every prescribed tree shape admits a graceful labeling. A tempting proof method using only a fixed short-color buffer is disproved by explicit arbitrarily large disconnected frozen-tail fibers: [fiber obstruction](graceful_fiber_obstructions_20261019.md).

## 6. Next proof questions, ranked by the structure now available

### A. Quantitative Pythagorean contact construction

The existence problem and adaptive index bounds are settled for these ports; for the prescribed plus-port search the first index has matching logarithmic order at the cusps. The next useful theorem bounds the output hypotenuse or bit complexity, including a prescribed spread window. Another refinement would determine whether the normalized first index has a limit or several accumulation values; matching order alone does not give an asymptotic equivalent. These questions require more than the number of elliptic additions.

### B. Denominator-controlled common-neighbor gluing

For a fixed AP core and target `k`, bound the necessary integer square dilation. Any fixed scale has only finitely many common neighbors, already by factoring the difference of two squares for two anchors. The infinite rational construction therefore needs a quantitative denominator theorem before it can help with dense intervals or economical integer realizations.

### C. Graceful exchange with a movable obstruction

Find a monotone invariant or repair move that escapes the frozen-tail fibers while retaining a prescribed global target. Three-edge path moves have strong finite evidence, but a proof must control all label orders and endpoint changes. Shape-preserving consequences require their own argument.

### D. Transfer to the major targets

A useful transfer must specify a map and what it preserves. Current square-sum `K_4` constructions cannot be transferred wholesale to planar unit-distance graphs. The selected-label cycles do not cover every consecutive integer. Angular norm series do not supply irrationality estimates. These are concrete missing steps for Hadwiger–Nelson, interval square-sum problems, and odd-zeta arithmetic, respectively.

No Collatz termination, Wilson-prime completeness, or planar Jacobian conclusion follows from these packets. The established bridge is the shared arithmetic object and its explicitly related coordinates; extending it to those targets remains a separate research problem.

## Publication and inheritance

All packets from this session are committed directly to `main`. Earlier work on the separate historical research branch has not been silently merged or deleted. The inherited tournament, sign-gauge, marked insertion, and square-sum results were used as starting points; this board replaces the initial question-only checkpoint with the current result map.
