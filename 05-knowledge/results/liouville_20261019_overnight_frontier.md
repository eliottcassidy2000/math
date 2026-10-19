# Overnight frontier: exact context transport, signed carriers, and cubic contact chains

Session: 2026-10-19 UTC. This continues the [earlier connection board](liouville_20261019_connection_board.md). The earlier sparse AP-triangle, fixed-contact, and graceful-exchange results remain in force. The present work concentrates on new proof mechanisms and the obstacles that survive them.

**No general Collatz, Gold Partition, 1/3–2/3, graceful-tree, or consecutive square-sum result is claimed.** All-size statements below have specified carriers and explicit proofs. Historical novelty is not established merely by finding a proof in this session.

## 1. Gold partition: exact weights now support genuine structural families

The [unbounded endpoint-flag theorem](gold_unbounded_flag_context_20261019.md) gives a fixed two-comparison Gold strategy and balance at least 2/5 for a nine-element core inside arbitrary terminal-upset/initial-ideal collars. Its proof retains the exact marginal weight

    (1+r·1{left port first})(1+s·1{right port last}), r,s≥0.

The strategy works on the entire nonnegative cone of the four flag measures, not only a finite sample of contexts. The [bounded four-collar strategy](gold_four_collar_20261019.md) has the stronger uniform normalized Gold slack 109/1174 when r,s≤1. The [terminal-chain classifier](gold_terminal_context_20261019.md) shows why a limiting equality is insufficient: an arbitrarily high-order near-uniform perturbation can leave a constant signed deficit forever. It also supplies an explicit repaired boundary selector.

The [symbolic two-chain family](gold_two_chain_fork_20261019.md) proves the exact parameter region m≤n and n²−mn−(m−1)(m−2)≤0 for a fixed strategy under all flag weights. It permits an arbitrary side poset attached through the first chain element. An eleven-element member has uniform normalized Gold slack 1/56 under arbitrarily large allowed contexts. The curved equality boundary has a complete Pell/Fibonacci recurrence; after an explicit affine change of coordinates its evolution is the fourth power of the Gold-budget matrix. The diagonal is a separate zero-margin boundary. This class is not all posets.

## 2. Poset balance: mixing is automatic; the inequality is the remaining difficulty

The [ideal transfer theorem](ideal_transfer_positivity_20261019.md) proves that maximum incomparability degree D forces every D-step ideal-layer matrix to be strictly positive. Positive diagonal scaling places every entry between 1 and binom(D,floor(D/2)). Hence the sharp uniform Hilbert contraction coefficient is

    tau_D=(binom(D,floor(D/2))-1)/(binom(D,floor(D/2))+1).

For D=7 it is exactly 17/18. Two disjoint D-chains attain the bound. This removes a previously plausible obstruction: no special bridge decomposition is needed merely to obtain positive mixing.

The [four-collar balance theorem](poset_four_collar_boundary_20261019.md) proves an infinite width-three, seven-thin family with a designated pair of balance at least 101/236. The [finite-cone framework](finite_cone_balance_route_20261019.md) states the exact coverage, invariance, and marked-event conditions a larger proof needs. A strict 1/100 witness margin survives arbitrary exterior weights after a conservative 742-rank buffer on each side at D=7.

The [algebraic boundary theorem](algebraic_boundary_certificates_20261019.md) repairs one-sided periodic contexts using exact algebraic Perron rays. Already a natural four-thin family requires a quartic field. Full-degree Perron eigenvalues force rational threshold equality to be entrywise equality for fixed rational marked rows; otherwise strict algebraic separation gives a finite tail certificate. Two-sided stationary equality and complex subdominant modes remain distinct obstructions.

**Missing global step:** a family of exact boundary cones and actual comparison witnesses covering all legal nonperiodic seven-thin presentations. Positive mixing alone does not supply that coverage. Fibonacci bulk comparisons remain below 1/3 while the finite boundary supplies the witness.

## 3. Gold budgets, the golden Collatz carrier, and signed Golay action

The [finite budget bridge](gold_budget_golden_golay_20261019.md) starts from the actual Gold obligation g(x)=(1-2x)/(1-x). On P¹(F_11) minus one fixed point, the coordinate z=(x-5)/(x-9) gives

    g : z→2z,       g² : z→4z.

The second is exactly the inherited eleven-point standard-Collatz period-five carrier. The one-step budget multiplier requires the signed affine group of order 110; the two-step subgroup has order 55.

The [signed relation/Golay construction](golden_signed_relation_geometry_20261019.md) gives the equivariant syndrome map. One budget step exchanges the reciprocal five-dimensional modules; on a single module it can be represented by a nonlinear canonical decoder. The largest linear subspace on which the decoder reflection is linear has dimension exactly two.

**Essential boundary:** every projective residue occurs as both a balanced and an unbalanced comparison probability in genuine finite posets. The finite conjugacy organizes a carrier; exact real count inequalities still need separate transport.

## 4. Integer Collatz: distinguish an affine identity from a point stabilizer

The [integer collision atlas](integer_stabilizer_20261019.md) classifies rational stabilizers and shows that an actual integer merge can remain in a nonidentity affine state forever. The pair (45,15) illustrates the dependence on the chosen affine presentation. Literal merging itself is presentation invariant.

For translated sources x+d and x at equal Terras depth N, unequal odd-step counts give isolated roots with

    x < d/2 + 2^N/12.

The coefficient 1/12 is asymptotically sharp uniformly in d, even for lawful positive integers. Above that finite head, any merge must be a whole-cylinder affine identity. The rational root denominator contains the exact period factor 3^j-1. Unequal clocks instead introduce determinants involving both powers of two and three; the 223/233 receipt has the factor 2^8-3^5=13.

**Missing global step:** source-wise coverage by parity-valid, positive, well-founded join receipts. Haar-almost-sure coalescence, odd-modulus agreement, and a finite-height atlas do not establish that coverage.

## 5. Square sums: an infinite rational chain with exact Pythagorean ancestry

The [boundary-start lemma](square_sum_growing_chain_bit_bound_20261019.md) supplies arbitrarily small contact roots from a prescribed seed with a quantitative signed elliptic-translation budget. It removes the need to start at the old seed coset or rebuild an atlas at every step.

The [cubic contact theorem](square_sum_cubic_contact_chains_20261019.md) alternates ++ and -- contacts. Their incoming and outgoing ports are opposite at every interior petal, so the AP-paired cycle transfer still works. Two cubic steps have cusp expansion h→h-4h³+O(h⁴), preserving a fixed chart for an infinite rational sequence.

The exact primitive arithmetic is

    C_(i+1)=C_i² R_i, gcd(C_i,R_i)=1,
    (7/8)C_i³ < C_(i+1) < (9/8)C_i³.

A Gaussian factorization identifies the new factor with an outer Euclidean/Pythagorean child or parent norm. Existing hypotenuse-prime valuations double and the new factors are coprime to the old hypotenuse. This is a map-level link to the ternary-tree work, not just a shared list of numbers.

One polynomial-in-seed-bit-length initialization avoids the seed forever. Every finite k-prefix gives an integer square-sum cycle through the square-scaled seed with O(3^k L³ log⁴(L+1)) label bits. Each prefix has its own nested scale; no single finite integer center serves the infinite sequence. For a fixed initialized orbit, bit growth is Theta(3^k), while spreads satisfy f_i sqrt(i)→1.

The [balanced contact-tree theorem](square_sum_balanced_contact_trees_20261019.md) uses two free block maps distinguished by cubic jets. Trees of logarithmic depth admit polynomial-bit embeddings after one subdivision per tree edge and a seed stem. Arbitrarily long chains remain exponentially costly in these constructions.

The [capacity theorem](square_sum_common_center_capacity_20261019.md) independently proves that polynomial numerical labels in the number of common-center petals are impossible. Polynomial bit length remains a meaningful open target.

## 6. The square-sum/graceful/Collatz analogy made exact

The [signed affine labeling note](signed_affine_tree_labelling_20261019.md) gives an exact switching and feasibility dictionary. Square-sum tree edges are reflections; a bipartite sign switch converts them to signed square differences. Consecutive coverage becomes a root-interval intersection with explicitly forbidden collision values. Graceful tree labeling becomes a distinct integer potential with span exactly the edge count and a permutation edge-difference palette.

Signed Collatz words obey the same affine cycle equation, with different parity, integrality, and positivity guards. The useful common target is an admissibility-preserving move or a well-founded defect, rather than an unproved isomorphism of the conjectures.

**Next packing obstacle:** sparse square-sum embeddings do not preserve a consecutive interval, and graceful exchange connectivity does not ensure every prescribed tree shape has a labeling. Future moves must retain the appropriate palette and control vertex-potential collisions.

## 7. Current proof targets

1. Extend robust poset witnesses to genuinely aperiodic multiport interfaces, including bridgeless incomparability graphs.
2. Extend the proved Gold injection schema beyond a single first-arrival attachment while preserving exact deficits.
3. Classify lower-degree rational contact maps; the cubic construction is effective, but degree minimality is not yet proved.
4. Exploit fixed-translation Collatz branch structure to obtain source-wise descent beyond a growing finite atlas.
5. Keep interval coverage and prescribed tree shape as explicit constraints in every proposed square-sum/graceful transfer.

Each target has a specified obstacle and carrier. None is inferred merely from numerical resonance or a common finite quotient.
