# Overnight frontier: realizable poset laws, elliptic seed dynamics, and the remaining coverage barriers

Session: 2026-10-19 UTC. This continues the [earlier connection board](liouville_20261019_connection_board.md). The earlier sparse AP-triangle, fixed-contact, and graceful-exchange results remain in force. The present work concentrates on new proof mechanisms and the obstacles that survive them.

**No general Collatz, Gold Partition, 1/3–2/3, graceful-tree, or consecutive square-sum result is claimed.** All-size statements below have specified carriers and explicit proofs. Historical novelty is not established merely by finding a proof in this session.

**Literature status.** Ordinary consecutive-interval square-sum existence is already established: Robert Gerbicz records paths for every n≥25 and cycles for every n≥32 in his January21,2018 [OEIS A090461 contribution](https://oeis.org/A090461), linking his deterministic construction. Any interval-coverage obstacle mentioned here concerns extending this particular structured construction or preserving additional prescribed data, not an open claim of ordinary existence.

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

The [capacity theorem](square_sum_common_center_capacity_20261019.md) independently proves that polynomial numerical labels in the number of common-center petals are impossible. Polynomial bit length remains an unresolved target for this common-center construction program.

## 6. The square-sum/graceful/Collatz analogy made exact

The [signed affine labeling note](signed_affine_tree_labelling_20261019.md) gives an exact switching and feasibility dictionary. Square-sum tree edges are reflections; a bipartite sign switch converts them to signed square differences. Consecutive coverage becomes a root-interval intersection with explicitly forbidden collision values. Graceful tree labeling becomes a distinct integer potential with span exactly the edge count and a permutation edge-difference palette.

Signed Collatz words obey the same affine cycle equation, with different parity, integrality, and positivity guards. The useful common target is an admissibility-preserving move or a well-founded defect, rather than an unproved isomorphism of the conjectures.

**Next packing obstacle:** sparse square-sum embeddings do not preserve a consecutive interval, and graceful exchange connectivity does not ensure every prescribed tree shape has a labeling. Future moves must retain the appropriate palette and control vertex-potential collisions.

## 7. Later overnight progress: aperiodic contexts and the exact Gold boundary

The [aperiodic core-word theorem](port_word_cone_aperiodic_balance_20261019.md) now permits arbitrary nonperiodic words in two cores and five seam types. The exact port-word cone proves balance at least5872/12875, exceeding9/20. It includes ordinal-indecomposable, bridgeless incomparability graphs of maximum degree7. The [Gold upgrade](gold_port_word_cone_20261019.md) supplies two-comparison strategies on the same class, with uniform normalized slack at least791/6015. These are specified infinite families, not a theorem for all seven-thin posets.

The [prefix-side-poset theorem](gold_prefix_side_poset_20261019.md) extends the Gold strategy to nested arbitrary side-poset attachments through a chain prefix. With r=m-j, its sharp universal region is

    m≤n, n²−(r+1)n−r(r−1)≤0.

The attachment-density ceiling tends to phi^(-2). Exact diagonal saturation forces j=a², m−1=b² with consecutive Fibonacci a,b of specified parity. The [Pell half-step theorem](gold_pell_collatz_halfstep_20261019.md) identifies the positive curved boundary with a norm−4 orbit. Modulo11, its positive step is multiplication by4 and its algebraic half-step multiplication by2, recovering the stated Collatz/Gold finite carrier. The half-step moves positive dimensions to formal negative ones; it is not an integer Collatz move or a proof-transfer mechanism by itself.

## 8. The cubic contact degree barrier is now proved

The [minimal-degree theorem](square_sum_cubic_contact_degree_minimality_20261019.md) classifies global rational contact identities: degree1 produces only the same unordered petal, degree2 is impossible, and the advancing cubic construction is therefore minimal. The exclusions have explicit polynomial identities and replayable standard-library certificates.

The [finite-catalogue height theorem](square_sum_finite_catalogue_height_obstruction_20261019.md) strengthens this: any distinct rational chain made from a fixed finite catalogue of advancing global contact maps has exponential height, even if its map choices are adaptive and nonperiodic. A finite rational change of coordinates cannot remove the degree barrier. Polynomial-bit common-center chains must use a genuinely different mechanism, such as unbounded catalogues or algebraic correspondences.

## 9. The requested seed has a non-torsion elliptic continuation

The [Pythagorean-seed interpolation](square_sum_ppt_seed_elliptic_interpolation_20261019.md) places482,3362,6242 at x(P),x(2P),x(3P) on

    y²=x³−482x²+6195836x−299007352,
    P=(482,51840).

The sequence continues158,57442,8198402/7921. Consecutive sums are chord-slope squares because the x² coefficient cancels x(P). Every primitive Pythagorean AP seed yields a nonsingular curve and a non-torsion section: the sixth multiple has an explicit nonintegral coordinate, so Nagell–Lutz supplies a uniform proof. A concrete outer Euclidean/Berggren parameter ray gives an infinite family with every multiple positive. This is an exact group-law/ternary-tree connection.

The [elementary sparse baseline](square_sum_elementary_sparse_path_cycle_bounds_20261019.md) is essential calibration. Arbitrary sparse square-sum paths and trees already admit greedy labels of O(N²) numerical size; a two-label closure gives cycles of the same order while retaining the exact seed. Elliptic interpolation adds arithmetic structure, not a stronger unrestricted existence or size bound. Square-scaling an elliptic prefix generally changes its initial integer labels.

The [elliptic translated-coordinate theorem](square_sum_elliptic_translation_obstruction_20261019.md) explains a boundary between these constructions. An affine combination A+B x(Q)+C x(Q+P), with nonzero B,C and non-two-torsion P, has an even divisor exactly in the chord-square case B=C and A=B x(P). The natural normalized-x Pythagorean contact construction misses that condition for every non-torsion shift. Its odd divisor survives every isogeny pullback. Faltings then implies only finitely many rational contacts for each fixed shift. Higher-degree parameter functions remain open.

## 10. Collatz source ranks: exact local control and exact recharge

The [fixed-translation descent theorem](fixed_translation_descent_20261019.md) improves the isolated-head bound to d/2+2^(N−c sqrt(N)+O_d(1)) and proves an actual smaller-source receipt when min(a,b)≤2x−d. In particular depths N≤2x−d suffice. The unresolved isolated configurations lie in a linearly growing small-source region, with too many odd events for that descent argument.

The [return-rank theorem](collatz_return_rank_obstruction_20261019.md) proves dual exact budgets for a word of length n with a odd steps and center p: forward repetition spends v_2(x−p) in increments n; inverse integer lifting spends v_3(y−p) in increments a. Every fixed or ultimately periodic inverse word code gives only finitely many distinct integer lifts of one fixed seed.

A lawful two-word bank can nevertheless recharge either local valuation rank arbitrarily on switching centers. Thus local ranks are not automatically one global rank. Real-height contraction works for that particular bank, but universal source entry remains unproved. Retaining positivity, parity cylinders, and the actual source is indispensable.

## 11. A genuine realizability constraint closes the rooted-five poset obstacle

The [rooted-five theorem](rooted_five_logconcave_balance_20261019.md) now handles arbitrary finite exteriors, with no degree or width bound. If I(x) union I(y) is exactly the five-point core with covers b<z<a, x<a, b<y, then one of x:y, x:z, z:y is balanced in the whole poset.

The decisive step is classical rank log-concavity applied to a new marker element. If R is the initial run of an ideal in a uniform internal exterior word and m=E(R), then its tail is submultiplicative and

    E binom(R,2) ≤ m².

The two exterior sides independently produce count parameters satisfying1≤c≤u≤c² and1≤d≤v≤d². Exact eleven-order counts then force a balanced comparison. A previously hostile mixture of individually legal collars violates v≤d², so it cannot arise from any actual exterior law. This closes the specific local obstruction by preserving a constraint that the arbitrary-mixture cone had discarded. It does not show this core is unavoidable in all seven-thin posets.

The [exact PPT positivity theorem](square_sum_ppt_exact_positivity_threshold_20261019.md) also replaces a sufficient parameter family with a complete criterion for the specified recurrence: a/b<beta, where beta is the upper root in(0,2−sqrt(3)) of t⁶+4t⁵+3t⁴−8t³+35t²−12t+1. Above beta, non-torsion real-circle dynamics forces negative labels infinitely often. The first-five-label golden inequality is weaker: the primitive triple427,1836,1885 has eight positive initial labels followed by a negative ninth.

The [guarded switching theorem](collatz_switching_height_certificate_20261019.md) supplies a common weighted-height potential for a finite graph precisely when all its directed multiplier cycles contract. It includes an explicitly expanding adjacent return J followed by a lawful descending reset on a specified residue cylinder. It also proves why a finite contracting graph cannot cover all lawful J-only prefixes. Guard coverage remains an independent mathematical obligation.

Finally, the [prescribed-edge absorber packet](square_sum_prescribed_edge_absorbers_20261019.md) gives the exact Pell insertion

    x+y=s²; u=24x+25y, v=25x+24y;
    u—v can be replaced by u—x—y—v.

A singleton version and Hall's matching condition give a precise simultaneous-absorption theorem, conditional on eligible hosts already occurring in the backbone. A two-interval construction supplies disjoint hosts and explicit cycles. Gerbicz's established interval construction is credited. Its bounded-endpoint affine lift does not automatically preserve the requested seed edges, so protected-edge states require a separate proof.

## 12. The latest checkpoints: all-size witness migration and the exact seed interval

The [all-size mixed-two-slave theorem](two_slave_all_common_balance_20261019.md) goes beyond the five-point class. For a critical pair with at most one exclusive neighbor on each side and c common incomparable neighbors, it proves

    sum_z [p(x<z)+p(z<y)] ≤ 4c/3

with arbitrary common-neighbor order and arbitrary exterior. Exact nonnegative prefix/suffix increments reduce the statement to three symbolic chain-threshold shuffle identities. Selecting a maximum-span pair in the transitive greater-than-two-thirds majority tournament then contradicts this inequality. The c=0 case reduces to a path incomparability graph, whose endpoint balance is a Fibonacci ratio. Thus a hypothetical counterexample's maximum-span roots must have at least two exclusive neighbors on one side; in the seven-thin case they have at most four common neighbors. Remaining asymmetric cases are not disposed of.

The [all-factorial-moment lemma](ideal_run_all_moments_20261019.md) extends the earlier ideal-run bound to E binom(R,k)≤(ER)^k, with sharp limiting examples from two-chain shuffles. The [core-only Gold obstruction](rooted_five_core_gold_obstruction_20261019.md) shows why a balance theorem cannot automatically be promoted to a Gold strategy on the same five endpoints: an actual eight-point poset defeats every such strategy, while exterior comparisons give a strict Gold certificate.

The [seed-preserving interval certificate](seeded_interval_cycle_6242_20261019/square_sum_seeded_interval_cycle_6242_20261019.md) now gives an explicit square-sum cycle on every integer1…6242 beginning482,3362,6242. A controlled split of49-lifted blocks first produces6247 vertices; a fixed24-move repair trace removes the excess without cutting the seed edges. The endpoint6242 is minimal simply because the seed contains6242. This is one exact finite witness, not a uniform theorem for all larger intervals. The [common-pair parity obstruction](square_sum_affine_lift_protected_edge_obstruction_20261019.md) explains why an unchanged intact-block recursion cannot simply preserve the first seed edge at every residue.

The [Berggren absorber forest](square_sum_berggren_absorber_bijection_20261019.md) supplies another exact link to the ternary-tree theme. The classical three matrices generate all positive norm-one triples a²+b²−c²=1 in three parity-sector trees, rooted at(1,1,1),(2,1,2),(1,2,2). Unique addresses give branch-preserving bijections with three copies of the primitive-Pythagorean tree. Non-boundary nodes are precisely the nondegenerate coefficient gadgets of the stated linear absorber form. The central branch gives the Pell coefficients(5,5,7),(29,29,41), while the seed's AAA address gives(8,31,32) in the even-first sector. No historical priority claim is attached to the norm-one tree.

The [odd adjacent-return theorem](collatz_odd_adjacent_return_20261019.md) disproves a proposed simplification: expanding adjacent returns can have an odd lower source with every prefix above its starting value. A fully guarded30-step example supplies arbitrarily long finite repetitions on varying source cylinders. A general expanding/reset lemma bounds the number of excluded initial lifts in a minimally contracting composite cylinder by ceil(c+p), independently of repetition count. The F/H reset families preserve adjacency and give actual smaller-source receipts; arbitrary source entry and insufficient-reset exits remain unresolved.

## 13. Complete conditional transformation and further poset exclusions

The [one-sided three-slave theorem](one_sided_three_slave_aggregate_20261019.md) proves the same aggregate4c/3 bound with no exclusive neighbors on one side and at most three on the other, including arbitrary orders and arbitrary exteriors. Its integer-boundary polynomial certificate is distinct from ordinary coefficientwise positivity. Tiny exact suffix certificates close the c=0 case. Four one-sided slaves already admit an actual-poset violation of this aggregate, although that poset still has a balanced pair; the method's boundary is explicit.

The [complete final split-lift](seeded_final_split_table_20261019/square_sum_complete_final_split_lift_20261019.md) now has verified certificates for all392 combinations of input parity, output residue, and cut-edge orientations. Given the precisely stated protected input nice pair, it produces an interval-spanning cycle with the seed for every residue24…72. The final transformation has no missing state. A uniform all-large-N theorem still requires the upstream protected input-state transformations and a complete consecutive base interval. The single6242 certificate does not replace those hypotheses.

## 14. Current proof targets

1. Extend maximum-span witness migration to the asymmetric multi-slave configurations left by the all-size mixed-two-slave theorem. Retain realizable boundary laws and actual comparison witnesses.
2. Extend the Gold prefix injection schema to more general attachment geometry while preserving exact count deficits.
3. Find port-compatible algebraic contact correspondences beyond the finite rational-map and normalized elliptic-coordinate barriers.
4. Prove source entry or a common descent potential for Collatz return banks despite center-relative rank recharge.
5. Complete the upstream protected input-state transformations and consecutive base interval for the uniform prescribed-seed construction. All392 final split states and the exact6242 witness are established; the remaining input-existence argument is separate.

Each target has a specified obstacle and carrier. None is inferred merely from numerical resonance or a common finite quotient.
