# A binary move with a frozen ternary carrier: what is and is not possible

## History braid

* The September [ternary/Berggren bridge](../ternary_bridge_20260925.md) keeps the signed even leg and odd denominator, not just a residue. Its sources include THM-3333/3756 and [Berggren edge transport](../collatz_mod6_20260921_berggren_edge_transport.md), where exact halving requires the numerator/denominator pair. [September arithmetic braids](../arithmetic_braids_20260917_collatz.md) separates a genuine 3-adic isometry from ordinary size descent; [reverse-tree pieces](../collatz_mod6_20260917_reverse_tree_pieces.md) and [inverse completion](../arithmetic_braids2_20260917_inverse_completion.md) give neighboring versions.
* October [binary–ternary guard fusion](../collatz_binary_ternary_guard_fusion_20261004.md) combines the [debt-word](../debt_word_join_search_20261004.md) and [joint-home address](../joint_home_address_20261004.md) experiments; exact binary words and ternary sibling addresses solve separate guard problems. [Join shields and lifts](../collatz_join_shields_and_lifts_20261004.md) and [inverse-ray ternary addresses](../inverse_ray_ternary_addresses_20261004.md) retain the carry and the actual node (source row, block height, parent). Residue compatibility never supplies an absolute lift by itself.
* [Ternary ports and golden defect](../ternary_ports_golden_defect_20261018j.md) connects graph-Pythagorean synthesis and earlier [Tait/Kempe triples](../procgen_kuratowski_20260925_tait_kempe_triples.md): nonempty boundary supports can be disjoint. [Ternary Frobenius reconstruction](../ternary_frobenius_reconstruction_20261019/TERNARY_FROBENIUS_RECONSTRUCTION.md) tracks cyclic order three via its norm operator; a formal residue defect is different from a true lost kernel. These point in opposite directions here: a constraint may destroy an entire route, or isolate the smaller space in which a real route exists.
* [Smooth interfaces](../adaptive_budgets_smooth_interfaces_20261020/README.md) and its [connections](../adaptive_budgets_smooth_interfaces_20261020/CONNECTIONS.md) follow [prime-batch triangles](../prime_batch_triangle_20261020/README.md) and the June THM-516 Burnside core. The gcd-response matrix factors as divisor zeta times diagonal totient times zeta transpose. Prime-only probes lose composite prime powers; Farey moves must retain content as well as normalized slope. This suggests studying the lattice of physical observations rather than separately studying their projections.
* Recent [Bott/Cayley slots](../chanpak_bott_cayley_periods_20261115/BOTT.md) supplies slot reciprocity and diagonal factorization. [Physical routes](../chanpak_physical_routes_20261115/PORT_ROUTES.md) supply literal readers, [equal-count routes](../chanpak_physical_routes_20261115/EQUAL_Q.md) and the [CRT firewall](../chanpak_physical_routes_20261115/CRT_FIREWALL.md) warn against assembling marginal coverage, and [Fibonacci repair](../chanpak_graceful_ternary_20261115/FIBONACCI_REPAIR.md) has 00/10/01 choices with a retained cross-block guard. The construction below moves a binary digit while freezing the larger ternary carrier, rather than merely comparing separate attainable residue sets.

## Refine the obstruction before looking for a witness

Let q_T(t) be the count of ordered t-slot covers of T by directed paths (empty slots allowed). Concatenate the slots to a permutation w. A cut between slots is mandatory at each backwards adjacent pair. If d(w) is the number of those pairs, its contribution is binom(t+n−1−d(w),n), by stars and bars with compulsory cuts. Reversing w replaces d by n−1−d, giving q_T(−t)=(−1)^n q_T(t). The universal integer-coefficient odd-cycle expansion (as in the Bott/Cayley note, via Grinberg–Stanley) also shows that q_T is monic and its other coefficients are even.

For positive a,b, q_T(a)≡(−1)^n q_T(b) mod (a+b). In particular **freezing q_T(3) allows at most two residues of q_T(2) modulo 5 even when order and parity vary**. A three-term arithmetic progression of distinct consecutive moves with stride coprime to 5 is impossible on that fiber. This excludes a uniform three-reading run at stride 3·2^k, more generally any 2/3-smooth stride, regardless of how many tournaments are allowed. More generally, three consecutive elements of an AP on a frozen q(b) fiber require that the additive order of the stride modulo a+b be at most two.

At fixed order, put ε=1 for odd order and ε=2 for even order. If both H=q(1) and q(3) are frozen, the difference polynomial is divisible in 2Z[t] by t^ε(t²−1)(t²−9), hence Δq(2) is a multiple of 60 or 120 respectively. If H and q(2) are frozen the same argument gives Δq(3) a multiple of 240 or 720 respectively. In both cases the proof is monic polynomial division after extracting the common factor 2. Thus a nontorsion horizontal dyadic clock must either tolerate a factor 5 or change its ternary carrier. Padding by a one-way transitive tail of size s multiplies q(2) by 2^s and q(3) by 3^s independently, so a defect in the latter becomes the larger-order term rather than an error that can be neglected.

## Physical separation by two-factor exchange

Use 8-vertex adjacency row codes with bit 8i+j indicating i→j. Exact profiles:

| Code | H | q(2) | q(3) |
|---:|---:|---:|---:|
| 4584474270165959682 | 77 | 2024 | 21789 |
| 9169081464381770758 | 15 | 960 | 14175 |
| 9160127196319778056 | 35 | 1352 | 16947 |
| 9169072685987728394 | 33 | 1440 | 18225 |

Take the ordinal sum D0 of rows 1–2 and D1 of rows 3–4. Optional-slot covers factor under ordinal sum. Both modules therefore have H=1155 and q(3)=308859075; their q(2) readings are 1943040 and 1946880. The 3840 shift is exactly 15·2^8. The equalities are genuinely simultaneous: the products of 77·15 and 35·33 agree, and the q(3) ratios on the two blocks exactly cancel (21789/16947=18225/14175=9/7). The q(2) ratios almost, but do not exactly, cancel. The modules differ by 11 specified internal arc reversals (7+4), treated together, with the interblock orientation unchanged. No claim is made about count control during individual internal flips.

Append a common transitive tail of size s to each module, and substitute it into the physical three-reader route from the prior session (skeletons 29130→29002→29000 at port 5). The companion script derives this route independently via external gap enumeration rather than using the tabulated polynomials as axioms. For each s≥0 the six strong tournaments fit an exact rectangle:
   
    H(s,j,c) = 4·308859075·3^s + 9·Q_j·2^s + 1155c,
    Q_0=1943040, Q_1=1946880,  c∈{−3,−1,1}.

Vertical steps are single exterior arc flips, each of size 2310. Horizontal steps are the same 11-edge simultaneous switch in every row, of size 34560·2^s=135·2^(s+8). These edits are on disjoint literal arcs; the mixed difference vanishes identically. The modules remain nonempty, and substituting a nonempty module into a strong skeleton preserves strong connectivity. Horizontal stride is Θ(H^log_3 2) with unbounded 2-adic valuation, not Θ(H) as in a pure dyadic tower. The common divisor of the two steps is 30; it is a divisor, not a covered fine interval.

## Many independent physical switches, but a curved image

There is an exact higher-dimensional lift, without any more reader ports. Ordinally concatenate m copies, independently choosing D0 or D1 in each copy, before adding the common tail. In every row the m physical switches commute on disjoint arcs and are all available from the start. For k selected copies, the total is

    4·308859075^m·3^s + 9·3840^m·506^(m−k)·507^k·2^s + 1155^m c.

This supplies m successive positive horizontal moves for arbitrarily large m, with H and q(3) held fixed *during that sweep*. But the 2^m states collapse to m+1 horizontal levels, and the steps form a geometric progression of ratio 507/506 rather than an equal grid. In particular, treating an independent Boolean choice as an independent additive digit gives a false window certificate. The mixed second difference under two independent horizontal switches is exactly 9·3840²·2^s times the product of the untouched q(2) readings. This gives an explicit positive error term to cancel in any attempt to flatten the strip. The reader flips continue to give a constant step 2·1155^m at each level.

One way to search for such pairs is to view each physical profile as a point in the multiplicative monoid of scan values. Seek a relation on the H and q(3) coordinates that does not persist in q(2); any such relation lifts to an actual ordinal product. Then seek other relations whose simultaneous legal states interact without destroying the reader. This differs from finding two independent residues that cannot coexist on the same host.

**Limits.** Each individual block gives one binary choice. The many-block lift gives arbitrarily long sweeps only by supplying more physical blocks in advance, with the steps no longer equal. Repeating the same switch returns to the original state. It gives neither a pure 2/3-smooth horizontal stride nor placement at a prescribed absolute target. A promising next exact obligation is to construct a *different* physical correction varying q(3) in coordination with q(2), or use two non-equivalent readers so that the forced mod-5 coupling is carried through the next state. The 2-versus-3 distinction here is structural: binary growth programs a digit, ternary growth transports the large baseline, and evaluation reciprocity forces a third prime into any repair that insists on keeping the two independent. This local obstruction and physical witness do not prove Chan–Pak.
