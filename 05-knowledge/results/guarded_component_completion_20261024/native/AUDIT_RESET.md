# Independent audit of the native minimum-reset compiler

Status: **VERIFIED proof review / FINITE-EXACT replay.** Reviewed
`../theory/NATIVE_MINIMUM_RESET.md`, `check_reset.py` and its frozen receipt.
Normal and optimized replay both pass30,886 checks. No substantive error found.

The source-specific dynamic program keeps exact selected count and the full
Q-charge. It charges each actual slot exactly once, using its supplied path
distance when its bit changes, including U/V choices in slots that are not
selected. Truncation at count three is safe because future increments are
nonnegative. Each binary assignment has a unique prefix decomposition, so
the minimum, minimizing multiplicity and lexicographically chosen bitmask
are correct. The independent comparison with all9,880 literal triples at
every requested charge checks all three outputs, not just attainability.

The fixed-route graph retains each FULL reservoir triple as its intermediate
vertex. This is necessary: pairwise fibre-distance minima alone need not
choose the same intermediate record. Minimizing the bottleneck first and
then the additive total with that cap fixed is sound. The checker also
enumerates all216 complete reservoir trajectories, independently confirming
the exact(peak,total) optimum(103,415), including the chosen tie-breaking path.
The greedy comparator is an actual complete path, with(103,427).

The native inputs match the twelve-control family's fixed six-arrival route
and literal endpoint records. The old core margin stays two, controlled
helper occupancy stays six, reservoir occupancy stays three, and total charge
stays134. Every optimized shared pin is evaluated on its embedded path owner,
not inferred from the bit or charge alone. The physical row replay includes
both source and target diagonal contributions. Kernel/resource completeness
and source/partner host gluing are explicitly inherited from the frozen native
packet, with exact hashes; they are not reproved by the min-plus program.

The result is correctly scoped to this fixed core route and the supplied full
menus, with the prescribed helper return. It does not improve the uniform120
bound to103 and does not claim optimality over other core routes. The stated
source-specific worst-debt figures concern the five encountered triples, not
all9,880 possible source triples. No old pin is relaxed and no competing
intermediate witness is silently concatenated.
