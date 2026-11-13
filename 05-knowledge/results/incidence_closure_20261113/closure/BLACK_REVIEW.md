# Independent root review of the black-core boundary

**PASS.** The root read the entire proof and checker separately from their
producer, after independently reconstructing the full black-support component
sizes from the frozen ordinary/merged matrices. The native census remains the
same declared n3..7 source; this is not a new independently generated n7 census.

The bound K=2(a+b)>=2|a-b| embeds every signed edge in black support; equality
requires one channel zero. Every chiral node has positive black degree because
all its darts are black and its Hamilton count is positive. Thus every contracted
signed component is black-connected. Projection and path-lifting prove the exact
component-count theorem; sum over both endpoints retains doubled internal loops
and the full source mass. The attachment corollary needs the stated nonempty
connected signed sector. The synthetic disconnected example is correctly marked
as an abstract graph, without a tournament-realization assertion.

The root's separate component traversal agreed with [1,1], [1,2], [1,1,1,7],
[1,1,32], [1,1,1,1,268] at n3..7, and independently verified that all listed
isolates have black degree zero. The n5 path 2--8--4 has two physical black
edges of multiplicity2 through a self-converse middle node; it refutes only
the converse implication stated in the note. No all-order connectivity follows.

The checker was source-reviewed for matrix normalization, component universe,
connected-fibre requirement, empty-graph convention, row-mass preservation and
SC attachment. Its normal and optimized receipts match with 69,050 always-active
checks. No repair was required. The status-only promotion is separately recorded
in BLACK_PROMOTION.json so earlier review-manifest hashes remain fixed.

Reviewed prepromotion proof: `712a9616ebb5425c35107a00c3bb4b11017d4b758afcc1bd6fb1325509774616`.

Promoted proof: `7106701edb9f2287442e9355105f9792e7770842819ece881352f794cd786a73`.
