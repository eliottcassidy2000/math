# Complete all-H census through eight vertices

**FINITE-EXACT, inherited census, newly compiled cumulative table.**

`ALL_H_THROUGH_8.csv` answers the finite-size question without restricting H
through73. It gives every odd H from1 through661 at each n=1,...,8, including
zero entries, and cumulative frequencies over positive sizes1,...,n.
All even H and all H>661 have zero frequency throughout these sizes.

Columns distinguish ordinary isomorphism classes, self-converse ordinary
classes, converse-merged classes, and labeled tournaments. Cumulative columns
add those same channels across sizes; they never include the auxiliary empty
tournament. Rows7 and21 explicitly vanish at every size.

The ordinary totals by order are1,1,2,4,12,56,456,6880. Across sizes1 through8,
there are7412 ordinary classes. These are complete finite counts, not estimates.

Reproduce with `python finite_census.py` or `python -O finite_census.py`.
The compiler reads only the self-contained frozen `finite_input.json`, whose
provenance links and source hash identify Session24's complete `extended_census`
input. No tournament enumeration or fresh validation of old claims is performed.
Constructor controls check ordinary totals, labeled totals2^(n choose2),
converse orbit arithmetic, and zero gaps. The normal and optimized receipts
match. CSV output uses LF line endings for identical API publication bytes.

This finite all-H table and the neighboring `build.py` answer different scopes:

- `finite_census.py`: **all H**, but sizes1 through8 only.
- `build.py`: **all sizes**, but only odd H1 through73 with its frozen inventory.

Order9 has known191536 ordinary classes, but its complete H-histogram was not
readily available in the frozen input used here. We do not invent its entries
or describe an order8 table as an order9 classification.
