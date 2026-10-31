# New small bases and the missing B51

The eleven positive base certificates in SMALL_BASES.json are independently
checkable finite witnesses. Their proof status does not depend on a search
algorithm or on search completeness.

B42–B50 and B52 came from the parent session's forward/reverse constrained DFS.
The final missing B51 was recovered by independently enumerating paths on 51
and storing their common-vertex positional-parity masks, then looking up masks
from paths on 52. The compiled discovery program reported a complete first
side traversal of 14,314,732 recursive nodes, 24,777 paths, and 12,006 distinct
profiles; the second traversal stopped at a witness after 3,206,372 nodes and
7,165 paths. These are discovery diagnostics, not a separate audited counting
theorem or a necessary dependency of the positive theorem.

The positive certificate is short: both paths begin
1,3,46,18,7,2,47,34,15,21,43,6,30,19,
then take different suffixes with the same position parity for every label
through 51. All vertices and square equations are replayed by both checkers.
The repair is to change the entire recursive base interval and cutoff, not to
append small seed files while retaining the old 2620 cutoff.
