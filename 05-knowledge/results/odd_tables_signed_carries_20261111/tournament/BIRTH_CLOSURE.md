# The remaining two paths are not supplied by one unrestricted birth

**FINITE-EXACT, 2026-11-11.** This controls one concrete source supplied
by the concurrent native patch investigation. It is not an obstruction
to other sources, more vertices, or the desired scalar count5993.

Start with parent rows

    (162,1020,1017,1009,993,960,897,768,513,1).

Starting at root0, append five vertices successively, each beaten only
by the preceding root, and retain the last vertex as root14. Reverse
the two arcs on unordered pairs `{0,10}` and `{0,11}`. The resulting
fifteen-vertex tournament is

    (2210,1020,1017,1009,993,960,897,768,513,1,
     3071,5118,10239,20479,8191).

Its Hamilton count is5991 and root14 still has indegree one. The
concurrent patch supplies its all-placeholder formula; this note uses
only the literal order15 specialization and asks for one more vertex
that raises its count to5993.

For every old subset A, endpoint DP computes paths ending at i and
starting at j. Define E_i,S_j on the full old set, and let W_ij count
ordered spanning two-path covers whose first path ends at i and second
path starts at j, both nonempty. For a new vertex beating exactly the
old set O, the exact count is

    sum_(i notin O) E_i + sum_(j in O) S_j
      + sum_(i notin O,j in O) W_ij.

Split a full path immediately before and after the new vertex: the
two empty-side cases give the first two sums and the nonempty-side
case gives the last. Thus a single subset-DP precomputation evaluates
all cuts without rerunning full path enumeration on each child.

The complete32768-cut test finds **no5993 child**. It does not impose
root-indegree retention, so the smaller class preserving the root's
indegree one is excluded too. The nearest count is5991, attained by
the global-sink cut O=empty. This closes only the proposed one-birth
completion of this particular5991 patch. Since there is no hit at
placeholder size1, no uniform placeholder identity is inferred.

Reproduce with

    python3 -B 05-knowledge/results/odd_tables_signed_carries_20261111/tournament/one_birth_close.py

The wrapper compiles `one_birth_close.cpp` in a temporary directory;
`one_birth_close_receipt.json` retains the source rows, root, universe,
target and result. No binary is stored in the repository.
