# Independent actual-path flag audit

2026-10-23, Hall/Bott recovery subagent. **PASS**, with the scope in README:
all-size local transport and deletion compatibility, plus finite exact
five-vertex carrier tables. No small-observer completeness theorem follows.

The proof's load-bearing coordinate is the actual path-edge forest F_S.
Every pair in S is globally retained precisely when it lies in F_S;
pairwise reversal therefore commutes with restriction. Deleting a vertex
restricts this forest and does not insert a new edge. The map preserves its
admissible domain and is an involution. Relabelling commutes because the
marked arcs are transported together with the tournament.

I checked the following implementation details independently:

- In `canonical(t,f,r)`, permutations with `p[0]==r` put the old marked
  root in the new position zero. All 4! orders of the other vertices occur.
  This is rooted canonicalization, not canonicalization that forgets r.
- The forest stores unordered pairs, but their directions remain in the
  tournament word. The permutation's reversal mask acts on tournament
  bits only; thus a marked arc remains the same directed arc after a
  relabelling. Forced forward arcs in the initial path gauge may legitimately
  become backward arcs in a canonical representative.
- A zero gap inserts one outside vertex; a retained gap inserts none.
  The selected vertices span positions from zero to at most eight.
  Prescribing their arcs never changes an unmarked global path edge,
  and marked ones are forced forward. Therefore all 5,184 local states
  really have a global Hamiltonian-path realization of size at most nine.
  No Hamiltonian-path assumption is made about the selected five-set.
- The deletion calculation restricts the original pair mask. Its six-bit
  complement at order four does not compress path gaps into shortcuts.
- Rooted-forest-free and cardinality-only witnesses are admitted local
  states and hence use the same global realization proof.

Both normal and optimized producer replay passed with 72,583 checks.
The independently written `audit_check.py` imports no producer functions.
It constructs directed adjacency sets, reverses the unmarked directed arcs
literally, and canonicalizes by scanning all 120 vertex permutations rather
than using precomputed bit-transform tables. It checks all 31,104 source/
target assignments (six canonicalizations per admissible state) and matches
both complete saved quotient maps. It obtains 676 unrooted flags, 3,360
rooted flags, and 24 fixed flags in each case.

Reproduction:

    python -B 05-knowledge/results/metagraph_flag_transport_20261023/check.py
    python -B -O 05-knowledge/results/metagraph_flag_transport_20261023/check.py
    python -B 05-knowledge/results/metagraph_flag_transport_20261023/audit_check.py

The producer source hash replayed was
`0a6bd3de16102b083506110c03497d865fadcb695efa4ad0ad6faa333e060cfd`.
The exact transport rule is useful for authenticating target motifs on a
known source/path port. Additional separation evidence is still needed
when using those motifs to infer a target class from a coarser record.
