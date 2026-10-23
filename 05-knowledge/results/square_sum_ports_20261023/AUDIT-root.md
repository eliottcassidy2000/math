# Independent publication review

23 October 2026. Reviewer: root. **PASS for the stated finite universe and
elementary all-size identities**, including the user's explicit 3--12 to
4--12 correction. The final review also repairs a false limitation: a
Hamiltonian cycle supplies both inclusion and omission for every single edge,
with unrestricted endpoints. See the correction below; joint requirements
and exact frequencies are not thereby classified.

The arrival decomposition has no missing factor of two: an interior newcomer
determines one unordered neighbor pair and one cover by two components;
fixing the lesser joining owner in the first subset makes the endpoint
convolution unique. Singleton path components use the stated endpoint
convention. The old-path and two-cover cases are disjoint and exhaustive.

For a proper S, deletion of S from a spanning path leaves
|S|+1-b_P(S)-e_P(S) components. Additional graph edges can only merge these,
which gives the inequality in the correct direction. Tight-cut searches need
no sizes above floor((n-1)/2). For the marked G25 relaxation, a selected
1--15 and endpoints18,t spend at most three units; sizes at least14 have
budget at least4. Therefore the declared search through size13 is complete
for that relaxation, not an unproved small-cut assumption.

The final exclusion certificate does not depend on trusting the full path
inventory: the explicit 21 cuts eliminate their endpoint choices, and the
remaining three deterministic degree propagations force a deficit or the
proper cycle 2--14--22--3--13--23--2. Each cycle edge is a literal square-sum
edge. A proper selected cycle cannot lie in a Hamiltonian path.

I independently replayed verify.py (28,245 checks) and guarded_forbidden.py
(56 checks). The first reconstructs full paths via edge subsets and compares
them with the vertex-order DFS output. The second needs no path enumeration.
The construction agent separately replayed the exhaustive C++ cut universes,
the subset-convolution recurrence, and normal/optimized Python consumers.

Publication found a byte-level reproducibility issue: CRLF source hashes would
change under the repository's LF checkout policy. The owned packet was repaired
to specify LF-normalized hashing and fixture comparison, retaining all
mathematical contents. This is a transport-format repair, not a change to a
proof or finite result. The full reproduction includes a CRLF round-trip
control; the manifest records the precise normalization.

**Mathematical correction before publication.** The draft said a Hamiltonian
cycle did not establish that every noncycle edge occurs in a spanning path.
For any chord uv, remove the two forward cycle edges at u,v and add uv. The
two resulting components join into one spanning path containing that edge.
Cycle edges are included by removing another cycle edge. Omission witnesses
come from ordinary cycle cuts. Thus every edge in a Hamiltonian undirected
graph is optional among unrestricted Hamiltonian paths. This was independently
confirmed by hall_bott_recovery. The construction agent added literal inclusion
witnesses for all46 G32 edges and all cycle/chord positions through order64.
The all-n>=32 square-sum consequence uses the inherited CITED cycle theorem,
also pinned in Robert Gerbicz's own 2018 OEIS A090461 comment; its full
construction was not re-audited here. The failed implication and scope are
recorded in 01-canon/MISTAKES.md.
