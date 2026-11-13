# Incoming two-clock constructions remain outside unit readiness under births

**PROVED, using the H5 hereditary criterion:** every member of the incoming
Pi and Psi arbitrary-module families, and every further induced extension,
has no unit-ready root. **FINITE-EXACT:** the two literal source profiles.
The incoming clock and one-arc identities remain credited to their source.

The incoming revision94cb3ea contains
`05-knowledge/results/clone_clocks_negative_collatz_20261113/CLOCKS.md`.
It was read directly from that Git blob before integration. Its full
commit, Git blob identity and SHA256 are pinned in `incoming_clocks.json`.
This is an immutable **historical** provenance pin, not a claim that every
hypothesis in that version remains current. Revision5e148f917a refuted its
stronger claim that the constant term can never be-1: products such as
(60*3^j-1)(2^k+1) are counterexamples. The corrected
[current CLOCKS section5](../../clone_clocks_negative_collatz_20261113/CLOCKS.md)
separates this **REFUTED** claim from its narrower still-open C2-line
hypothesis. No constant-term exclusion was used below; the Pi/Psi source
identities and their hereditary barriers are unaffected.
The source proves the actual arbitrary-module interfaces

    H(Pi[M,N])=60 P_M(3)P_N(2)+H(M)H(N),
    H(Psi[M,N])=160 P_M(3)P_N(2)+40P_M(3)H(N)-3H(M)H(N),

where P allows empty path bins. The source's complete run-skeleton
certificates justify these identities; this note does not rerun them.
For transitive positive module sizes j,k their counts are

    1+60*2^k*3^j,
    3^j(160*2^k+40)-3.

These are actual repeatable clone clocks. They do not by themselves
assert unit readiness at any root.

The independently counted base graphs are

    Pi =(462,508,120,496,481,449,384,260,4),
    Psi=(974,1020,568,1008,993,961,644,772,68,256).

Pi is also reconstructed independently from the incoming forward order
with exactly the four stated back arcs. Both bases are strong. Their
all-root profiles give

| base | H | minimum starts | minimum F-H |
|---|---:|---:|---:|
| Pi |361|6|29|
| Psi |1077|18|85|

Thus both satisfy the criterion in
[ROOTED_HEREDITY H5](ROOTED_HEREDITY.md). They have no unit-ready induced
extension at any count or order. Replacing their two module vertices by
**any nonempty tournaments** retains the base induced: choose one vertex
from each module. Therefore every member of both arbitrary-module
families, not merely their transitive specializations, inherits the same
all-count obstruction. Any further sequence of births retains that
forbidden induced base as well.

The incoming Pi transport reverses one old arc and changes the positive
module sizes(j,k) to(j+1,k-1), for k>=2, realizing H'=(3H-1)/2.
Both endpoints remain in the Pi family, hence retain this readiness
obstruction. This does not refute the physical arithmetic transport:
it identifies a separate invariant that the transport does not repair.
At k1 the asserted positive-block transport is unavailable, and this
note makes no claim about deleting the last block or changing the
quotient there. More general old-arc repairs can remove the induced
base and lie outside the birth-only obstruction.

This strengthens the incoming finite readiness screens for these two
specified families to an all-order, arbitrary-module statement. It does
not address every two-clock template, negate their Hamilton-count
coverage, or settle other interfaces. The proposed all-k indegree-one
unit-ready405 interface has separately been **REFUTED** by the complete
k17 exclusion reviewed in [K17_INCOMING_AUDIT.md](K17_INCOMING_AUDIT.md).

Reproduce `python -B .../geometry/incoming_clocks.py`, then with `-O`.
The referenced Git object must be present; no mutable working-tree copy
is used as the provenance source.
