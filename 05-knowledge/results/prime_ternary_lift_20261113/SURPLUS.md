# Small rooted surpluses: what dropping indegree one actually buys

2026-11-13, Codex. **PROVED** padding lemma; **FINITE-EXACT** census;
all-order small-count conclusion additionally uses the **CITED** Busch
minimum. The degree bound at the end is **OPEN**, not a dependency.

Let F(T,r) count ordered pairs of paths partitioning T, the first
nonempty and ending at r, the second possibly empty. Put delta=F-H.
Cutting a Hamilton path just after r injects its H paths into these
pairs. The remaining pairs are exactly those with nonempty second path
whose first vertex beats r. Thus delta>=0 counts broken joins, not all
two-path covers indiscriminately.

The inherited triangle context gives H'=H+2F=3H+2delta. It can implement
the prime-lane digits2 and4 precisely at delta1 and2. Unlike the earlier
paired-repair grammar, no indegree-one restriction is required here.
Explicitly, append a,b with r->a->b->r, both a and b beating every old
vertex other than r. This is the identified-root context of the inherited
[endpoint compiler](../endpoint_context_recursion_20261113/port_compiler.py).
The verifier constructs these literal arcs at all85 rooted controls.

## Exact head and an all-order padding argument

For an ordinal decomposition T=L -> C -> R, where C contains r,

    H(T)=H(L)H(C)H(R),
    F(T,r)=P_L(2) F(C,r) H(R).                        (1)

Here P_L(2) counts ordered two-slot covers allowing empty slots, with
P_empty(2)=1 and H(empty)=1. Each of the two paths takes its L portion
first; only the second can have an R portion, since the first ends at r.
Restricting to the three blocks and concatenating back proves the
bijection in(1). The blocks need not be strong for this identity.

Now take C to be the strong component containing r. If L is nonempty,
P_L(2)>=2H(L), so delta(T,r)>=H(T). For H>=3 this excludes delta1 or2.
Hence such a root must be in the first strong component. Then

    delta(T,r) = delta(C,r) H(R).

Since H(R) is positive odd, delta1 or2 forces H(R)=1. A tournament with
one Hamilton path is transitive (equivalently the inherited cyclic-triple
exclusion for H=1), so R is only transitive suffix padding. Consequently
for H>=3 the existence of delta1 or2 is determined entirely by a
strong core of that same count.

[THM-4744, interval chains and strong inventories](../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md)
records Busch's all-order minimum, including f(9)=75 and its increasing
continuation. Therefore, for **every attained odd h from3 through73**,
the complete order-eight census decides at all tournament orders whether
delta1 or2 exists. This is not an all-order classification for larger h.

The ordinary counts at orders2..8 are1,2,4,12,56,456,6880. Every vertex
occurrence was checked, without identifying equivalent roots. There are
307 occurrences with delta1 or2:

| delta | root indegree | occurrences |
|---:|---:|---:|
| 1 | 1 | 145 |
| 2 | 1 | 115 |
| 2 | 2 | 47 |

Thus discarding indegree one gives47 additional rooted occurrences,
all at digit4 in this universe. It does not remove all composition gaps:
neither digit exists at h=3,15,27,35,41,43,57,67. Some other counts have
only one digit. The complete lists and witnesses are in `surplus.json.gz`
and `surplus_audit.json`. For h=1, a transitive root of indegree k has
F=2^k, so delta1 exists and delta2 does not.

For example,15->47 cannot be produced by this direct triangle context
from any H15 tournament at any order. This excludes that operation,
not47 itself or a core-changing construction. The new full-response
catalog gives alternatives for some failed direct exits, including the
all-m mixed-mode family starting41->125. `triage.json` records repairs
in the prime/7p lanes from the normalized dyadic catalog, including their
exponent thresholds; it excludes49, which is already a finite atom seed.
Empty repair lists are only absences in that catalog.

## Niche question left open

Every tested root through order eight satisfies

    delta(T,r) >= 2^(indegree(r)) - 2.

The observed minima for indegrees0..7 are0,0,2,6,14,30,62,127. This suggests
an excursion or path-partition injection. The naive bound F>=2^k does
not prove the claim, because H must still be subtracted. No valid uniform
injection was obtained here; the assertion at arbitrary orders remains
**OPEN**. In particular the finite evidence does not prove that delta1
always forces indegree one.

Independent permutation controls extract F at all85 root occurrences
through order five. Literal ordinal joins check(1) with empty, transitive
and cyclic blocks. These controls supplement the proof rather than
replace it.
