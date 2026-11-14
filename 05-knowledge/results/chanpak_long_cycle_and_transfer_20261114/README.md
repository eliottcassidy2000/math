# Long Chan–Pak session: complete cycle encoding and sharper ordinary-path gates

**Chan–Pak remains open.** This session proves a complete neighboring
spectrum, constructs an exact escape for the known failed terminal,
and rules out two broad ways of trying to upgrade local control
into ordinary absolute coverage.

1. **[Every positive Hamiltonian cycle count, at logarithmic order](BINARY_CYCLE_GATE.md).**
   Two hubs over a transitive spine encode any \(K>0\) directly in
   binary, with explicit strong witnesses. At a fixed order \(n+4\)
   the construction fills \(1,\ldots,2^n\). The cycle count is the
   linear coefficient of the odd-order color polynomial divided
   by two, and is recoverable by a short signed list of ordinary
   path counts. No bounded fixed linear list of color evaluations
   (even with signs) can do this universally. This closes the
   encoding/readout issue, not the one-tournament path realization.
2. **[Exact ordinary path counts of those cycle encoders](CYCLE_DIGIT_GAP.md).**
   Removing one marked vertex gives a path-count difference
   made of three explicit digit sums. A leading-zero change cancels
   the nuisance terms in a signed combination of **four** literal
   ordinary counts, yielding exactly \(2K\). The whole ordinary
   count also has an exact dyadic padding ray; for cycle counts
   11, 13, and 26 these rays show three distinct width phases for
   divisibility by 7. The cycle-13 source gives the ordinary
   prime 397 and \(749=7\cdot107\) at consecutive paddings.
3. **[Changing the terminal does repair the exact old miss](MULTIPRIME_TERMINAL_LIFT.md).**
   The failure of every terminal length 7 at 180,131 does not
   preclude ready geometry itself. A literal 17-vertex ready
   construction with length 3 gives exactly 180,131. More broadly,
   eight actual roots and two terminal clocks realize every odd
   cofactor residue simultaneously at arbitrary powers of
   \(2,3,5,7\), with exact division by 7. This clock grammar still
   has only \(O((\log X)^2)\) distinct outputs through height \(X\).
4. **[Nonlinear apex controls do not cover on their own](APEX_RENEWAL.md).**
   The *entire* one-apex family over a transitive spine misses
   at least \(1-\log2\) of odd counts asymptotically. The width
   phases cancel in the leading term, not just at dyadic endpoints.
   A [complete literal inverse](HALF_INVERSE.md) for this family
   tests all orders at \(O(\sqrt N\log^2 N)\) arithmetic cost.
   This says nothing directly about how many primes it misses,
   or what higher-interface constructions supply.
5. **[A shared-root chain can carry arbitrary tournament blocks](SHARED_ROOT_CHAIN.md).**
   A reusable nonnegative three-state transfer handles arbitrarily
   many physical rooted modules. It gives exact ordinary repair
   rails and a dyadic terminal displacement. But with a fixed
   singleton reader, any continuation except one special singleton
   overshoots \(2H+1\); both binary digits cannot be implemented
   by appending blocks while preserving that reader. The obstruction
   is all-order, not a short word search.
6. **[Near-readiness has an exponential root penalty](ROOT_DEFECT.md).**
   The old outdegree-one theorem becomes a quantitative bound:
   \(R_+-H\ge H(A)(2^{d-1}-1)\) for root outdegree \(d\) and the
   earlier ordinal part \(A\). Explicit examples show exponential
   scale base two is unavoidable.

Each note points to its literal or exact arithmetic replay. The
remaining requirement is not another congruence: either physically
absorb the signed excess channels from the cycle encoder, or
supply a target-selectable ordinary return that changes its current
source or terminal reader. The fixed families whose limitations
are proved here cannot supply that return by themselves.

**Sharper follow-up.** [Two induced deletions of the same binary cycle
encoder](POSITIVE_CYCLE_GAP.md) already give an ordinary path gap
exactly (2K), with strong graphs of orders (lfloorlog_2 Kfloor+3)
and (+4), upper count at most (4K+1). No signed combination is
needed for the gap. The lower count is not freely selectable:
that distinction is still what blocks the Chan–Pak inference.
