# Independent Artin–Schreier supplement review

2026-10-20 UTC. **PASS**, including the final odd-degree restriction added to §2.

Reviewed complete source: `session10-artin/ARTIN_SCHREIER_AND_CRITICAL_PERIODS.md`, SHA-256 `952b8a3ae29f5737f2015ee920b70fad3635585c574b40e5667b5959020ae6f7`.

## Finite-field fibers and exceptional elements

For b of exact degree d, every root of a²+a=b satisfies sigma^d(a)+a=t, where t is b's minimal-field trace. If t=0, the roots lie in F_(2^d); their degrees are exactly d because F_2(b) is contained in F_2(a). If t=1, the roots have degree exactly 2d: their degree is a multiple of d, divides 2d, and does not divide d.

The two trace-zero roots cannot belong to one Galois orbit. A Frobenius power interchanging them would fix b, hence have exponent divisible by d, which already fixes both roots. For trace one, sigma^d interchanges the roots. This proves the orbit-indegree classification, including the distinction between ambient trace and minimal-field trace.

Zero has ordinary degree one and minimal polynomial X. Its predecessors are the distinct orbits of 0 and 1, so it is a branch. The element 1 has degree one and trace one; it is a source for odd ambient degree and has one predecessor orbit for even ambient degree. No degree-zero convention or omitted exceptional orbit is needed.

## Intrinsic stable subspace and restricted count identity

For m=e u with e a power of two and u odd, sigma^e fixes exactly L=F_(2^e). Thus im(sigma^e+I) has dimension m-e. It lies in the relative-trace kernel by telescoping. The relative trace is the identity on L, so is surjective and has kernel of the same dimension. Since D^e=sigma^e+I, the three stated descriptions of W coincide.

D preserves W and is injective there because 1 has relative trace 1. Hence it permutes W, and every element enters W in at most e derivatives. This supplies exactly the hypotheses needed for the previously audited routing theorem. Minimal-polynomial ordering is invariant under finite-field isomorphisms and gives a well-defined intrinsic tie-breaker; it need not agree with a normal-basis necklace ordering.

The final restriction is valid. Fix odd r and take m=2^a r. Within the preserved odd-degree class r, an ambient-trace-one source must have degree 2^a r: the repetition factor m/d must be odd. The trace-zero targets in that same class have degrees 2^j r for 0≤j≤a. Restricting the bijection therefore proves

    I_1(2^a r) = sum_(j=0)^a I_0(2^j r),

with I_0(1)=I_1(1)=1. It is a restriction of an actual bijection, not a deduction from equality of total counts alone.

## Exact critical period and characteristic boundary

The formula for P_n(c)=f_c^n(0) follows by induction. If P_n(c)=0, the identity P_n(c)²+P_n(c)=c^(2^n)+c forces d|n for c of exact degree d. At n=kd the sum is (k mod 2) times the minimal-field trace. These two facts establish the exact least return d or 2d, not merely an upper bound. They include c=0 with period one and c=1 with period two. The formal derivative of P_n is one, so the asserted separability and divisor-period decomposition are valid.

The note makes no inference about rational cycles in characteristic zero. In particular it does not confuse the trace-zero targets of several possible degrees with parameters all having one common exact period.

## Checker and citation scope

The checker was read in full. Its exhaustive nonzero-element test a^(2^m-1)=1 certifies that each supplied quotient is a field. It constructs the complete Galois-orbit graph, verifies the trace/indegree conditions and exact critical returns, and checks routing, inverses, and degree control. It uses explicit exceptions, so optimization does not remove checks.

The two frozen outputs are byte-identical and report 8,063 checks for m=1,...,8. They were inspected and compared, not rerun in this review. Finite evidence supports the implementation; the general conclusions above follow from the proof.

**CITED:** [Baker–Chen–Li–Qian, §2.2, Theorem 2.6](https://arxiv.org/html/2508.12924v2), checked directly on 2026-10-20, supports the normal-basis, trace-refined necklace correspondence and its attribution to Reutenauer. The real-dynamics part remains an external cited theorem; it was not independently reproved here. No novelty claim is made.

No Git or repository-checkout edits were made.
