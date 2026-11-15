# Joint annulus and residue windows for the hyperbolic reader

**Proved, not pointwise coverage.** Fix (rge1) and any nonempty open fractional interval ((a,b)subset(1,2)), including one entirely in ((3/2,2)). In every sufficiently high band ((a2^K,b2^K)), **every** odd residue modulo (2^r) is realized by the *same fixed order-five hyperbolic quotient*, with at least (c_{a,b,r}K) different counts per residue. A uniform form allows (2^rle c_{a,b}log K) for existence. Nonetheless the family realizes only (O((log X)^2)) absolute counts up to (X).

## Construction and mechanism

Use the first hyperbolic quotient of THM-5133 with its second port a singleton. Its output at first-port profile ((h,Q)=(H,Pi_2)) is (6Q-h). The profiles ((11,100)) and ((5,32)) have literal cores from the [fixed-kernel replay](HYPERBOLIC_MOVE.md). Take the ordinal join of (j) copies of the first, optionally one copy of the second ((ein{0,1})), and a transitive tail of length (tge1). Both (H) and (Pi_2) multiply in ordinal joins: an ordinary path splits into one path in each side; a two-color assignment of the vertices and the two ordered paths split independently in each side. Hence the literal count is
[
N(e,j,t)=6,32^e100^j2^t-5^e11^j.
]
For (tge r), the residue is just (-5^e11^j). For (rge3), the powers of 11 generate precisely the odd classes (1,3mod8) modulo (2^r). Indeed (v_2(11^2-1)=3), and each squaring adds one to this valuation. The order is (2^{r-2}), the size of the displayed subgroup. Multiplying by 5 covers its other coset. At (rle2) use the result for (r=3). A target odd residue therefore fixes (e) and one class of (jmod P), where (P=2^{max(r,3)-2}).

Put (mu=log_2 100), which is irrational by unique factorization, and (B=6,32^e100^j 2^t). Normalize (B) into ([2^K,2^{K+1})) by choosing (t=K-lfloorlog_2(6,32^e)+jmufloor). Restrict (0le jle K/(4mu)), so (tge3K/4-O(1)), and the subtractive term (5^e11^j=o(2^K)) uniformly. As (j) runs up its prescribed class modulo (P), normalized (B/2^K) is an irrational rotation at step (Pmu) on the logarithmic circle. Pick a closed interval strictly inside ((a,b)) and use the logarithm of its interior. The orbit of an irrational circle rotation is dense; finitely many translates of any open interval cover the circle by compactness. Thus any block of some fixed number of successive admissible indices visits that inner interval, independently of the starting phase. This gives (Omega_{a,b,r}(K)) hits in ((a2^K,b2^K)) for every sufficiently large (K); (e) is fixed per residue.

There are no collisions among these (j) in a given band for fixed (e). A collision would imply (11^{|j-j'|}equiv1mod2^{min(t,t')}). The valuation of the left side minus 1 is 1 for odd difference and (2+v_2(|j-j'|)) for even difference (factor 120 and square). This is at most (O(log K)), much smaller than the exponent on the right. The upper bound (O((log X)^2)) follows from (5^e11^jle B/12) when (tge1), and hence (j,t=O(log X)).

## Quantitative version with growing conductor

Take a circle interval of width (ell>0) strictly inside ((log_2a,log_2b)). By pigeonhole, for some (mle M=lceil4/ellceil+1) the step distance (delta=|mPmu|) satisfies (0<delta<ell/2). If (n) is the nearest integer to (mPmu), the distinct positive integers (100^{mP},2^n) differ by at least 1. By the mean-value theorem for the logarithm, (deltage C,100^{-MP}) for an absolute (C>0). Starting from any phase, successive multiples of the signed step of size (delta) hit the interval in at most (lceil1/deltaceil+1) steps (one turn of the circle, with gaps at most (delta)). So a visit is guaranteed within at most (C_ell,P,100^{MP}) indices on the original (j)-axis, including the initial step into the prescribed congruence class. This is less than (K/(4mu)) for sufficiently large (K) whenever (2^rle c_{a,b}log K), for a sufficiently small constant depending only on the inner interval. The choice has (tge r) and the same uniformly small subtractive term. This is a simultaneous real-and-dyadic placement theorem with growing conductor, not merely separate projection tests.

## Why the absolute step is still missing

Separate real-density and complete finite residue coverage can coexist even in one explicit, jointly controllable **sparse** family. To hit a prescribed (N) one needs the **exact** synchronization
[
N+5^e11^j=6,32^e100^j2^t,
]
not proximity to the normalized circle coordinate and not a residue modulo (2^r) where the absolute error can still be enormous. This is an instance of the same logical gap as rational transfer rank vs. its ranks at finite places: the coordinate that disappears in one projection may be the whole obstruction in the integer problem.

A productive next step is a replacement of the (11/100) block by a family of **paired** profiles that share the multiplicative positioning coordinate but have additive degrees of freedom at the prescribed scale. Independent repetitions of today’s two profiles cannot supply that freedom; their output bound proves it. Related: [the row-fold reader experiment](ROW_FOLDS.md) shows that changing the small differential without recentering the absolute baseline can also fail.
