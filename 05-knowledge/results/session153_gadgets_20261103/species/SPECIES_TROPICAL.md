# Strong-component generating series and tropical order

Status: exact reformulation and small independent checks. No novelty or literature-priority claim. The fixed-H formula and its eventual-polynomial consequences were already proved in `session24-multiplicities/PROOF.md`, equations (1), (4), (5); the min-plus interpretation is an immediate corollary. See also `math/05-knowledge/results/camion_busch_gaps_polyhedra_collatz_20261002.md`, §4, for the finite strong-core bound and the existing reduction. In particular, the series below does not solve the spectrum problem.

Let H(T) be the number of directed Hamiltonian paths and let the empty tournament have order 0 and H=1. Isomorphism classes are ordinary, not converse-merged or labelled. For d>1, define A_d(x)=sum x^|C| over nontrivial strongly connected tournament classes C with H(C)=d; write A_d=0 for even d. Camion gives |C|≤H(C), so these are finite polynomials. Let F_m(x)=sum over all T with H(T)=m of x^|T|. There is exactly one transitive class in each order, including 0, and F_1(x)=1/(1-x).

Work formally with arithmetic coefficients and Dirichlet multiplication. Define F(x,s)=sum_{m≥1} F_m(x)m^(−s), and A(x,s)=sum_{d>1} A_d(x)d^(−s). Coefficients of every fixed x^n involve finitely many classes; no analytic convergence is asserted. Then

  F(x,s) = 1 / (1 − x − A(x,s)).

Proof: Every tournament has a unique linearly ordered sequence of strongly connected components. A component is either a singleton, of weight x, or a nontrivial strong class, of weight x^|C| H(C)^(−s). All edges point forward between components, forcing each Hamiltonian path to traverse components in this order and making H multiplicative. Isomorphisms preserve the canonical linear component order, so no factorial or orbit correction is needed. Summing all component words proves the formula. The same identity with x^n/n! and labelled class counts requires the appropriate labelled product; it is not being claimed as this ordinary formula.

Extracting Dirichlet coefficients recovers precisely the prior fixed-H expression:

  F_m(x) = sum_{d_1⋯d_k=m, d_i>1} [∏ A_{d_i}(x)] / (1−x)^(k+1),

where factorizations are ordered and the empty factorization occurs only for m=1. Singletons form independent runs in the k+1 gaps. For each m this is finite because d_i≥3.

Define σ(d) as the least order of a nontrivial strong tournament with H=d, or infinity if there is none, and τ(m) as the least order of any tournament with H=m, or infinity if absent, with τ(1)=0. Taking the lowest power of x in the positive-coefficient series gives

  τ(m)= min_{d_1⋯d_k=m, d_i>1} Σ σ(d_i)
      = min_{d|m, d>1} [σ(d)+τ(m/d)]  (m>1).

There is no cancellation. For each m, the full order support is exactly {n≥τ(m)} if τ(m)<infinity, because singleton strong components may be prepended; otherwise it is empty. Consequently H=m at order n iff this finite divisor-lattice min-plus calculation gives τ(m)≤n. This is a precise finite criterion once the strong divisors' σ are known, not a shortcut for determining them. A finite acyclic transfer matrix on divisors d|m, with edge d→de weighted x^σ(e), gives the same min-plus shortest-path problem. This is a tropical encoding, not a separate valuation obstruction or a claim that the full tournament problem is context-free.

Example: the only nontrivial strong tournament with H=3 is the triangle: Camion bounds its order by 3. Thus A_3=x^3 and F_3=x^3/(1−x)^2, so exactly n−2 classes have H=3 for n≥3, and none for n<3. Similarly A_5=x^4 by the Camion bound plus checking strong classes through order 5, so τ(5)=4. Independent direct enumeration below shows two strong order-5 classes with H=9; no order <5 realizes H=9. Hence σ(9)=τ(9)=5, while treating 9 only as 3×3 costs 6 vertices. This illustrates why factor-by-factor minimization must include composite strong atoms.

Small independent check: `python check.py` enumerates labelled orientations through n=5, quotients via the lexicographically least encoding over vertex permutations, computes H by subset dynamic programming, and checks strong connectivity via transitive closure. It checks the H=3 formula and prints class counts for H=3,5,9; expected rows n=1..5 are (0,0,0),(0,0,0),(1,0,0),(2,1,0),(3,2,2), with strong-count row at n=5 (0,0,2). This is an illustrative finite check, not a proof of the general identities.
