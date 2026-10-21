# Renewal geometry of fixed-presentation tournament strata

Status: PROVED from the stated inherited finite-alphabet results; independent proof review passed. No claim about a uniformly random unrestricted tournament.

## Objects and inherited dependency
For an ordinary tournament isomorphism class let h(T)=H(T)/|Aut(T)|. Use the unique ordered strong-component decomposition and the finite actual strong alphabets proved in `fixed_degree_metagraph_rationality_20261021/session22-fixed-degree/PROOF.md`. Strong h=1 types are exactly K1 and C3. Fix a realizable odd h. Delete all K1/C3 components; call the resulting ordered list the skeleton. There are finitely many possible skeletons, each of bounded total order. Let r be the maximum skeleton length. For h=1 the empty skeleton has r=0. Write S for a skeleton, |S| for its vertex order, and let W be the finite set of length-r skeletons (individual strong isomorphism types, not just their h-values).

A class with skeleton S of length l is uniquely a sequence of l+1 words over K1,C3 separated by its l strong atoms. For a maximal skeleton, l=r. Let M(T) be the number of neutral C3 components. Let G_i be the number of vertices in neutral gap i, i=0,...,r, conditional on maximal skeleton length. Define rho as the positive solution rho+rho^3=1, a=rho^2, mu=a/(1+3a), and sigma²=a(1+a)/(1+3a)^3.

## Theorem A: maximal skeleton concentration and its Boltzmann weights
Choose an ordinary isomorphism type uniformly from order-n types with h(T)=h. Such types exist for every sufficiently large n by adding K1s to any witness. The probability its skeleton has length r is 1-O(1/n) (exactly1 if no shorter skeleton exists). Moreover, for each S in W,

P(skeleton=S) -> rho^|S| / sum_(S' in W) rho^|S'|.

Proof: write A(z)=1/(1-z-z³), c=1/[rho(1+3rho²)]. The generating function for one specific skeleton S of length l is z^|S| A(z)^(l+1). The root rho is simple and uniquely smallest in modulus; its principal part is c/(1-z/rho). Therefore its n-th coefficient is rho^|S| c^(l+1) n^l rho^(-n)/l! (1+O(1/n)). Finite summation gives the claims. No actual atom multiplicities are replaced by invented integer signatures.

## Theorem B: uniform-simplex neutral gap geometry
Conditional on any fixed maximal skeleton S, as n tends to infinity,

(G_0/n,...,G_r/n) converges in distribution to Dirichlet(1,...,1), the uniform probability measure on the simplex of r+1 nonnegative coordinates with sum1.

For r=0 the law is the point mass1. For r>=1: put N=n-|S|. If a_m=[z^m]A(z), then the probability of a particular gap-size vector g with sumN is product a_gi / [z^N]A^(r+1). Uniformly as m grows, a_m=c rho^(-m)(1+O(theta^m)) for some theta<1. Also a_m<=C rho^(-m) for all m. For fixed epsilon>0, all compositions with g_i>=epsilon N have asymptotically equal weights, uniformly. The total probability of any g_i<epsilon N is O(epsilon)+O(1/N), by the coefficient asymptotic and counting such compositions using the global bound. Riemann sums on the interior simplex followed by epsilon down to0 prove the assertion. Finite mixing over S preserves the same limit. Thus this limit is independent of the limiting skeleton type.

Equivalently the normalized positions of the r exceptional strong atoms converge to the order statistics of r independent uniform points in [0,1]; their bounded lengths disappear under normalization. This is not a statement that finite-n positions are independent.

## Theorem C: neutral triangles and automorphism entropy
With the same uniform fixed-h model,

(M(T)-mu*n)/sqrt(sigma²*n) converges in distribution to the standard normal law,
E M(T)=mu*n+O(1), Var M(T)=sigma²*n+O(1).

Numerically mu≈0.1942540040 and sigma²≈0.04956151754. These numerical values are not used in the proof.

Self-contained generating-function argument: mark each neutral C3 by u. The full series is a finite sum of z^|S| (1-z-u z³)^(-length(S)-1). For complex u near1, the implicit-function theorem gives a simple dominant root R(u), R(1)=rho, with the other roots separated uniformly. The top-length terms have nonzero analytic leading amplitude K(u) near1 (positive at1). Partial fractions give uniformly

[z^n]F_h(z,u)=K(u)*n^r*R(u)^(-n)*(1+O(1/n)).

The error can be differentiated on a smaller complex neighborhood by Cauchy's formula. Thus the logarithm of the moment generating function is

n*log(rho/R(exp(t))) + log(K(exp(t))/K(1)) + O(1/n).

Let y(t)=log R(exp(t)), A(t)=exp(t)R(exp(t))². Implicit differentiation of R+exp(t)R³=1 gives y'=-A/(1+3A), A'=A(1+A)/(1+3A). Hence -y'(0)=mu and -y''(0)=sigma²>0. The first two derivative identities yield mean and variance. Set t=i s/sqrt(sigma² n), subtract the centering term, and use Taylor's formula; the characteristic function tends to exp(-s²/2), proving the distributional limit. This argument uses only a finite rational function with a uniformly isolated simple moving root, not an assumption about independent components.

Since Aut(T) is the product of the component automorphism groups and each C3 contributes3, log |Aut(T)| = M(T) log3 + O_h(1). Consequently log |Aut(T)|/n tends in probability to mu log3, and its centered/scaled Gaussian law follows with variance constant sigma²(log3)². Nonrigid atoms are a finite fixed-h alphabet, so their contribution is uniformly bounded. The same conclusion holds for v_3(|Aut(T)|), with bounded additive difference from M(T). Since H(T)=h*|Aut(T)| at fixed h, exactly the same normalized limit and Gaussian fluctuation hold for log H(T). In particular H(T)^(1/n) tends in probability to 3^mu, despite the fixed number h of Hamilton-path automorphism orbits. The total SCC count L(T) equals n-2M(T)+O_h(1), so L(T)/n tends to 1-2mu and its Gaussian variance constant is 4sigma². The proportion of vertices in neutral triangles tends to3mu; those in neutral singletons tends to1-3mu. Nonrigid skeleton vertices contribute only O_h(1). These are orbit-versus-raw-count distinctions, not equalities of the two measures.

## Corollary: actual fixed-degree merged-node geometry
For any odd k with F_k nonzero, a uniformly chosen order-n node of fundamental black degree k is, with probability tending exponentially to1 up to a polynomial factor, a non-self-converse pair arising from an ordinary h=k type. This follows from Session22's p(z) versus p(z²) scales. For every converse-invariant statistic above (M, logAut, maximal skeleton length) the same limits hold for these nodes: each non-self-converse pair has exactly two ordinary types. To speak about the ordered gap vector, choose either representative with probability1/2. Reversal reverses the gap vector; the uniform-simplex law is reversal invariant.

For k=3^r the inherited P3=z5 makes the maximal skeleton uniquely r copies of R5. Thus the limiting positions have the explicit uniform-order-statistic law even though the finite strong-center polynomials are not known. r=0 includes degree1; r=1 gives a uniform single R5 position in the overwhelmingly dominant trivalent nodes.

## Boundaries and connection contract
Source: fixed-h actual SCC words. Target: finite weighted skeleton plus neutral renewal gaps. Map deletes neutral components but retains their ordered words and sizes. Preserved: exact class counts, neutral counts, component order, product of automorphism orders. Lost if only the skeleton is retained: gap locations, literal Hamilton paths, ports, metagraph wires. The gap words restore class identity, but none of these scalar counts recovers black-edge adjacency.

The conclusion is NOT a model of uniformly random labelled tournaments or unrestricted merged nodes, whose bounded-degree mass tends to0. Uniformity is over isomorphism types, not labels; automorphism-weighting would change the triangle fugacity and the constants. Maximal skeleton length is not automatically Omega(h): absent actual strong factors can lower it. Geometric order statistics do not imply graph connectivity. All claims are asymptotic at fixed h/k, not uniform in growing h/k.

Prior search: searched exact phrases Dirichlet gaps, uniform simplex, neutral Gaussian/central limit and quasi-powers in current theorem/results/reflection corpus; no matching tournament result located. This is not an external novelty claim. The generating-function proof is a direct consequence of the inherited finite-alphabet decomposition, not a new classification of the strong atoms.

## Theorem D: variance doubling on self-converse deficit slices
Fix an even deficit d>=0 for which the self-converse slice is nonempty. Restrict n to either parity on which the slice is nonempty for arbitrarily large n. For a uniformly chosen self-converse class of that order and deficit, the total number M of neutral C3 components satisfies

E M=mu*n+O_d(1), Var M=2*sigma²*n+O_d(1),
(M-mu*n)/sqrt(2*sigma²*n) -> Normal(0,1).

For d>0, the inherited center-divisor formula is a finite sum of F_x(z²) E_(x,d)(z), including x=1 and E_(1,d)=R_d. Each actual center C has fixed bounded size c and neutral-triangle contribution e_C equal1 if C=C3 and0 otherwise; its half U has h=x and order (n-c)/2. Unique SCC decomposition makes these disjoint families. For a fixed admissible parity, only centers with matching c occur. Every nonempty component family occurs at all sufficiently large n of that parity, since neutral K1 can be appended to U. Within one family,

M(T)=2 M(U)+e_C.

Theorem C at half-size gives mean mu*n+O(1), variance 2sigma²*n+O(1), and the displayed CLT. There are only finitely many families. Their means differ by bounded constants, so the total-variance identity adds only O(1), even when their mixing weights depend on n. A finite mixture of the same limiting normalized Gaussian is again that Gaussian. No assumption of independent mirror halves is made: they are identical after converse, which is exactly why the variance doubles. At d=0 use the inherited pure-blue center alphabet {empty,K1,C3,R5} and x=1, so the same argument applies.

Consequently every nonempty even fundamental degree k has this variance-doubling law on each supported parity; degree0 does too. Odd degree k with F_k=0 also lies wholly in this self-converse regime. If F_k is nonzero, Theorem C instead governs its overwhelmingly non-self-converse merged nodes. These are fixed-degree laws; neither mixture may be extrapolated to growing degree.

Hostile comparison: treating reflected halves as independent would incorrectly predict sigma²*n rather than2sigma²*n. The exact identity M(T)=2M(U)+e_C is the preserved dependence sidecar.

## External methodological context
The moving-pole/characteristic-function method is classical analytic combinatorics, not a novel probability theorem. See Flajolet and Sedgewick, Analytic Combinatorics, Chapters V and IX (author-hosted overview: https://sedgewick.io/books/analytic-combinatorics/). Our argument above supplies the particular root calculation and its proof directly; no external theorem is needed to import a tournament statement. The attempted full-book fetch was too large, so no claim is made here to have freshly audited the complete book text.

## Corollary E: a classification-free asymptotically uniform sampler
Fix r>=0. Consider the actual tournament family of ordered SCC words over {K1,C3,R5} containing exactly r R5 components. At total order n, put N=n-5r. For each 0<=m<=floor(N/3), let a=N-3m. The exact number of words with m triangles and a singletons is

w_m=(a+m+r)!/(a! m! r!).

For N>=0 an exact uniform sampler of this subfamily first draws m with probability w_m/sum_j w_j, then draws a uniformly random multiset word with a singleton symbols, m triangle symbols and r R5 symbols. Sequentially selecting the next symbol in proportion to its remaining multiplicity gives exactly that uniform multiset law. Every word determines a unique actual tournament isomorphism type; this is not labelled sampling and no isomorphism rejection is needed.

To sample non-self-converse merged nodes of degree 3^r, reject palindromic component words and identify each remaining word with its reversal. All three atoms are self-converse, so the rejection criterion is exact. Every surviving pair has two words and therefore the resulting distribution is exactly uniform on the resulting subfamily of merged nodes. This rejection sampler is asserted only for n for which a nonpalindromic word exists; such words exist for every sufficiently large n. Rejection probability tends exponentially to0 up to a polynomial factor, by the p versus p(z²) counting scales.

By Theorem A and the inherited fact that the unique maximal h=3^r skeleton is R5 repeated r times, this sampler's law has total-variation distance O_r(1/n) from the uniform law on all order-n fundamental degree-3^r nodes. The extra self-converse degree-3^r nodes have exponentially negligible proportion. The total-variation identity used here is exact: a uniform law on a subset A of a finite set B differs from the uniform law on B by 1-|A|/|B|. For r=0,1 no shorter realizable nonrigid skeleton contributes, so the error is exponentially small up to polynomial factors. No explicit finite-n error constant or effective onset is claimed without the remaining atom/center inventories.

This is a constructive transfer from generating-function asymptotics to sampling despite incomplete classification. It does not sample unrestricted tournaments, metagraph edges, or any requested degree growing with n.
