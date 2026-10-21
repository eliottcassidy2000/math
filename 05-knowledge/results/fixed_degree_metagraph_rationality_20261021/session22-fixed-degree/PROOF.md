# Every fixed fundamental black degree has a rational counting series

Status: PROVED from the explicitly inherited inputs; independent algebra/type review passed. This is a corollary of existing effective finite-atom and center laws, not a complete wiring classification. New Session22 synthesis; no historical lost computation is claimed restored.

## Conventions and inherited inputs

For an ordinary tournament isomorphism type T let h=H/|Aut(T)| (odd), and b be its number of blue Hamilton-word orbits. The fundamental converse-merged paired-black graph Q has degree h at a non-self-converse pair, and (h-b)/2 at a self-converse type; loops count twice. Literal doubled-edge metagraph degrees are NOT used.

Sources: effective_black_defect_order_bound_20261020/PROOF.txt proves every strong tournament of black deficit Delta=h-b<=d has order<16(d+2)^2(d+1). In particular there are finitely many strong atoms of any fixed h, and finitely many strong self-converse centers of any fixed positive deficit. The rigid strong atoms are K1,C3 (THM-2477 and current all-order supersession), and metagraph_classification_20261021a proves the only strong h3 atom is cyclic R5. Ordered SCC decomposition is unique, h multiplies, and converse reverses/opposes the SCC word. THM-2454 and deficiency_three_20261020/centers.md prove for T=U -> C -> U^op:

h(T)=x^2 a, b(T)=x b, Delta(T)=x(xa-b),

where x=h(U), a=h(C), b=b(C), and empty C is assigned (a,b)=(1,1). The existing general center-divisor formula is inherited, not newly discovered here. The conclusion below combines it with finite strong alphabets. It does not require the new stronger Session22 center cutoff.

Let p(z)=1-z-z^3, A(z)=1/p(z). Define:
P_j(z): strong ordinary atoms with h=j, by vertex order;
C_(a,b)(z): strong self-converse centers of signature(a,b);
C_a^+(z): all strong self-converse centers with h=a, plus the empty center when a=1;
F_j(z): all ordinary types with h=j (including empty type in F1).

All P_j and C_a^+ are finite polynomials, by the inherited bound and Delta<=h. P1=z+z^3, P3=z^5, C1^+=1+z+z^3, C3^+=z^5. These definitions count actual isomorphism types, never every algebraically admissible integer signature.

## 1. Ordered multiplicative factorization

For positive odd j, let Omega(j) count prime factors with multiplicity, Omega(1)=0. Unique SCC decomposition gives the finite expression

F_j(z) = sum over r=0..Omega(j) of A(z)^(r+1) * sum_{j1...jr=j, ji>1 odd} product_i P_ji(z).

The empty tuple contributes only for j=1. Explanation: delete rigid SCCs, leaving an ordered list of r nonrigid atoms whose h-values multiply to j; arbitrary K1/C3 words fill the r+1 gaps. Every list is counted exactly once. Each ji has at least one prime factor, hence r<=Omega(j). Therefore F_j has denominator dividing p^(Omega(j)+1). Set F_j=0 for positive even j.

This is an effective finiteness/rationality theorem, not a practical computation of the unknown P_j coefficients. It does not assert every odd j is realizable.

## 2. Self-converse counts at fixed h and fixed deficit

Let J_j count self-converse types with h=j. Unique SCC reflection gives

J_j(z)=sum_{x odd, x^2 divides j} F_x(z^2) C_(j/x^2)^+(z).

This is ordinary class counting: the first-half type determines the opposite second half; no factor1/2 occurs. Strong centers are mandatory, except the single empty possibility.

For d>0 even, let S_d count self-converse types of deficit d, and let R_d count strong self-converse centers of deficit d. For odd x>1 dividing d define E_(x,d)(z) to count empty/strong centers with x(xa-b)=d. The inherited inequalities b<=a and a>=1 yield x(x-1)<=d and a<=d/[x(x-1)]. Thus every E_(x,d) is polynomial. The exact inherited divisor identity is

S_d(z)=A(z^2)R_d(z)+sum_{x>1 odd, x|d, x(x-1)<=d} F_x(z^2)E_(x,d)(z).

This proves S_d rational with only p(z^2) as a possible denominator factor. This implication is not a computation of R_d or of E_(x,d).

## 3. All fixed-degree metagraph node counts

Let B_k(z) count nodes of degree k in Q. For k>=1,

B_k(z)=S_(2k)(z) + 1[k odd]*(F_k(z)-J_k(z))/2.

The terms are disjoint: self-converse classes versus two-element converse pairs. Converse acts freely on the non-self-converse part, so coefficients are integers.

Define L_k=max{Omega(x): x odd, x|2k, x(x-1)<=2k}; x=1 is included. For odd k, every x with x^2|k also occurs in this set. Hence

Denominator(B_k) divides p(z)^(Omega(k)+1) p(z^2)^(L_k+1) for odd k;
Denominator(B_k) divides p(z^2)^(L_k+1) for even k.

This proves rationality for EVERY fixed positive degree, not merely the trivalent layer. For degree0 the established pure-blue alphabet already gives
B0(z)=(1+z+z^3+z^5)/p(z^2).

The displayed denominator is an upper bound, not generally minimal; cancellations or empty strata may occur. Every fixed-degree sequence therefore satisfies an eventual integer constant-coefficient recurrence. Unknown finite center/atom coefficients affect the numerator and onset, so this is not an algorithm that outputs all coefficients without the finite classification inputs.

## 4. Two exponential scales and a fully explicit infinite family

Let rho in(0,1) solve p(rho)=0. It is the unique smallest-modulus root of p: if p(z)=0 and |z|<=rho, triangle equality in 1=z+z^3 forces |z|=rho and both terms positive real, hence z=rho. The roots of p(z^2) of smallest modulus are +/-sqrt(rho).

For every even k>=2,
[z^n]B_k = O(n^L_k rho^(-n/2)),
and the same exponential bound without polynomial factor holds for k0. Constants depend on the finite unknown polynomials. Zero sequences are allowed.

For k=3^r, r>=0, the term with the largest possible number r of nonrigid SCCs is uniquely all R5 atoms, since P3=z^5. Thus
F_(3^r)(z)=z^(5r) A(z)^(r+1) + terms with pole order at rho at most r.
The self-converse J and S terms have all poles of modulus at least sqrt(rho), strictly greater than rho, so they do not cancel the leading rho pole in B_(3^r). Consequently

[z^n]B_(3^r)(z) ~ rho^(4r-1) / [2 r! (1+3rho^2)^(r+1)] * n^r rho^(-n).

The constant is positive and explicit even though all higher center coefficients remain unknown. For r=0 this gives the inherited leaf leading scale; r=1 recovers the trivalent leading coefficient rho^3/[2(1+3rho^2)^2]. This is a generalization from the finite-atom mechanism, not a new proof of the h3 classification.

More generally, if the presentation count h=k is realizable by a tournament (equivalently F_k is not zero), then B_k has a dominant pole at rho: choose the maximum number r of realizable nonrigid factors in a strong-atom factorization of k; all corresponding polynomial weights are positive at rho, so no cancellation at the leading pole can occur. This r may be less than Omega(k). For k=1 the empty nonrigid factorization has r=0. Thus odd degrees k with F_k nonzero have rho^(-n) growth times a polynomial; this does not assert that existence of a self-converse degree-k node implies F_k is nonzero. If F_k=0, only the self-converse term remains. even fixed degrees have at most rho^(-n/2) growth times a polynomial. This concerns fixed-degree counts, not degree distributions at growing k.

## 5. Bounded-degree nodes are asymptotically negligible

For each fixed k, rationality with these pole sets implies [z^n]B_k/n! tends to zero. For each fixed K, the sum of counts of degrees0 throughK is at most C_K n^D_K rho^(-n), for constants depending on the finite unknown alphabets. Let M_n be the total number of converse-merged tournament isomorphism types at order n. Each merged type contains at most2n! labelled tournaments, so M_n>=2^(n(n-1)/2)/(2n!). It follows that

(sum_{k=0}^K [z^n]B_k)/M_n -> 0.

Thus a uniformly chosen merged node has degree tending to infinity in probability, in the precise sense that the probability of degree at most any fixedK tends to zero. This uses the fundamental graph degree convention. It is not a limiting distribution, rate-uniform statement inK, or connectivity theorem. The distinction matters: almost all nodes escaping bounded degree need not form one connected component.

## 6. Scope, hostile tests, and connection contract

The source is ordered SCC words and actual strong-type alphabets. The map records multiplicative h and the unique converse center; it preserves class counts and degree. It forgets literal Hamilton words, endpoints, nonpath-complement targets, edge multiplicities, bridges and connectivity. Those need explicit port sidecars. Rational node counts do not realize or classify wiring.

Hostile: a nonrigid pure-blue half R5 paired with its opposite is NOT pure-blue; h9,b3,degree3. This is exactly x=3 in S6. Quasi-polynomial growth is false even at degree0; rationality allows exponential poles. Shared denominator roots do not identify electrical currents or Collatz arrivals. The useful connection is finite alphabet plus unique composition, not similar numerical sequences.

Independent tests must check divisor bounds, exact word factorization, ordinary versus converse counts, known B1/B3 formulas, synthetic-alphabet multiplicities, and rational denominators. Synthetic tests are algebraic controls, not realizability of invented tournament atoms. No large-order center census is claimed.
