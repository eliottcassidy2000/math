# Fibonacci color addition: exact phase convolution, square targets and branch extremes

Status: proved for the Session41 canonical three-unit word and its classified occurrence-dependent branch family. The arithmetic carrier is positive integers with ordinary addition. No prime-weighted estimate, Goldbach theorem, shared-global-dictionary classification or tournament correspondence is claimed.

## 1. The colors are three rotation intervals

Put alpha=(3−sqrt5)/2=phi^(-2), b=1−2alpha and beta=1−alpha. Let

I_B=[0,b), I_C=[b,beta), I_A=[beta,1).

Then the canonical color of n>=1 is c exactly when {n alpha} belongs to I_c. Indeed Session41 gives

a(H)=floor((H+1)alpha), c(H)=floor((H+2)alpha).

The increment a(n)−a(n−1) is1 exactly when {n alpha}>=1−alpha. The increment c(n)−c(n−1) is1 exactly when {n alpha} lies in [1−2alpha,1−alpha). The remaining interval is B. No positive integer orbit point meets a boundary, by irrationality. Thus densities are lambda_A=lambda_C=alpha and lambda_B=b=phi^(-3).

This is a deterministic phase coding, not three interchangeable states. Complete Fibonacci addresses carry more information than these three intervals.

## 2. An exact phase function predicts every colored additive pair count

For N>=2 define the ORDERED count

R_cd(N)=#{1<=a<N: color(a)=c, color(N−a)=d}.

Define the circle convolution

g_cd(theta)=integral_0^1 1_(I_c)(t) 1_(I_d)(theta−t mod1) dt.

Then, uniformly in N and in c,d,

|R_cd(N)−(N−1)g_cd({N alpha})| <=16(1+floor(log_2(N−1))).       (1)

In particular R_cd(N)=(N−1)g_cd({N alpha})+O(log N). This is a parameter-free first-order formula with an explicit safe error bound. It is not merely the product of the two color densities at a fixed target: the target phase survives.

Proof. Since {(N−a)alpha}={Nalpha−{aalpha}}, the admissible phases of a form the intersection of a fixed interval with a reflected circle interval. It has at most two interval components, with total length g_cd(theta).

For completeness, an elementary discrepancy bound suffices. The convergent denominators of alpha are1,2,3,5,8,...; for a convergent p/q, |alpha−p/q|<1/q^2. In any q consecutive orbit positions, corresponding points of the shifted rational q-grid move by circular distance less than1/q. A membership change for an interval can occur only within1/q of either endpoint. Each of these two arcs contains at most three grid points. The rational-grid interval count differs from q times length by at most2. Hence the irrational block discrepancy is at most8, uniformly in the starting phase and interval.

Greedily subtract the largest available convergent denominator from a block length L. Consecutive denominators have ratio at most2, so each positive remainder is less than half the previous one. At most1+floor(log_2 L) blocks result. Applying the block bound to each of at most two interval components proves(1). Small q and a wrapped interval are harmless under the stated constant; no unproved optimal discrepancy constant is used.

### The requested sine/cosine connection is literal here

If I_c has length lambda_c and midpoint m_c, its nonzero Fourier coefficient is

e^(−2 pi i r m_c) sin(pi r lambda_c)/(pi r).

Consequently the exact continuous phase profile is the absolutely convergent series

g_cd(theta)=lambda_c lambda_d
 +2 sum_(r>=1) [sin(pi r lambda_c)sin(pi r lambda_d)/(pi^2 r^2)]
                 cos(2 pi r(theta−m_c−m_d)).                 (2)

The sine terms encode interval widths; cosine encodes the sum of their phase centers relative to the target. Formula(2) concerns a concrete additive convolution, not a free analogy with tournament spectra.

## 3. The forced-A tent and square targets

For A,A, direct interval intersection gives

ell(theta)=g_AA(theta)
 =max(0,min(1,theta+alpha)−max(beta,theta))
 =0                         for 0<=theta<=b,
 =theta−b                   for b<=theta<=beta,
 =1−theta                   for beta<=theta<1.             (3)

The support has length2alpha, its peak is alpha, and its integral is alpha^2. Outside this support R_AA(N)=0 EXACTLY, not just asymptotically.

Let U_A(N) count unordered pairs of DISTINCT positive A integers summing to N. Then

U_A(N)=[R_AA(N)−1_(N even and color(N/2)=A)]/2.              (4)

As K tends to infinity,

(1/K) sum_(k=1)^K R_cd(k^2)/k^2 -> lambda_c lambda_d,
(1/K) sum_(k=1)^K U_A(k^2)/k^2 -> alpha^2/2,               (5)

and the natural density among square roots of representable square targets is

lim_(K->infinity) #{1<=k<=K: U_A(k^2)>0}/K =2alpha.         (6)

Define the k=1 representation count as0; this initial convention does not affect a limit.

Proof of the analytic input. For every nonzero integer h and fixed shift r>=1, differencing the phase h alpha k^2 gives the linear phase 2hr alpha k+hr^2 alpha. Its normalized exponential sum tends to0 by the geometric-series formula, since 2hr alpha is irrational. The van der Corput inequality with a fixed number H of shifts therefore gives limsup of the squared normalized quadratic sum at most a constant times1/H. Let H tend to infinity. Every nonconstant Fourier coefficient vanishes, so trigonometric approximation gives equidistribution of {k^2 alpha}. This is the standard Weyl argument; an author-hosted reference is [Tao's notes, Proposition1, Lemma3 and Corollary6](https://terrytao.wordpress.com/2010/03/28/254b-notes-1-equidistribution-of-polynomial-sequences-in-torii/).

Now apply equidistribution to the continuous piecewise-linear g_cd. The error in(1), after division by k^2, tends to0 and has vanishing Cesaro mean. The integral of a convolution is the product of its interval lengths, proving(5), with(4) removing a negligible diagonal.

For(6), every phase outside(b,1) has no A pair. On any closed subinterval of(b,1), ell is bounded below by a positive constant, so(1) eventually gives R_AA(k^2)>1 and hence U_A(k^2)>0. Equidistribution supplies the lower density from such trimmed intervals and the upper density from the whole support. Let their endpoints approach b and1. This proves the exact density, without assigning a uniform finite arrival horizon to phases arbitrarily close to the boundary.

## 4. Every legal branch retains the forced square edges

Session41's occurrence-dependent family has A fixed, B allowed to become A or B, and, in the retained-zero family, C allowed to become B or C. Therefore every legal branch's A set contains the canonical A set. Equations(5) and(6) supply a branch-independent infinite family of square-sum monochromatic edges. The smallest forced-A edge is2+7=9.

Thus no member of this branch family is a PROPER coloring of the entire square-sum graph, even though the colors can be meaningful allegiance labels or additive statistics. The canonical branch attains exactly the density(6); other legal branches can add further A representations. Uniform A-pair averages are lower bounded by the canonical limits, not necessarily equal to them.

## 5. The opposite end: a legal A class that represents every integer from7 onward

Choose B->A at every occurrence and leave canonical C unchanged. This is a legal fixed selector in both occurrence-dependent families, independent of primality or the target N. Its A class is A_max=N_positive minus canonical C, of density1−alpha>1/2.

Write H=N−1. Inclusion-exclusion gives the exact identity

R_Amax,Amax(N)=H−2c(H)+R_CC(N),

so R_Amax,Amax(N)>=H−2floor((H+2)alpha).                    (7)

For H>=11, the right side exceeds (1−2alpha)H−4alpha, which already exceeds1 at H=11 and increases thereafter. At H=10 the integer bound is2. Thus every N>=11 has at least two ordered representations and therefore a distinct-summand pair. The remaining N=7,8,9,10 have explicit pairs2+5,3+5,2+7,2+8. The last failure is6: among A_max integers below6, only2,3,5 occur, and the only sum6 is3+3. Also5=2+3 works.

Hence A_max is an asymptotic additive basis of order2, with exact last DISTINCT-summand failure6. In particular every square target k^2 with k>=3 has an A_max square-sum edge, whereas canonical A succeeds for a proportion2alpha of square roots. Both branches have genuine but different mathematical roles; neither needs to be discarded.

The word 'basis' here concerns unrestricted positive integer summands in an explicitly defined set. It says nothing about requiring those summands to be prime.

## 6. A branch generating polynomial at each additive target

Among unordered distinct pairs {a,N−a}, count canonical types AA, AB and BB as r_AA,r_AB,r_BB. Any pair containing a canonical C can never become AA. Let

u=b(H)−r_AB−2r_BB,

where b(H) is the NUMBER of canonical B positions through H, not the interval endpoint b. These u B positions are in C-blocked pairs or at the excluded diagonal. Over all B->A/B selectors on the prefix, the exact generating polynomial for the number of unordered distinct A representations of N is

Z_N(z)=2^u z^(r_AA) (1+z)^(r_AB) (3+z)^(r_BB).             (8)

Proof. Reflection a<->N−a partitions the prefix into independent unordered pairs and a possible central singleton. An AA pair always contributes1. An AB pair contributes0 or1 with one choice each. A BB pair contributes1 only when both positions are promoted to A, yielding three choices of weight1 and one of weight z. The remaining B choices are unconstrained factors2. The optional retained-zero C->B/C choices multiply the entire polynomial by2^c(H) without changing A membership.

In particular every integer from r_AA through r_AA+r_AB+r_BB is attained, with the exact multiplicities in(8). This is an all-size formula for one target, not independent randomness of the colors across several targets. Different N share the same vertex selectors; their joint generating function must keep those shared variables rather than multiply the separate Z_N.

## 7. Scope and reproduction

The checker uses integer-square-root formulas and exact Q(alpha) comparisons for every color, convolution endpoint, tent value and error inequality. It independently regenerates the substitution prefix, checks all pair counts through1500, and checks square targets through250^2. Decimal sample averages are reported only as diagnostics, never as proofs of a limit. It also exhausts short-prefix selectors to verify(8) and the last-failure statement for A_max.

The discrepancy proof and square equidistribution explain the infinite results. Prime constraints would require additional estimates for prime-weighted, jointly constrained sums; ordinary equidistribution and finite Goldbach checks do not provide those estimates. The larger shared-global-expansion dictionary remains a different model from the occurrence-dependent family used here.
