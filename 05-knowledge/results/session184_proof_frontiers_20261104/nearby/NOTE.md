# Infinite bounded-degree hole modules and their exact atom profiles

Session 184, 2026-11-04. Status: elementary proofs supplied; executable exact controls supplied. This is a structural extension of THM-2500, not a resolution of a major open problem or a claim of literature priority.

## Sources and target selection

The main-branch portfolio ledger explicitly requests extension of the finite cage/ray decomposition to infinite hole modules. THM-2433 and THM-2500 define relative atoms using proper multiplicative pairs, with weak and strict diagonals treated separately. The latter bounds finite-hole remainders by M^3 and M^4, where M=max H. That coordinate is unavailable for infinite H. We replace it by the additive prime-factor degree Omega.

- https://github.com/eliottcassidy2000/math/blob/main/00-navigation/PROBLEM-LEDGER.md
- https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-2433-operation-fibre-deletion-incidence-and-startup-scar.md
- https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-2500-finite-hole-relative-atoms-prime-rays-and-cover-bounds.md
- Richard J. Mathar, *Series of Reciprocal Powers of k-almost Primes*, arXiv:0803.0900v3, https://arxiv.org/abs/0803.0900v3. The symmetric-function/prime-zeta expansion used below is established machinery, not new here. All required identities are rederived below.

Only targeted current-main files were consulted; no startup helper or correction archive was read.

## 1. Definitions

Let H be any subset of {2,3,...}, finite or infinite, and S its complement in that carrier. W_H consists of retained composite n having no factorization n=ab with a,b in S and 2<=a<=b. T_H consists of retained n having a nonempty ambient strict factor fibre but no retained factorization with 2<=a<b. Thus primes are excluded from both artificial sets and prime squares are additionally excluded from T_H.

Let Omega(n) count prime factors with multiplicity. Assume K is a positive integer and Omega(h)<=K for every h in H. No bound on the sizes of h, no density condition, and no multiplicative closure are assumed.

## 2. Sharp rank theorem, including the equality cases

**Theorem A.**

1. Every n in W_H satisfies Omega(n)<=2K+1.
2. Every n in T_H satisfies Omega(n)<=2K+2.
3. If Omega(n)=2K+1, then n is in W_H if and only if every divisor d of n with 1<=Omega(d)<=K belongs to H.
4. If Omega(n)=2K+2, then n is in T_H if and only if n=p^(2K+2) for a prime p and p,p^2,...,p^K all belong to H.

**Proof.** If Omega(n)>=2K+2, split its multiset of prime factors into a block of K+1 and a complementary block of at least K+1. Both products are outside H, giving a weak retained pair. This proves (1). If Omega(n)>=2K+3, take block sizes K+1 and Omega(n)-K-1>=K+2. The products cannot be equal, since equal integers have equal Omega, proving (2).

At rank 2K+1, each proper factor pair has exactly one member of rank at most K. Its other member has rank at least K+1 and is retained. All pairs are killed exactly when every such smaller-rank divisor is a hole. The target itself is retained automatically. This proves (3).

At rank 2K+2, all splits into two blocks of K+1 have both products retained. Strict atomicity requires that every such split have equal products. If two prime entries p_i and p_j in the factor multiset differ, choose a block of K+1 positions containing i but not j. Swapping i and j changes its product; not both block products can equal sqrt(n). Hence all prime entries must be equal, so n=p^(2K+2). Its proper strict splits are p^a*p^(2K+2-a) for 1<=a<=K. The second factor is retained, so all are killed exactly when p^a is a hole for every such a. This proves (4). QED.

**Sharpness.** For H={p,...,p^K}, the targets p^(2K+1) and p^(2K+2) attain the weak and strict bounds. Infinite sharp examples follow by taking every integer of ranks 1 through K as a hole.

This is a rank bound, not a height bound. Even H={2} gives unbounded-height atoms 2q. Nor may a bound on the number of DISTINCT prime factors replace Omega: H={p,p^2,...} has unbounded Omega although every hole has only one distinct prime factor. Indeed, for any other prime q and every j>=1, the retained target p^j q is a weak and strict artificial atom: in every proper factor pair, the factor not containing q is a positive power of p and hence a hole. These targets have unbounded Omega.

## 3. A complete infinite family

Let H_K={n:1<=Omega(n)<=K}. Then

W_(H_K)={n:K+1<=Omega(n)<=2K+1},

T_(H_K)=(W_(H_K) minus {p^2:p prime}) disjoint-union {p^(2K+2):p prime}.

The subtraction matters only at K=1, but is essential: prime squares have empty ambient strict fibre and are not artificial strict atoms.

Indeed any retained factor has rank at least K+1, so all retained targets below rank 2K+2 are weak atoms. After excluding ambient prime squares, all those weak atoms are also artificial strict atoms. At rank 2K+2 Theorem A classifies the additional strict atoms; above it there are none.

For Re(s)>1 put P(s)=sum_p p^(-s) and P_j(s)=sum_(Omega(n)=j)n^(-s), with P_0=1. Unique factorization gives the coefficientwise identity

sum_(j>=0) P_j(s)t^j = product_p (1-t p^(-s))^(-1)
 = exp(sum_(r>=1) P(rs)t^r/r).

Therefore

D_W(s)=sum_(j=K+1)^(2K+1)P_j(s),
D_T(s)=D_W(s)+P((2K+2)s)-[K=1]P(2s).

The elementary Euler-product expansion gives P(s)=log(1/(s-1))+O(1) as real s decreases to 1. Indeed log zeta(s)=P(s)+sum_(r>=2)P(rs)/r and the second term stays bounded, while integral comparison gives zeta(s) asymptotic to 1/(s-1). For fixed j the generating identity gives

P_j(s)=P(s)^j/j!+O(P(s)^(j-2)+1),

because every nonidentity cycle partition replaces at least two singleton cycles by bounded P(rs), r>=2. Consequently, writing L=log(1/(s-1)),

D_W(s) ~ D_T(s) ~ L^(2K+1)/(2K+1)!.

Thus infinite hole modules can have arbitrarily high odd logarithmic singularity order. They do NOT retain the finite-hole form "finite cage polynomial times a prime tail plus finite polynomial": that form has at most a first-order logarithm.

## 4. Arbitrary prime holes: exact support, no collision multiplicity

Now let H be any set of primes, and put Q(s)=sum_(p in H)p^(-s), Q_r=Q(rs). Let h_3(Q)=(Q(s)^3+3Q(s)Q_2+2Q_3)/6.

**Theorem B.** W_H is the disjoint union of:

- semiprimes pq, p<=q, with at least one of p,q in H;
- triple prime products pqr, p<=q<=r, all three primes in H.

T_H is obtained by deleting {p^2:p in H} from W_H and adding {p^4:p in H}.

**Proof.** Theorem A with K=1 bounds weak rank by three and strict rank by four. At rank two the only pair is killed exactly when a prime factor is a hole. At rank three every proper pair is prime times semiprime, and the semiprime is retained; every distinct prime divisor must therefore be a hole. Strict and weak coincide at odd rank. The rank-four strict equality case is exactly p^4 with p in H. Ambient prime squares are excluded from T_H. QED.

The exact Dirichlet series are

D_W=QP - Q^2/2 + Q_2/2 + (Q^3+3QQ_2+2Q_3)/6,
D_T=D_W-Q_2+Q_4.

To count semiprimes without collisions, separate a mixed pair with one prime in H and one outside H (Q(P-Q)) from an unordered pair entirely in H ((Q^2+Q_2)/2). The triple term follows from the cycle-index identity, or direct classification into three distinct primes, one repeated prime, and a triple repetition. Products across the rank-two and rank-three families cannot collide because Omega differs.

If H is nonempty, both reciprocal sums diverge and the absolute convergence abscissa is 1: fix p in H and use the distinct targets pq for all primes q other than p. If H is empty both sets are empty. This applies even to infinite prime sets with convergent reciprocal sum.

## 5. An exact conditional singularity transition

Assume H is nonempty and, as s decreases to 1,

Q(s) ~ c L^alpha,   c>0,   0<=alpha<=1.

This is a hypothesis on a weighted prime profile, NOT an asserted density theorem for arbitrary H. Theorem B implies, for both D_W and D_T:

- 0<=alpha<1/2: D(s) ~ c L^(1+alpha);
- alpha=1/2: D(s) ~ (c+c^3/6)L^(3/2);
- 1/2<alpha<=1: D(s) ~ (c^3/6)L^(3alpha).

Proof: QP has order L^(1+alpha), Q^3/6 has order L^(3alpha), Q^2 has order L^(2alpha), and Q_2,Q_3,Q_4 are bounded. The QQ_2 term is O(L^alpha); all listed remainders have smaller exponent than the dominant term. At alpha=0 Q tends to c and QP still supplies the first-order logarithm. At alpha=1 the cubic dominates even the quadratic terms. The strict correction is bounded. QED.

This gives a precise obstruction to transporting a finite-cage first-order logarithmic heuristic to arbitrary infinite prime holes: the crossover is Q approximately sqrt(P), not merely "H is infinite". No assertion of existence of every profile exponent is needed or made.

## 6. Verification and boundaries

Run `python3 verify.py` and `python3 -O verify.py` from this directory. All checks use explicit exceptions. The independent direct-fibre path uses trial divisors; classifiers use prime-factor multisets. The suite checks every subset of the first six primes; rank-cutoff holes and sharp endpoints; arbitrary bounded-rank small holes; exact integer Dirichlet coefficients including the diagonal corrections; and finite-hole, full-prime-hole, and sparse-hole hostile controls.

No twin-prime, Collatz, LRC, Jacobian, or analytic-continuation consequence is claimed. The prime-zeta formulas themselves are standard symmetric-function machinery. The contribution here is the sharp infinite-hole rank/equality classification, its full rank-ideal specialization, and the exact support-level prime-hole transition, linked to the repository's stated frontier.
