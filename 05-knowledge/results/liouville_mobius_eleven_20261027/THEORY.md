# A completely solvable model of joint incidence, transport, and guards

## Inheritance, portfolio and live board

Anchor: the operation-closed attachment carrier, now tested through arithmetic
observers. Niche: recovered occurrences of 11. Wildcard: explicit one-shot
circle-projector corners. The six live concepts are **support, sign transport,
joint guard, square dilation, Frobenius phase, native realization**.

Closest proved mechanisms: the [square-zero inventory, SQ1–SQ2](../squarefree_cycle_inventory_20261019/SQUAREFREE_CYCLE_BRIDGES.md)
and [native Boolean guard theorem A2–A3](../session70_attachment_guard_alignment_20261026/session70-vitali/REPORT.md).
Canonical hostile: the [81/27 C3-fibre pair, F2–F3](../tournament_composition_classification_20261026/fibre/RESULT.md).
Corrected near miss: [HYP-9029 / strong-interval-tiling-law](../../hypotheses/HYP-9029-strong-interval-tiling-law.md)
keeps global H-spectrum coverage open; the June reflection's claimed proof
does not follow from ruling out two low gaps. Least-used sidecar: indexed
Euler exponents behind the [June ghost atlas](../../../07-reflections/euler-product-ghosts-and-zeta-side-cancellation-gates.md).

The inheritance is mathematical: SQ1 already says to retain distinct prime
labels; A2 already realizes Boolean conjunctions in named re-marking. We add
an exact arithmetic continuation quotient and calculate new native readouts,
without presenting those older mechanisms as discoveries of this packet.

## D1. Two genuine dualities (classical, PROVED here elementarily)

Let Omega(n) count prime factors with multiplicity, lambda(n)=(-1)^Omega(n),
and mu(n)=lambda(n) if n is squarefree, otherwise zero. On positive integers,

    lambda(ab)=lambda(a)lambda(b),
    mu(ab)=mu(a)mu(b) [gcd(a,b)=1].

The second identity also holds when a or b is already nonsquarefree, since
both sides vanish. A sign and a guard are different coordinates.

The divisor-poset interpretation explains the guard intrinsically. For
`n=product p^e_p`, the interval[1,n] is the product of chains[0,e_p]. The
incidence inverse of the all-ones order relation on one chain has value1 at
length0, -1 at length1, and0 at every length>=2: the inverse recurrence gives
`-(1-1)=0` from length2 onward. Taking products gives arithmetic mu.
Liouville is the sign of the rank sum e_p. They coincide on Boolean divisor
intervals (squarefree n); a repeated prime makes a longer chain, whose
incidence cancellation is zero even though its rank still has a parity.
This identifies a precise source/target map to incidence algebra, beyond
merely observing that both functions take signs.

Write * for Dirichlet convolution, 1 for the constant-one function, epsilon
for its unit, s(n)=[n is a square], and u(d^2)=mu(d), with u zero off squares.
Then

    mu * 1 = epsilon,     lambda * |mu| = epsilon,
    lambda * 1 = s,       lambda = mu * s,       mu = lambda * u.

Proof: at a prime p the generating series of mu, lambda, s and u are
`1-z`, `1/(1+z)`, `1/(1-z^2)` and `1-z^2`; multiply the factors.
An equivalent genuine involution is `J(f)(n)=lambda(n)f(n)`:
complete multiplicativity gives `J(f*g)=J(f)*J(g)` and J^2=id.
It exchanges the inverse pair (1,mu) with (lambda,|mu|).

The analytic series are classical: for Re(s)>1, mu has Dirichlet series
1/zeta(s), lambda has zeta(2s)/zeta(s), and |mu| has zeta(s)/zeta(2s).
See [DLMF 27.4, equations 5, 7, 8](https://dlmf.nist.gov/27.4).
These identities are not conjectures or new literature claims.

## D2. Exact smallest multiplicative carriers (PROVED)

Fix k>=1 distinct primes S. The domain is the multiplicative monoid generated
by S. Two integers are equivalent for observer f when f(nx)=f(mx) for EVERY
x in that monoid. Repeated multiplication contexts reduce to one such x.

For mu, all nonsquarefree integers have the same zero continuation. Every
squarefree integer is represented by its support A subset S. Distinct A,B
are distinguishable: if p belongs to A but not B, multiplication by p gives
zero on A and a nonzero sign on B. A live state is distinguished from failure
by the empty continuation. This proves precisely 2^k+1 states, not just an
upper bound. Composition is disjoint union, with intersection sent to failure.

For lambda, equal signs remain equal after every continuation and opposite
signs differ already at x=1: precisely two states. For the JOINT observer
(mu,lambda), keep all squarefree supports and two nonsquarefree states
F_+,F_-. They remain nonsquarefree but transport lambda under multiplication.
All are reachable (p^2 and p^3 give the two failures) and distinguishable by
the preceding tests. Hence the joint minimum is 2^k+2.

With S={2,3,11}, the live integers are 1,2,3,6,11,22,33,66. Already 6 and22
have the same pair (+1,+1), but multiplication by3 gives mu(18)=0 and
mu(66)=-1. A summary of current output is not a continuation carrier.

The response matrix on live supports is

    R[A,B]=mu(n_A*n_B)
          =(-1)^(|A|+|B|) [A intersect B is empty]
          = K tensor ... tensor K,
    K = [[1,-1],[-1,0]],       K^-1=[[0,-1],[-1,-1]].

Thus R is unimodular and has rank 2^k over every field. The lambda response
is the outer product of the two sign vectors and has rank1. Over Q, adjoining
both outputs adds one new row direction for failure sign, giving rank2^k+1:
the live mu rows are independent and the two failure rows are opposite nonzero
lambda rows. These ranks omit an explicit constant mass coordinate. Appending
that coordinate gives ranks 2, 2^k+1 and 2^k+2 respectively over Q, enough to
recover distributions on all states. No exponential lower bound on the number
of deterministic memory *bits* follows; k support bits suffice on live states.

Generalization: for any nonzero multiplicative sign/weight character with
prime values c_p in a field, the guarded local matrix is
`[[1,c_p],[c_p,0]]`, determinant -c_p^2. Exponential rank arises from
independent exclusion guards, even when the transport factors are trivial.

## D3. Joint laws, square-zero algebra and the Vitali mechanism (PROVED)

For k=3 compare the uniform laws on the following squarefree integers:

    E: {1,6,11,66},          O: {2,3,22,33}.

Each prime has inclusion probability1/2 under both laws. Both mu and lambda
have the SAME full sign distribution: half+1, half-1. Nevertheless, after
multiplication by66, expected mu is -1/4 under E and0 under O. The only
surviving support is empty. A later joint guard detects hidden correlation.

For general k>=3 let the parity of the first k-1 support bits be even or odd,
and leave the last bit independent and uniform. Every marginal on <=k-2 bits
agrees between the two laws, and the sign distributions agree because of the
last free bit. Multiplication by the product of all k primes gives expected
mu=(-1)^k/2^(k-1) versus0. This proves that no bounded marginal order, even
augmented by the sign law, handles every later multiplicative guard.

The exact algebra is `A_S=Q[x_p:p in S]/(x_p^2:p in S)`. Its squarefree basis
x_A multiplies by disjoint union or zero. The linear functional
`phi(x_A)=(-1)^|A|` realizes R[A,B]=phi(x_A*x_B). The pairing is nondegenerate.
In the algebraic terminology A_S with phi is a Frobenius algebra; this word
does NOT identify its structure with finite-field Frobenius squaring.
Crucially phi is not an algebra homomorphism: phi(x_p^2)=0 whereas
phi(x_p)^2=1. This is the precise type warning already present in SQ1.

Relabel prime coordinates as named alignment bits. The disjointness tests
`[A intersect B=empty]` become the Boolean zeta tests after complementing the
state support. Diagonal signs then convert that matrix into R. This is a
literal invertible coordinate connection to native Vitali guard moments.
It preserves conjunction testing and rank; it does not map integer products
to tournament complements or make the two state spaces the same object.
Native tournament realization is supplied by the inherited A2 construction,
at n=7k+2, not by primality labels themselves. Unnamed averaging changes the
operation alphabet and may erase particular tests; its minimal carrier is
not decided by this named-guard theorem.

## D4. Square dilation, analytic scale and Euler ghosts (PROVED)

Let M(x)=sum_(n<=x)mu(n), L(x)=sum_(n<=x)lambda(n). D1 gives the exact identities

    L(x)=sum_(d<=sqrt(x)) M(x/d^2),
    M(x)=sum_(d<=sqrt(x)) mu(d)L(x/d^2).

For each fixed alpha>1/2, an O(x^alpha) bound for either summatory function
implies the same bound for the other, since sum d^(-2alpha) converges.
At alpha=1/2 this absolute-value proof loses log x. This identifies an exact
scale boundary of the comparison, without proving the required cancellation
bound or making a sign assertion about either summatory function.

Define the FORMAL ordinary Euler product
`P_b(q)=product_(n>=1)(1-q^n)^b(n)`. Each coefficient uses finitely many factors.
The same identities now give

    P_lambda(q)=product_(d>=1) P_mu(q^(d^2)),
    P_mu(q)=product_(d>=1) P_lambda(q^(d^2))^mu(d).

Proof: the exponent on (1-q^n) on the first right side is
sum_(d^2|n)mu(n/d^2)=lambda(n); the inverse is identical with u. No analytic
convergence is required. If G_b(n)=sum_(a|n)a*b(a), then

    G_lambda(n)=sum_(d^2|n) d^2 G_mu(n/d^2).

This connects the June distinction between exponent, ghost and coefficient
sparsity to an exact transform. Factor multiplicities survive; factor order
does not. Substituting lambda for mu as the inverse of divisor summation
instead returns the square-dilated original data, not the original data.

No fixed residue clock determines either function on all positive integers:
m and2m coincide modulo m but have opposite lambda; 1 and(m+1)^2 coincide
modulo m but have mu values1 and0. This is an elementary obstruction to a
fixed-modulus observer, not an undecidability claim. Both functions are
computable, and finite prime banks give the exact finite machines of D2.

## D5. Native tournament readouts and their limits (FINITE-EXACT)

Let H count directed Hamilton paths, with no identification by reversal.
Ordered join satisfies H(T->U)=H(T)H(U). Therefore lambda composed with H is
an exact join character, and mu composed with H is a character with the
coprimality guard gcd(H(T),H(U))=1. A five-tournament with outgoing bit rows
(6,28,24,17,1) has H=11. Joining it with C3 gives H=33; joining two copies
gives121. Joining two C3 gives9. The latter two have lambda=+1 but mu=0.
This realizes the {3,11} guard inside literal tournaments. Rédei's odd H
excludes the prime2 from H-values, so the full {2,3,11} arithmetic bank is
not asserted to be an H-factor bank.

For the twelve-vertex pair define tau_P by keeping the eleven directed edges
of P and reversing every other arc. We enumerate ALL legal P, quotient their
path-normalized presentations by source automorphisms, and multiply each
presentation's response by the free fibre size81 or27. The universe contains
462915 paths per source and5715 or17145 presentations. All black/blue types
are included; this is not a black-only statistic or a uniform-simple-edge walk.

Both sources have H=3^6*5*127. Their mu responses are -2916 and297, their
lambda responses -3321 and-1323. The normalized mu means are -4/635 and
11/17145; the visible numerator11 is an output of this calculation, with no
separate invariance theorem claimed. Their inherited second responses L^2(1)
are independently reproduced as202898835819 and202869738135.

This supplies a closure obstruction as well as a separator. Both sources
lie in the SAME joint arithmetic failure state F_+. Failure is absorbing
under ordered join, but native path complement can restore squarefree H and
can do so with different signed distributions from these two sources.
Consequently even D2's minimal JOINT multiplicative carrier is not closed
under the actual counted complement operator. Enlarging the operation alphabet
requires refining its states. Keeping prime support and sign solves one
composition problem; it does not replace the metagraph's path incidence data.
Indeed the full integer H is equal here, so even retaining its complete
factorization cannot determine these next-complement responses.

Every saved partner-H histogram bin is checked using a separate linear sieve;
the C++ producer instead trial-factors each H. A recursive Python path count
checks two sources and16 partners. The complete presentation traversal and
all partner H values are produced by one C++ implementation; the arithmetic
audit does not silently become a second full traversal.

## What unifies, and what remains to prove

The proposed research principle is **classify the response algebra before
discarding support or phase**. D2 is a solved special case, not a theorem that
all our objects have a small finite carrier.

| Source -> target | Map / preserved predicate | Lost information / required sidecar | Cheapest decisive test |
|---|---|---|---|
| integers -> two signs | lambda, multiplicative transport | prime identities and repeated-support guard | 6 vs22, multiply3 |
| prime supports -> square-zero basis | named primes become variables; no repetition | actual magnitude and primality remain external | repeated prime becomes zero |
| support laws -> Vitali coordinates | complement subsets, then diagonal signs | native operation names and transport maps | E/O laws above |
| tournaments -> arithmetic readouts | T maps to mu(H),lambda(H) | joint ports, source/target alignment, colours | actual81/27 complements |
| binary decorations -> finite field | chosen normal basis; rotation becomes squaring | phase gauge and full quotient automorphisms | degree11 vs modulus11 |
| divisor data -> Euler ghosts | indexed divisor sum | ordered history and unseen high coefficients | wrong inverse at121 |

For arbitrary contexts, observational equivalence must remain a congruence for
the allowed operations. The repo already supplies that abstract criterion.
The open task is its economical realization for unnamed complement/re-marking
with joint incidence, transport and realization guards. The present finite
arithmetic model shows why transport alone can collapse to rank1 while the
guards require exponentially many joint response coordinates. The 81/27
experiment demonstrates a useful new observable, not completion of that task.
