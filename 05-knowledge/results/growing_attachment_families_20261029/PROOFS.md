# Growing attachment families: rank, degree, conductor and actual roots

2026-10-29. **PROVED deductions from CITED classical inputs; same-agent
proof review; FINITE-EXACT controls.** No literature-priority claim, new
canon ID, or resolution of general complement/re-marking closure.

## P1. A finite polynomial bank, including every dependent query

Let q be odd, chi its quadratic character, and f_1,...,f_m nonzero
polynomials over F_q. Factor each as a nonzero constant times powers of
monic irreducibles. Let V be the binary matrix of their factor exponents,
and s the vector of quadratic-character bits of the constants. Put r=rank V.
Let Z contain every F_q-zero of every f_i, let z=|Z|<q, and let D be the
sum of the degrees of the irreducibles appearing in V with a nonzero column.
Factors occurring only to even powers still contribute to Z, though not D.
We may enlarge Z explicitly; its cost must remain in z.

On U=F_q minus Z the answer vector y, defined by chi(f_i(x))=(-1)^y_i,
belongs to the affine code C=s+im V, of size 2^r. For c in C write N_c
for its number of realizations. If r>0, then

    |N_c - (q-z)/2^r| <= (1-2^(-r)) ((D-1)sqrt(q)+z).       (1)

For r=0, there is exactly one answer vector, realized q-z times.
Vectors outside C never occur. V determines possible algebraic patterns;
it does not alone guarantee that all of them occur in a small finite field.

**Proof.** Choose a row-space basis of V and represent it by monic
squarefree polynomials g_1,...,g_r. The original answers are an injective
affine encoding of the r bits chi(g_j(x)). For each nonzero u in F2^r,
the squarefree representative h_u of product g_j^u_j is nonconstant,
has degree at most D, and is not a constant times a square. On U its
character agrees with that product. Weil's bound gives

    |sum_(x in U) chi(h_u(x))| <= (D-1)sqrt(q)+z.

Expand the indicator of any specified basis-bit vector by character
orthogonality. Its trivial term is (q-z)/2^r and its 2^r-1 other terms
obey this bound. This proves (1). The cited bound is stated in Jones,
*Infinite Paley graphs*, section 6, Remark 1, equation (3):
<https://arxiv.org/html/1912.02709v1#S6>. QED.

This proof is the usual character projector with its rank and deleted
points made explicit. It applies to every finite bank, with no assumption
that its displayed polynomials are independent. Infinitely many products,
repetitions, powers and nonzero constant multiples drawn from a fixed
finite factor support are handled by the same finite code. Arbitrary square
multiples are free only where their zeros have also been excluded: adding
new zeros can exhaust U. Rational functions have the same version after
retaining all zeros and poles and their parity valuations.

**Boolean closure without a union penalty.** Set

    rho=(2^r-1)((D-1)sqrt(q)+z)/(q-z).

Every nonempty subset B of C has

    |{x in U:y(x) in B}| = |B|(q-z)/2^r * (1+theta),
    |theta|<=rho.                                         (2)

Sum (1) over B; the relative error stays rho. Thus one estimate controls
all Boolean predicates in the displayed answers simultaneously. It does
not replace arbitrary operations on the underlying field points.

## P2. Joint colored attachments in an actual Paley tournament

Now q=3 mod4 and x->y means chi(y-x)=1. There are no ties for distinct
vertices; chi(-1)=-1 makes this an intrinsic tournament. Use a bank from
P1. Fix h distinct named roots outside Z. Prescribe a feasible bank color
for each of k ordered new vertices, every root-to-new incidence, and a
labelled tournament on the new vertices. Keep root-root incidences fixed.

At step j, 0<=j<k, there are t=h+j actual earlier vertices. The polynomials
X-a for these vertices are distinct and coprime to every bank factor:
a is outside Z. Their classes are independent modulo the bank row space.
The enlarged rank is exactly r+t, its support degree is D+t, and its
deleted set has exactly z+t points. With n_t=2^(r+t), put

    M_t=(q-z-t)/n_t,
    E_t=(1-1/n_t)((D+t-1)sqrt(q)+z+t)                      (3)

when r+t>0. If r+t=0 put E_t=0. Provided h+k<=q-z, the total
number L of these injective rooted colored copies obeys

    product_j max(0,M_(h+j)-E_(h+j)) <= L
      <= product_j (M_(h+j)+E_(h+j)).                     (4)

**Proof.** Apply P1 at each actual partial embedding. The same bounds
hold for every partial embedding, not just for an average root tuple.
Multiply lower and upper extension counts by induction. Empty partial
branches present no problem; a positive lower bound certifies extension.
Every root-new or new-new sign is a genuine Paley incidence, with signs
reversed when the requested arc points towards the earlier vertex. QED.

In particular, for sequences of banks and roots, uniformly over all their
placements and all feasible specified patterns,

    L = (1+o(1)) q^k / 2^(rk+hk+binom(k,2))               (5)

whenever k>=1, h+k<=q-z, and

    (D+z+h+k+1) 2^(r+h+k) / sqrt(q) -> 0.                (6)

Indeed the sum of the relative errors against q/n_t is bounded by a
constant times the expression in (6): use a geometric sum in 2^j, and
sum the deleted-point main-term corrections k(z+h+k)/q. Products of
1+o(1) factors then give (5). This deliberately simple sufficient budget
is not claimed sharp. It allows both the bank and the pattern to grow.

The same relative estimate holds after summing over ANY fixed subset of
the joint color/incidence/tournament answer space, because all complete
patterns have the same main term. No extra union bound over this answer
space is needed. Such predicates may vary with q; they must still be
predicates of these specified answers, not of unrecorded vertex identities.

**Linear special case and classical credit.** Without a polynomial bank,
for any t distinct roots and any signs, the one-vertex count satisfies

    |N - q/2^t| <= (t/2-1+2^(-t))sqrt(q)+t/2.            (7)

Expand the t linear factors: the singleton sums vanish and the degree-j
sums cost (j-1)sqrt(q). At each deleted root the normalized indicator
contributes at most 1/2. Equation (7) is classical, recorded explicitly
by Jones section 6, Remark 2, with credit to Bollobas--Thomason; the
tournament construction goes back to Graham--Spencer (1971),
<https://doi.org/10.4153/CMB-1971-007-1>. We do not claim it as new.
Sequential multiplication makes its joint counting consequence explicit:
all root-plus-new sizes h+k <= (1/2-eta) log_2(q), eta>0 fixed, satisfy
(5) with r=D=z=0. This is a sufficient window, not an optimal threshold.
Necessarily 2^t<=q-t if every root-incidence pattern is to have a witness.

## P3. An infinite multiplicative family and its transport

**Every finite binary code, with actual attachment generators.** Let
C be any affine binary code of length m and dimension r. Choose a full
column-rank matrix M:F2^r->F2^m and shift s with C=s+im M. Choose r
distinct field points a_l and a nonsquare nu. Put

    f_j(X)=nu^(s_j) product_l (X-a_l)^(M_jl).              (7a)

Then the characters of these polynomials are exactly the affine encoding
of the r actual root-to-x Paley incidence bits. P1 applies with D=z=r,
and every word of C is realized whenever

    q-r > (2^r-1)((r-1)sqrt(q)+r),                        (7b)

for r>0. For r=0 the constant code is realized everywhere. This realizes
every finite binary linear matroid's dependence data as well, by using
any binary matrix representation as the query exponent vectors. The
square-total kernel and the character-answer image are dual codes, not
the same code. Fano-type dependence therefore fits without an ordinary
graph-cycle representation. There is no assertion that an individual
product query is itself a new tournament vertex: the generators are
actual roots, while products are arithmetic observables of their joint
incidences. The static carrier uses exactly r bits once (7b) holds.

For any odd p and any finite set of pairwise distinct monic irreducibles
over F_p, retain their factors, parity matrix and zero set. Over F_(p^e)
their nonempty squarefree products remain nonsquares: splitting cannot
change simple roots into even multiplicities. Thus their original binary
independence persists even when individual irreducibles split. All
polynomial products supported on them are governed by P1. As odd e tends
to infinity, constant-character signs from F_p stay fixed. For p=3 mod4,
P2 simultaneously applies to the entire Paley family q=p^e, e odd.

For a growing example, take s distinct linear factors, all 2^s products
as queries, and any powers supported on those factors. Then r=s, D=z=s.
Any h roots outside those s points and any k new vertices satisfy the
uniform law if (2s+h+k+1)2^(s+h+k)=o(sqrt(q)). The number of displayed
queries may be exponential in s, or an unbounded redundant list; the
charged information is s bits together with its realization data.

Frobenius phi:x->x^p preserves every Paley arc. Transport coefficients,
roots, excluded points and vertices simultaneously: f^phi(phi(x))=
phi(f(x)) and chi(phi(u))=chi(u). This gives a bijection of each named
joint-response set with its transported response set. Fixing coefficients
in F_p simplifies this to root/vertex transport; it does not authorize
forgetting named-root phase. The characteristic-two F64/F256 carriers
from earlier work are different objects, not quadratic Paley tournaments.

## P4. Prime attachments with a growing conductor

This extends Session107's **fixed** modulus theorem, not its proof of the
large-sieve bound: [Theorems A/B](../session107_quasirandom_arithmetic_20261029/characters/CHARACTER_STRUCTURE.md).
Let a bank of nonzero signed integer bases vary with X. Let Q=8 times the
product of ALL odd primes dividing any base, including square factors.
Let r be the rank of the exponent-parity matrix on positive prime factors,
and C its affine sign code on primes p=3 mod4. Retain only primes not
dividing Q. For any fixed A>0 suppose Q<=(log X)^A.

Uniformly over all such banks and all c in C,

    #{p<=X:p=3 mod4, color(p)=c}
       = Li(X)/2^(r+1) * (1+O_(A,B)((log X)^(-B)))        (8)

for every fixed B>0, with ineffective constants. The prime count here
excludes divisors of Q. This is a deduction from the classical
Siegel--Walfisz theorem; its ineffectivity precludes extracting a certified
numeric starting height from this argument.

**Proof.** CRT gives exactly phi(Q)/2^(r+1) reduced classes for each
color. Siegel--Walfisz gives pi(X;Q,b)=Li(X)/phi(Q)+O_(A,L)(X/log^L X),
uniform in Q and b. Sum over those classes and divide by the main term.
The relative error is O(phi(Q) log^(1-L)X); choose L>A+B+1.
Tao's author notes, Exercise 64, explicitly state the uniform modulus
range and ineffective constants:
<https://terrytao.wordpress.com/2014/12/09/254a-notes-2-complex-analytic-multiplicative-number-theory/>.
The pi version follows by partial summation (apply the uniform statement
with a larger fixed exponent for intermediate heights). QED.

As in (2), (8) holds for arbitrary nonempty Boolean unions of colors,
with the same relative error. The length of a redundant bank is immaterial
once Q and r are charged. Multiplying by new square factors may enlarge Q:
its extra excluded primes cannot be discarded as a supposedly free bit.

Choose independent prime vertices in any k such nonempty Boolean fibers.
Their sizes are at least (1-o(1))Li(X)/2^(r+1). Inherit Session107 A:
edge error <= (B_X/2)sum_(i<j)1/sqrt(n_i n_j), where
B_X<<_epsilon X^(1/2+epsilon); collision error <= binom(k,2)/min n_i.
Consequently the probability of any labelled k-tournament, conditional
on distinctness, is (1+o(1))2^(-binom(k,2)), uniformly, if

    2^(r+binom(k,2)) k^2 log(X) X^(-1/2+epsilon) -> 0.   (9)

For Q<=log^A X we have 2^r<=phi(Q)/2<=log^A X. Thus every fixed eta>0
allows k <= (1-eta)sqrt(log_2 X), by choosing epsilon small enough.
This is a sufficient spectral window, not an optimal prime-pattern law.
For fully specified colors their joint probability also has the factor
2^(-rk); the product of k relative errors in (8) tends to zero here.

An explicit growing family avoids any unproved smooth-prime assumption:
take the largest initial segment of odd primes whose product times 8 is
at most log^A X, include 2, and permit every signed integer supported on
these primes. Its independent rank tends to infinity. Every finite list
of queries, and every Boolean predicate in its sign answers, is covered
at once. This says nothing about the distribution of Sylvester prime
terms or factors, which form a different, often thin source ensemble.

## P5. Recovered niche: an Apollonius discriminant family

The older [generic Apollonius squareclass proof](../apollonius_squareclass_code_20261019/APOLLONIUS_SQUARECLASS_CODE.md)
left numerical distribution separate from its generic parity group. Here
is a literal one-parameter discriminant specialization to which P1 applies.
Take labelled centers (0,0),(t,0),(0,1), and radii 1,2,3. The six defect
factors in that proof become

    F12+=t^2-1, F12-=t^2-9, F13+=-3, F13-=-15,
    F23+=t^2,   F23-=t^2-24.

Modulo squares its four discriminants are

    -3 A,  -15 A C,  -3 B C,  -15 B,
    A=t^2-1, B=t^2-9, C=t^2-24.                          (10)

If char(F_q) is outside {2,3,5,23}, A,B,C are squarefree and pairwise
coprime, so their squareclasses are independent. The displayed four rows
have rank 3 and sole nonzero relation their sum. The product of their
constants is a square, hence the affine answer code is precisely the
eight even-parity words. Remove Z={t:t A B C=0}, z<=7. P1 gives each word
count (q-z)/8 with error at most 7/8*(5sqrt(q)+z), uniformly for this
entire infinite family. Every odd-parity word is impossible exactly.

This is a proved distribution law for the four labelled discriminant
characters on this family. It is NOT a theorem counting eight distinct
geometric circles at every surviving t: quadratic-leading-coefficient,
zero-radius and solution-distinctness guards would also be needed for that
claim. The map and its lost information remain explicit. The excluded
characteristics have a mechanism: coincident or inseparable factors, or
zero constants. For example in characteristic 23, C=A and rank drops.

## P6. What the 81/27 control still forbids

The two inherited twelve-vertex tournaments share their labelled cyclic
triples but have automorphism orders 81 and 27. P2/P4 count their distinct
labelled patterns individually; they never identify their response carriers.
With no named roots and all moving colors equal, unlabelled counts multiply by 12!/|Aut|,
so the group-27 class has asymptotic frequency three times the group-81
class in these families. Root incidences can further break symmetry;
then use the automorphism group preserving the entire specified pattern.
Inherited exact source certificates remain in
[DEFECT_TRANSFER D6](../overnight_rank_return_20261029/DEFECT_TRANSFER.md#d6-full-three-run-matrices-of-the-8127-pair).
No new twelve-vertex response computation is claimed here.

**Relation to A000568.** Once every codeword is realized, its 2^r named
colors form a palette to which Session107's exact colored Burnside formula
[COLORED_RECURRENCE](../session107_quasirandom_arithmetic_20261029/colored/COLORED_RECURRENCE.md)
applies with palette size 2^r. P2/P4 supply arithmetic realization and
asymptotic weights of its labelled colored patterns. They neither replace
that recurrence nor compute new terms of A000568; the arithmetic source
measure and the class count are different quantities. Actual metagraph
complements still require the inherited full target and joint port data.

Rank is minimal for static answers only when every algebraic word is
realized: 2^r distinct outputs require at least r bits. P1 proves this
realization when its lower bound is positive. Outside that range rank
may overstate finite realized information. Under complement, re-marking,
addition, or exact path queries, this static minimality proves nothing
about the smallest closed carrier. The full relations, ordered joint
ports, source identities and transport must still be retained.

**Open next targets.** Improve the sufficient growth windows using more
than one-edge spectral telescoping; bound rank/degree growth under an
actual complement compiler; classify the extra Apollonius realization
divisors and their resulting finite-field circle counts. These are separate
obligations from the family laws proved above.

## P7. All finite Collatz valuation words supply guarded linear banks

Fix sigma in {+1,-1} and any finite legal positive valuation word
d_1,...,d_J for the odd-only map n->(3n+sigma)/2^v2(3n+sigma).
Let D_j=sum_(i<=j)d_i, D=D_J, and use its exact positive odd source
cylinder n_0=b_0+2^(D+1)t, t>=0, with the least positive residue b_0.
Inherit this cylinder representation from Session107,
[fixed-word source guards](../session107_quasirandom_arithmetic_20261029/orders/ORDERS_CARRIES_AND_CRT.md).
Along it every intermediate integer is

    n_j(t)=A_j t+B_j,  A_j=3^j 2^(D+1-D_j),              (11)

where B_j is the actual intermediate value starting at b_0. This follows
by the literal affine recurrence, retaining the exact valuation at each
step. For i<j, A_i B_j-A_j B_i is nonzero: the segment from i to j has
the form (3^(j-i)n_i+sigma C)/2^(D_j-D_i), with C>0.
Its determinant is sigma A_i C/2^(D_j-D_i).

Thus the following nonzero integer is an explicit exceptional-prime
certificate for this word:

    H=6 product_(0<=i<j<=J) (A_i B_j-A_j B_i).            (12)

For every prime p not dividing H the J+1 linear polynomials (11) have
nonzero slopes and distinct roots modulo p. Among t=0,...,p-1, exclude
the J+1 positions where an n_j(t) is divisible by p. For every prescribed
sequence of quadratic signs along this ACTUAL integer valuation word,
its count obeys (7) with q=p and t there replaced by J+1. Scalar signs
chi_p(A_j) merely shift the bit pattern. The count is therefore
p/2^(J+1)+O(J sqrt(p)+J), uniformly over the word whenever p avoids H.
For J=0 the same assertion is the elementary one-linear-factor count.

This is an infinite family for EVERY finite legal valuation word for
either fixed choice of sigma, not a bank of individually discovered
examples. An arbitrary finite collection of such words has a common
exceptional certificate, the product of its H's. Formula (7) bounds all
their patterns simultaneously, with the maximum word length charged.
The exceptional certificate can be enormous; no small least-prime bound
is claimed. It excludes primes for this chart, not exceptional Collatz
orbits or exceptional starting integers.

This also clarifies the earlier Fourier warning. The actual pairs
(n_i,n_j) lie on a deterministic affine line modulo p, retaining perfect
joint correlations in additive Fourier coordinates. Their quadratic
signs can nevertheless be asymptotically independent when the two roots
are distinct. The map from values to signs destroys information. Source
cylinders, intercepts, ordered valuation clocks and the prime p remain
sidecars; these finite-word statements imply no termination, temporal
equidistribution of one orbit, or cycle classification for either map.
