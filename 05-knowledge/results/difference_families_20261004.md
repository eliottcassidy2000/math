# Difference families and recursive Collatz certificates

2026-10-04. Owner-directed continuation of the tiling, four operations,
golden denominator and signed Collatz investigations.

**PROVED** below: elementary difference-pair arithmetic, triangular splice
law, a golden fundamental domain, the primitive-phase correspondence, and
the exact 19-branch phase tower; lattice-carry termination for supplied
golden-field inputs. **FINITE-EXACT**: the stated complete
denominator censuses and integer filters. **PROPOSED / OPEN**: the combined
family carrier and a normalizer covering every integer. The companion
[rational-anchor compiler](rational_anchor_returns_20261004.md) proves new
families of descents and completed signed certificates. Neither positive
Collatz nor completeness of the three known negative basins is proved.
No new canon theorem ID or mathematical priority claim is assigned.

The main result is a separation of three useful structures. Difference
pairs provide infinitely many representations of each integer. Golden
phases exhibit the requested recursion beyond 76, with exactly 19 child
cycles per parent. Guarded rational-anchor expressions actually certify
Collatz routes. Keeping all three, with their maps and guards, supplies a
more meaningful enriched numeral than an untyped collection of invariants.

## Inheritance and the working board

The closest proved mechanisms are:

- [THM-4528, parity in base phi](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md)
  and [the signed golden-carrier note](collatz_golden_carriers_20261004.md):
  exact phase, conserved denominator, ordered affine carry, and finite
  classifications at q=1,2,11,76.
- [THM-4512, coefficient-descent classes](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md):
  the extra valuation bit belongs in every route guard.
- [THM-4507, finite-valuation obstruction](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md):
  finitely many polynomial valuation coordinates, including primes 2,3,11,
  do not repair logarithmic height into a global nonincreasing potential.
- [The rational-shadow inventory](fibonacci_two_copies_collatz_20260926.md)
  and [shadow-error analysis](collatz_shadow_error_flp_20260927.md):
  -19/11 and its word (1,1,2) are already present.
- [THM-4139, rational three-cycle and order-six lift](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md),
  [THM-4146, three-cycle lift and fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md),
  [quadratic time reversal](arithmetic_braids_20260917_geometry.md), and
  [the three-to-six seam](arithmetic_seams_20260921_dynamics.md).
- Incoming commit 75daccc120:
  [reversible golden carries](golden_digit_carry_20261003.md),
  [prime clocks](golden_prime_clocks_20261003.md), and
  [mixed return compilation](collatz_mixed_return_compiler_20261003.md).
  Its distinction between a primitive coefficient pair and a ring unit is
  essential in the phase count below. The two notions are not identified.

Canonical hostiles: the rational cycle 1/13 sharing the norm-11 denominator
ideal with -5; arbitrary common shift above the -1 cycle; a legal growth
word with insufficient final division. Corrected near misses: the omitted
oddness bit, the golden upper boundary, the signed norm versus quotient
cardinality, and the attempted horizontal-fibre identification in
THM-4139/4146. The underused sidecar here is the conjugate real coordinate
of a+b*phi, which makes the phase count exact.

| Lane and live concept | Representation and operation | Decisive question |
|---|---|---|
| Anchor: certified families | Rational anchor, repeat count, exit guard | Does the endpoint pay the original source or reach a named root? |
| Difference as family | Nonnegative pair with fixed difference | What survives translating both endpoints? |
| Niche: golden phases | Conjugate pair in an L-shaped domain | Is a symbolic phase also an integer Collatz cycle? |
| Recursive scale | Primitive vectors modulo 4*19^k | How many child cycles lie above one parent? |
| Two operations | Triangular splice and radix subdivision | Are the two uses of multiplication correctly typed? |
| Wildcard: quadratic anchor | Difference, sum, chosen cycle phase | Which coordinate controls the next difference? |
| Golden constants | Trace, prime product and Higgs expressions | Is there a map, an exact identity, or only a numerical comparison? |

The assigned research target takes the anchor lane; LRC(14) remains open
and is not improved here. Used research cards: type the analogy; recover
the hidden coordinate; distinguish local coverage from universal coverage.
No new meta-pattern is promoted from this one session.

## Difference families with two useful projections

For n in Z define the infinite fibre

    D_n={(n_++k,n_-+k): k>=0},
    n_+=max(n,0), n_-=max(-n,0), pi(a,b)=a-b.          (1)

Every member maps to the same integer n. The exact observation behind the
7/11 and 8/12 comparison is

    (11,7),(12,8) in D_4.

They differ by translating both endpoints by one. Keeping the common shift
k retains their different locations; quotienting by simultaneous translation
forgets it. This is a precise family interpretation of the selected
expression n-r: if r is fixed, n-r determines n and does not by itself
create an infinite fibre. Both endpoints must vary together to keep a fixed
difference, or the anchor must be recorded as a separate coordinate.

There are natural commutative operations on the pairs:

    (a,b)+(c,d)=(a+c,b+d),
    (a,b) odot (c,d)=(ac+bd,ad+bc).                   (2)

The projection pi carries them to ordinary integer addition and multiplication.
Associativity and distributivity follow by identifying (a,b) with a+b*e,
e^2=1, before evaluating e=-1. This is classical difference arithmetic,
not a new construction of Z.

The second projection s(a,b)=a+b is equally useful: it is also additive
and multiplicative under (2). The pair is recovered from (s,n) by

    a=(s+n)/2, b=(s-n)/2,
    s>=|n|, s=n mod2.                                (3)

Thus a fixed difference n cuts out an infinite ray with totals
|n|,|n|+2,|n|+4,..., whereas a fixed total cuts out a finite diagonal.
These are two different projections of the same pair space. They should
not be silently interchanged when relating differences to a tournament size.

For the user's A<B split, set a=B+1, b=A+1. Then s=N=A+B+2 and n=B-A.
The finite size-N diagonal is the triangular split chart; moving along
D_n changes N by even increments while preserving the block imbalance n.
This gives a concrete meaning to a family of integers organized by one
difference.

There is also a useful categorical reading: regard (a,b) as the arrow
b -> a. The integer is its displacement, and D_n consists of all translates
of an arrow of displacement n. Sequential composition remembers the middle
endpoint:

    (a,b) compose (b,c)=(a,c).

Parallel addition is the first operation in (2). They satisfy the typed
interchange identity

    [(a,b) compose (b,c)]+[(d,e) compose (e,f)]
      =[(a,b)+(d,e)] compose [(b,c)+(e,f)]
      =(a+d,c+f).

This gives two genuine ways of composing a construction. Sequential
composition has one identity (b,b) at each object and is defined only when
endpoints match. There is no single common-unit operation on all arrows
to which the usual Eckmann-Hilton conclusion could be applied. Forgetting
the endpoints sends both compositions to addition of displacements and
loses precisely their different gluing instructions. This is a concrete
reason to keep the family members rather than only their integer labels.

## The triangular f law and a typed interpretation of g

With T(j)=j(j+1)/2 and N=A+B+2, the user's identity is exactly

    T(N-1)=T(A)+T(B)+(A+1)(B+1).                     (4)

It partitions the edges of K_N into the edges within blocks of sizes A+1
and B+1 and the cross rectangle. For tournaments, the cross rectangle
stores (A+1)(B+1) independent orientation bits. It must be retained if
one wants to reconstruct the tournament from the two smaller pieces.
A<B selects unequal unordered splits; equality is also valid in (4).

The splice A circle B=A+B+1 is associative, commutative, with formal unit
-1. Its cross term c(A,B)=(A+1)(B+1) satisfies

    c(A,B)+c(A circle B,C)
       =c(B,C)+c(A,B circle C).                     (5)

This cocycle identity says that three blocks have the same total number
of cross edges whichever two blocks are joined first. It is a precise
f-side law supplied by the clarification. The unparenthesized string
`(A+B)f1fAfB` still does not uniquely specify a scalar binary function f;
no arbitrary parse is imposed here.

The decimal split has the exact forms

    10^2=T(9)+T(10)=45+55,
    T(19)=T(9)+T(9)+10^2=190.                        (6)

For N=5 a fixed Hamiltonian path accounts for 4 of the 10 edge decisions,
leaving 6. The [tiling note](tiling_modular_atoms_20261003.md) owns the
full path presentation and the 15-cell mod-10 wedge; (4) supplies an
additional block decomposition, not a new count of isomorphism classes.

For g, the literal scalar reading of the clarified law is

    g(A,B)=AB*g(A,B).

Over the integers or reals this forces g(A,B)=0 whenever AB!=1. An
information-bearing alternative exists **at the family level**. For
positive A,B and c=AB, Euclidean division k=c*t+j gives a bijection

    D_n  <->  {0,...,c-1} x D_n,
    (n_++k,n_-+k) <-> (j,(n_++t,n_-+t)).             (7)

Every branch still projects to n. The right side is c disjoint labelled
copies of the family; it is not scalar multiplication of its integer value.
Iterating (7) records radix-c digits of the common shift. This realizes
the requested self-similarity without forcing the represented integer to
zero. It is a proposed typed interpretation of g, not a unique recovery
of its original definition. For finite k the digit record is eventually
zero; genuinely infinite digit addresses require a completion and need
not themselves be ordinary integer representatives.

The two operations have distinct roles: (4)-(5) build and reassociate
blocks; (7) subdivides a fibre while retaining the branch. This is why
an Eckmann-Hilton collapse cannot be inferred: a common carrier, compatible
units and the requisite interchange law have not been supplied. The
[paired-family note](paired_prime_collatz_20261003.md) already computes
the explicit interchange defect for addition and transported odd multiplication.

## A Collatz lift exists but the spare coordinate need not descend

The following map projects exactly to ordinary signed Collatz:

    if a-b is odd: (a,b) -> (3a+1,3b),
    if a-b is even: (a,b) -> (floor(a/2),floor(b/2)). (8)

In the second case a,b have the same parity, so their floor difference is
(a-b)/2. This proves the projection identity for every nonnegative pair.

The common part is not a Lyapunov function. Over the -1 cycle begin at
(k,k+1), k>=1. After its two ordinary steps,

    k -> floor((3k+1)/2),

which strictly increases: 1,2,3,5,8,12,18,27,... . The represented integer
returns to -1, while its extra coordinate grows. The failed implication is
that adding information automatically provides a decreasing rank. The
surviving structure is an exact semiconjugate lift. A useful certificate
needs a proved rewrite and a terminal or lower-source link.

This motivates the [rational-anchor compiler](rational_anchor_returns_20261004.md).
It proves such links rather than trying to make the arbitrary common shift
descend. It adds a -19/11 chart, a positive first-descent family starting
21223 mod32768, and symbolic power expressions reaching each signed root.

## An L-shaped fundamental domain for the golden phase

Put phi=(1+sqrt(5))/2, psi=1-phi, r=phi-1=-psi, O=Z[phi]. Retain the upper
boundary convention of the signed golden note:

    beta(x)=phi*x-d, d=1 if phi*x>1, d=0 otherwise.

This section concerns exact rational-field points x in Q(phi). If
x=(a+b*phi)/q in lowest scalar denominator, its conjugate is
y=(a+b*psi)/q. Multiplication by phi acts on residues by

    G=[[0,1],[1,1]], (a,b) -> (b,a+b) modq.           (9)

The relevant lattice and closed domain are

    Lambda={(m+n*phi,m+n*psi):m,n in Z},
    D=([0,1] x [-r,1]) union ([0,r] x [-phi,-r]).    (10)

**Fundamental-domain lemma.** D tiles the plane by Lambda, up to boundaries.
Every nonzero rational lattice coset has exactly one representative in D.

**Proof.** The lattice covolume is sqrt(5). The two rectangles have total
area phi+r=sqrt(5). If two interior points differ by
(u,v)=(m+n*phi,m+n*psi), then |u|<1 and |v|<phi^2. Hence
|n|< (1+phi^2)/sqrt(5)=phi, so n is -1,0,1. For n=0, m=0. For n=1 the
only possible m are -1,-2; m=-2 has v=-phi^2 and is excluded. The remaining
nonzero candidate is (r,-phi), or its negative. If (x+r,y-phi) were also
interior, x>0 would put its first coordinate above r, forcing its second
coordinate above -r, hence y>1, a contradiction. Thus interiors do not
overlap. Equal area gives full measure coverage on the torus. The locally
finite union of closed translates is closed, so a missed point would give
a positive-area hole; there are none.

For exact scalar denominator q>1, none of the horizontal or vertical
boundary coordinates is possible. Each boundary value and its conjugate
is an algebraic integer; equality of either coordinate would force
x in O and q=1. This removes boundary multiplicity from every nonzero
coset of q^(-1)O/O. QED.

The two-dimensional extension is

    F(x,y)=(phi*x-d,psi*y-d).

It maps D bijectively away from boundaries. Split the top rectangle at
x=r. Its left part maps to (0,1) x (-r,r^2); its right part maps to the
lower rectangle. The lower rectangle maps to (0,1) x (r^2,1). These
images partition D. The additional cut coordinates are algebraic integers,
so exact q>1 again avoids every cut.

**Primitive-phase theorem.** For every q>1, the purely periodic golden
points of exact scalar denominator q correspond bijectively to

    V_q={(a,b) modq:gcd(a,b,q)=1}.                    (11)

Their exact periods are the periods of G on these vectors.

**Proof.** The lemma supplies one representative of every primitive coset.
F preserves the finite set and projects to the invertible residue map G.
Uniqueness of representatives makes it a conjugacy, with exactly the same
periods. Conversely, for any periodic golden word extend the digits in
both directions. Since the first coordinate returns, its algebraic
conjugate returns and has the convergent expression

    y=-sum_(k>=0) psi^k*d_(-1-k).

Bounding the negative even powers and positive odd powers separately gives
-phi<=y<=1. If the preceding digit is 0, the lower bound improves to -r.
If it is 1, the next first coordinate is at most r. Thus the point lies
in D and is already counted. QED.

The q=1 boundary is genuinely exceptional: beta has the fixed point 0
and the period-two orbit 1,r, all with the trivial residue class. Do not
apply (11) to q=1. The inherited denominator-one classification treats
those boundary choices explicitly.

**Literature placement, CITED.** This is an elementary specialization of
the classical natural-extension and tiling framework for Pisot expansions;
see Kalle and Steiner, [Beta-expansions, natural extensions and multiple
tilings associated with Pisot units](https://arxiv.org/pdf/0907.2676),
version 3, Example 3.14 and section 4.3. The elementary domain and counting
proofs used here are supplied above. No novelty is claimed for that framework.

### A terminating difference normalization inside the golden field

There is a useful synthesis with the original difference proposal. For any
x in Q(phi) intersect (0,1), of exact q>1, let x_bar be the unique point
of D with the same residue. Their difference is an algebraic integer:

    x=x_bar+u+v*phi, u,v in Z.

Run beta on both. If their current digits are d and d_bar, the offset obeys

    (u,v) -> (v-(d-d_bar),u+v).                       (10a)

This is an exact carry equation: the finite periodic phase runs in V_q,
and an unbounded integer pair records the deviation from its periodic
representative. It retains precisely the data lost by reduction modulo q.

**PROVED termination on field inputs.** Every such offset eventually
becomes (0,0). The conjugate coordinate obeys y_next=psi*y-d, so it is
eventually bounded because |psi|<1. The first coordinate stays in [0,1].
At fixed denominator q, the two bounds leave finitely many possible
coefficient pairs. The deterministic orbit therefore becomes periodic.
By the primitive-phase theorem every periodic point is in D, where its
residue has the unique representative x_bar. Thus the offset vanishes.
This is the elementary quadratic Pisot mechanism in this coordinate system.

For the certified source 151,

    Theta(151)=(586-361*phi)/2
              =Theta(1)+(293-181*phi).

Its reference phase begins at Theta(1)=phi/2. The offset (293,-181)
becomes zero after the fifteen ordinary steps of 151 -> 227 -> 341 -> 1
(three odd steps and twelve divisions). Another control is
Theta(-9)=(30-12*phi)/11, with reference numerator (8,-1) and offset (2,-1).
Both phase and offset are verified against the actual parity digits.

The scope is crucial: termination is proved for **supplied field inputs**.
Constructing Theta(n) in Q(phi) from an arbitrary integer n without already
knowing its future remains the missing part. The formula does not establish
that every positive n has q=2. It gives a total normalizer on the proposed
algebraic carrier, while surjectivity of its integer-realizable sealed part
onto all positive integers is OPEN. This states exactly where a new
certificate-generating theorem would have to enter.

## Prime products and the precise role of 6 over pi squared

Inclusion-exclusion in (11) gives exactly

    |V_q|=J_2(q)=q^2*product_(p|q)(1-1/p^2).          (12)

A pair primitive in this sense is not necessarily a unit of O/(q).
Primitive means not both coefficients are divisible by any p|q; a unit
requires gcd(a^2+a*b-b^2,q)=1. For example the -5 numerator (-1,7) is
primitive modulo 11 but has norm -55 and is not a unit. The incoming
prime-clock note explains this distinction and the two split prime branches.

The all-prime version is the density of coprime pairs of ordinary positive
integers:

    lim_(X->infinity) #{1<=a,b<=X:gcd(a,b)=1}/X^2
       =product_p(1-1/p^2)=1/zeta(2)=6/pi^2.          (13)

For the density, use the sum of mu(d)*floor(X/d)^2; after division by X^2,
absolute convergence supplies the limit. The Euler product and
zeta(2)=pi^2/6 are classical identities, recorded in
[DLMF 25.2](https://dlmf.nist.gov/25.2) and
[DLMF 25.6.1](https://dlmf.nist.gov/25.6.E1).
Thus 6/pi^2 has a precise place in this project: it is the limiting
primitive-pair sieve behind the finite golden phase counts. It is not the
proportion J_2(q)/q^2 for any arbitrary fixed q.

Selected complete computations, each with integer realization separately
checked:

| Exact q | Primitive phases | Period : number of cycles | Nonzero integer Collatz roots |
|---:|---:|---|---|
| 2 | 3 | 3 : 1 | 1 |
| 3 | 8 | 8 : 1 | none |
| 4 | 12 | 6 : 2 | none |
| 5 | 24 | 4 : 1; 20 : 1 | none |
| 10 | 72 | 12 : 1; 60 : 1 | none |
| 11 | 120 | 5 : 2; 10 : 11 | -5 |
| 19 | 360 | 9 : 2; 18 : 19 | none |
| 29 | 840 | 7 : 4; 14 : 58 | none |
| 38 | 1080 | 9 : 6; 18 : 57 | none |
| 76 | 4320 | 18 : 240 | -17 |
| 100 | 7200 | 60 : 20; 300 : 20 | none |

Here 10*10=100 has a second useful reading: two modular coefficients give
100 addresses modulo 10, of which 72 are primitive. Modulo 100 there are
10000 addresses, of which 7200 are primitive. These are two-coordinate
residue counts, distinct from the multiplication-table wedge and (6).

The primes 2,3,11 retain several different roles: parity and ternary guards;
the inherited halving totals 2,3,11 of roots 1,-5,-17; and the inert or split
prime behaviour in O. Prime ordinal positions 1,2,5 do not determine these
maps. The phase theorem now supplies a common explicit residue carrier in
which the prime distinctions can actually be tested.

## Beyond 76 there is an exact tower with 19 branches

For k>=1 define q_k=4*19^k. Then every primitive golden phase at exact
denominator q_k has period

    L_k=18*19^(k-1),
    J_2(q_k)=4320*19^(2k-2),
    number of cycles=240*19^(k-1).                  (14)

Reduction modulo q_k sends each cycle at q_(k+1) onto a cycle at q_k.
**Every parent has exactly 19 child cycles**, each 19 times longer.

**Proof.** Modulo 4 every primitive vector has period 6: modulo 2 its
period is 3, while G^3-I=2G is nonzero on it modulo 4 and G^6=I modulo 4.
Modulo 19 the eigenvalues are 5 and 15, of orders 9 and 18. Every nonzero
vector therefore has period 9 or 18. CRT gives period 18 at q_1=76.

The exact integer identity

    G^18-I=76*G^9=19*U, U=4*G^9                    (15)

has U invertible modulo 19. The binomial expansion implies, for every
vector v nonzero modulo 19 and every s>=1,

    min valuation_19 of the coordinates of
       ((G^(18s)-I)*v) = 1+v19(s).                  (16)

For completeness, taking a nineteenth power raises the exact matrix
valuation by one and preserves an invertible leading matrix modulo 19;
taking a power coprime to 19 multiplies that leading matrix by a unit.
This proves (16), not merely the order of the matrix. Since every primitive
period must already be a multiple of 18, its exact lift is L_k. Formula
(12) gives the counts in (14). There are 19^2 lifts of each vector at the
next level. Over a parent cycle of length L_k these occupy
19^2*L_k/(19*L_k)=19 child cycles. QED.

| k | q_k | Common period | Number of cycles |
|---:|---:|---:|---:|
| 1 | 76 | 18 | 240 |
| 2 | 1444 | 342 | 4560 |
| 3 | 27436 | 6498 | 86640 |

This is an exact recursive continuation of the 76 phase layer. It is
the most concrete meaning found here for the proposed microcosm/macrocosm
boundary. The recursion lives in the golden residue dynamics; it does not
assert that all subsequent phases are integer Collatz cycles.

That distinction has a decisive finite test. Decode each of the 4560 cycles
at q=1444 into its raw parity word, with p odd letters, A divisions, and
ordered carry B. Its rational Collatz cycle starts at

    n=B/(2^A-3^p).                                   (17)

The complete computation finds **no integer realization** among those
cycles, covering all 1559520 primitive phases. At q=76 the same filter
retains the known -17 cycle. The extra 19 branches therefore do not simply
produce a fourth negative integer root. No all-k exclusion is claimed.

The computation also explains a stronger failure at this first lift.
Every q=1444 word has between 83 and 106 odd letters in its 342 raw steps.
A negative rational cycle would require 3^p>2^(342-p), hence p>=133.
Thus **every cycle at this exact denominator is a positive rational cycle**;
none can serve as a negative growth anchor. This is a complete finite
statement at k=2, not an extrapolation to all levels.
At q=76 there are 26 negative rational cycles, only one integer cycle.
Reduction between the two levels therefore does not preserve the sign of
the decoded Collatz root. The phase quotient retains period information
while losing the parity weight that determines that sign.

Another cheap extrapolation also fails: the magnitudes 1,5,17 continue
under a ->3a+2 to 53, but

    -53 -> -79 -> -59 -> -11 -> -1 -> -1

under U. This preserves a recursive integer sequence while losing the
claim that each term names a different terminal cycle.

## The quadratic parameters through the difference coordinate

For f_c(x)=x^2+c and a pair x,y, let delta=x-y and s=x+y. Then exactly

    delta_next=s*delta,
    s_next=(s^2+delta^2)/2+2c.                        (18)

This explains why a fixed difference needs a second coordinate: (2,1)
and (3,2) both have difference 1, but their squared differences are 3 and
5. If r_i is a marked cycle of f_c and x=r_i+delta, the anchored form is

    delta_next=delta*(2r_i+delta), phase i -> i+1.    (19)

For an affine Collatz word W with fixed rational anchor r, the corresponding
law is W(n)-r=(3^p/2^A)*(n-r). Both are exact displacement charts. The
quadratic term in (19) and the Collatz valuation guards prevent treating
them as a global conjugacy.

The inherited quadratic trace formulas, for chosen three-cycle trace sigma
and eta=2sigma+1, are

    c=-(sigma^2+sigma+2)=-(eta^2+7)/4,
    reversal: eta -> -eta-1,
    other cycle at the same c: eta -> -eta.           (20)

The first operation transports a chosen cycle by L(x)=-x-1/2 and reverses
its order. The second changes which cubic factor is selected; it need not
preserve rational splitting. Their composition is a unit translation of
eta, a valid family operation with a necessary field/splitting sidecar.

At eta=0 the parameter is -7/4 and the two traces meet; the inherited
conductor-seven construction yields a parabolic algebraic three-cycle.
At eta=-1/2 the reversal fixes the chosen trace, giving -29/16 and

    -7/4 -> 5/4 -> -1/4 -> -7/4.

This is the inherited unique rational arithmetic-progression three-cycle,
with multiplier 35/8. The same -7/4 appears here as a **point**, whereas
in the preceding sentence it was a **parameter**. Starting x^2-7/4 at zero
is a third, different orbit; its first values -7/4,21/16,-7/256 do not close.
The conductor-nine companion reverses to parameter -11/4 with multiplier
19. That 19 is a quadratic multiplier; a conjugacy with the tower (14) has
not been established.

There is a more direct arithmetic web around the new compiler anchor:

    rational Collatz: -19/11 -> -23/11 -> -29/11 -> -19/11,
    scaled 3n-11:       19   ->   23   ->   29   ->   19,
    final numerator: 3*29-11=76=4*19,
    golden phase: Theta(-19/11)=(-4+20*phi)/29.

Here scaling by -11 and the parity evaluation are explicit maps. The
golden denominator 29 follows from phi^7-1=7+13*phi, of norm -29. This
explains why a golden phase with no integer cycle can nevertheless be a
useful rational anchor for **integer** route families. The compiler retains
the denominator-11 integrality guard instead of discarding that phase.

The same last numerator 29 and total division factor 16 suggest comparing
with c=-29/16. An affine identification of the two three-point cycles
fails: the rational Collatz cycle has unequal successive sorted gaps
4/11 and 6/11, whereas the quadratic cycle is an arithmetic progression.
The numbers alone do not supply a dynamical conjugacy. The surviving
connection is the explicitly typed anchored-coordinate framework (18)-(20).

These are useful atlases of anchored dynamics. Their common content is
the exact displacement update with retained phase, not an assertion that
quadratic dynamics proves Collatz convergence.

## Four labels and the tournament question

The inherited golden denominator labels of the four known nonzero signed
cycles are

    root:          1, -1, -5, -17,
    denominator:   2,  1, 11,  76.

Their number alone does not supply six pairwise orientations of a tournament.
One legitimate observable is the ratio A/p of total halvings to odd steps:

    -1: 1,  -5: 3/2,  -17: 11/7,  1: 2.

Orient from smaller to larger ratio. There are no ties, and the result is
the transitive four-vertex tournament. The vertices are the four chosen
cycles, the observable is this clock ratio, and reversing the comparison
reverses all arcs. It preserves their clock order but loses the ordered
carry, phase and integer-realization condition. It therefore supplies no
classification of further cycles. Reachability among distinct terminal
cycles supplies no alternative tournament, because none reaches another.

The productive tournament connection here is instead the exact triangular
splice (4), where the missing cross-edge record has a defined type and size.
The earlier [route-tournament carrier](collatz_atom_route_memory_20261003.md)
already stores intrinsic route information in ordered components. Both
constructions retain the information their integer counts would discard.

## The two lambda formulas must remain distinct

The expression in the current request is

    lambda_current=(5*phi^4)^(-1/6)
                  =0.554855513344434017... .

The [September 30 note](collatz_circulant_20260930_circulants_lucas_cubic_monotile.md)
used a different expression:

    lambda_previous=1/[5*(2*phi^4)^(1/6)]
                   =0.129280563443397023... .

The exact golden trace phi^4+phi^(-4)=7 gives respectively

    25*lambda_current^12-35*lambda_current^6+1=0,
    976562500*lambda_previous^12-218750*lambda_previous^6+1=0. (21)

These algebraic identities are proved by putting t=lambda^6 and taking
the conjugate of phi^(-4). They organize the constants within the same
quadratic field used by the phase reader.

For a numerical comparison only, the Standard Model tree relation is
m_H^2=2*lambda*v^2; lambda is an input, not a predicted constant of the
Standard Model. See the [PDG 2025 Higgs review, equation 11.8](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-higgs-boson.pdf).
Using the same illustrative v=246.22 GeV gives 259.375097730384 GeV for
the current expression and 125.200177018301 GeV for the previous one.
The latter is the source of the earlier numerical Higgs comparison.
Neither (21), the phase tower nor the prime-pair density derives a physical
Higgs coupling; precision comparisons would also need scale and scheme.

## The proposed carrier and its exact research obligation

The most useful proposal is an **atlas of guarded difference families**.
Every state retains an integer projection and any known chart:

    integer representative (a,b), pi=a-b;
    anchor and phase, displacement from that anchor;
    repeat word and count, exact source guards;
    exit budget measured against the original source;
    ordered affine carry and timed golden polynomial;
    optional primitive golden phase and its lifting address;
    endpoint pointer and open or sealed proof status.

Every integer has an open state. A sealed state must contain a verified
finite route, a proved parameterized macro followed by a sealed suffix, or
a chain of positive lower-source certificates terminating at 1. Negative
states carry an explicit cycle and arrival phase. The representation may
be redundant, but each redundancy has a checkable map; arbitrary unrelated
annotations are not counted as a proof.

This is more informative than an integer in a useful sense: the integer
value forgets the chart and proof, while the structured numeral retains
how to certify its computation. A power expression can represent an enormous
integer together with a short repeated-block proof. The rational-anchor
companion supplies twelve actual signed examples and all-parameter rules.

| Source | Map and target | Preserved predicate | Lost data and required sidecar | Cheapest test |
|---|---|---|---|---|
| Pair | (a,b)->a-b | Integer addition and multiplication | Common shift; retain k | Growing lift above -1 |
| Block split | (A,B)->A+B+2 | Vertex count | Cross orientations; retain rectangle | Three-block cocycle |
| Family | k->(k modc,floor(k/c)) | Projected integer | Branch if forgotten; retain digit | Reconstruct k |
| Golden point | x->(a,b) modq | Period, for q>1 in D | Boundary at q=1; retain convention | Independent trap enumeration |
| Primitive phase | Reduction q_(k+1)->q_k | Parent cycle | Child branch and phase | 19 children per parent |
| Golden cycle | Word->B/(2^A-3^p) | Formal cyclic itinerary | Integrality and sign; retain guard | Complete q=1444 filter |
| Rational chart | Source->exit by a guarded repeated word | Actual Collatz path and proved return | Unproved suffix; retain pointer | Weakened-budget source 487 |
| Quadratic pair | (x,y)->(delta,s) | Exact coupled dynamics | Chosen cycle if using an anchor | Same difference, different sums |

The live board changes in three ways after the computation. The conjugate
coordinate turns a denominator census into an all-q phase theorem. The
19-adic lift gives a precise recursive structure but immediately encounters
the independent integer filter. The rational anchor supplies actual new
return rules, so it overtakes the purely symbolic tower as the most direct
route toward the certification target.

The outstanding question is now precise: can an orbit-independent normalizer
assign every positive input a certified lower-source rule, using an extensible
atlas of rational charts and finite heads? The analogous negative question
must also certify the terminal cycle. Neither a countable list of charts nor
arbitrarily large finite coverage is sufficient. A total construction or
a well-founded coverage argument is still required.

## Reproduction and audit scope

Run from the repository root:

    python3 04-computation/experiments/difference_families_20261004.py --large \
      > 05-knowledge/results/difference_families_20261004.out

The [script](../../04-computation/experiments/difference_families_20261004.py)
uses exact integer or rational arithmetic for all mathematical decisions.
Decimals occur only in the explicitly numerical lambda displays.
The [saved output](difference_families_20261004.out) records:

- All exact denominators 2..30,38,76,100: the domain points and every
  primitive residue agree with the earlier independently implemented trap
  enumeration. This is 32 complete denominator comparisons.
- Complete phase and integer-cycle counts for the eleven displayed q values,
  and all 4560 cycles / 1559520 phases at q=1444.
- Matrix-order controls at all eight levels k=1..8 of (14). The all-vector,
  all-k assertion rests on (16), not on matrix-order sampling.
- 28561 ordered pair arithmetic controls and 169 Collatz lifts; the growing
  common-coordinate hostile; 400 triangular splits, 8000 cocycle triples,
  and 63000 radix reconstruction instances.
- 722 rational quadratic-pair controls, the exact -29/16 cycle and its
  multiplier, the -53 hostile, and exact coefficient-pair checks of (21).
- Golden carry normalization for every exact-denominator point with
  q in {2,11,76}, -40<=a,b<=40 and 0<(a+b*phi)/q<1, plus the two signed
  itinerary controls 151 and -9. The output records its count and maximum
  entry time; the all-field termination proof is separate.

The rational compiler has a separate reproduction block and explicit
universes in its companion. Checks remain active under optimized Python.
This is a self-audited research note with independent computational paths,
not an assertion of an independent human or agent review. No Lean theorem
or universal Collatz result is claimed.

Repository-wide documentation validation also reports the pre-existing
maintained-prefix budget failure in `05-knowledge/hypotheses/INDEX.md`
(125 lines, budget 120). That file is unchanged in this session; it does
not affect the exact research replays above.
