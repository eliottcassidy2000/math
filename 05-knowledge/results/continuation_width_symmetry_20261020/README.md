# Continuation width, common symmetry support, and sparse defect words

20 October 2026. **PROVED:** the scoped continuation-width, symmetry, grammar,
and counting results below. **FINITE-EXACT:** the attached controls.
**OPEN:** polynomial complexity, a smaller quotient preserving structural
successor overlaps, and universal Collatz entry. The proof review and separate computational
representations were performed by the same agent; no independent-agent
audit or canon theorem-ID promotion is claimed.
This continues [the finite-horizon response recurrence](../tournament_future_boundary_20261020/README.md)
and [the fixed-feedback compiler](../defect_boundary_20261020/sparse_theory.md).
No historical-priority claim is made.

The anchor is uniform continuation width and symmetric sparse-defect grammars.
The niche is the difference between an individual permutation's support and
the common support of a group. The wildcard is a source-preserving affine
defect grammar for Collatz, including the obstruction to transferring a
tournament bound from edit count alone.

## Inheritance and live board

- Closest proved mechanisms: the future-response quotient; the individual
  `3k` support and `4k+1` canonical clone-core bounds; the exact affine carry
  cocycle and the local rotation-defect ideal.
- Canonical hostiles: independent prime coordinates; different histories with
  identical mass but different future gcd response; a periodic binary word
  with two ones and arbitrarily long zero gaps; an arbitrarily long lawful
  expanding Collatz word.
- Corrected near miss: the preceding session's scalar clique-separator
  quotient, and the current two-edit norm barrier, where a large resultant
  does not equal the actual simultaneous-root content.
- Least-used sidecar: the *whole* automorphism orbit partition, and the
  positive-definite continuation kernel suggested by incoming matrix-valued
  metagraph moments.

| Object | Coordinate retained | Cheapest decisive test |
|---|---|---|
| A000568 variable elimination | every still-admissible insertion response | independent prime-subset histories |
| Linear summaries | complete continuation matrix | exact rank of a principal submatrix |
| Sparse tournament symmetry | full group orbits and an optimal order | exhaustive small automorphism groups |
| Near-extremal symmetry | `delta=3f-s`, marked skeleton and ordered gaps | reconstruct all arcs and marked ports |
| Collatz defect words | source, exponent order, affine carry, integer clock | direct integer trajectory versus compressed evaluation |
| Metagraph compression | joint Gram matrix and forgotten directions | compare with incoming covariance obstruction |

## 1. A uniform bound and a superpolynomial obstruction

Use the descending odd-length elimination of the preceding packet, with
degree horizon N. At each stage the states retain committed mass d and
responses `R_j=sum_a gcd(a,j)` for unprocessed odd `j<=N-d`.
Let W(N) be its maximum number of occupied states over all stages.

**Theorem 1.** There is an absolute c>0 such that, for all sufficiently large N,

    exp(c sqrt(N/log N)) <= W(N)
                          <= (N+1) exp(pi sqrt(2N/3)).

The lower bound also holds, up to a polynomial factor, among histories with
the same mass and number of cycles. In particular this exact descending,
wordwise-continuation quotient has no polynomial state bound.
This does not rule out a better elimination order, symbolic summation, or an
algorithm preserving only the final scalar count.

**Upper bound.** There are at most `sum_(d<=N) p(d)` processed cycle
multisets, even if one allows all positive parts. For t>0,

    log product_(j>=1)(1-exp(-tj))^(-1)
      =sum_(r>=1) 1/[r(exp(tr)-1)] <= pi^2/(6t).

Thus `p(d)<=exp(tN+pi^2/(6t))` for d<=N. Optimize at
`t=pi/sqrt(6N)` and multiply by N+1. The case N=0 is separate.

**Lower bound.** Let L>=3 and let P be the odd primes in `(L,2L]`.
Write m=|P|, S=sum(P), and take N=4S. At the stage immediately before
processing the odd lengths at most 2L, all lengths `3p`, p in P, have
already been processed. For each subset A of P choose the history

    mu_A={3p:p in A}.

Every history has mass at most 3S, leaving at least S, so all the future
primes in P remain affordable. Since no prime in P is 3,

    R_1=|A|,    R_p=|A|+(p-1) 1[p in A].

The 2^m histories give different states. Pigeonholing by mass and
cardinality leaves at least `2^m/((3S+1)(m+1))` different histories in one
common mass/cardinality cell.

The only analytic input is the classical prime number theorem:
`pi(x)~x/log x`, [NIST DLMF 27.2.3](https://dlmf.nist.gov/27.2.E3).
It implies `m~L/log L` and `S=Theta(L^2/log L)`. Therefore
`m=Theta(sqrt(N/log N))`. For arbitrary large horizons choose
`L=floor(c0 sqrt(N log N))` with a sufficiently small fixed c0; then
`4S<=N` and the same lower estimate holds. Adding slack cannot remove
these distinguishable histories. QED.

### The obstruction survives linear combinations

**Theorem 2.** At the prime-subset stage above, an exact linear summary
preserving every permitted individual continuation has dimension at least
`2^m/((3S+1)(m+1))`.

Choose the common mass/cardinality cell F from the proof, with cardinality r.
For A,B in F, use future continuation `{p:p in B}`. It is affordable since
the total mass is at most 4S. The past-future edge-orbit interaction is

    K(A,B)=2^(r^2) product_(p in A intersect B) 2^(p-1).

Apart from positive row/column factors, this is the corresponding submatrix
of the actual weighted continuation operator. Put x_p=2^(p-1)>1.
The full subset kernel `H(A,B)=product_(p in A intersect B) x_p` is the
tensor product of the positive definite matrices

    [1  1]
    [1 x_p].

A principal submatrix is positive definite. Hence K on F has rank |F|.
A factorization through D linear coordinates has rank at most D, proving
the lower bound. This is a statement about a specified linear input/output
interface, not about all computational encodings or arithmetic complexity.
QED.

The incoming [metagraph moment theorem](../merged_metagraph_distance_moment_constraints_20261020/PROOF.md)
uses Gram positivity to detect losses after a fibre quotient. Here positivity
proves that the chosen future probes have no nontrivial linear null direction.
The carriers and predicates are different; the common mechanism is explicit.

### Reroute: eliminate largest-prime batches

The recurrence itself permits any order of odd cycle lengths. Order them
by decreasing largest prime factor, with length 1 last. Within each prime
batch any order is allowed. Retain the same degree-pruned future responses.
Let Wprime(N) denote the peak state count for this order.

**Theorem 2A (uniform improvement).** As N tends to infinity,

    Wprime(N) <= exp(O(sqrt(N) / log N)).

Thus this order has an asymptotically smaller upper bound than the proved
lower bound for the descending-length order. This concerns numbers of
states and, with polynomial overhead as justified below, bit-operation
cost. It is still superpolynomial as an upper bound and supplies no
polynomial algorithm. A version using only the elementary integer
Euler-product bound gives `exp(O(sqrt(N) log log N/log N))`;
the sharper displayed estimate uses the classical prime-count bound
`pi(t)=O(t/log t)`.

**Proof.** Consider a stage within the batch of an odd prime p. The
processed history splits into old cycles with largest prime factor greater
than p, and newly processed cycles with largest prime factor exactly p.
For an old length a, let c be its largest p-smooth divisor. Every future
length j is p-smooth, so `gcd(a,j)=gcd(c,j)`, and `a>p c`.
For a new length a, retain `b=a/p`, so `a=p b`.
The response vector is determined by the two multisets of c's and b's,
and their combined mass is at most N/p. Keep actual committed degree
as an additional coordinate; it must not be replaced by this kernel mass.

Writing M=floor(N/p), the number of these two-coloured partitions of mass
at most M is at most

    (M+1) exp(2pi sqrt(M/3)).

This follows from the two-colour Euler product by exactly the upper-bound
argument of Theorem 1. Hence at every stage of this batch,

    states <= (N+1)(M+1) exp(2pi sqrt(M/3)).          (4)

There is also a complementary bound. Write Psi(N,p) for the number of
p-smooth positive integers at most N, including 1. There are at most
Psi(N,p) future response coordinates, each an integer between 0 and N.
Consequently

    states <= (N+1)^(Psi(N,p)+1).                    (5)

Put v=log N, u=log log N,

    Y=(v/(512u))^2,    sigma=1/2-1/u.

For sufficiently large N, `1/4<=sigma<1/2` and Y>=2. Rankin's elementary
product bound and the sum over all integers give

    log Psi(N,Y)
      <= sigma v + log product_(prime q<=Y)(1-q^(-sigma))^(-1)
      <= sigma v + 16Y^(1-sigma).

Indeed `-log(1-t)<=t/(1-t)`, `2^(-sigma)<=2^(-1/4)<7/8`, and
`sum_(j=2)^floor(Y) j^(-sigma)<=Y^(1-sigma)/(1-sigma)<=2Y^(1-sigma)`.
Also `log Y/u<=2`, so `Y^(1/u)<=exp(2)<8`. Therefore

    log Psi(N,Y) <= v/2-v/u+16*8*v/(512u)
                 = v/2-3v/(4u).                    (6)

If p>=Y, (4) is `exp(O(sqrt(N) u/v))`, absorbing polynomial factors.
If p<Y, (5)-(6) have logarithm at most
`(sqrt(N) exp(-3v/(4u))+1) log(N+1)`, which is
`o(sqrt(N)u/v)`. Processing the final unit length uses at most
polynomially many degree/cycle-count states.
These alternatives cover every intermediate stage, not merely completed
prime batches. No prime-number input enters this first upper bound.

For the sharper claim, the prime number theorem cited in Theorem 1 gives
`pi(t)<=C t/log t` for all t>=2 after enlarging C. Uniformly for
`1/4<=sigma<=1/2`, partial summation gives

    sum_(prime q<=Y) q^(-sigma)
         <= C1 Y^(1-sigma)/log Y.                  (7)

One explicit justification is to bound primes up to sqrt(Y) by all integers,
giving `O(Y^((1-sigma)/2))`, and use `log t>=log Y/2` in the remaining
partial-summation integral. Since 1-sigma>=1/2, its constants are uniform.
The first contribution is absorbed into the displayed bound.

Now choose `Y=(a v)^2`, where a>0 is a sufficiently small fixed constant,
and keep `sigma=1/2-1/u`. The logarithm of the Euler product is at most
`8C1 Y^(1-sigma)/log Y <=8C1 exp(2) a v/u` for large N.
Choose `a<=1/(16C1 exp(2))`. Then

    Psi(N,Y) <= sqrt(N) exp(-v/(2u)).

For p>=Y, (4) has logarithm O(sqrt(N)/v). For p<Y, (5) is smaller.
This proves the sharper upper bound without invoking a smooth-number
asymptotic. QED.

The proof also applies to the signed clone-removal and fixed-point-free
core recurrences: their actual cycle mass is at most the committed degree,
and the extra substitution degree is retained. Cancellation can remove
occupied states, not add response profiles outside these bounds.

**Bit-cost corollary.** All coefficients through N of each of these three
series can be computed in `exp(O(sqrt(N)/log N))` bit operations
using elementary exact arithmetic, with a possibly different implied constant.
There are O(N) stages and at most O(N^2) transitions per state; response
vectors have O(N) entries of O(log N) bits. Each partial-history denominator
`z_mu` divides `|mu|!`, hence N!. Its orientation exponent is at most
`binom(N,2)`. Clone-substitution factors have absolute value at most
`2^(2N)`, and even the crude history count `(N+1)^(2N)` suffices to bound
the scaled numerator by `2^(O(N^2+N log N))`. Aggregate rational coefficients
therefore have polynomially many bits. Exact additions, multiplications,
gcd reductions, hashing or comparison of profiles, and the polynomially
many outgoing choices give only polynomial overhead. This is an
asymptotic algorithmic bound; the finite runtime table is a separate control.

## 2. One support bound for the entire automorphism group

Let f=f(T) be the minimum feedback-arc number. Let G=Aut(T), let
O_1,...,O_t be its non-singleton vertex orbits, and put

    ell_i=|O_i|,    s=sum_i ell_i.

Thus s is the size of the union of the supports of ALL automorphisms.

**Theorem 3.**

    f >= sum_i (ell_i^2-1)/8 >= s/3.

In particular all automorphisms are supported on one common G-invariant
set of at most 3f vertices. Every vertex outside that set is fixed by G.

**Proof.** The action on each O_i is transitive, so its induced tournament
is regular. Consequently ell_i=2h_i+1 is odd and at least three.
In any linear order restricted to O_i, the last h_i vertices have at
least `h_i,h_i-1,...,1` outgoing edges directed backward, counted from
the end. These are distinct arcs, giving `h_i(h_i+1)/2`.
The internal arcs of different orbits are disjoint. Finally
`(ell^2-1)/8>=ell/3` for odd ell>=3, with equality only at ell=3.
QED.

This strengthens the inherited elementwise statement. A bound for each
permutation alone would not in general justify a bound for their union;
the regular *whole-group orbits* supply the missing argument.

## 3. Exact extremizers and a near-extremal symmetric grammar

Define the nonnegative integer symmetry slack

    delta=3f(T)-s.

It is a new parameter, different from the fixed feedback budget.

**Theorem 4 (zero slack).** Delta=0 if and only if T is an ordered join of
singleton and directed-triangle components. With exactly k triangles,
`f=k`, `s=3k`, and `Aut(T)=C3^k`. Consequently their bivariate ordinary
class series, including the empty object, is

    1/(1-x-y x^3),

where x records vertices and y triangles/feedback. At fixed n,k their
number is `binom(n-2k,k)`, interpreted as zero when n<3k.
The fixed-point-free automorphism mass of this family is

    1/(1-(2/3)y x^3).

The latter is a symmetry-weighted series, not a count of classes.

**Theorem 5 (bounded exceptional skeleton).** Every tournament of symmetry
slack delta admits the following exact representation:

1. A marked vertex set B of size at most 2delta, invariant under Aut(T).
2. An order on B and arbitrary orientations on its internal pairs.
3. In each of its |B|+1 ordered gaps, a word over singleton and directed
   triangle atoms.
4. All inter-atom and atom-to-B arcs agree with the displayed order.

The three vertices of each unmarked triangle form an actual G-orbit and
a homogeneous module. Storing actual vertex identities and the phase
inside each triangle preserves any marked source vertices, path, and ports.
This is an existence/reconstruction grammar, not an injective class census:
arbitrary generated templates need not have the asserted actual slack.

**Proof of both theorems.** Choose an optimal order and let

    c=sum_i (ell_i^2-1)/8,    d=3c-s.

Let e be the number of backward arcs whose endpoints belong to different
G-orbits (singleton orbits included). Internal orbit costs give
`e<=f-c=(delta-d)/3`.

Put in B every non-triangular non-singleton orbit, and every orbit meeting
an endpoint of one of these e external backward arcs. Orbits of size at
least five contribute at least `4ell/5` to d, since

    3(ell^2-1)/8-ell >= 4ell/5   (odd ell>=5).

Their total size is at most 5d/4. All other orbits have size one or
three. Each external backward arc adds at most six further vertices.
Therefore

    |B| <= 5d/4+6e <= 2delta-3d/4 <= 2delta.       (1)

Every remaining triangle orbit A meets no external backward arc.
For any other orbit C, all cross arcs agree with the optimal order.
The number of out-neighbours in C is the same at each vertex of A,
because G acts transitively on A and preserves C. In this order that
number is the number of C-vertices after the A-vertex. Hence no C-vertex
can lie between two A-vertices. This holds for every C, including singleton
orbits: the three A-vertices are consecutive, and every outside vertex
sees them uniformly. They form one triangle atom.

All remaining singleton vertices meet no backward arc. Thus only B-B
pairs and the interiors of triangle atoms can disagree with the order.
This proves the grammar and exact reconstruction.

If delta=0, (1) makes B empty. The whole tournament is an ordered join of
the stated atoms. Conversely such a join has exactly one unit of feedback
and three active vertices per triangle. Its strong-component order is
unique, each triangle has C3 automorphisms, and singleton components are
fixed. This proves the classification, counting, and product-group claims.
An automorphism is fixed-point-free precisely when there are no singleton
atoms and each triangle rotation is nonidentity, a fraction (2/3)^k.
QED.

### Scope and a useful sharpness boundary

The canonical clone core need not be bounded for fixed delta: an ordered
join of arbitrarily many triangles has delta=0 and is pair free.
The bounded object is the exceptional skeleton; the unlimited part is
carried by a word. This resolves the apparent tension with the 4k+1 theorem.
Delta=1 and delta=2 are impossible: if d>0 then d>=4, while if d=0,
s is a multiple of three and so is delta. Delta=3 is realized by the rigid
four-vertex tournament `C3[2,1,1]`; delta=4 by the regular five-tournament.

The grammar's gap words have generating function `1/(1-x-yx^3)`,
the same atom recurrence appearing in the incoming rigid/pure-blue
classification. Here it follows from a different extremal predicate:
common automorphism support versus feedback. No equivalence of the two
classes is asserted without the appropriate additional hypothesis.

### The full symmetry action, including ports, factors

**Corollary 5C (bounded exceptional symmetry).** For the constructed B,
let t be the number of quiet triangle orbits and u the number of remaining
singleton vertices outside B. There is a permutation group H acting
faithfully on B, with |B|<=2delta, such that

    Aut(T) = H x C3^t

in its displayed action. In particular its full permutation cycle index is

    Z_Aut(T) = Z_H * ((p_1^3+2p_3)/3)^t * p_1^u.    (8)

If a collection of individually fixed ports meets h distinct quiet
triangles, then

    Aut(T; ports) = H_(ports intersect B) x C3^(t-h).

The singleton ports outside B impose no further restriction. Thus phase
and port information can be carried exactly through the unbounded words.

**Proof.** Each quiet triangle is a homogeneous module and an actual
G-orbit, so its internal rotations extend independently to automorphisms
of T, and no element of G exchanges two such modules. Any automorphism
restricts to a rotation on each triangle. Remove these independent
rotations: the remaining automorphism fixes every vertex outside B.
Call this group H. Its action on B is faithful; its support is disjoint
from the quiet rotations, which commute with it. This proves the direct
product and the cycle-index product. A nonidentity triangle rotation
fixes no vertex, so marking any of its ports removes exactly that factor.
QED.

The fixed-point-free fraction is zero when u>0. Otherwise it is
`(der(H)/|H|)(2/3)^t`, with the empty-action fraction interpreted as one.
This extends the zero-slack formula to every fixed-slack representation.
Every prime other than three dividing |Aut(T)| is confined to the
bounded group H and therefore to the odd part of `(2delta)!`.
The theorem does not claim that every H, skeleton, or gap allocation is
realizable with the prescribed actual orbit partition.

### A uniform growth law for every fixed symmetry-slack layer

**Corollary 5A.** Let b_D(n) count unlabelled tournaments with delta<=D,
where D is a fixed nonnegative integer. Coefficientwise,

    [x^n] 1/(1-x-x^3) <= b_D(n)
      <= [x^n] sum_(m=0)^(2D) 2^binom(m,2) x^m/(1-x-x^3)^(m+1).

Consequently

    b_D(n)=O_D(n^(2D) lambda^n),    lim b_D(n)^(1/n)=lambda,
    lambda^3=lambda^2+1,    lambda>1.

**Proof.** The lower bound is the injective zero-slack word classification.
For the upper bound, order the at most 2D marked vertices, choose their
internal tournament in at most `2^binom(m,2)` ways, and fill every ordered
gap with a word of singleton/triangle atoms. This is a surjective
description, so duplicate representations only enlarge the upper bound.
Internal triangle labels may be standardized for this *unmarked class*
bound. Actual marked ports would require their addresses to be retained.

The smallest-modulus zero rho of `1-x-x^3` is positive, simple, and unique
on its modulus circle: equality in `1=|x+x^3|<=rho+rho^3` forces both
terms to be positive real. Thus its reciprocal lambda is the unique
dominant growth rate. Poles of order at most 2D+1 give the polynomial
factor in the upper bound. The simple lower series has a positive
constant times lambda^n. QED.

This gives a finite-description growth theorem even when the clone core,
feedback number, and moved support all grow without bound.

### Feedback also pays for the size of the symmetry group

**Corollary 5B.**

    |Aut(T)| <= 3^f(T),

with equality exactly for the zero-slack singleton/triangle words.

**Proof.** A tournament has no involutive automorphism, since a transposed
pair would have its edge reversed. Pairing group elements with their
inverses shows that |G| is odd. Its faithful action on its orbits embeds
G in the product of symmetric groups, so |G| divides the product of
their factorials and is at most the odd part of that product.

For odd ell>=3,

    oddpart(ell!) <= 3^((ell^2-1)/8),

with equality only at ell=3. The cases 3 and 5 are `3=3` and `15<27`.
For the step ell -> ell+2 with ell>=5, the odd-part multiplier is at most
`(ell+1)(ell+2)/2 <= 3^((ell+1)/2)`. The latter holds at ell=5 and
propagates because its left-side successive ratio is less than three.
This proves the factorial inequality by induction. Apply Theorem 3.
Equality excludes every orbit of size at least five and forces f=s/3.
Theorem 4 gives the stated words, whose product group attains equality.
QED.

### The first two nonzero slack layers have exact normal forms

Put `A(x)=1/(1-x-x^3)`. For a tournament U, write `C3[1,U,1]` for
vertices x,y and U, with `x -> U -> y -> x`.

**Theorem 6.** The only possible symmetry slacks below six are 0,3,4.
Their exact ordinary class series are

    delta=0:  A(x),
    delta=3:  x^2 (A(x)-1-x) A(x)^2,
    delta=4:  x^5 A(x)^2.                           (9)

At delta=3 the unique nontrivial strong component is `C3[1,U,1]`,
where U is a singleton/triangle word on at least two vertices. Arbitrary
singleton/triangle words precede and follow it. At delta=4, replace that
component by the regular five-tournament R5.
Here “nontrivial” means different from singleton and triangle atoms.

**Proof for delta=3.** The orbit charge d in Theorem 5 is zero, so all
active orbits are triangles. There is at most one external backward arc.
Between two triangle orbits, the cross orientation has constant row and
column sums. A nonuniform 3-by-3 matrix of this kind is a permutation
matrix or its complement. There are three directed alternating four-cycles,
one for each pair of matching edges; no single edge belongs to all three.
Every total order therefore has at least two backward cross edges,
contradicting the budget. Thus each triangle orbit is a homogeneous module,
also relative to all fixed singleton vertices.

For arbitrary homogeneous modules, feedback equals the sum of internal
feedback numbers plus minimum weighted feedback in the module quotient.
The proof is the inherited representative-sampling argument: in an
arbitrary order choose one uniformly random representative from each
module. The expected weighted cross cost of the induced block order equals
the original cross cost, and placing optimal internal orders in the blocks
attains the converse inequality. Transitivity of the modules is not needed.

Here the module weights are 3 and 1, and the residual weighted feedback is
exactly one. A minimizing block order has just one backward arc, between
two singleton vertices x,y. All atoms between them form U; the outer
atoms form two ordered words. U cannot be empty, since that would leave
a transitive quotient. U cannot be a single singleton, since x,U,y would
then be a homogeneous triangle with independent rotations, contradicting
that x,y are actual fixed vertices.

Conversely let |U|=m>=2. In the resulting strong component, x has score m
and y score one. The maximum-score class contains x alone, or x and the
first singleton of U; in the latter case x is the source of the pair.
Similarly the minimum-score class contains y alone, or y and the final
singleton of U; y is the sink of that pair. Thus x,y are intrinsic and
fixed by all automorphisms, and U is recovered with its ordered atom word.
Its only automorphisms are the internal triangle rotations.
Feedback is one plus their number: the weighted quotient is nontransitive
and has the displayed order with one backward edge. This gives delta=3.
The strong-component order makes both outer words and the central
component unique, establishing the second series in (9).

**Proof for delta=4.** The only possibilities for d<=4 are d=0 or d=4.
The first would make delta divisible by three, impossible. Thus there is
one orbit of order five and all other active orbits have order three.
There are no external backward arcs. The contiguity argument of Theorem 5
applies to every orbit, so T is an ordered word of singletons, triangles,
and one regular five-tournament.

There is only one regular five-tournament up to isomorphism. Choose a
vertex v, order its out-neighbours a->b and in-neighbours c->d.
Regularity forces b->c,d, a->c, and d->a, giving the cyclic orientation
on (v,a,b,c,d). Its rotations act transitively, while an automorphism
fixing v fixes the two ordered neighbour pairs; hence its group is C5.
Its feedback number is three, by the orbit lower bound and the cyclic
displayed order. It has delta=4. Unique strong-component order proves
the third series.

Finally delta=1,2 cannot occur as noted above. At delta=5, either all
active orbits are triangles, giving delta=0 modulo 3, or there is exactly
one five-orbit, giving delta=1 modulo 3. Neither gives five. QED.

**A hostile exactly at six.** Take two directed triangles `(a_i)` and
`(b_i)`, with `a_i -> b_j` precisely when i=j. The simultaneous rotation
of both triangles is an automorphism. Their distinct score classes show
that every automorphism preserves the two triples; the matching locks their
phases, so the group is the diagonal C3, not C3 x C3.
The three alternating four-cycles above force at least two external
backward arcs; the internal triangles force two more. The displayed order
`(b_2,b_0,a_1,b_1,a_2,a_0)` has four, so f=4, s=6, delta=6.
Neither orbit is a module. This is the first permitted slack at which that
independence assumption can fail; the bounded exceptional group retains
the synchronization.

## 4. Collatz: what transfers and what fails

The tournament proofs use an order, disjointly charged violated edges, and
regular invariant orbits. A parity string with few changed digits has none
of these properties automatically. For example

    (1 0^L)^2

differs from the all-zero word in two positions but is preserved by a
rotation moving all `2L+2` positions. Thus a fixed edit budget cannot imply
a uniform rotation-support bound for binary words. The reference itself
has extensive symmetry.

Even an authenticated positive Collatz word can have one exceptional
exponent and an arbitrarily long quiet part: source `2^(2r+1)+1` has the
exact prefix `2^r 1`. After r steps it is `2*3^r+1`; every earlier
valuation is exactly two and the next is exactly one. A constant number
of defects therefore does not bound the integer clock. This family
descends substantially; it is a hostile to information deletion, not to
the Collatz conjecture.

The valid transfer is a representation principle: keep a bounded exceptional
part together with exact ordered lengths of the quiet pieces. For the
accelerated odd Collatz branch

    T_a(n)=(3n+1)/2^a,    a>=1,

set z=n-1. Then

    z -> alpha_a z+beta_a,
    alpha_a=3/2^a,    beta_a=4/2^a-1.

The quiet exponent 2 acts as `z -> (3/4)z`.
Every supplied exponent word with k non-2 exponents has the unique form

    2^r0 b1 2^r1 ... bk 2^rk,     bi!=2.

Store the bi, every gap ri, the actual source, and the total integer clock
`2 sum ri+sum bi`. The terminal displacement is exactly

    z_out=A z_in+C,
    A=(3/4)^(sum ri) product_i alpha_bi,
    C=sum_i beta_bi (product_(j>i) alpha_bj)
                       (3/4)^(ri+...+rk).           (2)

The exponent on the i-th carry starts with the gap AFTER bi.
All factors are rational; no real-clock approximation is used.
Splitting a quiet gap into smaller quiet gaps changes neither A nor C.
Actual integer admissibility still requires the stated valuation word;
the affine expression alone does not assert it.

Replacing each beta by its positive part gives Cplus>=C. Therefore

    A(n-1)+Cplus < n-1                              (3)

is a sufficient strict-descent certificate for an authenticated word.
Its margin survives every quiet-gap refinement. This elementary affine
certificate is a scoped implementation of the established carry cocycle,
not a new universal descent theorem.

### A completely explicit contracting block family

For a block `1^s 2^r a`, with s>=1,r>=0,a>=3, put
`u=(3/2)^s`, `t=(3/4)^r`. Its exact coefficients are

    A=3tu/2^a,    C=[6t(u-1)+4]/2^a-1.

If `2^a>=6t(u-1)+4`, then C<=0 and A<1. Thus every authenticated
positive odd source n>1 strictly descends under this block.
The strict slope comparison follows from
`6t(u-1)+4-3tu=4+3t(u-2)>=5/2`, since u>=3/2 and 0<t<=1.
This class includes every `1 2^r a` with a>=3.
For example source 19 has exponents (1,3) and reaches 11.
The inequality is sufficient, not necessary; failing it leaves the
source-specific test (3) or the exact signed carry test available.

The existing [two-edit endpoint classification](../collatz_two_swap_endpoint_exclusion_20261020/ENDPOINT_TWO_EDIT_CLASSIFICATION.md)
is much stronger in its cycle lane: it excludes entire all-length families
using simultaneous-root content and an audited logarithmic estimate.
The [local rotation-defect ideal](../collatz_rotation_defect_ideal_20261020/LOCAL_ROTATION_DEFECT_IDEAL.md)
retains clock sign, full prime-power content, and rotational phase.
Those results are inherited, not replaced by (2)-(3).

### Research targets, with hostile tests

- **OPEN, further width:** the largest-prime order evades the descending
  prime-subset cut and proves Theorem 2A. Determine its sharp width, or
  improve it by symbolic kernels. Theorem 2 is scoped to its particular cut.
- **OPEN, symmetric census:** normalize the fixed-delta skeleton grammar
  beyond the complete delta<6 normal forms, without losing its actual
  automorphism orbits or port identities. Delta=6 has a concrete synchronized
  triangle obstruction. A finite description is not already a deduplicated
  enumeration.
- **OPEN, Collatz entry:** prove that each positive source eventually admits
  a certified return using (2), a stronger existing controller, or an
  arithmetic content certificate. A bound on the number of exceptions
  without their placement, clock, and carry does not suffice.
- **Conditional useful target:** bound the discounted positive carries
  Cplus relative to `(1-A)(n-1)` at selected, authenticated return times.
  This is source-specific and refinement-stable, but establishing it for
  every source is still an unproved arrival/coverage obligation.

The inherited expanding word with exponents (1,2) has slope 9/8 and
positive carry. Its repetitions give arbitrarily long lawful expanding
prefixes at suitable positive integers. They remain the hostile to any
claim of bounded-depth universal descent from a fixed quiet alphabet.

## Validation and concurrent work

Exact checks, additional hostile controls, source hashes, and incoming
integration are recorded in the accompanying files. No independent-agent
audit is claimed. Finite controls do not establish the all-order proofs.

The [validation and proof audit](AUDIT.md) records the exact universes and
remaining limitations. [Incoming integration](incoming.md) distinguishes
proved dependencies, recovered side themes, and read-only context.
