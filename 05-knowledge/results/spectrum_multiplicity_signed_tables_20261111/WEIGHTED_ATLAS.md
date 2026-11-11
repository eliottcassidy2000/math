# Multiplicity is a weight on the signed multiplication-table atlas

Date: 2026-11-11. **PROVED:** the product law, injections, primitive-fibre
second differences, and growth interpretation below. **INHERITED:** the
strong-component factorization, exclusion of 7 and 21, the global two-class
inventory at H=63, and the two-prime-lane reduction. **OPEN:** unrestricted
odd-spectrum completeness.

## The operation behind 2, 5, 16

Use the signed source coordinates of
[PUNCTURED_ATLAS.md](signed_arcs/PUNCTURED_ATLAS.md): an odd modulus is
`m=4x+sigma`, where sigma is +1 or -1. For m>=3 this determines x>=1
uniquely. Treat m=1 as a separate multiplicative identity, with formal
coordinate (0,+1). A coordinate labels an edge of the signed atlas; it is
not the assertion that x itself is a Hamiltonian path count.

The ordinal sum of tournaments multiplies their Hamiltonian path totals.
Transporting that operation to these coordinates gives exactly

    (x,sigma) * (y,tau) = (4xy + tau*x + sigma*y, sigma*tau).       (1)

Indeed, `(4x+sigma)(4y+tau)` equals four times the displayed first
coordinate plus the displayed sign. Thus (1) is associative and commutative
on count labels, although ordinal sum of the underlying tournaments is
generally noncommutative. The lost information is the order of the strong
components, which must be retained when counting actual tournaments.

The directed triangle has H=3 and coordinate (1,-1). Its multiplication
action is therefore

    J(x,sigma) = (3x+sigma,-sigma).                            (2)

Starting with the forbidden label 7 gives

    labels:       7        21        63        189        567
    sources:    (2,-) -> (5,+) -> (16,-) -> (47,+) -> (142,-).

The first two coordinates label precisely the missing raw arrows 2->1
on the minus sheet and 5->16 on the plus sheet. The third labels the
minus halving arrow 16->8, whose modulus is the attained count 63.
This distinguishes the label operation J from a raw Collatz operation:
at even x, J does not halve x. In closed form, after a multiplications,

    m_a = 7*3^a,
    sigma_a = -(-1)^a,
    x_a = (7*3^a + (-1)^a)/4.                               (3)

This is an exact map linking the user's three numbers, rather than an
identification of their two different graph roles.

## Give each table edge its full multiplicity sequence

Let N_n(h) count ordinary tournament isomorphism classes of order n with
H=h, without identifying a class with its converse. Let L_n(h) count
tournaments on the fixed labelled set [n], and put P_n(h)=L_n(h)/n!.
Assign to the edge with modulus h the two weight sequences
`(N_n(h))_n` and `(P_n(h))_n`. Its support is exactly whether some term
is positive. The complete formulas and finite data are in
[multiplicity/PROOF.md](multiplicity/PROOF.md).

Appending a singleton sink is injective on isomorphism classes: recover
the old tournament by removing the final singleton in its unique ordered
strong-component decomposition. Hence every edge weight N_n(h) is
nondecreasing with n, and support through order n equals support at order n.
The analogous statement holds for normalized labelled mass P_n(h).

Appending a triangle is also injective, by removing the final strong
component. Consequently

    N_(n+3)(3h) >= N_n(h),
    P_(n+3)(3h) >= P_n(h)/3,
    L_(n+3)(3h) >= (n+3)! L_n(h)/(3 n!).                    (4)

The factor 3 is the triangle's automorphism order. Thus (2) transports
positive edge weight in a realizable direction, with a precise order and
symmetry cost. Its converse is false: the positive edge 63 has no positive
predecessor at 21. Removing a factor of three from a count is not the same
as removing a terminal triangle from an actual tournament.

Write a_r(q) for the number of nontrivial strong classes of order r and
count q. Splitting according to the final strong component gives the exact
growth identity

    N_n(h)-N_(n-1)(h)
      = sum_(q|h,q>1) sum_(r>=3) a_r(q) N_(n-r)(h/q).       (5)

Equation (5) holds for n>=1. Set N_0(1)=P_0(1)=1, both order-zero
counts at h!=1 to zero, and all negative-order counts to zero. The removed part on
the left consists exactly of tournaments ending in a singleton. A final
nontrivial component gives each term on the right uniquely. Replacing
a_r(q) by the sum of reciprocal automorphism orders gives the same
identity for P. Equation (5) describes how existing support grows; it
cannot create a new h without a strong atom supplying a new factor.

## The two holes leave an exact growth deficit

The inherited global H=63 inventory consists of two rigid strong classes,
both of order 8. Since every proper nonunit factorization of 63 contains
7 or 21, every H=63 tournament has exactly one of these components and
otherwise singleton components. Therefore, for n>=7,

    N_n(63) = P_n(63) = 2(n-7),
    sum_(j=1)^n N_j(63) = (n-7)(n-6).                     (6)

The source (16,-) is thus born at order 8 with weight 2 and grows linearly.
The preceding sources (2,-) and (5,+) have identically zero weights.

More generally, for fixed a>=2, the maximum-component words for
H=7*3^a consist of one of those two H=63 cores and a-2 triangles. Their
length is a-1; there are 2(a-1) ordered words, each of size 3a+2 and
automorphism order 3^(a-2). Their exact contribution is

    2(a-1) binom(n-2a-3,a-1),

with the term interpreted as zero before size 3a+2. Lower-component
words contribute lower-degree terms. In particular,

    N_n(7*3^a) ~ 2 n^(a-1)/(a-2)!,
    P_n(7*3^a) ~ 2 n^(a-1)/(3^(a-2)(a-2)!).              (7)

For comparison, the ray 3^a has ordinary degree a and leading coefficient
1/a!. The two absent initial labels have a measurable effect on the
maximum number of freely separated strong components. The ordinary
asymptotic is inherited from the fixed-H multiplicity packet; the weighted
version retains the automorphism sidecar.

## Recovering the indivisible strong inventory from growth

Call h>1 product-primitive relative to the attainable spectrum when it
cannot be written uv with attained u,v>1. For such h every tournament
has exactly one nontrivial strong component. Hence

    N_n(h) = sum_(r<=n) a_r(h)(n-r+1),
    Delta N_n(h) = sum_(r<=n) a_r(h),
    Delta^2 N_n(h) = a_n(h),                              (8)

where Delta f_n=f_n-f_(n-1), and values at negative orders vanish.
The same formulas for P recover automorphism-weighted strong counts.
These are exact finite differences at every order, not just asymptotics.
In particular N_n(h) is a nondecreasing discrete convex function, and its
eventual slope is the total number of strong classes with H=h.

The exclusion of 7 and 21 makes every attained member of the following
set product-primitive:

    p (odd prime p!=7),   7p (odd prime p!=3),   63,   343. (9)

Here 49 is already included by p=7 in the second lane. Proper
factorizations of 63 use 7 or 21, and those of 343 use 7. The inherited
THM-4094,
[`hamiltonian-matching-deficit-and-two-prime-lane-completeness`](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md),
shows that realizing the two infinite lanes in (9), together with the
known carry atoms 49,63,343, suffices for the entire desired spectrum.
Equation (8) translates that remaining existence problem into a precise
growth question: each required lane must have at least one positive
second difference. Polynomial growth on already attained lanes does not
imply positivity on a different lane.

For any proposed h the inherited strong-order bound limits its atom
search to r<=floor(1+3 log_5 h). Thus each individual positivity question
is finite, although a uniform positive answer over the infinite lanes
remains open. This is the manageable finite object for a fixed modulus.

## Connection contract and useful boundary

Source: an actual tournament and its ordered SCC word. Target: a signed
table edge with its multiplicity sequence. Map: H, then m=4x+sigma.
Preserved: count, support, multiplication, size filtration, and symmetry
weight when retained. Lost: the directed edges, endpoint response, and
which strong component supplied a factor. Needed sidecar for a birth
construction: an actual marked source and its endpoint/cover response.
Decisive test: 63 is positive although 21 is zero, so the inverse of J
cannot be inferred from scalar divisibility.

The signed atlas now has both a support problem and an exact measure on
each supported edge. The finite controller in PUNCTURED_ATLAS.md concerns
routes on a candidate support; (5) concerns growth given strong atoms;
the new return construction concerns an actual marked tournament.
Keeping these three maps explicit is what allows their positive results
to be combined without silently turning a route into a construction.
