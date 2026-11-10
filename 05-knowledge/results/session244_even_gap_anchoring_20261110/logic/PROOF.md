# Scalar abundance versus anchored construction

Status: **PROVED abstract countermodel**; arithmetic atom reduction **INHERITED**. Physical
inverse-tree invariant: **CONDITIONAL / OPEN**. Chan–Pak remains **OPEN**.
No claim of literature novelty. This packet does not construct a tournament
with an arbitrary prescribed full Hamilton-path count.

## Inheritance and the live concept board

The closest physical mechanism is Session199
`session199/physical-interpolation/PROOF.md`, §4: every even difference occurs
between infinitely many strong single-arc pairs. Its §1 matching obstruction
and §5 target-preserving-but-not-arbitrarily-placed rebases are the relevant
hostiles and near misses. Session221
`05-knowledge/results/session221_computation_to_construction_20261107/physical/PROOF.md`
provides the explicit half-gap selector and its uncontrolled digital baseline.
Session243
`05-knowledge/results/session243_prime_gap_bridges_20261110/constructions/GAP_UNIVERSALITY.md`
provides the compact positive pair with upper count at most 5N. The least-used
sidecar here is the endpoint coordinate, as distinguished from the full count
in THM-5001, `01-canon/theorems/THM-5001-every-positive-endpoint-quotient-has-a-small-strong-realizer.md`.
Arbitrarily long physical windows are inherited from THM-5000,
`01-canon/theorems/THM-5000-coprime-interface-composition-consecutive-odd-windows.md`.

Board: compact differences; scalar multiplicative closure; endpoint fibres;
full-count prime atoms; target-decreasing physical invariants. The first three
forget the target location; the last two retain exactly the information that
needs to be supplied. No subtraction or division closure is inherited.

## 1. A strong simultaneous scalar countermodel

Let

    S = {positive odd integers} \ ({7,21} ∪ {primes p >= 151}).

Then all of the following hold simultaneously.

1. S contains 1 and is multiplicatively closed.
2. Every positive even difference 2N has witnesses a,b in S with
   2N+1 <= a < b=a+2N <= 5N.
3. Every even difference occurs infinitely often in S.
4. S contains arbitrarily long intervals of consecutive odd integers.
5. S agrees with the conjectured spectrum below 151, yet contains no prime
   at least 151. Its relative density among odd integers is 1, by the standard
   elementary estimate pi(x)=o(x).

**Multiplication.** A product of two members cannot be 7, since 7 is prime.
It cannot be 21, whose only nonunit proper factorization is 3·7. If both
factors exceed 1 the product is composite and therefore cannot be an excluded
large prime. A unit factor leaves the other member unchanged.

**Compact pairs, N>=30.** Solve

    a = 1 (mod 2), a = 0 (mod 3), a = -2N (mod 5).

CRT gives one residue class modulo 30. The N consecutive integers in
[2N+1,3N] contain a representative. This a is an odd multiple of 3 greater
than 21. Its partner b=a+2N is an odd multiple of 5 greater than 5, with
b<=5N. Both are composite members of S. This is an actual bound on absolute
placement, not just existence of a residue class somewhere.

**Compact pairs, 1<=N<30.** Usually use a=2N+1, b=4N+1. The only failures
are N=3,5,10, resolved respectively by (9,15), (15,25), (23,43). These pairs
obey the same bounds. All entries are below 151. The range condition
2N+1<=3N holds even at N=1; the pair there is (3,5).

**Infinitely many pairs.** For any fixed N, the same CRT progression contains
arbitrarily large a. Beyond 21 both a and a+2N are appropriate odd composites.
Negative even differences follow by exchanging the endpoints; difference
zero is immediate.

**Long odd intervals.** For every L>=1 set M=(2L+3)!. The L successive odd
integers M+3,M+5,...,M+(2L+1) are all composite: the ith is divisible by
2i+3 and strictly larger than that divisor. They exceed 21.

The cutoff 151 can be replaced by any larger cutoff. Thus *any fixed finite
collection* of correct affirmative full-count evidence can be included while
still excluding every sufficiently large prime. This is stronger than
removing one isolated target from a scalar model.

This is a logical countermodel to the listed **scalar properties**, not a
counterexample to Chan–Pak, and not a model of all tournament identities.
In particular, exact graph relations, source recurrences, and full physical
state realizability have not been asserted for S. One cannot conclude that
those stronger data are useless from this countermodel.

### Uniform approximate anchoring still fails

For every N>=1 and every odd target anchor m>=51, the very same model has a
pair a,a+2N in S with

    m-28 <= a <= m.

Indeed the 15 odd integers in that window contain every odd residue modulo
30 once, including the CRT class above; a>=23 makes both endpoints proper
composites. Likewise every odd m>=23 has such a pair with m<=a<=m+28.
Thus the lower endpoints of every prescribed gap form a uniformly dense set,
with additive placement error at most 28, independently of both anchor and
gap. In particular every fixed relative-error bound holds at sufficiently
large anchors. Yet all sufficiently large prime targets remain absent.
Quantitative approximation must be coupled to an **exact physical repair**;
no asymptotic improvement that retains a nonzero error is automatically such
a repair. This statement does not claim that every smaller constant error
also fails, or that actual tournament gap endpoints obey this bound.

## 2. Why endpoint universality alone does not repair it

An additional scalar coordinate E ranging over every positive integer places
no restriction on the projection S without a coupling theorem. Even the
basic bounds 0<E<D and U>=D can be respected: for each N>=1 attach a formal
state

    (D,E,U,d) = (12N+69, N, 24N+138, 6N+33).

Both D and d are odd composites above 21 and hence lie in S. The basic physical
inequalities E<=d<=D, U>=D, U>=2d and oddness of d also hold. Its endpoint list
can formally be (N,D-N), with nonnegative entries summing to D. This gives
every endpoint N, including 7 and 21, while the full-count projection still
misses every large prime. These are formal tuples, **not physical tournament
states**; the tuple does not claim closure under O/Z/J or the full endpoint
identities. The precise logical point is that a marginal-surjectivity theorem
for E cannot supply arbitrary fibres (D,E) or a prescribed D.

Consequently, taking gcds of scalar increments, Bezout combinations of their
values, or subtracting one endpoint from another is not a physical constructor.
The countermodel makes those logical failures simultaneous with compactness,
multiplication, and very abundant values.

## 3. Inherited exact multiplicative atoms of the desired spectrum

**Inheritance correction.** This reduction is already proved and independently
audited in THM-4094,
`01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`,
§§5–7. That theorem also supplies the finite strong carry realizers 49,63,343.
The argument below is an explicit atom restatement and self-contained proof,
not a new completeness reduction. THM-4094 allows p=7 in its 7p lane, thereby
including 49 there; we separate 49 into the finite list instead.

Let A be all positive odd integers except 7 and 21. Its nonunit multiplicative
atoms are exactly

    odd primes p != 7;
    7p, where p is an odd prime other than 3 and 7;
    49, 63, 343.                                             (1)

Here an atom is an element not expressible as a product of two nonunit
elements of A. This is an arithmetic theorem about A, not a claim that those
numbers are currently all known tournament counts.

**Proof of irreducibility.** Primes are immediate. The only proper split of
7p has a forbidden factor 7. Proper splits of 49 use 7; those of 63 use
7·9 or 3·21; those of 343 use 7·49. Thus all listed values are atoms.

**Proof of generation.** Write an arbitrary allowed odd n as 7^a u with
7 not dividing u. If a=0, factor u into odd primes other than 7. If a>=2,
write a=2r+3s with r,s>=0 (every integer at least 2 has this form), and use
49^r343^s together with the prime factors of u. If a=1 and u has a prime
factor p other than 3, take the factor 7p and factor u/p. In the remaining
case u=3^b, allowedness forces b>=2, so n=63·3^(b-2).
This proves generation and completeness of the atom list.

For the actual tournament spectrum, multiplicative closure follows by an
ordered sum of two tournaments: all cross arcs go from the first to the
second, and each Hamilton path is a Hamilton path of the first followed by
one of the second. Thus the path counts multiply exactly. The output need
not be strong, which is immaterial to the unrestricted Chan–Pak target.

**Constructive payoff.** It suffices to realize the prime family p!=7, the
family 7p for p!=3,7, and the three finite atoms. Establishing just the prime
family does not suffice via multiplication: 7p is an additional infinite
family of atoms. In the countermodel of §1, all those composite atoms are
already present, so even knowing all composite targets would leave the
large-prime obstruction intact. A useful physical number-theoretic objective
is therefore a prime selector plus a guarded 7p lift, not an appeal to prime
factorization with the exceptional prime silently dropped.

## 4. A small anchored inverse tree that would suffice

This offers a precise sufficient **physical invariant to seek**, rather than
claiming that universal gaps already give it. Suppose there is a collection
C of realizable, fully specified rooted tournament states with these features:

- C contains a state with each count in B={1,13,15,41,43}.
- For every state in C of count n and each sign epsilon in {-1,+1}, whenever
  2n+epsilon is positive and not 7 or 21, there is an explicitly realizable
  state in C of count 2n+epsilon.

Then every positive odd value other than 7 and 21 is realized.

**Proof.** For a target m>1 define its unique odd parent

    P(m)=(m+1)/2 if m=1 (mod 4),
    P(m)=(m-1)/2 if m=3 (mod 4).

P(m) is positive odd and smaller than m. The only allowed targets whose
parent is 7 or 21 are 13,15,41,43. Otherwise repeated parent reduction ends
at 1 without leaving A. Reversing this descent gives the required physical
construction from a seed in B. In this fixed parent graph, the five seeds
are minimal: they are precisely the roots of the allowed-target forest.
This is **not** a claim of global minimality among all possible invariants.

The invariant must include every coordinate required to continue the chosen
construction. It need not hold on every tournament with count n; canonical
witnesses forming a stable family would be enough. Conversely, arbitrary
choice of a tournament with count n cannot replace the C-membership premise.
Session199's equal-scalar rebases and Session221's equal-count/different-gap
examples are exactly why this qualification matters.

The conditional rule is stronger than needed (different valid target-parent
maps could work), but it converts the missing anchor into a bounded pair of
exact scalar identities and explicit sidecar-closure obligations. Its binary
inverse terminates in O(log m) steps. No such physical stable family has been
proved here, and no universal O/Z/J implementation is asserted.

### A bounded repair version with concrete physical sidecars

A different sufficient condition can separate coarse selection from exact
repair. Suppose for each large odd target m one can construct a rooted source
s with scalar F(s)=D-E+U in [m-2r,m], together with a certified **same-level
sidecar fibre** admitting all corrections

    F(s_j)=F(s)+2j, 0<=j<=r.

Choose j=(m-F(s))/2. The coarse selector and the r+1 repair certificates
suffice; neither scalar approximate placement nor unanchored repair windows
alone suffices. The source s may depend on m, and the repair must belong to
that very source, not merely to another source of the same scalar count.

For actual O/Z/J states, a retained-root pre-suffix O/Z toggle has profile
change (Session199 §5)

    (Delta D,Delta E,Delta U,Delta d)
       = (U-E-(A+p)d, 0, -B d, 0), B=2^k.

Under F its exact change is

    U-E-(A+p+B)d.

Thus a candidate +2 repair requires the exact pre-toggle level equation
U-E=(A+p+B)d+2, and every further repair requires its own compatible
pre-toggle equation after all preceding operations. For the balanced
readout G=D+2U the equation instead is U-E=(A+p+2B)d+2.
These are actual sidecar equalities to test, rather than a scalar gcd claim.
They are necessary and sufficient for that specific toggle's scalar response,
conditional on the chosen physical source, suffix, and readout. Physical
strongness and retention of a usable source must be checked separately.
Session199's matching theorem already rules out obtaining an arbitrary
repair chain from single-letter changes in one long retained-root block.
A valid bounded repair scheme therefore needs different blocks, rebasing,
or other physical operations. No such complete scheme is proved here.

## 5. Alternative placement invariant for existing interval constructors

Suppose physical windows contain every odd integer in [A_j,B_j], with odd
endpoints, A_1<=23, B_j unbounded, and

    A_(j+1) <= B_j+2,  B_(j+1)>=B_j.

If the known finite values below 23 other than 7,21 are supplied, the windows
cover the entire target. This follows by induction on the right endpoint of
the union. THM-5000's widths and exponentially growing baselines do not
currently establish these overlap inequalities. Compact gap witnesses also
do not supply them. This is a second concrete sufficient anchored statement,
not an unproved assertion that current windows overlap.

## 6. Exact checks and scope

Run `python3 -I -B verify.py` in this directory and compare with `CHECK.json`.
Repeat under `-O`; all checks use explicit exceptions, not assert statements.
The standard-library checker tests compact pairs for N=1,...,10000,
multiplication for pairs of model members below 1000, the atom classification
against direct divisor testing through 19999, and factorial intervals of
lengths 1,...,100. It also checks 47,500 uniformly anchored CRT pairs and
the conditional inverse parent tree for targets below 100,000. The infinite results are proved above, not inferred from
these finite tests. No literature priority or proof of Chan–Pak is claimed.
