# Divisor fibres and the distinct roles of 7 15 20 24 25 26 and 33

**PROVED elementary collision formulas; inherited classifications credited;
FINITE-EXACT graph and group controls.** Shared integers are linked only
where the map and its preserved data are specified.

## Twenty and the binary divisor row

The prior [divisor-overlap theorem P1--P3](../aliquot_divisor_overlap_20261028/THEORY.md)
already settles the user's addition rule. With s(n)=sigma(n)-n,

    there exist positive integers a,b with ab=n and s(n)=1+a+b
        iff tau(n)=4 or s(n)=n+2.                        (1)

The two alternatives are disjoint. Here tau counts ALL positive divisors.
Exactly five proper divisors means tau(n)=6; under that condition (1)
holds only for n=20. Its actual divisor sum1+2+4+5+10=22 also equals
1+1+20, with (a,b)=(1,20). The repeated unit and the non-proper divisor20
are a numerical compression witness, not the actual proper-divisor list.

Thus20 really participates in both descriptions. Sixteen is a prime power
with four proper divisors and s(16)=15, but it fails (1). The categories
cannot be made exclusive by including16 in both sides without changing
the predicate. They should retain source, exponent shape and sum separately.

For an odd prime p,

    sigma(4p)=(1+2+4)(1+p)=7(p+1),
    s(4p)=3p+7,       s(4p)-4p=7-p.                      (2)

This inherited chart contains12->16,20->22,28->28, with excess4,2,0.
The coefficient7 is literally the binary divisor row1+2+4.

## Complete prime-guarded collisions between dyadic rows

Fix integers a>b>=1 and odd primes p,q. Put

    A=2^(a+1)-1, B=2^(b+1)-1, g=gcd(A,B), R=A/g, S=B/g.

Then

    sigma(2^a p)=sigma(2^b q)
       iff p=2S ell-1 and q=2R ell-1 for an integer ell>=1.  (3)

Both displayed outputs must remain prime. Their common divisor total is
2ell*lcm(A,B). This is a complete classification within the specified
shapes, with no assertion of infinitely many simultaneous prime values.

Proof: multiplicativity gives A(p+1)=B(q+1). Since R,S are coprime,
p+1=St and q+1=Rt. Both R and S are odd, while p,q are odd, so t=2ell.
The reverse implication follows by substitution. The familiar Euclidean
algorithm on powers gives g=2^gcd(a+1,b+1)-1. Thus every pair of dyadic
divisor rows has an explicit compatible chart, including its prime guards.

At a=2,b=1 this becomes

    p=6ell-1, q=14ell-1,
    n=4p=24ell-4, m=2q=28ell-2,
    sigma(n)=sigma(m)=42ell.                             (4)

The first pair is20 and26:

    20+22=42=26+16.

The next admissible pairs include68,82 and164,194. Furthermore,

    s(n)-n=8-6ell,      s(m)-m=4-14ell.

Consequently20/26 is the UNIQUE pair in this entire prime-guarded family
whose first member is abundant; the second member is always deficient.
This gives a real connection among20,26 and the1-2-4 divisor row. Equal
sigma does not preserve source, proper-divisor count, abundance, or the
next aliquot state. For example tau(20)=6 whereas tau(26)=4.

Incoming [Session141's aliquot classification](../session141_metric_aliquot_square_paths_20261102/aliquot_twenty/PROOFS.md)
also proves the complete fibre sigma^-1(42)={20,26,41}. Our chart specifies
two dyadic-prime source shapes; it does not claim that20 and26 exhaust the
fibre. That packet extends the excess-two classification through tau17.

## Why 16 and 33 can belong together under a different observer

For an odd prime p other than3,

    s(3p)=p+4.

Therefore

    s(2^a)=s(3p) iff p=2^a-5 is prime and p!=3.            (5)

The first valid exponent is a=4, giving16->15<-33. Later examples include
64->63<-177 and256->255<-753. The excluded a=3 would use p=3, but9 has
s(9)=4, not7: the coprimality guard cannot be dropped.

Equation (5) preserves the aliquot TARGET. Equation (3) preserves the total
DIVISOR SUM. They are different quotient maps. For instance sigma(16)=31
and sigma(33)=48, despite their shared aliquot target. These two examples
explain how seemingly conflicting lists can record different observations.

## Binary levels are a weighted lift of a divisor chain

The proper divisors of2^r form the r-element chain1,2,...,2^(r-1).
Replace the divisor2^j by all binary words of length j. This expands its
weight into2^j separate addresses, with heap label2^j+value(word). Hence

    s(2^r)=2^r-1 = number of binary words of length less than r.

At r=3 this is1+2+4=7; at r=4 it is1+2+4+8=15. The divisor chain itself
does NOT have7 or15 vertices. The lifted word address is the extra coordinate
needed for tree incidence and the geometric embedding in GAP_GEOMETRY.
Truncating a word gives tree ancestry; it is not ordinary Collatz iteration.

The number24 supplies a related two-factor example:
sigma(24)=(1+2+4+8)(1+3)=15*4=60, so s(24)=36. This retains the full
2^3*3 exponent rectangle, not just the count or the total of its cells.

## The square25 count anatomy

For the fixed square-sum guard on1..25, there are32 undirected edges and20
directed Hamilton paths, each prescribing24 arcs. There are four compatible
pairs, each prescribing26 union arcs, and no compatible triples. These are
the inherited [joint kernel](../joint_incidence_clocked_frames_20261102/PATH_KERNEL.md)
data, rechecked here. Consequently the guard-orientation generating function is

    2^32 +20*2^8(z-1)+4*2^6(z-1)^2
      =4294962432+4608z+256z^2.                           (6)

The remaining268 tournament edges are outside the guard, giving the user's
factor2^268. Thus25 is the vertex count,20 the path count,24 the length of
each path, and26 the union size of a compatible pair. They are tied by
literal shared-edge incidence, not by a guessed arithmetic isomorphism.

The [minimal complement carrier](../lower_frontier_complement_fibonacci_20261102/COMPLEMENT_CARRIER.md)
has12 isolated path types and four compatible pairs. A pair contributes
five closed response types: both paths, two leaves which return to both,
and two leaves which stay single. With one dead state,

    33=12+4*5+1;     24=20+4 initial active legality types.

One refinement gives32 active types. This33 counts observer states; it is
not an extra square-sum vertex or a physical honeycomb site. G25 is rigid,
so keeping the guard graph canonically recovers the named command frame.
Discarding the guard, adding arbitrary tournament paths or adding arithmetic
extension contexts changes the minimum carrier.

The incoming [integer-content refinement](../session141_metric_aliquot_square_paths_20261102/square_paths_content/PROOF.md)
shows that the full polynomial has exact coefficient gcd2^276, so (6) can
also be written2^276(16777197+18z+z^2). Integer translation preserves the
coefficient gcd; reversal and joint multiplicities explain its binary factors.
Its order26 extension is2^296(536870561+320z+29z^2+2z^3). These are incoming
results, not new enumerations here. The [common-block connection](../session141_metric_aliquot_square_paths_20261102/CONNECTIONS.md)
further records that reversing the25 pair changes its fixed-arrow successor
deletion cost from1 to2 while preserving its undirected cost1. A retained
terminal is thus a concrete additional coordinate, even at unchanged counts.

## The older geometric 15 20 24 and multiplicities 1 2 4

The corrected [cardinality atlas](../cardinality_bridges_20261018m.md)
recovers the June icosahedral work while preserving MISTAKE-059's warning:
a quotient of cardinalities does not prove uniform fibres. On five labels,
direct permutation enumeration gives

    A5: 1 identity, 15 elements of order2,
                   20 of order3, 24 of order5.

Its rotation-axis interpretation has15 axes with one nonidentity half-turn,
10 axes with two nonidentity third-turns, and6 axes with four nonidentity
fifth-turns. Thus15=15*1,20=10*2,24=6*4. Those multiplicities have their own
precise origin in cyclic stabilizers. They do not identify the20 rotations
with the20 square-sum paths.

The geometry also distinguishes the angle-defect test from the SIDE triangle
inequality: the reciprocal sums for(2,3,5),(2,3,6),(2,3,7) are respectively
above, equal to and below1. These are spherical, flat and hyperbolic triangle
types. All corresponding distance metrics still obey the ordinary triangle
inequality. Neither those group orders nor the carrier33 is thereby a proof
of Collatz convergence or a new Hadwiger-Nelson lower bound.
