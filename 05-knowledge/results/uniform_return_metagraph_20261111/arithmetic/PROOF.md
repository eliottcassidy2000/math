# Extremal aliquot sources, factor layers, and certified quadratic primes

**PROVED elementary scoped deductions; CITED Bertrand theorem; FINITE-EXACT
enumeration and certificate checks. Chan–Pak and Landau remain OPEN.**
No literature-priority, new-OEIS-term, or enumeration-record claim is made.
Here A(n) denotes A135244, not a tournament address: for n>=2 it is the
largest m with s(m)=sigma(m)-m=n, or0 if the fibre is empty.

## 1. Inheritance and the source information being retained

The [Session66 inverse report](../../session66_root27_aliquot_graph_interfaces_20261026/session66-aliquot/REPORT.md)
already supplies the sharp global bound m<=(n-1)^2, the prime-square
exception, a complete O(n)-candidate inverse, the even-source/square
parity split, and Booker's two-prime hyperbola. The prime/composite odd
square split is already in the cited Pomerance–Yang algorithm. These are
inherited, not discoveries of this packet. Booker's faster even-target
algorithm is published in [*Finite connected components of the aliquot
graph*](https://doi.org/10.1090/mcom/3299); our maximum-only prefix method
does not claim to supersede that individual-fibre algorithm.

The [S81 observer theorem](../../session81_counterexample_families_20261028/ALIQUOT_OBSERVERS.md)
already classifies every prime-power multiplication observer and gives
complete valuation-profile inversion for an arbitrary fixed multiplier.
Its multiplier40 collision is the canonical hostile to scalar compression.
The [aliquot124 correction](../../aliquot_power_124_20261026/aliquot/README.md)
keeps the distinct-Goldbach antecedent open and separates the unique
prime-square-to-perfect target6 from arbitrary-square assertions.

The [previous near-square note](../../chanpak_five_core_bridges_20261110/history/NEAR_SQUARE_AND_ADJACENT.md)
already identifies101 as the first one-apex failure and proves that the C3
dyadic slice6*2^m+5 contains only17,101,197 of the required square-plus-one
form. The [Fortunate/Landau note](../../fortunate_landau_bridge_20261110/README.md)
already proves that every fixed input residue class contains infinitely
many composite square-plus-one values. Neither route is revived here as
an all-height prime or tournament construction.

Concept board: exact inverse fibre / number of prime factors / balanced
extremal sources / guarded observer composition / local quadratic roots /
actual tournament chart. Anchor is source-return information, the niche
is maximum-source enumeration, and the wildcard is extremal equality as
an exact interface between two OEIS sequences. The least-used coordinate
is Omega(m), counting prime factors with multiplicity. The new practical
move is to certify the largest source without enumerating its whole fibre.

## 2. Every factor-count layer has an explicit source bound

Let m>1 have t=Omega(m), and write its prime factors with repetition as
p1<=...<=pt. For j=0,...,t-1, let d_j be the product of the j largest
factors (d_0=1). These are distinct proper divisors. Their geometric means
give d_j>=m^(j/t). Therefore

    s(m) >= F_t(m) := sum_{j=0}^{t-1} m^(j/t)
                   = (m-1)/(m^(1/t)-1).                 (1)

Equality holds exactly when m is a prime power p^t: equality in the
geometric-mean comparisons forces all p_i equal, and that case has exactly
the displayed proper divisors. For fixed m>1, F_t(m) strictly increases
with positive integer t, since its denominator strictly decreases.

Consequently, for every k>=2 and every source with Omega(m)>=k,

    n=s(m) >= F_k(m) > m^((k-1)/k),
    m <= (n-1)^(k/(k-1)).                              (2)

The second bound uses just the largest displayed proper divisor and1;
the first retains the entire chain and is sharper. Equality in n=F_k(m)
holds exactly for m=p^k. There is no claim of equality in every rounded
integer upper bound. Also n>=2^t-1, so t<=floor(log2(n+1)). Thus every
finite target fibre splits into finitely many factor-count layers, with
explicit shrinking bounds. This classifies their cost boundaries, not the
primality guards or every image of every exponent pattern.

For k=3 the sharp real bound can be inverted without numerical roots.
If r^2+r+1=n, then r=(sqrt(4n-3)-1)/2 and
r^3=(n*sqrt(4n-3)-3n+2)/2. Hence

    Omega(m)>=3 and s(m)=n  =>  m<=B3(n),
    B3(n)=(isqrt(n^2*(4n-3))-3n+2)//2.                 (3)

The floor is exact: flooring the inner square root before division by2
does not change the final floor. No floating-point arithmetic is used.
The prime-cube family gives the real equality boundary, for example
8->7,27->13,125->31. The n>=2 domain excludes the infinite prime fibre F(1).

## 3. Maximum-source certificates and an infinite classified family

The inherited sharp k2 bound gives

    A(n)=(n-1)^2  iff  n-1 is prime.                   (4)

If n-1 is not prime, every Omega2 source is pq with distinct primes and
p+q=n-1. For fixed sum, pq is largest when p and q are closest. Let G(n)
be this largest product, or0 if no distinct prime pair exists. Then

    G(n)>B3(n)  =>  A(n)=G(n).                        (5)

Proof: (4) excludes the prime-square branch, G exhausts distinct
semiprimes, and (3) excludes every other composite source. This is a
checkable sufficient certificate, not an assertion that G(n)>B3(n) for
all odd n. At n15, G=33 but B3=35; the sufficient certificate fails while
the true maximum is still33. At n7 there is no distinct pair and A(7)=8.

**Consecutive-prime theorem.** For every consecutive pair of odd primes
p<q,

    A(p+q+1)=pq.                                     (6)

Proof. Bertrand's theorem gives q<2p. This is the sole non-elementary
import here, in the elementary classical sense; see Erdős' original
[1932 proof, pp194–198](https://users.renyi.hu/~p_erdos/1932-01.pdf).
Put z=(pq)^(1/3)>2. Since1<q/p<2,

    p+q <= (3/sqrt(2))*sqrt(pq)
          < z^2+z,

where the strict inequality follows by squaring positive terms and using
(z-2)(z-1/2)>0. Thus F3(pq)>p+q+1, so every Omega>=3 source for that
target is strictly smaller than pq. The odd target cannot receive a prime
square: its required prime base p+q is even and greater than2. Finally,
any distinct-prime pair of sum p+q with smaller entry greater than p
would lie between the consecutive primes p and q. Thus p,q are the closest
pair and maximize the semiprime source. This proves (6). There are
infinitely many such pairs, so (6) is an unconditional infinite type.

The initial pairs produce

    A(9)=15, A(13)=35, A(19)=77, A(25)=143,
    A(31)=221, A(37)=323, A(43)=437.

The statement is an elementary deduction with a proof here, without a
priority claim. It does not classify every odd target, prove a Goldbach
statement, or constrain forward aliquot trajectories.

## 4. A complete adaptive prefix algorithm

The implementation uses the following exhaustive disjoint source shapes.

1. Prime squares are settled immediately by (4).
2. For even targets, every even source m>=6 satisfies m<=2n-6, because
   1,2,m/2 are distinct proper divisors. An odd source must be an odd
   square a^2. Prime a is already handled. For composite a, its least
   prime factor is at most sqrt(a), so n>=a^(3/2) and a<=n^(2/3).
3. For odd targets, find the closest distinct-prime pair. If (5) holds,
   its product is certified globally. Otherwise enumerate the remaining
   higher-factor sources through B3(n). All Omega2 sources were already
   considered. This fallback is always finite and complete; the program
   never assumes a Goldbach conjecture.

For the complete prefix2<=n<=1,000,000 the only fallback targets are5,7,15.
The resulting type counts are:

| Selected maximum type | Number of targets |
|---|---:|
| Empty fibre |150232|
| Prime square |78498|
| Closest distinct semiprime |499996|
| Even source, after the prior types |270688|
| Composite odd square |584|
| Residual exhaustive source |1|

The final row is7<-8. These are **FINITE-EXACT** categories determined by
the stated branch order; the closest-semiprime row includes15, whose
global certification needs the fallback. No asymptotic category density
is inferred. The computed normalized table SHA256 is
`c784facbd3430155430704b3882109b36f4dd55ee2010b4546c60cba2e54a033`.

`audit_maxima.py` independently uses an additive divisor sieve through
2,000,000, trial factorization of composite square roots, prime-list
predecessor search, and an equivalent squared rational inequality rather
than the integer-cap code. It agrees on all999,999 entries. A separate
complete source sieve through249001 audits targets2..500 with no source-
shape pruning. The full published OEIS table through10000 also agrees.
The algorithm improves this workspace's maximum-source enumeration;
no claim is made that the values or method are globally new or fastest.

## 5. An exact A002496/A135244 extremal interface

For every integer r>=1, (4) immediately gives

    r^2+1 is prime  iff  A(r^2+2)=(r^2+1)^2.           (7)

This is a genuine predicate-preserving map. It is reversible by testing
the sharp extremal equality, and the quartic source has aliquot target
r^2+2 exactly under that prime guard. It supplies neither prime infinitude
nor a faster primality algorithm. In the shared numerical domain r<=999,
all112 prime cases agree with the independently audited maximum table.

For A002496 directly, apart from r1 giving2, r must be even. An odd prime
ell dividing r^2+1 satisfies ell=1 modulo4. Such a prime has exactly two
roots of-1 modulo ell. If a is a quadratic nonresidue, a^((ell-1)/4) is
one root. The sieve enumerates these two progressions, retaining the
exception r^2+1=ell itself. For r<=R every composite value has a prime
factor at most R, because floor(sqrt(R^2+1))=R. Thus sieving all eligible
primes through R is complete. These elementary mechanisms are classical.

Every retained value P=r^2+1 also receives a separate Lucas certificate:
the complete factorization P-1=r^2 is known from the small integer r, and
the program finds a with

    a^(P-1)=1 (mod P),
    gcd(a^((P-1)/q)-1,P)=1 for every prime q|r.         (8)

To see why (8) proves primality, let ell be any prime divisor of P.
The order of a modulo ell divides P-1, and every prime-power part of P-1
must divide that order by the gcd tests. Therefore P-1 divides ell-1,
forcing ell>=P and hence ell=P. This short proof is the classical full-
factorization Lucas criterion; it is not a new primality theorem.

At R=1,000,000 the complete enumeration has54,110 primes, with54,109 Lucas
certificates beyond2 and445,891 proper-divisor certificates for rejected
even-root candidates. The largest successful Lucas base is83, at
P476985601. Independent literal trial division agrees for r<=10000, and
all10000 terms in the pinned OEIS b-file agree. This is reproducible
certified computation, not a new record or evidence sufficient for Landau.

## 6. Actual near-square tournament witnesses and a boundary

Intersecting the inherited five-core support receipt through2^24-1 with
r^2+1 shows that none of its322 first-three-core holes is a positive even
square plus1. Thus all2047 targets r^2+1, r=2,4,...,4094, already use the
K1,C3,S4 bank. This is an inherited-data consequence, not a new coverage
range. The fresh complete inverse records a witness minimizing min(A,B)
in that normalized grammar, using the nonnegative seven-coordinate law.

There are1974 axis witnesses. The axis-only proposal first fails at
2501=50^2+1; an actual C3 chart has cut1,A1,B157. The minimum shorter
flank is29 at8820901=2970^2+1, which in this bank is supplied only by C3:
cut1,A29,B24978. Other high requirements are18 at4857617 and17 at9884737.
Thus the observed support does not come from a uniform axis or a bounded
short flank below29. No universal fixed bound is refuted by this finite
data, and all-height square-plus-one coverage remains open.

The search visits672,392,310 parameter pairs, checks all2047 saved
witnesses with the inherited signed compiler, and verifies32 literal
graphs of order at most13 by Hamilton-path DP. Its source-size minimization
is subordinate to minimum short flank; no unrestricted minimum-order
tournament claim is made. The exact table is `near_square_charts.tsv`.

## 7. Composition and what the maximum selector loses

The inherited exact response for gcd(a,b)=1 is

    s(ab)=b*s(a)+sigma(a)*s(b),  sigma(ab)=sigma(a)*sigma(b).

In coordinates(s(a),sigma(a)), coprime multiplication by b is the matrix
[[b,s(b)],[0,sigma(b)]]. It is bI+s(b)E for the idempotent
E=[[0,1],[0,1]], so these guarded operators commute and diagonalize in
coordinates(a,sigma(a)). Repeated-prime factors require the valuation
profiles retained in the S81 observer theorem. The guard cannot be erased.

Maximum selection R(m)=A(s(m)) is defined for composite m, is inflationary
and idempotent within each finite fibre, but does not make its scalar
quotient a composition law. For instance s(6)=s(25)=6 while s(12)=16 and
s(50)=43. Source return to25 changes the response under multiplication
by2. This is the exact structural lesson for the tournament reset: a
canonical source is useful only together with its allowed next operations.
The comparison is observer descent, not a functor between the two subjects.

| Map | Preserved predicate | Lost data / necessary guard |
|---|---|---|
| Aliquot source -> factor layer | Exact finite upper bounds | Actual primes and their exponent pattern |
| Fibre -> largest source | The target and a canonical extremum | Other sources' future responses |
| Consecutive odd primes -> aliquot edge | Global maximum pq at p+q+1 | Does not cover arbitrary odd targets |
| Quadratic prime -> quartic source | Equivalence(7) | Prime guard remains the whole infinitude problem |
| Even-square target -> core chart | Actual Hamilton count in the finite grammar | Chart, addresses, and marked source must remain |

The next constructive question is whether the increasing short-flank
requirements admit a controlled descent returning an actual usable source.
The arithmetic layers suggest a certificate-first search architecture;
their bounds do not supply that missing tournament operation.
