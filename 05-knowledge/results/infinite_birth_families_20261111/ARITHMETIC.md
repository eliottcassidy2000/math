# Prime-factor strata, inverse divisor sums, and quadratic tournament rays

**PROVED elementary bounds and identities; FINITE-EXACT enumeration;
INHERITED seed construction; OPEN infinitude and universal seed questions.**
No algorithmic optimality or literature-priority claim.

## 1. Keep the three sequences' objects distinct

[A000568](https://oeis.org/A000568) counts ordinary unlabelled tournaments
by vertex count, starting at the empty tournament. It does not count paths.
[A135244](https://oeis.org/A135244) gives, for n>=2, the largest m with
s(m)=n, or zero if none, where s is the proper-divisor sum.
[A002496](https://oeis.org/A002496) lists prime values k^2+1.
The primary b-files were retrieved and retained with hashes.

The prior [aliquot board](../aliquot_power_124_20261026/CONNECTION_BOARD.md)
already records finite inverse fibres except at target one, prime-square
predecessors, and the loss of arithmetic labels in an unweighted shape.
The [prime-dilation response](../prime_dilation_response_20261030/PROOFS.md)
already retains prime exponents and valuations. The
[Fortunate/Landau note](../fortunate_landau_bridge_20261110/README.md)
already distinguishes quadratic values from their prime factors and proves
that a fixed finite congruence wheel cannot force primality.

## 2. A1: a complete three-stratum inverse enumerator

Let Omega(m) count prime factors with multiplicity. For n>=2, a predecessor
m cannot be prime. If Omega(m)>=r>=2, its smallest prime factor p obeys
p<=m^(1/r), and the proper divisor m/p contributes to s(m). Therefore

    s(m)=n and Omega(m)>=r  =>  m<=n^(r/(r-1)).       (A1)

In particular every predecessor is in exactly one of these families:

1. m=p^2, with p prime; then n=p+1.
2. m=pq, with distinct primes p<q; then n=p+q+1.
3. Omega(m)>=3; then m<=floor(sqrt(n^3)).

This gives the exact identity

    A135244(n)=max(square_max(n), distinct_semiprime_max(n),
                  high_Omega_max(n)),                (A2)

where a missing stratum contributes zero. For even n the distinct-semiprime
candidate must be 2(n-3), and exists only if n-3 is an odd prime. For odd n,
choose distinct primes p<q with p+q=n-1; maximizing pq is exactly choosing
the largest eligible p below (n-1)/2. The equal-prime case MUST be excluded:
s(p^2)=p+1, not 2p+1. For example 3^2 is not a predecessor of seven.

For all targets through N we can sieve s(m) and Omega(m) only through
floor(sqrt(N^3)), then add the square and semiprime strata explicitly.
This cuts the divisor-sieve range from quadratic to N^(3/2). Searching
semiprime candidates still has a cost; no near-linear runtime is asserted.
If a found semiprime product exceeds floor(sqrt(n^3)), it certifies that
the third stratum cannot improve it. This certificate is unconditional;
existence of suitable prime pairs for every odd n is not assumed.

**FINITE-EXACT:** at N=10,000 the high-Omega sieve ends at 1,000,000.
The older generic quadratic bound is 99,980,001. All 9,999 primary
A135244 entries agree, and trial-division inversion independently agrees
for every target 2..150 over the full generic source bound.
The winning-family counts are:

| family | targets |
|---|---:|
| no predecessor | 1,212 |
| prime square | 1,229 |
| distinct semiprime | 6,015 |
| at least three prime factors | 1,543 |

These finite frequencies imply no limiting density.

## 3. A2: the exact bridge between A135244 and A002496

For any n>=2, every predecessor satisfies m<=(n-1)^2, and equality
occurs precisely when n-1 is prime and m=(n-1)^2.

Proof. If m=p^2 is a prime square, n=p+1. Otherwise take the smallest
prime divisor p of m. The proper divisors 1,p,m/p are distinct, so
n>=1+p+m/p>=1+2sqrt(m). Thus the non-prime-square predecessors even
satisfy m<=(n-1)^2/4. Equality in the larger bound is impossible there.

Consequently, for every integer k>=1,

    k^2+1 is prime
      iff A135244(k^2+2)=(k^2+1)^2.                    (A3)

This is an actual equivalence, with the map p -> (target p+1, predecessor
p^2), not a numerical analogy. Computing an inverse fibre is generally
more work than testing p for primality; the equivalence is a structural
bridge and not advertised as a faster primality test. The first 99
quadratic inputs whose targets fit within 10,000 were checked directly.

## 4. The specified small numbers have exact inverse-fibre roles

| n | A135244(n) | factorization / role |
|---:|---:|---|
| 2 | 0 | no predecessor |
| 3 | 4 | 2^2 |
| 5 | 0 | no predecessor |
| 7 | 8 | 2^3 |
| 11 | 21 | 3*7 |
| 13 | 35 | 5*7 |
| 19 | 77 | 7*11 |
| 27 | 133 | 7*19 |
| 33 | 247 | 13*19 |
| 1039 | 268577 | 491*547 |

For 33, the balanced distinct-prime pair 13+19=32 gives 247, which
exceeds floor(sqrt(33^3))=189; the other strata cannot improve it.
Thus the appearance of 13 and 19 at 33 is explained by an exact
sum/product optimization. It is separate from 33/3=11 in the previous
tournament orbit calculation. No map equating those two mechanisms is
claimed. Likewise seven being a missing full tournament path count does
not make seven an untouchable aliquot target: eight is its predecessor.

## 5. An exact growing sieve for quadratic prime values

For each odd prime ell, the congruence k^2=-1 mod ell has two roots if
ell=1 mod4 and no roots if ell=3 mod4. The prime two eliminates odd k
apart from k=1, which gives the prime value two. For each ell=1 mod4,
a quadratic nonresidue g gives a root g^((ell-1)/4).

For all k<=K, sieve both root classes for every prime ell<=K, but retain
the special equality k^2+1=ell. This is a complete test: every composite
k^2+1 has a prime divisor <=sqrt(k^2+1), hence <=k. The cutoff grows
with K; a fixed wheel would leave composite survivors, as the inherited
Fortunate/Landau note demonstrates.

**FINITE-EXACT:** for 1<=k<=100,000 this yields 6,656 prime values,
all agreeing with the corresponding primary b-file prefix. The entire
100,000-input decision vector agrees with a separate deterministic
64-bit Miller-Rabin calculation, and inputs through 1,000 also agree
with trial division. Correctness of the full sieve is proved above;
the second algorithm is an additional error detector.

## 6. The genuine, limited tournament transfer

Write an even k uniquely as k=d*2^a with odd d and a>=1. Then

    k^2+1=1+d^2*2^(2a).                              (A4)

Thus quadratic values occupy the odd-square coefficient rays, at even
exponents, inside the unit-offset seed construction. The map preserves the
actual integer Hamiltonian-path target, but has no primality predicate.
The exceptional input k=1 gives the prime two, which cannot be a tournament's
Hamiltonian-path count because those counts are odd.
Its sidecars are d, the exponent parity, the seed base, and the selected
attachment. A composite quadratic value is still a valid tournament target.

We reuse the [complete nine-base seed certificates](../metagraph_aspect_sequences_20261111/README.md),
not discover new seeds. For each of the 32 odd d<=63, the coefficient d^2
has such a certificate. The inherited tail formula is

    H(closure(B,O,t))=H(B)-W(O)+2^t U(O), t>=1.

Here W=H(B)-1 and U=d^2*2^j, so the result is 1+d^2*2^(j+t).
The saved quadratic_tournament_rays.json retains all bases and masks,
rounds the start exponent up to an even number, and checks two consecutive
literal tails for each base. Every rounded start is at most eight.
Consequently every value k^2+1 with k=d*2^a, odd d<=63 and a>=4,
has a certified tournament construction. This covers infinite candidate
rays; it neither proves infinitely many prime values on a ray nor supplies
seeds for unbounded d.

There is also a useful symmetry boundary. For a prime p=k^2+1>2,
k is even and p=1 mod4. The quadratic residues contain -1, so orienting
x->y by whether y-x is a quadratic residue gives both directions or
neither, not a tournament. These primes support the undirected Paley
construction, while the earlier p=7,1039 quadratic-residue tournaments
use p=3 mod4. Other tournament constructions on p vertices are possible;
the quadratic-residue orientation specifically fails. Never transfer
Paley tournament claims across this modulus distinction.

The arithmetic sieve, the inverse-fibre equivalence, the seed ray map,
and the tournament class census therefore have explicit but different
interfaces. None proves Landau infinitude, Chan-Pak, Collatz convergence,
or a relation between those open problems.
