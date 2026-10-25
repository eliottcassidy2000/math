# Fixed binary precision requires only a bounded number of permutation cycles

2026-10-25. **PROVED:** a termwise 2-adic cycle bound, the even-order leading
bit and its equality classification, two normalized low-bit congruences,
and an exact bounded-cycle precision algorithm for A000568. **FINITE-EXACT:**
independent full rational sums and modular evaluation agree through n=50.
This is a low-precision filter of the inherited Burnside formula, not a
Burnside-free formula or a polynomial-time algorithm for the full count.

## 1. Inheritance and exact observable

Write a(n)=A000568(n). The April result
[THM-305, two-adic valuation of tournament counts](../../../../01-canon/theorems/THM-305-two-adic-valuation-tournament-count.md)
already proves v_2(a(n))=(n-1)/2 for odd n>=3. Its associated
`04-computation/qr_tournament_enumeration_s25.py` contains an even-order
exploratory table, without the even-order formula below. A bounded search
of the current theorem/result surfaces for that theorem, tournament
valuations and totient leading bits did not recover an existing statement
of the latter formula. No historical-priority claim is made.

The June [THM-514, unit-tail and prime-tail kernel](../../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md)
retains the full count while peeling chosen cycle lengths. Here the
observable is smaller: K binary digits after a proved common power of two
has been removed. The canonical hostile is reducing a numerator before
cancelling the power of two in its centralizer denominator. The less-used
sidecar is that denominator's *odd part*, which is invertible at the
requested precision. The concept board is cycle count / pairwise gcd /
multiplicity factorial / coefficient valuation / exact modular units.

For an odd partition lambda=(c_1,...,c_m) of n, let

    G(lambda)=sum_(i<j) gcd(c_i,c_j),
    e(lambda)=(n-m)/2+G(lambda),
    z(lambda)=product_d d^(m_d)*m_d!.

The inherited fixed-permutation formula is

    a(n)=sum_(lambda odd partition of n) 2^e(lambda)/z(lambda). (P1)

All valuations below are valuations of rational summands. They need not
be integers individually in the ordinary sense; their odd denominators
are legitimate units in the 2-adic integers.

## 2. A uniform quadratic bound in the number of cycles

Because all c_i are odd,

    v_2(2^e/z)=(n-m)/2+G(lambda)-sum_d v_2(m_d!).     (P2)

Every gcd in G is at least one, while the multinomial coefficient
m!/product_d m_d! is integral. Consequently

    v_2(2^e/z) >= (n-m)/2+binom(m,2)-v_2(m!).        (P3)

For odd n put B=(n-1)/2. Then m is odd and the relative lower bound is

    beta_odd(m)=(m-1)^2/2-v_2(m!).                   (P4)

It is zero at m=1, one at m=3, and positive thereafter. For example,
v_2(m!)<=m-1 gives beta_odd(m)>0 for m>=5. Thus the single n-cycle is
the unique term of valuation B; its normalized coefficient is 1/n, odd
2-adically. This reproves the inherited odd theorem by one uniform bound.
The bound is sharp as a termwise bound: the partition
(n-m+1,1,...,1) attains it for every admissible odd n>=m, with the all-ones
case interpreted using its actual multiplicity factorial.

For even n>=4 put B=n/2. For even m>=4 the relative bound is

    beta_even(m)=m(m-2)/2-v_2(m!).                   (P5)

It is at least one. The case m=2 needs its own argument; inserting m=2
into (P5) would give -1 and would not prove the claimed baseline.

If lambda=(a,b), a<b odd, its exact relative valuation is gcd(a,b)-1.
If lambda=(d,d), d odd, then n=2d and its exact valuation is n-2.
Here d>=3, since n>=4, so this is at least B+1. Hence every summand at
even n>=4 has valuation at least B, and the terms attaining B are exactly
the distinct coprime odd pairs a+b=n.

## 3. Even leading bit and exact equality classification

Each distinct coprime pair contributes the odd normalized unit 1/(ab).
All other terms are even after division by 2^B. Such unordered pairs
are counted by phi(n)/2: the map a->n-a pairs the units modulo the even
integer n, and has no fixed unit for n>=4. Therefore

    a(n)/2^(n/2) = phi(n)/2 mod 2,       even n>=4.    (P6)

In particular,

    v_2(a(n))=n/2
      iff n=4, or n=2*p^a with p=3 mod 4 prime and a>=1. (P7)

For completeness, write n=2^s product_i p_i^(a_i), s>=1. Then
v_2(phi(n))=max(s-1,0)+sum_i v_2(p_i-1). This is one precisely when
n=4, or s=1 with one odd prime factor p having v_2(p-1)=1. That proves
(P7), rather than inferring it from a finite valuation table.

Orders 0,1,2 are handled separately: their counts are all one. In
particular n=2 does not satisfy the even baseline n/2; its partition
(1,1) has valuation zero.

## 4. Two more low-bit consequences

For odd n>=3,

    a(n)/2^((n-1)/2) = n+2 mod 4.                   (P8)

Only m=1 and m=3 can matter. The m=1 term is 1/n=n modulo 4. Among
three-cycle partitions, relative valuation one requires a repeated part
equal to one: (1,1,n-2), including the special all-ones partition at n=3.
Its normalized term is twice an odd unit and therefore contributes two
modulo four. Distinct triples have relative valuation at least two;
all other repeated triples have larger valuation.

For even n>=4,

    a(n)/2^(n/2)
      = (-1)^(n/2+1)*phi(n)/2 + 2*[n in {4,6}] mod 4. (P9)

For n>=8 only the distinct coprime two-cycle terms matter modulo four.
Their products ab are -1 modulo four if n=0 modulo four, and +1 if
n=2 modulo four. This gives the first term of (P9). Noncoprime pairs
have relative valuation at least two; equal pairs contribute at relative
valuation one only for (3,3), n=6.

Here is the needed four-cycle check, beyond the coarse bound (P5).
Multiplicity patterns 3+1, 2+1+1, and 1+1+1+1 have relative valuations
at least 3,3,4 respectively. Pattern 2+2 has G=a+b+4 gcd(a,b)>=8 and
relative valuation at least four. Four equal parts d have relative
valuation 6d-5, so relative valuation one occurs only for four ones,
n=4. For m>=6, (P5) is at least eight. This proves precisely the two
corrections in (P9).

For example the normalized even leading bit vanishes at n=50, but (P9)
gives the normalized residue two modulo four, determining v_2(a(50))=26.

## 5. The fixed-precision cycle budget

Fix K>=1. To compute a(n) modulo 2^(B+K), equivalently its normalized
quotient modulo 2^K, retain only

    odd n:  odd m with beta_odd(m)<K;
    even n: m=2, and even m>=4 with beta_even(m)<K.    (P10)

Every discarded partition has valuation at least B+K by (P3)--(P5).
There is no cancellation assumption: those terms individually vanish
at the requested precision.

The number of cycles is uniformly bounded independently of n. Using
v_2(m!)<=m-1 gives

    beta_odd(m)>=(m-2)^2/2-1/2,
    beta_even(m)>=(m-2)^2/2-1.

Therefore every m needed by (P10) lies below

    2+sqrt(2K+2),

apart from harmless endpoint overestimates when rounding. The code uses
the safe cap 2+floor(sqrt(2K+2)) and then applies the exact beta test.
The cycle budget is O(sqrt(K)), independent of the target order.

For K=1,...,6 the budgets are particularly small:

| K | Odd-order cycle counts | Even-order cycle counts |
|---:|---|---|
| 1 | 1 | 2 |
| 2,3,4,5 | 1,3 | 2,4 |
| 6 | 1,3,5 | 2,4 |

Generate partitions only at these permitted lengths. A simple exact
recursion for nondecreasing odd parts chooses the first part d from
the odd integers between the current lower bound and floor(n/m), then
recurses on (n-d,m-1,d). It never generates the excluded high-cycle
partitions in the production path.

For a retained partition write z=2^t*z_odd. Compute

    alpha=e-t-B.

If alpha>=K, discard the term too. Otherwise its normalized contribution
is exactly

    2^alpha * inverse(z_odd modulo 2^K) modulo 2^K.  (P11)

Every inverse exists. Summing (P11) gives the required normalized residue;
multiplication by 2^B restores the absolute residue modulo 2^(B+K).
The production computation never forms a rational partition summand and
never inverts an even denominator.

## 6. Meaning, limits and reproduction

The map from an odd partition to (alpha,z_odd modulo 2^K) preserves its
contribution at the chosen precision. Deleting alpha>=K loses higher
binary digits deliberately. The missing coordinate when precision grows
is the discarded cycle/valuation tail; a K-bit result cannot be reused
at higher precision without adding that tail. This is the same exact
observable discipline as the recent finite-horizon moment work, with a
different observation map.

This method still ranges over arithmetic cycle partitions of the target
order. Their lengths are small at fixed K, but the part values are not
bounded independently of n. No polynomial-time full-count claim follows:
the precision needed to recover the entire integer grows with n, and so
does its permitted cycle budget. No literal tournament core or Hamiltonian
path is enumerated, and no port-transport conclusion is inferred.

Run `python check.py > normal.json` and `python -O check.py > optimized.json`.
The independent control sums all odd partitions using exact Fractions for
n=0,...,50 and compares all K=1,...,6 to the bounded-cycle modular path.
It checks each omitted term's valuation, the odd/even leading laws, both
second-bit congruences, the explicit equality classification, denominator
and small-order hostiles, and invalid inputs. There are 31,583 full
partition records and 306 modular comparisons. All checks pass normally
and under optimized Python.

At n=50 there are 3,658 full odd partitions. K=1 visits 13 two-cycle
partitions, of which ten affect the residue. K=6 visits 163 partitions
with two or four cycles, of which 61 affect the residue, yielding
a(50)/2^25=46 modulo 64. These are operation counts, not timing claims.
