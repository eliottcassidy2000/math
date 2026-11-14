# The exact 30,017 source under repeated primorial translation

This is a different mechanism from sieving the primality of exponential rails. Here the question is whether the *same prescribed H* can provide a physical clock Q after the target is translated by a primorial. At H=23 the answer can be calculated for every translate, with no size cutoff.

## 1. General constant-H return law

In a full digit bank of slope a, choose a fixed odd h prime to a. For the progression N=N₀+Pt, with P>0 even, the unique digit c is periodic in t with period d=2a/gcd(P,2a). Split into t=j+dℓ. On each lane the effective clock forced by the target is

    Q_eff = B_j + Δℓ,    B_j=(N₀+Pj−c_j h)/a,    Δ=Pd/a.

For fixed h, list the finitely many ordered strings of nonsingleton strong SCCs whose path counts multiply to h. Finiteness follows from Busch's strong bound and the fact that each such count is at least three. Their Q products are fixed integers q_i. All physical Q values are exactly q_i·2^s, s≥0, because the leftover SCCs are singletons. The selected digit phase r means Q_eff=q_i·2^u with u≥r. Thus each return is determined exactly by

    q_i·2^u ≡ B_j (mod Δ),   q_i·2^u ≥ B_j,   u≥r,

and ℓ=(q_i·2^u−B_j)/Δ. This is a small transient test at the two-adic part, followed by a finite subgroup/period test for powers of two at the odd part. It is not enough that the digit itself returns. A bounded h bank therefore produces at most logarithmically many exact returns below large t. More importantly, this calculation predicts *where* its candidate sources fail without pretending that all h can be held fixed.

## 2. The 30,017 return

Take a=99, h=23, N₀=30017, and P=30030. The [earlier exact witness](../session326_exact_source_profiles_20261114/PRIMORIAL_SWITCH.md) uses Q=296 and digit31. Since 23 is prime, its source has exactly one nonsingleton strong SCC. Busch excludes order seven because the order-seven minimum H is 25. Literal enumeration through six gives exactly the two possible strong clocks Q=272,296 at H=23. With singleton padding these give its entire clock spectrum, not merely a sample.

The translation step has period three in the digit ring. The exact phase data are:

| t mod 3 | digit c | phase r | effective clock on t=j+3ℓ |
|---:|---:|---:|---|
| 0 | 31 | 0 | 296+910ℓ |
| 1 | 97 | 0 | 584+910ℓ |
| 2 | 163 | 3 | 872+910ℓ |

For example 30030≡132 mod198 and 23⁻¹≡155, so c advances by 66 each time. The resulting clock increment within the first three lanes is (30030−66·23)/99=288; after a full three-translation cycle it is 910. In particular t=3 returns the digit but forces Q=1206, which is 2 mod4 and impossible. This is an exact source-phase switch rather than a prime sieve.

Divide the clock congruence by 2 and use the odd modulus 455=5·7·13. The residues required of 2^u in each lane are

| j | clock core 272 | clock core 296 |
|---:|---:|---:|
| 0 | 68 | 1 |
| 1 | 417 | 334 |
| 2 | 311 | 212 |

The powers of two have order 12 modulo455, by the orders 4,3,12 modulo 5,7,13. Of these six entries, only 1 occurs in that subgroup. Thus the whole answer for t≥0 is

    N is realizable in the 99 digit bank at fixed H=23
    iff t = 3·296(2^(12k)−1)/910 for an integer k≥0.

At k=0 this returns N₀. The first positive return is t=3996, N=120,029,897, at the actual effective clock 296·2^12. The integer 3996 is not fitted: 2^12−1=4095=9·455, so ℓ=(296/2)·9=1332 and t=3ℓ. Conversely each such translate has the literal 296 clock core with 12k singleton pads and the original digit31. All other translates fail at this fixed H, not necessarily with other sources or cuts.

This is one rigorous reason the numbers around 30,000 keep suggesting a recurrence: the primorial sets the modular phase of the digit, while the physical clock sets which turns of that phase are affordable. When the two constraints are kept separate, the same local source makes a precise prediction far beyond the original neighborhood.
