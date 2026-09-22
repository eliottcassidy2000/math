# Repeated digits, prime novelty, and the arithmetic of return periods

**Date:** 2026-09-21. **Status:** PROVED elementary statements; FINITE-EXACT
primality/factorization and primitive-part controls; CITED classical input
only where marked. No claim of a new proof of Zsigmondy's theorem, no
universal primitive-divisor theorem for the decimal family, and no Collatz
convergence consequence.

## Inheritance and live concepts

Read first: [the Zsigmondy triad note](collatz_mod6_20260917_zsigmondy_triad.md),
especially §§1 and 7; [THM-4139, rational three-cycle/order-six lift](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md),
as routed there; and the maintained guardrails/protocol. The inherited
mechanisms are the Mersenne exception 63, multiplicative order, and the
pairwise coprimality of Fermat numbers. Their proofs are not presented
as discoveries of this session. The hostile is **composite but all prime
factors new**, which distinguishes two different proposed failure patterns.
The least-used coordinate is the shifted first-hit index of a decimal
affine orbit, as opposed to an order alone.

| Lane | Object and operation | Retained information / decisive test |
|---|---|---|
| Anchor | Decimal threes ending in one; append a digit | Exact recurrence, first composite, first-hit offset |
| Niche | Shared factors of nearby terms | Index difference, multiplicative order, prime powers |
| Wildcard | Mersenne/Fermat comparison | Primitive prime versus primality; target 1, -1 or 7 |
| Boundary | Finite congruence sieves | Avoid a prescribed finite prime pool; no prime-infinitude conclusion |

## 1. The first composite is the ninth-digit number

Let R_k have k decimal threes followed by one, and include `R_0=1`.
Then

```text
R_k=(10^(k+1)-7)/3,           R_(k+1)=10R_k+21.              (P1)
```

The following is an exact bounded primality/factorization statement,
verified by complete trial division:

| k, the number of threes | R_k | Status |
|---:|---:|---|
|1|31|prime|
|2|331|prime|
|3|3331|prime|
|4|33331|prime|
|5|333331|prime|
|6|3333331|prime|
|7|33333331|prime|
|8|333333331|`17 * 19607843`, both factors prime|

Thus the first composite has **eight threes and nine digits**. In
particular, assigning its factorization to the preceding eight-digit
number would be an indexing error.

Prime 17 occurs in exactly the index class

```text
17 | R_k   iff   k=8 mod16,                                (P2)
```

because `ord_17(10)=16` and `10^9=7 mod17`. Later members of this
class are composite because they exceed 17. This supplies an infinite
composite progression, not a classification of all composite terms.

## 2. General append-digit families and an exact gcd law

Fix a base `b>=2` and a nonzero digit `1<=d<b`. Let A_k consist of k
copies of d followed by the digit one in base b. Set

```text
U_t(b)=(b^t-1)/(b-1),
A_k=1+d*b*U_k(b),
c=(d-1)*b+1.
```

Then `A_0=1`, `A_(k+1)=b*A_k+c`, and `gcd(A_k,bc)=1` for every k.
The latter follows from the last digit and induction using
`gcd(b,c)=1`. Therefore, for `m>=0,t>=1`,

```text
A_(m+t)=b^t*A_m+c*U_t(b),
gcd(A_m,A_(m+t))=gcd(A_m,U_t(b)).                            (P3)
```

This is the useful general theorem: factor reuse is controlled by an
index **difference** and a repunit return time. It does not say that
`A_m` divides `A_n` whenever m divides n. The decimal hostile is already
`31` not dividing `331`.

For a prime power q that divides some A_k, define its repunit return
rank `rho_q(b)` as the least positive t with `q|U_t(b)`. This rank exists:
`x -> b*x+1` is a permutation modulo q, so its orbit of zero returns.
Equation (P3) implies

```text
q | A_(k+t)   iff   rho_q(b) divides t.                     (P4)
```

All occurrences are therefore one residue class modulo this rank.
When `gcd(q,b-1)=1`, the rank is simply `ord_q(b)`. For a prime
`p|b-1`, its rank is p, since `U_t(b)=t modp`; such a prime can occur
in the digit family only when `p` does not divide d. This exception
must not be hidden by dividing by b-1 modulo p.

If p occurs and a_p is its first index, then `1<=a_p<rho_p(b)`.
Consequently a prime factor of A_k is primitive relative to
`A_0,...,A_(k-1)` exactly when `k<rho_p(b)`. The rank controls the
spacing, while the first-hit offset controls when the prime is new.

If d divides b-1, an especially simple affine coordinate is available.
Put `h=(b-1)/d`. Then

```text
h*A_k = b^(k+1)-(b-h).                                     (P5)
```

Digit d=1 gives repunits, with target b-h=1. Decimal d=3 gives h=3
and target seven. Changing that target changes the first-hit arithmetic;
it cannot be discarded while importing a primitive-divisor theorem.

## 3. Every fifteen consecutive decimal terms are pairwise coprime

For the decimal family, every R_k is coprime to 210. In particular
(P3) simplifies to

```text
gcd(R_m,R_(m+t))=gcd(R_m,10^t-1).                           (P6)
```

**PROVED, sharp uniform window.** For every `m>=0`, the fifteen terms
`R_m,...,R_(m+14)` are pairwise coprime. The window length cannot be
increased to sixteen uniformly: `gcd(R_1,R_16)=31`.

Suppose a prime p divides two terms separated by t. It is not among
2,3,5,7, and (P6) gives `10^t=1 modp`. Also `10^(m+1)=7 modp`, so
raising this to the tth power gives `7^t=1 modp`. Hence p divides

```text
G_t=gcd(10^t-1,7^t-1).
```

The complete exact table for the required differences is

```text
t:    1  2  3  4  5  6  7  8   9  10 11  12 13 14
G_t:  3  3  9  3  3  9  3  3 999  33  3 117  3  3.
```

The only candidate primes are 3,11,13,37. Prime three never divides
R_k. For the remaining candidates the powers of ten form the subgroups

```text
mod11: {1,10};
mod13: {1,10,9,12,3,4};
mod37: {1,10,26}.
```

None contains seven, so none ever divides R_k. This excludes every
gap t from one through fourteen, at **all** source indices. At gap
fifteen, `ord_31(10)=15` and `10^2=7 mod31`, proving sharpness.
The finite table is an exact certificate for a universal window theorem,
not a scan of a finite initial segment pretending to be that theorem.

**Consequence.** Every prime factor of each of `R_1,...,R_15` is
primitive. Thus the first composite R_8 introduces two new primes;
its loss of primality is not a primitive-prime failure.

## 4. Prime-power clocks retain a coordinate that prime support loses

For decimal terms, every occurring prime p is coprime to 30, so its
prime-power occurrence spacing is `ord_(p^a)(10)`. If

```text
h=ord_p(10),             v_p(10^h-1)=1,
```

the odd-prime lifting identity gives
`ord_(p^a)(10)=h*p^(a-1)`. Every first-hit class modulo h lifts to a
unique class modulo `h*p^(a-1)`: the powers of `10^h` generate the
principal units congruent to one modulo p. The condition on valuation
is required; no unconditional ordinary-lift assertion is made here.

Exact examples are:

| Prime power | First index k with `p^a|R_k` | Period |
|---|---:|---:|
|17|8|16|
|17²|248|272|
|17³|2152|4624|
|31|1|15|
|31²|226|465|
|31³|3481|14415|

The inherited 63 mechanism has this same **order-lifting operation**:
`ord_3(2)=2`, while `ord_9(2)=6`. A higher power of an existing prime
supplies the larger period. This is the relevant shared structure;
the first composite in a digit family and a missing primitive prime
in a Mersenne sequence are different predicates.

## 5. Mersenne, Fermat and repeated threes: the precise comparison

Here a primitive prime means a prime not dividing any earlier term of
the specified sequence. The three modular tests retain different targets:

| Sequence | Prime-divisor condition | Novelty mechanism |
|---|---|---|
|`M_n=2^n-1`|`2^n=1 modp`|primitive iff `ord_p(2)=n`|
|`F_n=2^(2^n)+1`|`2^(2^n)=-1 modp`|every prime factor has order `2^(n+1)`|
|`R_k=(10^(k+1)-7)/3`|`10^(k+1)=7 modp`, `p` coprime to 210|primitive iff k is the first-hit offset, equivalently `k<ord_p(10)`|

The Fermat condition forces that exact order because the order divides
`2^(n+1)` but does not divide `2^n`. It proves that distinct Fermat
numbers have disjoint prime supports. In particular,

```text
F_0,...,F_4 = 3,5,17,257,65537             are prime,
F_5=4294967297=641*6700417                 is composite,
```

and **both** factors of F_5 are new. The earlier Fermat identity
`F_n-2=product_(j<n)F_j` gives the same coprimality; it is inherited,
and the root lane separately connects it to its translated operations.

By contrast `M_6=63=3^2*7` has no new prime, because three already
divides M_2 and seven divides M_3. Equivalently no prime has order
six for base two. The statements about 63 and the quadratic critical
orbit at `-7/4` are already audited in the inherited triad note; nothing
here identifies those dynamical systems with the decimal recurrence.

**CITED classical scope.** The Bang–Zsigmondy theorem supplies a primitive
prime of `a^n-b^n` for every `n>6` when `a>b>0` are coprime. Combining
that with direct n=1,...,6 calculations gives exactly n=1 and n=6 as
the primitive-prime exceptions for `2^n-1`. This standard input is
recalled in the introduction of the research paper
[Voutier–Yabuta, *Primitive divisors of certain elliptic divisibility sequences*](https://www.impan.pl/shop/en/publication/transaction/download/product/82701),
Acta Arithmetica 151.2 (2012), pp.165–190. The original archival record is
[Zsigmondy, *Zur Theorie der Potenzreste* (1892)](https://zenodo.org/records/2131326).
The archive metadata was retrieved; its scanned PDF did not load in
this session, so the modern paper is the retrieved theorem statement.

The classical hypothesis is not met by replacing b^n with the constant
seven. Nor does merely satisfying the second-order recurrence
`R_(k+2)=11R_(k+1)-10R_k` establish the Lucas-pair hypotheses of
[Bilu–Hanrot–Voutier's primitive-divisor theorem](https://oskar-bordeaux.fr/handle/20.500.12278/166397).
That paper concerns its defined Lucas and Lehmer sequences. No claim
that it applies to R_k is made.

### 5a. An exact Fermat-prime bridge: ten is a generator

**PROVED.** If `p=2^(2^n)+1` is prime with `n>=2`, then

```text
ord_p(10)=p-1.                                             (P7)
```

Consequently every such Fermat prime divides the decimal family R_k,
in exactly one index class modulo p-1. This is a conditional statement
for every Fermat prime, not an assertion that further Fermat primes
exist.

Here is an elementary proof retaining the modular mechanism. Such a
p is `17 mod40`; write `p=40t+17`. In the list `10j modp`,
`1<=j<=(p-1)/2`, count the residues greater than p/2. They correspond
to the five intervals

```text
(2i+1)*p/20 < j < (2i+2)*p/20,       i=0,...,4.
```

Each contains exactly `2t+1` integers, so the count is `10t+5`, odd.
Replace every such residue by its negative. The resulting absolute
residues permute `1,...,(p-1)/2`; multiplying and cancelling their
nonzero factorial gives

```text
10^((p-1)/2)=-1 modp.
```

Since p-1 is a power of two, the order cannot be a proper divisor of
p-1. This proves (P7). As ten generates all nonzero residues, its powers
hit seven exactly once per period; division by three is legal here.

| Fermat prime | `ord_p(10)` | First k with `p|R_k` |
|---:|---:|---:|
|3|1|never: `R_k=1 mod3`; division by three cannot be discarded|
|5|undefined, ten is not a unit|never: `R_k=1 mod5`|
|17|16|8|
|257|256|226|
|65537|65536|29252|

This supplies a genuine link between the two families: the Fermat
prime's power-of-two unit group and the explicit target seven determine
an occurrence clock in the decimal sequence. It does not prove a
shared cause for their initial runs of prime terms or identify their
primitive-divisor behavior.

## 6. What the finite calculation does and does not establish

The exact experiment strips from R_k every prime already present in
`R_1,...,R_(k-1)` by repeated gcd division. The residual exceeds one
for **every `1<=k<=200`**. This proves that each of these 200 terms has
a primitive prime without requiring their complete factorizations.
An independent implementation strips against individual earlier terms,
rather than their product, and gives the same residuals.

This is FINITE-EXACT. A primitive prime at every subsequent index is
not proved here, and no classification of the prime R_k is claimed.

There is an elementary all-index survivor worth retaining. For **any**
base/digit family of §2, infinitely many distinct primes divide its
terms. Given a finite pool P of primes that have occurred, the map
`x -> bx+c` is a permutation modulo their product because none divides
b. Its orbit from `A_0=1` returns after a finite period H. Consequently
every `A_(tH)=1` modulo each prime in P. These positive terms exceed
one for t>0, so their prime factors lie outside P.

Thus no fixed finite collection of prime-divisibility progressions
covers all sufficiently large indices. This does **not** prove infinitely
many prime terms: the uncovered terms may be composite with other
factors. It also does not say that A_(tH) is primitive relative to every
intervening index, only that it avoids the chosen finite pool.

## Reproduction and remaining boundary

```text
python 04-computation/experiments/arithmetic_seams_20260921_primes.py
python -O 04-computation/experiments/arithmetic_seams_20260921_primes.py --output primes-optimized.json
```

The [source](../../04-computation/experiments/arithmetic_seams_20260921_primes.py)
and [JSON](../../04-computation/experiments/arithmetic_seams_20260921_primes.json)
record complete trial divisions for the first-composite claims, the
fourteen-entry window certificate, prime-power periods/offsets, a bounded
general-base gcd audit, and both primitive-part implementations. Checks
remain active under `python -O`. None of the primality verdicts relies
on probable-prime testing.

The source/target correspondence is multiplicative motion modulo a
prime; it preserves the hitting congruence and return period. It loses
the target and first-hit offset if reduced to an order alone, and loses
prime-power multiplicities if reduced to prime support. Retaining
those coordinates explains both the sharp window and the different
failure patterns. It does not force global prime production at every
index or a Collatz descent law.

**Independent audit:** the summand lane read the proof and code for the
general gcd/rank law, denominator exception, sharp fifteen-term window,
prime-power lift hypothesis, finite primitive-part census, and finite-pool
avoidance. A separate focused review checked §5a's Gauss intervals,
factorial cancellation, power-of-two order conclusion, and indexed
Fermat-prime controls. Both reviews passed without mathematical changes.
Ordinary and optimized runs produced identical JSON. These are elementary
proof reviews and exact controls, not a formal theorem-prover certificate.
