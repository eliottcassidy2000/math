# Decimal digits, ternary clocks, and the exact carry behind the 333...31 family

Status: **PROVED** for the elementary statements and proofs below; **FINITE-EXACT** for the explicitly bounded reproduction; **CITED** for external searches and large primality certificates; **OPEN** for circular-prime finiteness, repunit-prime infinitude, and a universal Collatz descent argument. This is a synthesis of inherited mechanisms, with no literature novelty claim. No theorem ID is reserved.

The useful positive result is concrete: for

\[
N_k=(10^k-7)/3,\qquad U(n)=(3n+1)/2^{v_2(3n+1)},
\]

every `k >= 4` satisfies `U^2(N_k) < N_k`. The reason is an exact binary carry, independent of whether `N_k` is prime. Separately, the digit-append family and the Collatz sibling ladder have compatible, isometric ternary residue clocks. These two statements preserve different information.

## Inheritance and concept board

Closest proved mechanism: **THM-4473**, [collatz-digit-chains-rotation-repunits](../../01-canon/theorems/THM-4473-collatz-digit-chains-rotation-repunits.md), especially digit rotation as multiplication modulo `10^k-1` and the exact rising run of `2^k-1`. Its long [digits/rotation/repunit note](procgen_repunit_20260924_digits_rotation_repunits.md) supplies the prior circular-prime computation. The sibling-ladder valuations and the absence of an automatic rank are also inherited from [creative_sibling_20260925](creative_sibling_20260925.md).

Canonical hostile: composite binary repunits `2^(2m)-1` have arbitrarily many consecutive growing odd steps. Corrected near misses: a prime run is not a circular-prime run (`331 -> 133=7*19`); an isomorphism of ternary clocks does not preserve the ordinary numerical unit guard without its coordinate dictionary. Least-used relevant sidecar: the leading decimal digit in the rotation carry, together with the binary valuation used by the actual Collatz step.

The live board is: (1) decimal append, (2) circular rotation, (3) ternary residue clocks, (4) binary Collatz carry, (5) multiplicative order of prime factors. Appending and rotating are distinct operations. A factor's multiplicative order gives a repeating congruence class of indices; it does not give a well-founded orbit rank.

## 1. Audit of the supplied arithmetic

For `k >= 1`, `N_k` has `k-1` copies of the digit 3 followed by 1. Thus `N_1=1`, and the displayed prime run begins at `k=2`.

| k | N_k | Exact status |
|---|---:|---|
| 2 | 31 | prime |
| 3 | 331 | prime |
| 4 | 3331 | prime |
| 5 | 33331 | prime |
| 6 | 333331 | prime |
| 7 | 3333331 | prime |
| 8 | 33333331 | prime |
| 9 | 333333331 | `17 * 19607843` |

The seven prime entries and both factors on the last line are verified by exhaustive trial division. Therefore “the last member of the initial prime run has eight digits and the next has nine” is correct. “Eight prime steps beginning at 31” is an off-by-one error: there are seven. No conclusion about all later indices follows from the initial run.

The Fermat arithmetic is correct with `F_j=2^(2^j)+1`:

\[
F_2=17,\qquad F_5=4294967297=641\cdot6700417.
\]

The factors of `F_5` are independently trial-division verified as prime. A short certificate for the divisor 641 uses

\[
641=5\cdot2^7+1=5^4+2^4.
\]

Modulo 641, the first equation to the fourth power gives `5^4 2^28 = 1`; the second then gives `-2^32 = 1`. This proves the divisibility without invoking a prime-generation principle.

For the cubic,

\[
a_n=2n^3+4n^2+n=n\bigl(2(n+1)^2-1\bigr).
\]

Its first five positive values are `7,34,93,196,355`; `a_2=2*17` and `a_4=14^2` are correct. For every `n>=2`, both factors in the displayed factorization exceed 1, so all those values are composite. The first twelve primes are

`2,3,5,7,11,13,17,19,23,29,31,37`,

whose sum is **197**, not 196.

## 2. A factor is a periodic index guard, not an eventual barrier

**PROVED.** For every integer `k>=1`,

\[
17\mid N_k\quad\Longleftrightarrow\quad k\equiv9\pmod {16}.
\]

Proof: 3 is invertible modulo 17, so divisibility is equivalent to `10^k=7 (mod17)`. Direct powers give `ord_17(10)=16` and `10^9=7 (mod17)`. Exactly one of the sixteen exponent classes works. The smallest positive such index is 9.

In comparison `ord_641(2)=64`, since `2^32=-1 (mod641)` and 64 is a power of two. The two factorizations thus have the same *type* of certificate, a multiplicative-order congruence, but different bases and different target residues. Identifying their displayed integers does not identify their dynamics.

The analogous repunit rule, for prime `p` not dividing 90, is

\[
p\mid R_k^{(10)}\quad\Longleftrightarrow\quad\operatorname{ord}_p(10)\mid k.
\]

For `p=3`, division by 9 is not invertible and the correct rule is `v_3(R_k^(10))=v_3(k)`. These exact periodic filters explain how a finite-looking prime pattern can fail. They do not say that the whole family is permanently composite after its first failure.

## 3. Exact ternary clock shared with a sibling ladder

Put

\[
S(n)=4n+1,\qquad A(x)=10x+10.
\]

The decimal append recurrence is `N_(k+1)=10N_k+21`; normalizing `x=(N-1)/3` turns it into `A`. In fact

\[
A^j(0)=10\frac{10^j-1}{9},\qquad 3A^j(0)+1=N_{j+1}.
\]

For all integers `n,x` and positive `t`,

\[
v_3(S^t(n)-n)=v_3(t),\qquad
v_3(A^t(x)-x)=v_3(t).
\]

For the decimal family itself,

\[
v_3(N_{k+t}-N_k)=1+v_3(t).
\]

Proof: the elementary lifting identity gives

`v_3(4^t-1)=1+v_3(t)` and `v_3(10^t-1)=2+v_3(t)`.

Use

\[
S^t(n)-n=(4^t-1)(3n+1)/3,
\]

\[
A^t(x)-x=(10^t-1)(9x+10)/9,
\]

and note that `3n+1` and `9x+10` are 3-adic units. The last formula follows directly by subtracting `N` values. The same proof applies to 3-adic inputs.

Consequently `S` and `A` each form one cycle of length exactly `3^r` modulo `3^r`, for every `r>=1`. Define the finite decoder by its orbit labels:

\[
H_r(S^j(1)\bmod3^r)=A^j(0)\bmod3^r.
\]

It is a well-defined bijection and satisfies `H_r S=A H_r`. Its reductions agree across `r`, because both orbit periods are exactly `3^r`. The compatible maps give an isometric homeomorphism of `Z_3` conjugating these two affine maps. Isometry follows because the valuation of the difference of any two orbit points is exactly the valuation of the difference of their indices, on both sides.

This is a precise ternary recursion isomorphism. Its scope is **affine sibling/append clocks**, not the accelerated Collatz map `U`. The finite clock is also not a finite bound on the integer magnitude: both positive integer append orbits grow indefinitely.

The guarded map is explicit but needs its dictionary. At the first level,

\[
H_1(n)=n-1\pmod3.
\]

Thus the original Collatz inverse-target guard `n not congruent 0 mod3` becomes `H_1(n) not congruent 2 mod3`; it is not preserved as the same numerical predicate. Under the further map `N=3H+1`, all values become `1 mod3`, so discarding `H` loses the original three classes. At depth `r`, keeping `H_r` transports any ternary congruence guard exactly. It does not determine `v_2(3n+1)`.

## 4. Rotation has a finite carry sidecar

For a fixed width `k`, let `d` be the leading digit. Left rotation obeys

\[
\rho_k(n)=10n-d(10^k-1).
\]

Leading zeros must retain that fixed width when applying the operation again. Width loss would incorrectly classify examples such as `101 -> 011 -> 110`.

Since

\[
v_3(10^k-1)=2+v_3(k),
\]

the following statement is exact: **for every width-k string**, rotation is multiplication by 10 modulo `3^r` if and only if `r<=2+v_3(k)`. Sufficiency follows from the displayed carry. For necessity, choose a string with leading digit 1: the carry has exactly that valuation.

Every decimal rotation preserves the residue modulo 9, since `10=1 (mod9)`. Above that precision there are two separate obstructions:

* Multiplication by 10 need not be the identity: `337 -> 373` changes the residue modulo 27, although at width 3 it equals multiplication by 10 modulo 27.
* The leading digit carry need not vanish: at width 2, `13 -> 31` is not multiplication by 10 modulo 27.

For an inverse-Collatz integrality guard `2^K n=c (mod3^r)`, rotating the target changes its right-hand side to

`10c - 2^K d(10^k-1) (mod3^r)`.

That transports the *integrality congruence*. It does not automatically establish positivity or the exact sequence of binary valuations. The source is a fixed-width digit state; the target is a ternary residue; the lost information is the high carry unless its digit is retained. Primality plays no role in this transport.

## 5. The actual Collatz carry in the decimal family

**PROVED.** For every `k>=2`, the first accelerated odd step of `N_k` grows and has halving exponent 1. For every `k>=4`, its second accelerated odd step is below `N_k`.

Indeed

\[
3N_k+1=10^k-6,
\]

which has valuation 1 for `k>=2`. Hence

\[
y=U(N_k)=(10^k-6)/2>N_k,
\qquad 3y+1=(3\cdot10^k-16)/2.
\]

At `k=4`, this gives `3331 -> 4997 -> 937`, with exponents `(1,4)`. For all `k>=5`, the second numerator has valuation 3, and

\[
U^2(N_k)=3\cdot10^k/16-1,
\]

\[
N_k-U^2(N_k)=(7\cdot10^k-64)/48>0.
\]

The two small counterexamples to this two-step claim are

`31 -> 47 -> 71` and `331 -> 497 -> 373`.

The claim concerns first descent below the starting value, not a universal root certificate for the landing value. It is a useful certificate of a bounded descent block, with clock two accelerated odd steps. Its ordinary Collatz clock is six steps for `k>=5` (`1+1` and `1+3` ordinary steps), and seven at `k=4`.

Primes `N_4,...,N_8` and composite `N_9` all obey this descent mechanism. In particular factor 17 is not the trigger: its first odd step still grows, and the same two-step descent already happened at prime indices. The binary congruence and the ternary order filter are independent coordinates.

**Hostile to a fixed growth horizon.** For every `m>=2`, the composite number

\[
n=2^{2m}-1
\]

is divisible by 3 and has the first `2m-1` accelerated odd steps

\[
U^j(n)=3^j2^{2m-j}-1>n\qquad(1\le j\le2m-1),
\]

each with halving exponent 1. This follows by direct induction while the remaining power of two is at least 2. Thus being composite, even having the fixed small factor 3, does not impose any uniform finite bound on consecutive growth. This does not refute eventual convergence; it identifies the invalid intermediate implication.

## 6. Circular-prime lengths and current source status

**PROVED inheritance.** A decimal circular prime with at least two digits uses only `1,3,7,9`. The statement that prime rotation fixed points are exactly prime repunits requires at least two digits: the single-digit primes `2,3,5,7` are fixed and are not repunits, as explicitly tested here. A prime repunit has prime length, because `R_d` divides `R_k` when `d` divides `k`. This prime-length restriction applies to repunit primes, not all circular primes: `1193` and `193939` have composite digit lengths 4 and 6.

**FINITE-EXACT here.** A sieve of all integers below `10^6`, checked independently by enumerating the alphabet `{1,3,7,9}`, gives 55 circular primes, at lengths 1 through 6. The respective counts are `4,9,12,8,10,12`; only `11` is a repunit among them. Fixed-width rotations are used on both paths.

**CITED current search.** The original research/search reports collected by [Patrick De Geest](https://www.worldofnumbers.com/circular.htm), accessed 2026-09-25, list no nonrepunit circular primes of lengths 7 through 25, with the length-24/25 exclusions credited to Saverio Castelli in October/November 2024. Its entries for lengths 26 and 27 remain unknown. This is a bounded computational exclusion, not a theorem excluding every larger length. Its discussion of finitely many nonrepunit circular primes is explicitly heuristic. Neither a complete finite family nor infinitude is established.

**CITED current primality status.** The known certified decimal repunit lengths currently listed are

`2,19,23,317,1031,49081,86453,109297`.

The first five are inherited classical values. The three later certifications are supported directly by [Underwood's R49081 proof announcement and certificate link](https://t5k.org/primes/page.php?id=133761), [Enge's R86453 research paper](https://arxiv.org/abs/2404.05506), and [Enge and Underwood's R109297 computation](https://enge.math.u-bordeaux.fr/blog/ecpp-109297.html). Enge hosts [independently checkable certificate files](https://www.multiprecision.org/cm/ecpp.html). This session did not rerun those large certificates. The current [repunit table](https://pzktupel.de/Primetables/TableRepunit.php) lists `270343,5794777,8177207` as probable-prime lengths, not certified additions. Repunit-prime infinitude remains conjectural, as stated by [PrimePages](https://t5k.org/glossary/page.php?sort=Repunit).

The inherited September-24 note's five-value list is explicitly restricted to lengths `<=1100` and remains correct; the newer giant certifications above update the wider context without retracting that bounded statement. A separate wording repair in that note is needed: **circular primality** is rotation-invariant by definition; **ordinary primality of a single orbit member** is not. Also, the at-least-two-digit restriction must accompany the prime-fixed-point statement wherever it is summarized. These scope/prose repairs were reported to the coordinator for the shared truth-surface log.

## 7. Reproduction and remaining obligation

Run from the repository root:

```text
python 04-computation/experiments/ternary_digits_20260925.py
```

The standard-library script writes [ternary_digits_20260925.out](ternary_digits_20260925.out). Universe and filters are printed before results. It checks all seven displayed prime candidates by exhaustive trial division; the small factorization certificates; all circular primes below `10^6` by two independent paths; all residues of the conjugacy modulo `3^r`, `1<=r<=8`; valuation identities for `1<=k,t<=60`; decimal Collatz blocks for `2<=k<=300`; and composite rising controls for `2<=m<=100`. The proofs, rather than those finite ranges, support the universal identities.

The positive bridge is a decoder of ternary residue clocks plus an independent binary descent block. To obtain a Collatz-to-Pythagorean root certificate, the next construction must intertwine the actual guarded Collatz moves, retain their binary clock, and transport a well-founded descent rank on the *reachable image*. An abstract isomorphism between infinite ternary trees, a finite residue cycle, or a factor congruence alone does not supply this last obligation.
