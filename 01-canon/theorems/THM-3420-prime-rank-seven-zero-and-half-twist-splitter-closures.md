---
id: THM-3420
title: "Prime rank-seven zero- and half-twist splitter closures"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  On a prime number p of
  sheets, a fixed-zero cover by at most seven strict danger blocks exists if
  and only if p=29.  Among primes p congruent to 13 modulo 14, a half-twist
  cover by at most seven blocks exists if and only if p=13.  Both positive
  cases have exact rank seven in their stated layer.  The theorem also gives
  the exact rank-seven capacity defects by order modulo 14 and reduced-
  numerator parity.  It does not classify the other half-twist prime classes,
  composite rank seven, arbitrary physical common times, or LRC(14).
source: root-rank7-reflection-2026-08-15
audit: independent reconstruction of both proofs and arithmetic CLEAN; standard-library exact companion; normal/optimized replay; Lucas certificate for the sole factor above 64 bits; hostile compatibility graphs
depends_on:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
related:
  - THM-3414-fixed-zero-six-owner-base-classification
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3407-hadamard-multi-toggle-response-and-sharp-trade-distance
script: 04-computation/lrc_prime_rank7_splitter_closures_thm3420.py
output: 05-knowledge/results/lrc_prime_rank7_splitter_closures_thm3420.out
script_sha256: f95aaf081ccd5c92cb7474104f242623bba9246bdc19eb04ce8da81f3c4e2af6
output_sha256: 22740907984a9e64e75208c0a01ce0222fd844dfe65826685ff11e64274ef959
semantic_sha256: 4dffcebc010848fc022824c46a916208d8b46c70ac3306168463b613ce58a7ed
hash_basis: LF-normalized bytes
---

# THM-3420 -- prime rank-seven zero- and half-twist splitter closures

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and exact scope

Retain the zero-mode-cochain layers of
[THM-3405](THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md).
For a prime `p`, a nonzero residue `r` defines the two strict sheet masks

```text
Z_(p,r)={ell in Z/pZ: ||r ell/p||<1/14},
H_(p,r)={ell in Z/pZ: ||r(2ell+1)/(2p)||<1/14}.        (1)
```

Then:

```text
some at-most-seven fixed-zero masks Z_(p,r) cover Z/pZ
iff p=29;                                                (2)

for prime p==13 (mod 14),
some at-most-seven half masks H_(p,r) cover Z/pZ
iff p=13.                                                (3)
```

The positive atoms are

```text
p=29, zero: residues (1,4,5,6,7,9,13),
  block sizes (5,5,5,5,5,5,5),
  multiplicity 7 at sheet 0 and 1 elsewhere;

p=13, half: residues (1,2,3,5,7,9,11),
  block sizes (2,1,2,2,2,2,2),
  multiplicity 1 on every sheet.                       (4)
```

The same capacity arguments show that neither atom can use six blocks in its
stated layer.  With THM-3416's already proved support through rank six, they
also give `rho_ZMC(13)=rho_ZMC(29)=7`.  No converse for all prime or all
composite sheet numbers is asserted.

## 2. Inheritance and loss ledger

THM-3405 reduces every primitive zero-complete-cochain cover to the literal
fixed-zero or Boolean half-twist layer.  THM-3416 succeeds at rank six because
one block must have density strictly greater than the limiting `1/7`; only
finitely many quotient orders survive.  Seven blocks sit exactly at that
limit, so the order cutoff disappears.

| field | exact connection |
|---|---|
| source | primitive fixed-zero and half-twist sheet masks |
| target | reflection-orbit capacities and multiplicative translate packings |
| map | retain quotient order and numerator parity, then identify a prime-sheet mask with a dilate of a short symmetric multiplier set |
| preserved | strict endpoints, block size, zero/fixed sheet, reflection, cover multiplicity, and primitive parity |
| destroyed by the first map | the multiplier ratio and hence exact translate intersection |
| required sidecar | the ratio set `A/A`, equivalently the splitter/factorization incidence |
| cheapest decisive controls | zero `p=29` versus `p=43`, and half `p=13` versus `p=41` |

The theorem restores exactly that multiplier sidecar.  It does not turn a
zero-mode-cochain cover into an arbitrary physical-time cover.

## 3. Exact rank-seven capacity and reflection defects

For a half-twist owner of quotient order `m`, write its reduced numerator as
`s`.  If `s` is odd, the block count is

```text
a(m)=2 ceil(floor((m-1)/7)/2).                         (5)
```

If `s` is even, then `m` must be odd and the count is

```text
z(m)=1+2 floor((m-1)/14).                              (6)
```

The fixed-zero count is also `z(m)`, independently of numerator parity.
Put

```text
delta_z(m)=7z(m)-m,       delta_a(m)=7a(m)-m.          (7)
```

As `r=m mod 14` runs from zero through thirteen, the exact tables are

```text
r:          0   1  2  3  4  5  6  7   8   9  10  11  12  13
delta_z:   -7   6  5  4  3  2  1  0  -1  -2  -3  -4  -5  -6
delta_a:    0  -1 -2 -3 -4 -5 -6 -7   6   5   4   3   2   1.  (8)
```

The half mask is invariant under `sigma(ell)=-1-ell`.  For odd `m`, an
even-reduced block contains the unique fixed point of `sigma`, while an
odd-reduced block avoids it.  For even `m`, `sigma` has no fixed point and
only the odd-reduced type occurs.  Thus `(5)--(8)` record the exact number of
reflection pairs as well as the scalar size.

For seven blocks of orders `m_i|Q`, types `t_i in {z,a}`, and pullback degrees
`d_i=Q/m_i`, the total overlap mass of a cover is

```text
Omega=sum_i |B_i|-Q=(1/7) sum_i d_i delta_(t_i)(m_i)>=0.  (9)
```

Equality holds exactly for a partition.  At zero twist all blocks contain
sheet zero, so in fact `Omega>=6`.  At half twist with odd `Q`, if `e` blocks
contain the fixed sheet then `Omega>=e-1`, and all remaining overlap occurs in
reflection pairs.  In particular

```text
Omega == e-1 (mod 2).                                  (10)
```

For even `Q`, there is no fixed sheet and `Omega` is even.

These exact gates still leave infinitely many formal profiles:

- for every `m==1 (mod 14)`, seven full-order zero blocks have total size
  `m+6`, exactly paying the forced common-zero tax;
- for every odd `m==13 (mod 14)`, six odd-reduced half blocks and one
  even-reduced half block have total size exactly `m` and the right fixed-point
  invoice.

All seven blocks can live on the same unbounded order, so a CRT product has no
independent second cylinder.  Reflection plus CRT therefore does **not**
restore the rank-six finite-order compiler.  Exact translate incidence is the
first missing coordinate.

## 4. Fixed-zero primes reduce to one multiplicative splitter

Let `p` be prime and put

```text
k=floor((p-1)/14),       I_k={0,+-1,...,+-k} in F_p.   (11)
```

Every nonzero fixed-zero mask is `r^(-1)I_k`.  All masks contain zero, so the
union of at most seven masks has size at most

```text
1+7(2k)=14k+1.                                        (12)
```

For `p<=7`, one has `k=0`, so every block is the singleton zero and no cover
exists.  For `p>7`, a prime residue modulo fourteen belongs to
`{1,3,5,9,11,13}`.  Formula `(12)` forces

```text
p=14k+1                                                (13)
```

and exactly seven blocks.  Equality in `(12)` then forces their nonzero parts
to be pairwise disjoint.  If the residues are `r_1,...,r_7`, set

```text
A={1,...,k},       C={+-r_i^(-1):1<=i<=7}.             (14)
```

The cover is equivalent to the unique multiplicative factorization

```text
F_p^* = A C,              |A|=k, |C|=14.              (15)
```

This is the precise splitter-set carrier that the capacity profile forgot.

### 4.1 Power sums leave only `29` and `211`

Write

```text
S_j(k)=sum_(a=1)^k a^j,       T_j(C)=sum_(c in C)c^j. (16)
```

For `1<=j<=14<p-1`, the unique factorization `(15)` gives

```text
0=sum_(x in F_p^*)x^j=S_j(k)T_j(C).                   (17)
```

Since `C=-C`, every odd `T_j` vanishes.  If all seven even power sums through
degree fourteen vanished as well, Newton identities would force every
elementary symmetric function of the fourteen nonzero elements of `C` to
vanish.  The degree-fourteen one is their nonzero product, a contradiction.
Therefore `S_j(k)=0 mod p` for some even `j<=14`.

From `(13)`, `k=-1/14 mod p`.  Exact Faulhaber evaluation gives the following
numerator factorizations; every denominator is a power of two times a power of
seven and hence is a `p`-unit:

```text
j= 2: 13
j= 4: 13*47
j= 6: 13*46567
j= 8: 13*17*131*821
j=10: 13*19*41*1045859
j=12: 211*32873822321
j=14: 13*29*16622136190229.                           (18)
```

Among the prime factors in `(18)`, only `29` and `211` are congruent to one
modulo fourteen.  Thus `(15)` forces

```text
p in {29,211}.                                        (19)
```

### 4.2 The `211` candidate is maximally hostile

At `p=211`, `k=15`.  Let `J={+-1,...,+-15}`.  For any
`lambda in F_211^*`, place

```text
0, lambda, 2lambda, ..., 15lambda
```

on the cyclic `211`-gon.  One of the sixteen circular gaps has length at most
`floor(211/16)=13`.  Consequently

```text
d lambda=+-a  (mod 211),       1<=d<=15, 1<=a<=13,    (20)
```

so `lambda in J/J`.  Hence `J/J=F_211^*`: every two fixed-zero blocks meet
away from zero.  Seven such blocks cannot realize the disjointness forced by
`(12)`.  This excludes `211` without an enumeration.

### 4.3 The `29` atom is a Paley splitter

For the residues in `(4)`, the set

```text
C={+-r^(-1)}={nonzero quadratic residues modulo 29}.  (21)
```

The element two is a nonsquare modulo `29`, so

```text
{1,2} C=F_29^*                                        (22)
```

disjointly.  This proves the positive cover conceptually and connects the
rank-seven atom to the Paley/Hadamard quadratic-residue carrier, rather than
merely replaying a seven-mask search.

## 5. The critical half-twist prime class

Now let

```text
p=14k+13                                               (23)
```

be prime.  An odd-residue half block has size `2k+2` and avoids the fixed
sheet; an even-residue block has size `2k+1` and contains it.  A cover needs
at least one even block.  Seven blocks have enough total capacity only when
there is exactly one even block and six odd blocks, in which case their total
size is exactly `p`.  Thus any cover is a partition.

Represent sheets by odd words `x=2ell+1 mod 2p` and remove the fixed word
`x=p`.  After multiplying all nonfixed words by one common unit, the even
block becomes

```text
E={+-1,...,+-k} subset F_p^*,                          (24)
```

while the six odd blocks are dilates `c_i O` of

```text
O={+-1,+-3,...,+-(2k+1)}.                              (25)
```

Therefore

```text
F_p^* = E disjoint_union c_1 O disjoint_union ...
                    disjoint_union c_6 O.             (26)
```

### 5.1 Newton obstruction

Put `y_i=c_i^2` and, for `n>=1`,

```text
S_(2n)(k)=sum_(a=1)^k a^(2n),
R_(2n)(k)=sum_(t=0)^k (2t+1)^(2n),
q_n=sum_i y_i^n.                                      (27)
```

For `p>13`, the exponents `2,...,14` are below `p-1`.  Power-summing the
partition `(26)` gives

```text
q_n=-S_(2n)(k)/R_(2n)(k),       1<=n<=7.              (28)
```

The first six values determine the elementary symmetric functions
`e_1,...,e_6` of the six elements `y_i` by Newton identities.  They must then
satisfy

```text
q_7=e_1q_6-e_2q_5+e_3q_4-e_4q_3+e_5q_2-e_6q_1.       (29)
```

Substitute `k=-13/14 mod p` from `(23)` and let `D` be the left side minus the
right side of `(29)`.  Exact rational simplification gives

```text
numerator(D)=
  7^10 * 13^2 * 277 * 2719 * 7779713652980688586832393,

denominator(D)=
  2^39 * 17^3 * 31 * 61 * 83 * 311 * 1487 * 4597^2
       * 84631 * 255443 * 49785481.                    (30)
```

The large numerator factor in `(30)` is prime: the exact companion verifies
the complete Lucas certificate

```text
N-1=2^3*7^2*13*191*619*643*4649*4319561459,
base=3.                                                (31)
```

Every factor in `(31)` is independently checked prime; `3^(N-1)=1 mod N`,
and the Lucas gcd is one for every distinct prime factor of `N-1`.

No numerator prime in `(30)` other than `13` is congruent to thirteen modulo
fourteen.  The only denominator primes with that residue are `83` and
`255443`.  At either exceptional prime, already

```text
R_14(k)=0 mod p,       S_14(k)!=0 mod p,               (32)
```

so the original power-sum equation `(28)` is impossible.  Therefore no prime
`p>13` in the class `(23)` supports a cover.  The `p=13` partition in `(4)`
completes the proof of `(3)`.

## 6. What remains at rank seven

The same fixed-point invoice sharply narrows, but does not close, the other
prime half-twist classes.  If `e` is the number of even-residue blocks, then

```text
p=14k+r, r in {1,3,5}:    Omega=e-r;
p=14k+r, r in {9,11,13}:  Omega=14-r-e.                (33)
```

Since the fixed sheet alone costs `e-1`, residues `3` and `5` are impossible.
At residue `1`, every nonfixed support must partition.  At residues `9,11,13`,

```text
e <= (15-r)/2.                                         (34)
```

Section 5 closes the equality endpoint `r=13,e=1`.  The residue-one mixed
splitters and the positive-overlap residue-nine/eleven cases remain open.

Composite rank seven also remains open.  Exact positive controls at
`Q=14,38,51,68,148` use nested quotient orders, so the coprime CRT-cylinder
product from the rank-six proof is not available.  No finite primitive basis,
support density, LRC row exclusion, or physical-time conclusion is inferred.

The compatibility graphs provide hostile controls rather than dependencies:

```text
fixed-zero clique number:
  p=29,43,71,211 -> 7,6,4,1;

half, maximum odd-block clique avoiding one even block:
  p=13,41,83 -> 6,5,3.                                (35)
```

## 7. Exact companion

Run

```bash
python3 04-computation/lrc_prime_rank7_splitter_closures_thm3420.py
python3 -O 04-computation/lrc_prime_rank7_splitter_closures_thm3420.py
```

The normal and optimized transcripts are byte-identical.  The companion uses
only the Python standard library.  It independently interpolates all power-sum
polynomials over `Fraction`, replays Newton identities, checks the factor
tables and Lucas certificate, verifies the Paley splitter, proves
`J/J=F_211^*` exactly, computes the hostile clique numbers, and replays the
rank-seven positive overlap profiles.
