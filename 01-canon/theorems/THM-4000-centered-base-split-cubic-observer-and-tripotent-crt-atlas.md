---
id: THM-4000
title: "Six-to-eleven integral observer compiler, split-cubic tripotents, and the centered Farey--Pell ray"
status: >
  PROVED elementary arithmetic + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED.  Consecutive integer samples
  reconstruct P(B) modulo the optimal falling-factorial modulus; the
  six-node B=10 specialization packages the divisibility observers 6 through
  11 and retains an exact-sample factor 12 beyond their residue-only lcm.
  The three-node specialization is the integral split cubic X^3-X: its
  evaluation lattice has cokernel Z/2, its roots modulo every n are classified
  prime-power by prime-power, and for n>=3 the centered roots exhaust exactly
  for odd prime powers and n=4.  The same centered packet gives one parabolic Farey
  ray and an odd-depth Pell square-factor tower.  Power-of-two cyclotomic
  quotients give a separate exact order-8 observer.  The E8 and OCF sections
  are typed applications of stabilization/evaluation, not a Bott, octonion,
  LRC, class-group, or cross-problem theorem.
source: root-2026-08-24-six-eleven-number-machinery
audit: >
  PASS (base_triple_nt + period_wrap_audit, 2026-08-24).  The sampler and
  determinantal-divisor Smith forms, split-cubic parity, tripotent prime-power
  classification including the 2-adic hostile, Farey/Berggren indexing, Pell
  factorization, cyclotomic exact orders, common-shell E8 normalization, OCF
  kernel, Moon--Busch optimizer, quantifiers, and scope firewalls were
  independently reconstructed.  Normal, optimized, and frozen outputs match.
depends_on:
  - THM-3555-catalan-thickening-universal-cubic-root-cover
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3380-hamiltonian-deletion-layer-monoid-semiring-and-small-order-boundaries
  - THM-1370-h-spectrum-omits-7-21-all-n
related:
  - THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase
  - THM-872-milgram-discriminant-residue
script: 04-computation/six_to_eleven_number_machinery_thm4000.py
output: 05-knowledge/results/six_to_eleven_number_machinery_thm4000.out
script_sha256: 4777f23c31f3dbe65f420b6a1b77f303677e94015c009ffd113bcbb0ce88cd73
output_sha256: 1c482f100f60496ab1d678f75687d52b3e33b6eaa6e4e9096948479c87d27072
hash_basis: raw LF bytes
---

# THM-4000 -- the six-to-eleven integral observer compiler

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The integers

```text
6,7,8,9,10,11
```

are not being identified by size or symbolism.  They form one exact packet
because they are the six distances from the decimal base `10` to the six
consecutive observer nodes

```text
4,3,2,1,0,-1.
```

The three central observers `1,0,-1` recover the base-neighbour packet
`9,10,11`; at base seven they recover `6,7,8`.  The governing integral
algebra is the split cubic `X^3-X`, already present as the completely split
fiber in
[THM-3555](THM-3555-catalan-thickening-universal-cubic-root-cover.md).
The new content here is its arithmetic observer, CRT, tripotent, Farey, and
Pell interfaces.

## 1. Consecutive samples give an optimal congruence compiler

Fix integers `a`, `m>=0`, and put

```text
A_(a,m)={a,a+1,...,a+m},
F_(a,m)(X)=product_(j=0)^m (X-a-j).                     (1)
```

For `P in Z[X]`, let `Delta` be the forward-difference operator and define

```text
R_(a,m)P(X)=sum_(j=0)^m Delta^j P(a) binom(X-a,j).      (2)
```

Here `binom(t,j)=t(t-1)...(t-j+1)/j!` is the integer-valued polynomial
binomial for every `t in Z`, including negative `t`.

Although the displayed binomial basis has rational coefficients, `(2)` is
integral: the falling factorials form a monic `Z`-basis of `Z[X]`, so
`Delta^j P(a)/j!` is an integer.  Thus `(2)` is the ordinary integral
remainder of `P` on division by the monic polynomial `F_(a,m)`.  It has degree
at most `m` and agrees with `P` at every point of `A_(a,m)`; hence

```text
P(X)-R_(a,m)P(X) in F_(a,m)(X) Z[X].                   (3)
```

Consequently, for every integer `B` outside `A_(a,m)`,

```text
P(B) == sum_(j=0)^m binom(B-a,j) Delta^j P(a)
        mod product_(j=0)^m (B-a-j).                   (4)
```

where the modulus means the absolute value of the displayed product.  At a
sample node the statement is instead exact equality.  The modulus in `(4)` is
optimal.  The ideal of values at `B` of integer
polynomials vanishing at all sample nodes is exactly

```text
F_(a,m)(B) Z,                                          (5)
```

because the kernel is the principal ideal `(F_(a,m))` and `P=F_(a,m)`
attains its generator.

### The factorial Smith sidecar

Evaluation on the consecutive nodes is not onto `Z^(m+1)`.  In the monic
falling-factorial basis

```text
f_k(X)=product_(r=0)^(k-1)(X-a-r),
```

one has

```text
f_k(a+j)=j^(under k)=k! binom(j,k).                     (6)
```

The Pascal matrix `(binom(j,k))` is integral unitriangular, and the change
from monomials to the `f_k` is also integral unitriangular.  Since
`0!|1!|...|m!`, the evaluation lattice therefore has Smith form

```text
SNF(ev_(a,m))=diag(0!,1!,2!,...,m!).                   (7)
```

These factorial invariant factors are the exact compatibility/carry
sidecars among integer samples.  They explain why exact sample lifts can
recover repeated prime powers which separate residue congruences lose.

## 2. One six-observer formula packages every modulus from 6 to 11

Take `a=-1`, `m=5`, `B=10`.  Then

```text
F_(-1,5)(10)=11*10*9*8*7*6=332640.                    (8)
```

Expanding `(4)` in the sample values gives the exact congruence

```text
P(10) == -252 P(-1)+1386 P(0)-3080 P(1)
          +3465 P(2)-1980 P(3)+462 P(4)     (mod 332640). (9)
```

Equivalently, the forward-difference weights are

```text
1,11,55,165,330,462 = binom(11,0),...,binom(11,5).     (10)
```

For each `r in {6,...,11}`, put `a_r=10-r`.  Polynomial congruence gives

```text
P(10)==P(a_r) (mod r).                                  (11)
```

Thus `(9)` simultaneously compiles the six ordinary divisibility observers

```text
r:       11  10   9   8   7   6
a_r:     -1   0   1   2   3   4.                       (12)
```

For a decimal digit polynomial, `P(1)` is the digit sum, `P(0)` the final
digit, and `P(-1)` the alternating digit sum; `P(2),P(3),P(4)` are the
additional small-base observers needed to include `6,7,8` in the same packet.

The six residue classes in `(11)` alone determine `P(10)` only modulo

```text
lcm(6,7,8,9,10,11)=27720.                              (13)
```

The exact samples in `(9)` determine it modulo their product `332640`, a
factor `12` more.  This is not a contradiction to CRT: `(9)` consumes exact
integer lifts satisfying the Smith constraints

```text
diag(1,1,2,6,24,120),                                  (14)
```

whereas `(11)` has discarded those lifts.

## 3. The centered three-node compiler is the split cubic

Set `a=-1,m=2`.  For

```text
u_-=P(-1),       u_0=P(0),       u_+=P(1),
```

endpoint parity gives `u_-==u_+ (mod 2)`, and the unique integral quadratic
remainder is

```text
r_P(X)=u_0+((u_+-u_-)/2)X+((u_++u_--2u_0)/2)X^2.       (15)
```

The evaluation matrix on `(1,X,X^2)` is

```text
[1 -1 1]
[1  0 0]
[1  1 1],                                               (16)
```

with determinant `2` and Smith form `(1,1,2)`.  Therefore

```text
Z[X]/(X^3-X) -> Z^3,
[P] |-> (P(-1),P(0),P(1))                              (17)
```

is injective with image

```text
{(u_-,u_0,u_+):u_-==u_+ (mod 2)}                       (18)
```

and cokernel `Z/2`.  Away from the prime `2`, the three split roots really
give a product of three copies; integrally, the endpoint collision leaves
one parity invoice.

Evaluating `(15)` at an integer base `b>=2` proves

```text
P(b) == [b(b+1)/2]P(1)+(1-b^2)P(0)+[b(b-1)/2]P(-1)
        (mod b^3-b).                                    (19)
```

Again the modulus is optimal, attained by `P=X^3-X`.

### The two requested centered packets

At `b=7`,

```text
(b-1,b,b+1)=(6,7,8),
P(7)==28P(1)-48P(0)+21P(-1) (mod 336).                 (20)
```

The separate neighbour residues recover only modulo

```text
lcm(6,7,8)=168;                                         (21)
```

An explicit compatible-residue compiler is

```text
P(7)==112P(1)+120P(0)+105P(-1) (mod 168).              (21a)
```

Its difference from `(20)` is
`84(P(1)+P(-1))+168P(0)`, so endpoint parity makes it vanish modulo `168`.
The exact lifts in `(20)` retain one further factor two.

At `b=10`,

```text
(b-1,b,b+1)=(9,10,11),
P(10)==55P(1)-99P(0)+45P(-1) (mod 990).                (22)
```

Here the three moduli are pairwise coprime, so ordinary CRT already reaches
the full product.  Its idempotent form is

```text
P(10)==550P(1)-99P(0)+540P(-1) (mod 990).              (23)
```

The difference between `(22)` and `(23)` is
`495(P(1)+P(-1))`, divisible by `990` by endpoint parity.

In general

```text
gcd(b,b-1)=gcd(b,b+1)=1,
gcd(b-1,b+1)=gcd(b-1,2),                               (24)
lcm(b-1,b,b+1)=
  b(b^2-1)       if b is even,
  b(b^2-1)/2     if b is odd.                          (25)
```

These are the specialization shadows of

```text
Res(X,X-1)=-1,  Res(X,X+1)=1,  Res(X-1,X+1)=2,
disc(X^3-X)=4.                                         (26)
```

The endpoints are also the first two cyclotomic values

```text
b-1=Phi_1(b),       b+1=Phi_2(b).                       (27)
```

Thus primes in the left endpoint see the base with order one, primes in the
right endpoint see `b==-1`; for odd endpoint primes this has exact
multiplicative order two, while modulo two it has order one.  Primes in the
center see the base as zero.  This is the exact meaning of the `1,0,-1`
observer; it is not digit symbolism.

## 4. Roots of `x^3=x`: the tripotent atlas

Let

```text
Trip(n)={x in Z/nZ:x^3=x}.                              (28)
```

CRT makes its cardinality multiplicative.  For an odd prime `p`, the three
roots `0,+1,-1` modulo `p` are simple, since the derivative `3X^2-1` is
`-1,2,2`; Hensel lifting therefore gives exactly three roots modulo every
`p^e`.  At the prime two,

```text
|Trip(2)|=2,       |Trip(4)|=3,       |Trip(2^e)|=5 for e>=3. (29)
```

For `e>=3`, an even tripotent must be zero; the odd ones are exactly the four
solutions of `x^2=1 (mod 2^e)`.

Writing `n=2^e product_(i=1)^r p_i^(a_i)` with distinct odd `p_i`, one gets

```text
|Trip(n)|=c_2(e) 3^r,
c_2(0)=1, c_2(1)=2, c_2(2)=3, c_2(e>=3)=5.             (30)
```

For `n>=3`, the centered three roots exhaust `Trip(n)` exactly when

```text
n is an odd prime power, or n=4.                        (31)
```

In any commutative ring a tripotent has an idempotent support and is a corner
involution.  If `x^3=x`, then

```text
e=x^2,       e^2=e,       ex=x,       x^2=e;            (32)
```

so `x` is an involutive unit in the corner ring `eR`.  Conversely `(32)`
implies `x^3=x`.  When two is invertible, this refines to the orthogonal
signed-idempotent decomposition

```text
e_+=(x^2+x)/2,       e_-=(x^2-x)/2,
e_+e_-=0,            x=e_+-e_-,       x^2=e_++e_-.     (32a)
```

Without that hypothesis, “signed idempotent” is false: `x=3 mod 8` is a
tripotent with support one, but it is neither `+1` nor `-1`.  The universal
coordinates are the support idempotent and its corner involution; local signs
are valid on odd-primary factors.

The requested moduli display every boundary:

```text
n   Trip(n)
6   {0,1,2,3,4,5}
7   {0,1,6}
8   {0,1,3,5,7}
9   {0,1,8}
10  {0,1,4,5,6,9}
11  {0,1,10}.                                           (33)
```

Thus `7,9,11` have the pure sign trichotomy; `8` has the four 2-adic unit
signs; and composite `10` has three additional mixed CRT states.  The
centered decimal triple `9,10,11` gives the diagonal `-1,0,+1`, not the full
tripotent algebra modulo ten.

There is also a distinct role for six:

```text
gcd_{z in Z}(z^3-z)=6.                                  (34)
```

Three consecutive integers make divisibility by `2` and `3` immediate, and
`2^3-2=6` proves maximality.  Hence

```text
x^3==x (mod n) for every residue x  iff  n divides 6.  (35)
```

The primes in `(34)` have different mechanisms.  The prime `2` is the
discriminant/resultant collision in `(26)`; the prime `3` occurs because the
three simple roots exhaust the three-element field.  They must not be
merged into one vague notion of ramification.

## 5. The same centered triple is a parabolic Farey edge

For an integer `b>=2`, define

```text
A_b=[b   b-1; b+1  b],          C=[0 1; -1 2].          (36)
```

Then

```text
det A_b=1,              C A_b=A_(b+1),                 (37)
(b-1)/b < b/(b+1),
b/(b+1)-(b-1)/b=1/[b(b+1)].                             (38)
```

The matrix `C` is parabolic and fixes the cusp `1`.  Therefore

```text
6/7<7/8       and       9/10<10/11                     (39)
```

are literally two edges of one parabolic cusp ray, separated by three
applications of `C`.  After reversing spinor coordinates, this is the
pre-Gaussian consecutive-spinor carrier of
[THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md),
not a new analogy.

The Gaussian-square map makes the connection operational:

```text
(b,b-1) |-> (b^2-(b-1)^2, 2b(b-1), b^2+(b-1)^2)
         =(2b-1,2b(b-1),2b(b-1)+1).                    (40)
```

Thus `b=7` gives `(13,84,85)` and `b=10` gives `(19,180,181)`, two
consecutive-coordinate states on the proved Berggren `U`-spine, three
`U`-steps apart.  The determinant/Farey carrier is
preserved; digit values, tripotent states, and divisibility tests are not.

## 6. Fixed-base powers produce a Pell square-factor tower

For an integer `b>=2`, put `d=b^2-1`, and let `T_n=T_n(b)`, `U_n=U_n(b)` be the Chebyshev
polynomials.  Since

```text
A_b=bI+[0 b-1; b+1 0],
[0 b-1; b+1 0]^2=dI,                                  (41)
```

one has for every `n>=1`

```text
A_b^n=[T_n, (b-1)U_(n-1); (b+1)U_(n-1), T_n],          (42)
T_n^2-d U_(n-1)^2=1.                                   (43)
```

The top-to-bottom projective ratios of the two columns bracket the quadratic irrational
`sqrt((b-1)/(b+1))`:

```text
(b-1)U_(n-1)/T_n
  < sqrt((b-1)/(b+1))
  < T_n/((b+1)U_(n-1)),                                (44)
```

with exact interval width

```text
1/[(b+1)T_n U_(n-1)].                                  (45)
```

For odd `n=2m+1`, the Chebyshev half-angle factorization gives

```text
T_n-1=(b-1)(U_m+U_(m-1))^2,
T_n+1=(b+1)(U_m-U_(m-1))^2,                            (46)
```

where `U_(-1)=0`.  Hence each centered base triple seeds infinitely many
common-center square-scaled endpoints.

For the decimal packet,

```text
X=T_(2m+1)(10):  10,3970,1580050,...,
X-1=9 A_m^2,              X+1=11 B_m^2.                (47)
```

For the base-seven packet the first centers are

```text
7,1351,262087,...,                                      (48)
```

with `X-1=6A^2`, `X+1=8B^2`.  Translation in the base
(`C A_b=A_(b+1)`) and powering at fixed base are different operations; the
theorem records both rather than conflating them.

## 7. A genuine arithmetic period-eight observer

For integers `b>=2` and `k>=1`,

```text
Phi_(2^k)(b)=b^(2^(k-1))+1.                             (49)
```

Modulo this integer, `b^(2^(k-1))=-1`; since the modulus is greater than
two, the order of `b` is exactly `2^k`.  Every odd prime divisor `q` also
has

```text
ord_q(b)=2^k,                 q==1 (mod 2^k).           (50)
```

For `k=3`, polynomial division by `X^4+1` gives the exact eight-phase fold

```text
sum_i a_i X^i |->
sum_(r=0)^3 [sum_(j>=0)(-1)^j a_(r+4j)] X^r.           (51)
```

Thus `P(b) mod b^4+1` is the alternating sum of four-digit base-`b` blocks.
The modulus is optimal, attained by `P=X^4+1`.  The two bases give

```text
Phi_8(7)=2402=2*1201,          ord_1201(7)=8,
Phi_8(10)=10001=73*137,        ord_73(10)=ord_137(10)=8. (52)
```

Equations `(49)`--`(50)` also give the elementary Euclid proof of infinitely
many primes `q==1 (mod 2^k)`: if a finite list existed, choose an even base
divisible by their product and take an odd prime divisor of `(49)`.

This is an arithmetic period-eight clock in a cyclotomic quotient.  It is
not Bott periodicity: varying `b` changes its modulus and prime splitting
without changing any Clifford class.

## 8. A second, independent arithmetic eight-shift: `L -> L direct_sum E8`

For a positive-definite even integral lattice `L`, use the half-squared-norm convention
`Theta_L(q)=sum_(v in L) q^((v,v)/2)`.  Orthogonal stabilization gives

```text
discform(L direct_sum E8)=discform(L),
Theta_(L direct_sum E8)=Theta_L Theta_E8=Theta_L E_4,   (53)
```

because `E8` is even unimodular and norms add.  This is a genuine typed
rank-`r` to rank-`r+8` number-theoretic operation.  Applied to the root
lattices `A_1,A_2,A_3`, it gives

```text
L          rank   |disc|  roots     L direct_sum E8 rank  roots
A1           1       2      2                9             242
A2           2       3      6               10             246
A3           3       4     12               11             252. (54)
```

So ranks `9,10,11` are exact `E8`-stabilizations of ranks `1,2,3`, with
their discriminant forms preserved.  They are not all copies of rank one,
and they are not the affine/hyperbolic objects informally called
`E9,E10,E11`.

[MISTAKE-484](../MISTAKES.md) records the decisive hostile: the historical
identity `theta_3^8=Theta_E8` was false.  In the common squared-norm
convention, `Theta_(Z^8)=theta_3(q)^8` but `Theta_E8=E_4(q^2)`.  Their norm-two
shells have sizes `112` and `240`, and their norm-four shells have sizes
`1136` and `2160`.  Stabilization multiplies theta series; it does not identify
them.

## 9. An unrelated repo application: the OCF evaluation kernel

The same discipline -- retain a vector before evaluating it at one scalar --
sharpens [THM-3380](THM-3380-hamiltonian-deletion-layer-monoid-semiring-and-small-order-boundaries.md).
For covered size `s`, group odd-cycle packing types by their number `ell` of
cycles, and let `L_s` be the attained lengths.  Fugacity specialization
`x=2` is the integer map

```text
E_(2,s):Z^(L_s)->Z,             e_ell |-> 2^ell.        (55)
```

If `ell_0<...<ell_t`, then

```text
im E_(2,s)=2^(ell_0) Z,
ker E_(2,s)=direct_sum_(i=1)^t
  Z(2^(ell_i-ell_0)e_(ell_0)-e_(ell_i)).               (56)
```

This follows by solving the first coordinate after dividing by
`2^(ell_0)`; the displayed vectors form a primitive basis.

For covered sizes six through eleven:

```text
s    odd-part packing types          L_s
6    3+3                             {2}
7    7                               {1}
8    3+5                             {2}
9    9, 3+3+3                        {1,3}
10   3+7, 5+5                        {2}
11   11, 3+3+5                       {1,3}.              (57)
```

The first scalarization kernel is therefore

```text
4e_1-e_3                                                   (58)
```

at covered size nine; it appears again at eleven, while ten has no
length-grading kernel.  THM-3380's exact order-nine witnesses realize `(58)`
as

```text
(a_9,b_333)=(161,9) <-> (157,10).                       (59)
```

Four Hamiltonian 9-cycles trade for one `3+3+3` packing under evaluation at
two.  It follows immediately from `(57)` that at order eleven

```text
(Dham,b_333,b_335) determines the full bivariate F_T(x,y). (60)
```

Necessity of `b_333` is witnessed (and persists under ordered-join singleton
padding).  Whether a realizable order-eleven pair makes `b_335` necessary is
**OPEN**.  This OCF result shares an evaluation-loss mechanism with the
integer observers; it shares no number-theoretic object and gives no LRC or
octonion transfer.

## 10. A second unrelated application: Moon--Busch is min-plus arithmetic

[THM-1370](THM-1370-h-spectrum-omits-7-21-all-n.md) imports the cited
Moon--Busch strong-tournament floor

```text
f(n)=min {3^a 5^b : a,b in Z_(>=0), 2a+3b=n-1}, n>=3. (61)
```

This is an exact optimization over the numerical semigroup generated by two
and three.  Replacing `(a,b)` by `(a-3,b+2)` preserves `2a+3b` and multiplies
the objective by `25/27<1`.  Hence the minimizer uses the largest

```text
b<=floor((n-1)/3) with b==n-1 (mod 2),                 (62)
```

and then `a=(n-1-3b)/2`.  In the requested window this gives

```text
n          6    7    8    9    10   11
(a,b)    (1,1)(0,2)(2,1)(1,2)(0,3)(2,2)
f(n)      15   25   45   75   125  225.                (63)
```

The input to this arithmetic observer is the order `n`; the target is a
weighted lattice point on `2a+3b=n-1`, and the preserved output is the minimum
Hamiltonian-path count.  The optimizer records parameters of an attaining
Moon construction.  The value formula alone does not classify every extremal
tournament or its block placement.  The historical quadratic and
Fibonacci-style fits are hostile controls: both agree for several terms and
then fail, as recorded in MISTAKE-055.

## 11. Verification, boundaries, and next machinery

The standard-library companion performs `148241` exact gates, as frozen in its
output.  It checks:

- consecutive-sampler cases on both sides of the node set and the Smith
  diagonals recovered from determinantal divisors;
- `10000` decimal digit-polynomial controls for `(9)`;
- the split-cubic congruence, CRT idempotents, and every tripotent modulo
  `2<=n<=512` against `(30)` and `(31)`;
- all Farey/Pell identities for `2<=b<=30`, `1<=n<=15`;
- power-of-two cyclotomic folds for `2<=b<=30`, orders through `32`;
- the `A_1,A_2,A_3` theta/root/discriminant stabilization rows;
- every odd-part packing type at covered sizes `6,...,11`; and
- the Moon--Busch min-plus optimizers at orders `6,...,11`.

The positive controls are bases seven and ten.  The hostiles are load-bearing:

```text
odd base b=7: neighbour residues lose one factor 2;
mod 8 and mod 10: centered tripotents do not exhaust;
covered size 9: x=2 loses 4e_1-e_3;
Z^8 versus E8: equal rank does not imply equal theta series. (64)
```

The natural next theorem is the higher-jet sampler.  Reconstructing modulo
`F_(a,m)(B)^k` requires Hasse derivatives through order `k-1` at every node;
the integral Smith/carry lattice of that confluent evaluation, rather than a
formal Taylor slogan, is the missing coordinate.  A second exact target is a
realizable order-eleven witness for the `b_335` kernel in `(60)`.

Nothing here upgrades the almost-complex structure on `S^6`, transfers class
group rank certificates, proves LRC(14), or identifies the cyclotomic clock
with Bott periodicity.  Those objects require their own maps and predicates.
