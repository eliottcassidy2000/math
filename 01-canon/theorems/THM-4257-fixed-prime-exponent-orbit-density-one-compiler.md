---
id: THM-4257
title: "Fixed-prime exponent-orbit density-one compiler"
status: >
  PROVED RELATIVE TO THM-4244/4250 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  For every odd core b>=3 and fixed odd integer P>1, suffix-good exponent
  classes have monotone orbit density tending to one with an explicit
  exponential defect bound. The theorem gives the exact orbit iff, the
  one/two-sheet lift, collar and coset criteria, and new prime-power factorial
  families for multipliers 18, 22, and 26. Finite suffix failure remains
  nonconclusive, and no converse to THM-4244 is asserted.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4244-even-multiplier-odd-core-complementary-pair-factorial-compiler
  - THM-4250-odd-core-suffix-lift-automaton-and-density-one-compiler
related:
  - THM-3474-factorial-binary-submask-polygon-and-prime-power-reset-families
  - THM-4237-multiplier-six-binary-adjacency-prime-power-factorial-closure
  - THM-4243-multiplier-ten-double-overlap-prime-power-factorial-closure
primary_script: 04-computation/factorial_fixed_prime_exponent_orbit_compiler_thm4257.py
primary_output: 05-knowledge/results/factorial_fixed_prime_exponent_orbit_compiler_thm4257.out
independent_audit_script: 04-computation/factorial_fixed_prime_exponent_orbit_independent_audit_thm4257.py
independent_audit_output: 05-knowledge/results/factorial_fixed_prime_exponent_orbit_independent_audit_thm4257.out
primary_script_sha256: 8d5df3cec17e38f9ef55fda3e133f5ceba6701bb6b239ed589b16ee8f5fe91ed
primary_output_sha256: 1979327118b98bd954def0e6cdd7f29a6bacf1b545fcc3fad1744e4bcabd3ac5
independent_audit_script_sha256: b9bf812914067c530116f779e4872bd93d5e8f1cd95e04be5bb2f79faab6cc53
independent_audit_output_sha256: af860179394b1ff9efa6ac3fcb7533854ed5acf7459dd09a6db44e8af74ccd00
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary and independent programs use distinct residue
  representations, contain no optimization-removable truth gates, and
  reproduce byte-identical normal/-O transcripts. The analytic fixed-orbit
  density proof identifies the orbit as one or two complete high-bit
  cylinders before applying THM-4250's disjoint reset blocks.
---

# THM-4257 -- fixed-prime exponent-orbit density-one compiler

**PROVED RELATIVE TO THM-4244/4250 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Objects, inherited map, and loss ledger

Fix an odd integer `b>=3` and `ell>=1`, and put

```text
M=2^ell,       [x]_ell=the least residue of x modulo M.           (1)
```

The proved THM-4250 suffix language is

```text
R_b(ell)={odd r mod M:
  [sr]_ell & [(b-s)r]_ell !=0 for every 1<=s<=(b-1)/2}.           (2)
```

Let `P` be any positive odd integer, write

```text
m_ell(P)=ord_M(P),                                                (3)
```

and define the exponent-orbit certificate set

```text
K_(b,ell)(P)
 ={k mod m_ell(P): [P^k]_ell in R_b(ell)}.                        (4)
```

Only the residue of `P` modulo `M` enters (3)--(4). In the factorial
application `P` is a fixed prime; in the finite atlas the symbol `u` denotes
an arbitrary odd unit residue and avoids pretending that every small residue
representative is itself prime.

The map in this packet is

```text
exponent class k  |-->  suffix [P^k]_ell.                         (5)
```

It preserves exactly the suffix predicate (2). It destroys every bit above
`ell-1`, so it does **not** classify full integer overlaps, actual common
factors, or bad factorial windows. Those losses are used explicitly in the
hostile controls below.

## 2. Exact orbit iff and infinite exponent progressions

### Theorem 1 (exact suffix-orbit compiler)

For every odd `b>=3`, every `ell>=1`, every positive odd integer `P`, and
every integer `k>=0`,

```text
[P^k]_ell in R_b(ell)
 iff
k mod m_ell(P) in K_(b,ell)(P).                                  (6)
```

Thus `K_(b,ell)(P)` is an exact iff classifier for the finite suffix
certificate, not merely a sufficient sublist.

Moreover, `0` never belongs to `K_(b,ell)(P)`. Hence every class represented
by `1<=kappa<m_ell(P)` in a nonempty `K` supplies the infinite positive
arithmetic progression

```text
k=kappa+t*m_ell(P),               t=0,1,2,... .                   (7)
```

#### Proof

By definition of multiplicative order,

```text
P^(k+m_ell(P))=P^k mod M.                                        (8)
```

Therefore (2) is constant on exponent classes modulo `m_ell(P)`, and (6) is
exactly (4) with well-definedness made explicit.

For the final assertion, `P^0=1`. In the `s=1` pair, `[1]_ell` is odd while
`[b-1]_ell` is even, so their bit-zero intersection is zero; if the latter
residue is zero, the whole intersection is still zero. Thus

```text
1 notin R_b(ell),                                                 (9)
```

and the zero exponent class is absent. Each remaining modular class has the
positive representatives (7). QED.

### Corollary 2 (fixed-prime factorial progressions)

Fix `q>=1`, put

```text
a=2^q b,                                                         (10)
```

and let `P>a` be prime. If the positive integer `k` has
`k mod m_ell(P)` in `K_(b,ell)(P)`, then the low-bit overlap in (2) is also
an overlap of the full products

```text
(sP^k)&((b-s)P^k),       1<=s<=(b-1)/2.                         (11)
```

Use THM-4244's notation

```text
L(x^r)=r!,       A_n^d(v)=L((d-x+v*x^2)^n).                     (11a)
```

That theorem therefore gives, with `N=aP^k`,

```text
gcd_Q(A_(N-1)^(N+1),A_N^(N+1))=1,                               (12)
```

and for every exact quadratic

```text
f=alpha+beta*x+gamma*x^2,       alpha*beta*gamma!=0,             (13)
```

at least one of

```text
L(f^(N-1)),             L(f^N),             L(f^(N+1))           (14)
```

is nonzero. Every class in `K` yields infinitely many such conclusions by
(7).

#### Proof

If two residues modulo `2^ell` share a one-bit, their full nonnegative integer
representatives share that same low bit. Hence (2) implies every overlap in
(11). The hypotheses `P>a>=6` put the row inside THM-4244's prime-power
factorial scope. Equations (12)--(14) are exactly its proved one-way
empty-barcode consequences. QED.

The converse is deliberately absent: failure of (2) can be repaired by a
higher bit, and even a nonempty THM-4244 candidate barcode need not contain an
actual factor.

### Theorem 2A (all-height union and AP completion)

For a fixed positive odd `P`, define

```text
E_ell(P)={k>=1:k mod m_ell(P) in K_(b,ell)(P)},

G_b(P)={k>=1:(sP^k)&((b-s)P^k)!=0
                 for every 1<=s<=(b-1)/2}.                       (14a)
```

Then

```text
E_ell(P) subset E_(ell+1)(P),
G_b(P)=union_(ell>=1) E_ell(P).                                  (14b)
```

Consequently every one-off full-height good exponent lies in an infinite
good arithmetic progression.

#### Proof

A common bit below level `ell` survives at every higher suffix level, giving
the first inclusion. Membership in any `E_ell` supplies full integer overlap
witnesses, so the union is contained in `G_b(P)`. Conversely, if `k` lies in
`G_b(P)`, choose one common occupied bit for each of the finitely many
complementary pairs and take `ell` above every chosen position. Their suffix
reductions retain all witnesses, so `k in E_ell(P)`. Theorem 1 then places
every exponent congruent to `k mod m_ell(P)` in the same good progression.
QED.

## 3. Self-contained binary unit coordinates and order formulas

The finite orbit can be constructed without assuming a black-box unit-group
classification.

### Lemma 3 (`C_2 x C_(2^(ell-2))` from first principles)

For `ell>=3`, every odd unit modulo `2^ell` has a unique expression

```text
u=(-1)^epsilon 5^t mod 2^ell,
epsilon in {0,1},             t mod 2^(ell-2).                    (15)
```

Consequently

```text
(Z/2^ell Z)^times = <-1> x <5>
                    ~= C_2 x C_(2^(ell-2)).                      (16)
```

#### Proof

For every `j>=0`,

```text
nu_2(5^(2^j)-1)=j+2.                                            (17)
```

At `j=0` this is `nu_2(4)=2`. If `X=5^(2^j)`, then `X=1 mod 4`, so
`nu_2(X+1)=1`; factoring `X^2-1` raises the valuation by exactly one. This
proves (17) inductively. Thus `5` has exact order `2^(ell-2)` modulo `2^ell`.

All its powers are `1 mod 4`. There are exactly `2^(ell-2)` odd residue
classes that are `1 mod 4`, so `<5>` is precisely that set. Every unit that is
`3 mod 4` becomes `1 mod 4` after multiplication by `-1`, proving existence
in (15). The two signs are disjoint modulo four, and (17) gives uniqueness of
`t`. QED.

For `ell=1` the unit group is trivial; for `ell=2` it is `C_2`. These small
levels are handled directly in the atlas.

### Lemma 4 (exact order of a fixed odd integer)

Let `P` be a positive odd integer. If `P=1` then `m_ell(P)=1` for every
`ell`; below assume `P>1`. For `ell=1`, `m_ell(P)=1`. For `ell>=2`:

```text
P=1 mod 4, alpha=nu_2(P-1):
  m_ell(P)=2^max(0,ell-alpha);                                  (18)

P=3 mod 4, beta=nu_2(P^2-1):
  m_ell(P)=2^max(1,ell-beta+1).                                (19)
```

#### Proof

When `P=1 mod 4`, the same factorization used in (17) gives

```text
nu_2(P^(2^j)-1)=alpha+j,       j>=0.                            (20)
```

When `P=3 mod 4`, start after one squaring; `P^2=1 mod 8`, and obtain

```text
nu_2(P^(2^j)-1)=beta+j-1,      j>=1.                            (21)
```

If an exponent is `2^j` times an odd number, factoring `X^t-1` with odd `t`
shows that the odd factor does not change the two-adic valuation. Therefore
the first exponent giving divisibility by `2^ell` is exactly the power of two
displayed in (18) or (19). QED.

### Corollary 5 (the `-1` orbit obstruction)

For `ell>=3`,

```text
-1 in <P mod 2^ell>       iff       P=-1 mod 2^ell.              (22)
```

#### Proof

Write `P=(-1)^epsilon 5^u` using (15). If `P^k=-1`, the sign equation forces
`epsilon=1` and `k` odd, while the cyclic equation is

```text
u*k=0 mod 2^(ell-2).                                            (23)
```

An odd `k` is invertible modulo `2^(ell-2)`, so `u=0` and `P=-1`. The reverse
direction is immediate. QED.

Thus the universal `-1` certificate meets `<P>` only in the single congruence
class stated in (22). This blocks any argument from that one residue alone;
the complete-cylinder argument below is what lawfully transfers density one
to a fixed orbit.

## 4. Exact exponent-class lift law

Fix the same odd integer `P` at two consecutive levels and abbreviate

```text
M=2^ell,       m=m_ell(P),       m'=m_(ell+1)(P),
d=m'/m.                                                            (24)
```

Reduction of an order modulo `2M` shows `m|m'`. Conversely, write
`P^m=1+cM`; then

```text
(P^m)^2=1 mod 2M.                                                (25)
```

Hence

```text
d in {1,2}.                                                       (26)
```

For `k mod m`, its reached exponent sheets modulo `m'` are

```text
E_k={k+j*m:0<=j<d}.                                              (27)
```

Let `r=[P^k]_ell`. Every `e in E_k` reaches one of the two residue children

```text
[P^e]_(ell+1)=r+epsilon_e M,       epsilon_e in {0,1}.            (28)
```

If `d=2`, the two sheets reach the two distinct children. If `d=1`, the orbit
reaches exactly one child.

For `1<=u<b`, write

```text
ur=x_u+M q_u,        x_u=[ur]_ell,        h_u=q_u mod 2.          (29)
```

Let `F(r)` be the set of complementary pairs failing below level `ell`.
For each failed pair let `e_s` and `o_s` denote its even and odd multiplier.
Define its allowed child-bit set by

```text
Epsilon(r)={0,1},                         if F(r)=empty;

Epsilon(r)={1 xor h_(o_s)},               if F(r)!=empty,
   every h_(e_s)=1, and every h_(o_s) agrees;

Epsilon(r)=empty,                         otherwise.             (30)
```

### Theorem 6 (one/two-sheet exponent lift)

For every `k mod m` and `e in E_k`,

```text
e in K_(b,ell+1)(P)       iff       epsilon_e in Epsilon(r).     (31)
```

In particular, a closed parent contributes all `d` reached sheets. A failed
parent contributes no sheet or one sheet; when `d=2`, it contributes one
exactly under the compatibility condition in the middle line of (30).

#### Proof

For a residue child `r+epsilon M`, direct multiplication gives

```text
[u(r+epsilon M)]_(ell+1)
 =x_u+M*(h_u xor (epsilon*(u mod 2))).                            (32)
```

An existing low-bit overlap persists in both children. For a failed pair the
only possible new overlap is the new bit: its even multiplier has fixed bit
`h_(e_s)`, while its odd multiplier has bit `h_(o_s) xor epsilon`. One global
child closes all failed pairs precisely under (30). Equation (28) restricts
this exact residue-child test to the sheets actually reached by the fixed
prime orbit, proving (31). QED.

If `a_k` is the number of allowed reached sheets over a failed parent class,
then the disjoint lift partition gives the exact recurrence

```text
|K_(b,ell+1)(P)|
 =d*|K_(b,ell)(P)| + sum_(k notin K_(b,ell)(P)) a_k.              (33)
```

Dividing by `m'=dm` proves a new monotonicity statement:

```text
|K_(b,ell+1)(P)|/m_(ell+1)(P)
 >= |K_(b,ell)(P)|/m_ell(P).                                    (34)
```

Thus fixed-orbit certificate density is monotone with suffix height.

### Theorem 7 (quantitative fixed-orbit density one)

Assume `P>1`, and put

```text
L=ceil(log_2(b-1))+1,

c_P=nu_2(P-1),       if P=1 mod 4,
c_P=nu_2(P^2-1),     if P=3 mod 4.                              (34a)
```

For every `ell>=c_P`,

```text
1-|K_(b,ell)(P)|/m_ell(P)
 <=(1-2^(-L))^floor((ell-c_P)/L).                              (34b)
```

Hence the left-hand orbit density tends to one exponentially. The full-height
good exponent set `G_b(P)` from (14a) has natural density one. In particular,
every fixed odd `P>1` has a suffix-good infinite exponent progression.

#### Proof

First identify the orbit exactly. If `P=1 mod 4`, put
`alpha=nu_2(P-1)=c_P`. The order formula (18) and containment give, for every
`ell>=alpha`,

```text
<P mod 2^ell>={x:x=1 mod 2^alpha}.                              (34c)
```

Indeed the right side contains the orbit and has exactly
`2^(ell-alpha)=m_ell(P)` elements.

If `P=3 mod 4`, put `beta=nu_2(P^2-1)=c_P`. Since
`nu_2(P-1)=1`, one has

```text
P=2^(beta-1)-1 mod 2^beta.                                     (34d)
```

Even and odd powers lie respectively in the two cylinders

```text
<P mod 2^ell>
 ={x:x=1 mod 2^beta}
  union {x:x=2^(beta-1)-1 mod 2^beta},       ell>=beta.          (34e)
```

Their total size is `2^(ell-beta+1)`, equal to (19), so containment is
equality.

Thus the bits in positions `c_P,...,ell-1` range freely and uniformly over
one complete cylinder in (34c), or over each of two complete cylinders in
(34e). Partition those `ell-c_P` free positions into
`floor((ell-c_P)/L)` disjoint blocks of length `L`. THM-4250 proves that any
one all-one block forces suffix closure. Therefore a nonclosing orbit word
must avoid the all-one pattern in every chosen block. The blocks are
independent and each is all one in a fraction `2^(-L)` of the cylinder,
giving (34b).

The periodic set `E_ell(P)` has natural density exactly
`|K_(b,ell)(P)|/m_ell(P)` and is contained in `G_b(P)` by Theorem 2A. Hence
the upper density of the complement of `G_b(P)` is bounded by the right side
of (34b) for every `ell`; letting `ell` grow proves natural density one.
QED.

There is also an explicit guaranteed AP level. Put

```text
ell_*=c_P+L.                                                     (34f)
```

The negative collar already meets the orbit there. If `P=1 mod 4`, use

```text
w=2^c_P-1.                                                       (34g)
```

If `P=3 mod 4`, the odd- and even-power cylinders admit respectively

```text
w=2^(c_P-1)+1,             w=2^c_P-1.                           (34h)
```

Each target `2^ell_*-w` lies in its displayed orbit cylinder, and
`w(b-1)<=2^(ell_*-1)` follows from `b-1<=2^(L-1)`. Thus THM-4250's negative
collar proves `K_(b,ell_*)(P)` nonempty without search. If in addition
`P>2^q b` is prime, Corollary 2 transfers a density-one set of exponents—and
in particular this explicit AP—to the THM-4244 factorial conclusions.

Density one is not cofiniteness: no claim is made that only finitely many
exponents fail.

## 5. Structural tests beyond enumeration

### 5.1 Negative-collar orbit criterion

THM-4250 proves the explicit sufficient collar

```text
C_b(ell)={-w mod M: w>=1 odd and w(b-1)<=M/2}
          subset R_b(ell).                                      (35)
```

For `ell>=3`, write the unique coordinates

```text
P=(-1)^epsilon 5^u,
c=(-1)^eta 5^v,                  c in C_b(ell),
n=2^(ell-2).                                                     (36)
```

Then

```text
c in <P>
 iff there exists k satisfying
     epsilon*k=eta mod 2,       u*k=v mod n.                    (37)
```

This is a compact necessary-and-sufficient collar-hit test. More explicitly:

- if `u=0`, one needs `v=0` and either `epsilon=1` or `eta=0`;
- if `u!=0`, put `g=gcd(u,n)`. One needs `g|v`; if `epsilon=0` one also needs
  `eta=0`, while if `epsilon=1` one needs `eta=v/g mod 2`.

Indeed, in the second case `u/g` is odd and invertible modulo the even power
of two `n/g`, so every solution of the cyclic congruence has parity `v/g`.
Every solution exponent obtained this way lies in `K`.

For a nontrivial example, at `(b,ell,P)=(13,7,43)` one has

```text
C_13(7)={127,125,123},
43=(-1)*5^29,       123=-5=(-1)*5^1 mod 128.                    (38)
```

The system is `k` odd and `29k=1 mod 32`, whose unique solution is

```text
k=21 mod 32.                                                     (39)
```

Thus fixed prime `P=43` gives an infinite exponent family for multiplier
`a=26`, even though its orbit is empty at the shorter level `ell=5`.

### 5.2 Subgroup/coset containment is exactly exponent-AP containment

Let `d|m_ell(P)` and fix `kappa mod d`. The image of the exponent congruence
class is exactly

```text
{P^k:k=kappa mod d}
 =P^kappa <P^d> mod M.                                          (40)
```

Therefore

```text
every k=kappa mod d lies in K_(b,ell)(P)
 iff
P^kappa <P^d> subset R_b(ell).                                 (41)
```

This gives a proof-level minimization procedure: inspect only divisors of the
power-of-two order and retain the inclusion-maximal accepted cosets. It is
strictly stronger than printing individual exponent classes.

## 6. Exact small-core atlas

Write `lambda(b)` for the least level with `R_b(ell)` nonempty, and
`ell_0(b)` for the least level satisfying `2^(ell-1)>=b`. Exact enumeration
gives

```text
b             3  5  7  9 11 13
lambda(b)     2  3  3  4  4  4
ell_0(b)      3  4  4  5  5  5.                                (42)
```

The controls are

```text
R_3(2)={3},          R_5(3)={5,7},          R_7(3)={7},
R_7(4)={3,5,7,15}.                                             (43)
```

The last row reproduces THM-4250's modulus-16 `b=7` exponent table. At the
least nonempty target level,

```text
R_9(4) ={9,15},
R_11(4)={3,15},
R_13(4)={5,15}.                                                 (44)
```

The non-`-1` seed masks, listed in complementary-pair order, are

```text
(b,r)=(9,9):    (8,2,2,4),
(b,r)=(11,3):   (2,2,8,4,2),
(b,r)=(13,5):   (4,2,2,4,8,2).                                 (45)
```

All entries are nonzero, so (44)'s new seeds have direct hand certificates.
At the least universal level the exact residue sets are

```text
R_9(5) ={9,15,19,25,29,31},
R_11(5)={3,7,9,15,19,21,27,31},
R_13(5)={5,15,21,31}.                                          (46)
```

Here is the complete nonempty fixed-unit orbit atlas modulo 32. Omitted unit
classes have empty `K`; `m` is their exact order.

```text
b=9:
u       3       5       9       11      13      15      19
m       8       8       4       8       8       2       8
K     2,5,6   2,3,6    1,3    2,3,6   2,5,6     1     1,2,6

u       21      25      27      29      31
m        8       4       8       8       2
K      2,6,7    1,3    2,6,7   1,2,6     1
empty u: 1,7,17,23.

b=11:
u       3        5      7  9       11       13     15
m       8        8      4  4        8        8      2
K     1,2,3,5   5,6     1  1     3,5,6,7   2,3     1

u       19       21     23 25      27       29     31
m        8        8      4  4       8        8      2
K      1,2,5,7   1,6     3  3    1,3,6,7   2,7     1
empty u: 1,17.

b=13:
u       5       13      15      21      29      31
m       8        8       2       8       8       2
K      1,5      3,7      1      1,5     3,7      1
empty u: 1,3,7,9,11,17,19,23,25,27.                            (47)
```

The scripts print every unit row, including empty-orbit hostiles, exact
membership period, and inclusion-maximal exponent cosets.

## 7. New fixed-prime infinite families

Take `q=1`, so the new multipliers are `a=18,22,26`. Direct orbit evaluation
inside (46), followed by the coset minimization (41), gives

```text
a=18, P=19, ell=5, m=8:
  K={1,2,6} = {k=1 mod 8} union {k=2 mod 4};

a=22, P=23, ell=5, m=4:
  K={3} = {k=3 mod 4};

a=26, P=29, ell=5, m=8:
  K={3,7} = {k=3 mod 4};

a=26, P=43, ell=7, m=32:
  {k=21 mod 32} subset K by the negative collar.                 (48)
```

The first three primes are the first primes above their respective
multipliers. The exponent moduli in (48) are minimal for the displayed
cosets: at the preceding divisor a rejected exponent enters the class.

For transparency, the relevant power cycles modulo 32 are

```text
<19>: 1,19,9,11,17,3,25,27;
<23>: 1,23,17,7;
<29>: 1,29,9,5,17,13,25,21.                                    (49)
```

Intersecting (49) with (46) gives the three exact `K` rows in (48). The
accepted suffixes are not the prior universal `-1` suffix: they include
`19,9,25`, `7`, `5,21`, and `-5`. Nor are these `b=7` lanes from THM-4250.
Thus (48) supplies genuinely new fixed-prime infinite exponent families.

For every exponent in (48), Corollary 2 proves the rational coprimality and
three-moment nonvanishing conclusions (12)--(14).

## 8. Empty-orbit and direction hostiles

1. **An exact empty suffix orbit can miss a true full-height closure.** At
   `(b,ell,P)=(9,4,23)`,

   ```text
   <23 mod 16>={1,7},             R_9(4)={9,15},
   K_(9,4)(23)=empty.                                             (50)
   ```

   At `k=1`, the full overlap masks are `(16,32,0,80)`, so one pair really
   fails. But at `k=3`, with `H=23^3`, the full integer masks are

   ```text
   (11264,19456,3072,44032),                                    (51)
   ```

   all nonzero. Thus THM-4244 closes the `a=18,P=23,k=3` factorial row even
   though the level-four exponent orbit is empty. Higher bits are essential.

2. **The same fixed orbit can enter later.** The `P=43,b=13` orbit is empty
   at `ell=5` but hits `-5` at `ell=7`, as (38)--(39) show. An empty finite
   observer is not a permanent orbit obstruction.

3. **One target is not the full density mechanism.** Equation (22) says
   exactly that the only cyclic orbit containing `-1` is the orbit generated
   by `-1`. The finite atlases still contain 36 empty unit-orbit rows at their
   declared levels. Theorem 7 transfers density one only after proving that a
   fixed orbit becomes one or two complete high-bit cylinders.

4. **Suffix iff is not factorial iff.** Equation (6) is an equivalence only
   for membership in the declared low-bit language. Nonmembership does not
   imply a nonempty full barcode, a factor, a common root, or three vanishing
   moments.

5. **Prime and multiplier inequalities remain load-bearing.** The pure orbit
   statements allow every odd unit. The factorial transfer requires an actual
   prime `P>2^q b`, `q>=1`, and the exact quadratic support inherited from
   THM-4244.

## 9. Verification universe and independence contract

The primary companion uses direct modular products, explicit cyclic power
orbits, a signed power-of-five chart, and direct exponent-sheet lifting. It
checks:

```text
all 2,044 odd unit cells through ell=11 for coordinate uniqueness,
the valuation order formula, and the exact -1 orbit criterion;
all odd b=3,...,31 through ell=9 for 437,100 orbit-iff cells;
218,790 exponent-parent lift cells through ell=8;
every negative-collar residue and 7,620 collar-orbit rows;
866,490 subgroup-coset / exponent-AP equivalences;
the complete control and target atlases in (42)--(47);
the four fixed-prime families in (48);
and the high-bit-repair hostile (50)--(51).                       (52)
```

The independent companion shares no primary code and uses occupied-bit sets,
cycle dictionaries, direct child enumeration, and a separately derived
signed-power coordinate solver. Its exact universe and both frozen hashes are
recorded in the packet manifest.

Both scripts use explicit runtime gates rather than optimization-removable
assertions. Normal and `-O` executions must byte-match their frozen outputs.

## 10. Scope firewalls

- The mathematical novelty is an exact finite-orbit compiler, its lift law,
  its coset/collar structure, and the new fixed-prime families. The factorial
  implication is exactly the proved THM-4244 transfer.
- No converse to THM-4244 or THM-3474 is asserted.
- No claim concerns `SFC(3)` or the three-variable Factorial Conjecture; this
  remains the exact three-slot quadratic slice isolated by THM-4244
  (MISTAKE-350 firewall).
- No factorial polynomial is reduced modulo a small prime, so the
  MISTAKE-363 denominator/leading-degree hazard does not arise.
- The finite unit atlas classifies residue orbits, not the distribution of
  primes in residue classes. The four displayed primes are explicit and need
  no prime-distribution theorem.

**QED.**
