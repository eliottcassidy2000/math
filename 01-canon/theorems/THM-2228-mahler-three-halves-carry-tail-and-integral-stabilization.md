---
id: THM-2228
title: "Mahler three-halves carry tail and integral stabilization"
status: >
  PROVED + VERIFIED-EXACT. Binary carry words are homeomorphic to the
  2-adic integers and conjugate the shift to the 2-adic extension of
  a -> ceil(3a/2). A word codes a positive Mahler Z-number exactly when
  every weighted suffix tail is strictly below one and its compatible
  canonical residues modulo 2^m eventually stabilize to an ordinary
  positive integer. The two conditions are independent, the equality-one
  tail boundary defeats every finite strict-prefix test, no positive-integer
  itinerary is eventually periodic, and explicit denominator-19 rational
  pseudo-Z-numbers survive every prescribed finite horizon. This is an exact
  reformulation and obstruction theorem, not a proof that Z-numbers exist
  or do not exist.
source: codex-2026-07-24-mahler-carry-tail
depends_on: []
related:
  - THM-2049-the-DC2-Ore-boundary-correction-complex-is-acyclic
  - THM-2163-radix-relation-carry-descent
  - THM-2172-frobenius-collapse-of-mahler-and-twisted-cyclic-packets
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
script: 04-computation/mahler_three_halves_carry_stabilization_thm2228.py
output: 05-knowledge/results/mahler_three_halves_carry_stabilization_thm2228.out
script_sha256: 5070dc34ebfa21ac609866014fac7f4c086eea5856d0788d6cf44dbd350df14a
output_sha256: 85ddadd5904935ca5ceee788ae529bba72b34a59761207887cb932bde6b7a39e
hash_basis: working-tree bytes (LF)
---

# THM-2228 -- Mahler carries need both a tail and termination

A positive real number `alpha` is a Mahler `Z`-number when

```text
0<=fractional_part(alpha(3/2)^n)<1/2       for every n>=0.       (1)
```

This theorem separates the two genuinely different requirements hidden in
(1): a real suffix-tail inequality and ordinary termination of a compatible
2-adic carry state.

## 1. Every binary word is a compatible 2-adic state

For a word

```text
c=(c_0,c_1,...) in {0,1}^N
```

define

```text
Y_n(c)=sum_(j>=0)c_(n+j)(2/3)^(j+1),                              (2)

C_m(c)=sum_(0<=j<m)c_j 2^j 3^(m-1-j),                            (3)

r_m(c)=-3^(-m)C_m(c) mod 2^m,       0<=r_m(c)<2^m.               (4)
```

The inverse in (4) exists because `3` is odd. The identity

```text
C_(m+1)=3C_m+c_m2^m                                      (5)
```

gives

```text
r_(m+1)=r_m mod 2^m.                                     (6)
```

Hence the residues define one point

```text
Phi(c) in Z_2.                                           (7)
```

In fact, `Phi` is a homeomorphism. For every finite word of length `m`,
(4) is the unique residue class `A mod 2^m` for which

```text
2^m a_m=3^m A+C_m                                       (8)
```

is integral. Thus the `2^m` binary cylinders map bijectively to the `2^m`
residue classes modulo `2^m`; passage to the inverse limit proves the claim.

Define the continuous map on `Z_2`

```text
T(a)=(3a+(a mod 2))/2.                                  (9)
```

On nonnegative integers this is

```text
T(a)=ceil(3a/2).                                        (10)
```

The cylinder bijection also gives the exact conjugacy

```text
Phi(shift(c))=T(Phi(c)).                                (11)
```

Thus a binary carry word always describes a legitimate **formal** 2-adic
orbit. It need not describe an orbit starting at an ordinary integer.

## 2. Exact characterization of Mahler Z-numbers

> **Carry-tail characterization.** A word `c in {0,1}^N` codes a positive
> Mahler `Z`-number if and only if
>
> ```text
> Y_n(c)<1                         for every n>=0,       (12)
> ```
>
> and the canonical integers `r_m(c)` eventually stabilize: there are
> `A>=1` and `m_0` such that
>
> ```text
> r_m(c)=A                         for every m>=m_0.     (13)
> ```
>
> When these conditions hold, the coded number is unique and equals
>
> ```text
> alpha=A+Y_0(c)/2.                                      (14)
> ```

### Forward direction

Suppose `alpha` satisfies (1), and write

```text
alpha(3/2)^n=a_n+x_n,
a_n in Z_(>=0),               0<=x_n<1/2.             (15)
```

Put

```text
c_n=2a_(n+1)-3a_n=3x_n-2x_(n+1).                      (16)
```

The right side lies strictly between `-1` and `3/2`; since it is an
integer,

```text
c_n in {0,1},                 c_n=a_n mod 2,           (17)
```

and therefore `a_(n+1)=T(a_n)`. Let `Z_n=2x_n`. Equation (16) becomes

```text
Z_(n+1)=(3/2)Z_n-c_n.                                  (18)
```

Iterating backward for `M` steps gives

```text
Z_n=sum_(0<=j<M)c_(n+j)(2/3)^(j+1)
       +(2/3)^M Z_(n+M).                               (19)
```

Because `0<=Z_n<1`, the remainder tends to zero. Hence `Z_n=Y_n(c)`,
which proves (12), and (14) holds with `A=a_0`. Iterating the integer
recurrence gives (8), so `A=r_m mod 2^m` for every `m`. Once `2^m>A`, the
canonical representative is literally `r_m=A`, proving (13). The case
`A=0` would force every `a_n` and `c_n` to vanish and then `alpha=0`, so
positivity gives `A>=1`.

### Reverse direction

Now assume (12)--(13). Compatibility (6) implies

```text
A=r_m mod 2^m                    for every m.           (20)
```

Therefore

```text
a_m=(3^m A+C_m)/2^m                                  (21)
```

is a positive integer. Equations (5) and (21) give

```text
2a_(m+1)=3a_m+c_m,          c_m=a_m mod 2,             (22)
```

so this is exactly the integer orbit of `A` under `T`. Equations (2) and
(18) imply inductively that the number (14) satisfies

```text
alpha(3/2)^n=a_n+Y_n(c)/2.                             (23)
```

Condition (12) makes the second term its fractional part and places it in
`[0,1/2)`. Thus `alpha` satisfies (1). Uniqueness follows from (13)--(14).

This proves the characterization.

### A finite certificate for each strict tail

For a fixed suffix and prefix length, put

```text
C_(n,m)=sum_(0<=j<m)c_(n+j)2^j3^(m-1-j),
Y_(n,m)=2C_(n,m)/3^m.                                  (24)
```

The omitted binary digits have the sharp uniform envelope

```text
Y_(n,m)<=Y_n<=Y_(n,m)+2(2/3)^m.                        (25)
```

Consequently the infinite strict inequality has the exact semidecision
certificate

```text
Y_n<1
 iff
there is an m>=1 with 2C_(n,m)+2^(m+1)<3^m.            (26)
```

The forward implication follows by taking `m` so large that the envelope
width is smaller than `1-Y_n`; the reverse implication is immediate from
(25). By contrast, the raw inequalities

```text
2C_(n,m)<3^m for every m                               (27)
```

imply only `Y_n<=1`. Equality at a finite stage is impossible because the
left side is even and the right side is odd, so this is a genuine infinite
boundary rather than a missed finite equality case.

## 3. Neither condition can be discarded

### Safe real tails without ordinary stabilization

For the pure period

```text
c=(100)^(infinity),
```

the three suffix phases are

```text
Y_n in {18/19,8/19,12/19}.                             (28)
```

They all satisfy (12). But one period has `C_3=9`, so its 2-adic state is

```text
Phi(c)=-9/(3^3-2^3)=-9/19.                             (29)
```

The canonical residues of this negative rational 2-adic integer never
stabilize to a nonnegative ordinary integer.

### Ordinary stabilization without safe real tails

Start the integer recurrence at `A=1`. Its carry word begins

```text
1011,
```

whose contribution to `Y_0` is already

```text
2/3+(2/3)^3+(2/3)^4=94/81>1.                          (30)
```

Its canonical residues stabilize to `1`, but (12) fails.

Together (28)--(30) show that the symbolic word or the real tail alone is
an invalid quotient.

## 4. The strict boundary is invisible to finite prefix tests

Set `z_0=1` and recursively take the greedy digits

```text
c_n=floor(3z_n/2),
z_(n+1)=3z_n/2-c_n.                                   (31)
```

Then `c_n in {0,1}`, and `0<z_n<1` for every `n>=1`. Indeed, (31) keeps
the orbit in `[0,1)`, while writing `z_n` over denominator `2^n` shows its
numerator remains odd and hence never vanishes.

For every `n,M`,

```text
z_n=sum_(0<=j<M)c_(n+j)(2/3)^(j+1)
       +(2/3)^M z_(n+M).                               (32)
```

Every finite suffix-prefix sum in (32) is therefore strictly below one.
Nevertheless,

```text
Y_0(c)=z_0=1.                                         (33)
```

Thus replacing the infinite strict condition (12) by all finite strict
prefix inequalities is false. The equality boundary carries information
that no fixed cylinder sees.

## 5. Periodic words and arbitrarily long pseudo-orbits

For a pure period `w` of length `L`, (8) over one period gives

```text
Phi(w^(infinity))=-C_L(w)/(3^L-2^L).                  (34)
```

If `w` is nonzero, this rational number is negative; if `w` is zero, it is
zero. Consequently:

> No carry itinerary of a positive ordinary integer under `T` is eventually
> periodic.

Indeed, an eventually periodic tail would make the positive integer
`T^N(A)` equal in `Z_2` to the nonpositive rational in (34), which is
impossible.

There is also an exact automata consequence. Let

```text
I_+={c(A): A is a positive ordinary integer}
```

be the language of infinite parity itineraries. Then `I_+` is not
omega-regular: it is nonempty but contains no ultimately periodic word,
whereas any nonempty language accepted by a finite-state Buechi automaton
contains an ultimately periodic word (repeat a loop between two visits to
one accepting state). This statement concerns ordinary-integer realization;
it does not assert that the unknown `Z`-number language or its finite-prefix
language is nonsofic.

There is nevertheless no finite-time obstruction to the Mahler inequality.
For every `M>=1`, choose the unique

```text
1<=k_M<=18,                 k_M 2^M=9 mod 19,          (35)
```

and put

```text
alpha_M=k_M 2^M/19.                                    (36)
```

Since `3/2=11 mod 19`, for every `0<=n<=M`,

```text
fractional_part(alpha_M(3/2)^n)
 =[9*11^n]_19/19
 in {9/19,4/19,6/19}
 subset [0,1/2).                                      (37)
```

Hence positive rational numbers survive the `Z`-number test for arbitrarily
long prescribed horizons. Their integer parts are

```text
A_M=(k_M2^M-9)/19
   =-9/19 mod 2^M,             0<=A_M<2^M,            (38)
```

so they are exactly the growing canonical prefixes of the safe formal state
`-9/19`. These finite points follow the same safe phase cycle for arbitrarily
long horizons and converge 2-adically to one formal state, yet that state
never terminates to an ordinary positive initial integer.

The denominator `19` is the first odd `q>1`, coprime to `6`, whose
multiplication-by-`3/2` permutation has a nonzero cycle entirely in the
lower half; the unique such cycle at `q=19` is `{4,6,9}`. This last
minimality statement is a finite exact check, not a claimed general
classification of rational cycles.

## 6. Scope and exact referee

This theorem does not decide whether a word satisfying both (12) and (13)
exists. It proves that:

1. finite carry cylinders are exactly dyadic residue classes;
2. the formal inverse limit is all of `Z_2`, much larger than the ordinary
   nonnegative integers required by Mahler's problem;
3. real tail admissibility and ordinary termination are independent; and
4. no universal fixed orbit-time cutoff can refute all candidates.

Run

```bash
python3 04-computation/mahler_three_halves_carry_stabilization_thm2228.py
python3 -O 04-computation/mahler_three_halves_carry_stabilization_thm2228.py
```

The exact companion verifies the cylinder bijection and conjugacy through
depth sixteen (`131070` nonempty words), the two hostile words, `256` greedy
boundary digits, all `8190` binary period blocks through length twelve, and
the open-tail certificate on `3586` periodic suffix phases through length
eight (`1793` safe phases, all certified). It also checks the odd-denominator
scan through `19` and `8384` rational pseudo-orbit points through horizon
`128`. It uses exact integers and `Fraction` arithmetic, explicit raising
checks, and produces the same frozen output normally and under optimization.

QED.
