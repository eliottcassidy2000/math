---
id: THM-3423
title: "Odd-interval ratio complement and dyadic clique law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A symmetric odd interval
  O_L in F_p has, beyond an explicit
  cubic threshold, an exact ratio-set complement consisting of the reduced
  opposite-parity fractions of height at most h.  The Cayley graph of
  pairwise-disjoint multiplicative dilates of O_L then has exact clique number
  1+floor(log_2(h-1)); its sharp coloring is the 2-adic valuation modulo that
  number.  The thresholds are sufficient rather than asserted minimal.
source: root-2608-odd-interval-ratio-law-2026-08-15
audit: independent proof reconstruction; multi-prime finite-field controls for every 2<=h<=13; normal/optimized/stored-output replay; hash and documentation audit clean
related:
  - THM-3420-prime-rank-seven-zero-and-half-twist-splitter-closures
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3416-zero-mode-cochain-global-rank-six-support
script: 04-computation/odd_interval_ratio_complement_thm3423.py
output: 05-knowledge/results/odd_interval_ratio_complement_thm3423.out
script_sha256: 230b9c90f341dc272b50392c6537c3fe0363ccf8e5f055a52c33c20df9b61a3f
output_sha256: 9227e53fe70e01793ff56f0f074f48d8da439f0f7138b539d462cd455a72aade
semantic_sha256: 05f52e47dd904880524b09744efd98c90e6e53b4b37976dbe805ff0d8199a6c2
hash_basis: LF-normalized bytes
---

# THM-3423 -- odd-interval ratio complement and dyadic clique law

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Exact statement

Fix an integer `h>=2`, an odd prime `p`, and a positive odd integer `L`.  Put

```text
c=p-hL
```

and assume

```text
0<c<2h,
p >= (h+1)(h(h+1)+c).                                    (1)
```

Thus `L` is the largest positive odd integer satisfying `hL<p`.  In `F_p`
define the symmetric odd interval

```text
O_L={+-1,+-3,...,+-L}.                                   (2)
```

Let

```text
D_h={epsilon*a/b in Q* :
       epsilon in {+-1}, a,b>=1, gcd(a,b)=1,
       a and b have opposite parity, a+b<=h}.             (3)
```

Reduce the fractions in `(3)` modulo `p`.  Then

```text
F_p^* minus (O_L/O_L) = D_h (mod p).                      (4)
```

If in addition

```text
p>2(h-1)^3,                                               (5)
```

the largest number of pairwise-disjoint multiplicative dilates of `O_L` is

```text
omega_h=1+floor(log_2(h-1)).                              (6)
```

Equivalently, the Cayley graph on `F_p^*` with

```text
u~v iff uO_L and vO_L are disjoint                        (7)
```

has exact clique number `(6)`.

The two thresholds are explicit sufficient thresholds, not asserted
minimal.  Formula `(4)`, unlike an asymptotic density estimate, is an exact
set identity.

## 2. Inheritance and connection ledger

The closest proved mechanism is THM-3420's multiplicative splitter reduction:
once scalar block sizes are tight, exact translate incidence is the missing
coordinate.  Its canonical hostile is the order-`211` fixed-zero interval,
whose ratio set is already all of `F_211^*`.  The corrected near miss is the
rank-six density compiler: at seven blocks density reaches exactly `1/7`, so
no finite order cutoff survives.  The least-used sidecar is numerator parity;
it turns a centered interval into the odd lattice coset `(1,1)+2Z^2`.

| item | exact content |
|---|---|
| source | multiplicative dilates of the odd interval `(2)` |
| target | the bounded rational connection set `(3)` |
| map | an eight-point/general `(h+1)`-point circular gap followed by an odd Bezout lift |
| preserved | exact disjointness and multiplicative ratios |
| destroyed | the locations of individual interval points and all non-pairwise overlap multiplicities |
| restored sidecar | `v_2` of the rational ratio |
| cheapest decisive test | `h=7`, where `D_7` has 24 elements and clique `{1,2,4}` |

This is a lawful graph: vertices are dilates, adjacency is the intrinsic
pairwise observable “disjoint,” there are no ties, and the target is exactly
the packing number.  No tournament orientation is introduced.

## 3. The circular-gap and parity-constrained lift

Fix `lambda in F_p^*`.  Place

```text
0,lambda,2lambda,...,h lambda
```

on the cyclic `p`-gon.  The `h+1` gaps sum to `p`, so two of the displayed
points give, after choosing orientation,

```text
b lambda=+-A (mod p),
1<=b<=h,       1<=A<=floor(p/(h+1)).                     (8)
```

Divide `A,b` by their gcd.  This is lawful because `p` is prime and the gcd
is less than `p`.  Continue to write the coprime pair as `A,b`.  Assumption
`(1)` implies `A,b<=L`.

If `A,b` are both odd, `(8)` already puts `lambda` in `O_L/O_L`.  Otherwise
they have opposite parity.  If `A+b<=h`, `(8)` puts `lambda` in `D_h`.

It remains to treat `A+b>=h+1`.  Since `gcd(A,b)=1`, the equation

```text
b x-A y=p                                                  (9)
```

has integer solutions.  Because `p` is odd and `A,b` have opposite parity,
its odd--odd solutions form one nonempty class

```text
(x,y)=(x_0,y_0)+2t(A,b),       t in Z.                   (10)
```

The real line `(9)` contains `(p/(A+b),-p/(A+b))`.  Choosing the nearest
integer `t` in the parity class `(10)` gives

```text
max(|x|,|y|)<=p/(A+b)+max(A,b).                           (11)
```

For all pairs in the current range,

```text
p/(A+b)+max(A,b)<=p/(h+1)+(h+1).                         (12)
```

Indeed, if `A<=b`, this is immediate from `b<=h`.  If `A>=b`, the left side
is the convex function `p/(A+b)+A` on

```text
max(b,h+1-b)<=A<=p/(h+1);
```

at its lower endpoint it is at most `p/(h+1)+h`, and at its upper endpoint
it is strictly less than `p/(h+1)+(h+1)`.  Finally `(1)` is exactly the
inequality

```text
p/(h+1)+(h+1)<=L.                                        (13)
```

Thus `(9)` supplies odd `x,y` with `|x|,|y|<=L`, and `(8)--(9)` give
`lambda=x/y in O_L/O_L`.

Conversely, suppose a reduced opposite-parity fraction
`epsilon*a/b` with `a+b<=h` belonged to `O_L/O_L`.  Odd `x,y` with
`|x|,|y|<=L` would make

```text
b x-epsilon*a y
```

an odd multiple of `p`.  It is nonzero but has absolute value at most
`(a+b)L<=hL<p`, a contradiction.  This proves `(4)`.

## 4. The exact dyadic clique law

Normalize a clique from `(7)` by multiplying every vertex by the inverse of
one vertex.  All remaining vertices lie in `D_h`.  If three rational
fractions from `D_h` satisfy one ratio relation modulo `p`, clearing
denominators produces an integer of absolute value at most `2(h-1)^3`.
Under `(5)` it must vanish over the integers.  Hence every normalized modular
clique is already a clique in the rational Cayley graph with connection set
`D_h`.

Put

```text
t=1+floor(log_2(h-1)).                                    (14)
```

Color `Q^*` by

```text
color(q)=v_2(q) mod t.                                    (15)
```

Every element of `D_h` has nonzero 2-adic valuation of absolute value at
most `t-1`: its reduced numerator and denominator have opposite parity and
are at most `h-1`.  Thus `(15)` is a proper `t`-coloring and every clique has
size at most `t`.

The bound is attained by

```text
{1,2,4,...,2^(t-1)}.                                     (16)
```

Every ratio in `(16)` is a signed power `2^j/1` or its reciprocal, and
`2^(t-1)+1<=h`, so it belongs to `D_h`.  This proves `(6)`.

## 5. The `h=7` half-clock specialization

For `h=7`,

```text
D_7={+-a/b : gcd(a,b)=1, a,b opposite parity, a+b<=7},
|D_7|=24,        omega_7=3.                              (17)
```

For prime half-clock orders `p=14k+r` with `r in {1,9,11,13}`, the odd block
has the form `(2)` with

```text
(r,L,c)=(1,2k-1,8),(9,2k+1,2),(11,2k+1,4),(13,2k+1,6).  (18)
```

Thus every such prime `p>=512` has at most three pairwise-disjoint odd
blocks.  This is the exact large-prime incidence input needed by the
rank-seven half-twist classification; scalar overlap accounting and the
finite primes below `512` are separate obligations.  This theorem alone does
not classify a cover, a composite modulus, an arbitrary physical time, or
LRC(14).

## 6. Equality and failure boundaries

- The obstruction set is precisely the opposite-parity short fractions.
  Dropping parity changes the theorem: odd/odd short fractions lie in the
  ratio set, not its complement.
- The clique lower bound `(16)` shows the logarithmic dyadic tariff is sharp.
- The proof needs primality when dividing the short circular-gap relation.
  A composite modulus sharing a factor with `A,b` can descend to a proper
  quotient; no composite analogue is claimed.
- Condition `(5)` prevents accidental modular identities among distinct
  small rational ratios.  Without it, the rational `v_2` coloring need not
  descend.
- Pairwise disjointness forgets higher overlap multiplicity.  A cover theorem
  must retain a separate orbit-overlap budget.

## 7. Exact companion contract

The standard-library companion checks `60,049` admissible rounding cells on
prime controls through `h=24`; constructs and colors the rational graphs
through `h=40`; verifies the sharp power-of-two cliques; finds prime
finite-field controls satisfying `(1),(5)` and compares `O_L/O_L` with
`D_h` exactly; and freezes the four `h=7` prime residue controls.  It uses
integer and `Fraction` arithmetic only, with no assertion-dependent gate,
float, external package, subprocess, or file write.

Reproduce with

```text
python3 04-computation/odd_interval_ratio_complement_thm3423.py
python3 -O 04-computation/odd_interval_ratio_complement_thm3423.py
```
