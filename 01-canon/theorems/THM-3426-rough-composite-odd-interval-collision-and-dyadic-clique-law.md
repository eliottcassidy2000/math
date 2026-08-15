---
id: THM-3426
title: "Rough-composite odd-interval collision and dyadic clique law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For an odd modulus n
  whose least prime factor exceeds the
  short-relation height h, the complement of the unit-multiplier collision
  set of the maximal odd interval is exactly the opposite-parity rational
  packet D_h.  Beyond the same cubic threshold as THM-3423, pairwise-disjoint
  unit dilates have sharp dyadic clique number 1+floor(log_2(h-1)).  This is a
  rough-modulus incidence theorem, not a mixed quotient-order cover theorem.
source: root-rough-composite-odd-interval-2026-08-15
audit: >
  Independent proof reconstruction; rough controls at h=2,7,10,16; exact
  n=969 short-gcd hostile; normal/optimized/isolated/stored replay; LF hash,
  AST safety, scope, and documentation audit clean.
depends_on:
  - THM-3423-odd-interval-ratio-complement-and-dyadic-clique-law
related:
  - THM-3421-prime-half-twist-rank-seven-classification
  - THM-3425-half-twist-rank-six-primitive-breaker-profile-closure
script: 04-computation/rough_composite_odd_interval_collision_thm3426.py
output: 05-knowledge/results/rough_composite_odd_interval_collision_thm3426.out
script_sha256: b96eeca089ac86e2d82f0d52fef37ca8b509d5ffd76615ed08524748791d3ced
output_sha256: 97fb799cf97d6485ef48aac93caa486c223d147207a4d5a72bc5f52ff2becb2c
semantic_sha256: f03aded979fb19d75378b50614b6df1dc100ed6e25d08cd482501810ae72a691
hash_basis: LF-normalized bytes
---

# THM-3426 -- rough-composite odd-interval collision and dyadic clique law

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Candidate statement

Let `h>=2`, let `n` be odd, and let `L` be a positive odd integer.  Put

```text
c=n-hL
```

and assume

```text
0<c<2h,
spf(n)>h,
n >= (h+1)(h(h+1)+c).                                  (1)
```

Here `spf(n)` is the least prime factor of `n`; the condition includes the
prime case and is new when `n` is composite.  In `Z/nZ` put

```text
O_L={+-1,+-3,...,+-L}.                                  (2)
```

The entries of `(2)` are **not** assumed to be units.  Let `U_n=(Z/nZ)^*`
act on `(2)` and define its unit-multiplier collision set

```text
C_n(O_L)={lambda in U_n: lambda O_L intersects O_L}.     (3)
```

Let

```text
D_h={epsilon a/b:
       epsilon in {+-1}, a,b>=1, gcd(a,b)=1,
       a,b have opposite parity, a+b<=h}.                (4)
```

Because `spf(n)>h`, every numerator and denominator in `(4)` is a unit
modulo `n`.  Reduce `(4)` in `U_n`.  Then

```text
U_n minus C_n(O_L)=D_h (mod n).                          (5)
```

If also

```text
n>2(h-1)^3,                                              (6)
```

the largest number of pairwise-disjoint unit dilates of `O_L` is exactly

```text
1+floor(log_2(h-1)).                                     (7)
```

The thresholds are explicit sufficient thresholds, not asserted minimal.
For prime `n`, every entry of `(2)` is a unit and `(5)` is exactly THM-3423's
ratio-complement law.  For composite `n`, writing `O_L/O_L` would generally
be ill-typed; `(3)` is the faithful carrier.

## 2. Inheritance and connection ledger

The closest proved mechanism is THM-3423.  Its circular gap, odd Bezout lift,
and dyadic coloring use primality only to cancel a short common divisor and
to interpret every interval entry as a quotient denominator.  The first use
survives under `spf(n)>h`; the second must be replaced by the action relation
`(3)`.

The canonical hostile is

```text
(h,n,L,c)=(9,969,107,6),       spf(n)=3.                 (8)
```

The size threshold in `(1)` holds, but even after discarding rational
fractions with nonunit numerator or denominator, the true complement has the
four additional unit multipliers

```text
{161,325,644,808}.                                       (9)
```

Thus the roughness assumption is an actual gcd-cancellation gate, not a
convenience.

| field | exact connection |
|---|---|
| source | the unit action of `U_n` on the possibly nonunit interval `O_L` |
| target | the rational Cayley graph with connection packet `D_h` |
| map | shortest circular gap, cancellation of its gcd, then the parity-constrained Bezout lift |
| preserved | exact pairwise collision/disjointness, strict interval endpoints, and multiplicative normalization |
| destroyed | the colliding point pair and the small-prime gcd arm |
| restoration sidecar | `spf(n)>h`; below it one needs a breaker profile rather than a scalar density |
| cheapest tests | rough composite `n=517,h=7`; small-prime hostile `(8)` |

The graph is intrinsic: vertices are unit dilates, adjacency is disjointness,
and there are no ties.  No tournament orientation is introduced.

## 3. The composite circular gap

Fix `lambda in U_n`.  The `h+1` residues

```text
0,lambda,2lambda,...,h lambda
```

are distinct: a collision would make `n` divide a nonzero integer of absolute
value at most `h`.  Their circular gaps sum to `n`, so two displayed points
give

```text
b lambda=epsilon A (mod n),
epsilon in {+-1}, 1<=b<=h, 1<=A<=floor(n/(h+1)).         (10)
```

Put `delta=gcd(A,b)`.  Since `delta<=h<spf(n)`, it is a unit modulo `n`, and
we may cancel it in `(10)`.  Renaming the resulting coprime pair as `(A,b)`,
condition `(1)` implies

```text
A,b<=L.                                                  (11)
```

If `A,b` are odd, `(10)` itself exhibits a collision in `(3)`.  Otherwise
they have opposite parity.  When `A+b<=h`, `(10)` places `lambda` in `(4)`.

It remains to handle `A+b>=h+1`.  The equation

```text
b x-A y=n                                                (12)
```

has odd--odd solutions.  Indeed, opposite parity of `A,b` and oddness of `n`
force one variable odd, and the ordinary solution step `(A,b)` lets us choose
the other odd; the odd--odd solutions then form one class with step
`2(A,b)`.  Choosing the nearest member of that class to the real point
`(n/(A+b),-n/(A+b))` gives

```text
max(|x|,|y|)<=n/(A+b)+max(A,b)
             <=n/(h+1)+(h+1)<=L.                        (13)
```

The middle inequality is the convex endpoint estimate from THM-3423, and
the last is equivalent to `(1)`.  Equations `(10)` and `(12)`, followed by
cancellation of the unit `b`, give

```text
lambda y=epsilon x (mod n),
```

so `lambda` lies in `(3)`.

Conversely, take `lambda=epsilon a/b` from `(4)`.  If it caused a collision,
odd `x,y` with `|x|,|y|<=L` would make

```text
b x-epsilon a y
```

a multiple of `n`.  It is nonzero by the opposite parity of `a,b`, but its
absolute value is at most `(a+b)L<=hL<n`, a contradiction.  This proves
`(5)`.

## 4. Modular descent and the sharp clique

Normalize a family of pairwise-disjoint unit dilates by one member.  Formula
`(5)` puts every other multiplier in `D_h`.  For three reduced fractions from
`D_h`, a modular ratio identity clears to an integer

```text
a d f-c b e
```

divisible by `n`, with absolute value at most `2(h-1)^3`.  Under `(6)` it
must vanish over the integers.  Hence every normalized modular clique is a
clique in the rational Cayley graph with connection set `D_h`.

Put

```text
t=1+floor(log_2(h-1)).                                   (14)
```

The coloring `q -> v_2(q) mod t` is proper on that graph: every connection
has nonzero 2-adic valuation of absolute value at most `t-1`.  The clique

```text
{1,2,4,...,2^(t-1)}                                     (15)
```

is sharp.  Its entries are units because they are at most `h-1<spf(n)`, and
all pair ratios belong to `D_h`.  This proves `(7)`.

## 5. Exact companion

The standard-library companion performs the following exact checks.

- It reconstructs all `9,160` rational connection cells through `h=40`, the
  dyadic coloring, and all `342` sharp-clique edges.
- For each `2<=h<=16`, it finds the first two composite moduli satisfying
  `(1),(6)` and `spf(n)>h`: `30` controls in total.  Every interval contains
  nonunits, so the test exercises the new carrier rather than a disguised
  field case.
- It directly computes the unit action and `(3)` in `1,032,666` checked action
  cells, verifies `(5)`, and independently compares `42,184` modular graph
  pairs with their rational lifts.
- It checks `22,201` parity-lift rounding cells and freezes hostile `(8)--(9)`.

Normal and optimized runs are byte-identical to the stored output.  The
companion is a hostile referee for the algebraic proof, not bounded evidence
for arbitrary moduli.

## 6. Boundaries and non-consequences

- The interval may contain nonunits; this is intentional.  Only the
  multiplier is required to lie in `U_n`.
- The least-prime-factor gate can fail exactly through short-gcd arms, as
  `(8)--(9)` demonstrates.  A profile that records those arms is additional
  data, not a corollary of `(5)`.
- At `h=2`, `D_h` is empty and the sharp packing number is one; the companion
  includes composite controls `27,33`.
- The numeric thresholds are sufficient and may not be minimal.
- At `h=7`, this theorem gives the same sharp three-disjoint-odd-block tariff
  on every sufficiently large rough modulus.  It does not combine different
  quotient orders, prove a Boolean cover impossible, or settle the reserved
  THM-3425 breaker-profile lane.
- No physical common-time or LRC(14) conclusion follows.
