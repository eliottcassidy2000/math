# Signed Collatz parameters: conserved content, primitive denominators, and nine cycles at `b=-5`

**Status: PROVED elementary parameter mechanisms + FINITE-EXACT bounded
cycle/trajectory censuses.** No complete classification of the cycles of
`3n+1`, `3n-1`, `3n+5`, or `3n-5` is claimed. No novelty claim is made for
the classical rational-cycle viewpoint or the standard cycle equation.

The result and its controls are in
[the script](../../04-computation/experiments/arithmetic_braids2_20260917_signed_cycles.py)
and [matching JSON](../../04-computation/experiments/arithmetic_braids2_20260917_signed_cycles.json).
The JSON distinguishes an exhaustive bounded-word universe from a separate
bounded-start trajectory universe; their cycle lists are not conflated.

## Inheritance, portfolio, and concept board

The closest proved mechanisms are the ordered carry equation and inverse
fibre in [the first-session Collatz note](arithmetic_braids_20260917_collatz.md),
and the prime-adic affine lift in
[the summand note](arithmetic_braids_20260917_summand.md).
The canonical hostile is the `5n+1` map, where local inverse-fibre fullness
coexists with multiple cycles. The corrected near miss is inferring a
complete cycle classification from a finite census. The least-used sidecar
is the gcd between the state and the additive parameter.

Anchor: understand the signs and `b=-5`. Niche: classify cycles by arithmetic
content rather than by their smallest visible number. Wildcard: compare
word-bounded and start-bounded searches; this exposed two period-17 cycles
outside the first search's universe.

| Live concept | Invariant / exact map | Lost coordinate or hostile |
|---|---|---|
| Parameter and state sign | `(b,n)->(-b,-n)` | This preserves forward time, not arrow reversal |
| State gcd with parameter | `gcd(n,b)` when `3` does not divide `b` | Forgetting it mixes scaled and primitive cycles |
| Ordered halving word | Rational denominator `q` and carry `B` | Equal length and total exponent do not determine a cycle |
| Inverse-fibre braid | `R_b(n)=4n+b` | Full triadic coverage does not distinguish basins |
| Positive/negative regions | A finite one-way sign-crossing set | A signed cycle need not be reachable from positive starts |
| Residues modulo the parameter | Principal ideal and a unit-subgroup coset | A quotient component does not determine an integer basin |

For external context, [Lagarias, *The 3x+1 Problem: An Overview*, section 4](https://arxiv.org/html/2111.02635v1)
explains the correspondence between rational `3x+1` cycles and integer
`3x+k` cycles after clearing odd denominators. The proof and all parameter
conventions used here are given explicitly below.

## 1. Reflection is conjugacy, not reversal

For odd nonzero `b` and an odd integer `n`, define

```text
U_b(n)=(3n+b)/2^v_2(3n+b),                               (1)
```

provided `3n+b!=0`; valuations use absolute values. For `3` not dividing
`b`, every odd signed integer is in the domain. If `3|b`, the exceptional
input `n=-b/3` would land at zero, outside the odd nonzero domain.

For every nonzero odd integer `d`,

```text
U_(db)(dn)=d U_b(n).                                    (2)
```

This follows because multiplication by odd `d` does not change the power
of two dividing the numerator. The case `d=-1` gives

```text
U_(-b)(-n)=-U_b(n).                                     (3)
```

Thus reflection sends each forward arrow at parameter `b` to a forward
arrow at parameter `-b`, preserving its halving exponent. For example,
`3->5` under `U_1` becomes `-3->-5` under `U_(-1)`.
It does not become `5->3`: indeed `U_(-1)(5)=7`.

The inverse **relation** instead consists of all solutions

```text
n=(2^k u-b)/3,          k>=1,          2^k u=b mod3.      (4)
```

It is generally infinitely multivalued. Under `U_1`, the target `1` has
predecessors `1,5,21,...`. Negation, a bijective change of coordinates,
cannot be confused with this reversal of the graph's arrows.

## 2. The gcd is a conserved arithmetic stratum

For every defined step, oddness of `b` gives the exact identity

```text
gcd(U_b(n),b)=gcd(3n,b).                                 (5)
```

Here and below gcd is positive. The numerator has gcd
`gcd(3n+b,b)=gcd(3n,b)` with `b`, and removing powers of two changes
nothing. In particular, if `3` does not divide `b`, then

```text
d=gcd(n,b) is constant along the whole orbit.            (6)
```

For general odd `b`, every prime component other than three is conserved,
while

```text
v_3(gcd(U_b(n),b))=min(v_3(n)+1,v_3(b)).                 (7)
```

Consequently, on a cycle the full power `3^v_3(b)` must already divide
every node. The gcd is constant there even when `3|b`.

**Cycle content theorem.** For any integer cycle, the following coincide:
`gcd(n_i,b)` at any node and the gcd of all its nodes. Call this number `d`.
It divides `b`, and dividing all nodes and the parameter by `d` gives a
cycle at parameter `b/d` with gcd one. Conversely odd dilation lifts every
such primitive cycle.

To prove equality with the node gcd, the conserved `d` divides every node.
Conversely any common divisor of the nodes divides
`b=2^k_i n_(i+1)-3n_i`, so it also divides `gcd(n_i,b)`.

For `b=+-5` there are exactly two possible strata: `d=5`, consisting
precisely of fivefold dilations of cycles at `b=+-1`, and `d=1`, containing
the genuinely primitive parameter-five cycles. This is an exact all-cycle
partition, even though neither piece has been completely enumerated.

## 3. Each ordered word carries its own minimal parameter

Fix a positive halving word `k_1,...,k_L`. Put `K_0=0`,
`K_i=k_1+...+k_i`, and

```text
B=sum_(i=0)^(L-1) 3^(L-1-i) 2^K_i,
Delta=2^K_L-3^L,
q=|Delta|/gcd(B,|Delta|).                                (8)
```

Then `B>0`, `Delta!=0`, and `q` is odd and coprime to three.
The fixed point of the composed affine maps is

```text
n_0=b B/Delta.                                          (9)
```

**PROVED iff.** This word is an exact signed odd integer cycle word at
parameter `b` if and only if `q|b`. Repeated words are allowed here;
their minimal period may be shorter than `L`.

Necessity is the reduced denominator in (9). For sufficiency, `q|b`
makes `n_0` odd integral. Under a cyclic rotation, the new carry is
`B'=(3B+Delta)/2^k_1`, an integer. Since `Delta` is odd and coprime to
three, `gcd(B',|Delta|)=gcd(B,|Delta|)`. Thus every rotated source is odd
integral and satisfies `3n_i+b=2^k_i n_(i+1)`. The odd successor makes
the valuation exactly `k_i`, proving sufficiency and closure.

This also proves that `q` is invariant under cyclic remarking. Cancelling
a repeated word does not change it because it is the reduced denominator
of the same rational orbit. At parameter `b`, the cycle content is

```text
d=|b|/q.                                               (10)
```

Indeed write `B/Delta=A/q` in lowest terms with signed `A`; then
`n=(b/q)A` and `gcd(A,q)=1`. Dividing states by `b` identifies every
cycle with a rational `U_1` cycle of exact denominator `q`.
The minimal positive parameter magnitude supporting that rational cycle
is therefore `q`. Parameter scaling, state content, and denominator are
three exact descriptions of the same invariant.

The rational clock condition is

```text
2^K_L=3^L mod q.                                       (11)
```

For primitive `b=+-5` cycles, `q=5`; since `3=2^3 mod5`, this becomes
`K_L=3L mod4`. This is necessary but not sufficient: the carry decides
whether the reduced denominator is actually five.

All members of a cycle have the same sign, namely the sign of `b/Delta`,
because every rotated carry is positive. Positive cycles at `b=-5`
therefore require `2^K<3^L`; negative cycles require the reverse inequality.

## 4. The inverse braid persists inside every stratum

When `3` does not divide `b`, a target has predecessors only if it is
coprime to three. Congruence (4) selects one parity of `k`, and consecutive
solutions satisfy

```text
R_b(n)=4n+b,
U_b(R_b(n))=U_b(n),
v_2(3R_b(n)+b)=v_2(3n+b)+2,
R_b^t(n)-n=(4^t-1)(3n+b)/3.                            (12)
```

Since `3n+b` is a three-adic unit,
`v_3(R_b^t(n)-n)=v_3(t)`. Every inverse fibre is again a full cycle on
odd residues modulo `2*3^s`, at every level `s`. Yet
`gcd(R_b(n),b)=gcd(n,b)`: it remains confined to its arithmetic stratum.
For `b=+-5`, inverse-fibre movement modulo five is simply `n->-n`.
Full coverage of the triadic coordinate does not erase this independent
mod-five information, nor does it choose a later basin.

## 5. Exact census: nine exhibited cycles at `b=-5`

The first universe is **all 1,744,435 ordered positive words** satisfying
`1<=L<=10` and `L<=K<=22`. There is no node-height or sign bound. Among
these, 37 words have `q=1` and 25 have `q=5`, including repeated and
cyclically remarked words. These give seven distinct `b=-5` cycles.

A separate direct iteration, using no cycle-word formula, tested all
10,000 signed odd starts `|n|<=9999` at each of `b=-5,-1,1,5`.
It found two more `b=-5` cycles, both of period 17 and total exponent 27,
outside the first universe. All 40,000 starts were resolved; none hit the
declared exploration cap. This is exactly why a bounded cycle count must
retain its period/height universe.

The full displayed `b=-5` list, with cyclic start chosen by least absolute
value, is:

| Cycle | `L` | `K` | Content `d` | Primitive parameter |
|---|---:|---:|---:|---:|
| `(5)` | 1 | 1 | 5 | -1 |
| `(25,35)` | 2 | 3 | 5 | -1 |
| `(85,125,185,275,205,305,455)` | 7 | 11 | 5 | -1 |
| `(-1)` | 1 | 3 | 1 | -5 |
| `(-5)` | 1 | 2 | 5 | -1 |
| `(-19,-31,-49)` | 3 | 5 | 1 | -5 |
| `(-23,-37,-29)` | 3 | 5 | 1 | -5 |
| The `-187` cycle below | 17 | 27 | 1 | -5 |
| The `-347` cycle below | 17 | 27 | 1 | -5 |

The complete longer witnesses are

```text
-187 -> -283 -> -427 -> -643 -> -967 -> -1453 -> -1091
 -> -1639 -> -2461 -> -1847 -> -2773 -> -2081 -> -781
 -> -587 -> -883 -> -1327 -> -1993 -> -187;

-347 -> -523 -> -787 -> -1183 -> -1777 -> -667 -> -1003
 -> -1507 -> -2263 -> -3397 -> -2549 -> -1913 -> -359
 -> -541 -> -407 -> -613 -> -461 -> -347.
```

Their ordered exponent words are respectively

```text
(1,1,1,1,1,2,1,1,2,1,2,3,2,1,1,1,5),
(1,1,1,1,3,1,1,1,1,2,2,4,1,2,1,2,2).
```

Both have `(L,K,q)=(17,27,5)` and
`Delta=2^27-3^17=5,077,565`, but their carries are `189,900,931` and
`352,383,011`. Substitution into (9) gives `-187` and `-347` exactly.
This is a sharp hostile to classifying cycles by period, total contraction,
and denominator alone: the ordered carry is essential.

Negation gives nine corresponding displayed cycles at `b=5`. For `b=-1`,
the same searches find positive cycles `(1)`, `(5,7)`,
`(17,25,37,55,41,61,91)` and negative `(-1)`; `b=1` is their reflection.
These counts are **FINITE-EXACT in the stated searches**, not universal
cycle counts. The elementary content theorem is universal; the catalogs
are not.

## 6. The one-way sign portal restricts which cycles are reachable

If `b>0`, the positive odd region is forward invariant. A negative input
can cross to positive only when `0<-n<b/3`. For `b<0` the negative region
is forward invariant, and a positive input can cross only when
`0<n<|b|/3`. Thus signs change at most once along an orbit.

For `b=-5`, the sole positive-to-negative portal is

```text
1 -> -1 -> -1.                                         (13)
```

No positive start can reach any other negative cycle. Multiples of five
cannot reach this portal because the gcd is invariant; they are conjugate
to positive `U_(-1)` orbits after division by five. Positive starts coprime
to five can reach `-1` only by first reaching `1`.

In the bounded-start control, all 4,000 positive odd starts at most 9999
and coprime to five do reach that portal; the remaining 1,000 positive
odd starts belong to the three displayed scaled positive cycles' basins.
The all-height extension is not proved. The negative primitive cycles
at `-19,-23,-187,-347` are real cycles, but inaccessible from any positive
start by the sign argument, independently of computation.

## 7. Annihilators and a stronger character invariant

Assume `gcd(b,6)=1`, write `M=|b|`, and work in `R=Z/MZ`. Every step is
multiplication by a unit in this ring:

```text
U_b(n)=3*2^(-k)*n mod M.                                (14)
```

Consequently the principal ideal `(n)` is invariant. Its integer encoding
is exactly the conserved content `d=gcd(n,b)`. A closing word with total
exponent `K` and length `L` satisfies

```text
n*(2^K-3^L)=0 mod M,
Delta in Ann_R(n)=(M/d),    equivalently M/d divides Delta. (15)
```

Indeed, `M|n*Delta` is equivalent to `M/d|Delta`, since `n/d` is coprime
to `M/d`. This is an exact bridge from a dynamical closing condition to a
zero-product condition in a finite ring. It is only a necessary closing
test: the carry still matters. At `b=5`, the word `(1,1)` has
`Delta=-5` and `B=5`, so the clock congruence modulo five holds, but its
only fixed point is `n=-5`, in content-five rather than the unit stratum.

The ideal does not retain all invariant modular information. Put `m=M/d`
and `a=n/d`, a unit modulo `m`. Let

```text
H_m=<2,3> inside (Z/mZ)^*.
```

The coset `a H_m` is invariant, because the next normalized state is
`3*2^(-k)*a mod m`. This is the complete component invariant of the
formal residue graph with an edge `a->3*2^(-k)*a` for every `k>=1`:
its strongly connected components are exactly the cosets of `H_m`.
To prove this, the multipliers for `k=1,2` have ratio two and generate
three as well. In a finite group their positive products include their
inverses, so they generate precisely `H_m`. Every such formal edge has
an actual odd integer lift: solve `n=a mod m` and
`3n+b/d=2^k mod 2^(k+1)` by the Chinese remainder theorem. This also
ensures the valuation is exactly `k`. Quotient connectivity does not
assert a cycle or a shared basin of the actual integer dynamics.
For `m=1` the assertion is the trivial one-element invariant.

At the parameter `b=23` (or `-23`), this becomes a familiar character.
Both `2=5^2 mod23` and `3=7^2 mod23` are quadratic residues, and two has
order eleven. Thus `H_23` is exactly the eleven nonzero squares, and
the Legendre symbol `(n/23)` is conserved on the coprime stratum.
For example, one and five have equal gcd content but opposite characters,
so their orbits cannot merge. Twenty-three is the smallest prime above
three where both two and three are squares: direct checks exclude
`5,7,11,13,17,19`. At modulus five, `H_5` is the entire unit group, so
this particular refinement adds nothing to the user's `b=-5` case.
It likewise gives no extra invariant for `b=1`, whose parameter ring is
trivial, and should not be transferred to the unrelated modulus-nine
output braid.

## Verification and next obligation

Run

```text
python 04-computation/experiments/arithmetic_braids2_20260917_signed_cycles.py
python -O 04-computation/experiments/arithmetic_braids2_20260917_signed_cycles.py
```

Both runs pass and produce byte-identical JSON. Checks remain active under
optimization. Beyond the exhaustive word universe and 40,000 direct starts,
the script performs 198 independent exact-fraction affine compositions,
25,992 gcd/dilation controls including parameters divisible by three, and
12,000 triadic-braid controls. The modular sidecar checks the generated
components for all odd moduli below 200 coprime to six, 19,132 character
steps for `b=+-23`, and 1,872 exact CRT edge lifts for moduli
`5,23,35,47,77` and exponents one through six. Every cycle includes explicit nodes,
exponents, carry, denominator, gcd content, and its primitive reduction.
The direct control caps only unresolved exploration segments at 2,000
steps and absolute height `10^80`; cached resolved paths retain actual
reachability. There were no censored starts.

The parameter-level mechanism now gives a better next target than counting
another few cycles: characterize which primitive denominator-five words
with `2^K<3^L` can exist, and whether every positive coprime-to-five orbit
must reach the unique sign portal. The congruence `K=3L mod4`, content
conservation, and triadic completeness are necessary structure, not a
resolution. The two period-17 cycles show why any proposed word quotient
must keep the carry or an equivalent ordered invariant.
