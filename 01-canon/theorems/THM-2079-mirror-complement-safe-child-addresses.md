---
id: THM-2079
title: "Mirror-complement safe-child addresses and equivariant outer-tail ownership"
status: >
  PROVED from THM-2073/2075. Reversal acts freely on every divisor-complete
  quotient safe set, so components occur in mirror pairs and their count is
  positive even. Every unique safe-child section is reversal-equivariant;
  over a depth-r tower, mirror components have complementary addresses
  a and 2^r-1-a. Terminal endpoint-owner sets mirror unchanged, while each
  original odd tail flips its nearest-integer owner parity. Thus the outer
  two-tail covering condition is equivariant, not contradicted by mirror
  parity. Address search may be halved, but LRC(14) is not proved.
source: codex-2026-07-21-LRC-mirror-address-correction
depends_on:
  - THM-2073
  - THM-2075
related:
  - THM-2078
  - MISTAKE-238
  - HYP-8845
---

# THM-2079 -- mirror-complement safe-child addresses

Put `delta=1/14` and let

```text
iota(t)=1-t=-t in R/Z.                                  (1)
```

Retain the THM-2073 tower

```text
C=Q_0,
Q_i=2Q_(i+1) union {h_i},    0<=i<r.                    (2)
```

Every `Q_i` is divisor-complete through `14`, and THM-2075 gives a unique
safe-child inverse section

```text
s_i:G_(Q_(i+1))->G_(Q_i)                                (3)
```

of doubling.

## 1. Free mirror pairing of components

Every safe set is reversal-invariant:

```text
iota(G_(Q_i))=G_(Q_i),                                  (4)
```

because `||q(1-t)||=||qt||`. The only fixed points of `iota` on the circle
are `0` and `1/2`. The first is never safe at positive threshold. Since
`Q_i` contains a multiple of `2`, it has an even speed `q`, and

```text
||q/2||=0,
```

so `1/2` is not safe either. Thus reversal acts freely on `G_(Q_i)`.

It also acts freely on its connected components. A component fixed setwise
would be either a fixed singleton or a compact interval on which an
involution has a fixed point; both contradict freeness on points. Therefore
the components occur in distinct pairs

```text
I <-> iota(I),                                          (5)
```

and

```text
#components(G_(Q_i)) is positive and even.              (6)
```

Positivity follows independently from settled lower-dimensional LRC. This is
the valid mirror-parity conclusion retained after MISTAKE-238.

## 2. Every safe-child section is equivariant

Doubling commutes with reversal:

```text
D(iota(theta))=iota(D(theta)).                           (7)
```

If `theta=s_i(sigma)` is the unique safe child over `sigma`, then
`iota(theta)` is safe by (4) and doubles to `iota(sigma)` by (7). Uniqueness
in (3) forces

```text
s_i(iota(sigma))=iota(s_i(sigma)).                      (8)
```

On a component `I`, write THM-2075's constant sheet bit as

```text
s_i(sigma)=(sigma+epsilon_I)/2,
epsilon_I in {0,1}.                                     (9)
```

For `0<sigma<1`, equations (8)--(9) give

```text
epsilon_(iota I)=1-epsilon_I.                           (10)
```

Indeed

```text
1-(sigma+epsilon_I)/2
 =(1-sigma+(1-epsilon_I))/2.
```

Thus mirror components take opposite child sheets at every level.

## 3. Complementary depth-r addresses

Let `I` be a terminal component of `G_(Q_r)`. THM-2075 supplies the composite
inverse

```text
s_I(sigma)=(sigma+a_I)/2^r,
0<=a_I<2^r.                                             (11)
```

Iterating (8), or applying reversal directly to (11), gives

```text
a_(iota I)=2^r-1-a_I.                                   (12)
```

The right side is the bitwise complement of the `r`-bit word for `a_I`.
Therefore the entire address assignment is determined by one component from
each mirror pair. A finite terminal address SAT or owner-word search can be
halved without losing any candidate.

Terminal endpoint-owner sets are paired without changing speed labels. If
`q in Q_r` owns an endpoint `sigma`, then

```text
||q iota(sigma)||=||q sigma||=1/14.                     (13)
```

After lifting, the inherited owner `2^r q` remains paired in the same way.
Guard ties may be added, as allowed by THM-2075, but the inherited owner set
is nonempty on both endpoints.

## 4. Why mirror parity does not close the outer tails

Let `z` be either original odd tail `x` or `y`. On `G_C`, strict seam failure
gives `||zt||<1/7`, so the nearest integer `N_z(t)` is unique. At the mirror
phase,

```text
N_z(1-t)=z-N_z(t).                                      (14)
```

Since `z` is odd,

```text
N_z(1-t)=1-N_z(t) (mod 2).                              (15)
```

Thus each tail kills the complementary outer lift on the mirror component.
If `x` and `y` have opposite owner parity at `t`, they still have opposite
parity at `1-t`; both bits flip together. The two-tail covering condition is
therefore exactly reversal-equivariant.

This is the precise correction to HYP-8920. Mirror parity does not transport
the empty full-row safe set to an empty terminal. Instead it says:

```text
terminal components and addresses come in complementary pairs,
and the outer tail owners cover those pairs equivariantly.                (16)
```

The symmetry halves the remaining address/owner search but supplies no
contradiction by itself.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that a free involution must obstruct a cover.
Here it organizes the cover: terminal components are the vertices, reversal
pairs them, child addresses complement, and odd-tail owner bits complement in
the same way. Quotienting by reversal preserves one representative component
and its address, but the complementary reconstruction rule (12) must be kept.

Orienting the two members of each mirror pair by smaller numerical address
gives disjoint two-vertex transitive tournaments, each with score histogram
`(0,1)` and one Hamiltonian path. Those fingerprints merely restate (12).
The faithful carrier is the mirror-paired component graph with address and
tail-owner bits. QED.
