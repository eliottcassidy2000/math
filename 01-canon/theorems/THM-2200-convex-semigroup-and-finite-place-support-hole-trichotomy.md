---
id: THM-2200
title: "Convex, semigroup, and finite-place support-hole trichotomy"
status: >
  PROVED, with exact LRC, tournament, and flow applications. For every
  nonnegative integer polynomial, an integer exponent outside the convex hull
  of its support remains absent on every dilation, while an exponent inside
  that hull appears in some power after a finite dilation. Nevertheless, for
  every prime p, a degree-one support hole has coefficient divisible by p at
  its p-fold dilation. Applied to THM-2192, the 17 linearly separated matching
  holes are persistent, the other 19 fill in Hafnian powers of degree at most
  20, and all prime dilates have the Frobenius divisibility. Applied to
  THM-509, the unrealized order-six tournament profile (8,10) fills after a
  two-fold ordered join as the realized order-twelve profile (16,20).
  THM-2177's K_3 flow obstruction lies on the convex-external side. These are
  operation-sensitive classifications, not transfers of the underlying
  objects and not a proof of LRC(14).
source: codex-2026-07-24-support-hole-spectrum
depends_on:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-509-the-two-layer-baby-hodge-dichotomy
  - THM-2177-planar-counterexample-to-goemans-unsplittable-cost-flow-conjecture
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
---

# THM-2200 -- the support-hole trichotomy

Let

```text
P(x)=sum_(a in A) c_a x^a,
A subset Z^d finite,                 c_a in Z_(>0),  (1)
```

be a Laurent polynomial with nonnegative integer coefficients. Write
`A^(+n)` for the `n`-fold Minkowski sum of `A`. Positivity gives the exact
support law

```text
supp(P^n)=A^(+n).                                     (2)
```

Fix an integer exponent `h in Z^d` absent from `A`.

## 1. Convex-external holes persist

If `h notin conv(A)`, rational hyperplane separation supplies a rational
linear functional `ell` and rational constant `C` such that

```text
ell(a)<=C for every a in A,          ell(h)>C.        (3)
```

Every exponent `b=a_1+...+a_n` in `A^(+n)` obeys

```text
ell(b)<=nC<n ell(h)=ell(nh).                          (4)
```

Therefore

```text
nh notin supp(P^n)                  for every n>=1.  (5)
```

This is a persistent support hole, certified by a linear inequality.

## 2. Convex-internal holes fill after finite dilation

If `h in conv(A)`, the rational polytope

```text
{(theta_a):theta_a>=0, sum_a theta_a=1,
             sum_a theta_a a=h}                      (6)
```

has a rational point. Choose a common denominator `n` and put
`m_a=n theta_a`. Then

```text
m_a in Z_(>=0),       sum_a m_a=n,
nh=sum_a m_a a.                                        (7)
```

Thus `nh in A^(+n)`, and positivity in (1) gives

```text
[x^(nh)]P^n>0.                                        (8)
```

Equations (5) and (8) prove the exact ordinary-support dichotomy:

> An integer ray through a degree-one support hole meets the support
> semigroup at some positive degree if and only if the hole lies in
> `conv(A)`.

The least successful `n` is an arithmetic denominator of the hole, invisible
to the convex body alone.

## 3. Finite-place survival is strictly finer than support

For every prime `p`, Frobenius in `F_p[x_1^(+-1),...,x_d^(+-1)]` gives

```text
P(x)^p=P(x_1^p,...,x_d^p)                 modulo p.  (9)
```

Because `h notin A`, the coefficient of `x^(ph)` on the right is zero.
Hence

```text
[x^(ph)]P^p=0                              modulo p. (10)
```

This remains true when `ph in supp(P^p)` and its integer coefficient is
positive. Ordinary semigroup saturation and finite-place survival are
therefore different predicates. Equation (10) is the abstract positive-model
shadow of THM-2022's whole-layer Frobenius mechanism; it does not identify
`P` with the NC2 balanced-face polynomial.

## 4. Scalar matching/Hafnian application

Take THM-2192's coloured perfect-matching polynomial

```text
P_S(x_1,...,x_6)
 =sum_(perfect matchings M of S)
       product_({u,v} in M) x_(ell(u,v)).             (11)
```

Its 36 missing guard-danger profiles split exactly as follows.

- Seventeen violate one of THM-2192 (21m). Summing the relevant endpoint
  inequality over `n` matchings gives its right side multiplied by `n`, so
  (5) applies and every dilation remains absent from `P_S^n`.
- The other nineteen have the rational certificates in THM-2192 (21o).
  Clearing their denominators gives successful powers

```text
(4,4,4,4,4,2,20,5,5,2,4,4,4,20,12,2,4,5,5),        (12)
```

  in the displayed hole order. Hence all nineteen fill by degree at most
  twenty.
- For every prime `p`, (10) says that the coefficient at the `p`-dilate of
  any of the 36 holes is divisible by `p`.

For the smallest illustration,

```text
2(12333)=11333+22333.                                 (13)
```

The two ordered products of the matching witnesses on the right make the
coefficient at `2(12333)` positive and even. This algebraic product uses two
independent matching carriers. THM-2198's actual multiplication-by-thirteen
image pump acts on one rooted fibre and must retain occupancy, owner labels,
and winding. Equation (13) does not supply those transition sidecars.

## 5. Tournament ordered-join application

THM-509 proves that at order six the profile

```text
h=(c_3,c_5)=(8,10)                                   (14)
```

is unrealized, while profiles

```text
a_-=(8,8),                 a_+=(8,12)                 (15)
```

are realized. Thus `2h=a_-+a_+`, the exact analogue of (13).

Choose realizing tournaments `T_-` and `T_+` and form their ordered join
`T_- vee T_+`, orienting every cross edge from the first block to the second.
No directed cycle can use a cross edge: after entering the second block
there is no edge back to the first. Consequently every directed cycle lies
inside one factor and, in particular,

```text
c_k(T_- vee T_+)=c_k(T_-)+c_k(T_+)       for k=3,5.  (16)
```

The order-twelve join therefore realizes

```text
(c_3,c_5)=(16,20)=2h.                                 (17)
```

This fills the profile ray after changing vertex order. It does not realize
`(8,10)` at order six. Vertex count and the ordered strongly-connected-block
word are the load-bearing operation sidecars; the profile alone forgets
them.

## 6. The unsplittable-flow obstruction is convex-external

In THM-2177, the three cheap-path selections have incompatibility graph
`K_3`. Every capacity-good integral routing satisfies the stable-set facet

```text
z_1+z_2+z_3<=1,                                      (18)
```

whereas the fractional marginals obey

```text
z_1+z_2+z_3=16/15>1.                                 (19)
```

Thus the flow obstruction is on the convex-external side of the taxonomy:
linear cost/facet duality separates it. Such a separator can address the
seventeen external scalar matching holes, but it cannot see the nineteen
internal holes, which require matching-semigroup or labelled-incidence data.

The three examples therefore share a support-hole spectrum, not a common
object:

```text
convex-external  -> persistent under independent products;
convex-internal  -> filled after an arithmetic dilation;
finite-place     -> coefficient valuation can obstruct after filling;
dynamic image    -> additionally needs the operation's transport sidecars.
                                                               (20)
```

Confusing these levels is exactly what makes convex relaxations, bare support,
or an unrelated multiplication law look stronger than they are. QED.
