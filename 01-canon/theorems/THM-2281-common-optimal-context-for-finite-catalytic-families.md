---
id: THM-2281
title: "Common optimal context for finite catalytic families"
status: >
  PROVED. In every commutative nonexpansive integer-valued metric monoid,
  any finite family of nonempty fixed-saving context loci has a common
  context: add one witness from each upper ideal. In particular one context
  simultaneously attains the catalytic root length of every member of any
  finite family. For knots, Schubert prime coordinates sharpen the sum to
  the coordinatewise maximum of the individual witness vectors. On their
  finite union prime alphabet, the common optimal locus contains a full
  translated orthant, has box density one, and has a finite antichain of
  minimal contexts. The theorem is noneffective, does not produce a positive
  catalyst, and does not extend to arbitrary infinite families or to
  nonattained real-valued infima.
source: codex-2026-07-25-common-optimal-context
depends_on:
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity
related:
  - THM-2272-persistent-interaction-packing-and-calibrated-catalytic-defect-spectrum
---

# THM-2281 -- one context realizes every optimum in a finite family

THM-2191 proves that catalytic response is antitone when a common summand is
added. Its immediate finite-intersection consequence is stronger than
choosing a separate optimal context for every object.

## 1. Abstract simultaneous-context theorem

Let `(M,+,0)` be a commutative monoid with an integer-valued metric `d`
satisfying joint nonexpansivity

```text
d(a+c,b+e)<=d(a,b)+d(c,e).                           (1)
```

Use THM-2191's notation

```text
ell(x)=d(x,0),
rho_x(y)=d(x+y,y),
ell_cat(x)=min_y rho_x(y).                           (2)
```

For an integer saving threshold `s>=0`, put

```text
I_s(x)={y:rho_x(y)<=ell(x)-s},                       (3)

I_opt(x)={y:rho_x(y)=ell_cat(x)}.                    (4)
```

THM-2191 proves

```text
y in I_s(x)  implies  y+z in I_s(x) for every z,    (5)
```

and proves that `I_opt(x)` is a nonempty upper ideal. The attainment in the
last assertion uses integrality of the metric.

> **Finite common-context theorem.** Let `x_1,...,x_m` be a finite
> labelled family. For each `h`, let `J_h` be either `I_opt(x_h)` or a
> nonempty fixed-saving locus `I_(s_h)(x_h)`. Then
>
> ```text
> intersection_(h=1)^m J_h is nonempty.              (6)
> ```

Choose `y_h in J_h` and set

```text
y=y_1+...+y_m.                                       (7)
```

For each `h`,

```text
y=y_h+sum_(j!=h)y_j.
```

Equation (5) therefore gives `y in J_h`. This holds for every `h`, proving
(6). No compatibility between the separately chosen witnesses is needed;
commutative addition supplies it.

In particular,

```text
rho_(x_h)(y)=ell_cat(x_h)             for every h.   (8)
```

Thus one common context simultaneously realizes all catalytic optima in a
finite family.

The same proof works for any finite collection of additive upper ideals in
a commutative monoid. The metric hypotheses are the mechanism which makes
the response loci upper ideals and, for optimal loci, makes them nonempty.

## 2. Knot corollary and the coordinatewise maximum

Specialize to the connected-sum monoid of oriented knots with Gordian
distance. Write

```text
rho_K(J)=d_G(K#J,J),
u_cat(K)=min_J rho_K(J).                             (9)
```

For any finite family `K_1,...,K_m`, choose individual optimal contexts
`J_h`. Equation (7) becomes connected sum:

```text
J=J_1#...#J_m.                                      (10)
```

Then

```text
d_G(K_h#J,J)=u_cat(K_h)               for every h.  (11)
```

Schubert prime decomposition, used in THM-2242, sharpens (10). Let
`P_1,...,P_r` be the finite union of the prime-knot types occurring in the
chosen `J_h`, and write

```text
J_h=P_1^(#a_(h,1))#...#P_r^(#a_(h,r)),
a^(h) in N^r.                                        (12)
```

Define the coordinatewise maximum

```text
b_i=max_h a_(h,i),                  1<=i<=r.          (13)
```

For every `h`, one has `b=a^(h)+c^(h)` for some `c^(h) in N^r`. Hence the
knot

```text
J_b=P_1^(#b_1)#...#P_r^(#b_r)                        (14)
```

lies above every `J_h` in the connected-sum divisibility order. The
upper-ideal law gives

```text
d_G(K_h#J_b,J_b)=u_cat(K_h)            for every h. (15)
```

This maximum may be strictly smaller than the sum vector in (10). It is a
faithful-section improvement specific to unique prime factorization; no
claim of global size minimality is made.

## 3. Orthant, density, and finite minimal bank

Retain the alphabet `P_1,...,P_r`. Equation (15) and upper closure show that
the common optimal locus contains

```text
b+N^r
 subset intersection_(h=1)^m I_opt(K_h).             (16)
```

Inside the box

```text
B_n={0,1,...,n}^r,
```

the complement of this translated orthant is contained in the union of the
coordinate slabs `{a:a_i<b_i}`. Therefore, for every `n>=max_i b_i`,

```text
#(B_n minus (b+N^r))
 <=sum_(i=1)^r b_i (n+1)^(r-1).                     (17)
```

Dividing by `#B_n=(n+1)^r` proves that the common optimal locus has box
density one on this finite prime alphabet.

It is itself an upper ideal: intersections preserve upper closure.
Dickson's lemma, proved directly in THM-2191, therefore gives a finite
antichain of minimal common optimal contexts in `N^r`. This is a finite
existence bank, not an effective computation of its entries.

## 4. Sharp boundaries and information ledger

Finiteness is load-bearing. In the additive monoid `N`, the upper ideals

```text
J_h={n:n>=h}
```

have nonempty every finite intersection, but

```text
intersection_(h>=0)J_h=empty.                        (18)
```

Attainment is also load-bearing for an optimal-context statement.
THM-2191's real-valued metric example has a diagonal response whose infimum
is not attained, so `I_opt` is empty even for one object. Fixed-saving
loci remain covered by (6) whenever each one is explicitly nonempty.

The typed connection is:

```text
source:
  one individually witnessed upper response ideal per labelled object;

target:
  one context satisfying every requested finite response bound;

map:
  add the witnesses, or take their coordinatewise maximum in a faithful
  finite prime-factor section;

preserved:
  every labelled saving threshold, including exact catalytic optimality;

destroyed:
  the first or smallest common context and any effective size bound;

sidecar:
  nonemptiness of each requested locus, and a finite prime alphabet for
  the maximum, density, and antichain refinements;

cheapest hostile tests:
  an empty positive-saving ideal, a nonattained real infimum, or the
  infinite descending intersection (18).                            (19)
```

The theorem produces no positive knot catalyst. If every individual optimum
is the root value, the common context merely realizes those trivial optima.
It also computes neither `u_cat(K)` nor a minimal common context. Its content
is the exact finite Helly-number-one behavior created by additive upper
ideals. QED.
