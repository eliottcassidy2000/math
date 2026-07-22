---
id: THM-2140
title: "Higher-order guard-pencil obligations force blocker depth"
status: >
  PROVED. In the remaining one-terminal-blocker guard pencils, clean 13^2-root
  fibres force integral higher-order generation. With six aligned
  nonblockers and one transverse terminal, the divided blocker lies in the
  integer span of the guard and every aligned triple. With seven aligned
  nonblockers, it lies in the span of the guard and every aligned five-set.
  Hence every corresponding determinant gcd divides the guard/blocker
  determinant. In nondecreasing valuation order, its 13-adic depth is at
  least a_4 for six aligned determinants and a_3 for seven. In particular
  the divided blocker is guard-aligned modulo thirteen, possibly with zero
  reduction.
  This is a recursive depth invoice, not an elimination of those two profiles
  or a proof of LRC(14).
source: codex-2026-07-22-LRC-higher-order-guard-pencil
depends_on:
  - THM-2125
  - THM-2139
related:
  - THM-2120
  - THM-2133
---

# THM-2140 -- higher-order guard-pencil depth invoices

Let `Gamma` be a rank-two character lattice, put

```text
K=Hom(Gamma,R/Z),                                     (1)
```

and let `g mod 13` be nonzero.  In both cases below there is exactly one true
terminal blocker

```text
c_*=13u,                                               (2)
```

and every displayed `c_i` is a nonzero multiple of `g mod 13`.  All remaining
terminal labels are transverse to `g mod 13`.  Suppose the terminal danger
bands cover the guard-safe region outside a null set.

We prove two exact statements.

1. In profile

   ```text
   (blockers,aligned,transverse)=(1,6,1),              (3)
   ```

   every triple `I subset {1,...,6}` satisfies

   ```text
   u in Zg+sum_(i in I) Zc_i.                          (4)
   ```

2. In profile

   ```text
   (blockers,aligned,transverse)=(1,7,0),              (5)
   ```

   every five-set `I subset {1,...,7}` satisfies the same generation law

   ```text
   u in Zg+sum_(i in I) Zc_i.                          (6)
   ```

Unlike THM-2139's minimal profile, these conclusions need not put all
characters on one rational line.  They instead retain a uniform hypergraph of
integral span obligations.

## 1. Six aligned labels: every triple captures the blocker

Fix a triple `I subset {1,...,6}`.  We first prove the phase containment

```text
{Y:||g.Y||<1/7 and ||c_i.Y||<1/14 for every i in I}
                  subset {Y:||u.Y||<=1/14}.           (7)
```

If the strict reverse of (7) were nonempty, choose a point `Y` in that open
set with a fully clean `13^2`-root fibre.  This clean choice is legitimate:
on a local inverse domain for `[13]:K->K`, the inverse images of the original
null exceptional set under the 169 root sheets have null union.

Choose guard/root coordinates as in THM-2125.  The blocker is constant and
safe on the fibre.  Since `||g.Y||<1/7`, the guard forbids exactly three of
the thirteen outer columns.  Each selected `c_i` has singleton quotient
phase and therefore forbids exactly one column.  Each of the other three
aligned labels forbids at most two.  The total outer-column capacity is

```text
3+3*1+3*2=12<13.                                     (8)
```

There is a column on which the guard and all six aligned labels are safe.
The sole transverse terminal hits at most two of its thirteen inner points,
so a clean uncovered root remains.  This contradicts the cover and proves
(7).

Let

```text
K_I=ker(g) intersection intersection_(i in I) ker(c_i). (9)
```

Every point of `K_I` lies in the source of (7), so the compact subgroup
`u(K_I)` is contained in the closed radius-`1/14` arc.  No nontrivial compact
subgroup of the circle is contained in such an arc: this includes every
finite cyclic subgroup as well as the full circle.  Hence `u(K_I)={0}`.
Exact annihilator duality gives

```text
K_I^perp=Zg+sum_(i in I)Zc_i,                         (10)
```

which is (4).  Notice that (10) is the integer span, not merely its rational
saturation.

## 2. Seven aligned labels: every five-set captures the blocker

Fix a five-set `I subset {1,...,7}`.  The same clean-root argument proves

```text
{Y:||g.Y||<1/7 and ||c_i.Y||<1/14 for every i in I}
                  subset {Y:||u.Y||<=1/14}.           (11)
```

Indeed, the guard, five singleton phases, and the other two aligned labels
have total outer-column capacity

```text
3+5*1+2*2=12<13.                                     (12)
```

There is no transverse terminal in this profile, so the surviving column
already contains uncovered clean roots.  Restricting (11) to the common
kernel of `g` and the five selected characters, and applying the same exact
annihilator argument, proves (6).

## 3. Determinant ideals and the recursive depth law

Choose an oriented basis of `Gamma` and put

```text
Delta_0=det(g,u),                 Delta_i=det(g,c_i). (13)
```

Taking the determinant with `g` in (4) or (6) gives the ideal inclusions

```text
Delta_0 in (Delta_i:i in I) subset Z.                (14)
```

Equivalently,

```text
gcd(Delta_i:i in I) divides Delta_0                  (15)
```

with the zero ideal interpreted literally.  This formulation also covers a
fully collinear selected set: if every `Delta_i=0`, equation (14) forces
`Delta_0=0`.

Use the convention `nu_13(0)=infinity`.

For six aligned labels, order

```text
a_1<=...<=a_6,             {a_i}={nu_13(Delta_i)}.   (16)
```

The triple of the three deepest determinants has minimum valuation `a_4`.
Applying (15) to it gives

```text
nu_13(Delta_0)>=a_4.                                 (17)
```

For seven aligned labels, order the seven valuations in the same way.  The
five deepest have minimum `a_3`, so

```text
nu_13(Delta_0)>=a_3.                                 (18)
```

Every aligned determinant is divisible by thirteen by hypothesis.  Thus
(17)--(18) imply

```text
13 divides det(g,u).                                  (19)
```

The divided blocker `u`, not only the original blocker `13u`, is therefore
guard-aligned modulo thirteen, possibly with zero reduction.  If its
reduction is nonzero, it is projectively parallel to the guard.  Equivalently,

```text
13^2 divides det(g,c_*).                              (20)
```

This is one additional thirteen-adic alignment digit.  After one division
step, `u` is either a new blocker (`u=0 mod 13`) or an aligned nonblocker.
The theorem does not by itself supply a fresh divided cover or license
iteration of the same fibre argument.

## 4. Scope and Tournament Analysis

The theorem does not turn (17) or (18) into rational collinearity.  A second
step must exploit the new divided alignment, lift the subset obligations to
the next digit, or combine them with the remaining multiple-blocker profiles.

The challenged assumption is that runner labels or the thirteen columns are
the final vertices.  Here the faithful vertices are the **proof obligations**:
the twenty aligned triples in profile (3), or the twenty-one aligned
five-sets in profile (5).  Their observable is the determinant ideal in
(14).  Orienting two obligations by larger minimum 13-adic depth, with label
order as a tie Hamiltonian path, produces a tournament whose score histogram
records search priority.  Its directed cycles, SCCs, edge flips, and path
counts do not determine which character subset generated the ideal.  The
preserved carrier is the uniform obligation hypergraph with its integral-span
and valuation sidecars.  QED.
