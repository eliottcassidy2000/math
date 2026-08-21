---
id: THM-3578
title: "Zariski-main boundary rank and the nonproper sheet-debt inequality"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  Let f:A2_C->Y be a dominant
  quasi-finite morphism of generic degree d to a normal affine surface.  In
  the finite normalization supplied
  by Zariski Main, the divisorial boundary has at least
  rank(Cl(Y)/torsion) irreducible components.  If f is etale and its
  divisorial nonproperness curves have generic affine fibre sizes n_C, then
  sum_C(d-n_C) is at least that rank.  The Chebyshev--Pell and three-arm S5
  Danielewski maps have exact sheet-debt vectors (d-1,d-1,1) and
  (4,2,2,3,3), respectively.  This is a necessary nonproperness-capacity
  law, not a planar Jacobian obstruction when the target is A2.
source: kps-s188, 2026-08-20
audit: >
  An independent reconstruction verified the localization/conorm rank bound,
  the generic-DVR sheet-length formula without a global flatness assumption,
  the boundary-to-Jelonek incidence, and both exact fibre/debt vectors.
depends_on:
  - THM-3404-factorized-danielewski-principal-parts-and-finite-cover-obstruction
  - THM-3566-chebyshev-pell-odd-keller-collision-tower
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
related:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3576-higher-exponent-belyi-keller-collision-tower
  - THM-3574-universal-reducible-target-graph-component-unit-no-go
---

# THM-3578 -- Zariski-main boundary rank and sheet debt

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem quantifies the
exact resource that separates the new Danielewski near-counterexamples from
a finite cover of affine space.
It does not prove `JC(2)`, and for target `A2` its class-group lower bound is
zero.

All varieties are over `C`.  The argument works over any algebraically
closed field after replacing geometric fibre cardinalities by separable
degrees and lengths.

## 1. Boundary-rank theorem

Let `Y` be a normal affine surface and let

```text
f:A2 -> Y                                                   (1)
```

be dominant and quasi-finite of generic degree `d`.  Let `K=C(A2)`, and let
`Xbar` be the normalization of `Y` in `K`.  Zariski's Main Theorem factors
`f` as

```text
A2 --j--> Xbar --pi--> Y,                               (2)
```

where `j` is an open immersion and `pi` is finite dominant of degree `d`.
Write

```text
B=Xbar\j(A2),
D_1,...,D_s = the codimension-one irreducible components of B. (3)
```

Then

```text
rank(Cl(Y)/Cl(Y)_tors) <= s.                              (4)
```

### Proof

The finite normalization `Xbar` is normal.  The localization sequence for
Weil divisor class groups is right exact:

```text
direct_sum_(i=1)^s Z[D_i] -> Cl(Xbar) -> Cl(A2) -> 0. (5)
```

Since `Cl(A2)=0`, the boundary classes generate `Cl(Xbar)`, and therefore

```text
rank Cl(Xbar) <= s.                                    (6)
```

Finite conorm and pushforward give

```text
pi_* pi^* = d : Cl(Y) -> Cl(Y).                          (7)
```

This is the valuation/norm identity proved in THM-3404, Section 6.2; it
does not require flatness.  Hence the kernel of
`pi^*:Cl(Y)->Cl(Xbar)` is killed by `d`.  Conorm is therefore injective on
the torsion-free quotient after tensoring with `Q`, so

```text
rank(Cl(Y)/torsion) <= rank Cl(Xbar).                  (8)
```

Combining `(6)` and `(8)` proves `(4)`.  Notice that isolated points of
`Xbar\A2` cannot pay this bill: deleting codimension two from a normal
surface does not change its class group.  Only divisorial escape counts.

## 2. The sheet-debt inequality

Assume in addition that `f` is etale.  Let `J(f)` be its nonproperness
locus, and let `C` range over the irreducible curve components of `J(f)`.
At a general complex point of `C`, let

```text
n_C = # f^(-1)(y).                                      (9)
```

Etaleness makes every affine-source point contribute length one.  The
finite map `pi` has total fibre length `d`, so define the missing-sheet debt

```text
delta_C=d-n_C.                                         (10)
```

Then `delta_C>0`, every `D_i` maps finitely onto one such `C`, and

```text
s <= sum_C delta_C.                                    (11)
```

Indeed, over the generic point of `C`, the boundary part of the finite
fibre has length `delta_C`.  Each boundary prime above `C` contributes a
positive ramification-index times residue-degree, hence at least one.  The
finite morphism `pi` preserves dimension, so a boundary curve cannot map to
a point.  Conversely, the finite-envelope criterion for properness says
that every divisorial component of `J(f)` receives at least one boundary
prime.  Summing over `C` proves `(11)`.

Together with `(4)`, this gives the nonproperness-capacity law

```text
rank(Cl(Y)/torsion)
   <= # {divisorial boundary primes in the finite envelope}
   <= sum_(C subset J(f)) (d-n_C).                     (12)
```

The middle term is deliberately retained.  Distinct boundary primes may
map to the same target curve, so class-group rank does not by itself lower
bound the number of components of `J(f)`.

## 3. Exact invoices for the two live Danielewski towers

### 3.1 Chebyshev--Pell odd tower

For THM-3566, `Y_beta:c^2e=b(b-beta)` has

```text
Cl(Y_beta)=Pic(Y_beta)=Z,                               (13)
```

and the map has odd generic degree `d`.  Its three nonproperness lines and
generic affine fibre sizes are

| curve | generic fibre size | debt |
|---|---:|---:|
| `C_0={b=0,e=0}` | `1` | `d-1` |
| `C_beta={b=beta,e=0}` | `1` | `d-1` |
| `L_beta={b=beta,c=0}` | `d-1` | `1` |

Thus

```text
debt vector=(d-1,d-1,1),       total=2d-1 >= 1.        (14)
```

The omitted point `(beta,0,0)` is not the whole boundary mechanism.  The
two additional target lines carry the sheets that escape while other sheets
still hit the same target points.

### 3.2 Symmetric three-arm S5 carrier

For THM-3572,

```text
Y_3:c^2e=b(3125b^2+256),
Cl(Y_3)=Pic(Y_3)=Z^2.                                  (15)
```

The degree-five map has five nonproperness lines.  In the order consisting
of the central `e=0` line, the two side `e=0` lines, and the two side `c=0`
lines, its generic affine fibre sizes and debts are

```text
n=(1,3,3,2,2),             delta=(4,2,2,3,3).          (16)
```

Hence

```text
sum delta_C=14 >= 2.                                   (17)
```

The five visible target curves force at least five finite-envelope boundary
primes, while `(4)` independently forces at least two.  The two lower bounds
measure different information: incidence with the Jelonek set versus the
free divisor-class lattice.

### 3.3 Higher-exponent Belyi tower

Incoming THM-3576 supplies, for every `n>=2`, an etale map of degree

```text
d_n=n(n-1)+1
```

to `Y_n:c^n e=b(b-beta_n)`, again with `Pic(Y_n)=Z`.  In the order central
`e=0`, side `e=0`, and side `c=0`, its generic affine fibre sizes are

```text
(1,(n-1)^2,n),
```

so its exact sheet-debt vector is

```text
(n(n-1), n, (n-1)^2),             total=2d_n-1.       (17a)
```

The total equals the global Euler defect

```text
d_n chi_c(Y_n)-chi_c(A2)=2d_n-1,                      (17b)
```

because `chi_c(Y_n)=2`.  Equality between the unweighted generic-curve debt
sum and the Euler defect is special to this three-arm passport; `(12)` does
not assert it for an arbitrary quasi-finite map.

## 4. Sharp controls and scope

1. **Finite cover.**  If `f` were finite, then `B` would be empty.  Formula
   `(4)` would force `Cl(Y)/torsion=0`, recovering the free-rank part of
   THM-3404's finite-factorial-cover obstruction.
2. **Codimension-two omission.**  A set-theoretic image can omit only points
   while the finite envelope still has divisorial boundary.  THM-3566 and
   THM-3572 are exact examples.  Image complement and nonproperness locus
   must not be conflated.
3. **Torsion.**  A degree-`d` finite cover can kill `d`-torsion.  Therefore
   `(4)` uses `Cl/torsion`; replacing it by the minimal number of generators
   of all of `Cl` is false.
4. **Normality.**  Without normality downstairs, integral elements in the
   shared fraction field can become regular upstairs, and conorm need not
   give the stated class-group injection.
5. **Target `A2`.**  Here `Cl(A2)=0`, so `(12)` is vacuous.  It is a design
   law for intermediate non-`A2` completions, not a proof of the planar
   Jacobian conjecture.

## 5. Concrete next target

For a putative Darboux pair `G:Y_Sigma->A2`, composition with one of the
etale maps `Phi:A2->Y_Sigma` would be a Keller counterexample.  The theorem
shows what must happen before composition: because `Cl(A2)=0`, `G` must
destroy the `h-1` arm classes through its own nonfinite boundary geometry.
The next useful invariant is therefore a **two-stage boundary matrix**:

```text
source boundary primes of Phi
 -> arm lattice Cl(Y_Sigma)
 -> boundary valuations of G.                           (18)
```

A scalar Euler characteristic or the identity `1={P,Q}` cannot see this
matrix.  The cheapest decisive experiment is to compute it for each
surviving `(2,4)` and `(3,3)` weight-support cell before attempting a full
coefficient elimination.
