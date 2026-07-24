---
id: THM-2148
title: "Three-target finite-kernel Fano obstruction"
status: >
  PROVED. If three radius-1/14 character bands cover the asymmetric source
  rectangle cut out by radii 1/7 and 1/14, then one target lies in the exact
  integer span of the two source characters (or on their saturated rational
  line in the dependent case). The only possible obstruction on the common
  finite kernel is the three-colour Klein-four cover. Translated-fibre
  coverage kills it: the three doubled target classes would need all nonzero
  parities of F_2^2, while the source polar diamond omits parity (1,1).
  Applied to the rank-eight profile (3,5,0), this gives a five-edge
  blocker-depth incidence invoice. It does not yet eliminate that profile.
source: codex-2026-07-24-LRC-three-target-fano-kernel
depends_on:
  - THM-2125
related:
  - THM-1166
  - THM-2139
  - THM-2140
  - THM-2141
script: 04-computation/lrc14_three_target_fano_kernel_referee_thm2148.py
output: 05-knowledge/results/lrc14_three_target_fano_kernel_referee_thm2148.out
script_sha256: 67ca620a3dedd75222c0b2a9958b9c28594c4421276f7e1b7f6a13e36c756a44
output_sha256: f92309860bf01db92e53bb01095ecd0e86e46eab8d656b041aad3b890f5a50b3
---

# THM-2148 -- three-target finite-kernel Fano obstruction

**PROVED.** The all-group proof below and its LRC application received an
independent hostile audit, including a separate endpoint sweep through
character-image order `99,999`.  The finite script is a control rather than
a dependency.

Let `Gamma` be a rank-two character lattice and

```text
K=Hom(Gamma,R/Z).
```

For a character `w`, put

```text
D_w={X in K:||w.X||<=1/14}.                           (1)
```

Take nonzero characters `a,b,w_1,w_2,w_3` and suppose

```text
{X:||a.X||<1/7 and ||b.X||<1/14}
                         subset D_(w_1) union D_(w_2) union D_(w_3). (2)
```

The theorem has two clauses.

1. If `a,b` are rationally independent, then for some `j`

   ```text
   w_j in Za+Zb.                                      (3)
   ```

2. If `a,b` are rationally dependent, then some `w_j` lies on their common
   saturated rational line.

The point of (3) is exact integral descent, not rational membership.  A
three-colour cover of the common finite kernel does exist abstractly: it is
the three-line cover of a Klein four-group.  The full source rectangle,
including its translated fibres, forbids exactly that exception.

## 1. Finite three-colour lemma

We first isolate the finite statement.

> **Finite lemma.** Let `F` be a finite abelian group and let
> `chi_1,chi_2,chi_3:F->R/Z` be characters.  If
>
> ```text
> F=union_(j=1)^3 {x:||chi_j(x)||<=1/14},              (4)
> ```
>
> then either one `chi_j` is trivial, or the three characters are exactly
> the three nonzero characters of a quotient `F->C_2 x C_2`.

### Proof

If a nontrivial character has image order `m`, its centered danger set has
relative size

```text
d(m)=[2 floor(m/14)+1]/m.                             (5)
```

For `m>=3`, one has `d(m)<=1/3`; equality is possible only at `m=3`.
If all three characters were nontrivial and none had order two, the sum of
their danger-set cardinalities would be at most `|F|`.  The identity belongs
to all three sets, so it is overcounted at least twice.  Their union would
have cardinality strictly below `|F|`, contradicting (4).  Thus, unless a
character is trivial, relabel so that `chi_1` has order two.

Put

```text
H=ker(chi_1).
```

The first danger set is exactly `H`, so the other two sets cover the other
coset `C=F\H`.  A translate of a radius-`1/14` arc contains at most

```text
floor(m/7)+1<=m/2                                    (6)
```

points of an `m`-grid; equality in the resulting density bound is possible
only for `m=2`.  If `chi_j|H` were trivial, then `chi_j` would factor through
`F/H=C_2`.  A nontrivial such character is identically `1/2` on `C`, so its
centered danger set misses `C`.  Hence both restrictions to `H` are
nontrivial.

Each of their danger sets occupies at most half of `C`.  Since the two sets
cover `C`, both occupy exactly half and they are disjoint there.  Equality
in (6) makes both restrictions order two.  Cosets of two distinct index-two
subgroups of `H` intersect, so disjointness forces

```text
chi_2|H=chi_3|H=:eta.                                 (7)
```

Let `J=ker(eta)`.  Then `F/J` has order four, and both `chi_2,chi_3` factor
through it.  This quotient cannot be cyclic: on `C/J`, a character of
`C_4` that is nontrivial on its order-two subgroup takes the values
`1/4,3/4`, both outside the radius-`1/14` arc.  Therefore

```text
F/J is isomorphic to C_2 x C_2.                       (8)
```

On this quotient, (7) and disjointness on `C/J` say that `chi_2,chi_3` are
the two distinct extensions of `eta`.  Together with `chi_1`, they are the
three nonzero characters.  This proves the finite lemma.  QED.

The endpoint convention is included in (5)--(6).  In particular, the
closed arc can contain both endpoints of an `m`-grid interval, but equality
with one half of the grid still occurs only at `m=2`.

## 2. The translated rectangle kills the Klein four-group

Assume now that `a,b` are rationally independent.  The map

```text
pi=(a,b):K->(R/Z)^2                                  (9)
```

is surjective and has finite kernel

```text
L=ker(a) intersection ker(b).                         (10)
```

At source phase `(0,0)`, containment (2) restricts to a cover of `L`.
If some `w_j|L` is trivial, exact annihilator duality gives

```text
L^perp=Za+Zb,
```

and proves (3).  Suppose instead that every restriction is nontrivial.  The
finite lemma says that they are the three nonzero order-two characters of
one Klein-four quotient of `L`.

Fix any phase

```text
(x,y),           ||x||<1/7, ||y||<1/14.               (11)
```

On its `L`-fibre, a target with nontrivial order-two restriction occupies
either no points or one coset of its index-two kernel.  Any two of the three
distinct Klein-four directions cover at most three quarters of the fibre.
Thus all three targets must be active on every source fibre.

Since `2w_j` kills `L`, write exactly

```text
2w_j=m_j a+n_j b,                   m_j,n_j in Z.      (12)
```

Activity of the `j`-th band on the fibre over `(x,y)` implies

```text
||m_j x+n_j y||<=1/7.                                (13)
```

This holds throughout the open rectangle (11).  Along a radial segment from
the origin the lifted value in (13) cannot jump across the gap between
different integer translates.  Taking the two coordinates to the signed
boundary therefore gives the polar inequality

```text
|m_j|/7+|n_j|/14<=1/7,
2|m_j|+|n_j|<=2.                                     (14)
```

There is now a parity obstruction.  Restriction gives an injective map

```text
(Gamma/(Za+Zb))[2] -> (Za+Zb)/2(Za+Zb),
[w] |-> 2w mod 2(Za+Zb).                              (15)
```

Injectivity uses that `Gamma` is torsion-free.  The three nonzero classes of
the Klein four-group must consequently have coefficient parities

```text
(1,0), (0,1), (1,1) in F_2^2.                         (16)
```

But the nonzero integer points of the polar diamond (14) are

```text
(m,n)=(+/-1,0),(0,+/-1),(0,+/-2).                     (17)
```

Their parities are only `(1,0),(0,1),(0,0)`.  The required `(1,1)` tooth is
missing.  This contradiction proves (3).

This is the exact loss ledger for the Fano/Klein-four view:

```text
source:       three target bands on the common finite kernel;
quotient:     the three nonzero colours of (C_2)^2;
preserved:    coverage of the zero fibre;
destroyed:    activity on translated source fibres;
sidecar:      doubled coefficient vectors in the source polar diamond;
decisive test: whether all three nonzero parity classes occur.
```

The zero-fibre quotient passes.  Its translated-fibre sidecar fails at the
third parity class.

## 3. The dependent case

If `a,b` are rationally dependent, let `alpha` generate their common
saturated rational line.  The connected circle `ker(alpha)` lies in the
source of (2).  A target not on the `alpha`-line restricts to a nontrivial
circle character, whose danger band has Haar measure `1/7`.  Three such
bands have union measure at most `3/7<1`.  Hence one target kills
`ker(alpha)` and lies in `Z alpha`.  This proves the second clause.

## 4. Application to the remaining three-blocker pencil

Return to the rank-eight profile

```text
(blockers,aligned,transverse)=(3,5,0).                (18)
```

Write the three true blockers and five aligned nonblockers as

```text
c_*j=13u_j,                  j=1,2,3,
c_1,...,c_5 in Gamma,                                      (19)
```

where every `c_i mod 13` is a nonzero multiple of `g mod 13`.  Suppose the
eight terminal danger bands cover the guard-safe region outside a null set.

For every aligned label `i`, the clean `13^2`-root ledger gives

```text
{Y:||g.Y||<1/7, ||c_i.Y||<1/14}
                    subset union_(j=1)^3 {Y:||u_j.Y||<=1/14}. (20)
```

Indeed, if the strict reverse were nonempty, choose a phase with all 169
root sheets outside the original null exceptional set.  The three blockers
are constant and safe on those roots.  The guard forbids three outer
columns, the selected aligned label one column, and the other four aligned
labels at most two each.  The total capacity is

```text
3+1+4*2=12<13.                                       (21)
```

A surviving column consists of clean roots safe for the guard, all five
aligned terminals, and all three blockers, contradicting the cover.

Choose an oriented basis and put

```text
Delta_i=det(g,c_i),             Delta_*j=det(g,u_j).  (22)
```

Applying the theorem to (20) gives, for every `i`, some `j` such that

```text
Delta_i divides Delta_*j.                             (23)
```

When `Delta_i=0`, the dependent clause says literally `Delta_*j=0`, so (23)
uses the zero ideal without saturation.  Otherwise (23) follows by taking
`det(g,-)` in `u_j in Zg+Zc_i`.

Thus the five aligned labels admit edges to three divided blockers, with
each edge carrying an exact determinant-ideal containment.  In particular,
with `nu_13(0)=infinity`,

```text
max_j nu_13(Delta_*j) >= max_i nu_13(Delta_i) >=1.    (24)
```

For the blocker selected by a deepest aligned determinant,

```text
nu_13(det(g,c_*j))
   =1+nu_13(Delta_*j)
   >=1+max_i nu_13(Delta_i).                          (25)
```

Also, after choosing one valid edge for each aligned label, one blocker has
at least two incident labels.  If those labels are independent of `g`, its
determinant is divisible by the least common multiple of their two aligned
determinants.

Equations (23)--(25) are the new blocker-depth invoice.  They remove the
finite-kernel Fano exception, but do not empty profile (18).  The remaining
possibilities include a genuine three-band cover after all three targets
descend to one source lattice, and line-rich scalar profiles with three deep
coefficients.  A height bound, a second translated-polar inequality, or a
new scalar `4+3`/`5+3` annulus law is still required.  No iteration of the
root argument is asserted.

## 5. Carrier and Tournament Analysis

The natural vertices are not the eight terminals alone.  They are the five
source rectangles and the three divided blockers, with a bipartite edge

```text
i--j iff u_j in Zg+Zc_i.                              (26)
```

The observable on an edge is the determinant ideal in (23); its
Archimedean and translated-fibre data must remain sidecars.  Forcing this
relation into a tournament would introduce cosmetic ties and forget which
source lattice annihilates which finite kernel.  The faithful carrier is
the labelled bipartite incidence graph, not a tournament.

The historical Fano/`chi_7` lens is useful here only as an obstruction
detector: three colours cover a Klein-four zero fibre, but the continuous
polar sidecar rejects the third colour.  It supplies no symmetry quotient
of the surviving determinant-incidence graph.

## 6. Exact hostile controls

The finite proof above is uniform and does not depend on enumeration.  The
matching exact referee nevertheless exhausts every invariant-factor group

```text
C_m x C_n,             m|n, n<=18,                    (27)
```

and every triple of nontrivial characters.  Across `58` groups it finds no
two-band cover and exactly `23` three-band covers.  Every cover is the
unique triple of nonzero order-two characters on a group with two-dimensional
two-torsion.  It also independently enumerates (17) and confirms that
`(1,1)` is the sole required Klein-four parity absent from the polar diamond.

Reproduce with

```text
python3 04-computation/lrc14_three_target_fano_kernel_referee_thm2148.py
```

The script prints the consequence objects--covering triples and polar
parities--rather than assumptions injected into the model.  The exhaustive
finite universe is a hostile endpoint/torsion control, not a bounded
substitute for the all-group proof.
