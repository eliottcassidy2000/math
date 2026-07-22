---
id: THM-2123
title: "Dependent 13-blocker closure and the complete transverse rank-eight invoice"
status: >
  PROVED. In a rank-two guard/eight-terminal torus cover, assume the guard is
  nonzero modulo 13, exactly one terminal is divisible by 13, that blocker's
  divided character is rationally dependent with the guard, and all seven
  residual terminals are transverse to the guard modulo 13. An almost-
  everywhere cover is impossible. THM-2120's pair-phase containment remains
  valid without guard--blocker independence. Character rigidity leaves either
  three terminal bands of total measure at most 3/7<5/7, or a rank-one span,
  contradicting the ambient rank-two hypothesis. With THM-2114, THM-2120,
  and THM-2122, every rank-eight cover must therefore have either a guard
  13-blocker or a nonblocker terminal projectively parallel to the guard mod 13.
source: codex-2026-07-22-LRC-dependent-thirteen-blocker
depends_on:
  - THM-2114
  - THM-2116
  - THM-2120
  - THM-2122
related:
  - THM-2098
  - THM-2124
---

# THM-2123 -- dependent 13-blocker closure

Let `Gamma` be a rank-two character lattice, let

```text
K=Hom(Gamma,R/Z),                                      (1)
```

and let the nonzero characters

```text
g,c_*,c_1,...,c_7 in Gamma                            (2)
```

span `Gamma tensor Q`. Assume

```text
c_*=13u,                           u!=0;
g mod 13 is nonzero;
g and u are Q-dependent;
c_i mod 13 is not proportional to g mod 13,  i=1,...,7. (3)
```

Suppose, toward a contradiction, that outside a null set `E` the terminal
danger bands cover the guard-safe region:

```text
{X:||g.X||>1/7}
 subset {X:||c_*.X||<1/14} union union_i {X:||c_i.X||<1/14}.   (4)
```

Put `epsilon=1/14`.

## 1. Pair-phase containment does not use independence

We first recover the load-bearing conclusion of THM-2120 under the weaker
hypotheses (3). Fix distinct residual labels `i,j` and suppose that

```text
O_ij={Y:||c_i.Y||<epsilon,
          ||c_j.Y||<epsilon,
          ||u.Y||>epsilon}                            (5)
```

is nonempty. It is open.

Multiplication by thirteen is a surjective finite covering `K->K`. For every
`Y`, its 169 roots differ by `K[13]`. Since `g mod 13` is nonzero, the guard
values across those roots contain a full translate of the thirteen-grid.
Some root is strictly guard-safe: a grid point lies within `1/26` of `1/2`,
and hence has norm at least `6/13>1/7`.

Choose such a root above a point of `O_ij`. Strictness and a local inverse
sheet give a nonempty open set of roots `P_ij` on which

```text
13X in O_ij,                         ||g.X||>1/7.      (6)
```

Let `v` span the kernel of `g` on `K[13]`. Discard the null set

```text
E_v=union_(r in F_13)(E-rv).                          (7)
```

Choose `X in P_ij\E_v`. On the thirteen-point orbit `X+rv`, the guard and
blocker are constant, and

```text
c_*.X=u.(13X),                     13(c_k.X)=c_k.(13X). (8)
```

Thus the blocker is safe. THM-2116's exact grid count says the strict
inequalities for `c_i.(13X)` and `c_j.(13X)` in (5) make their two orbit
danger sets singletons. Each other residual danger set has size at most two
by transversality. Total incidence is at most

```text
1+1+5*2=12<13,                                         (9)
```

contradicting the cover of all thirteen orbit points outside (7).

Therefore `O_ij` is empty for every pair, or equivalently

```text
{Y:||c_i.Y||<epsilon and ||c_j.Y||<epsilon}
 subset {Y:||u.Y||<=epsilon}.                         (10)
```

No step above used Q-independence of `g,u`. That hypothesis entered THM-2116
only to exhibit an easy safe base point, and entered THM-2120 only in its final
disposal of the all-line ledger.

## 2. Character rigidity leaves two ledgers

Apply THM-2120's exact two-character rigidity lemma to (10). For every
residual pair `c_i,c_j`:

```text
if c_i,c_j are independent, then u is +/-c_i or +/-c_j;
if c_i,c_j are dependent, then c_i,c_j,u are proportional.   (11)
```

If every residual is proportional to `u`, call this the all-line case.
Otherwise choose `c_1` off the `u`-line. For each `j>=2`, the pair `c_1,c_j`
cannot be dependent, since the second clause of (11) would put `c_1` on the
`u`-line. The first clause then forces

```text
c_j=+u or c_j=-u.                                     (12)
```

Hence exactly the same dichotomy as THM-2120 holds:

```text
(A) all seven residuals lie in Q u;
(B) exactly one lies off Q u and the other six equal +/-u.   (13)
```

## 3. Dependence makes both endings simpler

In case (B), the eight terminal danger bands reduce to at most three distinct
bands:

```text
the band for 13u,       the common band for +/-u,
and the exceptional residual band.                    (14)
```

Every nonzero character pulls back the radius-`1/14` arc to a set of Haar
measure `1/7`. Therefore the full terminal union has measure at most `3/7`.
The guard-safe set has measure

```text
1-2/7=5/7.                                             (15)
```

It cannot be covered even almost everywhere by (14).

In case (A), all terminals lie on the rational line `Q u`. The present
dependent hypothesis also places `g` on `Q u`. Thus all nine characters in
(2) span a space of rational dimension at most one, contradicting the stated
rank-two span.

Both cases are impossible, proving that (4) cannot hold. QED.

## 4. Complete transverse mod-13 invoice at rank eight

There is a coordinate-free corollary. In any rank-two rank-eight cover with
the span hypotheses above, THM-2114 says that at least one of the guard and
eight terminal characters vanishes modulo 13. Suppose the guard does not
vanish. Then at least one terminal is a 13-blocker.

If there are at least two terminal blockers, THM-2122 excludes the cover when
all remaining terminals are transverse to the guard. If there is exactly one,
write it `13u`: THM-2120 excludes `g,u` independent, while the theorem above
excludes `g,u` dependent. Consequently every such cover satisfies

```text
g=0 mod 13
or
there exists i with c_i!=0 mod 13 and
                    c_i in F_13^* g mod 13.            (16)
```

In determinant/content language, (16) is

```text
13|cont(g)
or
there exists a non-13-blocker c_i with 13|det(g,c_i).  (17)
```

This strictly sharpens the content invoice: after paying a terminal
13-blocker, a putative cover must also pay projective guard alignment with a
different terminal. The two labels cannot be conflated because the aligned
terminal is required to be nonzero modulo thirteen.

THM-2125 later strengthens the second branch to at least five nonblocker
terminals on the guard's projective mod-13 line. THM-2124 sharpens the guard-
blocker branch, when all terminals are nonblockers, to direction partition
`(8)` or `(7,1)`; under THM-2097's specialization hypotheses the latter has
seven rational guard-proportionals and one transverse terminal. Ranks nine
through twelve are not closed here.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption was that guard--blocker dependence invalidates the
global phase-rigidity mechanism. It invalidates only THM-2120's final
independent-coordinate choice. The thirteen-root guard sheet, two-singleton
capacity contradiction, and character-kernel lemma are unchanged; dependence
then turns the last all-line ledger into an immediate rank contradiction.

Candidate tournament vertices were character rays, terminal labels, residual
pairs, mod-13 projective points, and proof obligations. Orienting residual
pairs by which member of (11) captures `u`, with label order as the tie
Hamiltonian path, preserves the dichotomy search but destroys the common
character line and the rank-two span. The faithful carrier is the projective
mod-13 incidence `(g;c_*,c_i)` together with integral character kernels and
the ambient rational rank. Scores, cycles, SCCs, and edge flips do not recover
the determinant invoice (17).
