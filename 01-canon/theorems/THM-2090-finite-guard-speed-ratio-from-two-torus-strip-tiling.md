---
id: THM-2090
title: "Every rank-seven guard containment has a bounded reduced guard/speed ratio"
status: >
  PROVED by compactness and a connected two-torus strip-tiling obstruction.
  There is an absolute R_7 such that if the closed safe set of seven speeds
  lies in one open 1/7 guard comb, then some speed q has both entries of the
  reduced ratio (h/gcd(h,q),q/gcd(h,q)) at most R_7. The proof is qualitative
  and gives no numerical R_7. It strengthens the support-three relation
  alternative THM-2083 in arity, while THM-2085 remains stronger in giving
  the explicit coefficient height 57 for its three-term relation. This does
  not close the resulting finite ratio templates or LRC(14).
source: codex-2026-07-22-LRC14-two-torus-strip-tiling
depends_on:
  - THM-2080
  - THM-2081
related:
  - THM-2082
  - THM-2083
  - THM-2085
  - THM-2086
  - THM-2087
  - THM-2088
---

# THM-2090 -- a finite guard/speed-ratio alternative

Put

```text
D_q={t in R/Z:||qt||<1/14},
E_h={t in R/Z:||ht||<1/7},
C_h=(R/Z) minus E_h,
G_Q=intersection_(q in Q)((R/Z) minus D_q).             (1)
```

For positive integers `h,q`, define the two-term relation height

```text
rho(h,q)
 =min{max(|a|,|b|):(a,b) in Z^2 minus {(0,0)}, ah+bq=0}
 =max(h/g,q/g),                         g=gcd(h,q).       (2)
```

Thus `rho(h,q)` is exactly the height of the reduced rational ratio `q/h`.

> **Theorem.** There is an absolute positive integer `R_7` such that, for
> every seven-element set `Q` of positive integers,
>
> ```text
> G_Q subset E_h     implies     min_(q in Q)rho(h,q)<=R_7.  (3)
> ```

No parity, divisor-completeness, hereditary-primitivity, or terminal-height
assumption is needed. In the THM-2073 depth-four application those additional
coordinates remain available after (3).

The constant is not made explicit here. The point is structural: a putative
rank-seven obstruction cannot send **all seven** reduced guard/speed ratios
to increasing height. At least one speed belongs to a finite list of rational
multiples of the guard.

## 1. The compact connected-group limit

Suppose, for contradiction, that (3) has no uniform constant. Choose contained
packets

```text
G_(Q_n) subset E_(h_n),
Q_n={q_(n,1),...,q_(n,7)},
min_i rho(h_n,q_(n,i))->infinity.                       (4)
```

Label each set in increasing order and put

```text
v_n=(h_n,q_(n,1),...,q_(n,7)) in Z^8.                  (5)
```

Let `mu_n` be the pushforward of Haar measure on `R/Z` under

```text
t -> v_n t in (R/Z)^8.                                 (6)
```

For a character `m in Z^8`,

```text
mu_n_hat(m)=1 if m.v_n=0, and 0 otherwise.              (7)
```

Enumerate `Z^8` and take a diagonal subsequence on which every value in (7)
is eventually constant. Let

```text
L={m in Z^8:m.v_n=0 eventually},
K=L^perp subset (R/Z)^8.                                (8)
```

The limiting Fourier transform is `1_L`, so

```text
mu_n -> Haar measure m_K weakly.                        (9)
```

The subgroup `L` is saturated: if `k m in L` for a nonzero integer `k`, then
`k(m.v_n)=0` eventually, hence `m.v_n=0` eventually and `m in L`. Therefore
`Z^8/L` is torsion-free and `K` is connected.

Write `chi_0,...,chi_7` for the coordinate characters restricted to `K`.
Condition (4) says that `L` has no nonzero vector supported on `{0,i}`.
By duality, every map

```text
(chi_0,chi_i):K -> (R/Z)^2                              (10)
```

is surjective. In particular its pushforward of `m_K` is two-dimensional
Haar measure.

## 2. Containment becomes an exact seven-strip tiling

In `(R/Z)^8`, let

```text
C_0={x:||x_0||>=1/7},
B_i={x:||x_i||<1/14},
A_i=C_0 intersect B_i.                                  (11)
```

The set

```text
B=C_0 intersect intersection_i((R/Z)^8 minus B_i)       (12)
```

has `mu_n(B)=0` by (4). Its boundary lies in the finite union of coordinate
level sets `||x_0||=1/7` and `||x_i||=1/14`. Every coordinate marginal of
`m_K` is Haar, so this boundary has zero `m_K`-measure. Thus `B` is a
continuity set for (9), and

```text
m_K(B)=0.                                               (13)
```

Consequently the seven `A_i` cover `C_0` almost everywhere. Surjectivity in
(10) gives

```text
m_K(A_i)=(5/7)(1/7)=5/49,
sum_i m_K(A_i)=5/7=m_K(C_0).                            (14)
```

Equality of the sum and union measures forces multiplicity exactly one:

```text
1_(C_0)=sum_(i=1)^7 1_(A_i)       m_K-almost everywhere. (15)
```

In particular,

```text
m_K(C_0 intersect B_i intersect B_j)=0                 (16)
```

for every `i!=j`.

If the three characters `chi_0,chi_i,chi_j` were rationally independent,
their map from `K` onto the three-torus would be surjective. The measure in
(16) would then equal

```text
(5/7)(1/7)(1/7)=5/343>0,                               (17)
```

a contradiction. Hence every such triple is dependent. Fixing `i=1`, the
pair `chi_0,chi_1` is independent and every `chi_j` belongs to their rational
span. The eight coordinate restrictions generate the character group of
`K`, so

```text
rank(K)=2.                                              (18)
```

Thus any counterexample sequence to (3) would converge to an exact tiling of
the guard complement on one connected two-torus by seven centered character
strips.

## 3. Centered character strips cannot make that tiling

Let `Gamma=K_hat`, a free abelian group of rank two. Partition
`chi_1,...,chi_7` into rational rays in `Gamma tensor Q`, identifying opposite
rays because `||-x||=||x||`. For a ray `R`, choose its primitive character
`psi_R` and write

```text
chi_i=m_i psi_R,                    m_i in Z minus {0},
S_R(z)=sum_(i in R)1_{||m_i z||<1/14}.                  (19)
```

The integer-valued circle step function `S_R` is not constant. It equals
`|R|` near `z=0`, whereas multiplication by each nonzero `m_i` preserves Haar
measure, so

```text
integral S_R(z) dz=|R|/7.                              (20)
```

Therefore `S_R` has a boundary value `z_R` across which its net jump is
nonzero.

The characters `chi_0` and `psi_R` are rationally independent by (10), so

```text
(chi_0,psi_R):K -> (R/Z)^2                             (21)
```

is surjective. The boundary circle `psi_R=z_R` consequently meets the open
guard complement `||chi_0||>1/7` in open arcs. Boundary circles belonging to
other rational rays meet it in only finitely many points; the two guard
boundary circles do likewise. Choose a point on one of those open arcs away
from every such intersection.

In a small neighbourhood of that point, crossing `psi_R=z_R` changes `S_R`
by its nonzero jump, while every other ray contribution is locally constant
and both sides remain in the guard complement. The total multiplicity

```text
sum_R S_R(psi_R(x))=sum_i 1_(B_i)(x)                   (22)
```

therefore takes different constant integer values on the two adjacent open
cells. Equation (15) says it equals one almost everywhere on both cells, a
contradiction.

This rules out the sequence (4). The usual bad-sequence contradiction now
supplies one absolute `R_7`, proving (3). QED.

## 4. Frontier effect and scope

THM-2083 proves that some guard/two-speed triple has a uniformly bounded
support-three relation, and THM-2085 gives the explicit height `57` for that
statement. The present theorem exchanges effectivity for lower arity: some
single reduced ratio

```text
q/h=b/a,                    1<=a,b<=R_7                 (23)
```

comes from a finite list. This is exactly the ratio coordinate exposed by
THM-2082's translated-prime-grid branch, so the two results can now be joined
without pretending that rank-one code support remembers projective residues.

The theorem does not bound the common scale, select which of the seven speeds
has ratio (23), or eliminate the other six speeds. It therefore does not close
the depth-four terminal or LRC(14). A useful effective continuation is to
combine the degree-57 Selberg certificate with the boundary-ray proof to put
a numerical value on `R_7`, then run a symbolic ratio-template argument rather
than a height box.

The challenged assumption is that the relevant vertices are runners. In the
limit proof the faithful vertices are **rational character rays**; their
boundary circles preserve the exact outside-guard tiling predicate. Pairwise
outside-guard intersection is the observable, while rational proportionality
is the gauge. Orienting these rays as a tournament would discard the jump
sign that proves the contradiction. The faithful carrier is instead the
rank-two relation matroid together with its signed boundary arrangement.
