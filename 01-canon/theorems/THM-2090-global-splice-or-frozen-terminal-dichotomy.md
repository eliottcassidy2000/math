---
id: THM-2090
title: "Global guard-anchor splice or literally frozen terminal block"
status: >
  PROVED. In the rank-six persistent depth-four branch, the six cut-tree
  relations lift to height-114 support-three relations on the original
  thirteen speeds. Combine them with THM-2052's rank-eleven bounded relation
  code. Either the full row is finite, the last guard and one terminal speed
  anchor a height-91^6 star for all thirteen speeds, or the global template's
  terminal restriction has rank one. In the last branch terminal primitivity
  makes the normalized last-guard/seven-core block literally constant, and
  the only remaining motion is an affine lattice line in the five outer
  coordinates. THM-2092 subsequently bounds that line by the terminal needle,
  making the frozen branch finite. This is a localization, not LRC(14).
source: codex-2026-07-22-LRC14-global-terminal-splice
depends_on:
  - THM-2052
  - THM-2073
  - THM-2088
  - THM-2089
related:
  - THM-2062
  - THM-2065
  - THM-2079
  - THM-2082
  - THM-2086
  - THM-2092
---

# THM-2090 -- global splice or frozen terminal

Retain a hypothetical rank-seven terminal obstruction. Since
`|Q_r|=11-r=7`, its THM-2073 tower has depth exactly four. Write

```text
C=2^4 Q union {h_0,2h_1,4h_2,8h_3},
S=2C union {x,y},
h=h_3,
Q={q_1,...,q_7}.                                      (1)
```

Thus the eight original-row coordinates belonging to the terminal block and
last guard are

```text
H=16h,
V_i=32q_i,                    1<=i<=7.                 (2)
```

The other five coordinates are

```text
2h_0,4h_1,8h_2,x,y.                                   (3)
```

Assume the no-pair, rank-six persistent branch of THM-2088/2089. The theorem
below splices that terminal coefficient system to THM-2052's global bounded-
relation code.

## 1. Six local rows survive in the original coordinates

Choose any spanning tree of the THM-2087 complete cut. Each of its six
terminal relations is

```text
a h+b q_i+c q_j=0,
bc!=0,
max(|a|,|b|,|c|)<=57.                                 (4)
```

In the original coordinates (2), this becomes

```text
2a H+b V_i+c V_j=0.                                   (5)
```

Hence (5) is a support-three relation of height at most `114`. The six rows
remain independent by the leaf-elimination proof of THM-2088. Let `L` be
their six-dimensional span in `Q^13`.

Put

```text
Q_0=91^6=567869252041,
W=W_(Q_0,3)(S),                                       (6)
```

where `W_(Q_0,3)` is THM-2052's span of all support-two/three relations of
height at most `Q_0`. Equations (5)--(6) give

```text
L subset W,
dim W>=11.                                             (7)
```

Since every row in `W` annihilates the positive vector `S`, `dim W<=12`.

## 2. The rank-twelve branch is finite

If `dim W=12`, choose twelve independent defining relations from the bounded
generating packet of `W`. THM-2052's rank-twelve terminal gives

```text
max S<=3^6*(91^6)^12.                                  (8)
```

This is finite and will not be discussed further. Assume henceforth

```text
dim W=11,
U=W^perp,                    dim U=2.                  (9)
```

The six rows of `L` can be extended by five height-`Q_0` support-three rows
to a basis of `W`. Thus only finitely many such full-row coefficient
templates `U` occur.

## 3. Restrict the global template to the terminal block

Let

```text
pi_T:U -> Q^8                                           (10)
```

be restriction to `(H,V_1,...,V_7)`. Because `S in U` and its terminal
coordinates are positive,

```text
1<=rank(pi_T)<=2.                                      (11)
```

The image lies in the two-dimensional kernel of the six tree rows (5).

### Branch A: terminal restriction has rank two

Choose a basis of `U` and regard its thirteen coordinate columns as vectors
in `Q^2`. The column of `H` is nonzero. If `rank(pi_T)=2`, some terminal
speed column `V_j` is independent of it. For every other original coordinate
`S_k`, apply THM-2052's three-coordinate pigeonhole relation to the triple

```text
{H,V_j,S_k}.                                           (12)
```

It gives

```text
A_k H+B_k V_j+C_k S_k=0,
max(|A_k|,|B_k|,|C_k|)<=Q_0.                           (13)
```

The private coefficient `C_k` is nonzero: otherwise (13) would be a relation
between the two independent anchor columns. Therefore **all thirteen speeds**
lie in one two-anchor star whose anchors are specifically the last guard and
a terminal speed.

This is the desired global splice. It retains the last-guard coordinate,
unlike an arbitrary choice of the two THM-2052 anchors.

### Branch B: terminal restriction has rank one

Suppose `rank(pi_T)=1`. Undo the fixed dyadic scaling in (2) by the rational
isomorphism

```text
(H,V_1,...,V_7) ->(H/16,V_1/32,...,V_7/32).            (14)
```

The normalized restrictions of all vectors in `U` lie on one rational line.
The actual normalized terminal vector

```text
w=(h,q_1,...,q_7)                                     (15)
```

is primitive because `gcd(Q)=1`. Any other admissible integer point of the
same depth-four chart has normalized terminal vector `lambda w` for some
rational `lambda`. Integrality and primitivity of `w` force
`lambda in Z`; primitivity of its terminal core forces `|lambda|=1`; positivity
forces

```text
lambda=1.                                              (16)
```

Thus `(h,Q)` is not merely projectively fixed: it is **identical at every
admissible integer point of the template**.

The kernel of `pi_T:U->pi_T(U)` is one-dimensional and vanishes on all eight
terminal coordinates. Consequently the remaining integer points, after
saturation, lie on an affine lattice line

```text
S=S_0+n d,                    n in Z,                  (17)
```

where `d` is supported entirely on the five outer coordinates (3). Positivity,
distinctness, owner addresses, and hereditary primitivity cut (17) further.

Because the coefficient-row universe in (6) is finite, only finitely many
rank-one terminal restriction lines occur. Their primitive normalized
generators form a finite bank of frozen terminal blocks, even though each may
initially appear to carry an unbounded outer line (17). THM-2092 proves that
THM-2077's needle bounds that line absolutely once the block is fixed.

## 4. Exact full-row trichotomy

Every rank-six persistent depth-four obstruction therefore enters exactly one
of the following branches:

```text
I.   rank W=12: full original row lies in the finite box (8);
II.  rank W=11, rank(pi_T)=2:
     all thirteen speeds form a last-guard/terminal-anchor star (13);
III. rank W=11, rank(pi_T)=1:
     (h,Q) is literally frozen and only the five outer coordinates move.  (18)
```

This is the first exact splice of the THM-2087 cut to the original thirteen-
speed relation code. Branch II is ready for the determinant/CRT machinery of
THM-2062 once the coefficient rows are saturated. Branch III is ready for a
one-dimensional outer-line CRT/owner analysis, with terminal containment
precomputed once per frozen block.

THM-2092 subsequently applies the depth-four height transfer
`max(S)<=128 max(Q)/3`. It makes Branch III finite and lifts THM-2088's
bounded-terminal branch to the explicit original-row bound
`3900651012632789`. Consequently Branch II is the sole potentially unbounded
no-pair branch of this pipeline.

Every branch also inherits THM-2086's all-height filters. Writing
`B=max(Q)`, every live terminal satisfies

```text
7 does not divide h,
1<=#{q in Q:7 divides q}<=4,
sum_(q in Q minus {B})q+6h >=(17/1078)B.               (19)
```

In Branch II these are congruence and interval filters on the last-guard/
terminal-anchor lattice. In Branch III they are fixed properties of the one
frozen block and can be precomputed before the outer affine line is examined.

There is an important guardrail concerning THM-2065. The six local rows (5)
have height `114<2^20` and persist on both two-parameter branches. Therefore
THM-2065's existential persistent-circuit alternative is already satisfied
and cannot by itself collapse either branch to finite rays. The finish must
use the **multiplicity and location** of these six circuits, plus phase/owner
sidecars, not ask Fejer to supply one more unspecified relation.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that the terminal two-anchor template and the
global two-anchor template must be coordinated by a large elimination. Their
intersection is controlled by the rank of one restriction map. Rank two
splices the same anchors globally; rank one, together with primitivity, freezes
the terminal block exactly.

Candidate tournament vertices were the thirteen speeds, six local rows, five
outer quotient rows, terminal coordinates, and proof obligations. The faithful
binary split is rank of `pi_T`, which no pairwise orientation preserves. An
orientation by pivot order gives a transitive scheduling tournament with no
cycles or nontrivial SCCs; its Hamiltonian path only records Gaussian-
elimination order. The useful carrier is the restriction matroid of the
bounded relation code, with the dyadic scaling map and terminal gcd as
sidecars. QED.
