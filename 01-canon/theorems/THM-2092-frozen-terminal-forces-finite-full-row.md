---
id: THM-2092
title: "A frozen or bounded depth-four terminal block forces a finite full row"
status: >
  PROVED REDUCTION. THM-2077 gives max(S)<=128 max(Q)/3 at depth four.
  Therefore THM-2090's literally frozen terminal branch is finite, and
  THM-2088's rank-seven cut branch gives the explicit original-row bound
  max(S)<=3900651012632789. Consequently the only potentially unbounded
  no-pair branch left by THM-2087--2090 is the rank-eleven global template in
  which the last guard and one terminal speed anchor all thirteen speeds.
  The resulting finite banks have not been exhaustively decided, so this is
  not LRC(14).
source: codex-2026-07-22-LRC14-frozen-terminal-finite
depends_on:
  - THM-2077
  - THM-2088
  - THM-2090
related:
  - THM-2078
  - THM-2086
---

# THM-2092 -- bounded terminal means bounded full row

Retain the depth-four rank-seven tower

```text
S=2C union {x,y},
C=2^4 Q union {h_0,2h_1,4h_2,8h_3},
B=max(Q).                                               (1)
```

## 1. Exact terminal-to-full-row height transfer

THM-2077 proves at arbitrary depth `r` that

```text
max(S)<=2^(r+1)B*max(1,(12-r)/(r+2)).                  (2)
```

At `r=4`, this becomes

```text
max(S)<=32B*(8/6)=128B/3.                              (3)
```

Since `max(S)` is an integer,

```text
max(S)<=floor(128B/3).                                 (4)
```

This is a weak bound only because the lifted terminal speed `32B` may attain
its endpoint; the earlier-guard and tail inequalities used in THM-2077 are
strict. No endpoint convention is changed here.

## 2. THM-2090's frozen branch is finite

In THM-2090 Branch III, the normalized terminal block

```text
(h_3,Q)                                                 (5)
```

is literally identical at every admissible integer point of a fixed global
coefficient template. In particular its `B` is fixed. Equation (4) puts every
one of the thirteen original speeds in the finite box

```text
{1,...,floor(128B/3)}.                                 (6)
```

Thus the affine lattice line moving the three earlier guards and two original
tails has only finitely many positive distinct points satisfying the tower.

THM-2090 obtains its templates from a finite bounded-relation packet, so only
finitely many frozen primitive terminal blocks occur across the entire Branch
III ledger. Hence their union is finite, not merely finite template by
template. The bound is effective by enumerating the coefficient matrices and
taking their primitive terminal restriction generators; no such enumeration
is claimed here.

## 3. THM-2088's finite terminal lifts uniformly

In THM-2088's rank-seven cut-matrix branch,

```text
B<=91421508108581.                                     (7)
```

Combining (4) and (7) gives the explicit full-row bound

```text
max(S)
 <=floor(128*91421508108581/3)
 =3900651012632789.                                    (8)
```

Thus that entire branch is a finite set of original thirteen-speed rows,
including all earlier guards and original tails. This is considerably smaller
than the general THM-763 finite-height box, but it is still far beyond the
current exact census.

## 4. Revised unbounded frontier

Combine THM-2087--2090 with Sections 2--3. Apart from the bounded two-term
guard-ratio branch, every no-pair rank-seven obstruction is now either

```text
finite original row,
or
rank(W)=11 and rank(pi_T)=2,                            (9)
```

where the second branch is THM-2090's globally spliced star

```text
A_k(16h)+B_k(32q_j)+C_k S_k=0,
C_k!=0,
max(|A_k|,|B_k|,|C_k|)<=91^6                          (10)
```

for every original speed coordinate `S_k` outside the two anchors.

Therefore (10) is the sole potentially unbounded **no-pair** family produced
by the cut/holonomy/global-relation pipeline. It still carries six persistent
height-114 cut rows, THM-2086's modular/nonlacunary filters, THM-2079's owner
address, and THM-2077's relative height box. Classifying (10), plus the
separate bounded guard-ratio families, is the next decisive target.

Finite means mathematically reducible to exact computation, not already
computed. Neither (6) nor (8) proves those rows safe.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that an affine outer line remains an unbounded
problem after the terminal block freezes. The terminal needle controls every
outer speed by the same fixed `B`, so the line is cut by an absolute box.

Possible tournament vertices were terminal blocks, outer coordinates, height
inequalities, dyadic levels, and proof obligations. The transfer (2) is a
single scalar implication; orienting coordinates by their upper bounds yields
a transitive tournament with no cycles or nontrivial SCCs and adds nothing.
The faithful carrier is the rooted dyadic tower with its terminal radius and
the restriction-rank label from THM-2090. QED.
