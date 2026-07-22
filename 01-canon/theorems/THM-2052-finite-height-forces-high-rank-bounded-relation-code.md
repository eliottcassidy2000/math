---
id: THM-2052
title: Finite LRC height forces a high-rank bounded relation code
status: PROVED, conditional only through THM-763 on the settled twelve-speed LRC. Every primitive hypothetical LRC(14) counterexample has at least nine independent support-at-most-five relations of coefficient height at most 2^20, and at least eleven independent support-at-most-three relations of height at most 91^6. Hence it lies in a finite atlas of rational subspaces of dimension at most two. This is a dimension descent, not a classification of those subspaces or a proof of LRC(14).
source: codex-2026-07-21-LRC-unrelated-transfer
depends_on:
  - THM-763
related:
  - THM-935
  - THM-2051
  - HYP-8841
---

# THM-2052 -- finite height forces a bounded relation code

## 1. General pigeonhole--code lemma

Let `v=(v_1,...,v_n)` be positive integers with

```text
sum_i v_i<=B.
```

Fix integers `2<=s<=n` and `Q>=1` such that

```text
(Q+1)^s>QB+1.                                         (1)
```

Let `W_(Q,s)(v)` be the rational span of all integer coefficient vectors
`k in Z^n` satisfying

```text
k dot v=0,       2<=|supp(k)|<=s,       ||k||_infinity<=Q.   (2)
```

**Theorem.** Under (1),

```text
dim_Q W_(Q,s)(v)>=n-s+1.                              (3)
```

*Proof.* Fix any coordinate set `S` of size `s`. The `(Q+1)^s` coefficient
vectors `a in {0,...,Q}^S` produce integer sums

```text
sum_(i in S) a_i v_i in {0,...,Q sum_(i in S)v_i}
                       subseteq {0,...,QB}.
```

By (1), two coefficient vectors give the same sum. Their difference is a
nonzero vector `k` supported in `S`, of height at most `Q`, and orthogonal to
`v`. Positivity of the speeds rules out support one, so `k` belongs to (2).
Thus

```text
W_(Q,s)(v) intersect Q^S != {0}                       (4)
```

for every `s`-coordinate set `S`.

Put `U=W_(Q,s)(v)^perp` and choose a matrix whose rows form a basis of `U`.
Its coordinate columns have the property that every `s` of them are linearly
dependent: a dependence is exactly a nonzero vector in (4). If `dim U>=s`,
the matrix would contain `s` independent columns, a contradiction. Hence

```text
dim U<=s-1,
dim W_(Q,s)(v)=n-dim U>=n-s+1,
```

which proves (3). QED.

## 2. Exact LRC(14) specialization

Let `v` be a primitive hypothetical thirteen-speed counterexample. Since
`M(v)<1/14`, THM-763 with the settled twelve-speed LRC gives

```text
sum_i v_i<=91^12
 =322475487413604782665681=:B.                        (5)
```

First take

```text
Q=2^20=1048576,       s=5.
```

The exact comparison

```text
Q^4=1208925819614629174706176>B
```

implies `(Q+1)^5>QB+1`. Therefore

```text
dim W_(2^20,5)(v)>=9.                                 (6)
```

This uses exactly the height of Theorem A in THM-2051, but forces nine
independent short relation rows rather than merely one.

For a lower-dimensional atlas, take

```text
Q=91^6=567869252041,       s=3.
```

Now `B=Q^2`, so

```text
(Q+1)^3>Q^3+1=QB+1.
```

The theorem gives

```text
dim W_(91^6,3)(v)>=11.                                (7)
```

All possible coefficient rows in (7) belong to a finite set. Choosing eleven
independent ones from that set shows that every counterexample lies in one of
finitely many rational kernels of dimension at most two. This is an exact
finite atlas of two-parameter ambient subspaces, although each subspace still
contains infinitely many primitive positive integer points.

## 3. The rank-twelve terminal

There is a clean endpoint for any future relation-harvesting descent. Suppose
the same primitive `v in Z_{>0}^13` satisfies twelve independent relation rows,
each supported on at most `s` coordinates and of height at most `Q`. Put them
in a `12 x 13` matrix `K`. The signed maximal minors span the one-dimensional
kernel of `K`. Dividing their common gcd gives the primitive generator `v`.
Every row of each minor has Euclidean norm at most `sqrt(s)Q`, so Hadamard gives

```text
max_i v_i<=(sqrt(s)Q)^12=s^6 Q^12.                    (8)
```

At rank twelve the remaining LRC decision is therefore finite and exact by
the pair-sum ruler theorem. For the THM-2051 support ceiling `s=5`, this is
`max v_i<=5^6 Q^12`.

## 4. Consequence for the termination program

The useful state variable is not “a small relation was seen.” Repetition of
one relation does not descend. A valid HYP-8841 mechanism should track

```text
rank of the harvested bounded-relation code,
```

and at each unresolved step supply either a strict lonely exit or a new
independent row satisfied by the same speed vector. Equations (6)--(7) start
that rank at `9` or `11`; equation (8) makes rank `12` a finite terminal. The
remaining gap is therefore only one independent relation above the triple
atlas, provided the phase/owner data can force it.

This is coding theory rather than Tournament Analysis. The faithful vertices
are coordinate columns, and the preserved predicate is dependence of every
`s`-column set. Any binary orientation discards the signed coefficients and
relation height, so tournament scores, cycles, SCCs, and Hamiltonian paths do
not encode (1)--(8).
