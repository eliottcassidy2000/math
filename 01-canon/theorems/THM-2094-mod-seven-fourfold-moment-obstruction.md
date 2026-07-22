---
id: THM-2094
title: "Four seven-divisible terminal speeds are excluded by a conditional moment certificate"
status: >
  PROVED. If 7 does not divide the odd guard h and exactly four of seven
  terminal speeds are divisible by 7, the guard indicator is independent of
  the complete multiplicity distribution of those four danger combs. Combining
  that grid orthogonality with THM-2091's sharp three-charge ledger and the
  THM-1234 five-comb pair floor gives incompatible upper and lower bounds on
  one symmetrized quadratic moment, separated by 107/315315. Hence such a
  seven-comb safe set cannot lie in the guard. This removes the k=4 branch of
  THM-2086's all-height residual, but does not settle k=1,2,3, the remaining
  finite banks, or LRC(14).
source: codex-2026-07-22-LRC-mod7-conditional-moment
depends_on:
  - THM-1234
  - THM-2086
  - THM-2091
related:
  - THM-1122
  - THM-1166
  - THM-2080
script: 04-computation/lrc14_mod7_fourfold_moment_obstruction_codex_20260722.py
output: 05-knowledge/results/lrc14_mod7_fourfold_moment_obstruction_codex_20260722.out
script_sha256: c4c9b8096680483e99ab54f74706a2ba51f6b4cbe194276e243d1c908e795a73
output_sha256: 4844526e58f57db13744cf1be90c5f6563045cedaa76007f3eb8fa81ca0887ad
hash_basis: repository blobs with LF line endings
---

# THM-2094 -- the modulo-seven fourfold branch is impossible

Let `h` be odd with `7` not dividing `h`, and let

```text
Q=A disjoint_union B,
|A|=4,       |B|=3,
7|q for q in A,        7 does not divide q for q in B. (1)
```

As in THM-2091, write

```text
H=1_(E_h),          X_q=1_(D_q),
U=sum_(q in A)X_q,  V=sum_(q in B)X_q.                 (2)
```

Then

```text
measure(G_Q minus E_h)>0.                              (3)
```

In particular, `G_Q` is not contained in `E_h`.

## 1. Complete grid orthogonality

Every `X_q` with `q in A` is `1/7`-periodic. More generally, if `F` is any
integrable `1/7`-periodic function, then

```text
integral F H=(2/7) integral F.                         (4)
```

To prove (4), average over the seven shifts `t+j/7`. Since multiplication by
`h` permutes the nonzero residue classes modulo seven, away from the finite
guard-boundary set exactly two of

```text
h(t+j/7),       0<=j<=6,
```

lie within circle distance `1/7` of zero. Thus

```text
sum_(j=0)^6 H(t+j/7)=2                                 (5)
```

almost everywhere. The shifted values of `F` are all equal, so integration
gives (4).

Taking `F=1_(U=u)` proves the full finite-law independence

```text
integral (H-2/7) 1_(U=u)=0,          0<=u<=4.          (6)
```

This is stronger than the zero first-order charges in THM-2091: every
function of the four-speed multiplicity `U` is orthogonal to the guard.

The single-comb means are

```text
integral H=2/7,       integral U=4/7,
integral V=3/7.                                        (7)
```

## 2. The three remaining charges

THM-2091's three largest possible negative odd-guard charges have sum

```text
C_3=5/294+8/539+3/245=713/16170.                      (8)
```

The three speeds in `B` are distinct, so their reduced ratios are distinct
and that charge spectrum gives

```text
integral H V
 =sum_(q in B) measure(E_h intersect D_q)
 >=6/49-C_3
 =181/2310.                                            (9)
```

## 3. A five-comb lower moment

Choose uniformly three members of `A` and two members of `B`. THM-1234 says
the total global pair overlap in every resulting five-set is at least
`44/273`. Averaging those pair ledgers gives

```text
P(U,V)=(1/2) binom(U,2)+(1/2)UV+(1/3)binom(V,2),
integral P(U,V)>=44/273.                               (10)
```

Indeed a fixed `A-A`, `A-B`, or `B-B` pair is selected with probability
`1/2`, `1/2`, or `1/3`, respectively. This is a lossless symmetrization of
the five-subset pair obligations, not a sampled-family assertion.

## 4. The conditional quadratic majorant

Suppose for contradiction that `G_Q` is contained in `E_h` up to a null set.
Then outside the guard at least one danger comb is active:

```text
H=0  implies U+V>=1                                    (11)
```

almost everywhere.

For `0<=u<=4`, put

```text
alpha_u=(-91/12,-65/12,-121/20,-443/60,-113/12)_u.    (12)
```

On every state

```text
H in {0,1},       0<=u<=4,       0<=v<=3,
```

except the state `(H,u,v)=(0,0,0)` excluded by (11), the following exact
pointwise majorant holds:

```text
P(u,v)+(65/12)H-(65/42)u-(13/6)v
 +alpha_u(H-2/7)+(4/3)Hv <=0.                          (13)
```

After multiplying the negative of the left side by `420`, its nonnegative
slack tables, with columns `v=0,1,2,3`, are

```text
H=0:
u=0     X      0    770   1400
u=1     0    700   1260   1680
u=2   364    854   1204   1414
u=3   434    714    854    854
u=4   210    280    210      0

H=1:
u=0     0    350    560    630
u=1     0    140    140      0
u=2   630    560    350      0
u=3  1260    980    560      0
u=4  1890   1400    770      0.                       (14)
```

Here `X` marks only the forbidden state. Thus (14) is a complete exact proof
of (13), with no continuous relaxation or numerical optimization.

Integrate (13). Equation (6) kills the `alpha_U` term. Equations (7) and (9)
then give

```text
integral P(U,V)
 <=-(65/12)(2/7)+(65/42)(4/7)+(13/6)(3/7)
    -(4/3)(181/2310)
 =3901/24255.                                          (15)
```

But

```text
44/273-3901/24255=107/315315>0,                        (16)
```

contradicting (10). Therefore even almost-everywhere containment is
impossible, which proves the positive-measure conclusion (3). QED.

## 5. Transfer, scope, and assumption challenge

THM-2086 reduced its all-height relative-Hunter residual to

```text
7 does not divide h,
1<=#{q:7|q}<=4,
nonlacunary terminal packets.
```

This theorem removes the upper endpoint and leaves only one, two, or three
seven-divisible speeds in that analytic residual. THM-2092/2093 separately
put the no-pair relation-code branches in finite boxes; neither result decides
the remaining finite rows.

The mechanism came from combining two repo carriers that were previously
used separately. THM-1122's multiplicity LP supplies the finite-state
majorant, while THM-2086's vanishing seventh Fourier mode upgrades a zero
charge into the full conditional law (6). Keeping only the first moment of
`U` leaves a feasible hostile relaxation; the complete grid fiber is
load-bearing.

Candidate tournament vertices were the seven runners, the twelve `3+2`
five-subset obligations, the five multiplicity levels of `U`, and the forty
`(H,u,v)` states. No pairwise orientation preserves both (6) and (13).
The faithful carrier is the bipartition `A disjoint_union B` together with
the guard-conditioned multiplicity table. A tournament on runners would
discard exactly the residue class that proves the contradiction.

## Exact referee

The companion verifies the seven-shift count, the exact top-three charge
sum, all forty entries of the scaled slack table, the integrated upper bound,
and the strict rational gap (16). Runtime checks remain active under
optimization; normal and `python -O` transcripts byte-match the stored output
and end in `PASS`.
