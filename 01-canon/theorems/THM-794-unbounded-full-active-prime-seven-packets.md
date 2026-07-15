---
id: THM-794
title: Unbounded full-active prime-seven packets at fixed fastest ratio
status: PROVED (exact infinite family, wall order, coverage, holonomy, and extent refutation) + VERIFIED (Fraction replay for H=2 through 200)
source: codex-2026-07-14-S10
depends_on:
  - THM-783   # anchored simple-wall extension
related: [THM-779, THM-784, THM-786, THM-788, HYP-6840, HYP-6845, HYP-6850, MISTAKE-149]
verification:
  - 04-computation/lrc14_full_active_packet_holonomy_codex_S10.py
  - 05-knowledge/results/lrc14_full_active_packet_holonomy_codex_S10.out
---

# THM-794 — Unbounded full-active prime-seven packets

## Statement

For every integer `H>=2`, put

```text
F=49H+1,                 w_j=F-7j,  0<=j<=7,
f=w_0=F,                 g=w_1=F-7.
```

The eight owners are distinct, positive, and all congruent to `1 mod 7`.
For

```text
6H <= m <= 7H-2
```

define consecutive `f`-walls and one wall of each slower owner by

```text
x_m     = (m+1/2)/F,
y_(j,m) = (m+3/2-j)/(F-7j),       1<=j<=7.              (1)
```

Then the complete global wall order on `[x_m,x_(m+1)]` is

```text
x_m < y_(7,m) < y_(6,m) < ... < y_(1,m) < x_(m+1),     (2)
```

and **every wall in (2) is covered** in the prime-seven sheet problem.
Consequently a maximal blocking run contains the whole chain

```text
[x_(6H),x_(7H-1)],
```

which has

```text
H-1 consecutive active f-periods,
8H-7 covered global walls,
extent L_H=(H-1)/(49H+1).                               (3)
```

In particular:

1. the active-period count `A` in THM-788 is unbounded even while
   `R=ceil(f/g)=2` is fixed;
2. no bound on `A` can depend only on `R`, or even only on the labelled gap
   pattern `w_i-w_j=7(j-i)`;
3. for every `H>=5`,

   ```text
   L_H > 1/g + 2/f,                                     (4)
   ```

   so THM-786's proposed universal extent inequality is false.

This does **not** refute THM-788's proved counting inequalities.  It refutes
the proposed use of `A` as a compactness target.

## Proof

Write a wall of owner `u` with index `n` as

```text
z_(u,n)=(n+1/2)/u.
```

Thus `y_(j,m)` in (1) is the `w_j`-wall with index `m+1-j`.

### 1. Exact global event order

First compare `x_m` with `y_(j,m)`.  Since all denominators are positive,
cross-multiplication gives

```text
x_m < y_(j,m)
  iff 7j(m+1/2) > (j-1)F.                               (5)
```

The right side is strongest at `j=7`, and

```text
m+1/2 >= 6H+1/2 > 6H+6/49 = 6F/49.
```

So (5) holds for every `j`.

Put `a=m+3/2`.  For `2<=j<=7`, a second cross-multiplication gives

```text
y_(j,m) < y_(j-1,m)  iff  a < F/7.                      (6)
```

The comparison `y_(1,m)<x_(m+1)` gives exactly the same condition.  But

```text
a <= 7H-1/2 < 7H+1/7 = F/7,                             (7)
```

which proves all inequalities in (2).

Every slower wall mesh is

```text
1/w_j > 1/F=x_(m+1)-x_m.
```

Therefore a slower owner has at most one wall in this `f`-period.  Formula
(2) exhibits one from each of the seven slower owners, so these eight displayed
coordinates are exactly the full global event list.  Their inequalities are
strict, hence all events are simple.

### 2. The first fastest wall is covered

The wall preceding `y_(j,m)` is `y_(j,m)-1/w_j`.  Since

```text
0<y_(j,m)-x_m<1/F<1/w_j,
```

`x_m` lies strictly between these two consecutive `w_j`-walls.  The nearest
integer to `w_j x_m` is therefore

```text
n_j=m+1-j.                                               (8)
```

For a non-walling owner `u`, THM-779's sheet token is

```text
k_u(x)=-u^(-1) round(ux)  (mod 7).
```

All `w_j` have inverse residue one.  Hence at `x_m`,

```text
k_(w_j)(x_m)=-(m+1-j)=j-m-1  (mod 7).                   (9)
```

As `j` runs from one through seven, (9) runs through every element of
`F_7`.  The seven non-walling owners form a perfect sheet rainbow, so `x_m`
is covered.

### 3. Anchored extension covers the whole packet

The ordered wall-index word in (2), including both fastest endpoints, is

```text
m, m-6, m-5, m-4, m-3, m-2, m-1, m, m+1.               (10)
```

Modulo seven, each entry in (10) is its predecessor plus one.  Every owner
inverse is also one.  Thus each consecutive pair satisfies THM-783's exact
anchored extension congruence

```text
phi_next = phi_current + w_current^(-1)  (mod 7).
```

The first wall is covered by Step 2 and every next global coordinate is simple
by Step 1.  Induction through (10) proves that every wall in (2), including
`x_(m+1)`, is covered.

This argument applies consecutively for `m=6H,...,7H-2`.  There are `H-1`
such complete periods, each with seven visitors.  Starting with the first
`f`-wall and adding eight new walls per period gives

```text
1+8(H-1)=8H-7.
```

The endpoint difference is

```text
x_(7H-1)-x_(6H)=(H-1)/F.
```

This proves (3).

### 4. Fixed ratio and the failed extent conjecture

For every `H>=2`,

```text
1 < F/(F-7) < 2,
```

so `ceil(f/g)=2`.  Also

```text
w_i-w_j=7(j-i),
```

independent of `H`.  Yet (3) has `A>=H-1`, proving the first two
consequences.

Finally, (4) is equivalent to

```text
(H-3)g > F.                                              (11)
```

When `H>=5`, `H-3>=2`, while

```text
2g=98H-12 > 49H+1=F.
```

Thus (11) holds and the universal extent conjecture is refuted. ∎

## The structural correction: full-support holonomy is also gauge

THM-784 showed that counting raw fastest walls fails because arbitrarily many
**empty** fastest periods can refine a fixed slow rainbow chamber.  THM-788
contracted those empty periods and retained active periods.  The present
family shows the complementary failure mode: every retained period is active
in the strongest possible way, but it repeats a central gauge motion.

The visitor-incidence row in every period is

```text
(1,1,1,1,1,1,1),
```

and its inverse-residue sum is `7=0 mod 7`.  The owner word is the fixed cycle

```text
f -> w_7 -> w_6 -> ... -> w_1 -> f.                     (12)
```

Between matching sides of two consecutive `f`-walls, every owner crosses
exactly one wall.  Since every inverse residue is one, every token changes by
`-1`.  The whole collision state is therefore translated by the diagonal
element of `F_7`; after quotienting common sheet translation, (12) is a closed
eight-step holonomy cycle.  Repeating it creates no new normalized obstruction
state.

The next quotient must therefore contract at least two free dynamics:

1. empty same-owner refinement (THM-784); and
2. full-support packet cycles whose return map is a diagonal sheet
   translation (this theorem).

A promising finite coordinate is not the number of active packets but the
number of transitions between SCCs of the collision transducer **after**
diagonal deck translation is divided out, with metric lift data retained for
intersection with the core-safe component.  No bound on that corrected
coordinate is claimed here.

## Tournament Analysis and challenged vertex sets

Several tournament quotients expose exactly what is lost.

- **Packet occurrences as vertices.**  Chronological precedence gives the
  transitive tournament `T_(H-1)`: score histogram `0,...,H-2`, no directed
  cycles, singleton SCCs, and one Hamiltonian path.  It records repetition
  count but cannot see that every occurrence is the same normalized state.
- **Visitors inside one packet.**  Chronological order gives `T_7`, again with
  zero directed cycles and one Hamiltonian path.
- **Sheet/visitor labels modulo a marked cut.**  Orient `a->b` when
  `b-a in {1,2,3} mod 7`.  The cut-free tournament is the regular `R_7`:
  score histogram `{3:7}`, 14 directed triangles, one SCC, and 175 Hamiltonian
  paths.  Formula (9) changes only by a global rotation from period to period,
  so this tournament is constant.  Moving the marked cut by one turns its
  linear representative through six labelled edge flips.

Runner vertices preserve the fixed speed-gap pattern but forget event phase;
event vertices retain chronology but turn the repeated packet into a longer
transitive path; packet vertices retain full visitor support but forget the
diagonal return map.  The faithful object is the owner-labelled collision
cycle together with its deck holonomy and metric embedding.  This explicitly
challenges the assumption that either runners, raw walls, or merely nonempty
periods are canonical vertices for the exit problem.

## Exact replay

The verifier uses `fractions.Fraction` only.  For every `2<=H<=200` it
reconstructs the full event list in the certified interval, checks all strict
orders and nearest-integer tokens directly, verifies the rainbow condition at
every wall and in every intervening chamber, and checks (3), the fixed-ratio
claim, and every `H>=5` instance of (4).  It also computes the displayed
tournament fingerprints rather than treating them as labels.
