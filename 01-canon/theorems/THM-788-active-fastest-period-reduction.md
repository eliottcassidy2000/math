---
id: THM-788
title: Active-fastest-period reduction for prime-seven blocking runs
status: PROVED (elementary wall-mesh counting plus THM-783's period-sum law; independently audited, including 2,000 random tuples / 4,626 covered runs with zero inequality violations)
source: codex-2026-07-14-S10 with balanced-cascade parallel analysis; numbered after concurrent THM-785/786/787 claims
depends_on:
  - THM-779   # simple covered-wall criterion
  - THM-783   # period-sum and single-visitor laws, at corrected scope
related: [THM-778, THM-784, THM-786, HYP-6840, HYP-6845, HYP-6850, MISTAKE-147]
---

# THM-788 — active-fastest-period reduction

## Theorem

Let `W` be eight distinct positive owners, all nonzero modulo 7.  Let `f` be
the fastest owner, `g` the second fastest, and

`R = ceil(f/g)`.

Consider a prime-seven blocking run: consecutive coordinates in the **full
global event-coordinate order**, every coordinate simple and satisfying
THM-779's wall-rainbow condition. No intervening event coordinate is deleted.
Let

- `K` be its number of walls;
- `L` be the distance from its first to its last wall;
- a **complete f-period** be the interval between consecutive f-walls lying in
  the run;
- such a period be **active** if it contains a non-f wall; and
- `A` be the number of active complete f-periods.

Then

`K <= R(A+1) + 7A + 14`,                                      (1)

and

`L < (R(A+1)+1)/f < (A+1)/g + (A+2)/f`.                      (2)

Every active period contains between two and seven visitors, and their inverse
residues have sum zero modulo 7.

Consequently, an absolute bound `A <= C` would imply the scale-sensitive exit
bounds

`L < (C+1)/g + (C+2)/f`,

`K <= (C+1)ceil(f/g) + 7C + 14`.

Thus the corrected cascade target is the number of **active** fastest periods
after empty fastest-owner refinement has been contracted, not raw wall count.

## Proof

### 1. Same-owner blocks

A slower owner `c<f` has wall mesh `1/c > 1/f`.  Between two consecutive
c-walls there is therefore an f-wall.  Hence no slower owner can occur twice
consecutively in the global wall schedule.

Suppose a block contains `B` consecutive f-walls.  Its first and last wall are
separated by `(B-1)/f`, and there is no g-wall in the closed interval between
them.  Simultaneous endpoints are impossible in a blocking run.  Every closed
interval of length at least `1/g` meets the g-wall lattice, so

`(B-1)/f < 1/g`.

Therefore `B <= ceil(f/g)=R`.

### 2. Count the f-walls

Let `N_f` be the number of f-walls in the run.  If `N_f>=2`, the `N_f-1`
complete f-periods consist of `A` active periods and strings of empty periods.
There are at most `A+1` empty strings.  A string of `s` empty periods is a
block of `s+1` consecutive f-walls, so Step 1 gives `s<=R-1`.  Hence

`N_f-1 <= A + (A+1)(R-1) = R(A+1)-1`,

or

`N_f <= R(A+1)`.                                               (3)

The same inequality is trivial when `N_f` is 0 or 1 because `R>=2`.

### 3. Count the visitor walls

Each active f-period has length `1/f`, shorter than every slower wall mesh, so
it contains at most one wall from each of the seven slower owners: at most
`7A` interior visitors in total.

When the run contains an f-wall, the part before its first f-wall lies after
the preceding f-wall, and the part after its last f-wall lies before the next
f-wall; each edge interval has length less than `1/f`.  Each contains at most
one wall from each slower owner, giving at most 14 edge walls.  If the run has
no f-wall, it lies inside a single f-gap and has at most seven walls, which is
even stronger.  Combining this with (3) proves (1).

### 4. Metric extent

If the run has at least one f-wall, its two edge pieces have total length less
than `2/f`, while the span from its first to last f-wall is `(N_f-1)/f`.
Thus

`L < (N_f+1)/f <= (R(A+1)+1)/f`.

With no f-wall, `L<1/f`, so the same displayed bound holds.  Finally
`R=ceil(f/g)<f/g+1`, and substitution gives the second strict inequality in
(2).

### 5. Zero-sum content

Every active complete period contains at least one visitor and at most one
from each of seven slower owners.  THM-783's period-sum law says that their
inverse residues sum to zero.  A singleton inverse residue is nonzero modulo
7, so the visitor count is at least two and at most seven.  ∎

## Why this is the right quotient

THM-784 has `A=0`: its `21N` fast walls are a long empty refinement inside one
fixed slow rainbow chamber.  Formula (1) correctly allows raw wall count to
grow like `f/g`, while (2) keeps metric extent at the slow scale.

For Tournament Analysis, take wall events as vertices, chronological order as
the pairwise observable, and time orientation as the gauge.  The event
tournament is transitive: score histogram `{0,...,K-1}`, no directed cycles,
singleton SCCs, and a unique chronological Hamiltonian path.  Those
fingerprints do not distinguish an empty refinement from genuine interaction.
For fixed `K`, the bare tournament cannot distinguish empty refinement from
genuine interaction (although its vertex count of course records `K`). The
exact normal form for `N_f>=1` is

`E_0, V_1, E_1, ..., V_A, E_A`,

where each `E_j` is an f-block of 1 through `R` f-walls, absorbing at most
`R-1` empty periods, and each `V_j` is an ordered two- through seven-owner
zero-sum visitor packet. There are also two edge packets of at most seven walls
each. Contract the `E_j` blocks while decorating them by absorbed length; the
`V_j` packets are fibre/hyperedge data not supplied by the transitive
tournament.

The retained proof object is therefore the owner-labelled tie-free Hamiltonian
path together with one absolute slow scale (`g`, or `1/g`), the ratio `f/g`,
the decorated f-block word, and each active packet's zero-sum visitor set.

This challenges the assumption that event vertices at the finest available
resolution are canonical. Decorated period/block vertices preserve the
LRC-relevant visitor predicate and slow metric; undecorated period vertices do
not. The raw event tournament obscures both by treating arbitrarily fine
same-owner subdivision as new complexity.
