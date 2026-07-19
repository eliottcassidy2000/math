---
id: THM-1216
title: THE COMMON-PERIOD BEAT-HOLE LAW — every nontrivial common reduced period has a forced phase-aware residue hole, with B(Q)=5A(Q)-3 and a sharp distinct-mask/run-length refinement
status: PROVED (exact cyclic-mask argument; all Q>=2; deterministic referee; Lean finite-cardinality and arithmetic core)
source: codex-2026-07-18-S81
depends_on: [THM-1192, THM-1204]
related: [THM-1178, THM-1179, THM-1182, THM-1193, HYP-7845]
script: 04-computation/lrc14_common_period_beat_hole_referee_codex_S81.py
output: 05-knowledge/results/lrc14_common_period_beat_hole_referee_codex_S81.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCCommonPeriodBeatHole.lean
script_sha256: bc6d52ff05f3303001a0f8414e86e1b7aec53f5291cae50911fac312d3360288
output_sha256: 96b01b0daeda0c214902e3e2aa0de22b0a2fe09c37a4b8f8ca7b4d03c5848f63
formalization_sha256: b148513ed97aaff67edbb77f1a82fe233a6ecf2c80a0c74cb3dd43940632a4ad
---

# THM-1216 — the common-period beat-hole law

## Statement

Let

```text
a < b1 < b2 < b3 < b4 < b5 < b6
```

be positive integers, put `q=b6-b5=dQ`, and assume

```text
Q >= 2,
gcd(bi,q)=d  for i=1,...,6.                              (1)
```

For a closed slow-carrier safe gap

```text
G_k(a)=[(14k+1)/(14a),(14k+13)/(14a)]
```

let `P_q(k)` be its consecutive block of beat numerators:

```text
P_q(k)={p in Z: p/q in G_k(a)},   N=#P_q(k).             (2)
```

Define

```text
h(Q)=ceil(Q/14),
A(Q)=2h(Q)-1,
B(Q)=5A(Q)-3=10h(Q)-8.                                  (3)
```

If `N>=B(Q)`, there is `p in P_q(k)` such that

```text
||a p/q||, ||b1 p/q||, ..., ||b6 p/q|| >= 1/14.         (4)
```

Consequently the six faster dangerous combs cannot cover that slow gap.  A
simple phase-uniform geometric supplier is

```text
6q/(7a) >= B(Q),
equivalently q/a >= 7B(Q)/6.                             (5)
```

Under (5), every slow gap has such a witness.  This is an all-scale theorem
on the stated common-period stalk.  It does **not** prove that an arbitrary
six-comb slow gap has a common reduced period, and it does not close LRC(14).

## The reduced masks and the forced hole

Write `bi=d ri`.  Assumption (1) says `gcd(ri,Q)=1`.  At a beat point
`t=p/q`, the strict-danger predicate becomes

```text
||bi p/q|| < 1/14
iff 14 min(ri p mod Q,-ri p mod Q) < Q.                  (6)
```

Thus every fast comb gives a mask `M_i` on `Z/QZ`.  Multiplication by the
unit `ri` permutes residues, so THM-1192's exact count gives

```text
#M_i=A(Q)=2ceil(Q/14)-1.                                 (7)
```

Because `b6=b5+q`, the defining masks agree: `M_6=M_5`.  A point outside

```text
U=M_1 union M_2 union M_3 union M_4 union M_5            (8)
```

makes `b5,b6` and all four complementary fast speeds safe.  Every mask also
contains residue zero.  This common point is the phase information discarded
by the scalar union bound.  Hence

```text
#U <= 1+5(A(Q)-1)=5A(Q)-4=B(Q)-1.                        (9)
```

The arithmetic threshold fits within one period.  If `2<=Q<=14`, then
`A=1`, `B=2<=Q`.  If `h=ceil(Q/14)>=2`, then

```text
Q >= 14h-13,
B=10h-8,
Q-B >= 4h-5 > 0.                                        (10)
```

Therefore `B` consecutive integers have distinct residues modulo `Q`, but
only `B-1` of all residues can lie in `U`; one lies outside.  Membership in
the closed gap makes `a` safe, proving (4).

The proof isolates the missing feature in the phase-free law: all five rows
share the same zero residue.  That four-unit overlap credit makes **every**
`Q>=2` close, including the scalar-law exceptional periods.

## The faithful sharp invariant: longest dangerous run

For the actual five masks define

```text
L(U)=maximum length of a cyclic consecutive run contained in U.   (11)
```

This is the exact phase-aware obstruction.  Any `L(U)+1` consecutive beat
numerators contain a residue outside `U`, so the sharp mask-specific supplier
is

```text
N >= L(U)+1,
or sufficiently  q/a >= 7(L(U)+1)/6.                    (12)
```

Equation (9) only says `L(U)+1<=B(Q)`.  Composite periods can be much better
because their unit masks cannot occupy long consecutive runs.  This run
length, rather than mask density, is the right finite statistic to compute
on a concrete stalk.

There is also an intermediate distinct-mask refinement.  If the five masks
in (8) have only `r` distinct values, `1<=r<=5`, then

```text
#U <= 1+r(A(Q)-1),
B_r(Q)=2+r(A(Q)-1),                                      (13)
```

and `N>=B_r(Q)` suffices.  For `Q>=15`, the masks are indexed by unit sign
classes.  More precisely, let

```text
S_Q={x mod Q:14 min(x,-x)<Q}.
```

For `Q<=14`, `S_Q={0}` and all unit masks agree.  For `Q>=15`, its unit
stabilizer is exactly `{+1,-1}`.  To prove this, note that `S_Q` is a cyclic
interval of length `A` with `2A<=Q`.  If a unit `u` stabilizes it, multiplication
by `u` carries `S_Q intersect (S_Q+1)` to `S_Q intersect (S_Q+u)`.  The first
intersection has size `A-1`; a cyclic interval shorter than a semicircle has
an intersection this large only at shifts `+1` and `-1`.  Thus two masks
agree exactly when their normalized unit residues agree up to sign.

In particular,

```text
r <= min(5,phi(Q)/2)  for Q>=3,                          (14)
```

and if all normalized speeds are congruent up to sign, then `r=1`,

```text
N>=A(Q)+1=2ceil(Q/14),
q/a>=7ceil(Q/14)/3                                      (15)
```

suffices.

## The sharpened q=14 corollary

At `Q=14`, every mask is the singleton `{0}`.  Equations (13)--(15) give

```text
N>=2,
q/a>=7/3.                                               (16)
```

This strictly strengthens THM-1204, which assumed `q/a>=7`, obtained `N>=6`,
and then invoked the phase-free five-cap contradiction.  THM-1204 is now an
immediate weaker corollary of (16).

The equality-scale row

```text
(a;b1,...,b6)=(6;9,11,13,15,17,31),   q=14
```

has exactly two beat numerators in each slow gap and an exact lonely beat
point in each.  The two-point threshold is phase-uniformly sharp: for
`(a,d,k)=(13,2,6)`, `q=28` and the middle beat block is the singleton `{14}`,
the common dangerous residue.  More generally, in the family `a=6d+1` with
even `d`, the odd carrier centers the middle slow gap at `t=1/2` and its beat
block is the singleton `{7d}`.  The ratios `q/a=14d/(6d+1)` approach `7/3`
from below.

## Exact classification of the phase-free scalar tail

The weaker THM-1192 cap is

```text
U(N,Q)=floor(N/Q)A+min(N mod Q,A).                       (17)
```

Its common-period necessary law is `N<=5U(N,Q)`.  Put

```text
beta=Q-5A.
```

An eventual scalar contradiction exists exactly when `beta>0`, namely

```text
Q in {6,7,...,14} union {16,17,18,...}.                  (18)
```

For `Q=1,...,5` and `Q=15`, every multiple `N=mQ` passes the scalar law, so
there is no phase-free tail.  When `beta>0`, write

```text
m*=floor(4A/beta),
N_last=m*Q+min(Q-1,5A-m*beta),
N_pf=N_last+1.                                          (19)
```

Then `N_last` is the exact last passing size and every `N>=N_pf` contradicts
the scalar law.  Indeed, for `N=mQ+r`,

```text
N-5U = m beta-4r        if 0<=r<=A,
N-5U = m beta+r-5A      if A<r<Q.                       (20)
```

Formula (19) follows by locating the last nonpositive value.  For example,
the scalar thresholds at `Q=6,14,16,17,28,29` are respectively
`26,6,196,106,16,151`, while the exact uniform thresholds (3) are
`2,2,12,12,12,22`.  The disparity measures phase information, not a small
rounding loss.

## Why Q=1 is the unique exact obstruction

At `Q=1`, the reduced clock has only residue zero and every beat numerator is
dangerous for every fast speed.  No beat point can witness noncoverage.  This
is a counterexample to this **method**, not a covering construction and not a
counterexample to LRC.  Equations (9)--(10) show that it is the only common
reduced period without a forced exact residue hole.

## Mixed-period forward bridge

The finite-union proof is deliberately stated independently of equal mask
periods.  On a master clock `L`, if five lifted masks have sizes `C_i` and a
common point, their union has size at most

```text
1+sum_i(C_i-1).                                         (21)
```

For arbitrary `g_i=gcd(b_i,q)`, `Q_i=q/g_i`, and
`L=lcm_i Q_i`, one has `C_i=L A(Q_i)/Q_i`; the defining pair still shares one
mask.  Thus (21) is the phase-preserving consumer needed for a future
mixed-gcd theorem.  THM-1216 keeps the common-`Q` hypothesis so that its
geometric threshold is transparent; it makes no unproved mixed-period claim.

## Structural and tournament audit

The useful object is not a tournament on runners.  It is the nerve of the
five reduced danger masks, marked by the universal residue-zero vertex and
decorated by the consecutive block phase.  This quotient preserves the exact
beat-danger predicate and forced overlap.  Cardinality forgets the placement
of holes; `L(U)` restores precisely the placement needed by a consecutive
block.  Passing only to the phase-free caps destroys both overlap and run
length.

For the required Tournament Analysis audit, take common-zero intersection as
the pairwise observable.  It is symmetric, so increasing speed is merely a
tie gauge with Hamiltonian path

```text
b1 -> b2 -> b3 -> b4 -> b5 -> b6.
```

The gauge tournament is transitive: score histogram `{0,1,2,3,4,5}`, zero
directed cycles, six singleton SCCs, one Hamiltonian path, and fifteen edge
flips under the reverse gauge.  None of those fingerprints proves the
theorem.  The proof lives in the common intersection of five masks, an
incidence feature that a binary orientation erases.

Vertices were challenged as runners, slow gaps, fixed circle sections,
section boundaries, wall events, beat numerators, residue masks, cover arcs,
and proof obligations.  Residue masks plus the block-phase sidecar are the
faithful choice.  The quotient preserves exact danger at beat points but
destroys behavior away from that rational grid, so it cannot by itself close
an arbitrary continuous slow-gap branch.

## Formalization boundary

`LRCCommonPeriodBeatHole.lean` kernel-checks:

```text
Q>=2 -> B(Q)<=Q,
B(14)=2,
|M0 union ... union M4|+4 <= sum_i |Mi|
  when all Mi contain one common point,
|Mi|<=A(Q) and |block|>=B(Q) -> block is not covered.
```

The typed consumer leaves the real-gap block supplier and the speed-gcd unit
reduction as explicit hypotheses.  It contains no hidden geometric premise,
no `sorry`, and no floating arithmetic.

## Exact replay

```text
python3 04-computation/lrc14_common_period_beat_hole_referee_codex_S81.py
python3 -O 04-computation/lrc14_common_period_beat_hole_referee_codex_S81.py
```

The normal and optimized outputs are byte-identical.  The frozen referee
checks the mask count and budget through `Q=10000`, all unions of up to five
distinct masks through `Q=35`, the stabilizer classification through `Q=100`,
the sharp scalar tail formula through `Q=300`, and the exact q=14 examples.
