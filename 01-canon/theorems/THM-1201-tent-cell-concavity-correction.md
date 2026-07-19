---
id: THM-1201
title: TENT-CELL CONCAVITY CORRECTION — between consecutive zeros of the combined runner arrangement, the lower envelope g(t)=min_v ||vt|| is concave, so the triangle through a cell maximum lies below g and area is at least H L / 2, not at most H L / 2; the proposed THM-1195 integral certificate is therefore withdrawn
status: PROVED exact structural correction; the old cell inequality is refuted inside the actual 13-speed tight family {1,...,11,13,24}; no implication from integral g >= 1/28 to LRC(14) is established
source: codex-S78 (audit while closing the LRC(14) slow-gap frontier)
depends_on: [THM-1195 (withdrawn criterion), THM-1200 (its tent/independence paragraph is likewise withdrawn), MISTAKE-174]
scripts: 04-computation/lrc14_tent_cell_direction_audit_codex_S78.py -> 05-knowledge/results/lrc14_tent_cell_direction_audit_codex_S78.out
---

# THM-1201 — the zero-cell envelope is concave

Let

\[
g(t)=\min_{v\in V}\lVert vt\rVert
\]

and let \([L,R]\) be a cell between consecutive points of the combined zero
set \(\{k/v:v\in V,\ k\in\mathbb Z\}\).  No runner has a zero in the open
cell.  Therefore each restriction \(t\mapsto\lVert vt\rVert\) is one
triangular arch on \([L,R]\), hence is concave.  A pointwise minimum of
finitely many concave functions is concave: its hypograph is the intersection
of their convex hypographs.  Thus \(g|_{[L,R]}\) is concave.

Both endpoints are combined zeros, so \(g(L)=g(R)=0\).  If \(g(x)=H\) is a
cell maximum, concavity puts the two chords

\[
(L,0)\longrightarrow(x,H)\longrightarrow(R,0)
\]

**below** the graph.  Integrating those chords gives the correct universal
direction

\[
\boxed{\int_L^R g(t)\,dt\ \geq\ {H(R-L)\over2}.}
\]

Equality holds exactly when the lower envelope itself is the two-sided affine
tent (up to harmless collinear subdivisions).  There is no corresponding
upper bound \(HL/2\); the trivial upper bound is \(HL\).

## Exact 13-speed counterexample to the old direction

Take the tight LRC(14) family

\[
V=\{1,2,\ldots,11,13,24\}
\]

and its consecutive-zero cell \([1/24,1/13]\).  Exact arrangement reduction
gives

\[
g(t)=
\begin{cases}
24t-1,&1/24\leq t\leq1/23,\\
t,&1/23\leq t\leq1/14,\\
1-13t,&1/14\leq t\leq1/13.
\end{cases}
\]

The slopes \(24,1,-13\) are nonincreasing, visibly exhibiting concavity.  Its
height, length, area, and the former claimed upper bound are

\[
H={1\over14},\qquad R-L={11\over312},\qquad
\int_L^R g={185\over100464},\qquad
{H(R-L)\over2}={11\over8736}.
\]

Consequently

\[
\int_L^R g-{H(R-L)\over2}={3\over5152}>0,
\qquad
{\int_L^R g\over H(R-L)/2}={370\over253}.
\]

This is a counterexample to the proof premise in the same 13-runner setting
where THM-1195 was asserted, not a lower-dimensional analogy.

## What survives

The exact integral routine in the THM-1195 script still computes
\(\int_0^1g\) correctly.  Its labels “CERTIFIES” and the interpretation of
\(1/28\) as a sufficient threshold do not survive.  The numerical identity
between the independence heuristic \(1/(2n)\) and the *proposed* tent
threshold is algebraically true but no longer compares two valid estimates.

The corrected concavity law supplies a different statistic:

\[
\int_0^1g(t)\,dt\ \geq\ {1\over2}
  \sum_C |C|\,H_C.
\]

It is a lower bound on the length-weighted local-max profile.  It does not turn
a large integral into a \(1/14\) pointwise witness.  Any future integral route
must add genuine control of cell shapes or the distribution of the heights;
the zero partition alone has the opposite curvature from that used in
THM-1195.

## Tournament / carrier audit

The natural vertices here are **zero cells**, not runners.  A pairwise runner
crossing only determines where the active carrier of the lower envelope may
switch.  Quotienting to an active-carrier tournament preserves switch order
and slope descent, but destroys cell lengths and the height/area data needed
for an integral implication.  On the counterexample cell the switch path is

\[
24\longrightarrow1\longrightarrow13,
\]

with slope word \((24,1,-13)\).  There is no directed cycle and no meaningful
tie Hamiltonian multiplicity on this three-vertex path.  The challenged
assumption is precisely that runner vertices, or their pairwise crossings,
retain enough metric information to control area: they do not unless the zero
cell and its lengths are carried along.
