---
id: THM-2122
title: "Multiple terminal 13-blockers force a guard-kernel orbit escape"
status: >
  PROVED. Let a rank-two torus have one guard and eight terminal characters.
  Assume the guard is nonzero modulo 13, at least two terminals are divisible
  by 13 in the saturated character lattice, and every remaining terminal is
  transverse to the guard modulo 13. Then the terminal danger bands cannot
  cover the guard-safe region even almost everywhere. Settled LRC for at most
  eight speeds produces an open quotient phase where all divided blockers are
  safe. On a guard-safe thirteenth-root sheet, the guard-kernel orbit has 13
  points, while each of the at most six residual terminals hits at most two;
  hence 2(8-b)<13. A finite-sheet pigeonhole preserves positive measure. No
  independence among the blocker directions is required.
source: codex-2026-07-22-LRC-multiple-thirteen-blocker-capacity
depends_on:
  - THM-2114
related:
  - THM-2116
  - THM-2120
external: Settled Lonely Runner Conjecture for at most eight integer speeds (LRC(9) in the total-runner convention).
---

# THM-2122 -- multiple 13-blockers force orbit escape

Let `Gamma` be a saturated rank-two character lattice and

```text
K=Hom(Gamma,R/Z).                                      (1)
```

Take a nonzero guard `g in Gamma` and eight nonzero terminal characters
`c_1,...,c_8`. Put

```text
C={X in K:||g.X||>1/7}.                               (2)
```

Assume:

```text
g mod 13 is nonzero in Gamma/13Gamma;                  (3)

B={i:c_i in 13Gamma},             b=|B|>=2;            (4)

c_i mod 13 is not proportional to g mod 13
                                      for every i notin B.  (5)
```

We prove that

```text
C subset union_(i=1)^8 {X:||c_i.X||<1/14}             (6)
```

is impossible even after allowing a null exceptional set.

## 1. Make every divided blocker safe at once

For `i in B`, write

```text
c_i=13u_i,                    u_i in Gamma\{0}.        (7)
```

Discard duplicate sign classes among the `u_i`. Choose a cocharacter
`d in Hom(R/Z,K)` outside the finitely many rational hyperplanes

```text
u_i.d=0,              (u_i-u_j).d=0,              (u_i+u_j).d=0             (8)
```

for distinct surviving sign classes. Such an integral `d` exists because a
finite union of proper rational lines cannot exhaust the rank-two
cocharacter lattice. The positive integers

```text
s_i=|u_i.d|                                             (9)
```

are nonzero and distinct.

There are at most `b<=8` of them. The settled Lonely Runner Conjecture for at
most eight integer speeds supplies `t in R/Z` such that

```text
||s_i t||>=1/(number of sign classes+1)>=1/9>1/14.    (10)
```

Put `Y_0=t d in K`. Equation (10), including the discarded duplicates, says

```text
||u_i.Y_0||>1/14                    for every i in B.  (11)
```

The inequalities have uniform positive margin. Hence there is a nonempty
open neighborhood `U` of `Y_0` on which every divided blocker remains safe.

This is the only use of the settled small-runner theorem. It is applied to
the blocker **directions**, not to the original terminals or an arbitrary
shifted LRC instance.

## 2. Choose a guard-safe thirteenth-root sheet

Multiplication by thirteen,

```text
[13]:K->K,                   X |->13X,                 (12)
```

is a finite covering with kernel `K[13]` of order `13^2`. Fix a root
`X'_0` of `Y_0`. As `z` ranges over `K[13]`, condition (3) makes `g.z` range
surjectively over the thirteen-torsion grid, each value thirteen times.

Any translate of the thirteen-grid has a point of circle distance greater
than `1/7`: its open radius-`1/7` danger arc has length `2/7<4/13` and hence
contains at most four grid points. Choose `z_0` with

```text
||g.(X'_0+z_0)||>1/7.                                 (13)
```

The covering (12) has a local inverse sheet `sigma` through `X'_0+z_0`.
After shrinking `U`, we may assume

```text
13 sigma(Y)=Y,                 ||g.sigma(Y)||>1/7      (14)
```

for every `Y in U`.

Let `v` be a nonzero element of the one-dimensional kernel

```text
L=ker(g:K[13]->(R/Z)[13]).                             (15)
```

Then the thirteen local root sheets

```text
X_j(Y)=sigma(Y)+jv,                   j in F_13,       (16)
```

all satisfy `13X_j(Y)=Y`, and their guard value is the same strictly safe
value from (14).

For every blocker `i in B`, equations (7), (12), and (16) give

```text
c_i.X_j(Y)=13u_i.X_j(Y)=u_i.Y.                         (17)
```

Thus every blocker is strictly safe on every sheet (16), for every `Y in U`.

## 3. Residual orbit capacity is twelve, not thirteen

Fix `i notin B`. Over `F_13`, the annihilator of the line `L` in (15) is the
line spanned by `g mod 13`. Assumption (5) therefore gives

```text
c_i.v !=0 in (R/Z)[13].                                (18)
```

Consequently, as `j` ranges over `F_13`, the values

```text
c_i.X_j(Y)=c_i.sigma(Y)+j(c_i.v)                       (19)
```

form a complete translate of the thirteen-grid. A radius-`1/14` open danger
arc has length `1/7<2/13`, so it contains at most two of these values. Hence

```text
#{j:||c_i.X_j(Y)||<1/14}<=2.                           (20)
```

There are `8-b<=6` residual terminals. The union of all their dangerous
indices on the orbit (16) therefore has size at most

```text
2(8-b)<=12<13.                                         (21)
```

For every `Y in U`, at least one sheet index `j` is residual-safe. By (14)
and (17), that same point is also guard-safe and blocker-safe.

## 4. The finite-sheet null-set sidecar

The index selected after (21) may depend on `Y`; this causes no loss. Define

```text
A_j={Y in U: ||c_i.X_j(Y)||>=1/14 for every i notin B}. (22)
```

These thirteen measurable sets cover the positive-measure open set `U`.
Therefore some `A_j` has positive Haar measure. The local sheet map
`Y |-> X_j(Y)` is a diffeomorphism onto its image and preserves Haar-null
sets up to a constant Jacobian. Its image of `A_j` has positive measure and,
by (14), (17), and (22), lies in `C` but outside every open terminal-danger
band. Removing the finitely many residual endpoint hypersurfaces leaves a
positive-measure set on which all terminal inequalities are strict.

Thus the uncovered subset of `C` has positive measure, contradicting (6)
even in its almost-everywhere form. QED.

## 5. Frontier effect

THM-2114 forces at least one guard/terminal 13-content blocker in every
rank-eight cover. Under `g mod 13!=0` and residual transversality:

```text
exactly one terminal blocker, independent of g   -> empty by THM-2120;
at least two terminal blockers                    -> empty by THM-2122.      (23)
```

The blockers in the second line may be mutually dependent, mutually
independent, repeated up to sign, or individually proportional to the guard;
none of that affects the proof. THM-2123 later closes the unique dependent
blocker as well. The live rank-eight branches are now a guard 13-blocker or a
residual terminal nontransverse to the guard modulo 13. Ranks nine through
twelve remain outside the `2(8-b)<13` count as stated.

## 6. Assumption challenge and Tournament Analysis

The challenged assumption was that several 13-blockers make the finite orbit
harder because each blocker covers all thirteen points when active. Dividing
all blockers first reverses the burden: settled small-runner LRC supplies an
open phase on which **none** is active, and every additional blocker removes a
residual two-point toothpick from the orbit ledger.

Candidate tournament vertices were terminals, blocker sign classes,
thirteenth-root sheets, guard-kernel orbit points, residual toothpicks, and
coverage obligations. The faithful finite vertices in Section 3 are the
thirteen root sheets. Orienting sheets by the first residual label that hits
one and using cyclic `F_13` order as the tie Hamiltonian path records a search
schedule, but destroys the uniform `<=2` label capacities and the
positive-measure sets `A_j`. The preserved carrier is the bipartite
sheet--residual incidence graph together with the quotient-safe open set and
local-sheet null-set sidecar; no tournament quotient alone proves (21)--(22).
