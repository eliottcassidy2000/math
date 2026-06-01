---
source: oracle-2026-06-01-S523
status: result (computational; round=A000016 verified 5 terms, equivalences exhaustive m<=5)
tags: [LRC, tournaments, round-tournament, locally-transitive, A000016, A000568, necklace, THM-381, roots-of-unity]
---

# LRC(n) is a tournament question — and the realizable iso classes are a necklace (A000016), not all of A000568

**Prompt (user):** see how LRC at a particular n is really a question modeled by a
tournament of some size, and understand the structure of the set of its possible
isomorphism classes.

## 1. The tournament, and its size

This is exactly THM-381. For a speed system `S = {0, v_1, …, v_{n-1}}` (stationary
observer `0`, threshold `1/n`) define the **observer-source tournament** on the
`n` vertices `{0,…,n-1}`:
- observer–runner: `0 -> i  iff  ||v_i t|| >= 1/n`,  else `i -> 0`;
- runner–runner: the half-turn circular comparator (fixed Hamiltonian tie path).

**THM-381:** the observer is lonely at time `t` iff vertex `0` is a **source** of
`T_S(t)`. So **LRC(n) ⇔ for every primitive speed system, the clock-movie
`t ↦ T_S(t)` (a tournament on `n` vertices) visits a class where the marked
vertex is a source.** The source-marked targets are counted by `A000568(n-1)`
(once `0` is a source, the rest is a runner tournament on `n-1` vertices).

So the size is `n` (full movie) / `n-1` (the runner sub-tournament that the source
mark exposes). That much was known. The new content is *which* classes occur.

## 2. The realizable classes are ROUND tournaments = a necklace = A000016

The runner sub-tournament is the half-turn tournament of points on a circle, so
**every vertex's out-set is a contiguous clockwise ARC**. Tournaments with that
property are exactly the **round** (equivalently **locally transitive**)
tournaments. Generically (open times) the movie can realize *only* these.

Computed (`lrc_realizable_isoclasses_s523.py`):

```
 m              3    4    5    6    7
 ROUND classes  2    2    4    6   10      <-- realizable by the open movie
 A000568(m)     2    4   12   56  456      <-- all tournaments
 fraction      1.0  0.5 0.33 0.11 0.02
```
Verified **round == locally-transitive == half-turn-realizable** by exhaustive
enumeration for m=3,4,5. Cross-checks S512's independent note ("the open clock at
total n=5 sees only four of twelve classes" = round(5)=4).

**The sequence 2,2,4,6,10 is `A000016`**, with the closed form
```
   round(m) = A000016(m) = (1/2m) · Σ_{d | m, d odd} φ(d) · 2^{m/d}
```
(predicts round(8)=16, round(9)=30). `A000016` is a **necklace** count, and the
formula sums over **odd divisors only** — the same odd-cycle Burnside structure
that produces `A000568` itself (CLAUDE.md: `Fix(σ)=0` for even cycles). So:

> **The realizable LRC tournaments are the cyclic/necklace refinement of all
> tournaments.** A round tournament *is* a binary necklace — the cyclic pattern,
> around the circle, of which gaps the half-turn window covers. The complement
> `T ↔ T^op` is the necklace's reversal/antipodal flip. This is, literally, the
> "complement necklace" the twin-Goldbach thread (S521) was circling.

## 3. The structure has two layers, and the conjecture lives on the seam

- **OPEN layer (round, A000016):** what the movie sees at generic times. A
  *vanishing* fraction of `A000568` (2% at m=7). Its most symmetric element is the
  regular `m`-gon `R_m` = the `m`-th roots of unity (S522; `R_7`, H=175).
- **BOUNDARY layer (compactification):** at wall times two runners are antipodal
  (gap = 1/2) and the Hamiltonian tie path resolves the tie, producing classes
  that need **not** be round. This is why the S512 boundary count jumps (n=5: 4
  open → 11 with boundary) and why the S520 lonely-**menu** (1,2,6,6,≥12 for
  n=4..8) can exceed `round(n-1)` (e.g. menu 6 > round(5)=4 at n=6): the extra
  lonely classes are boundary/tie-resolved, the *tight* LRC cases.

So the iso-class structure of LRC(n) is a **two-storey object**:
```
   A000568(n-1)                      all tournaments
      ⊋  boundary-compactified set   (tie-resolved wall classes; the hard tail)
            ⊋  ROUND = A000016(n-1)  (open movie; a necklace; the generic body)
                  ⊋  lonely menu     (source-reachable; S520)  ∩ both layers
                        ∋  R_{n-1}   (regular polygon; roots of unity; S522)
```
The generic body is a tiny, exactly-counted necklace `A000016`; **all the
difficulty of LRC is pushed onto the boundary seam**, where non-round
(tie-resolved) classes appear and where the tight extremal speed systems live.
That matches the S520 "min-near-runner histogram {0: many, 1: 1}": almost every
set is lonely at an *open* (round) time; the rare hard set is lonely only at a
*boundary*.

## 4. Why this sharpens the program

- "LRC lives in `A000568(n-1)`" is a 10x overstatement of the generic picture:
  generically it lives in `A000016(n-1)` — a *necklace*, 2% the size. A proof
  should treat the round body (closed-form, cyclic) and the boundary seam
  separately; only the seam carries the obstruction.
- Both the ambient `A000568` and the realizable `A000016` are **odd-divisor
  Burnside** objects. The round ones are the rotationally-coherent (single-cycle)
  refinement — which is exactly the roots-of-unity / regular-polygon story of
  S522 made enumerative.
- Open question (→ HYP-1998): is the *boundary-compactified* realizable set also a
  named subsequence between `A000016` and `A000568`? Count it for m=4..7
  (S512 gives 11 at m=5); identify the tie-resolution map round ⊕ Hamiltonian-path
  → extra classes, and whether the lonely menu is exactly
  `round(n-1) ∪ (boundary source classes)`.

## Anchor
`04-computation/lrc_realizable_isoclasses_s523.py` (+ `.out`): round=2,2,4,6,10
=A000016; round==loc-transitive==half-turn (m<=5). Builds on THM-381, HYP-1987
(S512/S520), HYP-1995 (S522).
