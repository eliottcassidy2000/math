# The column denominator is an Euler-characteristic ceiling

*monad-explorer-2026-06-07 (deep-research / analytic lane, 9th session). Builds on THM-438
ADDENDUM-6 (column rationality conjectured) → ADDENDUM-7 (denominator PROVED).*

## The thing that wanted explaining

The 8th session found, empirically, that each cycle-rank column of the `(★★)` triangle is rational
with a very specific denominator:

```
   T_m(x) = Σ_k t(k,m) x^k = P_m(x) · x^m / (1−x)^{2m−1}.
```

Why `2m−1`? Not `m`, not `3m−3`, not `2m`. The number `2m−1` sat there as a fitted constant. A
fitted denominator is a fact; it is not yet understanding. This session it turned into a one-line
consequence of Euler's handshake — and the *reason it is `2m−1` and not the naive `3m−3`* is the
whole story.

## The ceiling

The pole order of a rational counting series is the **maximum number of independent geometric
parameters** in the objects being counted. Here the objects are even-series patterns of fixed cycle
rank `m`; the geometric parameters are the **lengths of the lines** (series classes). Each line can
be even-subdivided independently — one factor `x/(1−x)` apiece. So:

> pole order of `T_m` = max number of lines over rank-`m` even-series patterns.

Contract every line to an edge and you get the **core** — a connected multigraph of cycle rank `m`,
series-reduced (so min degree ≥ 3, bar the single-vertex rosette). Maximizing the number of edges
`e = V_core + m − 1` means maximizing the number of vertices `V_core`. Two ceilings press down on
it from opposite sides:

- **Series-reduction pushes vertices UP**: more branch vertices = more edges. A trivalent core
  (all degree 3) is the greediest: `V_core = 2m−2`, `e = 3m−3`.
- **The Eulerian-trail condition pushes vertices DOWN**: the pattern is a single open walk, so the
  core may have **at most two odd-degree vertices**. A trivalent core has *every* vertex odd —
  `2m−2` of them — and so admits no single trail. It is geometrically real but *physically
  unreachable* by one walk.

The handshake resolves the tension exactly. With at most two odd vertices and the rest of even
degree ≥ 4, `2e = Σ deg ≥ 4V_core − 2`, while `2e = 2V_core + 2m − 2`. Subtract: `V_core ≤ m`. So

```
   #lines = V_core + m − 1 ≤ 2m − 1,
```

attained by cores with exactly two degree-3 vertices (the walk's two ends) and `m−2` degree-4
vertices. **`2m−1` is the largest an Eulerian, rank-`m`, series-reduced skeleton can be.** The
denominator is an Euler-characteristic ceiling.

## Why this is the right kind of explanation

It is the same move that runs through the whole project, in a new costume. "Everything is the
triangle": invariants are controlled by a topological skeleton, and the constants that look
arithmetic are really combinatorial-topological. Here `2m−1` *looked* like it might be arithmetic
(it is the dimension of the cycle space of a maximal core, `e − V_core + 1 = m`, plus the
`V_core−1 = m−1` tree edges — `m + (m−1) = 2m−1`). It is. The pole order counts **all** the edges
of the maximal Eulerian core: `m` chords + `m−1` tree edges.

And the *exclusion* is as informative as the bound. The trivalent `3m−3` cores are exactly the ones
genus-blindness (ADDENDUM-5) and the loop equation already told us not to expect: they have too
many odd vertices to be a single trail. The Eulerian constraint — that there is **one** walk, not a
disjoint union — is the same "connectedness" that makes the cluster expansion linked, that makes
`S_k` a free *cumulant* rather than a moment. The pole order `2m−1` and the Catalan law are two
shadows of the same single-walk discipline.

## A correction that sharpens, not weakens

The tempting next step — "then `R_s(m,e)`, the coefficient of `(x/(1−x))^e`, must be the number of
reduced patterns with `e` lines" — is false. The brute census gives `13,14,4`; the actual
coefficients are `13,33,20`. They agree only on the diagonal. The gap is trail-ordering symmetry:
the lines of a core are visited in a definite order along the one walk, so permuting equal-length
lines is *not* a symmetry of the pattern, and the naive product formula mis-counts. This is the
fourth time on this single thread (MISTAKE-060, -061, -062, now this) that a clean final number was
reached by a per-class story that did not survive contact with the census. The pattern is itself the
lesson: **on this object, the totals are trustworthy and the local stories are not.** The signed sum
is Catalan; the unsigned count is uncatalogued; the column is rational; the per-line count is not
what it looks like. Structure lives in the cancellations and the poles — the global invariants —
never in the term-by-term picture. That is exactly what a free-probabilistic / topological-recursion
object should feel like.

## The next ceiling to climb

`deg P_m = m−2` is now one clean statement: `Σ_e (−1)^e R_s(m,e) = 0`, i.e. the line-count
alternating sum of (s-expansion) weights vanishes above rank 1 — `V(−1,y) = −y`. That is a
*within-column* sign-reversing involution flipping #lines by one. It is humbler than the headline
involution (handoff #2, shifting cycle rank `m` to prove `(★★)` itself), but it lives one level down
and may be where the real involution first becomes visible. When the denominator is a ceiling, the
numerator is what is left after the ceiling is subtracted — and that residue, degree `m−2` with
leading coefficient `2^m−1` (one per nonempty subset of the `m` cycles?), is the next thing asking
to be understood.
