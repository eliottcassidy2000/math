# Reflection: LRC@14 via permutohedron geometry — honest assessment

**Session:** opus-2026-06-01-S526
**Date:** 2026-06-01

## What the permutohedron geometry reveals

The n=14 problem lives on the 13-torus T^13. The runner positions
(v₁t,...,v₁₃t) mod 1 trace a LINE. The lonely set L = [1/14, 13/14]^13
is a BOX of volume (6/7)^13 ≈ 13.5%. LRC says the line hits the box.

The permutohedron Π₁₃ decomposes the torus into cells (one per circular
ordering). The LRC walk moves through cells by adjacent transpositions.
Within each cell, the lonely condition is LINEAR.

## What DID work

1. **All 13 runners are essential** — removing any single runner opens a gap.
   No redundancy in the slab covering.

2. **The CRT factoring is powerful.** At δ = {7t} near 0, each CRT class
   constraint has width ~12/14 and the product is (6/7)^6 ≈ 37.6% — a huge
   target. The line in the 6-torus (θ₁,...,θ₆) = (t,2t,...,6t) must hit
   a 37.6%-volume box. For most speed sets this is easy.

3. **Equidistribution counts confirm LRC exhaustively.**
   - Initial segment: 6/360,360 lonely lattice points (the φ(14)=6 wall times)
   - Primes: 210/2310 ≈ 9.1%
   - Spread: 5,258/51,051 ≈ 10.3%

4. **The initial segment is the extreme case.** It has only 6 lonely points
   out of 360,360 — ratio 0.0017%. All other tested speed sets have
   MUCH more generous lonely measures.

## What DID NOT work

**The Bonferroni claim was WRONG.** The script claimed "Bonferroni proves
μ(lonely) ≥ 0.906" — this is an error. The second-order Bonferroni gives
a LOWER bound on μ(∪ slabs), not an upper bound. The correct statement:

μ(∪ slabs) ≥ 13/7 - Σ pairwise = 0.094

This says the union covers at LEAST 9.4%, but we need it to cover LESS
than 100%. The Bonferroni doesn't prove this.

**MISTAKE LOGGED.** Higher-order inclusion-exclusion (triple overlaps, etc.)
is needed to get an UPPER bound on the union. This is computationally
intensive and the convergence is slow.

## The proof gap — stated precisely

LRC@14 is equivalent to: for every primitive (v₁,...,v₁₃) with gcd=1,
the Weyl sum S(N) = (1/N) Σ_{k=0}^{N-1} 1_L(v₁k/N,...,v₁₃k/N) > 0
where N = lcm(v₁,...,v₁₃).

The equidistribution theorem gives S(N) → (6/7)^13 ≈ 13.5% as
max(v_i) → ∞. So for LARGE speeds, LRC@14 is easy.

The hard case is SMALL speeds — specifically the initial segment {1,...,13}
where S(N) = 6/360360 ≈ 0.0017%. The 6 lonely times are the φ(14)=6
regular polygon configurations t = k/14 for gcd(k,14)=1.

## What the permutohedron adds

The face lattice of Π₁₃ encodes WHICH orderings are adjacent. The LRC
walk visits orderings in a specific sequence. The lonely orderings (where
both observer-adjacent runners are far enough) form a specific set of faces.

The key question: is the SET OF LONELY FACES always on the walk's path?

For the initial segment: the walk passes through the IDENTITY ORDERING
(runners in speed order, all equally spaced) at t=k/14. This ordering
is the CENTER of the lonely face. The walk hits it exactly 6 times.

For general speed sets: the walk visits DIFFERENT orderings. The
permutohedron guarantees that every ordering is reachable from every
other in at most O(n) steps. But the LRC walk is DIRECTIONAL (THM-387),
not free — so reachability needs the directional constraint.

## The honest bottom line

LRC@14 is NOT proved by any computation in this session. The proof
requires either:
A. A rigorous upper bound on μ(∪ slabs) via higher-order inclusion-exclusion
B. A Weyl sum bound proving S(N) > 0 for all primitive speed sets
C. A permutohedron argument proving the walk must cross a lonely face
D. A CRT-factored argument proving the 7-class coupon collector completes

All four approaches are within reach but require more work. The
computational evidence is overwhelming (498/498 at n=14, all speed sets
lonely), but a formal proof remains open.
