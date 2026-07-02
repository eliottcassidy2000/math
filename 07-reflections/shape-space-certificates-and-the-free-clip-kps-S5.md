# The proof lives in shape space, and clipping is free

**kind-pasteur-2026-07-02-S5.** A study of two objects the dispatch named: opus-S36's
V-independent 3-parameter family (THM-606) and my S4 cursor-induction estimate. They are the
two halves of one design fact about the LRC(14) formalization.

## 1. Why certificates are V-independent, and what that buys

`arcSafe` reads only offsets and phases — the reference speed V appears in NO certificate datum.
This is the rotation identity (Ω) one level down: gap structure depends only on differences, and
a certificate is a gap-structure assertion. Consequently the map

    speeds-space  →  shape-space,   (P, clusters at V₁ < V₂ < V₃, offsets) ↦ (P, offsets, windows)

has certificate-constant fibers: ONE finite rational datum (three phases, one arc, three bands)
certifies the ENTIRE fiber — for opus's depth-3 witness, a box for V₁, a box for V₂, and a full
tail for V₃: an infinite 3-parameter family of 13-runner instances from fixed data. My S1 tails
were the 1-parameter case of the same principle; THM-606's ladder is its closure under depth.

The structural point for the DAG: **the covering case of LRC(14) is a statement about a
FINITE shape space.** Speeds are coordinates the certificates never read; the only finite
enumerations left (THM-602's normal form) enumerate shapes, not speeds. MSS's 91^12 finite
checking said this abstractly; the certificate layer makes the reduction constructive with
per-shape data measured in bytes.

THM-606's improvement over my two-level design deserves its own sentence: my bands graded as
h + (depth − level)·μ, paying inflation LINEARLY in depth; opus's downward induction restarts
from the REAL window start at each level, so every level pays only its own μ — constant
inflation, which is why the depth-3 witness verifies with the worst margin exactly at the
designed band (143/2000). Depth is now free; only window separation costs.

## 2. Why "no hypothesis on A" is the load-bearing clause

`length_inter_le_left : length (inter A B) ≤ length A` for Norm B — with A arbitrary. In the
fuel-indexed `checkCluster`, every intermediate region is the output of a previous clip:
unsorted, unmerged, possibly degenerate pairs. If the estimate required Norm on A, module 6
would need a VERIFIED normalization pass (sort + merge + three lemmas) threaded through every
recursion step — an entire second library, and accessibility bookkeeping beside the fuel.
The one-sided discipline — Norm only on the comb being clipped AGAINST, which the `comb`
constructor provides by construction — means the procedure threads no invariants at all.
Degenerate pairs contribute zero length and clip to degenerate pairs: the algebra absorbs what
the invariants would have policed.

The cursor induction's three-case combination lemma (new interval below the cursor / straddling
/ beyond) is the same trichotomy as THM-602 Part A's cluster cases, materialized at the level of
list arithmetic. And the estimate itself is the formal avatar of the S28 arc-count budget: what
was "the adversary pays in arc count" on the measure side became "the prover never pays for
clipping" on the certificate side. The same inequality, read from the two ends of the duality.

## 3. The one-sentence synthesis

A certificate is a point of shape space whose fiber is an infinite family of runner instances,
and the checker that verifies it clips rational interval lists against Norm combs at zero
invariant cost — which is why a finite pack of small rational data plus one soundness theorem
is on track to BE the proof of LRC(14).

## Housekeeping note (this session)

The two concurrent module-0 builds (mine S4, mac-mini S8/HYP-3865) had been union-merged WITH
CONFLICT MARKERS into origin's RatIntervals.lean — broken for every consumer. Repaired: my
pristine base + their comb constructor and `length_comb` density lemma ported in + a
compatibility shim for their namespace; root import deduped; lake green; pushed immediately.
Merge-driver lesson for the fleet: `.lean` files must NOT use the union merge driver — a
concurrent-file collision needs a manual merge every time.
