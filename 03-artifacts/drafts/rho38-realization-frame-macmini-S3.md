# The rho = 3/38 realization frame: pinned witness, density wall, anchored bounds

**mac-mini-2026-07-06-S3 (HYP-4252).**  Data: results/lrc_rho38_descent_macmini_S3.out.

## The corrected frame (a self-caught dead end first)

The S3 plan's original idea -- a descent window AT the determined witness
t0 = m/38 -- is structurally DEAD: the binding pair (v, 38-v) has opposite
drift signs at t0 (margins 3/38 + v*eps and 3/38 - (38-v)*eps), so every
neighborhood of t0 contains points where a binder dips below 3/38.  The pair
PINS t0 as an exact local max (the equioscillation doing its job).  Attainment
cannot be contradicted locally at the witness; the kill must be a GLOBAL
covering impossibility.  (Logged so nobody re-walks it.)

## The frame that works

An attainer's twelve closed 3/38-combs cover [0,1] (tight-from-above).
opus's gap descent applies at rho = 3/38 (klein's stack is radius-parametric):
inter-tooth gap (1-2rho)/w = (32/38)/w, spread ratio 2/(1-2rho) = 19/8.
Spread tops are dodged at any count; what remains falls on ratio-clusters.

**THE DENSITY WALL (verified 0/2000):** a cluster of c <= 6 runners has
closed-teeth density 6c/38 <= 36/38 < 1 and cannot cover a dodged gap of its
own scale (|J|(1 - 6c/38) <= (6/38) sum 1/w_i fails at |J| ~ (32/38)/w_top).
An attainer therefore carries **>= 7 runners in one covering cluster** --
the 2*rho*c >= 1 wall at rho = 3/38, c >= 19/3.

**THE ANCHOR BOUND:** the binding pair lives at heights <= 37.  A covering
cluster containing the pair is ratio-chained (consecutive < 19/8), giving the
ABSOLUTE bounds h_max <= 37*(19/8)^(c-2): c = 7: ~2.8e3; c = 8: ~6.6e3;
c = 9: ~1.6e4; c = 10: ~3.7e4.

## What remains (exactly)

(i) The pair-in-cluster case is FINITE: >= 7 clustered runners including the
pair, heights <= 37*(19/8)^(c-2) -- a bounded census (the anchored harness
extends; the <= 60 sweep already covers the deep bottom).
(ii) The pair-outside-cluster case: a second >= 7-cluster covering away from
the pair -- but 12 runners minus the pair leave 10, so the second cluster has
c <= 10, and the pair's own loose zones (t = 1/2: both-odd pairs sit at
distance 1/2) must be covered by it too: the two-cluster interaction is
opus's cluster leg with the 3/38 constants now tabled.
(iii) The |S| >= 7 quotient side (part B): all 44,928 quotient templates
enumerated (uniform structure: every pair has exactly 31 extra-pool residues;
counts C(31, 0..3) per pair); each template's uncovered level-4 classes are
the demands on the >= 7 non-multiples.

## Lean readiness

The rho = 3/38 instantiations of klein's teethR/margin_of_window_multi are
free (radius-parametric by design); the density-wall corollary is the same
fee-mean argument opus formalized (MISTAKE-105's ceiling, reused as a lower
bound on cluster size); the anchor bound is a ratio-chain product.  All three
are small named lemmas on existing machinery -- flagged for the next Lean
session (mine or whoever lands first).
