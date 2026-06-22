# LRC Realizability Is Scaling-Invariant Exact-Tiling, Not Combinatorial Forcing

*mac-mini-2026-06-22-S39. The owner asked for realizability arguments "similar to tournament analysis,
but maybe the LRC's missing structure is slightly different." It is — and pinning exactly how sharpened
the finish. Companion to HYP-2888, kps HYP-2885.*

## The shared shape: an obstruction that cannot be realized

**Tournament.** H(T) = I(Ω(T),2). H=7 is forbidden not by counting but by REALIZABILITY: the conflict
graph Ω cannot equal K_3, because three pairwise-vertex-sharing triangles force a common vertex → a
fifth vertex carrying a directed C_5 → a fourth odd cycle (kps THM-200). The forbidden value is the
shadow of a non-realizable combinatorial configuration. The forcing is **discrete and combinatorial**.

**LRC(14).** Restated as coverability: LRC ⟺ the safe set {t : ‖st‖ ≥ 1/14 ∀s ∈ S} is non-empty ⟺
the open arc-systems U_s = {‖st‖ < 1/14} (s arcs of width 1/(7s), total 13/7 ≈ 1.857) do not cover
the closed circle. A counterexample is a 13-set whose arcs OVER-cover (kill even the boundary). The
LRC asks: is over-covering realizable? Same shape as Ω=K_3 — an obstruction we want to show unrealizable.

## How the LRC's structure is "slightly different"

The forcing is **continuous, arithmetic, and scaling-invariant** — not combinatorial. Verified (S39):

- **Exact coverage (meas safe = 0) is achieved ONLY by the consecutive-multiples** S = d·{1,…,13}.
  These tile [0,1) by their arcs to measure exactly 1. (0 of 5668 random/perturbed sets reached it;
  max non-multiple was 0.998.)
- **Additive energy is the WRONG invariant.** A(E) is translation-invariant — every length-13 AP has
  A = 1469 (the max). But coverage is SCALING-invariant: {2,…,14} and {3,…,15} have A = 1469 yet
  positive safe measure (0.061, 0.098). So coverage ≠ additive energy; they have opposite symmetries.
  The exact-coverage extremal is a strictly smaller, scaling-invariant class than the max-A APs.
- **The tight case is explicitly safe.** Each d·{1,…,13} has the witness t = 1/(14d):
  ‖jd/(14d)‖ = ‖j/14‖ = min(j,14−j)/14 ≥ 1/14 for j = 1..13. So the unique tight family carries a
  closed-form safe boundary point — no measure, no limit needed.

## The finishing structure (and the honest crux)

LRC(14) ⟸ **(2)** every d·{1,…,13} is safe at t = 1/(14d) [RIGOROUS] **+ (1)** every 13-set that is
NOT d·{1,…,13} has positive safe measure [the open crux, 0/5668 counterexamples].

Crux (1) is now stated cleanly and scaling-invariantly: *the only exact tilers of the circle by the
arc-systems U_s are the consecutive-multiples.* This is, in essence, the uniqueness of the LRC
extremal — so the reframing **sharpens** the hard core rather than removing it, but it does two useful
things: it isolates the hardness into a tiling-rigidity statement (no measure/extremality LP), and it
disposes of the tight family in closed form. The complementary route (HYP-2876/2878) is
extremality-free: a bounded-denominator rational witness for every set, where the crux is instead the
over-determination of the covering system.

## The apex-7 thread, both sides

Tournament: H=7 = I(K_3,2), forbidden at the prime 7. LRC(14): threshold 1/14 = 1/(2·7), tight on the
consecutive-multiples whose boundary witness sits at t = 1/(14d) — the D=14 = 2·7 denominator. Both
obstructions are non-realizable structures pinned at the apex prime 7; the tournament's is combinatorial
forcing (Ω=K_3 ⇒ C_5), the LRC's is geometric tiling-rigidity (only {1,…,13}·d tiles at 1/14).

Related: [[seven-is-the-unique-forbidden-clique-value]], [[the-even-graph-is-the-tournaments-cycle-half]],
HYP-2888 (exact coverage), HYP-2885 (kps additive-energy route), HYP-2876 (rational-witness route).
