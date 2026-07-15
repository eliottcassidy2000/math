---
id: THM-855  # renumbered from THM-810 (collision with the LRC four-replacement coset interface, first-pushed 21:04 vs 23:57 on 2026-07-14; opus-S312 hygiene)
title: The wiggly dictionary formulas — degree = m·fiber; the stationary law (H/|Aut| is the equilibrium measure of the metagraph walk); the per-flip drop law Δx = 4(d_loser − d_winner) + 8 (upsets climb, expected results descend, margin-2 is level); the unit-current cut law; and the level-mean status after the S308 decision
status: PROVED (F1–F4, each with a one-line proof; referee-verified exhaustively n = 4..6) + the level-mean law status recorded per the n=8 verdict
source: opus-2026-07-14-S309 (owner directive: prove the level-mean law, more node/edge/tiling statements, small formulas everywhere)
depends_on:
  - THM-787/790/793   # the flow, the leg law, the tower
related: [THM-801 (codex: the A+B+C−D−E−F+G Čech descent), THM-588, HYP-6870, HYP-6875]
verification: 05-knowledge/results/small_formulas_referee_opus_S309.out (exhaustive n = 4..6)
---

# THM-855 — the wiggly dictionary formulas

> (Formerly numbered THM-810; renumbered by opus-S312 — the LRC THM-810 was pushed first.)

Four exact formulas, each a one-line proof, each a tournament ↔ metagraph
dictionary entry (per the owner's principle: every metagraph fact has a
tournament reading, and vice versa).

## F1 — the degree law

> **The weighted wiggly degree of a class (self-loops doubled) is exactly
> m·fiber(class) = m·H/|Aut|.** *Proof:* every tiling has exactly m flips;
> project the fiber. ∎
> Tournament reading: a class's total connectivity in the metagraph is its
> Hamiltonian-path count, rescaled by its symmetry.

## F2 — the stationary law (the equilibrium meaning of Rédei's count)

> **The wiggly random walk on classes has stationary measure
> π(class) = fiber/2^m = H/(|Aut|·2^m).** *Proof:* the tiling-level walk on
> Q_m is regular, so uniform is stationary; class-lumping sends uniform to
> fiber-proportional, and F1 makes the class walk reversible w.r.t. it. ∎
> Tournament reading: **H per automorphism IS the long-run frequency with
> which random arc-perturbation visits an isomorphism class** — Rédei's
> quantity is the metagraph's equilibrium law. (Rédei's theorem, H odd, says
> every class has POSITIVE stationary mass of a definite parity; the
> transitive's H = 1 is the minimal possible equilibrium weight.)

## F3 — the per-flip drop law (the upset dictionary)

> Flipping one arc whose current winner w beats loser l changes the axis by
> **Δx = 4(d_l − d_w) + 8** (centered scores d = 2s − (n−1)). Hence:
> - **flipping an UPSET (d_w ≤ d_l) strictly climbs toward the transitive**
>   (Δx ≥ +8);
> - **flipping an expected result with margin d_w − d_l > 2 strictly descends
>   toward the circulant**;
> - **margin exactly 2 ⟹ a LEVEL flip (Δx = 0)** — the level edges of the
>   metagraph are precisely the flips of minimum-nontrivial-margin results.
> *Proof:* d_w ↦ d_w − 2, d_l ↦ d_l + 2 in x = Σd². ∎
> This refines THM-801's high/gap/low tri-partition with an exact per-edge
> current, and gives the wiggly layer's Δx spectrum {4k + 8 : k = d_l − d_w}.

## F4 — the unit-current cut law (the flow of transitivity, exactly)

> In the Smith network of G_n (unit current, transitive → distributed rail),
> **the net current through EVERY axis-threshold cut equals exactly 1.**
> *Proof:* Kirchhoff conservation: every threshold cut separates source from
> sink. ∎ This is the proved skeleton of the flow picture: whatever the local
> discordances (S308), one full unit of transitivity crosses every level
> boundary, at every n.

## F5 — the margin identity and the exact drift law (PROVED; sampled-exact to n = 20)

> **Σ over arcs of (s_winner − s_loser) = x/2** for every tournament — the
> total score margin is half the axis. Hence the uniform arc-flip walk has the
> EXACT Ornstein–Uhlenbeck drift
> **E[Δx | one uniform arc flip] = 8 − 4x/C(n,2)**,
> whose equilibrium x* = n(n−1) is exactly the stationary mean of x under
> uniform tournaments (E[Σd²] = 4n·Var(Binomial(n−1,½)) = n(n−1)) — uniform is
> arc-flip-stationary, and the drift vanishes precisely at its mean.
> *Proof:* Σ_arcs s_w = Σ_v s_v² and Σ_arcs s_l = (n−1)Σs − Σs²; subtract and
> use x = 4Σs² − n(n−1)². ∎ (The wiggly version carries a bounded base-path
> correction p(t) = Σ_path(d_w − d_l): E[Δx] = 8 − 4(x − p(t))/m.)
> Tournament reading: the mean-reversion of score variance under random arc
> perturbation is EXACT and linear — the mechanism behind the axis's
> mean-field behaviour (F4, the level-mean law below).

## F6 — the second-moment closure and the exact fluctuation–dissipation law (opus-S312; referee-verified exhaustively n = 4..6)

> With N = C(n,2), for EVERY tournament T:
> **E[Δ² | one uniform arc flip] = 16(n−4)·x/N + 64** — exact, pointwise in T.
> *Proof:* Δ = 4(d_l − d_w) + 8; Σ_arcs (d_w − d_l) = x (F5, doubled to
> d-units); and the orientation-free quadratic identity
> **Σ_pairs (d_u − d_v)² = n·x** (expand; Σd_v = 0). ∎
> Consequences (all exact):
> - **(E[x_t], E[x_t²]) evolve autonomously**: the moment hierarchy closes at
>   order 2. E[x_t] − 2N = (1 − 4/N)^t·(E[x_0] − 2N): explicit relaxation
>   rate 4/N.
> - **Stationary variance Var(x) = 2n(n−1)(n−2) = 4N(n−2)** (48, 120, 240 at
>   n = 4, 5, 6 — matches direct uniform-tournament computation, where
>   Cov(d_u², d_v²) = 0 exactly).
> - **The fluctuation–dissipation relation holds EXACTLY, not asymptotically**:
>   Var = E_stat[Δ²]·N/8 = (32(n−4) + 64)·N/8 = 4N(n−2). The axis walk is a
>   bona fide integer OU process.
> - **The closure boundary is sharp (F6d):** E[Δ³ | T] is NOT a function of x
>   — it splits already at n = 5 (x = 8, 16, 32) via the new odd invariant
>   y₃ = Σ_arcs (d_w − d_l)³. Even powers reduce (orientation-free); odd
>   powers generate new invariants y₅, y₇, ...
>
> **LRC transport (the seven-comb wall).** The collar overlap energy
> K = Σ_j C(d_j−1, 2) of the LRC14 frontier carries step-eight increments
> (THM-787) and the endpoint-defect flux Δ(−K) = d_dep − d_entry − 1
> (THM-785) — the same drift arithmetic as F3/F5. The seven-comb wall is a
> FIRST-moment stall (mean danger density 14/13 ≥ 1): F6's lesson is that the
> walk still has exact SECOND-moment structure. Candidate replacement
> potentials: Φ = (K − K*)² or degree-weighted energies Σ f(d_j) with f
> convex, whose event-drift E[ΔΦ] can be computed by the same
> size-biased-sampling algebra once THM-815's event-sampling law is fixed.
> The metagraph proof gap (drift → level-mean monotonicity, below) and the
> comb wall's replacement potential are THE SAME PROBLEM SHAPE: converting an
> exact linear drift into monotone hitting statistics.

## The level-mean law (status per the S308/S309 verdicts)

Exact concordance is FALSE (n=7: 132 real adjacent overlaps; n=8: spreading
to distance 32 — S308's decision). **Level-MEAN and level-MEDIAN monotonicity
hold at n = 5, 6, 7 AND 8** (the n=8 next-scale check PASSED even as the
pairwise overlap deepened to 5.1 mean-gaps) — the mean law survives where the
pointwise law died. Proved content: F4 (unit current through every cut), F5
(the exact mean-reverting drift — the mechanism), and the max-principle
endpoints. **The level-mean monotonicity conjecture** (now n ≤ 8 verified,
drift-mechanism named): the harmonic potential's level averages descend
strictly along the axis at every n. The proof gap: converting the walk's
exact drift into hitting-probability level-averages — a potential-theory
lemma over the F5 drift, the sharpest named target this thread has produced.
