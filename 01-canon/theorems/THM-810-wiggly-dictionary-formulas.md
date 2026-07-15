---
id: THM-810
title: The wiggly dictionary formulas — degree = m·fiber; the stationary law (H/|Aut| is the equilibrium measure of the metagraph walk); the per-flip drop law Δx = 4(d_loser − d_winner) + 8 (upsets climb, expected results descend, margin-2 is level); the unit-current cut law; and the level-mean status after the S308 decision
status: PROVED (F1–F4, each with a one-line proof; referee-verified exhaustively n = 4..6) + the level-mean law status recorded per the n=8 verdict
source: opus-2026-07-14-S309 (owner directive: prove the level-mean law, more node/edge/tiling statements, small formulas everywhere)
depends_on:
  - THM-787/790/793   # the flow, the leg law, the tower
related: [THM-801 (codex: the A+B+C−D−E−F+G Čech descent), THM-588, HYP-6870, HYP-6875]
verification: 05-knowledge/results/small_formulas_referee_opus_S309.out (exhaustive n = 4..6)
---

# THM-810 — the wiggly dictionary formulas

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

## The level-mean law (status per the S308/S309 verdicts)

Exact concordance is FALSE (n=7: 132 real adjacent overlaps; n=8: spreading
to distance 32 — S308's decision). Level-MEAN and level-MEDIAN monotonicity
hold at n = 5, 6, 7; the n = 8 check: see
05-knowledge/results/level_mean_n8_opus_S308.out (recorded in the session
log). The provable content extracted so far is F4 (cuts) plus the max-
principle endpoints; full mean-monotonicity, if the n=8 verdict sustains it,
is the named conjecture — with the S308 lesson applied: no proof attempt
until the next-scale check passes.
