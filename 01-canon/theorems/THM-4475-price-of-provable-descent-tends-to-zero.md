---
id: THM-4475
title: "The price of provable descent tends to zero: explicit L-step-descent trees in the Collatz pairing family at flip density <= 2 rho_L <= 2^(1-(1-h)L), h = h(log_3 2) (HYP-9136 proved), with a lower bound 2^(-0.774 L)"
status: >
  PROVED + INDEPENDENTLY AUDITED.
  In the pairing family of THM-4470 (a bit per pair {2i-1, 2i} chooses which
  member goes up; Collatz is the all-zero word), let P_L be the members in
  which every n >= 3 falls below itself within L steps.
  (A) For every L >= 8 there is an explicit member G_L of P_L. It is a tree
  whose only cycle is {1,2}, and its flip density is at most 2 rho_L, where
  rho_L = |Bad_L|/2^L <= 2^(-(1-h)L) is the density of classes on which
  Collatz is undecided, and 1 - h = 0.050044. Hence
  delta_L <= 2^(1-(1-h)L) -> 0: HYP-9136 holds, at the conjectured rate as
  an upper bound.
  (C) Every member of P_L has flip density at least rho_L / W_L, where
  W_L ~ c lambda^L and lambda = (3 + sqrt 13)/4. So
  delta_L >= 2^(-(0.7737+o(1))L).
  UPDATE 2026-09-26: THM-4478 proves the sharp lower bound
  2^(-(1-h)L-O(sqrt(L log L))); HYP-9137 is PROVED. The construction here
  supplies its matching upper bound. Collatz remains OPEN.
  SHEET: the construction transfers to the 3n-1 pairing, giving provable
  trees that destroy all three 3n-1 cycles.
source: collatz-procgen-20260922 session, price lane (2026-09-25), proving the session's HYP-9136 (candidate P1 of the brackets lane) with the coordinator's conjectured exponent; audited and promoted by the session orchestrator 2026-09-25
depends_on:
  - 01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md
related:
  - 05-knowledge/hypotheses/HYP-9136-provable-pairing-price-tends-to-zero.md (now SETTLED)
  - 05-knowledge/hypotheses/HYP-9137-sharp-provability-price-exponent.md
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md (the periodic analogue, HYP-9138)
  - 05-knowledge/results/collatz_procgen_20260922_choice_ladder.md (|Bad_L| = 2^(hL+O(log L)))
script: 04-computation/experiments/procgen_price_20260925_greedy.c
script_audit: 04-computation/experiments/procgen_price_20260925_orchestrator_check.py
output: 05-knowledge/results/procgen_price_20260925.out
output_audit: 05-knowledge/results/procgen_price_20260925_orchestrator_check.out
script_sha256: 3bae6f46160b8ce6872e35085070caa461f9ed220d7ba488dbc029fed6f47b50
script_audit_sha256: 896239b3251b63818a1ef1a6deba9e068f3882b3b55f7548caf99ffa59f23705
output_sha256: 98420c6bc474e1c564b2443edd783f18e7b4d11929a94f282da354be07b2add8
output_audit_sha256: ac0a8821c563f1fc3f726d697816b2f2a5cc06e2f998536b6dab5bd835b4c2be
hash_basis: raw LF bytes
audit: >
  The orchestrator read the proof (G1-G9 of the lane note) line by line
  and checked:
  * G3 (partner isolation via the multiple-of-3 entry argument);
  * G5 (freeness of pair(T(n)));
  * G7's forced parity steps, e.g. (27n+23)/32 < n and
    (243n+319)/256 < n for n >= 25, giving n = 507 (mod 512);
  * G8's late-visit case, including the endpoint formulas (81s+65)/16,
    (81s+38)/16, (243s+211)/32 and (243s+130)/32, and the
    non-divisibilities 81 ∤ 59, 81 ∤ 167, 243 ∤ 103, 243 ∤ 59;
  * G9's partner bound y_9 = (729n+1085)/512 < (3n+3)/2 for n >= 9;
  * the lower-bound recurrence g_k <= 1.5 g_(k-1) + 0.25 g_(k-2), giving
    lambda = (3+sqrt13)/4 and exponent 0.7737.
  An independent re-implementation of G_L, written from the note's prose
  only (procgen_price_20260925_orchestrator_check.py), finds, for
  L = 8, 12, 16 up to 10^6:
  * no stuck n;
  * every n <= 5*10^5 descends within L;
  * flip densities 0.06806 and 0.02850 at L = 8 and 16, equal to the lane's
    values;
  * B-rescues only at n = 507 (mod 512).
  The lane's pipeline was re-run: byte-identical output.
---

# THM-4475 -- the price of provable descent tends to zero

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_price_20260925_provability_price](../../05-knowledge/results/procgen_price_20260925_provability_price.md).

## 1. Setting

* **The pairing family (THM-4470).** Pairs are `{2i-1, 2i}` with one bit `eps_i` each. A number `n` in pair `i` goes up by `i` iff (`n` odd) XOR `eps_i`. Flipping pair `i` sends `2i-1` down to `i-1` and its partner `2i` up to `3i`. Collatz is `eps = 0`.
* **Provable members.** `P_L` is the set of members in which every `n >= 3` falls below itself within `L` steps. By strong induction such members have every orbit entering a finite set. `delta_L` is the infimum over `P_L` of the flip density.
* **The undecided density.** `rho_L = |Bad_L|/2^L` is the density of classes mod `2^L` on which Collatz has no `L`-step descent. `rho_L <= 2^(-(1-h)L)` (Chernoff), `rho_L >= 2^(-(1-h)L)/poly(L)` (cycle lemma), with `h = h(log_3 2) = 0.949956`.

## 2. The construction

* **Processing.** Process `n = 3, 4, 5, ...` in order. The default path of `n` uses frozen bits and reads free bits as 0 (Collatz). If it descends within `L` steps, freeze the pairs it visits.
* **Rescue.** Otherwise flip one free pair and freeze the new certificate. Take the first option that works, in the order A, F, B when `n = 3 (mod 8)`, and F, B, A otherwise:
  * **A:** `n`'s own pair;
  * **F:** the pair of the first point `w < 2n` of the path that is the image of an up-move and satisfies `w = 3 (mod 4)`;
  * **B:** the pair of `T(n)`.

## 3. Why it works

* **Partner isolation.** An F/B flip point `v = T(p)` is `2 (mod 3)`. So its partner `u = v+1` is a multiple of 3, entered only from `2u` or from the partner of `p`, whose pair is frozen at 0. Partners are therefore visited before descent only by themselves.
* **Rescues exist.** Every rescued `n` is odd and `3 (mod 4)` and lies in `Bad_L`. Its pair `pair(T(n))` is always free, so B is always available.
* **Where B is needed.** B is needed only at `n = 507 (mod 512)`, the 2-adic neighbourhood of the cycle `-5 -> -7 -> -10`. This follows from the forced parity word `110110110` and a late-visit exclusion by divisibility.
* **Partners come down.** Partners of A/F flips descend in 2 steps. Partners of B flips descend within 8 steps, re-entering the owner's own orbit at `y_7`.
* **Tree.** `eps_1 = 0`, so the only cycle is `{1,2}`.
* **Density.** Each flip is owned by a rescued `n` in `Bad_L` and has pair index in `(n/2, n]`. So the flip density is at most `2 rho_L`. ∎

## 4. Lower bound

**Current bound, 2026-09-26.** [THM-4478](THM-4478-critical-tube-affine-capacity-sharp-provability-price.md)
proves `delta_L >= 2^(-(1-h)L-O(sqrt(L log L)))` by a critical growth band
and actual-integer capacity cut. Thus HYP-9137 is proved. The bounds below
remain valid historical estimates, and the construction above is unchanged.

**Update 2026-09-25.** THM-4477 improves the lower bound to `delta_L >= rho_L^2/M2(L) >= 2^(-(0.1445+o(1))L)` (Cauchy–Schwarz with a certified second-moment majorant). It also shows that `2(1-h) = 0.1001` is the limit of any distribution-only argument.


* If `n in Bad_L`, its first `L` points must meet a flipped pair.
* Weighting by `1/n` and bounding how many bad orbits a flip at `v` can serve by the backward-tree weight `W_L`, which satisfies `g_k <= 1.5 g_(k-1) + 0.25 g_(k-2)`, gives `delta >= rho_L/W_L = 2^(-(0.7737+o(1))L)`.

## 5. Meaning

* **Collatz sits between two dense regions of the pairing cube.**
  * Density-zero flips make it divergent (THM-4470 §4).
  * Flips of density `2^(-0.05L)` make it provably a tree by `L`-step descent (this theorem).
* **Consistency with Theorem A of THM-4474.** Bounded-lookahead provability is bought exactly on the Collatz-undecided classes, whose density is governed by the `0.95` exponent. Here the construction pays for them one flip each.
* **What it does not do.** It does not prove Collatz. Density-zero and small-density modifications change the truth value (DEFECT), so approximating Collatz by provable trees says nothing about Collatz itself.
