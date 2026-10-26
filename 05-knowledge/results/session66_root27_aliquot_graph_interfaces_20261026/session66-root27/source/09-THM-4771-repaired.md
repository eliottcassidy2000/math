---
id: THM-4771
title: "Odd-orbit and countdown forms of BB(6) Cryptids. (a) Antihydra: the odd first-return map is the twisted Syracuse map y -> 3^(v-1)(3y+1)/2^v on y = -o < 0 (v = v_2(3y+1)), and the machine halts iff 6 + sum_{l<i}(2v_l - 3) = 0 for some i along the odd orbit of 27; equivalently the counter at o_i is 6 + sum_{l=1..i} (2 v_3(o_l) - 1) (klein's ANTIHYDRA.md of 2026-10-25 is the primary and fuller reference). (b) A4 (1RB0LF_1LC1RA_0RE0RD_---1LE_1LF1RC_1LC1LA): C(a,3) -> C(floor(3a/2)+2, 1 or 3 as a is even or odd), i.e. Antihydra's map floor(3n/2) at n = a+4, run while a is odd, inside a doubling cascade a_i = 2^(i+2) - 2; the machine halts iff some visit C(a,1) has a = 2^(m+2) - m - 4. (c) A11 (1RB1LA_1LC0RE_1LF1LD_1RE0LC_1LB1RA_---0LA): the orbit is a chain of countdowns; with x = a0 + b + 18 one round is exactly the Hydra map x -> floor(3x/2) with a0' = x mod 2, iterated floor(c/2) times from a seed x0 in {18, 19, 27, 30}; at a d = 3 end after K rounds the machine halts iff x_{K-1} is even and x_K = 1 - (c mod 2) (mod 2) (d = 6 ends reach such a check within one round); every non-halting d = 3 end restarts with c = floor(b/2) + O(1), so countdown lengths grow tetrationally. (d) A2, A5, A10: every non-end rule preserves the parity of the fuel (w for A2/A5, c for A10); halting needs the dangerous parity (w even, c odd), which only the end case o_small creates for A2/A5, and only b0_small and the c = 0 rules c0_a, c0_b (b = 0) create for A10 (the latter output c = 1 and reach b0_small when b = 0); each dangerous episode ends at the next cascade end of the other kind."
status: >
  PROVED from the rule systems: (b)-(d) from Dos Santos et al. A.4, A.11, A.2, A.5, A.10 (simplified forms;
  Lean-verified per the paper, CITED; re-validated against TM simulation in validate.out); (a) from the public
  Antihydra rule.
  (a) duplicates, and is weaker than, klein's
  05-knowledge/results/fixed_twists_conductors_reverse_potential_20261025/antihydra/ANTIHYDRA.md (committed
  2026-10-25 19:13): Mahler carry-map conjugacy, cylinder theorem, the ballot barrier sum v_i <= (3k-6)/2, the model
  probability phi^-6, 2e6 exact steps. klein asked opus to absorb or cite it; it is the primary reference for (a).
  The v_3 form of the counter (sum from l = 1; an off-by-one in the first draft was caught by the audit) is the only
  addition here.
  Checks:
  - (a): 2466 odd steps.
  - (b): exact against the rule system for a < 3000 (a4_hydra_reduction.out).
  - (c): the blank-tape orbit's four countdown ends; end statistics in a11_countdown_P9.out and a11_rule_check.out.
    The Hydra conjugacy x = a0 + b + 18 is the audit's; re-checked on 1e5 random states (session).
    Prior art for (c): Enan456/bbchallenge-proofs (GitHub, public 2026-10-22), dossier-7 "Kaskaskia" = A11, has a
    kernel-checked conditional closure Safe -> not Halts (Safe open, no machine decided) and describes a unary ladder
    drained and reloaded at growth ~3/2. The bbchallenge wiki page lists A11 as a potential Cryptid, probvious
    behaviour not analysed. The countdown/Hydra form and the parity check at countdown ends are the session's and the
    audit's.
  - (d): labelled steps equal the rule systems on 2e5 macro steps; episode statistics on 6e7 macro steps
    (parity_episodes_*.out).
  Independently audited (05-knowledge/results/bb6_cryptids_odd_orbits_20261025/AUDIT.md; MISTAKE-669): (b) and (d) confirmed; (a)
  index and (c) scope corrected; the A4 halting set re-checked for a < 20,000; the A2/A5/A10 parity lists enumerated
  complete.
session: opus-2026-10-25-S56
source: 05-knowledge/results/bb6_cryptids_odd_orbits_20261025/README.md
scripts:
  - 04-computation/bb6_cryptids_s56/a4_hydra_reduction_20261025.py (+ .out/.json)
  - 04-computation/bb6_cryptids_s56/a11_countdown_20261025.py (+ a11_countdown_P9.out), a11_rule_check_20261025.py (+ .out)
  - 04-computation/bb6_cryptids_s56/audit_s56/ (independent referee scripts: a11_check.py, a11_seed_parities.py, antihydra_a4_check.py)
  - 04-computation/bb6_cryptids_s56/parity_episodes_20261025.py (+ parity_episodes_A2/A5/A10.out)
  - 04-computation/bb6_cryptids_s56/rs_accel_20261025.py (labelled rule systems)
related:
  - 05-knowledge/results/fixed_twists_conductors_reverse_potential_20261025/antihydra/ANTIHYDRA.md (klein, primary for (a))
  - 01-canon/theorems/THM-4770-bb6-cryptid-a6-phase-map-halting-tops-capture-family-certified-acceleration.md
  - 05-knowledge/hypotheses/HYP-9380-classification-of-the-eleven-deepmind-bb6-cryptids.md
  - 01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-strict-real-tails-and-ordinary-stabilization.md
---

# THM-4771 — odd-orbit and countdown forms

## (a) Antihydra

H(n) = ⌊3n/2⌋ from n = 8; the counter gets +2 on even n and −1 on odd n, and the machine halts at an odd n with
counter 0.

**Odd orbit.**
- For odd o, (3o − 1)/2 = 2^(v−1) r with v = v_2(3o − 1). The v − 1 even steps give
  o' = 3^(v−1)(3o − 1)/2^v; with y = −o this is y' = 3^(v−1)(3y+1)/2^v.
- The counter changes by 2v − 3 per odd step. 8 → 12 → 18 → 27 gives counter 6 at o_0 = 27, so the machine halts iff
  6 + Σ_{l<i}(2v_l − 3) = 0 for some i.

**3-adic form.** (3o − 1)/2^v is prime to 3, so v_3(o_{l+1}) = v_l − 1. The counter then reads
6 + Σ_{l=1}^{i}(2 v_3(o_l) − 1) at o_i. The second sum starts at the first *returned* odd value o_1, not at
o_0 = 27. At i = 0 it is empty and gives 6, and the counter at o_1 = 135 is 6 + (2·3 − 1) = 11. Positivity at
return i is exactly Σ_{1≤l≤i} v_3(o_l) ≥ (i − 5)/2. Survival needs this at every i, not just on average.

The first draft summed from l = 0. The index was repaired independently by the root27-machine session (MISTAKES
2026-10-26 entry, commit 02065cd2cc) and by this session's audit (MISTAKE-669).

The odd terms must carry on average at least half a factor of 3 (a random integer carries exactly 1/2).

**klein's results.** klein's note proves much more: shift conjugacy via THM-2228, cylinders, the barrier
Σ v_i ≤ (3k − 6)/2, and P(halt | i.i.d. model) = φ^(−6).

## (b) A4

The rules are C(a,b) → C(2a+2, b−2a−2) (b > 2a+2); C(a,1) → C(2,2a+2); C(a,2) → HALT; and
C(a,2k) → C(x div 4, x mod 4 + 1) with x = 6a + 2k + 2, resp. x = 6a + 2k + 6 for b = 2k + 1.

- **Mod rules.** x is even, so b' ∈ {1, 3}.
- **b = 3** (k = 1, odd branch, a ≥ 1). x = 6a + 8.
  - a = 2m: x ≡ 0 (mod 4), a' = 3m + 2.
  - a = 2m + 1: x ≡ 2 (mod 4), a' = 3m + 3.
  - In both cases a' = ⌊3a/2⌋ + 2, and with n = a + 4 this is ⌊3n/2⌋.
- **Visits.** From C(a,1) → C(2, 2a+2) the cascade consumes Σ_{i<m}(2^(i+3) − 2) = 2^(m+3) − 2m − 8. The remainder
  is r = 2a + 2 − 2^(m+3) + 2m + 8. HALT ⟺ r = 2 ⟺ a = 2^(m+2) − m − 4.

## (c) A11

C(a,b,c,d) with d ∈ {3,6}.
- **a-loop.** C(a+2,b,c,d) → C(a,b+3,c,d), so a becomes a0 = a mod 2 and b becomes b + 3⌊a/2⌋.
- **Round.** For a0 ≤ 1 and c ≥ 2: C(a0,b,c+2,d) → C(a0+b+4, 3, c, d). After its a-loop this is
  ((a0+b) mod 2, 3 + 3⌊(a0+b+4)/2⌋).
- **End at c ∈ {0,1}, d = 3.** The rules C(0,2k+1,0,3) → HALT and C(0,2k,1,3) → HALT give the stated condition.
- **End with d = 6.** At c = 0 the identity C(a,b,0,6) = C(a,b+3,0,3) makes the check use b + 3. At c = 1 the rule
  C(a,b,1,6) → C(a+b+4,0,3,3) leads, after one round, to a d = 3 check with c = 1.
- **Hydra form (audit).** Put x = a0 + b + 18. One round sends a = a0 + b + 4 to a0' = a mod 2 and
  b' = 3 + 3⌊a/2⌋. So x' = a0' + b' + 18 = ⌊3x/2⌋ and a0' = x mod 2.
  - Restart seeds: (a0, b) = (0 or 1, 0) gives x0 = 18 or 19. C(6,3,·) collapses to (0, 12), x0 = 30. C(6,0,3,3)
    collapses to (0, 9), x0 = 27.
  - A d = 3 end after K rounds halts iff x_{K−1} is even and x_K ≡ 1 − (c mod 2) (mod 2).
  - So A11 asks for Hydra parities at tetrationally growing indices. The 5th check concerns H^K(30) with K ≈ 2^336.
  - Seed 18 lies on Antihydra's own sequence 8, 12, 18, 27, ….
- **Restarts.**
  - C(a,2k,0,3) → C(a,0,k+3,3) and C(1,2k+1,0,3) → C(6,3,k,3) start countdowns with c = k + 3 or c = k, where k ≈ b/2.
  - C(a,2k+1,1,3) → C(a,0,k+3,6), C(1,0,1,3) → C(6,0,3,3) and C(1,2k+2,1,3) → C(6,3,k,6) do the same with d = 6.
  - With d = 6: C(a,b,0,6) = C(a,b+3,0,3), and C(a,b,1,6) → C(a+b+4,0,3,3) reaches a d = 3 check after one round.
- **Growth.** b grows by a factor ≈ 3/2 per round, so after one countdown the new c is ≈ (3/2)^(c_old/2).

## (d) A2, A5, A10

Branches as in `rs_accel_20261025.py`.

- **A2.**
  - o_big, e_big, loop_a and loop_w change w by an even amount.
  - a_small and e_small output odd w: a_small gives (a−1, 1, b+2+2w) with b odd; e_small gives (b/2, 1, 2w−1).
  - o_small, on (0, b odd, w ≤ 1), outputs (0, 1, b+1+2w), with even w.
  - HALT occurs on (0, b even, 0), which needs w even.
  - The same holds for A5 with z in place of b.
- **A10.**
  - loop_b, loop_c and b0_big preserve the parity of c. b0_mid, b1, b2, b3, c0_b3 and c0_b with b ∈ {1, 2} output
    even c.
  - b0_small outputs (1, 0, 2a+2c+1) with c odd.
  - c0_a gives (a−1, b, 1), and c0_b with b = 0 gives (1, 0, 1). Both have c = 1. With b = 0 the next step is
    b0_small. With b ≥ 1 the next step (b1, b2 or b3) outputs even c, with no halting check in between.
  - HALT occurs at b = 0, c = 3.
  - In 6·10^7 macro steps every dangerous episode began with b0_small. The c0 rules never fired; the audit notes that
    c = 0 is not reached from the start.
