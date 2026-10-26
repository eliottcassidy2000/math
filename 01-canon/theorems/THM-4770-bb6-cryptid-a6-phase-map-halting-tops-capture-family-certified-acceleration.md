---
id: THM-4770
title: "BB(6) Cryptid A6 = 1RB0LF_1RC1RA_0RD0RB_1LE---_0LA1LE_1RD1LA (Dos Santos et al., arXiv:2610.22603v2, A.6). (a) Phase map: on phases P(v0,w0) = J(0,v0,w0) the rule system acts by a doubling cascade v_j = A 2^j - 6, w_j = w0 - 2A(2^j - 1) + 3j (A = v0 + 6), ending at the least j with w0 + 2A + 3j + 3 < 4A 2^j, followed by an explicit 13-case end table; the quantity Q = 2v0 + w0 increases by 3 per cascade level and by 7-20 in every end case, so S(W) > W for every returning sweep. (b) Halting: a phase halts iff w_j = 2v_j + 8; a sweep start J(0,0,W) halts in its first phase iff W = T_d = 24 2^d - 3d - 16, the largest sweep start of depth d. (c) Capture family: P(v0, 4v0+29) (v0 even >= 4) -> P(v0+1, 4v0+37) and P(v0, 4v0+33) (v0 odd >= 5) -> P(v0+3, 4v0+41), both at depth 1, so every configuration of the family runs forever (a translation in macro-parameter space; on the tape its period grows quadratically, so it is not a bbchallenge translated cycler); the sweep start N_d = 18 2^d - 3d - 4 (d >= 2) enters it in its second phase. (P(2,37) -> P(3,45) -> return J(0,0,62): v0 = 2 is not in the family.) (d) Parity protection: sweep starts with eps = (w-3) mod 2 = 0 are never targets (T_d, N_d have eps = 1); at d = 1 (mod 4) such starts with u odd form runs W -> W + 3d + 13 across the window. (e) Three-point lemma: for a piecewise-affine rule system whose labelled branches are affine on (polyhedron) cap (lattice coset), three exact evaluations certify a translation chain in parameter space on a whole range, and an evaluation at k = 2^E (E >= 3 maxbits + 300, maxbits over all states the chain passes; for A6 sweep chains the sweep starts suffice by Q-monotonicity) certifies an infinite one. (f) FINITE-EXACT: every one of the first 1.33e12 sweeps of the blank-tape orbit (J(0,0,0) after 35 TM steps) returns to J(0,0,.): no halt, no capture, into depth 44 (sweep starts < 2^49; about 75% of depth 44, the top T_44 not yet passed)."
status: >
  PROVED (a)-(e) from the paper's simplified rule system for A.6. The paper states that system is Lean-verified
  (CITED; Lean files in github.com/google-deepmind/cryptids, BB6, listed there as Holdout). This session re-validated it against direct TM simulation: 60/60 predicted configurations in order
  (validate.out). The phase map (a) also equals the rule system for every W < 20000 (a6_phase_map_20261025.py).
  FINITE-EXACT (f), via the certified two-level accelerator a6_accel_20261025.py. Its accelerated orbit equals the
  plain orbit exactly at sweeps 3e6 (W = 449995092) and 2e7 (W = 3481675260, a6_accel_crosscheck_2e7.out).
  Exhaustive halting and capture sets for W < 2^18 (a6_target_sets_E18.json): each depth has 1-4 halting starts
  (always including T_d) and, for d >= 2, 1-8 captured starts (always including N_d); scope d <= 13 only. Every capture
  certificate is the family (c).
  Prior art (2026-10-25): the bbchallenge wiki page for this machine (Category: Potential Cryptids, edited 2026-10-21)
  reproduces the paper's rules and states that its probvious behaviour has not been analysed; nothing beyond the rule
  system in the paper. The jump technique is a variant of standard inductive macro proofs in BB simulators
  (Ligocki's Quick_Sim inductive rules, Barsdell's proof machine); the three-point formulation on affine label domains
  is the session's.
  Independently audited (05-knowledge/results/bb6_cryptids_odd_orbits_20261025/AUDIT.md; MISTAKE-669).
  The audit re-derived (a), (b), (d) with its own code: 30,000 random phases up to 2^120, T_d and N_d for d <= 69.
  It corrected (c) (v0 = 2 is excluded; tape period quadratic) and supplied the Q-monotonicity bound in (e).
  It matched the orbit's depth entries for d = 1..31 exactly.
session: opus-2026-10-25-S56
source: 05-knowledge/results/bb6_cryptids_odd_orbits_20261025/README.md
scripts:
  - 04-computation/bb6_cryptids_s56/deepmind_cryptids_20261025.py (+ validate.out): rule systems, TM validation
  - 04-computation/bb6_cryptids_s56/a6_phase_map_20261025.py: phase map vs rules, W < 20000
  - 04-computation/bb6_cryptids_s56/a6_accel_20261025.py: certified accelerator, self-tests
  - 04-computation/bb6_cryptids_s56/a6_orbit_deep_20261025.py (+ a6_orbit_deep_3000s.out, a6_orbit_deep.json): orbit to depth 44
  - 04-computation/bb6_cryptids_s56/a6_target_sets_20261025.py (+ a6_target_sets_E18.json): exhaustive targets
related:
  - 01-canon/theorems/THM-4771-odd-orbit-and-countdown-forms-of-bb6-cryptids-antihydra-a4-a11-a2-a5-a10.md
  - 05-knowledge/hypotheses/HYP-9380-classification-of-the-eleven-deepmind-bb6-cryptids.md
  - 05-knowledge/hypotheses/HYP-9381-a6-funnels-capture-windows-and-probviously-decidable.md
---

# THM-4770 — A6: phase map, halting tops, capture family, certified acceleration

The rule system (Dos Santos et al., A.6, simplified form) acts on J(u, v, w) = 0^∞ 1^5 10^u 1 A> 11 01^v 1^w 01 0^∞.

## (a) Phase map

Write P(v0, w0) = J(0, v0, w0) and A = v0 + 6.

**Cascade.** For w ≥ 2v + 9 the rules give J(0,v,w) → J(v+3,0,w−3), then v+3 loop steps J(u,v,w) → J(u−1,v+2,w−2),
ending at J(0, 2v+6, w−2v−9). So v_j + 6 = A·2^j and w_j = w0 − 2A(2^j − 1) + 3j. The cascade stops at the least j
with w_j < 2v_j + 9, which is the same as w0 + 2A + 3j + 3 < 4A·2^j. Write (v, w) for its last level.

**End case w ≥ 3.** Rule J(0,v,w) → J(v+3,0,w−3), then t = ⌊(w−3)/2⌋ loop steps, give J(u, 2t, eps) with
u = v + 3 − t ≥ 1 and eps = (w − 3) mod 2.

- eps = 1:
  - u = 1: HALT.
  - u = 2: return J(0,0,4t+10).
  - u even ≥ 4: J(2,u,4t+5) → P(u+4, 4t+1).
  - u odd ≥ 3: P(u, 4t+5).
- eps = 0:
  - u = 1: return 4t+12.
  - u even: J(2,4,2u+4t+5) → P(8, 2u+4t+1).
  - u odd ≥ 3: return 2u+4t+10.

**End case w ≤ 2.**

- v even: return 2v+12 (w = 0), P(8, 2v+5) (w = 1), return 2v+18 (w = 2).
- v odd: P(8, 2v+3) (w = 0), return 2v+14 (w = 1), return 2v+22 (w = 2).

Every arrow is one of the paper's rules or a composition with ≤ 2 loop steps.

## (b) Halting

HALT requires eps = 1 and u = 1, i.e. t = v + 2, i.e. w = 2v + 8.

For v0 = 0, write A = 6 and W = w0. Then w_d = 2v_d + 8 at the last level is the condition W = 24·2^d − 3d − 16 = T_d.
The depth condition for W = T_d holds exactly at d. T_d is the largest W of depth d, because the depth-d condition is
W + 15 + 3d < 24·2^d.

## (c) Capture family

**Even case.** Let v0 be even and w0 = 4A + 5.
- Depth 1, since 4A + 5 + 2A + 3 ≥ 4A and 6A + 11 < 8A.
- The last level has v = 2A − 6 and w = 2A + 8. So t = A + 2, eps = 1, and u = A − 5, which is odd and ≥ 3 for
  v0 ≥ 2.
- Output P(A − 5, 4A + 13) = P(v0 + 1, 4(v0+1) + 33).
- The family needs v0 even ≥ 4. Then v0 + 1 ≥ 5 is in the odd case. For v0 = 2, P(3, 45) has u = 2 and returns
  J(0,0,62) (audit).

**Odd case.** Let v0 be odd and w0 = 4A + 9.
- Depth 1 for A ≥ 9.
- The last level has w = 2A + 12. So t = A + 4, eps = 1, and u = A − 7, which is even and ≥ 4 for v0 ≥ 5.
- Output P(A − 3, 4A + 17) = P(v0 + 3, 4(v0+3) + 29).

The labels are (1, e1uodd) and (1, e1ueven). Neither is HALT or a return, so by induction the machine runs forever from
the family: {P(v0, 4v0+29) : v0 even ≥ 4} ∪ {P(v0, 4v0+33) : v0 odd ≥ 5}.

This is a translation in macro-parameter space, not a bbchallenge translated cycler. The audit ran the TM from
J(0,4,45) and found the return periods 3052, 5028, 7516, 10516, 14028, … growing quadratically (second difference
512), with the tape growing 24 cells per period.

**N_d.** The sweep from N_d (d ≥ 2) has its first phase at depth d. It ends with w = 6·2^d + 8, t = 3·2^d + 2 and
u = 3·2^d − 5 (odd, eps = 1). So it outputs P(3·2^d − 5, 12·2^d + 13), and 12·2^d + 13 = 4(3·2^d − 5) + 33: this is
the odd member.

## (d) Parity protection

For a sweep start W of depth d let w = W − 12(2^d − 1) + 3d. Then w(T_d) = 12·2^d − 4 and w(N_d) = 6·2^d + 8 are
even, so both targets have eps = 1.

- *Odd d.* The step 3d + 13 is even. A start with eps = 0 and u odd returns W + 3d + 13, again with eps = 0, and its
  u drops by (3d + 13)/2.
- *d ≡ 1 (mod 4).* (3d + 13)/2 is even, so u keeps its parity. The run continues until u < 3 or the depth changes,
  and no point of it is a target.

## (e) Three-point lemma and infinite certificates

Let D_L be the domain of a label sequence L: a polyhedron intersected with a lattice coset, on which the composite
G_L(x) = Mx + g is affine. Suppose b, b + D, b + KD ∈ D_L with G_L(b) = b + D and G_L(b+D) = b + 2D. Then:

- b and b + D share the coset, so D ∈ Λ.
- Convexity puts b + kD (0 ≤ k ≤ K) in the polyhedron.
- MD = D, so G_L(b + kD) = b + (k+1)D.

Along such a chain every intermediate state is integer-affine in k, so each constraint value has the form α + βk with
integers α, β and |α| < 2^(2·maxbits + O(1)), where maxbits ranges over all states at k = 0. If the chain also holds
at K = 2^E with E ≥ max(1000, 3·maxbits + 300), then β ≥ 0 for every constraint, so the chain is infinite.

**Q-monotonicity (A6, from the audit).** Let Q = 2v0 + w0.
- Q rises by exactly 3 per cascade level.
- In the 13 end cases it rises by +8 (u = 2), +11, +7, +13, +20, +13 (w ≥ 3), and by +12, +20, +16, +19, +13, +20
  (w ≤ 2).
- So S(W) > W for every returning sweep, and every intermediate state of a sweep has 2v0 + w0 < W'.
- Hence, for sweep-level (level-2) chains, maxbits over the sweep starts (plus one) bounds every state.

**Limits.** The generic engine `rs_accel_20261025.py` has no such bound proved at level 2, but it never issued a
level-2 certificate. The K returned by the search is a verified lower bound and need not be maximal; this does not
affect validity.

## (f) The orbit

`a6_orbit_deep_20261025.py` uses level-1 (phase) and level-2 (sweep) jumps, each justified by (e).

- It ran 3000 s and stopped inside depth 44 after 1.33·10^12 sweeps, at W = 369779745580288 (< 2^49). That is about
  75% of the way through depth 44; T_44 = 422212465065836 has not been passed.
- Every sweep returned to J(0,0,·). The run would have aborted on a halt, on a capture certificate, or on a sweep
  longer than 2·10^6 phases; none occurred.
- Entry offsets at d ≡ 3 (mod 4) are long jumps, e_d ≈ (4/9)·2^d. The values for d = 15, 19, ..., 39 are listed in
  the source note.
- At every odd depth 29–43, three of the four nearest approaches (both sides of N_d, and the last start below T_d)
  are eps = 0 sweeps of depth d. By (d) they cannot be targets.
- The fourth, the first start past T_d, is a depth-(d+1) start, eps = 1 in its own depth. It cannot equal T_d either.
- At even depths the orbit passes exactly one target in the eps = 1 class: T_d when d ≡ 0 (mod 4), N_d when
  d ≡ 2 (mod 4). It misses each by a positive gap; the gaps are in `a6_orbit_deep.json`.
