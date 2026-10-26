# BB(6) Cryptids through odd orbits: Antihydra, the eleven DeepMind accelerations, a certified accelerator (opus-2026-10-25-S56)

Owner prompt (2026-10-25), with Mazur's *Explicit positive-density Collatz convergence in logarithmic time*
(v2.1, Sept 2026) and Dos Santos et al., *At the limits of computability* (arXiv:2610.22603v2) attached:
"think in terms of odd orbits ... we must solve a Collatz-like problem in order to solve BB(6) ... read recent results
carefully to draw inspiration from them."

Typing: PROVED / CONDITIONAL / FINITE-EXACT / VERIFIED / CITED / OBSERVATION / HEURISTIC / READING.
"Rule system" always means the paper's simplified Appendix A system. The paper states these are Lean-verified (CITED).
This session re-validated each system against a direct TM simulation (VERIFIED, `validate.out`).
Everything below is a statement about those rule systems.

## 0. Summary

1. **Certified translated-cycle acceleration (PROVED, §3).** Each labelled branch of these rule systems is an affine map
   on (polyhedron) ∩ (lattice coset). Hence three exact evaluations certify a translated cycle for a whole range, and a
   fourth evaluation at k = 2^E certifies an *infinite* one, which is a proof of non-halting. Two engines:
   `a6_accel_20261025.py` (A6 phase level and sweep level) and `rs_accel_20261025.py` (A1, A2, A5, A10).
   Both reproduce plain simulation exactly where plain simulation can reach.
2. **A6 (1RB0LF_1RC1RA_0RD0RB_1LE---_0LA1LE_1RD1LA), PROVED structure (§4.1–4.3):**
   - a closed-form phase map;
   - the halting sweep starts are exactly the "tops" T_d = 24·2^d − 3d − 16 (one per depth) plus O(1) sporadic ones;
   - a **capture family**: the depth-1 translated cycler P(v0, 4v0+29) ↔ P(v0, 4v0+33). Every configuration in it runs
     forever, and every sweep start N_d = 18·2^d − 3d − 4 (d ≥ 2) enters it.
   - Every capture found anywhere (exhaustively for W < 2^18, and in about 10^4 window experiments) is this family.
3. **A6 blank-tape orbit, FINITE-EXACT (§4.4).** Computed exactly into depth 44: 1.33·10^12 sweeps, sweep starts
   up to 2^49. Every sweep returned to J(0,0,·): no halt, no capture, no uncertified long sweep.
   - At odd depths every near miss of a target comes from the "odd" (eps = 0) class, and such sweeps can never equal
     a target (PROVED parity protection, §4.5).
   - At even depths exactly one target is exposed: the halting top T_d at d ≡ 0 (mod 4), the capture point N_d at
     d ≡ 2 (mod 4).
4. **A6 funnels and verdict (OBSERVATION/HEURISTIC, §4.5–4.6).**
   - Each depth window (size 12·2^d) funnels into ~25–45 exits. From the natural entries there is a single exit at
     d ≡ 0, 2 (mod 4) for 18 ≤ d ≤ 30, so those 400 entries are effectively one sample per depth.
   - Capture combs (residues mod 6d + 50, density ≈ 0.057) surround the long-jump landings at d = 15, 27 and 39. They
     are absent at d = 19, 23, 31. This period-12 law was pre-registered as P7/P8 and both passed.
   - Heuristic verdict: **probviously decidable, not probviously non-halting.** Three channels (halting at T_d,
     capture at N_d, capture at the combs) carry a divergent total hazard. A translated cycler is favoured,
     P ≈ 0.75–0.9 against P(halt) ≈ 0.1–0.25. The decision comes around depth 10^2–10^3, far out of reach.
5. **A11 is probviously HALTING, at tetrational times (§6, HEURISTIC; structure PROVED).**
   - Its orbit is a sequence of countdowns, each a Hydra-type map b ↦ 3 + 3⌊(b+a0+4)/2⌋ iterated c/2 times.
   - At each countdown end a parity check halts with probability 1/4. End parities were measured uniform over random
     countdowns (0.262, 0.233), and directly on the rule system (0.253–0.272 for every (c mod 2, d)).
   - The blank-tape orbit has survived 4 checks. The 5th comes after about 2^336 rounds.
6. **A2, A5, A10 are doubly parity-protected (PROVED lemma, measured rates, §5).**
   - Halting needs a rare parity-flipping end case and then an exact hit inside that one episode.
   - Measured: about one dangerous episode per doubling of scale. Heuristic residual ≈ c/L ≈ 10^-7 at the current
     L ≈ 2^26.
   - A1, A3, A4, A7, A8, A9 have exponentially growing parameters with exact-coincidence halting conditions, so their
     residual halting probability is ≲ 2^-(current bits), between 2^-2000 and 2^-276000.
7. **Odd orbits (PROVED, §2).**
   - Antihydra's odd orbit is the *twisted Syracuse map* y ↦ 3^(v−1)(3y+1)/2^v on y = −o < 0. Its counter is the
     valuation walk 6 + Σ(2v_l − 3), equivalently 6 + Σ(2 v_3(o_l) − 1).
   - A4's visit map contains exactly Antihydra's parity dynamics (f(a) = ⌊3a/2⌋ + 2, conjugate to ⌊3n/2⌋ by n = a+4),
     inside a doubling cascade. It halts iff a visit equals 2^(m+2) − m − 4.

## 1. What the two papers (and the community) supplied

- **Mazur v2.1 (READING).**
  - §2 works entirely with odd orbits: the Syracuse map S(x) = (3x+1)/2^a(x), inverse histories as valuation words w
    with weight 3^d·2^(−A(w)), and injective affine maps F_w : Z/3^t → Z/3^(t+d) whose weighted transfer preserves
    means up to 2^(−A(w)) (his Lemma 2.1).
  - §3 selects words by the *first crossing of a barrier* by the valuation partial sums. It controls them by Tao's
    fine-scale mixing (his eq. (2.3)) uniformly in the seed, which gives z(M) > 0 for one seed M = (4^j − 1)/3.
    His Lemma 4.1 is that these seeds permute residues mod 3^q, by lifting the exponent.
  - The Antihydra halting event is exactly a barrier first-crossing for the odd-orbit valuation sum (§2.1). Mazur's
    machinery counts the *density* of such crossings over starting values. It says nothing about one fixed orbit. That
    is the same wall as for Collatz, and for the forward-growing Antihydra map the density statement is elementary
    (the first k parities are a function of n mod 2^(k+1)).
- **Dos Santos et al. (CITED).** Eleven new Collatz-like accelerations (Appendix A). The status of all eleven is
  left open. The paper's lesson, which we follow: find exact macro rules first, verify them, then reason on the macro
  level. The Appendix A.11 transitions were independently checked by another agent today (MISTAKES 2026-10-25 entry,
  "576 stated native controls"). A p.5 Antihydra clock display is excluded from dependencies there and here.
- **bbchallenge wiki (READING, 2026-10-25).**
  - Classification vocabulary: probviously halting / non-halting / decidable.
  - Racheline's tetrational Cryptid 1RB1LE_0LC0LB_1RD1LC_1RD1RA_1RF0LA_---1RE, rated "3/5 translated cycler, 2/5 halting",
    is the precedent for §4.6 and §6.
  - None of A1–A11 appears on the wiki's Cryptid list. No analysis of them was found (web search for the machine
    strings; wiki Cryptid and Space Needle pages). A11 shares its first three states with Space Needle but differs in
    D, E, F.

## 2. Odd-orbit forms

### 2.1 Antihydra (PROVED, elementary)

Antihydra iterates H(n) = ⌊3n/2⌋ from n = 8 with a counter: +2 on even n, −1 on odd n. It halts when an odd n is met
at counter 0. Let o be odd. Then H(o) = (3o − 1)/2 = 2^(v−1) r with r odd and v = v_2(3o − 1). Each of the following
v − 1 even steps multiplies by 3/2, so the next odd value is

  o' = 3^(v−1) (3o − 1) / 2^v,   and with y = −o:   y' = 3^(v−1) (3y + 1)/2^v,  v = v_2(3y + 1),

a Syracuse step twisted by the unit 3^(v−1). The counter changes by 2(v − 1) − 1 = 2v − 3 per odd step. From 8 → 12 →
18 → 27 the first odd value is o_0 = 27 with counter 6. So **Antihydra halts iff 6 + Σ_{l<i} (2v_l − 3) = 0 for some
i**. This was checked against the direct iteration for 2466 odd steps.

Since (3o − 1)/2^v ≢ 0 (mod 3), v_3(o_(l+1)) = v_l − 1, and the counter equals 6 + Σ_{l≤i}(2 v_3(o_l) − 1). In words,
the odd terms must carry on average at least half a factor of 3. A random integer carries exactly 1/2 (Σ 3^(−k)).
Under the 2-adic Haar heuristic Antihydra's odd terms carry 1, since P(v_3 ≥ k) = 2^(−k). The drift +1 per odd
step is this difference.

### 2.2 A4 (PROVED from the rules; checked by `a4_hydra_reduction_20261025.py`)

Rules (A.4): C(a,b) → C(2a+2, b−2a−2) if b > 2a+2; C(a,1) → C(2,2a+2); C(a,2) → Halt; and two mod-4 rules.
- (i) The mod rules always output b' ∈ {1, 3}, since 6a+2k+2 and 6a+2k+6 are even.
- (ii) For a ≥ 1, C(a,3) → C(f(a), 1 if a even else 3) with f(a) = ⌊3a/2⌋ + 2. With n = a + 4 this is ⌊3n/2⌋, and
  a is odd exactly when n is.
- (iii) From a visit C(a,1) the cascade a_i = 2^(i+2) − 2 consumes 2^(i+3) − 2 per level from 2a+2. The remainder r
  decides: r = 1 is the next visit, r = 2 halts, and otherwise a mod rule leads to a visit (b' = 1) or to a Hydra run
  (b' = 3: iterate f while a is odd).
- (iv) Halting: r = 2 ⟺ a = 2^(m+2) − m − 4. This is exact against the rule system for all a < 3000.

Over 3000 visits: 0.76 bits gained per visit, half of the visits enter a Hydra run, and run lengths are geometric. The
heuristic residual after the first visit is Σ 2^(−bits) ≈ 0.24, and it is dominated by the first visits.

## 3. The certified accelerator (PROVED)

**Labels.**
- *A6.* Labels are the phase end cases of §4.1 together with the depth j.
- *Generic engine.* Labels are the rule branches. A collapsed self-loop is split by which bound stops it, and by the
  parity entering a floor, e.g. A10's loop gives labels `loop_b` (2b+1 ≤ c) and `('loop_c', c mod 2)`.
- Every label is an affine map on D_l = (polyhedron) ∩ (lattice coset): all conditions are linear inequalities with
  integer coefficients, plus parities.

**Lemma (three points).** Fix a label sequence L = (l_1..l_p). Its domain D_L is again a polyhedron intersected with
a coset, and its composite G_L is affine. Suppose b, b + D, b + KD ∈ D_L with G_L(b) = b + D and G_L(b + D) = b + 2D.
Then b + kD ∈ D_L and G_L(b + kD) = b + (k+1)D for every 0 ≤ k ≤ K.

*Proof.*
1. b and b + D lie in one coset, so D lies in its lattice and b + kD is in the coset for every k.
2. The polyhedron is convex, so it contains b + kD for 0 ≤ k ≤ K.
3. G_L(x) = Mx + g, and MD = G_L(b+D) − G_L(b) = D. Hence G_L(b + kD) = b + D + kD. ∎

**Infinite chains.** Along the chain, every intermediate state is integer-affine in k: x_i(k) = x_i(0) + k·δ_i with
δ_i ∈ Z^n. Every constraint value is therefore an integer-affine function of k, with |value at k = 0| < 2^(2·maxbits+O(1)).
If the chain also holds at K = 2^E with E ≥ 3·maxbits + 300 (and E ≥ 1000), no constraint can decrease along the ray.
So the chain holds for all k ≥ 0. When no label of L is HALT or a return, this is a **certificate of non-halting**.

**Jumps.**
- The maximal K is found by exponential plus binary search, using only verified evaluations.
- Level 1 works inside a sweep or segment, with fresh history.
- Level 2 works on whole sweeps. The token sequence of a sweep is a function of its start, it fixes the underlying
  label sequence, and the lemma applies to that sequence.

**Checks (VERIFIED).**
- *A6 engine.*
  - It equals the phase map for every W < 20000 (the phase map in turn equals the rule system there).
  - Its accelerated orbit equals the plain orbit at sweep 3·10^6 (W = 449995092). It computed 234,034 phases and
    took 1,156 level-2 jumps.
- *Generic engine.*
  - Labelled steps equal the rule systems for 2·10^5 macro steps.
  - A2, A5 and A10 match plain stepping.
  - A1 reproduces all 51 anchor states that a 1.5·10^8-step plain run reaches.

## 4. A6

### 4.1 Phase map (PROVED from the rules)

Write P(v0, w0) = J(0, v0, w0) and A = v0 + 6. While w ≥ 2v + 9, the rules give

  J(0,v,w) → J(v+3, 0, w−3) → (loop, v+3 times) → J(0, 2v+6, w−2v−9).

So v_j = A·2^j − 6 and w_j = w0 − 2A(2^j − 1) + 3j. The phase ends at the least j with w0 + 2A + 3j + 3 < 4A·2^j.

**End at J(0, v, w), w ≥ 3.** Put t = ⌊(w−3)/2⌋, u = v + 3 − t ≥ 1 and eps = (w − 3) mod 2.

| eps | u = 1 | u = 2 | u even | u odd |
|---|---|---|---|---|
| 1 | HALT | return 4t+10 | P(u+4, 4t+1) | P(u, 4t+5) |
| 0 | return 4t+12 | → P(8, 2u+4t+1) | P(8, 2u+4t+1) | return 2u+4t+10 |

The u = 2 entry in the eps = 0 row is the even case (P(8, 2u+4t+1)).

**End at w ≤ 2.**

| | w = 0 | w = 1 | w = 2 |
|---|---|---|---|
| v even | return 2v+12 | P(8, 2v+5) | return 2v+18 |
| v odd | P(8, 2v+3) | return 2v+14 | return 2v+22 |

"Return W'" means J(0,0,W'). Every entry is one of the paper's rules composed with at most two loop steps; for
example J(2, x, y) → J(0, x+4, y−4). The table equals the rule system for all W < 20000 (`a6_phase_map_20261025.py`).

**Halting.** A phase halts iff w_j = 2v_j + 8. For a sweep start (v0 = 0) the first phase halts iff W = T_j = 24·2^j − 3j − 16,
which is the largest W of depth j.

### 4.2 The capture family (PROVED)

For v0 even ≥ 2 and v0 odd ≥ 5:

  P(v0, 4v0 + 29) → P(v0 + 1, 4v0 + 37) = P(v0+1, 4(v0+1) + 33),
  P(v0, 4v0 + 33) → P(v0 + 3, 4v0 + 41) = P(v0+3, 4(v0+3) + 29).

*Proof.*
- **Case v0 even.** Write w0 = 4A + 5 with A even.
  - Depth 1: w0 + 2A + 3 ≥ 4A, and w0 + 2A + 6 < 8A ⟺ A > 5.5.
  - Then v = 2A − 6 and w = 2A + 8, so t = A + 2, eps = 1, and u = A − 5 is odd and ≥ 3 when v0 ≥ 2.
  - Output P(u, 4t + 5) = P(A − 5, 4A + 13).
- **Case v0 odd.** Write w0 = 4A + 9 with A odd.
  - Depth 1 holds for A ≥ 9.
  - Then w = 2A + 12, so t = A + 4, eps = 1, and u = A − 7 is even and ≥ 4 when v0 ≥ 5.
  - Output P(u + 4, 4t + 1) = P(A − 3, 4A + 17).
- The family recurs with v0 increasing by 4 every two phases. Neither label is HALT or a return, so the machine runs
  forever from any configuration in the family. ∎

**Main capture chain.** For d ≥ 2 the first phase of the sweep from J(0,0,N_d), N_d = 18·2^d − 3d − 4, has depth d. It
ends with w = 6·2^d + 8 and u = 3·2^d − 5 (odd, eps = 1), so it passes to P(3·2^d − 5, 12·2^d + 13). That is the odd
member of the family.

### 4.3 Exhaustive target sets, W < 2^18 (FINITE-EXACT, `a6_target_sets_20261025.py`)

| depth d | halting W − 24·2^d | captured W − 18·2^d |
|---|---|---|
| 2 | −29, −22 | −32, −21, −10 |
| 4 | −90, −28 | −16 |
| 6 | −34 | −22, −14, 23, 46, 66, 79, 122, 135 |
| 7 | −1363, −1352, −63, −37 | −25 |
| 8 | −1254, −40 | −28, 1036 |
| 10 | −46 | −34, 2886, 2899, 2978, 2991 |
| 11 | −20587, −20576, −49 | −37, 1781 |
| 12 | −19698, −52 | −40 |
| other d ≤ 13 | −3d − 16 only | −3d − 4 only |

- Every halting set contains the top (−3d − 16), and every capture set contains N_d (−3d − 4).
- Every capture certificate found (here and in all later experiments) is the §4.2 family:
  period 2, translation (4, 16), labels (1, e1uodd), (1, e1ueven).

### 4.4 The blank-tape orbit (FINITE-EXACT; `a6_orbit_deep_20261025.py`, `a6_orbit_deep*.out/json`)

- The orbit runs 3000 s and stops inside depth 44 after 1.33·10^12 sweeps, at W = 369779745580288 (< 2^49); resume
  from there.
- Every sweep returned to J(0,0,·). A sweep exceeding the 2·10^6-phase cap would have aborted the run; none did.
- At depth d one sweep applies the paper's loop rule J(u,v,w) → J(u−1,v+2,w−2) about W/2 times, about 10^14 rule
  applications at depth 44.
- Cross-check: the accelerated orbit equals the plain phase-map orbit exactly at 3·10^6 and at 2·10^7 sweeps
  (W = 3481675260; `a6_accel_crosscheck_2e7.out`). The latter run passes a sweep of 2,945,221 phases.
- Entry offsets e_d (first W of depth d minus T_{d−1}):
  - d ≡ 3 (mod 4): a big jump, e_d ≈ (4/9)·2^d. Values: e_15 = 14694, e_19 = 233174, e_23 = 3728426,
    e_27 = 59652526, e_31 = 954437418, e_35 = 15270995110, e_39 = 244335917530.
  - Otherwise e_d is 8–252.
- Near misses: the gaps (below, above) to N_d and T_d are in `a6_orbit_deep.json` (the depth table is reprinted in
  `a6_orbit_deep_3000s.out`).

### 4.5 Parity protection, funnels, basins

**Parity lemma (PROVED).**
- A sweep start W of depth d is in class eps = (w − 3) mod 2 with w = W − 12(2^d − 1) + 3d. Both targets have
  eps = 1, since w(T_d) = 12·2^d − 4 and w(N_d) = 6·2^d + 8. A point at odd distance from a target therefore has
  eps = 0.
- For odd d the step 3d + 13 is even, so eps = 0 sweeps with u odd return W + 3d + 13 in the same class.
- For d ≡ 1 (mod 4), u also keeps its parity: the run W → W + 3d + 13 crosses the whole window, and no point of it is
  a target.
- Observed: at every odd depth 29–43 all four nearest-approach gaps are odd. This was pre-registered as P5.

**Coalescence (FINITE-EXACT, `a6_basins_20261025.py`).** Fates of all W of the window, for 4 ≤ d ≤ 17:
- S(W) > W always holds inside a window.
- The whole window, 12·2^d − 3 starting points, funnels into 24–45 distinct exits.
- Halting basins: 1 element (the top) at d ≡ 1, 2 (mod 4). At d ≡ 0 (mod 4): 109 (d = 8), 4700 (d = 12), 45 (d = 16).
- Capture basins: 1–8 elements, except 396 at d = 8 and 2454 at d = 15.

**Transfer maps from the natural entries (FINITE-EXACT, `a6_transfer_6_30_400.out`).** Entries T_{d−1} + o,
o = 1..400:
- Exactly one exit at d = 12, 18, 20, 22, 23, 24, 26, 28, 30, and it is always the orbit's actual exit.
- 10–25 exits at d ≡ 1 (mod 4).
- Halts only at d = 7 (11/400).
- Captures at d = 6 (1), 8 (72), 15 (46), 27 (15).

**Uniform basin densities (sampling, `a6_basin_sample_12_26.out`).**
- Halting: 0.102 ± 0.007 at d = 12 (exhaustive value 0.096), 0.0015 at d = 15, and 0/2000 at every other d in 13–26.
- At d = 24, 1988 of 2000 uniform samples share the main-stem exit 108: almost the whole window is one basin.
- Captures: 0.006 at d = 15, otherwise 0/2000 (d = 12–26).

**Landing neighbourhoods (`a6_landing_hazard_R150.out`).** Entries T_{d−1} + e_d + δ, |δ| ≤ 150:

| d | d mod 12 | captures /301 | halts | actual landing |
|---|---|---|---|---|
| 15 | 3 | 17: residues {1,12,43,50,87,94,97,108} mod 140 = 6d+50 (8/140) | 0 | exit 65 (δ = 1 is captured) |
| 19 | 7 | 0 | 0 | exit 129 |
| 23 | 11 | 0 | 0 | exit 137 |
| 27 | 3 | 17: 12 residues mod 212 = 6d+50 (12/212) | 0 | exit 189 |
| 31 | 7 | 0 (pre-registered P7: passed) | 0 | exit 121 |
| 39 | 3 | 2/81 for |δ| ≤ 40, at δ = −39, −28 (pre-registered P8 ≥ 2%: passed) | 0 | exit 49 |

The d = 39 row comes from `a6_landing_hazard_d39_R40.json`. Both combs are exactly periodic with period 6d + 50 and
density ≈ 0.057, and their gaps alternate 7/11.

### 4.6 Hazard and verdict (HEURISTIC)

**Model.** The window dynamics coalesces, so the fate of a traversal is a function of its entry. The 400 natural entries
per depth are therefore effectively *one* sample: their zero halts at d = 8..30 is weak evidence. Hazards must be read
off the orbit's actual path, channel by channel, using the parity classes of its near misses (depth table in
`a6_orbit_deep.json`).

| depth class | target the orbit passes | class of the two bracketing sweeps | sum of the two gaps |
|---|---|---|---|
| d ≡ 1 (mod 4) | T_d and N_d | eps = 0 (cannot hit) | 3d + 13 (one e0uodd run step) |
| d ≡ 3 (mod 4) | T_d and N_d | eps = 0 (cannot hit) | 6d + 30 |
| d ≡ 0 (mod 4) | T_d (halting); N_d is jumped over | eps = 1 (exposed) | 120, 172, 128, 228, 272, 300, 200, 236, 596 (d = 8..40) |
| d ≡ 2 (mod 4) | N_d (capture); T_d is jumped over | eps = 1 (exposed) | 88, 180, 252, 276, 324, 268, 420, 468, 516 (d = 10..42) |

**Channels.** Assume the offset of the passing orbit is uniform among the eps = 1 positions.

| channel | probability per event | history (expected hits) |
|---|---|---|
| Halting at T_d (d ≡ 0 mod 4) | ≈ 2/(gap sum) | 0.09 over d = 8..40 |
| Capture at N_d (d ≡ 2 mod 4) | ≈ 2/(gap sum) | 0.08 over d = 10..42 |
| Capture at a landing comb (d ≡ 3 mod 12) | ≈ 0.057 | 0.17 over d = 15, 27, 39 |

None of these events was observed, which is consistent with the totals.

**Forward.**
- The gap sums grow roughly like 7–15 d, so the first two channels contribute O(1/d) per depth: about 0.05/d for
  halting and 0.04/d for N_d capture.
- If the combs persist, the comb channel contributes about 0.0047 per depth, a constant.
- The total hazard diverges in either case.

**Verdict: probviously decidable** (not probviously non-halting). Outcome odds:

| scenario | P(translated cycler) | P(halt) | when decided |
|---|---|---|---|
| Combs persist | ≈ 0.9 | ≈ 0.1 | typically by depth a few hundred |
| Combs do not persist | ≈ 0.5–0.75 | ≈ 0.25–0.5 | depth ~10^3–10^6 |

All of this is far beyond exact simulation, whose cost grows like 2^(0.47d). The translated cycler, if it comes, is
already certified by THM-4770(c).

## 5. A2, A5, A10 (double parity protection) and A1

**Lemma (PROVED from the rules).**
- *A2, A5.* Every rule except the end cases changes the fuel w by an even amount. The end cases set w odd (a_small,
  e_small, a1_small) or even (only o_small: a = 0, b odd, w ≤ 1, giving (0, 1, b+1+2w)). HALT needs w = 0 at an a = 0,
  even-b state, so it can only occur inside an even-w episode. Each episode lasts until the next a_small/e_small,
  typically one cascade.
- *A10.* The same holds for c with "odd" and "even" swapped. HALT needs b = 0 and c = 3. Odd c is created only by:
  - b0_small (b = 0, c ≤ 2);
  - the c = 0 rules c0_a and c0_b with b = 0. These output c = 1. With b = 0 the next step is b0_small; with b ≥ 1 the
    next step resets c to even, with no halting check in between.

**Measured (`parity_episodes_*.out`, 6·10^7 macro steps each).**

| machine | dangerous episodes | log2 L at their start | Σ 1/L |
|---|---|---|---|
| A2 | 11 (all from o_small) | 3, 6, 7, 8, 10, 14, 16, 19, 21, 22, 24 | 0.129 |
| A5 | 5 | 5, 9, 10, 12, 18 | 0.019 |
| A10 | 10 (all from b0_small) | 3, 8, 11, 12, 13, 13, 15, 15, 18, 23 | 0.128 |

- That is about one episode per doubling of scale.
- Inside episodes the smallest halting margins are 2 (at L = 9), then values of order L/8.
- With P(halt | episode) ≈ c/L, the residual is ≈ 2c/L ≈ 10^-7 at the current L ≈ 2^25.5–2^26.3: probviously
  non-halting.
- The generic engine gives these machines no speed-up, because their sweeps are chaotic: 2a + w is nearly conserved
  and the position on that circle follows an expanding map.

**A1.**
- The only translated cycle is (a1b, c): A(1, b, 0, d) → A(1, b+2, 0, d−1), run d times. The rest of each segment is a
  doubling cascade (db5, loop_a0).
- The certified engine reaches 6,711 bits after 17,241 segments in 2400 s (`rs_deep_A1.json`), against 2^28 for plain simulation
  in 60 s. Neither HALT nor CAPTURE occurs.
- Halting needs d = b + 2 exactly at a cascade end, a range of about 1.5b. Residual ≲ 2^(−4000): probviously
  non-halting.

## 6. A11 is probviously halting (structure PROVED, verdict HEURISTIC)

**Structure (from the rules).**
- After its a-loop, C(a,b,c,d) becomes C(a mod 2, b + 3⌊a/2⌋, c, d).
- With a0 = a mod 2 ≤ 1 and c ≥ 2, one round is C(a0, b, c+2, d) → C(a0+b+4, 3, c, d). So the countdown iterates
  (a0, b) ↦ ((a0+b) mod 2, 3 + 3⌊(a0+b+4)/2⌋) exactly c/2 times. That is a Hydra-type map, b ×1.5 per round.
- At the end (c ∈ {0,1}, d = 3), HALT occurs iff a0 = 0 and (c = 0 and b odd, or c = 1 and b even).
- Every other end case restarts a countdown with c ≈ b/2: C(a,2k,0,3) → C(a,0,k+3,3), C(1,2k+1,0,3) → C(6,3,k,3), and
  so on. The d = 6 ends reach a d = 3 check within one round.
- So the countdown lengths grow tetrationally.

**Measured.**
- End parities are uniform. Over random starts (`a11_countdown_20261025.py`) the halting fraction is 0.257/0.243
  (seed 5) and 0.262/0.233 (seed 77, c/2 ∈ [400, 1200); pre-registered P9).
- Directly on the rule system (`a11_rule_check_20261025.py`), from 6000 random configurations C(a0, b, c, d) with
  c/2 ∈ [100, 400): the first countdown end halts with frequency 0.253, 0.272, 0.268, 0.256 for (c mod 2, d) =
  (0,3), (0,6), (1,3), (1,6).
- The blank-tape orbit's ends (macro step, a0, b mod 2, c) are (1,0,0,0), (8,0,0,0), (33,1,1,0) and (1173,1,1,0).
- All four ends have a0 ≡ b (mod 2). This is not an invariant: the round map gives a0' ≡ b' exactly when
  a0 + b ≡ 1, 2 (mod 4).
  None halts; each had probability 3/4 of not halting.
- The current countdown has c ≈ 2^337. Its end, the 5th check, comes after about 2^336 rounds.

**Verdict.** A11 halts with heuristic probability 1, at the first check that fails. That happens at a time around a
tower of height ~5. This is far below the BB(6) champion but forever out of reach of direct simulation. The precedent
is Racheline's tetrational Cryptid.

## 7. Classification of the eleven (HEURISTIC; see HYP-9380)

| machine | growth (60 s plain; certified engine) | halting condition | residual / verdict |
|---|---|---|---|
| A1 | exponential per segment (6,711 bits) | d = b + 2 at cascade end | ≲ 2^-6000; non-halting |
| A2 | slow, L ≈ 2^25.5 | double parity | ≈ 10^-7; non-halting |
| A3 | exponential (2^199945) | K = M, small (K, M) | ≲ 2^-199945; non-halting |
| A4 | exponential, 0.76 bits/visit | a = 2^(m+2) − m − 4 | ≲ 2^-2000; non-halting |
| A5 | slow, L ≈ 2^26.3 | double parity | ≈ 10^-7; non-halting |
| A6 | slow, W ≈ 2^49 at depth 44 | exact tops T_d; capture family | probviously decidable: translated cycler ≈ 0.75–0.9, halt ≈ 0.1–0.25 |
| A7 | exponential (2^186527) | coincidence, then a parity | ≲ 2^-186000; non-halting |
| A8 | exponential (2^276185) | a − k ∈ {0,1}, k − a ∈ {4,5} | ≲ 2^-276000; non-halting |
| A9 | exponential (2^272084) | b − c ∈ {0,1} | ≲ 2^-272000; non-halting |
| A10 | slow, L ≈ 2^25.8 | double parity | ≈ 10^-7; non-halting |
| A11 | tetrational countdowns | parity check, P = 1/4 per countdown end | **halting**, tetrational time |

Every "non-halting" in this table means probviously non-halting.

## 8. What this suggests for the repo's open problems

1. **Odd orbits plus barriers.**
   - Antihydra halting is a first crossing of a barrier by an odd-orbit valuation sum. Mazur's §3 is a machine for
     exactly such first crossings, and it averages over seeds.
   - The fixed-orbit question is untouched by it, as for Collatz. The repo's Collatz stance (refuel bill as
     codelength, minimum-map renewal, two-sheet calculus) is the right frame: what is provable is density or ensemble
     statements.
   - The odd-orbit twisted-Syracuse form is the bridge between Antihydra and that machinery.
2. **Coalescence is the common mechanism.**
   - A6 windows funnel 12·2^d starting points into O(d) exits. This is the Cryptid face of the repo's Collatz
     coalescence results (rank-one coalescence S22; coalescence synthesis 2026-10-04).
   - Hazard must be computed as basin mass at the actual entries. The naive "near-miss distance" hazard overstates it,
     because a near miss in the wrong parity class carries zero risk.
3. **Certified acceleration is cheap and general.**
   - The three-point lemma turns any piecewise-affine rule system with linear branch conditions into one with
     certified jumps and certified translated cyclers.
   - It applies verbatim to the repo's base-p Collatz maps, 5x+1-type maps, and the Erdős ternary-type Cryptids, wherever
     a translated cycle is the mechanism.
4. **Classify by mechanism, not by "Collatz-like".**
   - Exact coincidences of exponentially growing parameters are safe.
   - Parity-protected slow growth is polynomially safe.
   - Funnels with exposed targets of both kinds are probviously decidable.
   - A constant-probability check at tetrational times is probviously halting.
   - Only the last two can ever change the BB(6) picture.

## 9. Pre-registration outcomes (PREREGISTRATION.md)

- P1: no halts in the A1/A2/A5/A10 runs. Holds; A2/A5/A10 were stopped at 540 s because the engine gave them no
  speed-up.
- P2: A6 depth ≥ 43 within 3000 s, no halt or capture. Holds: depth 44.
- P3: d = 28 and d = 30 have 1 exit, d = 29 has 25 exits. Holds.
- P4: no halts at d = 28–30. Holds.
- P5: parity protection at odd depths 29–43. Holds.
- P6: every capture is the §4.2 family. Holds.
- P7: holds. d = 31 has 0/301 captures and 0 halts.
- P8: holds, narrowly. d = 39 has 2/81 captures (2.5%, threshold 2%) and 0 halts.
- P9: holds (0.262, 0.233).

## 10. Files

All scripts are in `04-computation/bb6_cryptids_s56/`, with an `.out`/`.json` next to each.

- **Rule systems and validation:** `deepmind_cryptids_20261025.py`.
- **A6:** `a6_phase_map`, `a6_fast`, `a6_accel` (engine), `a6_orbit_deep`, `a6_target_sets`, `a6_basins`, `a6_transfer`,
  `a6_basin_sample`, `a6_landing_hazard`, `a6_window_hazard`, `a6_sweep_map`.
- **Generic engine:** `rs_accel`, `rs_deep`.
- **Other machines:** `parity_episodes`, `a11_countdown`, `a4_hydra_reduction`, `uncapped_residues` (an earlier capped
  version, `halting_opportunities`, had capped channels whose residue claims are void), `slow_machines_sweeps`,
  `long_run`.
