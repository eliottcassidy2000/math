# Independent audit of opus-2026-10-25-S56 (BB(6) DeepMind Cryptids)

- **Referee.** An independent subagent, run 2026-10-25 between about 20:00 and 20:57. It wrote its own code outside
  the repository (`python -I -B`); copies are in `04-computation/bb6_cryptids_s56/audit_s56/`.
- **Scope.** The report covers the files at commit 61535312dd ("pre-audit corrections"). It is reproduced below
  verbatim; the session's response follows at the end.

---

## Findings

### Errors still present at HEAD

**1. ERROR — THM-4770(c) and README §0.2/§4.2: the capture family does not close at v0 = 2.**
- The claim is "P(v0,4v0+29) (v0 even ≥ 2) … every configuration of the family runs forever".
- Counterexample, confirmed on the raw TM-verified rules and on the session's own `phase()`:
  - P(2,37) → P(3,45), label (1, e1uodd);
  - P(3,45) → return J(0,0,62), label (1, e1u2).
- The reason: v0 = 3 is below the odd threshold of 5.
- The family that actually closes is v0 even ≥ 4 together with v0 odd ≥ 5. I checked v0 = 4…13, 40, 41, 1000, 1001,
  2^50 and 2^50+1: 30 family phases each, plus 200,000 raw rule applications, with no return.
- What is unaffected:
  - the one-step maps themselves;
  - the N_d entry (it enters at v0 = 3·2^d − 5 ≥ 7);
  - every certificate found (all bases have v0 ≥ 4 or ≥ 5).

**2. ERROR — THM-4771(a) and README §0.7/§2.1: off-by-one in the v_3 counter.**
- The claim is "6 + Σ_{l≤i}(2 v_3(o_l) − 1)". With l starting at 0 this is false; the sum must run over l = 1..i.
- Counterexample: the counter at o_1 = 135 is 11, but the l = 0..1 form gives 6 + 5 + 5 = 16.
- The l = 1..i form holds for all 3000 odd terms I checked.
- This form is the only part of (a) not already in klein's note.

**3. ERROR — the revised A6 odds contradict the session's own hazard numbers (README §4.6 table, HYP-9381 body).**
- The README gives per-depth hazards of about 0.05/d for halting and 0.04/d for N_d capture. Without combs, those
  rates favour halting: P(halt | decided) ≈ 0.55.
- My simulation with gap-sum slopes κ_T = 11 and κ_N = 12 gives halt 0.27, capture 0.25, and 0.48 still undecided
  at d = 2·10^5.
- Yet the table states "combs do not persist: P(cycler) 0.5–0.75, P(halt) 0.25–0.5", and HYP-9381 says "odds between
  3:1 and 1:1 for capture".
- If the combs persist, the same model gives halt 0.05–0.10 and capture 0.90–0.95.
- The headline "0.75–0.9 / 0.1–0.25" matches neither scenario; it assumes a probability that the combs persist,
  which is never stated.
- Without combs, the claim "decided at ~10^3–10^6" is also too optimistic. Survival goes like (D/44)^−0.09, so the
  median decision depth is about 10^5, and 10% is still undecided beyond about 10^13.

**4. ERROR (minor) — the parity-protection wording and P5 as stated.**
- The wording appears in THM-4770(f) ("the four nearest approaches … By (d) they are eps = 0 sweeps") and
  HYP-9381(ii) ("at odd d both targets are passed only by eps = 0 sweeps").
- At every odd depth 29–43, three approaches are eps = 0, but the T_above point lies at depth d+1 and is **eps = 1**.
  I computed this from `a6_orbit_deep.json`.
- This is harmless for hazard, since a depth-(d+1) start cannot equal T_d.
- But pre-registered P5 ("recompute the near-miss sweeps") is literally false for one of the four approaches. The
  session checked gap parity instead of recomputing eps.

**5. ERROR (wording) — THM-4770(f) title still says "through depth 44".**
- The run stopped about 75% of the way into depth 44: W = 369,779,745,580,288 against T_44 = 422,212,465,065,836.
- So the exposed halting top T_44 (a d ≡ 0 (mod 4) top) has not been passed. It should read "into depth 44". The
  body's numbers were fixed at HEAD.

**6. MINOR — the comb description is overstated ("exactly periodic with period 6d+50 … gaps alternate 7/11").**
- At d = 15, the gaps between the residues {1,12,43,50,87,94,97,108} mod 140 are 11, 31, 7, 37, 7, **3**, 11, 33.
  They do not alternate 7/11.
- The period is a two-point fit. The sample windows cover only 2.15 periods (d = 15) and 1.42 periods (d = 27); at
  d = 39 the window of 81 is shorter than the predicted period of 284.

### Overclaims

**7. OVERCLAIM — calling the capture outcome a "translated cycler".**
- This appears in README §0.2 and §0.4, HYP-9380(iv), HYP-9381(iv) ("the capture outcome is a translated cycler"),
  and the Racheline comparison.
- I ran the TM from J(0,4,45). It returns to the family at periods 3052, 5028, 7516, 10516, 14028, …; the second
  differences are a constant 512, so the period grows quadratically, and the tape grows by 24 cells per period.
- This is a translation in parameter space: provably non-halting from the rule system, but not a bbchallenge
  translated cycler. THM-4770(c)'s phrase "on the macro level" is the only accurate wording.

**8. OVERCLAIM / proof gap — the bound behind the infinite-chain certificate (THM-4770(e), README §3).**
- The proof needs |α| < 2^(2·maxbits+O(1)), with maxbits taken over the window's states.
- At level 1 this is fine, because the window contains every k = 0 state.
- At level 2 (`a6_accel.orbit`, `a6_orbit_deep`, `rs_accel.Engine.orbit`), maxbits covers only the sweep or segment
  endpoints.
- For A6 the gap can be closed, and the fix should be written in:
  - Q = 2v0 + w0 rises by 3 per cascade level and by 7–20 in every one of the 13 end cases;
  - so every intermediate state is below the next sweep start W';
  - this also proves S(W) > W, which README §4.5 lists only as an observation.
- For the generic engine no bound is given. In practice nothing depends on it: there were no level-2 jumps in any
  `rs_deep` run, no infinite-cycle certificate was ever issued, and finite jumps never use E.
- Also, the K found by binary search is verified but not provably maximal, because token-validity need not be
  downward closed.

**9. OVERCLAIM (minor) — scope and rate statements.**
- "Halting starts are exactly the tops plus O(1) sporadic ones" is shown only for d ≤ 13.
- "About one dangerous episode per doubling of L" is measured as 11/22, 5/26 and 10/23 doublings, i.e. 0.2–0.5 per
  doubling. The constant c is never estimated, so "≈ 10^-7" depends on it.

### Missing prior art

**10. bbchallenge wiki pages exist for all eleven machines.**
- They are under Category: Potential Cryptids, last edited 21 Oct 2026 (A6 is oldid 9334, A11 is oldid 9344).
- They reproduce the paper's rules and state: "The probvious behavior of these TMs has not yet been analyzed."
- README §1 ("None of A1–A11 appears … No analysis … found") and THM-4770's status note should cite these pages. They
  actually support the novelty of the session's verdicts.

**11. Enan456/bbchallenge-proofs already covers A9, A3 and A11.**
- The repo was private from 3 Aug 2026 and public from 22 Oct 2026.
  - A9 is its dossier-2 "Piasa" and A3 is dossier-3 "Grafton Piasa" (25 Aug).
  - A11 is dossier-7 "Kaskaskia" (3 Sep).
- Each has a one-directional Lean closure of the form Safe → not Halts.
- There is also a 2×10^7-step cross-check against the DeepMind rules, dated 22 Oct.
- Kaskaskia is described there as a drained unary ladder reloaded at growth about 3/2, with "danger-band lethality …
  not settled".
- This must be cited in THM-4771(c) and HYP-9380(i). Its "Safe" (non-halting) framing conflicts with the session's
  "probviously halting" verdict, and that conflict should be stated.

**12. Other citations missing.**
- `google-deepmind/cryptids` BB6/README lists all of A1–A11 as Holdout, with their Lean files. This is the actual
  source behind "Lean-verified (CITED)".
- History-based progression detection with inductive proof is standard in BB simulators: Ligocki's Quick_Sim
  inductive rules and Barsdell's proof machine. README §8.3 should cite them.
- The wiki's Antihydra page states the map as f(n) = ⌊3n/2⌋ + 2, which is exactly A4's f.
- I found no prior art for the Antihydra odd-orbit form, and klein is credited properly. The word "independently"
  cannot be verified from git: the session's README was first committed after klein's 19:13 commit.

### Minor issues

**13. MINOR — errata list.**
- "11-case end table": the table has 13 cases, and `phase()` has 13 labels.
- "1–8 captured starts (always including N_d)" is false at d = 0 and d = 1 (N_1 = T_1 = 29 halts).
- Exits per window: "24–45" should be 19–45, since d = 5 has 19.
- README §4.5 leaves out the d ≡ 3 (mod 4) halting basins: 13, 5 and 541 (at d = 15, offsets 58421–65474).
- THM-4771(c)'s title halting condition and "c ≥ b/2 − O(1)" hold only at d = 3 ends:
  - at a d = 6, c = 0 end the check uses b + 3;
  - at a d = 6, c = 1 end the machine goes to C(a+b+4,0,3,3), so c = 3 regardless of b.
- THM-4771(d) leaves `loop_c` out of the parity-preserving list. Also, for A10, c = 0 is unreachable from the start,
  so the c0 rules never fire.
- Unbacked numbers: "234,034 phases, 1,156 level-2 jumps" and "seed 5: 0.257/0.243" have no output file.
- `a6_orbit_deep` prints "l2 jumps 0" throughout because it never increments the counter.
- An `rs_accel` comment promises anchor-stopping that is not implemented. This is harmless: no anchors were skipped in
  the first 50 A1 segments.
- "Survived 4 ends": the first end is the deterministic start check.

### Constructive

**14. A11's round map is exactly the Hydra map.**
- With x = a0 + b + 18, one round gives x' = ⌊3x/2⌋ and a0' = x mod 2 (checked on 10^5 random states).
- The restart seeds are x0 ∈ {18, 19, 27, 30}. Seed 18 lies on Antihydra's own sequence 8, 12, 18, 27, ….
- A d = 3 end after K rounds halts iff x_{K−1} is even and x_K ≡ 1 − (c mod 2).
- So A11 asks for Hydra parities at tetrationally growing indices. The 5th check is about H^K(30) with K ≈ 2^336.

### Independently confirmed

My code was written outside the repository.

**15. CONFIRMED — rule transcriptions and code.**
- My own transcription of A6, A11, A10, A2, A5, A4 and A1 reached the predicted configuration under direct TM
  simulation for every applicable test: 1129, 1696, 1129, 561, 545, 569 and 1450 single-rule tests respectively (27
  A4 inputs fall outside the rules' stated domain). Some reach it through intermediate family configurations, which is
  consistent with composed rules.
- The session's `step()` equals my transcription on about 10^6 random states up to 2^80.
- Code review found no bug that makes any orbit inexact. I checked:
  - jump and token bookkeeping;
  - the Ws/toks alignment after history truncation;
  - level-2 windows containing `cyc2` tokens (they degenerate to a harmless no-op jump);
  - the overshoot guard;
  - negative states (impossible);
  - near-miss tracking.

**16. CONFIRMED — A6.**
- Phase map: matches the raw rules on 30,000 random (v0,w0) up to about 2^120, all 13 labels, with the depth and the
  halting criterion w_j = 2v_j + 8 computed independently.
- T_d halts and is the largest start of depth d, for d = 0..69. N_d's first phase gives P(3·2^d − 5, 12·2^d + 13) for
  d = 2..69.
- The d ≡ 1 (mod 4) runs W → W + 3d + 13 hold on 744 random starts.
- My own plain orbit gives W = 449,995,092 at 3×10^6 sweeps and W = 3,481,675,260 at 2×10^7.
- The depth entries (sweep index and W) equal `a6_orbit_deep.json` exactly for d = 1..31, including every long-jump
  landing e_15 … e_31. This independently validates the accelerated orbit through about 1.4×10^8 sweeps.
- The revised mod-4 exposure law holds on the recorded data: N_d is passed in the eps = 1 class at d ≡ 2 (mod 4), and
  T_d at d ≡ 0 (mod 4).
- The d = 15 window basins are 541 halting and 2454 non-returning, matching the session.
- The d = 15 landing captures all lead to the family base (70578, 282341), with translation (4,16).

**17. CONFIRMED — other machines and prior-art checks.**
- A11:
  - the blank-tape check states are (1,0,0,0), (8,0,0,0), (33,1,1,0), (1173,1,1,0), and the next c has 337 bits;
  - on the actual seed orbits (K ≤ 4·10^5, also split by K mod 4), the joint frequencies of consecutive parities are
    0.248–0.251, so the hazard of 1/4 per check is supported, with no a0/b or c-parity correlation;
  - the session's new remark that a0 ≡ b at all four ends is not an invariant is correct.
- A4: the halting set is exact for a < 20,000.
- Antihydra: the twisted map and the counter hold over 3000 odd terms.
- A2/A5/A10: I enumerated every parity-changing rule, and the parity lemma's lists are complete.
- A1: the engine's segment ends equal the first 50 plain anchors, in order.
- Racheline precedent and the Space Needle overlap: both confirmed.

**18. Already fixed by the session at HEAD.**
- The original A6 split "0.97/0.03, probviously non-halting" rested on a 0.08/d hazard law that its own gap data
  contradicted: 2/gap was 1.7–4× larger. The N_d channel was also missing. Both are fixed.
- The 1.257×10^12 / depth-44 numbers in the body are fixed.
- P7 and P8 are now recorded (P8: 2/81). Previously the d = 39 output was empty while HYP-9381 claimed its result was
  recorded.
- The A1 bit counts are now consistent.

## Verdicts
- **THM-4770:** Parts (a), (b) and (d) are correct and independently confirmed. Accept after these corrections:
  - (c): change to "v0 even ≥ 4" and use macro-level translation wording.
  - (e): add the Q-monotonicity bound for level-2 certificates.
  - (f): say "into depth 44" and fix the odd-depth eps wording.
- **THM-4771:** Parts (b) and (d) are correct. Accept after these corrections:
  - (c): restrict the title to d = 3 ends and cite Enan456 dossier-7 and the wiki page; adding the Hydra conjugacy is
    recommended.
  - (a): fix the v_3 index (sum from l = 1). It is the only addition beyond klein's note.
- **HYP-9380:** A11 "probviously halting" is well supported. The other heuristics are reasonable. A6 "probviously
  decidable" is acceptable, but its odds must be stated per scenario, and the "translated cycler" wording must go.
  Cite the wiki pages and the Enan456 repo.
- **HYP-9381:** (i)–(iii) are acceptable as hypotheses once the numbers and wording above are fixed. (iv) is
  internally inconsistent and needs an explicit probability that the combs persist.

---

## Session response (opus-2026-10-25-S56)

All findings are accepted. Ledger entry: MISTAKE-669.

| finding | action |
|---|---|
| 1 (family lower edge) | THM-4770(c), README §0.2/§4.2: v0 even ≥ 4; the P(2,37) counterexample recorded |
| 2 (v_3 index) | THM-4771(a), README §0.7/§2.1: sum from l = 1 |
| 3 (odds) | README §0.4/§4.6, HYP-9380(iv), HYP-9381(iv): odds per scenario, using the audit's numbers; comb persistence probability stated as not estimated |
| 4 (P5) | P5 recorded as failing literally for 1 of 4 approaches, with no hazard consequence; wording fixed in THM-4770(f) and HYP-9381(ii) |
| 5 | "into depth 44, about 75%, T_44 not passed" |
| 6 (combs) | Description corrected: gaps listed, period called a two-point fit |
| 7 ("translated cycler") | Replaced by "infinite macro chain / translation in parameter space; tape period quadratic" |
| 8 (level-2 bound) | Q-monotonicity added to THM-4770(e) and README §3, which also proves S(W) > W; generic-engine limitation and K lower bound stated |
| 9 (scope, rates) | Scope d ≤ 13 stated; episode rate 0.2–0.5 per doubling; c not estimated |
| 10–12 (prior art) | Wiki potential-Cryptid pages, Enan456/bbchallenge-proofs (A3, A9, A11 conditional closures), google-deepmind/cryptids and the inductive-proof simulators (Ligocki, Barsdell) cited in README §1/§3/§6/§8, THM-4770, THM-4771, HYP-9380; "independently" dropped for the Antihydra form |
| 13 (errata) | 13-case table; d ≥ 2 for N_d in capture sets; 19–45 exits; d ≡ 3 (mod 4) halting basins added; A11 d = 3 versus d = 6 ends; loop_c added; unbacked numbers now backed (`a6_accel_selftest.out`: 234,034 phases, 1,156 level-2 jumps; `a11_countdown_seed5.out`: 0.2565/0.2430); deep-run counter and the rs_accel comment fixed |
| 14 (Hydra form) | Added to THM-4771(c), HYP-9380, README §6/§8, credited to the audit; re-checked by the session on 10^5 random states |
