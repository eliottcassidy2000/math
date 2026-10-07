        # Message: mac-mini-2026-10-07-twoanchor: residual first-reset-2 rule made uniform in both run lengths (THM-4600/4601: anchors, +1 barrier D<=8000, universal state (3,1-27), ladder compiler at +1); THM-4602 chirotopes/Legendre friezes; HYP-9240/9241; audited (MISTAKE-586)

        **From:** mac-mini-2026-10-07-S?
        **To:** all
        **Sent:** 2026-10-07 16:40

        ---

        mac-mini-2026-10-07-twoanchor close-out.

Owner prompt: "work towards a variable-depth rule for the residual first-reset-2 branch, think cluster algebras, total positivity and how friezes and conway-coxeter connect with our tournament adjacent work". Four OpenAI PDFs were attached and treated as trusted.

Note: 05-knowledge/results/twoanchor_reset2_friezes_20261007.md. It complements tiling-modular-atoms' reset-two notes (anchor −1, grounding). Independent audits A and B were applied; MISTAKE-586.

COLLATZ IS STILL OPEN.

The residual branch has a second unbounded parameter: the two-run length J after the reset. On the Mersenne line, J = 1 + floor(v2(K−1)/2). This session makes the rule uniform in J.

- THM-4600 (PROVED). Anchored pair-chain debts are run-transparent.
  - If u − c_a = 3^k(v − c_a), with c_a = 1/(2^a−3), then a-runs pass unchanged and end together.
  - a = 1 is THM-4555 (−1). a = 2 is the trivial cycle +1.
  - Group form: the states form the centralizer torus in the Borel group. The anchor moves as a virtual third orbit.
- THM-4601: the two-anchor reduction of the residual branch.
  - (i) PROVED. The branch is the lag-one pair (x, x−1) with x ≡ 1 mod 8.
  - (ii) +1 barrier: PROVED for D ≤ 8000 (exact limit chains); HYP-9240 for all D. No 2-adic deletion certificate completes by Terras time j after the run end, i.e. within the two-run. 3-adic certificates are not bound.
  - (iii) PROVED. For J ≥ 3, h_3 and h_4 end every two-run in the universal state (3, 1−27), whatever K, t and J are. The proof uses the clearing identity F_222(x) − 1 = 27(F_411(y) − 1).
  - (iv) PROVED. Ladder compiler at +1: head pairs with |v| = |u|+3 and a 4z+1 ladder condition give certificates uniform in both run lengths K and J. Examples, both used on the Mersenne line:
    - (10)~(3,1,1,3), post-run depth 11: 2^1889−1 ~> 2^1886−1 and 2^5249−1 ~> 2^5246−1;
    - (1,10)~(1,1,2,2,3), post-run depth 12: 2^6129−1 ~> 2^6126−1.
  - (v) PROVED. End-of-run states depend only on (D, J) and are constant from J_0(D). Certification is a fixed chain problem in the post-run bits, and Haar coverage tends to 1.
  - NUMERICAL: the reset child collapses as J grows (0.45 → 0.003), while h_3/h_4 stay flat at about 0.38. Deletion children with D ≤ 8 cover 0.65 of the branch within 400 steps.
  - (vi) PROVED. Every Mersenne shell v2(K−1) ≥ 4 has the same certified density.
- THM-4602 (PROVED; statements KNOWN: Babai–Cameron 2000, Gunderson–Semeraro 2017/2022).
  - A tournament is a Gr(2,n) sign pattern iff all 4-Pfaffians are ±1 iff it is locally transitive.
  - Positive friezes are a slice of the transitive (totally positive) part; Conway–Coxeter friezes are its integer points.
  - Over F_p, the Legendre pattern of P^1(F_p) is Paley + sink up to switching. Singer strips close as friezes (new to the repo).
- HYP-9241 (NUMERICAL). Pair independence for Collatz translation partners: q_2 = q01 q12 q02 within 3%, exponent R(R+1)/4. Karlin–McGregor repulsion (constant π/4) is rejected, so there is no TP structure here.
- Catalan: THM-438's C_k is a signed Möbius sum (MISTAKE-060/061), so its match with frieze counts is a coincidence, not a dictionary.
- PDFs, ANALOGY only.
  - Sextic torus packets: long diagonal segments and inactive cuts correspond to runs and anchored states.
  - Anticanonical paper: "finite generation attains the limit" corresponds to finitely many heads per depth over free run parameters.

LIMITS.
- Positive pairs that reach 1 with unequal odd counts never merge.
- Orphans (2% of odd exponents in [10^3, 6·10^3]) need 3-adic, non-Mersenne or depth ∝ K certificates.
- The i = 1 head grammar is incomplete: i ≥ 2 and source ladders also occur.

NEXT.
1. HYP-9240 for all D.
2. Complete ladder grammars per (D, J) state, with the exact Haar coverage of the branch.
3. A positive-integer absorption law.
4. Prove pair independence for R = 2.
5. The −5 and −17 anchors; the tiling session's G_1 and G_5 are their inverse branches.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
