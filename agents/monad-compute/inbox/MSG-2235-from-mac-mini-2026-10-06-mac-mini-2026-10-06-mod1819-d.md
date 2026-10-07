        # Message: mac-mini-2026-10-06-mod1819: D62 answered (THM-4554), uniform switches = collisions at -1 (THM-4555), Mersenne debt chain + odd least shift (THM-4556), Aut(T_k) = F_21 (THM-4557, new proof of Hanaki 2020), HYP-9213/9214, openai/math read

        **From:** mac-mini-2026-10-06-S?
        **To:** all
        **Sent:** 2026-10-06 21:14

        ---

        mac-mini-2026-10-06-mod1819, close-out. The owner sent three prompts this session:
(1) merge mod 18 / mod 19 with 7 / 63 and fractal recursion, and keep pushing toward Collatz;
(2) push the next step;
(3) "{7,21} versus odd shift distances and 37 of 60, then explore github.com/openai/math for information-compression ideas, especially the integer-multiplication paper". Two OpenAI preprints were attached: joint Dickman, and Snaky in 21.

Records:
- Notes:
  - 05-knowledge/results/mod18_mod19_seven_sixtythree_fractal_20261006.md (parts 1–2)
  - 05-knowledge/results/seven_twentyone_mersenne_openai_math_20261006.md (part 3)
- Canon: THM-4554, THM-4555, THM-4556, THM-4557. All four were independently audited; corrections are in MISTAKE-569..572.
- Hypotheses:
  - new: HYP-9213 (Mersenne coalescence), HYP-9214 (reset-2 debt resolves with probability -> 1);
  - HYP-9162 RESOLVED.
- Scripts: mod1819_20261006_{backward_sieve, clock_tower, trailing_ones_switch}.py and sevens_20261006_{mersenne_plateaus, drt_doubling_aut, debt_resolution_trend}.py, each with .out.

What changed:

1. THM-4554 (D62 answered).
   - The minimal counterexample's backward 3-adic sieve keeps a positive proportion: s_inf ∈ [0.28820, 0.29912]. The first moment alone gives >= 0.22425.
   - Moran duality: rho(θ) = 3^(θ−1)/(2^θ − 1) and min rho = 3^-(1-h). Backward first-passage words and forward glide words are mirror families across 2^K = 3^d.
   - Independent check at N = 1e7: the backward-minimal fraction is 0.2970.

2. THM-4555.
   - Run-length-uniform switches are exactly the collisions f_u(-1) = f_u'(-1) of the maps (3x+1)/2^c. The partner is (n+1)/2^D − 1 (delete D trailing ones).
   - The reset switch is the root collision (prior art: Ahmed 2016).
   - Reset-2 sources have no downward equal-length root switch, and there is no uniform 3-multiple switch.

3. THM-4556: the Mersenne line is a chain of debt states.
   - U^a(2^a − 1) = oddpart((3^a − 1)/2), i.e. the binary repunit passes through the ternary repunit.
   - The reset pairing glues {2^(2k−1) − 1, 2^(2k) − 1}. The least shift D is odd, and the nearest partner is 2^(2k) − 1 = 3(4^k − 1)/3 (the 63 = 3·21 pattern).
   - Switches are periodic mod 2^(K−2), the 2-adic clock of 3. Certified density is 15719/131072 at K = 20; opus S18 independently reached 0.1556 at K = 27 and a parity law for every source.
   - sigma(2^a − 1) (A390816) takes 104 values for a <= 6000, and 98% of odd a in [1e3, 6e3] share theirs with a smaller exponent.

4. HYP-9214: reset-2 debt resolution grows with size.
   - It is 0.19 at 16 bits, 0.82 at 2048 and 0.91 at 8192, with merges mostly late.
   - A log-density proof would make the rewrite program's obstruction a density-zero set.

5. THM-4557.
   - Doubling a doubly regular tournament on >= 7 vertices keeps Aut, so Aut(T_k) = F_21 for all k >= 3 and HYP-9162 is closed.
   - The audit found the theorem in Hanaki 2020 (arXiv:2011.06141, Thm 3.4). Ours is a new proof via the kernel of an odd skew matrix. MISTAKE-571: search association-scheme literature first.

6. {7,21}.
   - 7 = 111_2 and 21 = 10101_2 approximate −1 and −1/3 (the append maps 2x+1 and 4x+1, per S18).
   - Catalan's equation 2^3 + 1 = 3^2 makes 63 = 7·9 Zsigmondy's exception.
   - 7 is the only Mersenne prime with <2> = QR_p.

7. openai/math, read from LaTeX by six parallel readers.
   - The multiplication paper's negacyclic ring Z[y]/(y^r+1) has a 3-adic twin: the clock tower splits into 2^(3^n) − 1 and 2^(3^n) + 1 halves, and 19·27 = 2^9 + 1.
   - Its power saving is a subcritical recursion, i.e. a similarity dimension (cf. h).
   - Collisions are vanishing {2,3}-unit sums: S-unit finiteness holds per term count, and 73% of degenerate collisions come from Catalan's 9 = 8 + 1. They are also the 3-adic analogue of a dyadic rational's two expansions, and cost <= log2(n+1) bits.
   - Conditional on the OpenAI paper: (P+(n), P+(U(n))) are independent Dickman over odd n.
   - Negative cycles block every information-only exclusion of divergence.
   - The Artin, Hadamard and Diophantine clusters gave no usable inputs.

Proof status: no Collatz theorem.

Next obligations:
- HYP-9214 in log density (coalescing-walk model: the (L, Delta) sibling process).
- Push the certified Mersenne density by refining only the uncertified classes.
- The 3-adic Thorp covering A_cov(j) = O(j).
- Sharpen D62's lower end with a survival-conditioned first moment.
- Audit the readers' unaudited derivations (Dickman divisibility classes; entropy at −1) before any canon use.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
