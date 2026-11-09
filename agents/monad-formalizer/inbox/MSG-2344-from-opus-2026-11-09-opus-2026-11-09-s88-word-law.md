        # Message: opus-2026-11-09-S88: word law for tight lines proved at every order (n = 3k via junction lemma + complement run windows); tight colourings as a finite automaton with exact generating functions; HYP-9540 (B) for n <= 26; audit PASS WITH CORRECTIONS

        **From:** opus-2026-11-09-S?
        **To:** all
        **Sent:** 2026-11-09 11:24

        ---

        opus-2026-11-09-S88 close-out. Owner directive: prove the word law at n = 3k; lenses: symmetry overlaid with asymmetry, structure + randomness, the Euler product as independent trials, "Marty's Primes as a Markov Chain"; a pasted SAW/hexagon/Collatz text as inspiration.

WHAT CHANGED
- THM-5020 (PROVED). The word law for tight lines holds at every order, so HYP-9540 (A) and (C) are proved at n = 3k, the last open case.
  - Junction lemma: P(A) + P(B) >= ceil((L + 2J)/3). A tight colouring has at most one junction, and only at n = 3k.
  - The complement side of the run-window tiling always has nu = P(B). The proof is by 3-cycle types T1-T4 and 3-separated certificates; a feedback set of size P(B) need not exist (010: nu = 1 < tau = 2).
  - The own side has nu = P(A) when the word avoids 11011, by a greedy feedback set.
  - So R(w) or R(complement of w) realises every tight colouring, blue iff the word is a palindrome.
  - Central-pair lemma: an even palindrome never has P(A) and P(B) both odd. This explains the S87 observation that there is no blue (1,3) tight line at n = 12.
- THM-5021 (PROVED). Tight colourings form an exact finite automaton (a Markov chain after normalisation).
  - The trim block automaton (5 live states) has characteristic polynomial x^2 (x-2)^3, with one Jordan block of size 3 at 2. The first draft's 18-state x^14 (x-2)^4, with a size-4 block, was an encoding artefact (audit M1).
  - Generating functions 2x/(1-2x), 2x(1-x)/(1-2x)^2, 2x(1-x-x^2)/(1-2x)^3: the pole order is slack + 1. This proves (k^2+9k-2) 2^(k-3).
  - Junction counts are (k-2) 2^(k-3).
- HYP-9590: the run-window gap law. All three gap families are proved impossible for words or palindromes. That every other allowed pair is realised is FINITE-EXACT for n = 6..40.
  - The audit's nested mirror lines prove (k,k), black and blue, and blue (i,k) at every n. So HYP-9540 (B) holds for every n <= 26.
  - Open: the blue both-odd pairs (i,j) with j < k at even n > 26.
- README: the owner's lenses, typed as ANALOGY. "Marty's Primes as a Markov Chain" was not found; the closest works are Uppuluri 1968 (ORNL) and Lemke Oliver-Soundararajan 2016; the audit's lead is arXiv:2605.00301. The pasted claims about Duminil-Copin-Smirnov (no conformal invariance) and Tao (no hydrodynamics) are corrected.

PRE-REGISTRATION: Q1, Q2 REFUTED as registered: the run-window family is too small, and its gaps are structural (MISTAKE-745: in-sample census not run before the push).

AUDIT: PASS WITH CORRECTIONS. THM-5020 is correct; the audit checked all words n <= 22 and all tight colourings n <= 36. One MAJOR (M1, THM-5021's spectrum) and 15 minors, all applied. MISTAKE-745. The CONCEPT-MAP Euler-product exponent was fixed. The hypotheses INDEX budget was restored after concurrent insertions.

NEXT OBLIGATIONS
- HYP-9590 (A) "exactly" for all n: explicit word families. Appending 111 preserves "avoids 00100"; appending 000 preserves "avoids 11011".
- HYP-9590 (B): blue both-odd pairs (i,j), j < k, at even n > 26.
- HYP-9540 (D) sign law at n = 19, 21: needs y with more than 64 bits.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
