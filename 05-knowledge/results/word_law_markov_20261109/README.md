# The word law for tight lines at n = 3k, and tight colourings as a Markov chain

opus-2026-11-09-S88. Canon: [THM-5020](../../../01-canon/theorems/THM-5020-word-law-for-tight-lines-at-every-order.md),
[THM-5021](../../../01-canon/theorems/THM-5021-tight-colourings-as-a-markov-chain.md),
[HYP-9590](../../hypotheses/HYP-9590-run-window-gap-law-and-residual-pairs.md); updates
[HYP-9540](../../hypotheses/HYP-9540-packing-number-line-laws-tight-words-allowed-pairs.md) ((A) and (C) PROVED at
every n; (B) for every n <= 26). Independent audit: [audit/AUDIT.md](audit/AUDIT.md), PASS WITH CORRECTIONS, all
applied (MISTAKE-745).

**Owner directive (2026-11-09).** "Start another session proving the word law at n = 3k, consider blue and black
lines being symmetry overlayed with asymmetry and primes being structure combined with randomness and seeing how the
deep structure behind the euler product formula as a chain of independent trials and marty's 'Primes as a Markov
Chain' intuition", with a pasted text as inspiration, archived in [INSPIRATION_PASTE.md](INSPIRATION_PASTE.md).

## 1. The word law at n = 3k (THM-5020, PROVED)

THM-4951 proved, at n != 0 (mod 3), that the triple words of tight lines (nu(t) + nu(t + 1) = floor(n/3)) are
exactly the tight colourings, and that the blue ones are exactly the palindromes. The case n = 3k was open because the
run-window tiling fails there on words like 1110111. The proof has three ingredients.

**Junction lemma.** A junction is an isolated letter between two runs of the other letter of length >= 2: the centre
of 11011 or 00100.
- Deleting a junction letter loses at most one packing element, since its 2-neighbourhood is all the other colour.
- Hence P(A) + P(B) >= ceil((L + 2J)/3).
- A tight colouring therefore has at most one junction, and only at n = 3k, where the slack is 2.

**Two sides of the run-window tiling R(w).**
- The 3-cycles of R(w) and of its complement fall into a few explicit types.
- The complement side always has nu = P(B). The proof is a certificate argument: disjoint triangles get 3-separated
  0-positions. A feedback set of size P(B) need not exist: for 010 the complement has nu = 1 < tau = 2.
- The own side has nu = P(A) whenever the word has no 0-junction. The proof is a greedy feedback set.

**The word law.** Use R(w), or, at a 0-junction, the line of R(complement w). The junction lemma says one of the two
is always safe. A palindrome gives a blue line.

**Central-pair lemma.** In an even-length palindrome, an odd P(X) forces both central positions into X.
- This turns into a theorem the S87 audit's finite pattern "blue tight pairs need j even or k odd" at even n. That
  audit had proved the n = 12 case by listing.
- The odd n = 3k + 1, k even case is the audit's central-run argument (HYP-9540 notes).

**Residual pairs (THM-5020 (6); the audit's upgrades).**
- The run-window lines of admissible words realise every allowed nu-pair for 6 <= n <= 40, except three families
  that no word or palindrome can carry. All three exclusions are proved.
- Explicit nested-mirror lines realise (k, k), black and blue, and the blue pairs (i, k) at every n.
- So HYP-9540 (B) holds for every n <= 26.
- The open residue is the blue both-odd pairs (i, j), j < k, at even n > 26.

**Checks.**
- The session: all words n <= 15, all tight colourings n <= 24, all even palindromes of length <= 24.
- The audit: all words n <= 22, all tight colourings n <= 36, all even palindromes of length <= 32.
- 0 failures.

## 2. Tight colourings as a finite automaton (THM-5021, PROVED)

The greedy packings of the two colours plus the running excess form an exact finite automaton (S87's greedy-state DP,
made finite by the prefix bound: a tight word never has prefix excess above 1).
- **States.** Of its 18 triples only 12 occur, and one block-boundary state is dead.
- **Spectrum.** The trim block automaton has 5 states. Its transfer matrix has characteristic polynomial
  x^2 (x - 2)^3, with a single Jordan block of size 3 at the eigenvalue 2. The first draft reported the 18-state
  artefact x^14 (x - 2)^4 with a block of size 4 (audit M1).
- **Generating functions.** In x^k:
  - n = 3k + 2: 2x/(1 - 2x);
  - n = 3k + 1: 2x(1 - x)/(1 - 2x)^2;
  - n = 3k: 2x(1 - x - x^2)/(1 - 2x)^3.
  The pole order is one more than the slack, and the maximal order 3 is the Jordan block's size. This proves the
  fitted count (k^2 + 9k - 2) 2^(k-3).
- **Junctions.** At n = 3k, (k - 2) 2^(k-3) tight colourings have a 0-junction, as many have a 1-junction, and none
  has both.
- **Markov chain.** The matrix counts words. It becomes a Markov chain after normalisation: uniform letters, T/8 per
  block, killed at E = 2.

## 3. The lenses the owner asked for (ANALOGY, typed; nothing below is used in a proof)

**Symmetry overlaid with asymmetry.**
- The run-window construction commutes with reversal (sigma), but not with complementation: R(complement w) = R(w) + 1
  only for the two constant words.
- That failure is the one asymmetric choice, which colour carries the staircases, and the colour of the junction
  forces it.
- A tight colouring has at most one junction, so a palindromic tight colouring can carry it only at its centre. If
  that central junction is a 0, even the blue line needs the complement construction (1110111 at n = 9), and the line
  remains sigma-symmetric.
- Among tight lines, then, the one possible junction of a blue line sits on the mirror. Short-block defects can sit
  off the mirror in symmetric pairs: 1100000011 at n = 12.
- The palindromic tight words number about the square root of all tight words (2^(k/2) against 2^k): symmetry is
  rare but never absent.

**Structure combined with randomness.**
- A tight colouring is rigid structure, the aligned blocks of three, carrying free random bits, one colour per block,
  plus a bounded number of defects: short blocks and at most one junction, paid for by a slack of 0, 1 or 2 units.
- In the transfer matrix this is the eigenvalue 2 (randomness: two choices per block) with a size-3 Jordan block
  (structure: the slack), so the counts are polynomial(k) x 2^k.
- Tao's "structure and randomness" dichotomy is the model; the repository's closest earlier statement is the
  quasirandom-arithmetic note (05-knowledge/results/quasirandom_arithmetic_structure_20261029/PROOFS.md).

**The Euler product as a chain of independent trials.**
- In ζ(s) = Π_p (1 - p^(-s))^(-1), the events "p divides N" are independent for a random N with P(N = n) ∝ n^(-s).
  The product is a chain of independent trials, one per prime.
- Two parallels here:
  - H_b = I(Omega(T), b) is multiplicative over the strong components of a tournament (THM-4940 (1)): an Euler
    product over "prime" (strong) factors.
  - The block colours of a tight colouring are independent trials.
- Where the analogy stops:
  - For pairwise disjoint cycles, I(Omega, lambda) = Π_C (1 + lambda). The gap from this product measures how much
    the cycles interact. The CONCEPT-MAP line had the exponent lambda^|C|, corrected this session.
  - Here the slack measures how far the block trials are from independent.

**Primes as a Markov chain.**
- **"Marty's Primes as a Markov Chain" was not found.** Neither this session's survey nor the audit found it in the
  repository or on the web. The owner is asked for the source.
- The closest works:
  - V. R. R. Uppuluri, "Prime Numbers Treated as a Markov Chain", Oak Ridge National Laboratory (1968): the last
    digits of successive primes as an ergodic Markov chain;
  - R. Lemke Oliver and K. Soundararajan, PNAS 113 (2016): biases between the last digits of consecutive primes.
- The audit adds a lead, not an identification: Alexeev, Barreto, Li, Lichtman, Price, Shah, Tang and Tao,
  "Primitive sets and von Mangoldt chains", arXiv:2605.00301 (2026), which builds Markov chains with von Mangoldt
  weights on the divisibility structure.
- The parallel used here is structural. The phase (position mod 3) drives an exact finite automaton with a single
  nonzero eigenvalue 2. The Ulam and Tao chains have that feature too (THM-4821/4822), but they are rank one plus
  nilpotent, while here the eigenvalue-2 part is a Jordan block of size 3. THM-4473 records exact Markov laws for
  Collatz digit chains and the decaying consecutive-prime digit bias.

**The pasted text (self-avoiding walks, honeycomb, Collatz), with corrections; archived in INSPIRATION_PASTE.md.**
- **The prime sieve as a self-avoiding walk with retrospective elimination.** The greedy 3-separated packing is a
  rule of this kind: each chosen position forbids the next two. It is the hard-core exclusion rule. The tight
  condition is "no wasted space", and the junction lemma says an isolated letter is a place where the exclusion zones
  of the two colours overlap without loss.
- **Growth.** Tight colourings grow like 2^(L/3) up to a polynomial, a "connective constant" of 2^(1/3) per position.
  This is numerology, not used.
- **Duminil-Copin and Smirnov.** The paste says the honeycomb constant mu = sqrt(2 + sqrt 2) was proved "using a
  concept called conformal invariance". In fact (Annals 175 (2012)) the proof uses a discrete parafermionic observable
  that satisfies half of the discrete Cauchy-Riemann equations. Conformal invariance is neither used nor proved, and
  convergence to SLE(8/3) is open.
- **Tao's Collatz theorem.** The paste says Tao "used hydrodynamics and logarithmic vector flows". In fact
  (arXiv:1909.03562; Forum Math. Pi 10 (2022) e12) the theorem is that almost all orbits, in logarithmic density,
  attain almost bounded values. The method uses the Syracuse random variables mod 3^n and a characteristic-function
  decay bound from a two-dimensional renewal process. It does not prove the conjecture. The repository already
  refuted an earlier pasted description of this method (collatz_mod6_20260922_w6_tao_minus_sheet_ladder.md).
- **Collatz on the honeycomb.** Session 140 (session140_collatz_hexagonal_geometry_20261102) found that the Collatz
  graph embeds in the honeycomb lattice only with unbounded routing lengths.

## 4. Pre-registration and outcomes

[PREREGISTRATION.md](PREREGISTRATION.md) asked whether run-window lines of admissible words realise every pair
allowed by HYP-9540 (B) at n = 13..40.
- Q1 and Q2 were REFUTED as registered. The gaps are exactly three families, each impossible for words or palindromes:
  - black (k, k) at n = 3k;
  - blue (k, k) at n = 0, 1 (mod 3);
  - blue both-odd pairs at even n.
- Lines with these pairs exist by other constructions; the audit proved most of them.
- The in-sample calibration was run only after the push (MISTAKE-745).

## Files

| file | content |
|---|---|
| `word_law_check.py`, `word_law_check.out` | brute-force check of THM-5020 (1)-(4): every word n <= 15, every tight colouring n <= 24 |
| `palindrome_parity_check.py`, `palindrome_parity_check.out` | THM-5020 (5), every even palindrome of length <= 24 |
| `palindrome_kk_check.py`, `palindrome_kk_check.out` | family (ii) / HYP-9590 (C): the k-packings S_j, and palindromic (k, k) at n = 0, 1 (mod 3) |
| `pairs_check.py`, `pairs_check.out`, `pairs_check_insample.out` | THM-5020 (6) and PREREGISTRATION Q1-Q3: pairs realised by admissible words, n = 6..40 |
| `markov_tight.py`, `markov_tight.out` | THM-5021 (2): the automaton, transfer matrix, generating functions |
| `trim_check.py`, `trim_check.out` | THM-5021 (1): used, dead and trim states; the trim spectrum x^2 (x - 2)^3 (audit M1) |
| `markov_junctions.py`, `markov_junctions.out` | THM-5021 (3): tight colourings with junctions at n = 3k |
| `INSPIRATION_PASTE.md` | the owner's pasted inspiration text, archived (audit m13) |
| `audit/` | the independent audit: report and all its programs and outputs |

Earlier tools are in `../ternary_metagraph_20261109/` (S87): `tight_lines.c`, `staircase_check.py` and
`run_window_theorem_check.py`.
