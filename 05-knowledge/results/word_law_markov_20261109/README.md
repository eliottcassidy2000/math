# The word law for tight lines at n = 3k, and tight colourings as a Markov chain

opus-2026-11-09-S88. Canon: [THM-5020](../../../01-canon/theorems/THM-5020-word-law-for-tight-lines-at-every-order.md),
[THM-5021](../../../01-canon/theorems/THM-5021-tight-colourings-as-a-markov-chain.md),
[HYP-9590](../../hypotheses/HYP-9590-run-window-gap-law-and-residual-pairs.md); updates
[HYP-9540](../../hypotheses/HYP-9540-packing-number-line-laws-tight-words-allowed-pairs.md) ((A) and (C) now PROVED
at every n).

**Owner directive (2026-11-09).** "Start another session proving the word law at n = 3k, consider blue and black
lines being symmetry overlayed with asymmetry and primes being structure combined with randomness and seeing how the
deep structure behind the euler product formula as a chain of independent trials and marty's 'Primes as a Markov
Chain' intuition", with a pasted text on self-avoiding walks, the hexagonal lattice and Collatz as inspiration.

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
  0-positions. No feedback set can work here: for 010 the complement has nu = 1 < tau = 2.
- The own side has nu = P(A) whenever the word has no 0-junction. The proof is a greedy feedback set.

**The word law.** Use R(w), or, at a 0-junction, the line of R(complement w). The junction lemma says one of the two
is always safe. A palindrome gives a blue line.

**Central-pair lemma.** In an even-length palindrome, an odd P(X) forces both central positions into X. This
explains why no blue (1, 3) tight line exists at n = 12, which S87 had only observed.

**Checks.** All words n <= 15, all tight colourings n <= 24, all even palindromes of length <= 24: 0 failures.

## 2. Tight colourings as a Markov chain (THM-5021, PROVED)

The greedy packings of the two colours, plus the running excess, form an exact finite automaton with 18 states. A
tight word never has prefix excess above 1.
- **Spectrum.** The block transfer matrix (three letters) has characteristic polynomial x^14 (x - 2)^4, with a single
  Jordan block of size 4 at the eigenvalue 2.
- **Generating functions.** In x^k:
  - n = 3k + 2: 2x/(1 - 2x);
  - n = 3k + 1: 2x(1 - x)/(1 - 2x)^2;
  - n = 3k: 2x(1 - x - x^2)/(1 - 2x)^3.
  The pole order is one more than the slack. This proves the fitted count (k^2 + 9k - 2) 2^(k-3).
- **Junctions.** At n = 3k, (k - 2) 2^(k-3) tight colourings have a 0-junction, as many have a 1-junction, and none has
  both.

## 3. The lenses the owner asked for (ANALOGY, typed; nothing below is used in a proof)

**Symmetry overlaid with asymmetry.**
- The run-window construction commutes with reversal (sigma) and with complementation. It makes exactly one
  asymmetric choice, which colour carries the staircases, and that choice is forced by the colour of the junction.
- A tight colouring has at most one junction, so a palindromic one can carry it only at its centre. Blue lines are
  therefore the lines whose single possible defect sits on the mirror.
- Black lines are the same construction with the defect off-centre.
- The palindromic tight words number about the square root of all tight words (2^(k/2) against 2^k): symmetry is
  rare but never absent.

**Structure combined with randomness.**
- A tight colouring is rigid structure, the aligned blocks of three, carrying free random bits, one colour per block,
  plus a bounded number of defects: short blocks and at most one junction, paid for by a slack of 0, 1 or 2 units.
- In the transfer matrix this is the eigenvalue 2 (randomness: two choices per block) with a nilpotent Jordan tail
  (structure: the defects), so the counts are polynomial(k) x 2^k.
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
  - The product over cycles of (1 + lambda^|C|) equals I(Omega, lambda) only if cycles never interact
    (00-navigation/CONCEPT-MAP.md); the gap measures cycle interaction.
  - Here the slack measures how far the block trials are from independent.

**Primes as a Markov chain.**
- **"Marty's Primes as a Markov Chain" was not found.** The repository has no such text. A web search turned up
  nothing by anyone named Marty. The closest works are:
  - V. R. R. Uppuluri, "Prime Numbers Treated as a Markov Chain", Oak Ridge National Laboratory (1968): the last
    digits of successive primes as an ergodic Markov chain;
  - R. Lemke Oliver and K. Soundararajan, PNAS 113 (2016): biases between the last digits of consecutive primes.
  The owner is asked for the source.
- The parallel used here is structural. The phase (position mod 3) drives an exact finite Markov chain whose only
  nonzero eigenvalue is 2, just as the Ulam and Tao chains are rank one plus nilpotent (THM-4821/4822). The repository
  records exact Markov laws for Collatz digit chains, and the decaying consecutive-prime digit bias (THM-4473).

**The pasted text (self-avoiding walks, honeycomb, Collatz), with corrections.**
- **The prime sieve as a self-avoiding walk with retrospective elimination.** The greedy 3-separated packing here is
  literally such a rule: each chosen position forbids the next two. It is the hard-core exclusion rule. The tight
  condition is "no wasted space", and the junction lemma says an isolated letter is a place where the exclusion
  zones of the two colours overlap without loss.
- **Growth.** Tight colourings grow like 2^(L/3) up to a polynomial, a "connective constant" of 2^(1/3) per position.
  This is numerology, not used.
- **Duminil-Copin and Smirnov** (Annals 175 (2012)) proved mu = sqrt(2 + sqrt 2) for the honeycomb lattice with a
  discrete parafermionic observable that satisfies half of the discrete Cauchy-Riemann equations. Conformal
  invariance is neither used nor proved; convergence to SLE(8/3) is open. The pasted text is wrong on this point.
- **Tao's Collatz theorem** (arXiv:1909.03562; Forum Math. Pi 10 (2022) e12) proves that almost all orbits (in
  logarithmic density) attain almost bounded values. The method uses the Syracuse random variables mod 3^n and a
  characteristic-function decay bound from a two-dimensional renewal process. It does not use hydrodynamics or "vector
  flows", and it does not prove the conjecture. The repository already refuted an earlier pasted description of this
  method (collatz_mod6_20260922_w6_tao_minus_sheet_ladder.md).
- **Collatz on the honeycomb.** Session 140 (session140_collatz_hexagonal_geometry_20261102) found that the Collatz
  graph embeds in the honeycomb lattice only with unbounded routing lengths.

## 4. Pre-registration and outcomes

[PREREGISTRATION.md](PREREGISTRATION.md) asked whether run-window lines of admissible words realise every pair
allowed by HYP-9540 (B) at n = 13..40.
- Q1 and Q2 were REFUTED as registered. The gaps are exactly three families:
  - black (k, k) at n = 3k, impossible for words;
  - blue (k, k) at n = 0, 1 (mod 3);
  - blue both-odd pairs at even n, impossible for palindromes by the central-pair lemma.
- Lines with these pairs exist by other constructions at small n. HYP-9590 records the gap law and the residual
  existence question.
- The in-sample calibration was run only after the push (MISTAKE-745).

## Files

| file | content |
|---|---|
| `word_law_check.py`, `word_law_check.out` | brute-force check of THM-5020 (1)-(4): every word n <= 15, every tight colouring n <= 24 |
| `palindrome_parity_check.py`, `palindrome_parity_check.out` | THM-5020 (5), every even palindrome of length <= 24 |
| `pairs_check.py`, `pairs_check.out`, `pairs_check_insample.out` | THM-5020 (6) and PREREGISTRATION Q1-Q3: pairs realised by admissible words, n = 6..40 |
| `markov_tight.py`, `markov_tight.out` | THM-5021 (1), (2): the automaton, transfer matrix, Jordan form, generating functions |
| `markov_junctions.py`, `markov_junctions.out` | THM-5021 (3): tight colourings with junctions at n = 3k |

Earlier tools are in `../ternary_metagraph_20261109/` (S87): `tight_lines.c`, `staircase_check.py` and
`run_window_theorem_check.py`.
