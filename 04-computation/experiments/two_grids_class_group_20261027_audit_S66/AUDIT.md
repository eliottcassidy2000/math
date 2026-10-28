# AUDIT: opus-2026-10-27-S66 (THM-4825, THM-4826, HYP-9429, the results note, three programs)

*Saved by the author from the auditor's returned report text. The auditor's scripts and outputs are in this folder,
except the compiled `.exe` files and two arXiv full-text extractions (copyrighted papers, used only for reading), which
are not committed.*

This is an independent, adversarial audit of the worktree `math-wt-s66-eisenstein` at the pre-audit commit
**a29feb4ac6**, carried out on 2026-10-27.

**How it was done.**
- The worktree was treated as read only; nothing was edited, committed or pushed.
- Every claim was re-derived with the auditor's own code: Python 3.13 with numpy 2.3 and sympy 1.14, and C compiled
  with gcc 15 plus OpenMP on 12 threads. python-flint 0.9 is installed but was not needed.
- The author's programs were run only for the reproduction (item 6) and for one exact cross-count.

**The worktree changed during the audit, but not because of the auditor.**
- HEAD advanced to fd63e18127 through other sessions' commits. None of them touches the S66 deliverables.
- At 19:27–19:28 uncommitted edits appeared in four files: THM-4826, the note, and two_grids_class_group_20261027.py/.out.
  They were the author's own revisions.
- Everything below refers to a29feb4ac6. Where an edit already addresses a correction, it is marked "(done in the 19:27
  working copy)".
- The reproduction ran at 18:46, before those edits.

## Headline findings

**1. 6 and 54 are sums of four Eisenstein fourth powers.** This settles two of HYP-9429's seven open cases.

- **6.** 6 = (12 − 149ω)⁴ + (150 + 541ω)⁴ + (123 + 427ω)⁴ + (297 − 309ω)⁴.
  - Verified exactly in two ways: integer arithmetic in the basis (1, ω), and sympy with ω = e^{2πi/3}.
  - The bases fall in exactly the classes that THM-4826 (6) forces for n ≡ 6 (mod 16): fourth powers ≡ ω, 9ω, ω², 9ω².
  - The least box is |a|, |b| ≤ 541. The exhaustive search finds nothing with bases ≤ 540 and exactly this
    representation at 541.
  - With bases ≤ 800 the only representations, up to signs and order, are this one and its complex conjugate.
- **54.** 54 = 9·6 = (√−3)⁴·6 with √−3 = 1 + 2ω, so 54 = (310 + 173ω)⁴ + (−932 − 241ω)⁴ + (−731 − 181ω)⁴ +
  (915 + 903ω)⁴.
  - Verified exactly.
  - There is no representation with bases ≤ 800.
- **Consequences.**
  - Prediction (1) ("none for 6 with bases ≤ 150") stays CONFIRMED.
  - The open list shrinks to 60, 75, 78, 124 and 140. None of them has a representation with bases ≤ 650.

**2. Missed prior art inside a paper the theorem already cites.** Dechakulkamjorn–Rungtanapirom, arXiv:2208.02454
(2022), is cited by THM-4826 only for (3). It already proves:
- the Z[ω] case of THM-4826 (1), with the same unit-cycling argument (their Thm 4.1);
- the d = 7 case of (2) (their Thm 3.2(c));
- the Heegner d ≡ 3 (mod 8) cases of (2) (their Lemma 3.3, Lemma 3.5 and Thm 3.6).

**3. Everything else checks**, apart from a list of small errors and over-statements given below. All three programs
reproduce their committed outputs exactly.

## Summary table

| Item | Verdict | One-line reason |
|---|---|---|
| 1. THM-4825 (1)–(5) | SOUND WITH CORRECTIONS | Every statement and number was re-derived. (3) is false at d = 1, which the setting allows, and "the second grid G_d" misnames the union of the two grids. |
| 2. THM-4826 (1)–(6) | SOUND WITH CORRECTIONS | All mathematics and numbers check. (1)–(2) is largely prior art in the cited paper; the proof of (2) gives the ring-class-field degree "over K" where it means over H; the Candido/Proth attribution is inverted. |
| 3. HYP-9429 | SOUND WITH CORRECTIONS; 6 and 54 now SETTLED | Searches reproduced by three independent methods and extended. Pre-registrations precede their runs. `box` is exact and its class restriction valid. "No local obstruction" is true at every place (argument supplied). But 6 and 54 are representable; "six of the nine" multiples of 3 is seven; "47 needs large bases" is false; the local argument stops at 7; and the 2⁷ densities are not limits for nine n (the table is unaffected). |
| 4. Typing and novelty | CORRECTIONS NEEDED | THM-4826 (1)–(2) is PRIOR ART. Proth/Haldeman/Candido needs fixing. The Mersenne corollary is folklore. The design reading rests on a classical correspondence. No repository prior art exists for THM-4825, THM-4826 (4)–(6) or HYP-9429. |
| 5. Numbers | SOUND (one miscount) | Every quoted count, share, residue count, pair count, class size, density and rank was reproduced. The only wrong number is "six of the nine" (it is seven). |
| 6. Reproduction | REPRODUCED | Outputs are identical to the committed .out files apart from timing and line endings. `all 120` needs about 10 GB, and the obstruction script peaks at about 30–36 GB. |
| Results note | SOUND WITH CORRECTIONS | Inherits all of the above, plus four over-statements: "the Mersenne primes are exactly …", "there is no Gaussian quartic identity", "the first general bound … came from the two-grid ring", and "the second grid repairs every congruence obstruction". It also drops "squarefree" in the idoneal parenthesis. |

## 1. THM-4825

**What was checked** (`a1_thm4825.py`, `a1_classno.c`):

- **(1)–(2).** Every odd squarefree d ≤ 999 (204 with d ≡ 3, 200 with d ≡ 1 mod 4), in exact rational arithmetic.
  - G_d = Z + Zτ is closed under multiplication iff τ² ∈ G_d, iff d ≡ 3 (mod 4).
  - For d ≡ 1 (mod 4), the Z-span of {2, 2√−d, 1 + √−d, (1 + √−d)√−d} has Hermite basis {(1, 1), (0, 2)}, which is the
    checkerboard.
  - The products of that basis span 2Z[√−d], so 𝔭₂² = (2).
  - x² + dy² = 2 is solvable only for d = 1, and τ² ∉ G_d.
  - Symbolically, 2N(x + yτ) = 2x² + 2xy + ((1 + d)/2)y².
  - The proofs are correct.
- **(3).** h(−4d) for all 40,520 squarefree d ≡ 1 (mod 4) below 2·10⁵ (d = 1 included in that count).
  - Method: a complete count of reduced primitive forms, with no early exit.
  - Cross-check: Dirichlet's formula h(D) = Σ_{a ≤ |D|/2} χ_D(a) / (2 − χ_D(2)) for all 6,075 values 1 < d < 30,000.
  - Result: h = 2 exactly at d = 5, 13, 37. h = 1 occurs only at d = 1, and no odd h > 1 occurs, as [𝔭₂] forces.
  - The fundamental discriminants with |D| ≤ 20,000 and h(D) = 2 are exactly the 18 terms of OEIS A014603 (fetched).
    Only −20, −52 and −148 have the form −4d with d ≡ 1 (mod 4).
  - The reduced forms are confirmed, and so are the listed class numbers (17: 4, 29: 6, 41: 8, 53: 6, 61: 6, 85: 4,
    21: 4, 33: 4, 65: 8).
  - Baker (Ann. of Math. 94, 1971) and Stark (Math. Comp. 29 (1975) 289–302) are the right citations.
- **(4).** For d = 5, 13, 37, all split odd primes q < 2·10⁵ (8,962 + 8,939 + 8,938 = 26,839).
  - Each satisfies q = x² + dy² ⟺ q ≡ 1 (mod 4) ⟺ q = a² + b².
  - Otherwise q = 2x² + 2xy + ((1 + d)/2)y², tested independently as 2q = X² + dy² with X ≡ y (mod 2). No prime is
    represented by both forms.
  - The converse holds: for every squarefree d ≡ 1 (mod 4) with d ≤ 3,000 there is a split q ≡ 1 (mod 4) below 2·10⁵
    that is not x² + dy², except exactly for d ∈ {1, 5, 13, 37}.
  - All ten listed least counterexamples match.
  - The Bauer/Chebotarev argument is correct: K(i) lies in the genus field, hence in H, and equal split sets force
    H = K(i).
- **(5).** For d = 5, 13, 37, every n = 2^a d^b ≤ 10⁸ is represented by exactly one of the two forms, the first iff a is
  even. The owner's chain is as stated (5, 25, 100 on the first grid; 10, 50 on the second), and so are 13, 52 / 26,
  37, 148 / 74 and the table's witnesses.

**Findings:**

- **F1 (error, minor).** (3) is false at d = 1, which the setting allows. At d = 1, Cl = {[L₁]} equals
  {[L₁], [𝔭₂]} as a set, but h(−4) = 1. The same slip appears in the THM title and in the note's §2 Theorem.
- **F2 (wording).** "For d > 1, the second grid G_d = ½𝔭₂ represents an ideal class": G_d is the union of the two grids.
  The centred copy alone is a coset, not an ideal.
- **F3 (citation precision).** H(K) = K(i) comes from THM-4566's PROVED identification H(K) = F(i),
  F = Q(√p : p | m), read at the primes m = 5, 13, 37. The "owner-directive" completeness and the χ = 2 list (which
  includes 85, where H ≠ K(i)) are irrelevant to this step.
- **F4 (note §2).** {1, 5, 13, 37, 85} are the *squarefree* idoneal m ≡ 1 (mod 4) with no prime factor ≡ 3 (mod 4).
  Without "squarefree", the idoneal number 25 qualifies too (h(−100) = 2, forms (1, 0, 25) and (2, 2, 13)).
- **F5 (note §2, the d = 12 analogy).** The (½, ½) copy of Z[√−12] is not an ideal of Z[√−12]. As a lattice it is a
  proper ideal of its own multiplier ring, the conductor-8 order Z[√−48] (primitive form (4, −4, 13), discriminant
  −192).

## 2. THM-4826

**What was checked** (`a2_thm4826.py`):

- **(1).**
  - Z[ω], n ≤ 20,000: r(n) = 6Σχ₋₃. On odd shells the first grid carries r/3 = 2Σχ₋₃ and the centred grid 2r/3. Even
    shells have no centred points.
  - Hurwitz order, n ≤ 400: for odd n, Lipschitz 8σ(n) and centred 16σ(n); for even n, Lipschitz 24σ(odd part) and
    centred 0.
  - Both units of order 3 satisfy N(ζ − 1) = 3.
  - (1 + i)·(Hurwitz units) is exactly the set of 24 roots of D₄, under both left and right multiplication.
- **(2).**
  - **d ≡ 7 (mod 8).** For the 12 values d = 7, 15, …, 103 and odd n ≤ 3000, no centred element occurs.
  - **d ≡ 3 (mod 8) with h(−d) = 1** (d = 11, 19, 43, 67, 163), split p < 3·10⁵: every shell consists of 4 elements on
    one grid, the first iff p = x² + dy². Counts 4310/12963, 4321/12959, 4350/13011, 4313/12927, exactly the author's;
    for d = 163, 4278/13017.
  - **The ten d with h(−d) = 2** (35, 51, 91, 115, 123, 187, 235, 267, 403, 427, all ≡ 3 mod 8), split p < 10⁵: about
    half the prime shells are empty; the non-empty ones lie on one grid, with first-grid share 0.326–0.336 of the
    non-empty shells but about 0.16 of all split primes.
  - **Class-number ratio.** h(−4d)/h(−d) = 1, 3, 1 for d = 3, d ≡ 3 (mod 8), d ≡ 7 (mod 8): all 42 squarefree d ≡ 3
    (mod 4) below 200 and all 4,055 below 20,000.
- **(3).** All 39,231 primes p ≡ 1 (mod 3) below 10⁶, with primary primes computed independently by Euclid's algorithm in
  Z[ω] on gcd(p, r − ω).
  - b even ⟺ 2 is a cubic residue ⟺ p = x² + 27y²: 0 mismatches; 13,032 residues (0.3322).
  - The stronger form also holds for every p: (2/π)₃ equals π mod 2 in {1, ω, ω²}.
  - The two conjugate primary primes have b of the same parity, and their cubic symbols are conjugate, so the grid of
    "the" primary prime is well defined.
  - The p < 50 list is right.
  - arXiv:2208.02454 Thm 4.2 states exactly (1) ⟺ (2) ⟺ (3), with proof "[Ireland–Rosen, Prop. 9.6.1–2]" (read in the
    PDF).
  - OEIS A040028 cites Cox Thm 4.15 for the x² + 27y² criterion, conjectured by Euler and proved by Gauss.
- **Mersenne.** The ten listed representations (q = 5, …, 127) are exact. 2 = (2^m)³ with 3m ≡ 1 (mod q) is verified.
  For q ≤ 61 the primary prime recomputed by Euclid gives the same (x, y). q = 3 is correctly excluded.
- **(4).** Exact design strengths, from a moment test and a Gegenbauer test that agree:
  - {±1}: 1; {±ω, ±ω²}: 1; hexagon: 5; 16-cell: 3; tesseract: 3; 24-cell: 5; square: 3. The table is right.
  - Candido's identity, Liouville's identity, Σ_hexagon⟨u, x⟩⁴ = (9/4)|x|⁴ and Σ_24-cell⟨u, x⟩⁴ = 3|x|⁴ hold
    symbolically.
- **(5)–(6).** Over all α mod 16, α⁴ takes exactly the values {0} ∪ μ₃ ∪ 9μ₃, with α⁴ ≡ u or 9u where u = α⁴ (mod 2).
  μ₆⁴ = μ₃. The proof of (6) is correct, and a machine enumeration of all 210 class multisets mod 16 confirms that none
  sums to n ≡ 5 (mod 8).

**Findings:**

- **F6 (prior art, significant).** arXiv:2208.02454 (Dechakulkamjorn–Rungtanapirom, 2022) contains:
  - **Thm 4.1**, the Z[ω] case of (1): x² + 3y² counted as 2Σ(c/3) for odd n and 6Σ(c/3) for even n, with the same
    proof ({1, ω, ω²} is a reduced residue system mod 2, and exactly one of z, zω, zω² is ≡ 1);
  - **Thm 3.2(c)**, the d = 7 case: odd-norm elements lie in Z[√−7];
  - **Lemma 3.3**, for Heegner d ≥ 11: a split p is x² + dy² iff its prime factors are ≡ 1 (mod 2);
  - **Lemma 3.5**, for every squarefree d ≡ 3 (mod 8): even shells lie on the first grid;
  - **Thm 3.6**, the full count.

  The Hurwitz case of (1) is the classical step from Hurwitz's 24σ(n) to Jacobi's 8σ(n). "NEW TO THE REPOSITORY" is
  literally true, but the paper must be credited. This repeats the standing lesson about reading cited papers in full.
- **F7 (proof of (2)).** "Its degree over K is h(−4d)/h(−d) = 3" is wrong as stated. The ring class field L of Z[√−d]
  has [L : K] = h(−4d) = 3h(−d); the ratio h(−4d)/h(−d) is [L : H]. The step is correct only under the statement's
  hypothesis h(−d) = 1, which the proof does not restate, and the note's §3 drops that hypothesis altogether.
- **F8 (citation).** "Candido's identity (Candido 1951; the special case is older, Proth 1878)" has the history
  backwards. M. Ibrahim (arXiv:2001.08545, §8), citing Dickson's *History*, vol. II, pp. 650 and 657: the identity
  x⁴ + y⁴ + (x + y)⁴ = 2(x² + xy + y²)² itself was used by Proth in 1878 and is credited to C. B. Haldeman (1905).
  Candido (Scripta Math. 17, 1951) applied it to Fibonacci numbers, which is where the name comes from.
- **F9 (clarity).** "the primary prime over p": there are two primary primes over p, and both lie on the same grid.
  The proof of (3) shows only "first grid ⟹ x² + 27y²"; the converse is missing. Both done in the 19:27 working copy.
- **F10 (typing of (4)).** The correspondence between quartic Hilbert identities and antipodal 4-designs is classical
  (e.g. Reznick, Mem. AMS 96 (1992) no. 463).

## 3. HYP-9429

**Independent searches** (`a3_search.c`, `a3_search2.c`, `a3_brute.py`) used three methods:
- **`win`:** all pair sums in the box sorted in the lexicographic order of Z² (a translation-invariant total order); one
  two-pointer sweep covers the whole target interval, with no y-buckets or per-element binary search.
- **`deep`:** one target at a time. It enumerates the admissible class multisets mod 16 by machine (for n ≡ 6 (mod 16)
  exactly one, {ω, 9ω, ω², 9ω²}), splits each into two pairs, partitions the pair sums by residue mod 15 or mod 105, and
  matches the partitions class by class by hashing, radix sort and merge. Every hit is re-verified in 128-bit
  arithmetic.
- **numpy single-key brute force**, for B ≤ 30.

The three methods agree on every n in [0, 200] at B = 20 and at B = 30.

**Results:**
- **Window 0..50 at B = 20.** The misses are {5, 6, 13, 21, 29, 37, 45}, the same as the author's. Every represented
  n ≤ 50 has a representation with |a|, |b| ≤ 5 (30, 38, 44 need 5; 47 needs 4):
  47 = (4 + 2ω)⁴ + (2 − 2ω)⁴ + (4 + ω)⁴ + (3 − ω)⁴.
- **51..200.** At B = 60, 27 of 150 values have no representation (the 19 values n ≡ 5 (mod 8), plus 54, 60, 75, 78,
  102, 124, 140, 186), matching the author. At B = 120, 25 of 150 have none. The author's witnesses for 102 and 186 are
  recovered, and `deep` finds second ones:
  - 102 = (92 + 41ω)⁴ + (14 − 99ω)⁴ + (73 − 31ω)⁴ + (93 − 13ω)⁴;
  - 186 = (18 − 36ω)⁴ + (24 + 40ω)⁴ + (81 + 44ω)⁴ + (63 − 14ω)⁴.

  All quoted representations were verified symbolically.
- **The seven persistent misses.**

  | n | bases searched | result |
  |---|---|---|
  | 6 | ≤ 150, 300, 500, 540 | none |
  | 6 | 541 | representable (rep above) |
  | 6 | 800 | 2 hits (that rep and its conjugate) |
  | 54 | ≤ 800 | none |
  | 54 | 932 | representable (by scaling) |
  | 60, 75, 78, 124, 140 | ≤ 300, 500, 650 | none |

  The B = 650 and B = 800 runs each process 3.4·10¹⁰ to 5.6·10¹⁰ pair sums.
- **Cross-count with the author's `box`.** Identical representation counts for n = 22, 38, 70, 102, 118 at B = 150 (9,
  3, 2, 2, 5).
- **Calibration** (new).
  - At B = 200 the counts for the class-mates are small. In the class of 6: 22: 9, 38: 3, 70: 3, 86: 1, 102: 2, 118: 5,
    134: 1, 150: 2, 166: 1, 182: 2, 198: 11. The other three classes range from 1 to 30.
  - The counts barely grow from B = 200 to B = 500: 22 goes 9 → 10, 38 stays 3, 86 stays 1, 134 goes 1 → 2, 166 stays 1.
  - Non-obstructed misses in 51..200 fall from 17 (B = 20) to 15 (30), 8 (60), 6 (120), 6 (300) and 6 (500).
  - Persistence to large B is therefore weak evidence; 6 and 54 support the HYP's "sparse representation" reading.

**The author's C logic:**
- **`box B n`** searches exactly the box |a|, |b| ≤ B in the four classes: bases with (a even, b odd) and α⁴ ≡ ω or 9ω
  (mod 16); bases with (a, b odd) and α⁴ ≡ ω² or 9ω².
  - Signs are deduplicated canonically inside the symmetric box.
  - The restriction is a valid necessary condition: for n ≡ 6 (mod 16) the only admissible multiset is
    {ω, 9ω, ω², 9ω²}, and the parity filters follow from the residues. So 0 hits is exhaustive for the box.
  - Class sizes match (225/240/225/225 at B = 30; 5625/5700/5625/5625 at B = 150), and so do the pair counts.
  - Negative control: with the multisets for 14 or 12 (mod 16), 22 gets 0 hits instead of 3.
- **`all`** is exhaustive in the box. Robustness issue: `wit` is declared `long`, which is 32-bit on Windows. Pair
  indices ≥ 2³¹ (reached around B ≳ 175) would wrap and could produce a false "none". This is harmless for B ≤ 120 but
  should be `long long`.
- **`pat`** covers the conjugate box, as its own comment says. No final claim depends on it.

**Local obstructions** (`a4_local.py`):
- Residue counts match the author's at every modulus quoted: 7, 28, 100, 388, 1543 (2⁴–2⁸); 247, 2215, 19930 (3³–3⁵);
  151, 3751 (5², 5³); 484 (7²). Also computed: 6172 (2⁹), 28 (3²), 21904 (7³), 16 and 1600 (13, 13²), 31 and 3631
  (11, 11²).
- The only congruence obstruction for n ≤ 200 is n ≡ 5 (mod 8), at every 2-power modulus.
- **No obstruction at any place:**
  - at 2, every n ≢ 5 (mod 8) has a solution mod 2⁵ with a unit coordinate, which Hensel lifts;
  - at 3, every n has a unit-coordinate solution mod (1 − ω), which lifts;
  - at each prime 5 ≤ p ≤ 211, every n ≤ 200 has a solution mod p with a nonzero coordinate, which lifts;
  - for p > 211, Weil's bound gives N ≥ q³ − 81q^{3/2} > 0 because q ≥ p > 18.7, so a solution exists and lifts;
  - the archimedean place imposes no condition.
- **So the claim is TRUE but under-argued.** "6 and the eight are locally representable everywhere" holds, but the HYP
  and the note justify it only by tests at 2, 3, 5 and 7 plus the Chinese remainder theorem, which says nothing about
  primes ≥ 11.

**Local densities:**
- Every table value and product, the median 2.6667, the minimum 0.4444, the ranks 11, 12, 19, 42, 50, 66, 109 and 53,
  117, the lowest class {4, 20, 28, …, 196} and "6 of 7 below the median" reproduce exactly.
- But the 2⁷ values are not limits for nine n: 16, 32, 48, 64, 128, 144, 160, 176, 192 (e.g. σ₂(16) = 11, not 10;
  σ₂(128) = 12, not 12.25). The exact recursion is σ₂(16m) = σ₂ᵘ(16m) + σ₂(m), which needs modulus 2⁸ for n ≤ 200.
  Similarly σ₃(9m) = σ₃ᵘ(9m) + σ₃(m), so modulus 3⁵ is exact for v₃(n) ≤ 4. With the stable values, the median and every
  rank are unchanged.

**Pre-registration integrity:**
- 06ca4c5d48 (17:02:49) registers (1) and (2); their text is unchanged afterwards (diffed).
- d0a9ad6895 (17:09:34) appends the (1)/(2) results, registers (3), and first commits the C program.
- a29feb4ac6 (17:19:38) appends the (3) result. The .out containing the `all 120` run (102 s) appears only here.
- So each prediction precedes the commit that reports its run.
- Protocol notes: the registration names a "class-restricted" search, but (2) and (3) were run with the stronger
  exhaustive mode `all`, and the note was first committed after the runs; the front matter (title and status) was never
  updated after the runs.

**Further findings:**
- **F11.** "Six of the nine misses (6, 54, 60, 75, 78, 102)" are multiples of 3: 186 = 2·3·31 is one too, so it is seven
  of nine. The base rate among the candidates is 45/132 = 0.34; the hypergeometric tail is P(≥ 7) = 0.0075 for the nine,
  and 0.045 for 5 of the 7 persistent misses. This is post hoc. Multiples of 3 have the larger 3-adic density, and
  54 = 9·6 is a forced pair.
- **F12.** "Some need large cancelling bases: 47 …" is false: bases ≤ 4 suffice for 47.
- **F13.** The local argument is incomplete, as described above.
- **F14.** "6, 54, … are open" is now outdated.

## 4. Typing and novelty

**Repository grep** (01-canon, 05-knowledge, 07-reflections, excluding the S66 files):
- "x² + 27y²" and "2 is a cubic residue" appear only in the S40 note, which quotes Gauss's criterion (corrected under
  MISTAKE-626). Classical typing is right; a cross-reference would help.
- No prior statement of THM-4825, THM-4826 or HYP-9429 turns up for: Candido, Liouville's identity, g(4) ≤ 53; Hurwitz
  order, Lipschitz, 24-cell (only tournament numerology), spherical design; class number 2, (2, 1 + √−d); Eisenstein
  fourth powers.
- The concurrent THM-4824 (F) is credited.

**Status of each claim:**
- THM-4826 (1)–(2): PRIOR ART (F6).
- Mersenne corollary: elementary and folklore. B. Allombert (pari-users, 16 Nov 2024) states that 2 is a cube mod M_q
  for q > 3 and that M_q = 4x² + 27y².
- THM-4826 (4): a reading of a classical fact (F10).
- THM-4825: typing is right. (2, 1 + √−5) is the textbook non-principal ideal.
- THM-4826 (6) and HYP-9429: "new to the repository" is right. Two web searches found no literature.

**Attributions checked:** Baker 1971 ✓; Stark 1975 ✓; Candido 1951 ✓, but the Proth 1878 attribution refers to the
identity itself (F8); Liouville 1859 and g(4) ≤ 53 ✓; Balasubramanian–Deshouillers–Dress 1986 ✓;
Delsarte–Goethals–Seidel 1977 ✓; Jacobi ✓; Hurwitz 24σ(n) ✓; Cox Thm 4.15 ✓ and Thm 3.22 ✓ (Thms 7.24 and 9.4
consistent with the book as the auditor knows it); Ireland–Rosen Prop. 9.6.1–2 ✓.

**Over-statements in the note:**
- "the Mersenne primes are exactly the Eisenstein primes where 'doubling is a cube'" is false: 2 is a cube mod 43, which
  is not a Mersenne prime. (Done in the 19:27 working copy.)
- "the first general bound for fourth powers came from the four-dimensional two-grid ring" is anachronistic: Liouville's
  identity dates from 1859, the Hurwitz order from 1896.
- "there is no Gaussian quartic identity" is false. The square, together with the normalised centred points
  (±1 ± i)/√2, forms the octagon, a 7-design, which gives 6(x² + y²)² = 4x⁴ + 4y⁴ + (x + y)⁴ + (x − y)⁴ (verified).
- "the second grid repairs every congruence obstruction" overstates the second grid's role: the first grid alone already
  reaches 9–12 (mod 16). (Done in the 19:27 working copy.)

## 5. Numbers

Reproduced exactly:
- THM-4825: 40,520 squarefree d; 39 and 42 values below 200; h(−4d) = 2 only at 5, 13, 37; the 12 listed class
  numbers; the reduced forms; 8,962 / 8,939 / 8,938 / 26,839 split primes; all ten counterexamples.
- THM-4826: the shell laws; 4310/12963 (0.3325), 4321/12959 (0.3334), 4350/13011 (0.3343), 4313/12927 (0.3336);
  39,231 / 13,032 / 0.3322; the Mersenne table, the design table and the residues mod 16.
- HYP-9429: 841 / 7,321 / 29,041 fourth powers and 354,061 / 26,802,181 / 421,704,361 pair sums; 7 of 51, 27 of 150 and
  25 of 150 misses; 19 obstructed values; class sizes and the 32,062,500 × 31,640,625 pair-sum product; all densities,
  ranks and the median.

Wrong: "six of the nine" (it is seven). Not limits: the nine σ₂ values at 2⁷ listed in §3.

## 6. Reproduction

The C program was compiled (`gcc -O2`) and run in the modes `all 20 0 50`, `pat 30 22 1 9`, `box 30 22`, `box 30 38`,
`box 100 6`, `box 150 6`, `all 60 51 200` and `all 120 51 200`; both Python scripts were run from the worktree. After
normalising CRLF and timing lines, all three outputs are byte-identical to the committed .out files (374, 60 and 33
lines). Resource notes for the reproduction section: `all 120` allocates about 10 GB; the obstruction script peaks at
about 30–36 GB in the mod-243 step of `rational_residues`; the note's list of modes omits `box 30 38` and `box 100 6`.

## 7. Negative controls

- THM-4825: a perturbed second form (c = (3 + d)/2) breaks exclusivity for d = 5 at 20 of 34 split q; d = 17 fails the
  agreement at q = 13; the two class-number methods must agree on 6,075 values.
- THM-4826: Candido's identity with 3 in place of 2 and Liouville's with 5 in place of 6 are rejected; the 24-cell minus
  one point has strength 0; adding a false residue 5 to the fourth-power set makes 5 and 13 reachable.
- HYP-9429: the wrong class patterns give 0 hits for 22; a perturbed representation of 102 does not evaluate to 102; 13
  has no admissible multiset; with integer residues only, 5–15 become obstructed.

## Exact corrections (sentence → replacement)

*(Items 1–36 as returned by the auditor; the author's application of each is recorded in the note's §8 and in
MISTAKE-690.)*

THM-4825:
1. Title: insert "for d > 1" before "the two grids are the whole class group".
2. Statement (2): "the second grid G_d = ½𝔭₂ represents an ideal class of order exactly 2" → "the two-grid set G_d = ½𝔭₂
   (the first grid with its centred copy) is a fractional ideal whose class [𝔭₂] has order exactly 2; the centred copy
   alone is a coset, not an ideal".
3. Statement (3): add "d > 1" and "(At d = 1 the two classes coincide and h(−4) = 1.)".
4. Status: cite THM-4566's PROVED identification H(K) = F(i), F = Q(√p : p | m), at the primes m = 5, 13, 37.

THM-4826:
5. Status, novelty: credit arXiv:2208.02454 (Thm 4.1; Thm 3.2(c); Lemma 3.3, Lemma 3.5, Thm 3.6) as PRIOR ART for
   (1)–(2), the Hurwitz case as classical (24σ to 8σ), and the design reading to the classical correspondence (Reznick,
   Mem. AMS 463 (1992)).
6. Status, Candido: Candido (Scripta Math. 17 (1951)) applied it to Fibonacci numbers; the identity itself was used by
   Proth (1878) and is credited to C. B. Haldeman (1905) in Dickson's History, vol. II, pp. 650, 657.
7. Status, HYP pointer: 6 and 54 are representable (6 with bases ≤ 541; 54 = 9·6 by scaling with 1 + 2ω = √−3).
8. Statement (3): "a primary prime over p (either of the two conjugates; they lie on the same grid)". (Done.)
9. Proof of (2): degree h(−4d) = 3h(−d) over K, i.e. 3 over the Hilbert class field; density 1/3 under h(−d) = 1.
10. Proof of (3): add the converse. (Done.)

Results note:
11. Status: 6 and 54 representable (audit); 60, 75, 78, 124, 140 open (none with bases ≤ 650).
12. §2 Theorem: add "d > 1".
13. §2: cite THM-4566's PROVED identification at the primes m = 5, 13, 37.
14. §2: add "squarefree" to the idoneal parenthesis (else 25 qualifies).
15. §2, d = 12: "neither a ring nor an ideal of Z[√−12] (as a lattice it is a proper ideal only of its multiplier ring,
    the conductor-8 order Z[√−48])".
16. §3: ring class field of degree h(−4d) = 3h(−d) over K (3 over H); share 1/3 among split primes when h(−d) = 1, and
    in general among the non-empty prime shells (about 1/6 of all split primes when h(−d) = 2).
17. §4, Mersenne: "among", not "exactly". (Done.)
18. §5, Candido: the Proth/Haldeman/Candido history as in 6.
19. §5: replace the history sentence and "there is no Gaussian quartic identity" by: Liouville's identity (1859) is the
    4-design identity of the 24 roots of D₄ = (1 + i)·(Hurwitz units), proved algebraically before Hurwitz's order
    (1896); the square's units with the normalised centred points (±1 ± i)/√2 form the octagon, a 7-design, giving
    6(x² + y²)² = 4x⁴ + 4y⁴ + (x + y)⁴ + (x − y)⁴ (Liouville's identity at x₃ = x₄ = 0).
20. §5, mod-16 repair: first grid adds 9..12, centred grid adds 6, 7, 8, 14, 15. (Done.)
21. §5: "No other local obstruction", with the place-by-place argument (2, 3, 5 ≤ p ≤ 211, Weil above, C).
22. §5, 47: small bases suffice: 47 = (4 + 2ω)⁴ + (2 − 2ω)⁴ + (4 + ω)⁴ + (3 − ω)⁴.
23. §5, 6: no representation with bases ≤ 150; none ≤ 540; exactly one (up to symmetry) ≤ 800.
24. §5: 6 is representable in every completion of Q(ω).
25. §5, densities: exact except for n = 16, 32, 48, 64, 128, 144, 160, 176, 192 (2-adic factor needs 2⁸); median, ranks
    and lowest class unchanged.
26. §6 table: class-group row "(d > 1)"; split row PRIOR ART (arXiv:2208.02454); Mersenne row folklore (Allombert,
    pari-users 2024); designs row classical (Reznick 1992; Proth/Haldeman/Candido); search row updated.
27. §7: memory notes (about 10 GB, about 30 GB) and the modes `box 30 38`, `box 100 6`.

HYP-9429:
28. Front matter: status and title updated to OPEN for 60, 75, 78, 124, 140 with the registration commits and results.
29. 47 example as in 22.
30. Local argument for 6 as in 21/24.
31. Local argument for the eight as in 21.
32. Multiples of 3: seven of nine (6, 54, 60, 75, 78, 102, 186); base rate 45/132; P(≥ 7) = 0.0075 (5 of 7: 0.045);
    post hoc; multiples of 3 have the larger 3-adic density; 54 = 9·6.
33. Stability: the 2-adic density at 16m needs modulus 2⁸; the mod-128 values are not limits for the nine listed n;
    median and ranks unchanged.
34. "computed exactly" → exact except for the nine n.
35. Open question: 60, 75, 78, 124, 140 (6 and 54 are representable).
36. Protocol note: (2) and (3) were run with the stronger exhaustive mode `all`; the source note was first committed
    after the runs.

Also propagate items 5, 11, 26 and 35 to the INDEX entries.

## Optional strengthenings

1. **Local criterion as a theorem.** A rational n is a sum of four fourth powers in every completion of Q(ω) iff
   n ≢ 5 (mod 8): at 2 the unit-coordinate check mod 2⁵ covers all residues; at 3 always; at p ∤ n with q ≥ 19 Weil; at
   p | n with q > 27 Weil applied to x₁⁴ + x₂⁴ + x₃⁴ = −1; the remaining cases (p = 5, 7, 13, 19) by the direct check.
2. **Scaling.** Representability is preserved by multiplying by rational fourth powers of Z[ω] (16 = 2⁴,
   9 = (1 + 2ω)⁴), so misses are closed under division by 9 and by 16; this is how 54 follows from 6.
3. **THM-4825 (5)** holds for every prime d ≡ 1 (mod 4) (checked d = 17, 29, 41 to 10⁸).
4. **THM-4826 (2) for h(−d) > 1.** Split-prime shells are empty or on one grid; first-grid share 1/3 among the non-empty
   ones (checked for the ten d with h(−d) = 2).
5. **THM-4826 (3) with values:** (2/π)₃ ≡ π (mod 2) for all 39,231 primes. (Done.)
6. **Critical-dimension recursion:** σ₂(16m) = σ₂ᵘ(16m) + σ₂(m) and σ₃(9m) = σ₃ᵘ(9m) + σ₃(m) (Jacobian 1), the local
   face of the c(n) log B heuristic.
7. **Calibrate future pre-registrations** on class-mates and on the fill-in curve; prediction (3) failed because the
   survivors are selected for small c(n).
8. **Mersenne:** x is even, so M_q = 4x′² + 27y².
9. **Robustness:** change `wit` to `long long` in the author's C program.

## Files (this folder)

THM-4825: `a1_classno.c`, `a1_classno_sweep.out`, `a1_classno_ratio200.out`, `a1_classno_ratio20000.out`,
`a1_classno_listed.out`, `a1_thm4825.py` (+ .out). THM-4826: `a2_thm4826.py` (+ .out). HYP-9429 searches:
`a3_search.c`, `a3_search2.c`, `a3_brute.py` (+ .out), `a3_reps_check.py` (+ .out), `a3_six_verify.py` (+ .out),
`a3_54_verify.py` (+ .out), `a3_win20.out`, `a3_win60.out`, `a3_win120.out`, `a3_deep_controls.out`, `a3_deep150.out`,
`a3_deep300_a.out`, `a3_deep300_b.out`, `a3_deep500.out`, `a3_deep650.out`, `a3_deep800.out`, `a3_six_minbox.out`,
`a3_calib200.out`, `a3_calib500.out`, `a3_crosscount150.out`, `a3_negative_control.out`, `a3_classes_all.out`. Local
obstructions and densities: `a4_local.py` (+ .out), `a4b_recursion.py` (+ .out), `a5_mult3.py` (+ .out). Reproduction:
`run_repro.sh`, `compare_repro.sh`, `repro_e4.out`, `repro_two_grids.out`, `repro_obstructions.out`, `repro_summary.out`,
`repro_compare.out`. Batch scripts: `run_deep500.sh`, `run_calib.sh`, `run_calib500.sh`, `run_deep800.sh`,
`run_six_minbox.sh`, `run_deep650.sh`. Pre-registration snapshots: `hyp_v1.md`, `hyp_v2.md`. Other: `x0_pdftext.py`,
`omp_test.c`.
