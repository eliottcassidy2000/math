# Independent audit: THM-4823, PROBLEM_CLASSIFICATION §5, and the Σ⁰₄ / four-vertex reflection

*2026-10-27. An adversarial audit of the mac-mini-2026-10-27-slab pre-audit checkpoint by an independent subagent,
working from the statements rather than the producer's scripts. The auditor could not write report files from its
sandbox, so the producer saved this report verbatim from its final message. The scripts and outputs in this directory
are the auditor's own. Targets audited: THM-4823, PROBLEM_CLASSIFICATION §5 as of db9fb1c110, and the reflection,
including its later "First data" paragraph.*

## Verdict summary

| target | verdict |
|---|---|
| THM-4823 (A) census identities | **PASS**: proved, and verified on every class at orders 4–10. GAPs A1–A3 are wording. |
| THM-4823 (B) S4-free / TT4-free | **PASS** |
| THM-4823 (C) W4-free theorem and recursion | **PASS**: the proof is sound; the structure is checked class by class at orders ≤ 10 and on labelled tournaments at orders ≤ 7. GAPs C1–C2. |
| THM-4823 OEIS claims | **PASS**. Suggest adding W(n) = A133670(n) − n. |
| THM-4823 typed reading | PASS, with GAPs R1–R3 |
| PROBLEM_CLASSIFICATION §5 | **4 FAIL** (ETH level; full-Collatz horizon; JC(2) horizon; the Reading contradicts its own rows), about 15 GAPs |
| Reflection (Σ⁰₄ / four vertices) | **6 FAIL**, all local. The reverse-mathematics claims PASS (checked on the web) and the new five-vertex data PASS (reproduced). |

## 1. THM-4823

### (A) Census identities: PASS

**Proof re-derived.**
- A 4-set has at most one source. Sources occur only in TT4 and W4, so Σ C(s,3) = #TT4 + #W4. Sinks give the second
  identity the same way.
- In an oriented graph a closed 4-walk cannot repeat a vertex, so tr A⁴ = 4c₄.
- A directed 4-cycle spans a strong 4-set, which is S4. S4 has exactly one Hamilton cycle (1 of the 6 possible directed
  cycles).
- The solved formulas are correct.

**Own code (`census4.c`).**
- Each 4-set is classified twice, by sorted local scores and by the source/sink bits, and its 4-cycles are counted
  directly.
- Classes checked: 4, 12, 56, 456, 6,880, 191,536 and 9,733,056 at orders 4–10.
- Result: no failures of the identities or the solved formulas. tr A⁴ was computed from an explicit integer A².
- Score sequences on which #S4 is not constant: 0/4, 2/9, 11/22, 41/59, 136/167 at orders 4–8, exactly as claimed.
  Extension: 437/490 at order 9 and 1,397/1,486 at order 10.
- Witnesses at order 5: scores (1,1,2,3,3) occur with #S4 = 2 and with #S4 = 3 (`1101111101`, `1101111110`). Scores
  (1,2,2,2,3) occur with #S4 = 3 and 4.

**Gaps.**
- **A1.** "From order 5 on" is only computed to order 10. One line closes it for all n: add a sink tail to two order-5
  witnesses. The score sequences stay equal and #S4 is unchanged, because any 4-set containing a tail vertex has a sink.
- **A2.** "Exactly one matrix product" is only an upper bound. One product gives A², and then
  tr A⁴ = Σ (A²)ᵢⱼ(A²)ⱼᵢ costs O(n²). Showing that the scores do not determine #S4 does not prove a product is
  necessary. Note that c₃(v) = Σ_{u∈N⁺(v)} s_u − C(s_v,2) is O(n²). Write "at most one matrix product; no lower bound
  claimed".
- **A3.** The THM-1780 sentence ("matches … H spectral at orders 4, 5") is untyped inside a PROVED section. It compares
  different thresholds (5 vs 6) and different levels (scores vs spectrum); #S4 is spectral at every order. Type it
  ANALOGY or delete it.

### (B) PASS

**S4-free.**
- Every S4-free class at orders ≤ 10 has all strong components of size 1 or 3.
- Every class with a component of size ≥ 4 contains S4.
- Counts: 3, 4, 6, 9, 13 at orders 4–8, then 19 and 28 at orders 9–10 (A000930).
- Labelled counts follow n!·a(n) with a(n) = a(n−1) + a(n−3)/3: 40, 240, 1,760, 15,120, equal to brute force at
  orders 4–7.

**TT4-free.**
- Counts: 3, 3, 1, 1, 0 at orders 4–8.
- The unique order-7 class is doubly regular and isomorphic to QR₇. The unique order-6 class is P₇ − v.
- Labelled counts: 40, 184, 240, 240 at orders 4–7.

### (C) Theorem C: PASS

**Step 1 is sound.** N⁺(a) contains {b, x, y} with b → x (because C ⊆ N⁻(x)) and x → y. Transitivity of N⁺(a) forbids
y → b. The same argument in N⁺(b) gives c → y.

**Step 2 is sound.** Step 1 uses only C ⊆ N⁻(x) and W4-freeness, so it propagates along every directed path.

**Step 3 is sound.**
- R is nonempty.
- R is out-closed, so a sink of T[R] would be a sink of T.
- T[R] is transitive and nonempty, so it has a sink: contradiction.

**Equivalent form (not a gap).** W4-freeness makes every strong component after the first a singleton. So the theorem
is the same as "strong and W4-free implies locally transitive".

**Gaps.**
- **C1.** The count needs the decomposition to be unique. It is, because a tournament has at most one sink, so peeling
  sinks is canonical. Add one sentence.
- **C2.** The closed form also needs LT(n) = A000016(n) for all n. That is CITED (Brouwer 1980), and FINITE-EXACT only to
  order 10. The recursion itself is PROVED.

**Computation.**
- W4-free counts 3, 6, 11, 20, 35, 64, 115 and locally transitive counts 2, 4, 6, 10, 16, 30, 52 at orders 4–10. This
  reproduces `w4free.out`.
- The two definitions of W4-free (no W4 in the census; every out-neighbourhood transitive) agree on every class.
- Peeled-core check on every W4-free class:
  - the core is K₁ exactly once (TT_n);
  - otherwise it is sinkless, strong and locally transitive (0 violations);
  - the cores of order k number LT(k) − 1 = 1, 1, 3, 5, 9, 15, 29, 51 for k = 3..10.
- Every locally transitive class with a sink is transitive.
- Labelled check, without gentourng:
  - locally transitive counts are 48, 384, 3,840, 46,080 = 2^(n−1)(n−1)!;
  - W4-free counts are 56, 544, 6,384, 85,728, equal to W_lab(n) = n! + Σ_{k=3..n} C(n,k)(n−k)!(2^(k−1)(k−1)! − k!).

**OEIS.**
- A000016 (offset 0) matches 2, 2, 4, 6, 10, 16, 30, 52 at n = 3..10. The entry links Brouwer's "The Enumeration of
  Locally Transitive Tournaments" (Math. Centr. ZW138, 1980); no comment line states it, the link does.
- W(n) has no OEIS match, in the full window or the shorter ones.
- Add: W(n) = A133670(n) − n, where A133670 is the partial sums of A000016 (checked for n = 1..10). The labelled
  sequence 56, 544, 6,384, 85,728 is not in OEIS either.

**Side finding (THM-402).**
- The `is_round` test in `04-computation/lrc_round_2dichromatic_s592.py` checks only out-neighbourhoods, i.e.
  W4-freeness.
- Its "544 round tournaments at m = 5" is exactly W_lab(5) = 120 + 40 + 120 + 264. The locally transitive count at
  m = 5 is 384.
- THM-402's χ ≤ 2 still holds: on L ⇒ TT_j, append the tail to either transitive part.
- Its Verification line mislabels the tested set. That is outside the targets; correct it in THM-402, and THM-4823 can
  cite this.

### Typed reading

- **R1.** The NO TRANSFER bullet is about a #SAT-type counter, which counts labelled tournaments, but THM-4823 proves
  class counts. Add the labelled closed forms, which the auditor verified: W_lab above, S4-free = n!·a(n), and TT4-free
  (empty from order 8).
- **R2.** Cite the preprint with its tier, "T3, PREPRINT CLAIM", as TRUSTED-ADVANCES' one-tier rule requires.
- **R3.** "Non-recursive counts" literally means uncomputable, and "begin with 5-vertex patterns" is not proved. Rephrase
  it as a proposal.

**Producer code.** `pattern4.py` and `w4free.c` are correct.

## 2. PROBLEM_CLASSIFICATION §5

### 2.1 ETH level: FAIL

The same error appears in the reflection's title, §1 and §3, and upstream in klein's reflection §1.

- **Standard form.** Let Corr(A, φ, t) mean "A halts on φ within t steps with the correct answer"; it is decidable. Then
  ETH (s₃ > 0) is
  ∃δ∈ℚ>0 ∀A ∀c ∃φ ¬Corr(A, φ, c·|φ|^c·2^(δ·n(φ))),
  which is **Σ⁰₃**.
- **The formula as written.** The hypothesis "A total" is redundant: any algorithm meeting a computable time bound on
  every input is total. "∃x, A slow on x" has no constant c; it must be ∀c ∃x.
- **Negation.** ¬ETH is **Π⁰₃**: one algorithm for each δ. It needs no uniformity in δ and no proofs.
- **The "no 2^(o(n)) algorithm" form.** This form is Π⁰₄, and its negation (a single subexponential algorithm) is Σ⁰₄.
- **The preprint.** Its claim is one explicit algorithm with an explicit computable bound, which is **Σ⁰₂**.

No standard reading puts ETH at Σ⁰₄. As a sanity check from memory, P ≠ NP is Π⁰₂ by the same count.

### 2.2 Definitions (GAPs)

- **D1.** "One per residue class" contradicts b = O(states × symbols). FRACTRAN's modulus is the lcm of the denominators,
  which is exponential in the number of fractions; only k + 1 distinct affine rules appear. Define b as the number of
  distinct rules.
- **D2.** Klein's step 1 asked for the FRACTRAN constants, which are not computed here. The levels are also extended
  beyond klein's (I2 gains linear forms in logarithms; I3 gains analytic, geometric, ordinal). Say "adapted from klein".
- **D3.** The Sloane/Erdős digit tests, weak Collatz and Antihydra have no residue stop class. Say that for natural rows
  b counts the branches of the dynamics.
- **D4.** The Sierpiński/Riesel cell and Rule 30's "8 local rules" are not branch counts. Use "generic" or "not a loop".

### 2.3 Rows

| row | verdict | reason / correction |
|---|---|---|
| Goldbach | PASS | The I1 claim is true: every even N is a sum of two units mod m for every m (CRT argument; checked m ≤ 300). Label Helfgott and "density one" as from memory. |
| RH | GAP | Horizon is "—". Move "zeros on the line to a large height" into h (Platt–Trudgian 2021: 3·10¹², from memory). |
| Sloane / Erdős #406 | PASS | 2 is a primitive root mod 3^k (checked k ≤ 12). Zero-free units number 2^k, a share of (2/3)^(k−1), so density zero holds. Saye's 2·3⁴⁵ ≈ 5.9·10²¹ covers both conjectures. |
| weak Collatz | GAP | THM-4633 records more than 1.375·10¹¹ odd elements per nontrivial cycle (Hercher, with Barina's 2⁷¹). THM-4591's period 301,993 is far weaker. |
| full Collatz | **FAIL** | "≥ 2⁶⁸ (from memory)" contradicts canon THM-4633: 2⁷¹ (Barina, J. Supercomput. 81 (2025)). Listing Tao 2019 only at I1 understates it: it is a log-density result. THM-4821 and THM-4822 are pre-audit; say so. |
| Sierpiński / Riesel | GAP | "Open candidates are exactly those without a covering" is false as written. Intended: every open candidate lacks a covering. |
| FLT | GAP | Kummer's result is algebraic, not I1. |
| PA / ZF / ZFC | GAP | Gentzen's ε₀ consistency proof for PA is an ordinal (I3) partial result. |
| JC(2) / DC(2) | **FAIL** | Horizon "degrees 18, 22" ignores what the repo already records: degree ≤ 100 (Moh; HYP-9343, THM-3549) and ≤ 108 (Guccione et al.; MISTAKES.md, MISTAKE-606). |
| Mahler 3/2 | PASS (wording) | The block law is a forbidden-factor (I1) condition; "≤ ⌊ℓ/3⌋+1 ones" is its I2 consequence. |
| Hadwiger–Nelson | GAP | "Σ⁰₁ witnesses" is not a level word. |
| magic square of squares | GAP | "No congruence obstruction" has no in-repo source and is uninformative, since degenerate points exist. THM-4765 and THM-4773 say nothing about integer squares. |
| H ≥ disc | GAP | THM-4051 (audited) settles the strong base at order 7, so H ≥ disc is exhaustive through order 7. |
| black core | GAP | The connectivity pointer should be THM-4811 (E), not (D2). THM-4819 (D) and (E) are audit-pending. |
| ABC | GAP | Stewart–Yu (linear forms in logarithms, from memory) is I2 by §5's own definition. |
| Kakeya / AMM | GAP | A level cell "—" breaks the column rule. |
| LRC(14), twin primes, odd perfect, Antihydra, Rule 30 (level) | PASS | |
| new tournament row | PASS (wording) | Π⁰₁ and the horizon match THM-4819 (C). But L = ½Σ\|s₍ᵢ₎ − i\| and R = ½Σ\|s₍ᵢ₎ − m\| are ℓ¹ distances, not linear functionals; they become linear only on the polyhedral slab. Type "Parikh vector" as an analogy. |

### 2.4 Reading and coverage

- **FAIL.** "All b = 2 rows have I1 and I2 results" contradicts the weak-Collatz and Antihydra rows, which list I2 only.
- **GAP.** "None at I3" is contestable under §5's own definition (I3 = analytic): Tao 2019, Krasikov–Lagarias x^0.84, and
  Flatto–Lagarias–Pollington for Mahler (all from memory). Suggest "none at I3a, and no I3 result decides an orbit".
- **GAP.** The census-rows sentence misfits odd perfect (its row says structural) and RH (blank horizon).
- **GAP.** "No row can be resolved that way" conflicts with Sun 2-4-6-8, which §3 itself records as refuted by an
  explicit Σ⁰₁ witness.
- **GAP.** Rows are missing: GMC(2)/NC2 and Sun 2-4-6-8 are omitted while FLT is kept.
- **GAP.** Literature claims in the level column carry no "CITED from memory" label.

## 3. Reflection (Σ⁰₄ / four vertices / galactic algorithms)

**FAILs.**
- **ETH at Σ⁰₄; negation Π⁰₄ "uniform family with proofs".** See §2.1.
- **"Smallest nontrivial strong tournament".** That is C3, whose count is a function of the scores. Call S4 "the strong
  4-tournament".
- **"The mechanism is Lipschitz".**
  - The Lipschitz argument (4 per arc) gives only about q^(3/2)/4 edits.
  - Lipschitz plus Frobenius gives q²/8 (THM-4814 (C)).
  - The q²/4 bound needs trace duality: ⟨A,B⟩ = q(q−1) − 4e ≤ ‖A‖_op·‖B‖_*, with ‖A‖_op = √q (Session75 RESULT (3.2)).
  - At q = 100,003 the three bounds are 7.5·10⁶, 1.1·10⁹ and 2.4·10⁹.
- **"Rank is not Lipschitz".** False: one arc reversal changes rank by at most 2.
- **"The Klein-four grading is the first 2-adic digit".** H(W4) = H(L4), so H cannot carry the V4 labels.
  c₃ = (H−1)/2 is the coefficient i₁, not a binary digit (H(S4) = 5 = 101₂). H mod 4 sees only c₃ mod 2.
- **§3 "no 2^(δn) algorithm counts H … Σ⁰₄".** It is Π⁰₂ for fixed δ, Σ⁰₃ with ∃δ.
- **"A Collatz termination proof would have to supply an ordinal-type measure".** An untyped non-sequitur; analytic
  routes exist.

**GAPs.**
- **The fifteen links.** There are 14 bullets. THM-133 is the Z₇-circulant formula, not a P₇-only one. "Every one of
  the fifteen links lives at one of those two levels" is unsupported: only 6 are sorted, and the R(TT4)/P₇ extremality
  and the nuclear-norm maximizers are not census levels.
- **"H is #P-hard".** Not verified for tournaments (a web search found no reference). Label it from memory.
- **Labelled counts.** "Counted in polynomial time" holds for labelled counts too (formula above), but THM-4823 states
  class counts only. The pattern system is a 6-CNF with 8·C(n,4) clauses, which is not sparse; the preprint as recorded
  covers all k-CNF, so it still applies.
- **"One-line proof".** That is true only for counting occurrences in a given tournament. Counting pattern-free
  tournaments needs Theorem C.
- **"IDENTITY of cost".** It is an upper bound only.
- **"Two distinct metrics".** The metric is the same Hamming metric; the target sets differ (low-rank matrices vs parse
  tournaments).
- **Typing.** The EM "IDENTITY" should be an ANALOGY.
- **I3 split.** I3a/I3b drops "analytic, geometric".
- **"I2 decides only cycles".** It should say "bounds cycle lengths".

**PASSes.**
- **Nuclear norm.** The nuclear norm of 2(E_ij − E_ji) is 4 (singular values 2, 2).
- **Paley matrices.** The skew matrix of P_q is circulant and skew, with singular values √q (multiplicity q − 1) and
  one 0. Its nuclear norm is exactly (q−1)√q, and SSᵀ = qI − J. Checked for q = 7, 11, 19, 23, 31, 43, 47, 59; this is
  maximal at odd order (THM-4814 (C)).
- **Other §2 statements.**
  - Parse tournaments have nuclear norm O(n log n).
  - The q²/4 − o(q²) bound is right (Session75 spectral, independently audited).
  - The Dvir–Liu, Alman–Williams and Babai–Luks citations carry from-memory labels.
  - ω < 2.3714 is right; it is 2.371339 (from memory).
  - The skew-Hadamard doubling tower checks out to order 128 (nuclear norm N√(N−1)).
  - Odd order of automorphism groups needs Cauchy's theorem. THM-4812 does use Feit–Thompson.
  - Lemma B's failure list (m ∈ {5, 10, 13, 17, 25, …}) and (★) match THM-4819.
- **Hamilton paths at four vertices.** H = 1 + 2c₃ = 1, 3, 3, 5, checked on all 64 labelled 4-tournaments.
- **Counts.** A000568, and regular classes 1, 1, 3, 15, 1223 (gentourng).
- **Logical levels.** Δ⁰₀, Π⁰₁ and Π⁰₃/Σ⁰₃ are correct.
- **Reverse mathematics** (checked on the web).
  - RCA₀ ⊢ RT²₂ ↔ EM + ADS (Bovykin–Weiermann). Their original split used CAC, and Montalbán showed ADS suffices; some
    sources co-credit Hirschfeldt–Shore.
  - The Lerman–Solomon–Towsner abstract: EM is not equivalent to RT²₂.
  - EM is Π¹₂.
- **Paris–Harrington at fixed exponent** (from memory) and the **Cantor-space** remark are fine.
- **New "First data" paragraph.** Reproduced independently. At order 6 only TT5 and the (1,2,2,2,3) type with c₅ = 3,
  H = 15, |Aut| = 3 are determined. At orders 7–8 none are. At order 8 the pattern-free counts run from 139 to 2,971 (R5).
  - It holds in a stronger form: even sorted scores plus the full spectrum (all of tr A³…tr Aⁿ) determine no type at
    orders 7–8.
  - tr A³ in the list is redundant, since it is a function of the scores.

## 4. Corrections, consolidated

1. **ETH.** ETH is Σ⁰₃ and its negation Π⁰₃. The 2^(o(n)) form is Π⁰₄ / Σ⁰₄, and the preprint's claim is Σ⁰₂.
   "No 2^(δn) algorithm counts H" is Π⁰₂ for fixed δ, Σ⁰₃ with ∃δ. Flag klein's reflection §1 upstream.
2. **THM-4823.**
   - "Exactly one" → "at most one" matrix product.
   - Type the THM-1780 sentence as ANALOGY.
   - Add the padding line and the canonical-peeling sentence.
   - Status: the closed form is "PROVED modulo the CITED Brouwer count".
   - Add W(n) = A133670(n) − n and the labelled formula.
   - Add a note on THM-402's 544.
   - Name tier T3 at the preprint citation; reword "non-recursive".
3. **§5 horizons.**
   - Full Collatz: 2⁷¹ (THM-4633).
   - Weak Collatz: 1.375·10¹¹ odd elements (THM-4633).
   - JC(2): degree ≤ 100 / ≤ 108.
   - H ≥ disc: exhaustive through order 7 (THM-4051).
   - Black core: THM-4811 (E).
   - RH: put the zero height in h.
4. **§5 definitions and Reading.**
   - Definitions: b counts distinct rules; describe how natural rows stop; replace the two non-count cells; say
     "adapted from klein".
   - Reading: fix the I1 claim; replace "none at I3" with "none at I3a".
   - Re-level ABC, PA, Kakeya/AMM, Kummer, HN and the magic square.
   - Add the missing resolved rows; extend the "from memory" labels to the level column.
5. **Reflection.**
   - Rigidity: the mechanism is trace duality; delete "rank is not Lipschitz".
   - Rename "smallest nontrivial strong tournament" to "the strong 4-tournament".
   - Fix the Klein / 2-adic-digit sentence.
   - Type EM as ANALOGY; delete or type the ordinal-measure inference.
   - Label "#P-hard" from memory; weaken the fifteen-links claim.

## 5. Reproduction

Run inside this directory:
- `cc -O3 -o census4 census4.c`, then `gentourng -q n | ./census4 n` for n = 4..10 (about 35 s).
- `cc -O3 -o labeled4 labeled4.c && ./labeled4 7`
- `python3 checks_misc.py`
- `python3 pattern5_check.py 8`

Sources:
- Lerman–Solomon–Towsner, arXiv:1302.0828
- arXiv:1601.00050; arXiv:1601.04428; arXiv:1509.05340
- Patey, *Dominating the Erdős–Moser theorem*, arXiv:1505.03425
- Saye, arXiv:2202.13256
- #P-hardness search: arXiv:1506.00699 (no tournament hardness reference found)
- OEIS A000016, A133670
