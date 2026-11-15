# Span 3 with reachable cones: the full literal alphabet

The open span-3 statement is HYP-9742 for arbitrary literal letters. A literal letter is the set of gaps 1, 2, 3 at
which a new vertex beats an older one. The claim is that every bottleneck-free span-3 word has H ≥ C·ρ^n with
ρ ≈ 1.77, above the floor 5^(n/3). Bottleneck-free means rank ≥ 2 in the pinned minimal 62-dimensional transfer.
- *What is already proved.* THM-5212 (E) covers words in the 13 syllables. The syllable bridge covers syllable runs
  inside literal words.
- *Status here.* The literal case is **not** settled. This note proves a reduction, extends the periodic evidence to
  periods of up to 12 letters, and pins down why a uniform-memory reachable-cone certificate cannot work.

## 1. The literal language (FINITE-EXACT)

**Letters.** Ranks of letter powers in the 62-dimensional transfer (`lit1.py`), and growth of letter runs
(`lit2.py`):

| letter (back arcs) | ranks of powers 1–6 | growth of a^k |
|---|---|---|
| ∅ | 13, 3, 1, 1, 1, 1 | 1 |
| {1} | 13, 4, 3, 3, 3, 3 | 1.83929 (tribonacci) |
| {2} | 17, 13, 10, 9, 8, 8 | 2.28829 |
| {1, 2} | 17, 13, 12, 12, 12, 12 | 2.56297 |
| {3} | 48, 38, 36, 35, 34, 34 | 2.68076 |
| {1, 3} | 48, 45, 45, 45, 45, 45 | 2.89281 |
| {2, 3} | 57, 56, 56, 56, 56, 56 | 3.11359 |
| {1, 2, 3} | 57, 54, 54, 54, 54, 54 | 3.32260 |

- *Only the empty letter is a bottleneck.* M_∅³ has rank ≤ 1 over ℚ (exact 2×2 minors), while M_∅² does not.
- *Runs of {1} are tribonacci.* H({1}^k) is exactly the tribonacci number T(k+1), for k = 3–16. This is the same
  growth as the sole-exit tower of earlier work.
- *No cheap single-letter runs.* Every nonempty letter run grows faster than ρ* = 1.7700005, so none is a
  counterexample.

**Reduction (PROVED).** A bottleneck-free word never contains 000. So, after a trailing run of at most two empty
letters, it splits uniquely into the 21-word prefix code {0^i a : i ≤ 2, a ∈ {1, …, 7}}. Literal span-3 words are
therefore words over a 21-letter code. The cheap behaviour comes from the floor syllables, which in code words read
like 02|05 or 01|04 (`verify_span3_code.py`).

## 2. Periodic words up to 12 letters: every cheap one is a bottleneck (VERIFIED modulo two primes)

Census of all 52,094 cyclic code words of at most 4 code units, which is all periodic literal words with at most
four nonempty letters per period, periods up to 12:
- 31,875 grow more slowly than 1.77 per vertex.
- For every one of them, the 64th power of the period matrix has rank ≤ 1 modulo both primes 33,554,393 and
  33,554,383 (`codecycles_exact.py`). **No cheap cycle in this range is bottleneck-free.**
- This extends THM-5202 (D), which checked all periods of at most 6 letters.
- The cheapest cycles carry exactly geometric counts, the signature of factorization, for example:

  | code word | letters | H(word^k), k = 1, 2, 3, … |
  |---|---|---|
  | 2\|05 | {2}, ∅, {1, 3} | 1, 5, 25, 125, … |
  | 1\|04 | {1}, ∅, {3} | 3·5^(k−1) |
  | 1\|1\|7\|4 | {1}, {1}, {1, 2, 3}, {3} | 15·9^(k−1) |

**A warning for numerics.** Floating-point eigenvalue ratios flagged 18 cheap "bottleneck-free" cycles, with
|λ₂/λ₁| between 10⁻⁶ and 10⁻³. All 18 are bottlenecks exactly; the ratios came from nilpotent blocks of these
62×62 integer matrices. Classify by modular or exact rank, never by floating eigenvalues.

## 3. Why uniform-memory reachable cones stall, and what would work

The reachable (monotone) cones of the 13-unit certificate are images of the nonnegative orthant of the 119-state
transfer under the last few blocks ([CONE13_REPLICATION.md](../chanpak_reverse_complement_collisions_20261115/CONE13_REPLICATION.md)).
They carry over to the literal code verbatim, because the code letters are products of literal matrices.
- *The obstruction.* A certificate at rate ρ needs every cycle of its automaton to grow at least ρ. With uniform
  memory of m code units, the automaton still contains floor-syllable powers whose short factors are
  bottleneck-free. A syllable needs up to 4 repetitions, up to 12 letters, before it becomes rank one. So
  uniform-memory cones can certify no more than the floor 5^(1/3) unless m ≥ 8 code units, and there are too many
  such states (about 10⁷).
- *What would work.* Memory that is adaptive along the 13 syllables, as in the concurrent adaptive certificate. The
  states would be (syllable-parse state, last two code units): long memory only while the word follows a floor
  syllable, short memory otherwise.
  - The census of §2 says the target rate 1.77 is consistent: no cheap periodic word escapes.
  - The 21-word code makes the parse state finite.
  - This construction and its LP are the precise next step. They are not executed here.
