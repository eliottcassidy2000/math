# Independent adversarial audit of opus-2026-10-26-S60 (aliquot table merge note, THM-4805, HYP-9415)

Auditor: independent subagent, own code only (no session script was executed). Worktree audited: `/tmp/math-wt-s60`
at HEAD `54617edb71` (the commit `41713d3aa2` under audit plus the follow-up `54617edb71`, which landed during the
audit and added the `10^6` twin test; both are covered). One uncommitted working-tree edit to the note (Theorem D'(ii),
last sentence) was also read; see Minor 13.

## Verdict: PASS WITH CORRECTIONS

Every number, set, list and proof in the note, THM-4805 and HYP-9415 was recomputed independently and reproduces
(table below), including a full brute force of `s^-1(8128)` to `8128^2` without the parity lemma and the complete
fibre of `33550336`. The pre-registration story is supported by the git history (predictions in the draft commit
authored 14:09:47 with the C++ source but no C++ output; outcomes in the commit authored 14:18:37; the `10^6` range
was in HEAD `41713d3aa2` before the `54617edb71` run). The hypotheses INDEX has 120 lines before and after.

What must change is typing/attribution and wording, not mathematics: Theorem D'(ii) is the textbook abc example and
must be typed classical; Theorem A'(ii) is contained in S57 (and aliquot124 explicitly declines to treat it as new);
two `related:` paths in THM-4805 are broken; the "tau <= 3" label is a misnomer for "three proper divisors";
C'(5) `D(1) = P(2)` is off by the element 1.

## Majors (type or attribution changes; no statement is mathematically false)

**M1. Theorem D'(ii) is classical and must be typed so (note sec. 2 and sec. 5; THM-4805 title/status; results INDEX
row; PROBLEM-LEDGER).** Quoted: "(ii) `(1, 2^k - 1, 2^k)` is an abc hit (`rad(abc) < c`) iff `2^k - 1` is not squarefree
... if squarefree, `rad = 2(2^k - 1) >= 2^k`; otherwise `rad(2^k - 1) <= (2^k - 1)/3 < 2^(k-1)`" and the typing table row
"Mersenne abc hits = non-squarefree Mersenne numbers; the list to 64 | PROVED (D') + FINITE-EXACT", under the heading
"What this note adds". The family `a = 1, b = 2^(6n) - 1, c = 2^(6n)` with `9 | b` and `rad(abc) < (2/3) c` is the
standard textbook demonstration that `rad(abc) < c` happens infinitely often (it is the worked example in the
Wikipedia article "abc conjecture", section on triples with small radical, fetched 2026-10-26; the note's bound
`rad(2^k - 1) <= (2^k - 1)/3` is that computation verbatim), and `p^2 | 2^(p(p-1)) - 1` gives the same for every odd
prime. The "only if" direction and the `ord_(p^2)(2) | k` rewording are immediate. Correct typing: "CLASSICAL
(textbook family `1 + (2^(6n) - 1) = 2^(6n)`; the iff is immediate); FINITE-EXACT: the list `k <= 64`". Remove D'(ii)
from "what this note adds" / "Added:" in the results INDEX and the ledger, and from the sentence "No earlier note
contained Theorems A'-D'" (sec. 6) add "D'(ii) is classical outside the repository". This is exactly the
MISTAKE-614/W1 pattern the repository already logs (a classical-looking equivalence typed as new).

**M2. Theorem A'(ii) is already in S57 and must be CITED, not "added" (note sec. 2, sec. 6; results INDEX "Added: ...
3 and 7 are the only Mersenne values with a single source"; PROBLEM-LEDGER "3 and 7 the only single-source Mersenne
values").** S57 sec. 3 states "Below `10^8` exactly five odd numbers have a single source of any kind: `3 = s(4)`,
`7 = s(8)`, `9 = s(15)`, `11 = s(21)` and `39 = s(217)`" and, one paragraph earlier, the mechanism "Odd `k >= 9` is
`s(pq) = 1 + p + q` whenever `k - 1` is a sum of two distinct primes" (conditional on distinct Goldbach). A'(ii) is the
restriction of these two sentences to `k = 2^k - 1`; the note's own proof reads "(ii) The table ... and S57's list.
QED". The same-day packet aliquot124 (`aliquot/README.md` line 23) writes of these facts: "its MISTAKE-672 correction
preserves the small single-source exceptions 3 and 7. We do not treat those statements as new." Correct typing:
"CITED (S57 sec. 3; MISTAKE-672) + corollary: among Mersenne values the single-source ones are exactly `3, 7`
(FINITE-EXACT to `k = 26` from S57's list, extended here to a distinct-prime partition of `2^k - 2` for every
`4 <= k <= 80`; conditional beyond)". The sentence "No earlier note contained Theorems A'-D'" is false for A'(ii).

## Minors

1. **"tau <= 3" is a misnomer** (note title, status line, Theorem A' heading, typing table; THM-4805 title and
   heading; PROBLEM-LEDGER "the tau <= 3 regime below 15"). The note defines `t = tau - 1`; the shapes `pq` and `p^3`
   have `tau = 4`. Write "`t <= 3` (at most three proper divisors)" or "`tau <= 4`".
2. **A'(i) quantifier**: "Every value `v <= 14` has only sources of the three shapes" is false for `v = 0` (source 1)
   and `v = 1` (the primes). Write `2 <= v <= 14`. (THM-4805 title has the same phrase.)
3. **C'(5) is off by the element 1**: "`D(1) = P(2)`" is false under the convention `k >= 1` that the rest of C'
   needs (`D(1) = {1, 2, 4, ...}`, `P(2) = {2, 4, 8, ...}`); with `k >= 0` every column would contain 1 and C'(1)
   would fail. Recomputed: the only row equal to a column for `n, m <= 200` is `D(2) = P(2)`, and
   `D(1) = P(2) ∪ {1}`. Also "2 is the only prime on a row" needs "on the row of a smaller number" (every prime lies
   on its own row). Fix: "the row of 2, i.e. the row of 1 without its root, is the only row that is a column".
4. **B'(i) title phrasing** (THM-4805 title; note sec. 2): "the sources of `P + 1` of the row form `2^a q` ... are
   exactly the twin source `2(P - 2)`" — add "when `P - 2` is prime, and there are none otherwise" (the body says this).
5. **Two broken YAML paths in THM-4805 `related:`**:
   `01-canon/theorems/THM-4790-presentation-count-fibres-7-and-21-are-absent-for-every-tournament-on-at-most-11-vertices.md`
   (actual: `THM-4790-presentation-fibres-7-and-21-absent-through-eleven-vertices.md`) and
   `01-canon/theorems/THM-3833-abc-conditional-cube-radical-growth-and-hyperbolic-power-finiteness.md`
   (actual: `THM-3833-abc-conditional-cube-radical-and-hyperbolic-power-finiteness.md`). All 20 markdown links in the
   note, all 6 in HYP-9415 and the INDEX/ledger links resolve.
6. **"Booker's `4 x 10^18` verification"** (note sec. 1, `{2, 5}` row): the Goldbach verification to `4·10^18` is
   Oliveira e Silva–Herzog–Pardi (2014), as S57 cites; session69 cites Booker (arXiv:1610.07471, footnote 2) only for
   the distinct-prime reading of that computation. Write "the Oliveira e Silva–Herzog–Pardi verification, read for
   distinct primes through Booker's footnote".
7. **abc census boundary** (note sec. 2 D'(iii), sec. 4; THM-4805): 121 and 419 are the counts for `c <= 10^4` and
   `c <= 10^5`; the standard tabulation (OEIS A216370 "Number of ABC triples with quality q > 1 and c < 10^n":
   1, 6, 31, 120, 418, 1268, ...) uses strict inequality and reads 120 and 418. Recomputed: the boundary hits are
   `2197 + 7803 = 10^4` (`13^3 + 3^3·17^2`) and `19 + 99981 = 10^5` (`19 + 3^3·7·23^2`). State the convention and cite
   A216370 so the counts do not look inconsistent with the literature.
8. **OEIS citations missing**: the source-count list for `v = 2..60` is OEIS A048138 term by term (recomputed and
   compared: identical), `a(v) = 0` iff untouchable (A005114), `a(v) = 1` iff `v` is a "hermit number" (A057709);
   `28 ∈ A057709`, so "28 is its own only aliquot source" is an OEIS-level fact and HYP-9415 (2)/(4) should be phrased
   "28 is the only perfect hermit number (A057709)".
9. **HYP-9415 (2)**: "Falsified by a perfect number with no companion, or any non-perfect fixed point" — the second
   clause is vacuous (a fixed point of `s` is perfect by definition). Delete it.
10. **Quotation source** (note sec. 2 B'(ii)): the quoted sentence "another perfect target such as 28 is not excluded
    from having other source shapes" is from the session69 arithmetic-corollary audit
    (`session69-arithmetic-corollary-audit/ROOT_ADDENDUM.md`, line 40: "... by this theorem alone"), not from the
    theorem text; name the addendum.
11. **Section pointer** (note sec. 1, "extra rare" row): the twin identity and the complete small fibres are in
    session69 aliquot's "Immediate conclusions" items 1-3 (and INCOMING_S69.md of aliquot124), the bound in sec. 2;
    "sec. 2-3" is imprecise.
12. **Lemma 1, `t >= 6`**: `s >= 22` is valid but the true minimum is `s(24) = 36`; optionally sharpen (the argument
    "six distinct proper divisors sum to at least 21, equality forcing `60 | n`" is correct).
13. **Theorem D'(ii), last sentence, HEAD `54617edb71`**: "session69's square-later carrier has its first excess at the
    same triple `1 + 8 = 9` read from `n = 72`, and the row of 1 reaches that triple only through `s(8) = 7`, `7 + 1 = 8`,
    `8 + 1 = 9`" is muddled: `(1, 8, 9)` is not a row-of-1 triple (`9` is not a power of two). The uncommitted
    working-tree edit replaces it with "this is a different family ... a triple with an odd `c`", which is correct;
    commit it.
14. **Draft-only, no action**: the draft HYP-9415 (commit `e68e0ea0e7`) wrote "`[0.56, 2.68]` below `10^4`", mixing the
    exactly-two share below `10^4` (115/205 = 0.561) with the mean over `(10^3, 10^4]` (2.68; the mean over all 205
    values is 2.634). The final HYP no longer contains the parenthetical.
15. **Theorem A'(ii) edge case**: with `k = 0`, `s(1) = 0` also has the single source `2^0`; say `k >= 1` (or `k >= 2`).

## Recomputation table (auditor's own code; claimed / recomputed / match)

| item | claimed | recomputed | match |
|---|---|---|---|
| bound `n < v^2` for composite sources | stated | proved (least prime `p <= sqrt n`, `n/p` proper) and checked to `10^6`, 0 violations | yes |
| source table `v = 0..60` | session `.out` | identical for all 61 values; counts `v = 2..60` = OEIS A048138 | yes |
| explicit fibres 3, 6, 7, 15, 21, 28, 31, 63 | {4}, {6,25}, {8}, {16,33}, {18,51,91}, {28}, {32,125,161,209,221}, {64,177,817} | same | yes |
| sources of 127 | eleven | 11: 128, 1469, 1853, 2033, 2369, 2813, 3293, 3569, 3713, 3869, 3953 | yes |
| untouchables `<= 1000` | 89; primes {2, 5}; odd {5} | 89 (list identical to A005114); {2, 5}; {5} | yes |
| Lemma 1 minima | `t=4: 15` (16), `t=5: 16` (12), `t>=6: >= 22` | 15 at 16; 16 at 12; `t>=6` minimum 36 at 24 (bound 22 valid) | yes |
| shapes by `tau` | `p^2; p^3, pq; p^4; p^5, p^2 q` | sympy check to `10^5`, no failure | yes |
| A'(i) | all sources of `2..14` have `tau in {3,4}`; first `t=4` value 15, first `t=5` value 16 | true; 15 <- 16; 16 <- 12 (first `t>=6`: 36 <- 24) | yes |
| A'(ii) | single source exactly `k = 2, 3` | exhaustive `k <= 9`; distinct Goldbach partition of `2^k - 2` for all `4 <= k <= 80`; single-source odd `v < 1000` = {3, 7, 9, 11, 39} (= S57) | yes |
| B'(i) | row-type sources of `P + 1` = `{2(P-2)}` iff `P-2` prime | algebra re-derived; brute force all primes `P < 5000`, no exception | yes |
| B'(ii) | `{6,25}, {28}, {496,652}, {8128,10682,10828}` | full brute force over all `m <= v^2` (to 66,064,384), no parity lemma: identical | yes |
| `s^-1(33550336)` | `{33550336, 44655764}` | sieve to 67,100,670 + odd squares `q < N` via spf: identical; `44655764 = 2^2·1669·6689` | yes |
| companions | `652 = 4·163`, `10828 = 4·2707`, `10682 = 2·7^2·109` | same | yes |
| `(N-7)/3` prime | true for 28, 496, 8128; false for next four | 7, 163, 2707 prime; 11183443, 2863289683, 45812897107, 768614336046650707 composite | yes |
| B'(iii) `<= 10^4` | 205 twin values, 115 exactly two, first third at 74 <- {70,142,5329} | 205; 115; 74 (max 6, mean 2.634) | yes |
| B'(iii) `(10^4, 10^5]` | 1019; 494 (0.4848); mean 2.816; max 11 | 1019; 494 (0.4848); 2.816; 11 | yes |
| `10^6` test (54617edb71) | 0.4427 in `(10^4,10^6]`; 0.4457 all; mean 2.943; max 13 | 7964/3526 = 0.4427; 8169/3641 = 0.4457; 2.943; 13 | yes |
| even untouchables `< 10^5`, `< 10^6` | 13,862; 150,231 | 13,862; 150,231 (+5: 13,863 = A057978(5), 150,232 = A057978(6); OEIS fetched) | yes |
| C'(1),(2),(4) | corollaries | own bookkeeping `n <= 200`, bound `10^9`: 0 violations; proofs checked | yes |
| later-doubling list `n <= 40` | 6, 10, 14, 18, 20, 22, 26, 28, 30, 34, 36, 38 | same; first 6 | yes |
| 9, 25, 27 row hits | 36; 100, 400; 216, 1728, 13824 | same (also 49 -> 196, 784; 81 -> 324, ..., 5308416) | yes |
| C'(5) | `D(1) = P(2)` only row-column | only `D(2) = P(2)`; `D(1) = P(2) ∪ {1}` | **no (off by 1)** |
| D'(ii) hits `k <= 64` | 6, 12, 18, 20, 21, 24, 30, 36, 40, 42, 48, 54, 60, 63 | same (sympy factorisation, repeated primes 3, 5, 7 only) = multiples of 6, 20, 21 | yes |
| orders | 6, 20, 21, 110, 156, 253 | `ord_(p^2)(2)` for 3, 5, 7, 11, 13, 23: same; only 3, 5, 7 have order `<= 64` (p < 500) | yes |
| Wieferich orders | 364, 1755 | 364, 1755 (= `ord_(p^2)`) | yes |
| abc `c <= 10^4` | 121; 94; 48; 30 | 121 (120 for `c < 10^4` = A216370); 94; 48; 30; one even term always | yes |
| abc `c <= 10^5` | 419; 323; 119; 63; best 1.5679 at 1 + 4374 = 4375 | 419 (418 strict = A216370); 323; 119; 63; 1.5679 at (1, 4374, 4375) | yes |
| `H(P_7)`, `H(P_11)` | 189; 95095 | own subset DP: 189; 95095 | yes |
| numerology arithmetic | `s(18)=s(51)=s(91)=21`; square-divisor part of 64 = 21; `189 = 3·63` | same | yes |
| hypotheses INDEX | 120 lines | 120 before and after 54617edb71 | yes |

## Attribution table

| note's claim | cited to | found in source? | same quantifiers/status? |
|---|---|---|---|
| `s(p^3) = Phi_3(p)`; `21 = Phi_3(4)` does not give `s(64) = 21`; `Phi_3(x^2) = Phi_3(x)Phi_6(x)` gives `21 = 7·3` | session69 aliquot sec. 4 | yes | yes |
| `{1,2,4} = <2>` Paley-7 difference set, doubling an automorphism, Frobenius on `F_8^x`, `H = 189`, `|Aut| = 21`, no proper module | session69 power sec. 7 | yes | yes |
| `{1,2,4}` has rank 2 in `F_8`, not a normal basis | aliquot124 result 5 | yes | yes |
| 2 and 5 untouchable PROVED; "only" conditional on distinct Goldbach; UNPROVED/CONDITIONAL; primes to `10^6` checked | session69 aliquot sec. 2-3 | yes ("UNPROVED / CONDITIONAL", 78,498 primes) | yes |
| S57's dropped antecedent repaired | aliquot124 CORRECTIONS | yes | yes |
| Pomerance–Yang; "Booker's `4·10^18` verification" | session69 sec. 3 | Booker cited for the distinctness footnote; the computation is Oliveira e Silva–Herzog–Pardi (S57) | **misattributed (Minor 6)** |
| THM-1370: `H` never 7 or 21, odd-spectrum completeness OPEN; Busch floor 3, 5, 9, 15, 25, 45, 75; "`21 = 3·7` recurses into the first gap"; `H = 63` attained, first `7·3^k` | THM-1370 | yes (75 is THM-1370's prediction, Busch's theorem) | yes |
| THM-4790: `tc` never 7 or 21 for `n <= 11` | THM-4790 | yes (file name differs: Minor 5) | yes |
| 7, 21 touchable; 5 an `H`-value (24 of 64 labelled 4-tournaments); `H != 7, 21` does not exclude `tc`; no predicate-preserving map; "similar small-number patterns remain an analogy" | session69 aliquot sec. 6 | yes, verbatim | yes |
| descent mod 7 fails (`2 ≡ 23 = s(57)`), `tc(P_7) = 9` | session69 CONNECTIONS | yes | yes |
| ray tableau: new doubling ray iff odd, new squaring ray iff non-square; "4 adds nothing, 6, 8, 10 add only their squares, 9 adds only its doubles, 25 caps the block" | lead_reflection2 sec. 3 | yes, verbatim | yes |
| doubling/squaring forests | arithmetic_seams sec. 2 | section exists | yes |
| `S_n`, `Q_n`, `P_r D_k = D_(rk) P_r`, `n = 2^a b^e`, `S_n ∩ S_m = S_L`, 36 has five addresses | session69 power sec. 2-5 | yes | yes |
| Theorems P1-P3 | aliquot124 chains/RESULT.md | yes | yes |
| earliest-seed theorem: new iff odd non-perfect-power; 15 new; 9 not a new seed of 36 | session69 power sec. 3, README | yes | yes |
| square-later predicate `a >= 1, u > 2^a`; six least; `n_min(a) = 2^(2a) + 2^a`; density 1/2; 4859 to `10^4`; even roots make 4 first | session69 power sec. 6 | yes | yes |
| 3 first with odd bases; 6 first even base for every `k >= 2`; `k_min = 2 + floor(log_u 2^a)` | aliquot124 README result 2, chains sec. 3 | yes | yes |
| `F(6) = {6,25}` first finite multiple fibre, `F(1)` infinite; twin identity; `2p = p + 3` iff `p = 3` | session69 aliquot (Immediate conclusions 1-3) | yes (section pointer imprecise: Minor 11) | yes |
| unique perfect prime-square target `(6, 5)`, both audited | RARE_SIX; aliquot124 result 1 | yes | yes |
| "another perfect target such as 28 is not excluded from having other source shapes" | "the audited uniqueness theorem" | in ROOT_ADDENDUM.md line 40 (the audit of RARE_SIX), with "by this theorem alone" | yes (name the addendum: Minor 10) |
| primitive triple `(n, s(n), sigma(n))/g`, perfect -> `1 + 1 = 2`; fixed-core bound `C <= (U+V) rad(ABC)` | session69 aliquot sec. 5 | yes (5.1, 5.2) | yes |
| square-later bijection with `(u - 2^a) + 2^a = u`; first excess at 72 (`1 + 8 = 9`) | session69 power sec. 6a | yes | yes |
| `6^2 = 4·9` multiplicative, not a triple | session69 sec. 5.1 / CONNECTIONS sec. 3 | yes | yes |
| `F_1` reading: `s^-1(0), s^-1(1)` agree with `tau - 1`; "2 with nothing" first place prime weights matter; `s(p^a) = |P^(a-1)(F_p)|` | S57 sec. 4 | yes | yes |
| "failure would need both shapes to fail at once" | S57 sec. 3 | yes, verbatim | yes |
| five single-source odd values below `10^8` | S57 sec. 3 | yes: 3, 7, 9, 11, 39 | yes, but A'(ii) is thereby S57's (Major M2) |
| row formula `s(2^a q) = (2^a - 1) q + 2^(a+1) - 1` = guarded formula with `C = 2^a` | aliquot_guarded_cycles_20261014e | the guarded identity `s(Cp) = (sigma(C) - C) p + sigma(C)` is there (and restated in perfect_core_twin_edges) | yes |
| `13,863 = A057978(5)`, `150,232 = A057978(6)` | OEIS | A057978 = "Number of untouchable numbers <= 10^n": 2, 5, 89, 1212, 13863, 150232, ... | yes |
| concurrent commits 2f483a9201, 365ab05300, 52fde1af0d | git | exist, 13:48 / 14:06 / 13:56, before the draft's author time 14:09:47 | yes |

## Prior-art findings (from the auditor's own knowledge, with what was verified online)

1. **`s(n) = 7` only for `n = 8`**: trivial folklore; it is `A048138(7) = 1`, and 7 is a "hermit number" (A057709).
   Fine as an elementary PROVED fact; it is not a discovery.
2. **"3 and 7 are the only Mersenne numbers with a unique aliquot source"**: not a named literature result, but a
   two-line consequence of the distinct-Goldbach mechanism, and already in the repository (S57 sec. 3; Major M2).
   Elementary; type CITED + corollary.
3. **`(1, 2^k - 1, 2^k)` is an abc hit iff `2^k - 1` is not squarefree**: textbook. The family
   `1 + (2^(6n) - 1) = 2^(6n)` with `9 | 2^(6n) - 1` and `rad(abc) < (2/3) c` is the standard example (Wikipedia "abc
   conjecture", verified by fetch; also in expository sources such as Granville–Tucker), and `p^2 | 2^(p(p-1)) - 1`
   gives it for every odd prime. Must be typed classical (Major M1).
4. **Almost perfect numbers give `1 + s(n) = n`**: this is the definition (`sigma(n) = 2n - 1`); "the powers of 2 are
   the only known ones" is classical (Guy, UPINT B2). Fine as stated, no novelty.
5. **28 is its own only source**: `28 ∈ A057709` (hermit numbers); elementary OEIS-level data. The conjecture "28 is
   the only perfect hermit number" (HYP-9415 (4)) is fine as OPEN; cite A057709.
6. **abc census counts**: 120 / 418 for `c < 10^4`, `c < 10^5` are OEIS A216370 (de Smit's tables); the note's
   121 / 419 are the `<=` counts (Minor 7). The best quality `1.5679` at `1 + 4374 = 4375` is a well-known high-quality
   triple.
7. **Untouchable counts** 89, 1212, 13863, 150232: A057978, as the note says.
8. A'(i) (the three-shape regime below 15) and B'(i) (row-source uniqueness): elementary, no prior statement found
   in the repository (greps for "proper divisors", "tau <= 3", "single source", "self-only", "squarefree",
   "almost perfect", "Wieferich", "ord_" over 05-knowledge, 01-canon, 07-reflections: only unrelated hits, plus
   THM-4644's `ord_(p^2)(2)` remark for Fermat/Wieferich, which is a different statement).

## Typing and pre-registration checks

- HYP-9415 reports the refuted range (power-of-two share 0.2840 vs `[0.30, 0.50]`) in its title, status and body;
  the note sec. 4 and both INDEX rows and the ledger say "one range refuted". Honest.
- Pre-registration: draft commit `e68e0ea0e7` (author 14:09:47; committer date rewritten by the rebase) contains
  HYP-9415 with the three ranges and no outcomes, the `.cpp`, and no C++ `.out`; `41713d3aa2` (author 14:18:37) adds
  the `.out` and the outcomes; the note's section 4 numbers equal the `.out` (0.4848, 2.816, 44655764, 0.7709, 0.2840,
  0.1504, 13,862). The `10^6` range `[0.42, 0.52]` is in HYP-9415 at `41713d3aa2`, its outcome 0.4427 in
  `54617edb71`; `_v1e6.out` matches my recomputation. Consistent with the stated story (git cannot prove the output was
  unread before 14:09, but nothing contradicts it).
- Theorem title/status vs note: consistent except the labels in Majors/Minors above ("tau <= 3", "every value
  v <= 14", B'(i) phrasing, `D(1) = P(2)`). THM-4805 "AUDITED: pending" matches the note's sec. 7.
- Section 5 typing table: the rows "Mersenne abc hits ... PROVED (D')" and "Mersenne single sources exactly 3, 7 |
  PROVED (A')" must become CLASSICAL and CITED (S57) respectively; all other rows are typed correctly.

## Commands run (auditor's scratch directory `.../scratchpad/s60_audit/`)

```
python3 recompute_a.py <session .out>      # sieve to 1e6: bound check, table v=0..60 vs session .out, untouchables <= 1000,
                                           # Lemma 1 minima and shapes, A'(i), Mersenne fibres k <= 9, Goldbach partitions k <= 80
                                           # (stopped after section A'(ii): my first B'(i) brute force was written too slowly)
python3 recompute_a2.py                    # B'(i) brute force P < 5000 (even sources n <= 2v-2), B'(ii) via parity lemma,
                                           # (N-7)/3 primality, C' bookkeeping n <= 200 (bound 1e9), D' factorisations k <= 64,
                                           # orders mod p^2, abc census c <= 1e4 (numpy), H(P_3), H(P_7), H(P_11) (own DP)
g++ -O2 -std=c++17 -o recompute_b.exe recompute_b.cpp && ./recompute_b.exe   # sigma sieve to 67,100,670 (7 s):
                                           # FULL brute force s^-1(6/28/496/8128) to v^2; s^-1(33550336) (+ spf sieve for odd squares);
                                           # twin counts <= 1e4 and (1e4, 1e5]; even untouchables <= 1e5; abc census c <= 1e5 (9 s total)
g++ -O2 -std=c++17 -o recompute_c.exe recompute_c.cpp && ./recompute_c.exe   # 1e6 twin test and untouchables <= 1e6
python3 linkcheck.py /tmp/math-wt-s60      # every relative link / YAML path in the note, THM-4805, HYP-9415, INDEX rows, ledger
git log/show on e68e0ea0e7, 41713d3aa2, 54617edb71 (author vs committer dates, draft HYP/note contents, diffs); git diff (working tree)
grep/ripgrep sweeps: "squarefree", "almost perfect", "self-only", "single source", "tau <= 3", "proper divisors", "ord_", "Wieferich",
                     "2^k - 1" near abc/rad, "not excluded from having other source shapes"
WebFetch: oeis.org/A057978, oeis.org/A048138, OEIS search 6,31,120,418,1268 (A216370), en.wikipedia.org/wiki/Abc_conjecture
```

Outputs: `recompute_a.out`, `recompute_a2.out`, `recompute_b.out`, `recompute_c.out`, `linkcheck.out` in the same directory.
No repository file was modified.
