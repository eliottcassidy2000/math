# Audit of opus-2026-10-14-S33 (number_tournament_bridges_20261014)

Independent adversarial referee: a general-purpose subagent, 2026-10-14, 98 tool uses. The harness did not let the referee write files, so the author saved this report verbatim from the referee's reply. The referee's scripts and outputs are in this directory. No other repository file was changed by the referee.

Audited versions (git blob hashes at audit time):

| file | hash |
|---|---|
| note | `41ef0d6276` |
| THM-4644 | `648c43686e` |
| HYP-9273 | `d43f638829` |
| results INDEX | `65a9dc869a` |
| hypotheses INDEX | `475010e625` |

The author kept editing during the audit; sections 5b and 5c, (B′) and `q = 31` were audited in their edited form.

**Verdict.**
- THM-4644 is mathematically correct; all of its numbers reproduce.
- All computed tables reproduce exactly, including `q = 31 → 41 (mod 256)` by a different algorithm (2.3 GB).
- 6 ERRORs: one novelty claim, one statement of Dickson's theorem, four labels or small numbers. 16 IMPRECISE items.

## ERRORS

- **E1. Novelty of THM-4644 part 1 overclaimed.**
  - Locations: note §1b ("the converse direction is new"), §0.3 ("New, exact"); THM-4644 "Relation"; hypotheses INDEX ("the converse is new").
  - Prior art: E. Akin, *Rock, Paper, Scissors, Etc — Topics in the Theory of Regular Tournaments*, arXiv:1806.11241 (2018; v4 2020).
    - Theorem 3.13: for odd `2n+1`, the multipliers act freely on the game subsets (connection sets) iff `2n+1` is a squarefree product of Fermat primes.
    - Theorem 3.14: every nontrivial odd-order multiplier subgroup fixes some connection set, and all invariant ones are counted (`2^m`).
    - After 3.14 Akin writes "We do not know the answer for 17 or the other known Fermat primes", i.e. whether `Aut = Z_p` for every rotational tournament.
  - What is actually added: the passage from multipliers to the full Aut by Burnside/Alspach (as in THM-871), which settles Akin's question for every Fermat prime.
  - Fix: cite Akin and restate the novelty accordingly.
- **E2. The census counts an odd prime power, not an odd prime.**
  - Location: note §0.2 "Data" and the §3b table header ("an odd prime divides every `σ(n_i)` but not every member").
  - The script counts `odd(Q) ∤ G`: 229/231, 9/10, 4/4, 2/2, 0/1, 0/1, 0/1, all reproduced. Under the stated criterion (an odd prime dividing `Q` but not `G`), the amicable count is **221**.
  - In 8 pairs the excess comes only from a higher power of a prime that already divides `G`:
    - (69615, 87633): `G = 819`, `Q = 157248`;
    - (947835, 1125765): `G = 135`, `Q = 2073600`;
    - (8619765, 9627915), (11498355, 12024045), (13671735, 15877065), (17257695, 17578785), (22508145, 23111055), (46521405, 53011395).
  - Fix: say "odd prime power (`odd(Q) ∤ G`)" and optionally add the prime-level column 221, 9, 4, 2, 0, 0, 0.
- **E3. σ-divisibility percentages mislabelled.**
  - Location: note §3b, bullet 3 of "What the table shows", and the script label.
  - 98.5% (`3 | σ`) and 93.6% (`9 | σ`) are over all 535 members of non-trivial cycles (lengths 2, 4, 5, 28). For the 502 members of 2- and 4-cycles they are **99.6%** and **97.2%**. The baselines 85.1% and 57.4% are correct.
- **E4. Factor in "The mechanism" wrong.**
  - Location: note §1, "grows by a factor `p` in the exponent".
  - The unit-orbit count grows by `p`; the exponent `E` grows by `(p+1)/2` (`12 → 6564 = 547·12`; `2 → 3512 = 1756·2`).
- **E5. Dickson misstated.**
  - Location: note §0.2, "a sociable cycle of odd length `k > 1` cannot be regular, i.e. the gcd … cannot be a unitary divisor of each member".
  - KPP Theorem E assumes `a = gcd > 1`. With gcd 1 the gcd is trivially unitary, and no theorem excludes coprime odd cycles; KPP's own introduction uses the same loose sentence. §3b states it correctly.
- **E6. Pattern count wrong.**
  - Location: note §3, "AADD (8 up to rotation)".
  - The correct count is **7**: 7 + 1 + 2 = 10 four-cycles.

## IMPRECISE

1. **(I1)** §3b Lemma: Pollack's "`Q | G` unless `2 ∥ Q`" needs `k > 1` (an even perfect number has `Q = 2G` with `4 | Q`).
2. **(I2)** The codex remark on `I + P` treating lengths 3 and 5 alike is in `aliquot_dyadic_frontier_20261014e.md` §6, not in the guarded-core packet linked.
3. **(I3)** §1 conflates THM-4631's `<2>`-orbit count (1093: `3 → 3279`) with THM-4644's `<g>`-orbit count (`12 → 13116`). They jump by the same factor `p` and coincide only at 3511.
4. **(I4)** "doubling-invariant tournaments on `Z_(p^2)`" (Contents, hypotheses INDEX, commit message) is false at 1093: `−1 = 2^182 ∈ <2>` mod `p^2`, so none exist. Say "invariant under the odd part of doubling".
5. **(I5)** THM-4644 proof part 1: Burnside gives `Aut ≤ AGL(1,p)` only up to permutation isomorphism. Missing step: `Aut` has a normal Sylow `p`-subgroup of order `p`, which must be the translation group, so `Aut ≤ N(Z_p) = AGL(1,p)` in the given coordinates. The front matter also says THM-871 is "via Burnside"; THM-871 used Alspach.
6. **(I6)** Stale sociable data. Moews' list (djm.cc/sociable.txt, header 11-IX-2026): 5446 cycles of length > 2, of which **5433** of length 4 and **5** of length 8. The exhaustive search (member preceding the largest) reaches **10^14** (Bodyagin, Oct 2016), not 4.27·10^11; for odd smallest members, 10^15 and 10^16.
7. **(I7)** The `x/(log x)^1.5 ≈ 10^7` comparison drops the `o(1)` in the Iannucci–Luca exponent. It is indicative only; the count is dominated by the primes `2b^2 − 1`, about `√x/log x`.
8. **(I8)** HYP-9273 (A) lists the double-flip census as evidence, but double flips test only `k = 2`. For `k >= 4`, the two members in `□ ∪ 2□` of a mixed cycle need not be consecutive.
9. **(I9)** Heuristic stated as mechanism.
   - "must avoid every common odd divisor" should read "every common odd divisor of the `σ(n_i)` must divide every member".
   - "This is why / explains why odd lengths are rare": KPP and Pollack prove only upper bounds for odd order.
10. **(I10)** §3 corollary 5 yields only an abundant member. The deficient member follows from corollary 4: the minimum member is abundant and the maximum deficient.
11. **(I11)** References:
    - Pollack, "A remark on sociable numbers of odd order", *J. Number Theory* 130 (2010) 1732–1736;
    - Cohen, *Math. Comp.* 24 (1970) 423–429 is missing (source of "regular");
    - the Moews papers (*Math. Comp.* 57 (1991), 61 (1993)) and the list URL are missing;
    - Hagis pages 963–968;
    - "Moon (1966)" in the status block against the 1968 book in the references;
    - Bratley–Lunnon–McKay `5·10^7` is attested only via Wikipedia (AMS returned 403).
12. **(I12)** §0.6 calls the "rigid radius-2 hinge … a small certificate". Mac-mini has SOLVER verdicts only (no DRAT) on roughly 1,000-vertex graphs, and its note says the other cross-edges are essential. Drop "small certificate".
13. **(I13)** Sylvester's 2, 3, 7 and 43 are all Heegner numbers, not only 7 and 43.
14. **(I14)** Latent code issues, not triggered:
    - `find_cycles` in both aliquot scripts caps walks at 1000 steps and marks truncated paths done; an uncapped census still gives 248;
    - the palindromic script header says 1 GB, but it needs 3.2 GB at `q = 47`.
15. **(I15)** §5b: "(B) ⇔ `F ≡ (1−P)/2 (mod 4)`" needs only `P` odd. (B′) is needed only for "`F` odd iff `q ≡ 7 (mod 8)`".
16. **(I16)** The INDEX lines inherit E1 and I4 (hypotheses INDEX) and I9 (results INDEX).

## VERIFIED OK, with method

### A. THM-4644
- Proofs of parts 1–3 checked line by line (gap I5 noted).
- `E` for every odd prime `p < 4000` by explicit orbit decomposition of `Z_(p^2) \ {0}`: jumps exactly at 1093 (6564 vs 12) and 3511 (3512 vs 2). All six table rows match.
- Brute force over all connection sets:
  - `Z_9`: 4; `Z_25`: 16; `Z_49`: 4 (base 2);
  - the Wieferich-regime analogue `Z_49` with `g = 18` (`18^3 ≡ 1 mod 49`): 256, matching the Wieferich formula.
- Backtracking `|Aut|` over all circulant tournaments for `p = 3…17`: rigid exactly at 3, 5, 17. Stabilizers trivial exactly at 3, 5, 17 for `p <= 31`.
- Orders: `ord_3511(2) = 1755`, `ord_1093(2) = 364`, `ord_1093(3) = 7`.
- Rotational classes: 1, 1, 16, `2^120`, `2^32752`.
- THM-4631 numbers 3 → 3279 and 2 → 7022 confirmed.

### B. Paley
- `H(P_q)` by a full DP over all start vertices: 3, 189, 95095, 1172695746915, 15760206976379349. Backtracking `|Aut| = q(q−1)/2`; orbit counts all `≡ 1 (mod 8)`.
- `P_27`, rebuilt on `F_3[y]/(y^3 + 2y + 1)`: `|Aut| = 1053`; `H ≡ 95 (mod 256)`; `H/1053 ≡ 3` and `H/351 ≡ 1 (mod 8)`. In fact `H/351 ≡ 1 (mod 256)`.
- `Z_19` circulant residues `{1:128, 3:126, 5:96, 7:162}` with `|Aut| = p|M|` by backtracking for all 512; `Z_13` and `Z_17` mixed; `Z_5` all 3.
- DR15: `H = 198335025`, `|Aut| = 21`, orbit count `9444525 ≡ 5 (mod 8)`.
- `q = 31`: 41 (mod 256) by an independent layered colex-rank DP (2.3 GB, 40 s), which reproduces 9, 193, 41, 185 at `q = 7, 11, 19, 23` and 1 at `q = 27`. The pre-registration is genuine: commit `e399906e95` names `q = 31` as the first test.
- The §5b split is proved correctly (the only involutions in the coset are the `ρ_b`).
- `pal` = 1, 9, 185, 573057, 63871533, 2108041228485; brute force confirms 9 and 185. Residues 143157 (mod 172032) and 36777 (mod 188416), all `≡ 1 (mod 4)`. `P mod 8` and `F` match, and `F ≡ 31 (mod 128)` at 31.

### C. Aliquot
- Two independent σ-sieves agree on `1..10^8` (a linear multiplicative sieve and a segmented divisor-pair sieve to 388,251,360).
- 248 cycles: 5, 231, 10, 1, 1. No 3-cycle; all parity-homogeneous; none meets `□ ∪ 2□`. 54 all-odd: 53 pairs plus the 28158165 four-cycle.
- 236 amicable pairs with smaller member `< 10^8`, matching OEIS A066873; 5 of these have the larger member above `10^8`, leaving 231.
- Cycle identities hold; every non-trivial cycle has an abundant member and a deficient member.
- Double flips to `10^9` by a segmented trial-division sieve:
  - 4476 in total, 4458 of them prime, matching an independent count of primes `2b^2 − 1`;
  - 18 composite; all decade counts and all 18 factorizations match;
  - 11 of 18 have a prime factor with `p + 1 = 2b^2`;
  - maximum ratio 0.8299 at 949818597;
  - none closes a pair.
- Overflow bound: product bound 3.0565, observed maximum 3.0058.
- `p = 71` conic: `t = 205, 212, 213` give `q = 409, 3391, 85625 = 5^4·137`, and a direct prime scan finds exactly 409 and 3391.

### D. Frustration
- The lemma holds on all 255 cycles. Pollack's Lemma 5 and KPP Theorem E read in the papers.
- `det(I+P)`, the odd-`k` inverse and the even-`k` kernels checked symbolically for `k = 1..9`.
- Census:
  - regular counts 189 / 8 / 3 / 2 / 0 / 0 / 0;
  - excess counts 229 / 9 / 4 / 2 / 0 / 0 / 0;
  - `G, Q` = 8, 8 (5-cycle) and 2, 2 (9-cycle);
  - median excess 13041, and 3 divides it in 93.9% of pairs.
- The seven published cycles close, with lengths 6, 6, 6, 6, 8, 8, 9.

### E. Citations and facts
- Gioia–Vaidya, D'Abramo, Hagis (`> 10^60`), Iannucci–Luca (`x/(log x)^(3/2+o(1))`, and Pollack's twice-square remark), Pollack 2011, KPP and Dickson: confirmed against the papers or AMS pages.
- Search records: Chernykh complete below `10^20` (4197267 pairs, A066873); PrimeGrid Wieferich search to `2^64` (Dec 2022); the Wikipedia table matches the note but is stale (I6).
- Ramsey: `R(5,5) <= 46` (Angeltveit–McKay, arXiv:2409.15709); Wikipedia attributes the uniqueness of `P_17` to Evans–Pulham–Sheehan, matching the note.
- Recomputed:
  - Paley clique numbers 2, 3, 3, 4, 4, 5, 7 (own branch and bound plus networkx);
  - pseudoperfect numbers 2, 6, 42, 1806, 47058;
  - `2^32 − 1` as the product of the five Fermat primes, and the collector denominators;
  - the `(2,3,k)` Euler characteristics;
  - `1807 = 13·139`, `139 = 3^7 − 2^11`;
  - no `τ = 4` number is abundant;
  - the mod-21 residues.
- Canon THM-200, THM-029, THM-1370, THM-307, THM-871 and THM-582 are consistent with the note.

### F. Typing
- Nothing merely computed is typed PROVED.
- The answers in brief are honest on the main points: the problem is OPEN, the `10^18` search is inside Chernykh's range, and Dickson does not single out 3.
- Codex is credited correctly except the I2 link.

**Note for the author.** `|Aff|` is odd, so (B) is equivalent to `H(P_q) ≡ q(q−1)/2 (mod 8)`. Only `H mod 8` is needed, and the layered DP brings `q = 31` down to about 2.3 GB.

## Files in this directory

`audit_A_thm4644.py`, `audit_B_paley.py`, `audit_B2_q31_palindromic.py`, `audit_C_aliquot_cycles.py` (+ `audit_C_cycles.json`), `audit_C_double_flips.py` and `audit_E_misc.py`, each with its `.out`.

Sources the referee used:
- Pollack, socodd;
- KPP;
- Pollack, *Powerful amicable numbers*;
- Akin, arXiv:1806.11241;
- Moews, djm.cc/sociable.txt;
- Wikipedia: Sociable number, Paley graph;
- OEIS A066873, A001220;
- Chernykh, sech.me/ap;
- D'Abramo, arXiv:math/0501402;
- Hagis (AMS);
- Angeltveit–McKay, arXiv:2409.15709.

## Disposition (author, same session)

All 6 ERRORs and 16 IMPRECISE items were applied: note, THM-4644, HYP-9273, both INDEX lines, and scripts (the cycle-walk cap removed; prime-level column and population labels added). Recorded as MISTAKE-608.
