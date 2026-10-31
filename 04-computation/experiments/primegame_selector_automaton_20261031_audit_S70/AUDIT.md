# Independent audit: THM-4852 / HYP-9448 (opus-2026-10-31-S70): PRIMEGAME phase law, selector automaton, n = 5 metagraph

Auditor: independent referee agent (machine `opus`), 2026-10-31. Worktree `math-wt-s70-fractran`, commit `0bfb7da871`
(pre-registration commit `adfb264520`). All code in this folder was written from scratch; nothing was imported from the
session's scripts. No file outside this folder was edited; nothing was committed.

**VERDICT: SOUND WITH CORRECTIONS**

Every number in the deliverables reproduces with independent code. The phase law and
the agreement of the two lists hold for both lists to n = 2500 with a stronger check, and the step counts match OEIS. The
digit formulas are proved and stress-tested. But HYP-9448 is the published proof of Conway's Theorem 1, so it is not open:
THM-4852 (1) holds unconditionally for all n. The owner's list is Conway's own 1987 PRIMEGAME, not a variant. The selector
automaton is prior repository work and prior literature. The composite map is not a cellular automaton in the strict
sense. The metagraph negative test is vacuous.

Files audited:
- `05-knowledge/results/primegame_variant_selector_automaton_metagraph_20261031.md` (below "note")
- `01-canon/theorems/THM-4852-primegame-phase-law-...-base-29-primorial.md` (below "THM")
- `05-knowledge/hypotheses/HYP-9448-primegame-phase-law-...-minus-one.md` (below "HYP")
- `04-computation/experiments/primegame_variant_phase_law_20261031.c` (+ `.out`) and `primegame_selector_automaton_metagraph_20261031.py` (+ `.out`)
- the two INDEX lines added in `0bfb7da871` (`05-knowledge/hypotheses/INDEX.md` line 2, `05-knowledge/results/INDEX.md` line 3)

## 1. Findings

### MAJOR-1: HYP-9448 is a published theorem (Conway 1987, section 5, Figure 2); THM-4852 (1) is unconditional for all n

**Where.**
- HYP: line 4 (`status: OPEN for all n`), lines 27-28 (the file treats its all-n form as a hypothesis to be checked),
  lines 33-41 (pre-registered predictions).
- THM: line 3 (title, `verified n <= 2500; HYP-9448`), line 4 (`the all-n equivalence conditional on HYP-9448`),
  lines 29-33, novelty lines 63-66 ("New to the repository: the phase-law explanation ...; the role of 1/7").
- Note: lines 6-7, 11, 46-47, and the role list at lines 50-54.
- `05-knowledge/results/INDEX.md` line 3 (`all-n list equivalence conditional on HYP-9448`) and
  `05-knowledge/hypotheses/INDEX.md` line 2 (status list `... / OPEN`).

**Evidence.**
- Conway, "FRACTRAN: a simple universal programming language for arithmetic" (1987), is the paper all three files cite.
  The relevant pages were read here from a scan of the Springer chapter: pp. 4-6, 9-10 and 25-26.
  - Theorem 1 (p. 4) states PRIMEGAME with exactly the owner's fourteen fractions, ending 15/2, 1/7, 55/1.
  - Section 5, "Your PRIMEGAME Guarantee!" (pp. 9-10), labels them A-N, with L = 15/2, M = 1/7, N = 55/1, and proves
    Theorem 1 by Figure 2.
- **What Figure 2 shows.** Take 0 < d < n and the configuration 5^n 7^d 13.
  - If d does not divide n, the orbit returns to 5^n 7^(d-1) 13.
  - If d divides n, rule I produces the state-free number 2^n 7^(d-1). Then the word L^n M^(d-1) N gives 3^n 5^(n+1) 11,
    and (EF)^n K gives 5^(n+1) 7^n 13.
  - The text concludes that the game tests n-1, n-2, ..., 1 until it first meets a divisor of n.
- That is HYP-9448's law. D(n) is the first divisor met, and the owner's reset route quoted at note line 44 is Conway's
  `L^n M^(d-1) N` verbatim.
- **Secondary sources.**
  - Wikipedia "FRACTRAN" states the same law for input 2^n 7^m with 0 <= m < n: the output is 2^(n+1) 7^(k-1), k the
    largest proper divisor of n+1. It refers to Havil (2007) for a step-by-step account.
  - OEIS A007546 carries A. P. Heinz's 2009 closed form for the steps between consecutive phases:
    f(n) = n - 1 + (6n+2)(n - b) + 2 sum_{d=b}^{n-1} floor(n/d), with b the largest proper divisor. This formula presupposes
    the law. The audit verified it exactly at every phase k <= 2500 for Guy's list, and its shift by D(k-1) - 1 for the
    owner's list (section 2).
- **Both lists were already known to emit exactly the primes.** Conway 1987 Theorem 1 covers the owner's list. Guy 1983
  presents the 15/14 list as Conway's prime-producing machine; its text was not accessed here, but Figure 2's A-K
  analysis covers it verbatim (section 4). The repository already typed "Conway's theorem" CITED in
  `procgen_atlas_20260924` L1. So "the two programs emit the same primes forever" was never open.

**Correction (exact).**
- **HYP-9448 status:** `CLOSED - KNOWN (Conway 1987, section 5, Figure 2 = the proof of Theorem 1, for Conway's list
  = the owner's; Figure 2's A-K analysis applies verbatim to Guy's 1983 list, whose reset (15/14)^(D-1) (15/2)^(n-D+1)
  (55/1) is valid since D(n) - 1 < n). The pre-registered run is a FINITE-EXACT replication, not evidence for an open
  statement.`
- **HYP lines 27-28:** replace with "The law is Conway's own proof of Theorem 1 (1987, section 5, Figure 2); Wikipedia
  states it; OEIS A007546's step-count formula presupposes it."
- **THM title:** `(verified n <= 2500; HYP-9448)` becomes `(Conway 1987, Figure 2; replicated n <= 2500)`.
- **THM status:** `PROVED for all n (elementary corollary of the CITED phase law, Conway 1987 Figure 2) + FINITE-EXACT
  + CITED`. Delete "conditional on HYP-9448" in the THM, the note (lines 6-7, 11, 46-47) and results INDEX line 3, and
  replace OPEN with KNOWN in hypotheses INDEX line 2.
- **THM novelty:** "CITED, not new: the phase law, the trial-division semantics, and the roles of L = 15/2, M = 1/7,
  N = 55/1 and I = 1/17 (Conway 1987, Figure 2). New to the repository: the explicit two-line corollary that the two
  published lists pass through identical configurations at every phase boundary."
- **MISTAKE ledger:** add an entry. The cited paper's own proof section was the claimed hypothesis, the same pattern as
  MISTAKE-690/692. A complete all-n proof for both lists is in section 4 below.

### MAJOR-2: the selector-automaton formulation is not new; repository and literature credit are missing

**Where.** THM lines 45-52 and novelty lines 63-65 ("New to the repository: ... the selector automaton formulation");
note lines 63-69 and 125-135.

**Evidence.**
- **Repository prior work.** `05-knowledge/results/collatz_oscillation_gilbreath_20260926.md` section 3 (line 119;
  PROVED, brute-force checked, audited) already states the same structure for Collatz. In base 6, T is a radius-1
  cellular automaton on the digit string with one global bit (the parity of the lowest digit) selecting the branch.
  - It credits Cloney-Goles-Vichniac 1987, Korec 1992 and Kari 2012.
  - `MISTAKES.md` line 1497 records the earlier attribution slip for it.
  - The note cites instead `collatz_one_out_edge_systems_20260930.md` section 4, which is only the base-2 non-locality
    remark.
- **Literature.**
  - Kari-Kopra, "Cellular automata and powers of p/q", RAIRO-ITA 51(4) (2017) 191-204, arXiv:1710.05737: a CA that
    multiplies by p/q in base pq. This is the per-fraction digit locality.
  - Conway 1987, section 11 (p. 26): a FRACTRAN game is a Collatz function g with g(N)/N periodic, i.e. the branch
    depends only on N mod D. Conway also remarks there on Collatz functions whose segments are selected by the last
    decimal digit.
  - Cloney-Goles-Vichniac, Complex Systems 1 (1987) 349-360: "quasi cellular automaton" for 3x+1.
- THM-4852 (3) is therefore the routine combination of these two ingredients.
  - Each reduced fraction p/q with p, q | B acts as a radius-1 CA.
  - The branch is a function of x mod lcm(q_k), hence of x_0 when lcm(q_k) | B.
- No paper presenting general FRACTRAN programs this way was found, so the general statement is at most new to the
  repository as an extension.

**Correction (exact).** THM novelty and note section 3 read: "New to the repository: the extension of the base-6 Collatz
selector automaton (`collatz_oscillation_gilbreath_20260926.md` section 3) to every FRACTRAN program in lowest terms whose
numerators and denominators divide B (B = P# for primes <= P). CITED:
- the per-fraction locality: Kari-Kopra 2017 (the p/q automaton in base pq), Kari 2012 (3/2 in base 6);
- the branch-by-residue structure: Conway 1987, section 11;
- the quasi-CA notion: Cloney-Goles-Vichniac 1987."

Add the gilbreath note to THM `related`.

### MINOR-1: provenance and naming of the two lists are inverted in places

**Where.**
- HYP title (line 3) and statement (lines 14-18): the 15/14 list is called "Conway's PRIMEGAME", the owner's list "the
  owner's variant with 15/14 split into 15/2, 1/7".
- Simulator header (`.c` line 1; "standard" = Guy's list, "variant" = the owner's).
- Note lines 53-54: "1/7 ... exists only because the owner's list splits Guy's 15/14 into 15/2 and 1/7".
- Note lines 25-26 ("Conway's PRIMEGAME in its Wikipedia/OEIS form") and THM lines 18-20 ("the Wikipedia/OEIS form")
  never say that this is the list Conway himself published in 1987.

**Evidence.**
- Conway 1987 Theorem 1 (p. 4) and section 5 (p. 9) print the owner's list. Rule M = 1/7 is part of Conway's published
  PRIMEGAME.
- Guy's 1983 Math. Mag. list ends 15/14, 15/2, 55/1. Sources: the OEIS A202138 comment; Wikipedia's variant footnote
  citing Guy 1983 p. 26 and Conway-Guy 1996 p. 147; OEIS A007546 "fast version" = the 15/14 list.
- OEIS A007547 "slow version" = the 1/7 list. The audit's step counts reproduce A007546 from the 15/14 list and A007547
  from the owner's list (34/34 terms each). The labels "Guy 1983" and "Wikipedia/OEIS" used by the session and by
  `procgen_atlas_20260924_primegame.py` are therefore correct as far as they go.
- But the owner's list is Conway's 1987 original, not a split of Guy's. Historically the order is reversed: Guy 1983
  published 15/14 first, and Conway 1987 printed 15/2, 1/7.

**Correction.**
- **HYP title:** `PRIMEGAME phase law: Conway's 1987 PRIMEGAME (...15/2, 1/7, 55/1; the owner's list) and Guy's 1983
  list (...15/14, 15/2, 55/1) pass, at each entry into the state-free set, through 2^n 7^(D(n)-1) ...`.
- **HYP lines 14-18:** name the lists the same way.
- **THM lines 18-20 and note lines 25-26:** "the owner's list is Conway's PRIMEGAME as printed in Conway 1987, Theorem 1
  (also Wikipedia, OEIS A007542/A202138/A203363); Guy's 1983 list (also Conway-Guy 1996, p. 147; OEIS A007546 'fast
  version') ends 15/14, 15/2, 55/1".
- **Note lines 53-54:** "1/7 (Conway's rule M) erases the counter 7^(D(n)-1); Guy's earlier list instead consumes it
  with 15/14".
- **HYP naming paragraph (lines 51-54):** add "the simulator's 'standard' = Guy 1983, 'variant' = Conway 1987".

### MINOR-2: the phase law is mis-stated as a statement about every state-free configuration

**Where.**
- Note line 39: "Whenever no state prime is present, the configuration is 2^n * 7^(D(n)-1)".
- THM line 29: "the k-th state-free configuration is 2^k 7^(D(k)-1)".
- HYP lines 15-16: "Every time the configuration contains none of the state primes ..., it equals 2^n 7^(D(n)-1)".
- Results INDEX line 3: "Every state-free PRIMEGAME configuration is 2^n 7^(D(n)-1)".

**Evidence.** As written this is false. Each state-free stretch contains further state-free configurations. For example,
the owner's list sends 112 = 2^4 * 7 to 840 = 2^3 * 3 * 5 * 7, and Guy's sends 112 to 120 = 2^3 * 3 * 5. Neither is of
the form 2^m 7^j. All the code (session and audit) tests the first state-free configuration after a stateful one.

**Correction.** "Each time the orbit enters the state-free set (the first state-free configuration after a stateful
one; the start 2 is entry 1), the configuration is 2^n 7^(D(n)-1) for n = 1, 2, 3, ... in turn". In the THM and the
INDEX, "the k-th state-free configuration" becomes "the k-th entry into the state-free set".

### MINOR-3: a step of the equivalence proof is missing

**Where.** THM line 22 / note lines 84-94. The text proves only one direction: without a state prime, only the last
three rules act. It then asserts "From then on the lists are identical until the next state-free configuration".

**Evidence.** The assertion also needs the converse. In every stateful configuration, the first applicable rule must be
among A-K, the rules common to both lists. This is true: each state prime 29, 23, 19, 17, 13, 11 is by itself the
denominator of F, G, H, I, J, K respectively, and all of these precede the last three. The audit checked it at every
stateful step to n = 2500 (check (c)).

**Correction.** Insert after "only the last three rules can act": "Conversely, in a configuration with a state prime the
first applicable rule is among A-K, because 29, 23, 19, 17, 13 and 11 are the whole denominators of F, G, H, I, J, K;
so the two lists act identically between phase boundaries, and neither ever halts."

### MINOR-4: "radius-1 cellular automaton" over-states (3); the composite is not a CA

**Where.** THM title (line 3: "PRIMEGAME is a radius-1 cellular automaton on base-29# digits ..."), THM lines 45-47,
note lines 63-65.

**Evidence (`audit_selector_ca.out`, T5).**
- **Not shift-commuting.** F(29# * 2) = (17/91) * 2 * 29# = 2417248020, but 29# * F(2) = 15 * 29# = 97045398450.
- **Unbounded dependence on the units digit.** Take x = c * 29#^6 + 2 and x' = c * 29#^6 + 1 with c = 123456789. They
  differ only in digit 0, yet digit 6 of the image is 925925917 against 320430165, because the units digit selects 15/2
  against 55/1. So digit i of F(x) depends on x_0 at distance i.
- The mathematics of (3) is correct: the digit formulas and the selector-by-x_0 hold (tests T1-T4 in
  `audit_selector_ca.out`). What is wrong is the type. Each fraction is a genuine radius-1 CA, and the program is a finite family of them with a global selector.

**Correction.**
- **Title:** "PRIMEGAME is a selector (quasi-)cellular automaton on base-29# digits: fourteen radius-1 local rules, one
  chosen globally by the units digit".
- **THM (3) first sentence:** "... is computed digit-locally on base-B digit strings: each fraction acts as a radius-1
  cellular automaton, and the fraction to apply is a function of the units digit alone. The composite map is not itself
  a cellular automaton (it does not commute with x -> Bx); it is a quasi cellular automaton in the sense of
  Cloney-Goles-Vichniac."

### MINOR-5: THM-4852 (4), "no adjacency-preserving bijection", is vacuous and depends on the chosen graph

**Where.** THM lines 54-61 (stated as a theorem clause); note lines 70-72 and 146-153 ("in either reading").

**Evidence (`audit_anatomy_flow_metagraph.out`).**
- **The session's own numbers are right.**
  - Flow graph of the owner's list: 24 arcs, 19 pairs, no loops; simple degrees 3,3,3,3,3,4,4,5,5,5.
  - Metagraph, parsed directly from THM-4729's n = 5 table: 8 blue + 24 black lines, 18 simple pairs, black loops at
    nodes 5 and 8, pure-blue nodes 0, 3, 9; simple degrees 1,1,1,3,3,5,5,5,6,6.
- **The verdict is decided by the edge count alone (19 against 18).** In the weaker, homomorphic sense it is decided by
  the degrees. The metagraph has two degree-6 nodes and the flow graph has maximum degree 5.
- **Non-isomorphism is the generic outcome.** 0 of 2000 random G(10,18) graphs are isomorphic to the metagraph's simple
  graph.
- **The conclusion depends on the graph.** Take the prime co-occurrence graph (p ~ q when both divide a*b for one
  fraction a/b; 28 edges). There are 70 bijections nodes -> primes that carry every metagraph pair to an adjacent prime
  pair. None of them sends {TT5, 1->C3->1, T5} onto {7, 17, 5 or 11}, the owner's 7, 17, 55.
- **Guy's list has a different flow graph:** 26 arcs, 21 pairs. The THM's numbers are for the owner's list only.
- **A graph on the fractions.** The owner's {7, 17, 55} are fraction parts, so a graph on the fractions is the more
  natural carrier. The rule-succession graph of the orbit (Conway's flowchart; the arc set is the same after 40 and
  after 120 phases) has 14 nodes, 23 arcs and 18 simple pairs, with self-loops at L = 15/2 and M = 1/7. That is the same
  "18 pairs + loops at two nodes" as the metagraph, but on 14 nodes, which shows how cheaply such count matches arise.
  - Merge Conway's loop pairs AB, DG, EF and the inverse pair JK. This gives 10 macro-rules, of which exactly three
    (I, M, N) are pure. It has 13 pairs and is neither isomorphic to the metagraph nor embeddable in either direction.
- The owner's triple {7, 17, 55} consists of fraction parts, not registers, so "primes <-> nodes" is not even the
  owner's correspondence. 55 = 5 * 11 is of shape pq, with defect -2, so the triple also fails the {p, p^3, p^2qr}
  reading shape by shape. Only 7 and 17 are balanced.

**Correction.** Demote THM (4) to a remark, typed FINITE-EXACT and uninformative:
"PRIMEGAME's register-flow graph is not isomorphic to the merged metagraph's simple graph (19 against 18 edges; degree
6 against at most 5). This is the generic outcome, and it depends on the chosen graph: for the prime co-occurrence graph
there are 70 adjacency-preserving bijections, none matching the owner's {7, 17, 55} to the pure-blue nodes. No invariant
is transported in any reading tried."

In note line 153, "in either reading" becomes "for this flow graph, in either the simple or the multigraph reading".

### MINOR-6: "lcm = 29#, a Conway-Collatz map of modulus 29#" is uncredited

**Where.** THM line 36; note lines 59-60.

**Evidence.**
- Conway 1987, section 11, p. 26: g(N)/N is periodic, so a FRACTRAN game is a Collatz function with b_N = 0.
- The repository's L1 printed this exact lcm on 2026-09-24 (`05-knowledge/results/procgen_atlas_20260924.out` line 716:
  "the fraction used depends only on n mod lcm(denominators) = 6469693230").

**Correction.** Append "(CITED: Conway 1987 section 11; computed in procgen_atlas_20260924 L1)".

### MINOR-7: wrong source for the F = U + S convention

**Where.** THM line 38: "Under the repository's F = U + S (THM-4733's divisor counts)".

**Evidence.**
- THM-4733 defines no F, U, S. Its status explicitly does not claim "any arithmetic consequence for F = U + S".
- The convention F = tau - 2, U = omega - [prime], S = 2^omega - 1 - [squarefree] is written out in
  `05-knowledge/results/aliquot_66_five_fifteen_20261028/PROOFS.md` line 50. The session's formulas agree with it, and
  brute force to 2*10^5 confirms that the balanced shapes are exactly p, p^3, p^2qr.

**Correction.** Change the citation to "(the repository convention F = tau-2, U = omega-[prime], S =
2^omega-1-[squarefree]; aliquot_66_five_fifteen_20261028/PROOFS.md)".

### MINOR-8: the session's C checker does not test "for the next n ... in turn"

**Where.** `primegame_variant_phase_law_20261031.c` lines 38-44. The code reads n from r2 and tests only r3 = r5 = 0 and
r7 = D(r2) - 1. It never compares r2 with the phase index. HYP P1 ("for the next n = 1, 2, 3, ... in turn") was
therefore checked only indirectly, through the count `phases=2500/2500` at stop.

**Evidence.** The Python driver does check r2 = index, but only to n = 300. The audit simulator checks r2 = k exactly at
every entry, for both lists, to k = 2500. A negative control proves the check is live: `audit_primegame_sim 200 control`
falsifies the expected law at k = 97 only and fails exactly there, at step 1274952 (owner) / 1273490 (Guy). These are
A007547(25) and A007546(25). There is no effect on the result.

**Correction.** In the HYP results, add "(the session's simulator tests r3 = r5 = 0, r7 = D(r2) - 1 and the phase count;
the audit simulator additionally tests r2 = k at every entry, to k = 2500)".

### MINOR-9: pre-audit status claims and a leftover reservation

**Where.** THM line 4 and note line 13 say "independently audited (see the note's section 8)", while section 8 (note
line 226) says "(pending)". `01-canon/theorems/THM-4852-RESERVED-opus-s70.md` is still in commit `0bfb7da871`; it is
deleted only in the uncommitted working tree.

**Correction.**
- Replace the audit claims with a pointer to this file. Fill note section 8 with the verdict and the corrections applied.
- Commit the deletion of the reservation stub.

### NOTES

1. **Hypotheses of (3).** "Squarefree" and "P#" are not needed: the formulas hold for any p, q | B. T6 tested 2000
   non-squarefree cases in B = 900, 27000, 2520, 44100, e.g. 4/9 in base 900. Any B divisible by every numerator and
   denominator works. The supersingular "base 2*3*...*71" item (note section 6) is one instance among all finite prime
   sets.
   - What *is* needed is lowest terms. FRACTRAN multiplies by the rational number. An unreduced 6/10 fires on 35
     (giving 21), while the selector "10 | x_0" would not.
   - Say "each fraction in lowest terms".
2. **The 3x-1 cycle remark** (note line 191; typed NUMEROLOGY, correctly). It is at baseline. 10 of the 38 odd numbers in
   [17, 91] are PRIMEGAME parts. A random 7-set of them contains >= 3 parts with probability 0.257 (expected 1.84 hits).
3. **"7 = T*_5"** (note line 141) checks out. Node 7 is the end of TT5's blue line (THM-4729 L4: R_5, path plus all
   upward tiles), and R_5 is isomorphic to THM-4732's T*_5, with (H, |Aut|) = (9, 1).
4. **Scope of the session's CA check.** The 20,000-step check (driver section C, which the note calls
   `selector_automaton_check`; no function of that name exists) only meets configurations of about three base-29#
   digits. The audit extended it to 20 large random 29-smooth starts per list x 3000 steps, plus 400 random squarefree
   programs (134,111 steps) and 4000 + 3000 random formula and locality tests.
5. **Pre-registration timeline: verified.** The simulator in `adfb264520` (11:12:32) is byte-identical to the final one.
   The `.out` first appears in `0bfb7da871` (11:33:40), and the 9 min 40 s run fits in between. The protocol was
   followed, but it protects nothing here, since the statement was already a published theorem (MAJOR-1).
6. **Supersingular, Ogg, Monster and Fermat facts in note section 6: correct and correctly typed.** Ogg's fifteen; unused
   {31, 41, 47, 59, 71}; first non-supersingular output 2^37; Fermat ∩ Ogg = {3, 5, 17}; the product of the fifteen is
   31# * 41 * 47 * 59 * 71 = 1618964990108856390. The section correctly says the content is arithmetic of registers,
   not Monster structure.
7. **Repository direction D52** (`collatz_one_out_edge_systems_20260930.md` line 237: a digit system in which T is a
   sliding-block code). In the "selector" sense this is already answered by the 2026-09-26 base-6 automaton.
   THM-4852 (3) neither adds to it nor conflicts with it.
8. **Fraction anatomy.** The table (note lines 104-119), the defects, the net change {2: -1, 5: +2, 11: +2, 17: -2},
   both products, the lcm and the pure fractions are all exactly right (section 2).
9. **Clean checks requested by the brief.**
   - In a state-free configuration no rule before the last three can fire: A = 17/91 needs 13 as well as 7; B, C, I
     need 17; D, H need 19; E, K need 11; F needs 29; G needs 23; J needs 13. The audit simulator confirmed that only
     rules L, M, N fire in state-free configurations at every step to k = 2500.
   - b <= a holds at every entry: b = D(n) - 1 <= n - 1 < n = a (n = 1: b = 0). Equivalently, Guy's (15/2)-count
     a - b is >= 1 in every stretch (check (b)).
10. **Driver docstring.** The session's driver docstring (`primegame_selector_automaton_metagraph_20261031.py` line 3)
    says "the first 40 primes"; the code checks the 62 primes below 300.

## 2. Recomputed numbers

| quantity | session | audit (independent code) | status |
|---|---|---|---|
| first phase entries (r2,r3,r5,r7) | (1,0,0,0) ... (10,0,0,4) | same, k = 1..12 | OK |
| phase law at every entry, both lists | n <= 2500 (r2 = index not checked in C) | k <= 1000 and k <= 2500, r2 = k checked; negative control fires | OK |
| rule word of each state-free stretch | stated | owner (15/2)^a (1/7)^b (55/1); Guy (15/14)^b (15/2)^(a-b) (55/1); then 3^a 5^(a+1) 11; all k <= 2500 | OK |
| stateful configurations | (implicit) | exactly one state prime, exponent 1; fired rule in A-K; all steps to k = 2500 | OK |
| primes emitted | 367 below 2500, identical | 367, identical, both lists | OK |
| steps to entry 1000 | - | owner 1,344,064,699; Guy 1,343,900,918 | new |
| steps to entry 2500 | - | owner 20,961,560,628; Guy 20,960,532,292 | new |
| per-phase steps | - | Guy = Heinz's f(k); owner = f(k) + D(k-1) - 1, every k in 2..2500 | new |
| extra steps of the 1/7 list | - | sum_{k<1000} (D(k)-1) = 163,781; sum_{k<2500} = 1,028,336 | new |
| OEIS step counts | - | Guy = A007546 (fast), owner = A007547 (slow), 34/34 terms each | naming OK |
| exponent changes, 14 rules | table | identical | OK |
| defects F-U-S | 0/-2, -3/-2, ... | identical; balanced shapes p, p^3, p^2qr (brute force N < 2*10^5) | OK |
| owner net change / product | {2:-1,5:+2,11:+2,17:-2}; 55^2/(2*17^2) = 3025/578 | identical | OK |
| Guy product | 15/2 times owner = 45375/1156 = 3*5^3*11^2/(2^2*17^2) | identical | OK |
| lcm of denominators | 29# = 6469693230 (both) | identical (also procgen L1, 2026-09-24) | OK, credit |
| pure fractions | owner {1/17, 1/7, 55/1}; Guy {1/17, 55/1} | identical | OK |
| flow graph (owner) | 24 arcs, 19 pairs, no loops; (out,in) degrees | identical | OK |
| flow graph (Guy) | - | 26 arcs, 21 pairs | new |
| metagraph n = 5 | 8 blue, 24 black, 18 pairs, loops at 5 and 8, pure-blue 0, 3, 9 | identical (parsed from THM-4729) | OK |
| tilings per node | 1,1,2,3,5,9,9,10,11,13 | identical (= multigraph degrees, sum 64) | OK |
| co-occurrence embeddings | - | 70 adjacency-preserving bijections; 0 match {7, 17, 55} | new |
| rule-succession graph (owner) | - | 14 rules, 23 arcs, 18 pairs, loops at L, M; 10 macro-rules: 13 pairs, not comparable to S | new |
| CA digit formulas | proof sketch + 20,000 steps per list | 4000 formula cases (B up to 20 primes, x up to 25 digits); 3000 locality perturbations (0 violations at distance >= 2; radius 0 fails in 2694); 400 random programs; large PRIMEGAME starts | OK |
| CA as a global map | "radius-1 CA" | not shift-commuting; digit 6 depends on x_0 | MINOR-4 |

## 3. Prior-art log

| source | what it contains | consequence |
|---|---|---|
| Conway 1987, FRACTRAN (Open Problems in Communication and Computation, pp. 4-26; scan read, pp. 4-6, 9-10, 25-26) | Theorem 1 = the owner's list; section 5 / Figure 2 = the phase law with the reset L^n M^(d-1) N; section 11 = FRACTRAN games are Collatz functions (g(N)/N periodic), plus a remark on segments selected by the last decimal digit | MAJOR-1, MINOR-1, MINOR-6 |
| Guy 1983, Math. Mag. 56, 26-33 | the 15/14, 15/2, 55/1 list; its proof that it emits the primes (JSTOR, not accessed; list confirmed via OEIS A202138 and Wikipedia's footnote) | MINOR-1 |
| Conway-Guy 1996, The Book of Numbers, p. 147 | the 15/14 list (per Wikipedia) | MINOR-1 |
| Wikipedia "FRACTRAN" (fetched 2026-10-31) | PRIMEGAME = 1/7 list, attributed to Conway 1987; the 2^n 7^m -> 2^(n+1) 7^(k-1) law; the 15/14 variant "a little faster" (Guy 1983); Havil 2007 for a step-by-step explanation | MAJOR-1 |
| OEIS A007542, A202138, A203363 | 1/7 list; Guy's version differs from term 132 on; A202138 comment: Guy 1983 ends 15/14, 15/2, 55 | MINOR-1 |
| OEIS A007546 / A007547 | step counts, fast = 15/14 list / slow = 1/7 list; Heinz's 2009 closed form (encodes the phase law) | MAJOR-1; reproduced exactly |
| Kari 2012 (DLT, LNCS 7410) | 3x/2 as a CA in base 6 | cited by the session |
| Kari-Kopra 2017, RAIRO-ITA 51(4) 191-204, arXiv:1710.05737 | CA multiplying by p/q in base pq | MAJOR-2 |
| Cloney-Goles-Vichniac 1987, Complex Systems 1, 349-360 | 3x+1 as a quasi cellular automaton (base 2) | MAJOR-2, MINOR-4 |
| Korec 1992 (Math. Slovaca); Bruschi 2005 (arXiv:nlin/0502061) | CA forms of the 3x+1 map | background |
| web search "FRACTRAN" + "cellular automaton" | no paper presenting general FRACTRAN programs as (selector) CAs | extension unpublished but routine |
| repo `collatz_oscillation_gilbreath_20260926.md` section 3 | base-6 Collatz = radius-1 CA with one global bit from the lowest digit (PROVED, audited) | MAJOR-2 |
| repo `procgen_atlas_20260924.out` L1 (line 716) and `_primegame.py` | lcm = 6469693230; both lists to 15 primes, step counts = A007546/7; naming Guy 1983 / Wikipedia-OEIS; "Conway's theorem" CITED | MINOR-6; naming confirmed |
| repo `aliquot_66_five_fifteen_20261028/PROOFS.md` line 50 | the F = U + S convention | MINOR-7 |
| repo THM-4729 (n = 5 table), THM-4732 (T*_5), THM-4733 (no F = U + S claim) | metagraph data, node labels | MINOR-5, MINOR-7, NOTE 3 |

## 4. The phase law for both lists, all n (adapted from Conway 1987, Figure 2)

**Rule priority.** Rules A-K are common to both lists. In a configuration with exactly one state prime s, the first
applicable rule is the one below. Rules that need an absent state prime cannot fire, and among the rules with s in the
denominator the data guards decide, in list order.

| state | first applicable rule |
|---|---|
| 11 | E (29/33) if r3 > 0, else K |
| 29 | F |
| 13 | A (17/91) if r7 > 0, else J |
| 17 | B (78/85) if r5 > 0, else C (19/51) if r3 > 0, else I (1/17) |
| 19 | D (23/38) if r2 > 0, else H |
| 23 | G |

**One state prime at a time.**
- Rules with a state prime in the denominator have exactly one in the numerator, except I, which exits to the state-free
  set. N creates state 11 from the state-free set.
- So stateful configurations carry exactly one state prime, with exponent 1.
- With no state prime, only the last three rules can fire.

**Loops.**
- (i) In state 11 with r3 = c: (EF)^c K reaches state 13 with r7 += c, r3 = 0.
- (ii) In state 13 with r7 = d, r3 = 0, r5 = R >= d: (AB)^d J reaches state 11 with r7 = 0, r3 = d, r5 = R - d,
  r2 += d. Then (i) restores r7 = d.
- (iii) In state 13 with r7 = d, r3 = 0, r5 = r < d: (AB)^r A reaches state 17 with r5 = 0, r3 = r, r7 = d - r - 1,
  r2 += r.
  - If r = 0, I fires and gives the state-free configuration 2^(r2) 7^(d-1).
  - If r > 0, C, then (DG)^(r2) H, then (EF)^(r-1) K reach state 13 with r5 = r2, r2 = 0, r7 = d - 1, r3 = 0.

**One trial division.** Start in state 13 with (r2, r3, r5, r7) = (0, 0, N, d), 0 < d < N.
- After floor(N/d) rounds of (ii) + (i), loop (iii) applies with r = N mod d.
- At that point r2 = N, so the outcome is (0, 0, N, d - 1) in state 13 if d does not divide N, and 2^N 7^(d-1) if it
  does.
- No state-free configuration occurs before the exit.

**Reset.** From 2^n 7^b with b < n, the last three rules give 3^n 5^(n+1) 11:
- Conway / owner: (15/2)^n (1/7)^b (55/1);
- Guy: (15/14)^b (15/2)^(n-b) (55/1), which needs b < n.

Then (i) gives state 13 with (0, 0, n+1, n).

**Induction.** Entry 1 is 2 = 2^1 7^0. If entry n is 2^n 7^(D(n)-1), then D(n) - 1 < n, so both resets apply. The trial
divisions d = n, n-1, ... of N = n+1 stop at the largest proper divisor D(n+1). Entry n+1 is therefore
2^(n+1) 7^(D(n+1)-1), with no state-free configuration in between. ∎

Consequences:
- The two lists coincide at every phase boundary for all n.
- They emit exactly 2^p, p prime, in order: D(n) = 1 iff n is prime.
- Heinz's per-phase count follows. For N = n + 1:
  - Guy's reset costs N steps and (EF)^(N-1) K costs 2N - 1.
  - A failing divisor d costs 6N + 2 + 2 floor(N/d) steps.
  - The successful divisor D = D(N) costs 4N + 2N/D + 2.
  - Summing gives f(N) = N - 1 + (6N+2)(N - D) + 2 sum_{d=D}^{N-1} floor(N/d). The owner's reset costs D(n) - 1 more
    steps.

## 5. Rerun

```text
cd 04-computation/experiments/primegame_selector_automaton_20261031_audit_S70
gcc -O2 -o audit_primegame_sim audit_primegame_sim.c
./audit_primegame_sim 1000               # ~25 s  -> audit_primegame_sim_1000.out (ALL CHECKS PASSED)
./audit_primegame_sim 2500               # ~6 min -> audit_primegame_sim_2500.out (ALL CHECKS PASSED)
./audit_primegame_sim 200 control        # negative control -> audit_primegame_sim_control.out (must FAIL at phase 97)
python3 audit_selector_ca.py             # ~10 s  -> audit_selector_ca.out (ALL CA CHECKS PASSED)
python3 audit_anatomy_flow_metagraph.py  # ~8 s, needs networkx; parses THM-4729 -> audit_anatomy_flow_metagraph.out
```

The Conway 1987 scan was read from the public Gwern mirror of the Springer chapter (`gwern.net/doc/cs/computable/1987-conway.pdf`);
the page images were rendered outside the repository and are not committed.
