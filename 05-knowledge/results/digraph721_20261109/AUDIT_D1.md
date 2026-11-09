# Audit D1: mac-mini-2026-11-09-digraph721 (THM-5055, HYP-9625)

VERDICT: PASS WITH CORRECTIONS

Independent adversarial audit, 2026-11-09 (auditor: an independent Claude subagent on mac-mini, working only in this
worktree; no commits made).

**Versions read.**
- I read PROOF.md, README.md, THM-5055, HYP-9625, PREREGISTRATION.md and all scripts and outputs at the pre-audit
  checkpoint 153409f791.
- I re-read them at checkpoint 2 (c11abf3a22, 13:18), which landed during the audit. It added Part II (THM-5056) to
  PROOF.md, changed HYP-9625's status to "PROVED by THM-5056", and gave `wcheck.c` an odd-parity set "O".
- Part I of PROOF.md (sections 0-7) and the THM-5055 canon file are unchanged between the two checkpoints, apart
  from the PROOF.md header.
- **All line numbers below refer to c11abf3a22.**

**State of the worktree at the time of writing.** HEAD was 5dfac0a1ba, plus uncommitted edits by the session owner.
Those edits apply corrections from opus's audit of THM-5040 (a41a48937f) and from D2. Each finding says what they
already address.

**Scope.** THM-5056 (Part II) is audit D2's object; this audit touches it only where HYP-9625 depends on it.

**Method.**
- Every lemma was re-derived by hand.
- Every computation was redone with programs written for this audit (`04-computation/experiments/digraph721_20261109/audit/`;
  standard library only, no third-party code). They share no code with the session's.
- The session's own pipeline was also rerun and compared byte for byte.

## Overall

- **The proof of THM-5055 is correct; no gap was found.**
  - Lemma 4's evenness argument is sound. The hub is the unique vertex of maximum degree, so it is fixed. The
    excluded vertex u is then fixed. The twin permutation of D is even. The restriction to V(G) is an automorphism
    of G (in fact of Aut(G)_u).
  - R = PH(pendant, 1) meets Lemma 4.2's l = 1 hypothesis for every good graph: Seed2 on >= 3 vertices implies no
    universal vertex.
  - Both R and PH(3, 2) return good graphs, and the offset bookkeeping is right.
  - Lemma 3's co-component argument is right (for n >= 1).
  - Proposition 3.3's classification is complete.
  - In Lemma 8, every case A1/A2/A3/B1/B3 is right. The cases are exhaustive and all bounds hold: 1 <= i <= J in A3;
    i <= J or i = J + 3 in B3; J - 1 >= L(i) + 8; J >= L(J + 3) + 8.
  - The four ranges of section 5 partition the odd k <= M(n).
  - The interleaved induction is well founded. Lemma 8 for j <= n + 1 calls the theorem only at orders
    <= L(n + 1) - 1 < n.
- **Every computational input of the proof is independently certified**, without relying on the census:
  - F(n) is contained in E(n) for n <= 12 (4,028 pairs (n, k));
  - Lemma 7 at n = 8, 9;
  - the Lemma 8 base for all 2,048 odd j < 4096.

  All of these are witnessed by explicit graphs that my own verifier checks. My own closure, started from my own
  brute-force base on <= 8 vertices, reproduces E(n) = F(n) for every n <= 26.
- **The exclusions hold.** Section 3 was confirmed by exhaustive enumeration of all graphs on <= 9 vertices,
  Proposition 3.3 over its whole candidate class for odd n <= 15, and Proposition 3.2 by an independent census of
  small d-values.
- **HYP-9625's proved reduction (items 1-4) is correct**, and all its named numbers check out.
- **No MAJOR findings.** There are nine MINOR findings, none touching the theorem or its proof:
  - an invalid sentence in the justification of Remark 7.1 (the conclusion is true);
  - a false interpretive claim ("top-band values need prime graphs") and an inaccurate "exactly" in the canon's
    Reading;
  - one inaccurate NEW item relative to THM-5040;
  - an inherited misattribution of all-graphs completeness (Alikhani 2013 instead of Brouwer-Csorba-Schrijver 2009);
  - HYP-9625's PROVED status, which rests on theorems whose audits are not yet applied;
  - wording, credits and reproducibility details.

## Findings

### 1. MINOR: Remark 7.1's justification uses a wrong bound; the conclusion is correct

**Location.**
- PROOF.md lines 293-300 (section 7, remark 1, "They cannot suffice for all n").
- Relied on by THM-5055 statement 5 (canon line 93), the canon title, and README Results 5 (line 55).

**Issue.**
- The sentence "The other moves produce, at order n, only values below 3 * 2^(n-2) or at least 2^(n-1) + 1" does not
  support "only through R". The target 2^(n-1) - 81 is itself below 3 * 2^(n-2) for every n.
- The step showing that, for even n >= 84, 2^(n-1) - 81 can enter E_n only through the pendant set is missing.
- Statement 5 is labelled FINITE-EXACT, but "must fail at even n >= 84" is an argument about the closure for every n.
  It rests on one finite census fact.

**Fix.** Replace the two bullets with the following argument, writing v_n = 2^(n-1) - 81.
- **(i) Pendant bound.** A graph on m vertices with a pendant p (neighbour q) has d <= 3 * 2^(m-2), because every set
  avoiding {p, q} fails to dominate.
  - So in the closure, padding into P_n or S_n gives values <= 3 * 2^(n-3). That is below v_n for n >= 10.
  - The cone over G u K_1 gives >= 2^(n-1) + 1.
  - Hence, for n >= 10, v_n is in S_n iff v_(n-2) is in S_(n-2).
- **(ii) Entry into E_n.** For even n >= 84, v_n > M(n-1). Cones give >= 2^(n-1) + 1, and E_2-joins give
  >= 3 * 2^(n-2) - 1. So v_n is in E_n iff v_n is in W_n iff v_n is in S_n.
- **(iii) The chain.** For even n the chain ends at v_8 = 47, which is not in S_8: 47 is prime, and its three
  connected even realisers have no pendant.
- **Typing.** Type statement 5's second half as "PROVED from the FINITE-EXACT census fact (iii)".

**Evidence.**
- `audit/aw_four.py` re-implements the four-move closure from my own base.
  - It is exact for n <= 26.
  - The pendant bound holds at every order.
  - v_n is outside S_n for every even n in [8, 26], and inside S_n for every odd n in [9, 26].
- `audit/aw_smalld.out` confirms the census fact (item H below).

**Working tree.** Not yet addressed.

### 2. MINOR: "top-band values need prime graphs" is false; the canon's "exactly" list is inaccurate

**Location.**
- PROOF.md lines 307-308 (remark 2, last bullet).
- README lines 53-54 (Results 4).
- THM-5055 canon lines 101-103 (Reading).

**Issue.**
- **(a) Unions reach some top-band values.** Values just below a power of 2 can be disjoint unions. Each of the
  following is P_3 u H, with H the session's own verified witness. All four are even, on their least orders
  (`audit/aw_union.py`, `aw_union_topband.out`):
  - 505 = 2^9 - 7 = 5 * 101, on 10 vertices;
  - 1015 = 2^10 - 9 = 5 * 203, on 11 vertices;
  - 2045 = 2^11 - 3 = 5 * 409, on 12 vertices;
  - 4085 = 2^12 - 11 = 5 * 817, on 13 vertices.

  The same works for every top-band k with 5 | k and k/5 in F(n - 3). The opus audit of THM-5040 (M2) found 505
  independently.
- **(b) The "exactly" list is wrong under either reading.** The canon says the hub moves "supply exactly the values
  that no cone or join can reach at any order: those just below a power of 2, and the shifts 2^(n-1) + {3, 7, 13}
  and 3 * 2^(n-2) + {1, 5, 11}".
  - 3 * 2^(n-2) + {1, 5, 11} are cones: over the Lemma 7 graphs, exactly as section 5, range 4 builds them.
  - Under the reading "cone over a padded graph", the list omits 2^(n-1) + x with x in the top band of order n - 2.
    These are cones over Lemma 8 graphs (section 5, range 3, second bullet).
- **(c) Off-by-one bound.** "A join d >= 3 * 2^(m-2)" is off by one: an E_2-join has d >= 3 * 2^(m-2) - 1. For
  example, C_4 = E_2 v E_2 has d = 11.

**Fix.**
- Say "these values are never cones over padded graphs or E_2-joins over realisers of lower order".
- Drop "exactly", or list all families.
- Use ">= 3 * 2^(m-2) - 1".

**Working tree.** (a) and (c) are fixed in PROOF.md and README, and the canon Reading gained the union caveat.
(b) is still open: the canon keeps "exactly ... 3 * 2^(n-2) + {1, 5, 11}".

### 3. MINOR: the NEW list overstates the novelty of the hub move relative to THM-5040

**Location.**
- THM-5055 canon lines 28-29 ("NEW here, relative to THM-5040: the hub-move lemma for general r and l, with its
  evenness criterion").
- PROOF.md section 0, lines 20-27 (the "shared" list).

**Issue.** THM-5040 already has the hub move with an evenness criterion, in two cases.
- **Jump family.** THM-5040's Lemma 8 jump family J(H, z) is exactly PH(comp(H), z, k) with l = r = k, since
  r = |N_comp(H)[z]| = rho - h = k.
  - Indeed d(J) = d(comp H) + (2^k - 1) 2^rho = 2^(m-1) - D(H).
  - It uses the same unique-maximum-degree argument.
  - Its evenness criterion is an iff through Aut(H)_z, equivalently through the evenness of Aut(comp H)_z by B2.5.
    That is sharper than Lemma 4.2's sufficient condition "G even".
- **SR.** THM-5040's SR is the case (r, l) = (2, 1).

What is new here:
- l independent of r, in particular PH(3, 2);
- the class of good graphs closed under R and PH(3, 2), with the parity-switching shift 3;
- the seed lemma;
- the R-only obstruction;
- the HYP-9625 reduction.

**Fix.**
- Reword the NEW item: "the hub move with r and l independent (THM-5040 has (2, 1) = SR and l = r = its jump
  family, with an iff criterion via Aut(H)_z)".
- Add the jump family to the shared skeleton in section 0.
- Optionally, note in Lemma 4.2 that only Aut(G)_u needs to be even, since the proof uses only
  b|V(G) in Aut(G)_u.
- The value formula is an instance of THM-4996 Lemma B3's vertex-addition count. Label it elementary, as THM-5040
  does for its Lemmas 7-8.

**Working tree.** Not yet addressed.

### 4. MINOR (credit): completeness of d over all graphs is misattributed to Alikhani 2013

**Location.**
- THM-5055 canon line 35 ("Prior art: completeness of d over all graphs is S. Alikhani, Graphs Combin. 29 (2013)").
- README line 94 (credits).

**Issue.**
- Opus's audit of THM-5040 (a41a48937f, `05-knowledge/results/egcp_20261109/audit/AUDIT.md`, Finding 1) records the
  source as Brouwer-Csorba-Schrijver 2009, Proposition 1.2: every odd m with 0 < m < 2^n is d of a graph on n
  vertices. That covers completeness and conciseness.
- It also records that Alikhani 2013 concerns D(G, -1).
- THM-5055 inherited the Alikhani attribution from THM-5040's draft.
- I could not check the external sources myself (no web access), but the recorded evidence is specific: a
  PDF glyph check, and a restatement as Alikhani-Fatehi-Klavžar, arXiv:1512.07514, Prop. 2.1(iii).
- THM-5055 claims no novelty for this result, so it is a misattribution, not an overclaim. The sharp orders for even
  graphs remain new (THM-5040 and THM-5055).

**Fix.** Credit Brouwer-Csorba-Schrijver 2009 Prop. 1.2, citing MISTAKE-747.

**Working tree.** Applied in the canon and README. PROOF.md has no such line.

### 5. MINOR: HYP-9625's PROVED status rests on theorems whose audits are not yet applied

**Location.** HYP-9625 status, line 5, and the hypotheses INDEX row, line 2.

**Issue.**
- The status says "PROVED by THM-5056"; the working tree adds "and independently by THM-5042".
- THM-5056's audit D2 reports PASS WITH CORRECTIONS (`AUDIT_D2.md`), but it is uncommitted and not yet applied at
  the time of writing.
- THM-5042's own status still says "independent audit pending".

This audit checked:
- HYP-9625's reduction, items 1-4, which are correct;
- the named values;
- THM-5056's four seeds and their PH(w, 2) images: d = 13 -> 69, 33 -> 257, 101 -> 549, 165 -> 1061, all odd
  (`audit/aw_hyp_verify.out`);
- independently, that odd graphs on <= 9 vertices realise every odd k < 512 except 1, 5 and 37 (item I below).

It did not check THM-5056's or THM-5042's proofs.

**Fix.** Write "PROVED by THM-5056 (audit D2: PASS WITH CORRECTIONS) and THM-5042 (audit pending)", or the
equivalent, until both audits are applied. Mirror this in the INDEX row.

### 6. MINOR: HYP-9625 wording

**Location.** HYP-9625 items 2, 4 and 5 (lines 32, 43, 46).

**Issues and fixes.**
- **Item 2.** "For every graph G on m vertices, K_2 v G ... d = 3 * 2^m + d(G)" needs m >= 1: K_2 v (empty graph) =
  K_2 has d = 3, not 4. The value 3 = d(K_2) is covered separately. Part II's O1(c) in the working tree already says
  m >= 1.
- **Item 4.** "pad K_2 u G_115 to 10 vertices, where G_115 is any graph with d = 115": G_115 must have <= 8 vertices,
  for example the 7-vertex threshold graph.
  - Verified: 1369 = d(K_1 v (K_2 u T_115 u K_1)), and the graph is odd.
  - 25 = d(K_2 v E_3), 185 = d(K_1 v ((K_2 v P_4) u K_1)) and 517 = d(paw u T_47) are also verified odd.
- **Item 5.** "the chain avoids O" should read "avoids values realised by odd graphs within the available order".
  The cone over a padded realiser of x needs a realiser on <= L(p) - 1 vertices. This is now superseded by THM-5056;
  it is a wording fix only.
- **Arithmetic.** All of HYP-9625's arithmetic is right:
  - 1061, 2053, 4133 and 65537 are prime;
  - 517 = 11 * 47, 1025 = 5^2 * 41, 4097 = 17 * 241, 4613 = 7 * 659;
  - 37^2 = 2^10 + 3 * 115.

**Working tree.** Not yet addressed in HYP-9625.

### 7. MINOR: Lemma 3 needs n >= 1; value-formula credits are missing from PROOF.md

**Location.** PROOF.md line 62 (Lemma 3) and lines 59-60 (Lemma 2).

**Issue.**
- For n = 0, Lemma 3's formula gives 2, but d(E_2) = 1.
- The value formulas of Lemmas 2 and 3 are the x = 1 case of the Alikhani-Peng join formula
  d(A v B) = (2^|A| - 1)(2^|B| - 1) + d(A) + d(B). Lemma 3 is also THM-5040 Lemma 5 with t = 1.
- PROOF.md does not cite Alikhani-Peng (the canon and README do).

**Fix.**
- "for every graph G on n >= 1 vertices".
- Add the credit in PROOF.md.

**Working tree.** Not yet addressed.

### 8. MINOR: small precision and typing points

- **PROOF.md line 36.** "k is in F(L + 1), since 2^L <= M(L + 1)" fails for L = 1, where M(2) = 1. It is harmless,
  since then k = 1 is in F(1). Add "(L >= 2)".
- **PROOF.md line 89.** "(and n + l >= 3)" is automatic when l >= 2. Drop it.
- **THM-5055 canon line 10.** "Its base is a finite closure: E(n) = F(n) for n <= 26 ..." overstates what the
  induction uses.
  - The proof needs only three inputs: F(n) contained in E(n) for n <= 12; Lemma 7 at n = 8, 9; and the Lemma 8 base
    (L(j) <= 12, orders <= 21).
  - All three are certified by the 6,147 explicit, independently verified witnesses, whatever the census's
    exhaustiveness.
  - The census enters only the exclusions (Proposition 3.2) and Remark 7.1.
  - The closure to 26 is supplementary. Stating this makes the dependency structure, and its strength, explicit.

### 9. MINOR: reproducibility of section 6

- **Negative controls.** "Negative controls reject K_2, P_4 and a wrong value" (line 286) has no committed script or
  output. I reproduced it: the session's `wcheck` and my `aw_verify` both reject the three controls in
  `audit/negctl.txt`. Commit the controls and their output.
- **Compile commands.** The table lacks the compile commands for `evenclos`, `closure4` and `wcheck`, whose binaries
  are gitignored.
- **witness.out.** The committed `witness.out` contains the stderr line "base representatives: 2436". The documented
  command, which redirects stdout only, does not produce it. Use `2>&1` or drop the line.
- **Independence.** "Reproduced by an independent numpy implementation" (canon line 13; README line 73) is fair for
  the code. But `witness.py` uses the same census-derived base and the same move set. Say so; this audit supplies a
  closure that is independent in both base and code.

## What was independently reproduced

All commands are in `audit/run_audit.sh`; hashes are in `audit/sha256_audit.txt`.

**A. The session's pipeline (PROOF.md section 6), rerun and compared byte for byte.**
- `dom_even 511 16`: 21 s. Its output is identical to `dom_K511.out`, with SHA-256 6307dfa1... as in `sha256.txt`.
  The class and parity counts match `dom_K511.err`.
- The following are byte-identical:
  - `mkbase.py` gives `base9f.txt` and `base9.txt`;
  - `evenclos` gives `evenclos_full.out`, `evenclos_noR.out` and `gamma_dump.out`;
  - `closure4` gives `closure4.out`.
- `seedcases.py` gives the same counts: A1 7,140, A2 9, A3 19, B1 7,145, B3 23.
- `witness.py` gives `witnesses.txt` with an identical SHA-256 (5013c326...). `wcheck` reports ok = 6,147.
- `lemmacheck.py` runs 706 tests and `chaincheck.py` 217, each with 0 failures.
- The Euler transform of the connected even counts (1, 0, 1, 2, 8, 43, 398, 6,413, 184,596) is
  1, 1, 2, 4, 12, 56, 456, 6,880, 191,536 (A000568).

**B. An independent verifier, `aw_verify.c`.**
- How it works:
  - d by inclusion-exclusion over undominated sets;
  - parity from a strong generating set of Aut(G), via a stabiliser chain whose base comes from equitable
    refinement, with complete individualisation-refinement searches;
  - flags from degrees.
- Self-tests, all with 0 disagreements:
  - against brute force over all n! permutations and direct subset enumeration, on 30,000 random graphs (<= 8
    vertices, with twin constructions);
  - on 200 structured symmetric graphs (<= 11 vertices: cocktail parties, K_{r,r}, paths, cycles, Petersen, Q_3,
    hub images and others, each also randomly relabelled).
- **All 6,147 witnesses pass:** d, evenness and flags.
- **Coverage, checked independently of the tags (`aw_coverage.py`):**
  - every (n, k) with n <= 12 and k in F(n) (4,028 pairs) has a verified even witness on <= n vertices;
  - all 12 Lemma 7 graphs are present and good;
  - every odd j < 4096 has good witnesses at two consecutive orders n0, n0 + 1, with n0 <= L(j) + 8 (maximum 8).
- So the proof's computational inputs are certified.

**C. An exhaustive enumeration, `aw_enum.c`.**
- It covers all labelled graphs on <= 9 vertices with non-increasing degree sequence (every isomorphism class
  appears), by Gray code:
  - n = 8: 268,435,456 labelled graphs, 602,544 representatives;
  - n = 9: 68,719,476,736 labelled graphs, 37,693,734 representatives, about 5 minutes.
- E(n) = F(n) exactly for every n <= 9:
  - the maximum is M(n);
  - X(5) = 25, X(7) = 119 and X(9) = 501 are not realised;
  - 3, 7 and 13 are not realised.
- The complete set of flag tuples of even graphs equals `base9f.txt` (1,398 tuples), so the base is independently
  reproduced.

**D. The auditor's own closure, `aw_closure.py`.**
- It uses big-int bitsets and the section 2 and section 4 moves, started only from the auditor's base on <= 8
  vertices (no 9-vertex census).
- Results:
  - E(n) = F(n) for every n <= 26 (0 missing, 0 extra);
  - Gamma(+c) for c in {1, 3, 5, 7, 11, 13} is identical to the session's;
  - the Lemma 8 base holds for all 2,048 odd j < 4096, with maximum n0 - L(j) = 8.
- A sample of 989 claims was reconstructed as explicit graphs and verified by brute force, with 0 bad. It covers:
  every top-band value, every X(n-1), every 2^(n-1) + x with x delicate, and every 3 * 2^(n-2) + {1, 5, 11} at each
  order, plus random values, the offsets, and seeds, up to 26 vertices.

**E. Constructions exactly as the proof prescribes (`aw_construct.py`).**
- It built 436 graphs on 13-22 vertices and asserted every stated hypothesis along the way:
  - for each order n = 13..20, every value in ranges 2, 3b, 3c and 4b of section 5, plus 25 random values;
  - Lemma 8's inductive cases at J = 13: every A2, A3 and B3 case, plus 6 A1 and 6 B1 cases, at orders 21 and 22,
    with R / PH(3, 2) chains.
- All 436 pass the brute-force check of d, evenness and flags.

**F. Arithmetic (`aw_arith.py`).**
- The Lemma 8 case analysis holds for 13 <= J <= 400, and the section 5 ranges for 13 <= n <= 400.
- The four inequalities hold up to 10^6.
- There are 0 failures, and the case counts for J = 13-15 match `seedcases.out`.

**G. Proposition 3.3 (`aw_prop33.py`).**
- For odd n = 5, 7, ..., 15 I took every graph whose complement has maximum degree <= 2 (the only candidates, by the
  counting step): 11, 29, 70, 156, 334 and 686 graphs.
- The realisers of X(n) number 4, 6, 11, 21, 31 and 50, and all are odd.
- Those with <= 1 universal vertex have exactly the four listed shapes (three at n = 5).

**H. Small-d census (`aw_dump.c`, `aw_smalld.py`).**
- The census is complete by Ore's bound: a connected graph on n vertices has d >= 2^(n/2 + 1) - 1, so connected
  realisers of d <= 62 have <= 9 vertices.
- Results:
  - 3, 7 and 13 are realised only by K_2, K_3 and K_4 - e, all odd;
  - 5 only by P_3;
  - 37 by 4 connected graphs, all even;
  - 11, 23 and 47 have 1, 2 and 3 connected even realisers, none with a pendant;
  - 25 has 5 connected realisers, all odd.
- This matches `dom_K511.out` class by class.

**I. HYP-9625.**
- `aw_oddvals.c` is exhaustive over n <= 9. Odd graphs on n vertices realise every odd k < 2^n except 1, 5 and 37 for
  6 <= n <= 9 (except 1 and 5 for n = 3..5). So below 512, O misses exactly 1, 5 and 37.
- The explicit odd realisers of 25, 185, 517 and 1369 and the arithmetic facts were verified (Finding 6).

**Not checked.**
- THM-5056's and THM-5042's proofs (audit D2 and opus's pending audit).
- External bibliographic details: no web access was used, and Finding 4 relies on the evidence recorded by opus's
  audit.
- THM-5055 does not claim priority over THM-5040. It calls itself a second, independent proof, cites c46a3504dc and
  MSG-1988, and credits THM-5040 for the sharp formulation and the exclusions. Its credits for von Brömssen,
  RPGFD 2023, Brouwer-Csorba-Schrijver, Alikhani-Peng and Chan-Pak Conjecture 6.2 are present in the canon and
  README. Findings 3 and 4 are the only credit corrections.
