# Merged ledger: A000568, the tiling metagraph, the {7,21} Hamilton spectrum, and the module-reversal frontier

kind-pasteur-2026-11-04-S9, 4 November 2026. This note is a ledger, not a new theorem.
- **Statuses** are quoted from the cited files, read at origin/main 89c0723a76 and 92aaef65f8. The newest status wins;
  superseded items are marked.
- **New result:** THM-4864, `reversal_cube_mod4_anchor_20261104/`.
- **Residual:** HYP-9456.

## 0. One map, four questions

**Later November4 integration:** [the recovered even-graph programme](even_graph_response_merge_20261104/MERGE.md)
records already-constructed SCC/Lyndon, primitive-support, prime-local and
prime-circulant type bijections, each with its preserved data. The remaining
natural-map question must not erase these constructions or imply labelled
equivariance. Its marked Fourier specialization keeps the complete cubic and
quartic observer. Separately, reviewed Session189 commits `c5725b3a7f` and
`c87ca88cff` advance the original-module XOR frontier below to **prime order
at least nine**; see [the proof](session189_response_frontier_20261104/frontier_nine/PROOF.md)
and [independent review](session189_response_frontier_20261104/nine_integration_audit/REVIEW.md).
The older rows retain their dated input scopes, not a competing current boundary.

Send a named tournament T to three things: its isomorphism class, its Hamilton count H(T), and its two-port response
F_T. The response refines the count: F_T(0,0) = H(T). The four threads ask about this map.

| thread | question about the map | current boundary |
|---|---|---|
| A000568 | How many classes are there? | Closed formulas and specified type bijections are PROVED. A general natural-map criterion remains OPEN; labelled equivariant orbit bijections are excluded from five. |
| Tiling metagraph | How are classes wired by complement lines of tilings, with fibres H/\|Aut\|? | Parity and line counts are PROVED. Connectivity and a compressed recurrence are OPEN. |
| {7,21} spectrum | Which odd values does H take? | 7 and 21 are excluded (PROVED). Completeness is OPEN (Chan–Pak Conj. 6.2). |
| Response frontier | When does F separate T from an original-module XOR target? | Highest changed nodes must now be prime of order ≥ 9 (reviewed S189 integration); earlier transitive-child and singleton exclusions retain their stated scopes. |

## 1. Counting layer (A000568 and companions)

| item | file | status as recorded | content |
|---|---|---|---|
| THM-305 | `01-canon/theorems/THM-305-*` | PROVED | v_2(A000568(n)) = (n−1)/2 for odd n ≥ 3 |
| THM-514, THM-516 | canon | PROVED | Davis/Burnside formula over odd partitions; checked against A000568 through n = 100 |
| THM-283 | canon | PROVED | self-converse counts, A002785: 2, 2, 8, 12, 88, 176, 2752, 8784 for n = 3..10 |
| HYP-3799 | hypotheses | CONFIRMED (n ≤ 5) | imports the Royle–Praeger–Glasby–Freedman–Devillers equinumerosity (arXiv:2204.01947); natural bijection OPEN |
| THM-4721 | canon | PROVED (pre-audit) | Euler(connected even graphs) = INVERT(strong tournaments); with A051337 = 1, 0, 1, 1, 6, 35, 353, 6008 this gives A000568 |
| THM-4722, THM-4723 | canon | PROVED (pre-audit) | no S_n-equivariant type bijection for n ≥ 5 |
| HYP-9340 | hypotheses | (a) REFUTED at n = 10 (C15: 32 vs 31); proved for every odd p-group | the refined equinumerosity |
| THM-4728, HYP-9341 | canon / hyp. | failures only C15 on 5+3+1^f up to n = 12; decorated cancellation OPEN | |
| trijection S3 note | `trijection_four_nodes_three_edges_20261020.md` | PROVED | marked line, cycle or star ⇔ Tour(n−1) modulo 1, Z_{n−1} or S_{n−1}; (8,4,4), (64,16,8), (1024,208,24) |
| THM-4835, THM-4836 | canon | PROVED; audited (MISTAKE-689) | Euler-characteristic trijection criterion; pair-or-triple law exactly for n ≤ 5; t_8 = 6880, s_8 = 6008 |
| THM-4845 | canon | audited | A000568 as a Gaussian average; Davis values to n = 80 |
| THM-4803 | canon | CITED + PROVED + FINITE-EXACT; all-order Aut/isotype invariance REFUTED (corrected) | A000568 from the prime cycle-index census. Boussaïri–Ille–Lopez–Thomassé: a prime tournament is determined up to duality by its 3-cycles. |
| THM-4814, THM-4816 | canon | FINITE-EXACT (4814, unaudited); PROVED (A) + FINITE-EXACT (4816) | The prime cycle-index tower: the nine order-10 types of 7,906,161 prime tournaments, and orders 11–13 from non-rigid primes. |

**Open in this layer:** a natural non-equivariant bijection with even graphs; a decorated defect cancellation; a
compressed grammar or closed evaluation.

## 2. Metagraph layer

Definitions (THM-643, THM-644, THM-1430):
- A tiling is a vector in GF(2)^m with m = C(n−1,2). Tilings are exactly switching classes.
- σ is the grid transpose and acts as converse. Blue tilings are the σ-fixed ones and lie over self-converse classes.
- A line is {t, t + all-ones}.
- The fibre of a class is H/|Aut|, and sum_classes H/|Aut| = 2^{C(n−1,2)}.

| item | status as recorded | content |
|---|---|---|
| HYP-2250 | OPEN; naive simple-graph claim refuted through n = 6 | |
| THM-643 | PROVED (T1–T5); C1 REFUTED at n = 8, 9; pure-blue-leaf prose REFUTED at n = 5 | fibres odd; pure-black ⇔ not self-converse; blue lines 2^{e−1}, black lines 2^{m−1} − 2^{e−1}, e = ⌊(n−1)²/4⌋ |
| THM-644, 646, 647, 648 | PROVED | blue ⇔ i ↦ n+1−i is an anti-automorphism; score complement; anti-Rédei; no blue self-loops at odd n |
| THM-649 | REFUTED at n = 10 | the even-n self-loop law 2^{n/2−2} |
| THM-1430 (tiling-class file; the ID is shared by two other files) | MIXED, per part | the dictionary; 272 merged nodes at n = 7 |
| THM-2444 | REFUTED (pure-blue(9) = 6, not 5) + PROVED (rigidity lemma) | superseded in formula by THM-2454 |
| THM-2454 | PROVED (2026-10-20) | strong pure-blue alphabet {1, C3, regular C5}; pure-blue(n) via A000930: 2, 1, 3, 2, 4, 3, 6, 4, 9 for n = 3..11 |
| THM-4729, THM-4731 | L1–L5 PROVED pre-audit, L4 conditional; THM-793 conjecture REFUTED at n = 9 | |
| HYP-9345, HYP-9351 | OPEN | connectivity verified for n ≤ 9; second-order port data does not determine the metagraph at n = 6 |
| session50 | exact aggregate reconstruction | predicts 6880 and 191,536 |
| session67, adaptive_metagraph_reformulation_20261103 | count-free reconstruction (exact, exponential); one-temporary-vertex native compiler PROVED | |
| session68 METAGRAPH_BOUNDARY | flux law c3(U) − c3(T) = s(first) − s(last) − 1 | composition families not closed under lines |
| THM-4857 | audited | a GF(2) automaton classifies merged nodes exactly for n ≤ 5 and fails for every n ≥ 6 |
| HYP-9455, THM-4863 (opus-S72, today) | HYP OPEN, with P1–P3 confirmed at orders 8, 9 and P4 refuted at 8; THM-4863 checkpoint pre-audit | The path-weighted deletion tower determines merged nodes at orders 8 and 9. THM-796's face row is the Hamiltonian-endpoint-weighted deck. |

**Open in this layer:** all-order connectivity; blue self-line counts (n = 12 is next); self-line closed forms; a
compressed recurrence for A000568 on the metagraph.

## 3. The {7,21} Hamilton spectrum

- **Exclusions (PROVED).**
  - Main file: `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`. Its arguments:
    - multiplicativity over strong components;
    - h(T) ≥ h(T−v);
    - Moon 1966;
    - exhaustive strong minima 9, 15, 25, 45 for n = 5..8;
    - Busch 2006 floors.
  - H = 7: THM-343. H = 21: THM-115, with its historical body.
  - The ID THM-1370 is shared with the elliptic-Jacobian theorem. MISTAKE-209 already records this. Cite the full
    filename.
  - Fibre analogue: THM-4812 proves H/|Aut| ∉ {7, 21} at every order.
- **Completeness (OPEN).** It is exactly Chan–Pak's Conjecture 6.2.
  - CITED COMPUTATIONAL EVIDENCE: Royle's MathOverflow answer (updated January 2022) reports that every odd value
    below 80557 except 7 and 21 occurs by order 11.
  - Status note: `spectrum_basins_20261020e/spectrum.md`.
- **The two-prime-lane reduction (PROVED, already in canon).** THM-4094 is independently hostile-audited. Order joins
  multiply H, 7 and 21 are omitted, and the finite strong carry atoms are 49, 63 and 343. With these, completeness is
  equivalent to strong realization of two infinite lanes: p for odd primes p ≠ 7, and 7p for odd primes p ≥ 5.
  - This session re-derived the same reduction independently and checked the closure identity
    Odd∖{7,21} = ⟨primes ≠ 7, 7p (p ≥ 5), 63, 343⟩ up to 200,000. That is a reproduction, not a new result.
- **Progress along the lanes.**
  - THM-4097 (PROVED + VERIFIED-EXACT): the order-9 strong spectrum contains every odd value from 85 to 2881.
  - THM-4102 and THM-4104 move the lane targets to 80407 and 7·11527 = 80689.
  - Global completeness remains OPEN.
- **The one-apex family is not a route.** `spectrum_basins_20261020e` §5 already shows that it misses the prime 13, so its product closure misses 13 too. Exploratory reproduction from this session: Take a transitive tournament with one apex, so
  that T − v is transitive. Then H = 1 + 2m, where m counts even subsets Y (|Y| ≥ 2) with min Y an out-neighbour of
  v and max Y not one.
  - Over at most 16 transitive vertices it misses 176 of the 401 values m ≤ 400.
  - The first misses are m = 3, 6, 10, 14, 17, …, i.e. H = 7, 13, 21, 29, 35, ….
  - Completeness therefore needs genuinely strong building blocks, not products of one-apex blocks.

## 4. The module-reversal response frontier

| item | status | content |
|---|---|---|
| S188 (`00-navigation/S188.md`) | reviewed snapshot | Prime quotients of order 3, 5, 6 are closed in every ambient. Every highest changed quotient of a surviving collision is prime of order ≥ 7. |
| `parity_profile_transfer` (S188) | reviewed | Endpoint parity p_i(Q) = S_i + E_i is a one-query separator. For all-nontrivial children the whole coarse contrast is 2Σ p_i(a_i + c_i) mod 4. Mod-4 core separators lift to every ambient. |
| E3/E4, 76b3b1a | proved | Every tournament differs from its named converse mod 4. A single module converse is never hidden. |
| `simultaneous_module_reversal_20261104` | PROVED S1–S3 | With the quotient fixed, r disjoint nontrivial module reversals are decoded by r actual queries. This is the reversal cube without a quotient converse. |
| **S189 collision (your statement)** | reviewed PASS | A **transitive** module of order ≥ 5 has its physical order fixed, up to whole converse, by its internal shared-name quartics in every fixed ambient. So a collision leaves such a child unchanged or reversed. The threshold 5 is sharp for this battery. |
| Session189 checkpoints 04–10 (today, concurrent) | reviewed, except checkpoint07 integration pending | Prime-seven and prime-eight arbitrary-child ambient converse separation (finite complete certificates, including singleton positions). Frontier order ≥ 8 (reviewed), ≥ 9 (pending). The Seidel-energy dyadic separator, and the source-relative separation of every XOR family of a dyadic TT_N thickening (checkpoint09, not the uninflated source). An all-even modulo-16 transfer. |
| **THM-4864 (this session)** | PROVED; independently audited, SOUND WITH CORRECTIONS, corrections applied (MISTAKE-703); pre-registered P1–P5 | (A) A singleton child gives a coarse mod-4 separator, for any Q and any same-size target interiors. (B) For m ≥ 2, a child with a vertex that is a source or sink of both its source and target interiors gives an anchored mod-4 separator. In particular **unchanged and reversed transitive children are both separated**, at every order and in every ambient. |
| **HYP-9456** | OPEN | Residual frontier at every highest changed prime node (order ≥ 7, p(Q) = 0): no singleton child, and no child with a common end vertex. So every child is non-transitive or a crossed TT4. Conjecture: one-child separation **modulo 8**. Modulo 4 already fails for Paley-7 with C3 children kept unchanged, and the exact gap there is 4 mod 8. |

**What separates the unchanged and reversed cases.** Contract the child minus its source b_1 to one vertex β. Then:
- **Reversal** is visible at b_1 through a weighted version of the E1 identity.
- **Unchanged** is visible through the τ-invariant functional L_{b1} + L_β + [a_{b1} a_β c_β], whose value is
  1 + p_j(Q). When p_j(Q) = 1, coarse queries already separate.

The S188 prime-eight endpoint-obstruction quotient, which is blind to every coarse mod-4 query, is separated this way
at every transitive child (P4 of the pre-registration).

## 5. Cross-thread observations (typed)

- **PROVED, shared input.** One parity fact drives three results: by Rédei, every induced subtournament has odd H.
  - It gives odd fibres (THM-643).
  - It gives the mod-4 response identities (E1, Proposition 4 of THM-4864).
  - It is why the spectrum question is posed among odd values only.
  - Nothing further is claimed to follow.
- **PROVED.** H(T) = F_T(0,0). The spectrum gaps 7 and 21 are therefore gaps in the scalar coordinate of the response.
  The response-injectivity program studies fibres of the full response over module-XOR orbits. They are different
  questions about one map.
- **ANALOGY only.** The v_2 law for A000568 (THM-305) and the mod-4 response residues (q_2 ≡ 2, q_3 ≡ 2[n odd] mod 4)
  are both 2-adic statements about tournament counts. No implication is claimed.

## 6. Ranked open queue

1. **HYP-9456: non-transitive children.**
   - **Exploratory finding.** Modulo 4 is blind on Paley-7 with C3 children kept unchanged, but the exact gap is
     4 mod 8, already on a coarse coefficient. The next layer is a modulo-8 one-child law. Template: the S188
     five-type coarse_mod8_affine law.
   - Treat TT4 crossed reorderings with the full four-name weighted model; the evidence so far has 0 blind cases.
2. **Chan–Pak 6.2.** By THM-4094 this is exactly the two prime lanes, p and 7p.
   - Needed: a strong family realizing every prime p > 80407 and every 7p > 7·11527.
   - THM-4097's one-vertex ears are the current engine. The one-apex family alone is insufficient: it misses the prime 13 (spectrum_basins §5).
3. **The metagraph recurrence for A000568**, using HYP-9455's path-weighted deck as the candidate state.
4. **A natural bijection between tournaments and even graphs** (HYP-3799). The refined version fails at n = 10.
