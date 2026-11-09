# Audit G1: genesis (THM-4985, THM-4986, HYP-9565)

Independent adversarial audit of session mac-mini-2026-11-09-genesis, at checkpoint 561b0319b6 plus the uncommitted
working tree of 2026-11-09 (PROOF.md, README.md, THM-4985, THM-4986, HYP-9565 and the HYP-9532 status edit), with the
pre-registered outcomes and the new structural census included. All code is my own, in [`audit/`](audit/). None of
the session's programs were run.

## Verdict: PASS WITH CORRECTIONS

Every computational claim I could test reproduces exactly with independent code. This covers:
- the order-12 failure of the five-operation grammar;
- the 12 extension classes and the 193 classes in all;
- the defects, creator depths and structural census;
- both band certificates and k*(15) >= 4;
- the first-creator shares;
- the tau table and the pre-registered outcomes.

The proofs of Lemmas W, B1, B2, B3 and T1, the (S-inf) theorem, the (M_eps) theorem and its corollary are correct.

Three MAJOR corrections are needed:
- **M1.** The value "max H(12) = 545,835 (OEIS A038375)" has no source. OEIS A038375 ends at a(11) = 95,095. It is used
  as FINITE-EXACT input in canon (R_H(13) = 0.359) and in several evidence statements. Relatedly, max H(13) is not the
  circulant 3,711,175: a 13-tournament with 3,719,831 paths already exists in the repo.
- **M2.** The band tournaments are THM-4880's constant-letter ribbon tournaments. THM-4880 already gives their H-series
  for s = 2, 3. This is not credited.
- **M3.** The bounded-defect question behind HYP-9565 (A) was posed by Session230 before this session. This is not
  credited.

There are 23 MINOR findings.

## What I reproduced, and how

All files are in `05-knowledge/results/genesis_20261109/audit/`. The C tools are built with `cc -O2 -o X X.c`; the binaries were removed after the runs. `aud_spectra.py` regenerates the `S_n.txt` files that the other scripts read, and `chk_census_wit.py` documents how to regenerate its census extracts.

### Validation of my own tools

| tool | validation |
| --- | --- |
| `aud_core.c`: H by subset DP; lambda by all 2^n - 2 cuts; primality and module defect by enumerating every vertex set | the DP equals brute-force permutation counts on 2,400 random tournaments with n <= 8. My decoding plus DP reproduces H for 19,353 census witness lines (W and S lines) at orders 9-11 (`chk_census_wit.py`), which confirms the code convention ('1' = i -> j, upper triangle, row-major) |
| `aud_iso.py`: colour refinement plus exact backtracking | 200 random relabellings with n <= 12; |Aut| against brute force for n <= 7; 662 random pairs with equal score sequences against brute-force isomorphism. 0 failures |
| `aud_dec.c`: fewest reversals to a non-strong or non-prime tournament, by exhaustive reversal sets of size <= 3, independent of the module-defect formula | agrees with min(lambda, formula) on 1,600 random tournaments (n = 5..8) and on 62 highly connected ones (QR7, QR11, 60 random tournaments with minimum semidegree 3 at n = 8, 9). Module defect <= lambda in every case |
| `aud_spectra.py` | S(n), n <= 11, read only from census T-lines |

### Results

| claim | my check | result |
| --- | --- | --- |
| THM-4985 (1) | S17's 7 parents decoded from the integer codes: H = 15,213, 15,561, 11,475, 15,273, 15,213, 15,273, 15,561; all lambda = 3 and not prime. All 7 × 2^11 = 14,336 extensions classified (`aud_core ext`, `ext12_audit.out`) | all 13 label categories match `ext12.out` count for count (8,808 FC+FCV, ..., 36 escapers). The 36 escapers form **12 classes, 3 each**; 6 converse pairs, none self-converse; H = 48,997, 49,215, 49,935, 50,221, 51,451, 52,303; all |Aut| = 1, lambda = 3, twin distance 1, module defect 1 (`aud_classes.out`) |
| THM-4942 (3) exhaustive through 11 | the logic: lambda <= 2 means SUM or SUM2; lambda >= 3 means FC at orders <= 10; at order 11, FC or one of the 7 SUB classes | correct. The order-11 list rests on S17's `lam11.c` and on S17's independent audit (`aud_o11`: all 6,538,516 order-10 classes with minimum degree >= 2, the same 97 pairs and 11 classes). I did not regenerate the 9.7M order-10 classes. So **order 12 is the first failure**: confirmed |
| THM-4985 (2) | all 193 classes in `escapers12_all_classes.txt` | all escapers, all lambda = 3, pairwise non-isomorphic. They contain the 12 classes of (1); 17 are self-converse and 99 have their converse in the list. H runs over 32,897..58,017 with 150 distinct values, 43 of them doubled. The triples are (3,1,1): 190, (3,2,2): 1 (H = 57,643), (3,3,3): 2 (H = 48,537 and 48,849) |
| decomposability defect | brute force over reversal sets (`aud_dec_codes.out`) | H = 57,643: no single reversal decomposes it, a pair does. H = 48,537 and 48,849: **no set of <= 2 reversals** gives a non-strong or non-prime tournament, a set of 3 does. Both are self-converse with |Aut| = 3 and scores 3^6 8^6 |
| creator depth (`aud_depth.out`) | the H of all 66 order-10 and all 220 order-9 subtournaments, against S(9) and S(8) | 192 classes have depth 2 (4..27 first-creator 10-subtournaments each). H = 33,021 has none at order 10 and **54** at order 9, so depth 3 |
| new structural census (`aud_struct.out`) | own code: 3-back-arc cuts; maximal modules from the pair-closure relation, a different algorithm from `subfc.py` | **106 / 47 / 13 / 53** reproduced. All FC-FC cuts are 6+6, with the same H-pair frequencies as `fcpairs.out` ((43,43) in 15 classes). The defect-3 classes are two first creators of 45, with the 3 back arcs a matching. The depth-3 class is two first creators of 33. The doubled-vertex count is 0 (see m2) |
| bands (`aud_bands.out`, `aud_misc.out`) | own DP | H(band_3(13)) = 237,817, with deletions 20,485..72,073 (7 values), all in S(11) by the census. H(band_4(16)) = 64,053,917, with 8 deletion values. All 10 witnesses recomputed: order 12 (237,817), eight of order 14, order 15 (64,053,917). Both bands are prime, with lambda 3 and 4 and twin distance 1, so **both are escapers**. band_4(15): H = 15,927,743, lambda = 4, prime, certified in S(14), so **k*(15) >= 4**. H(band_4(14)) = 4,013,709 |
| Lemmas B1, B2, W (`aud_bandlemmas.out`) | B2: max-flow for s <= 6, 2s+1 <= n <= 24, and full cut enumeration for n <= 20. B1: closure for n <= 24 and full subset enumeration for n <= 18. {0, 2s+1} is a module at n = 2s+2. W: 22,075 instances; at n = a+b-1 all 398 pairs fail, so the threshold is sharp. B2's construction: a coprime bijection found by search, the D_r inside the band, pairwise arc-disjoint and strong, for s <= 8 and 3s <= n <= 40 | 0 failures |
| Lemma B3 (`aud_b3.out`) | every Hamiltonian path of 300 random tournaments with back-span <= s (n <= 9): run count, head window, injectivity of the encoding | 0 failures. The bound s+1 on the number of runs is attained |
| shares (`aud_spectra.out`) | census T-lines | 0.5000, 0.2500, 0.5000, 0.5893, 0.7193, 0.7932, 0.8545, 0.8884, **0.9200** (831,487,872 / 903,753,248). Solid prefixes of S(10) and S(11): 14,683 and 80,451 |
| R_H, R_tau, THM-589 (`aud_ratios.out`) | exact arithmetic | R_H(5..12) as stated, **not monotone** (0.571 -> 0.600, 0.466 -> 0.474, 0.404 -> 0.407). R_tau formula and its >= 2^(m-1)/(e·m) bound hold for m <= 59. W(n) for n <= 8 equals THM-589. CV^2 = W(n)/n! - 1 by brute force for n <= 5. n·CV^2 = 1.9998 at n = 160 |
| C3[R5, R5, TT_2] | own construction | H = **229,725**, lambda = 4 (scores 4..7). C3[R5, R5, 1] = 11,475. minH witnesses re-checked: f_3 values 621, 2,055, 5,753 and 15,213 at orders 8-11; 232,001 at order 12 and 926,401 at order 13, each with minimum semidegree 4 |
| tau table (`aud_tau_le8.out`; all 268M labelled graphs at m = 8, Bareiss, checked against brute-force tree counts) | own code | every entry for m = 4..8 matches. m = 8: 251,548,592 connected graphs, 3,700 values, 3,314 new, labelled share 0.5719. 2 never occurs. **Lemma T1** holds on all 1,026,505 labelled 2-connected graphs with m = 3..7, with equality only for cycles. The vertex-insertion row holds on 138,381 (G, non-cut v) pairs. alpha(N) agrees with OEIS A182258 for every N <= 60 (all N first attained at m <= 8; the rest have alpha >= 9 there too), and alpha(10), alpha(13), alpha(22) >= 9 |
| born at the Busch bound (`aud_busch_born.out`) | census plus Busch's f | {3, 5, 9, 11, 13, 17, 19, 23, 35, 39, 63} |
| pre-registration | read-only git | the predictions (lines 1-43) are byte-identical in fb31d2c38a and now; only the Outcomes header and body changed. fb31d2c38a reached origin/main at 01:06:14; `sa13` and `tau_share` were built at 01:07:43 and 01:08:08. Both sa13 logs end "done, found 0". The Q3 numbers reproduce. **All three outcomes are stated honestly** |

### Proofs checked line by line

- Lemma W.
- Lemma B1: the distinguisher table and Claims 1 and 2, including the d = 2s + 1 case.
- Lemma B2: degrees, THM-4811 (B) for n <= 4s + 1 (whose proof I re-derived), and the Pomerance-Selfridge range n >= 3s.
- Lemma B3: run heads lie in [v - s, v - 1], so there are at most s + 1 runs; there are s + 3 choices per side; the
  encoding is injective.
- The (S-inf) theorem: g can be taken nondecreasing; the insertion lemma, oddness and the 7/21 exclusion give
  H(U) in S(|U| - 1).
- The (M_eps) theorem: by Chebyshev, P(not FC) <= P(H <= (1 - eps) E H) <= CV^2/eps^2.
- The down-universality corollary.
- The module-defect formula. For fixed M only the arcs between M and V \ M matter, independently for each outside x, and
  any tournament with module M' needs at least cost(M') reversals.
- Lemmas T1 and T3, and Theorem T2.

### Literature checked (web)

- OEIS A038375: DATA ends at a(11) = 95,095 (G. Royle, 2022).
- OEIS A182290: |S_tau(n)|, 1, 1, 2, 5, 16, 65, 386, 3700, ... (J. Azarija).
- OEIS A182258: alpha(n). It records the fixed points 3, 4, 5, 6, 7, 10, 13, 22 and the o(log n) conjecture.
- Azarija-Škrekovski, Math. Bohem. 138 (2013) 121-131, arXiv:1210.6335: the (n + 4)/3, (n + 7)/3 iff statement, as
  cited. Their paper credits the fixed points to Nebeský and completes the case 22.
- Stong, Australas. J. Combin. 82(2) (2022) 182-196: O((log n)^(3/2)/log log n). Stong attributes the o(log n)
  conjecture to Azarija and Škrekovski.
- Chan-Pak, arXiv:2308.10214: §6.5 is "Hamiltonian paths in tournaments", and Conjecture 6.2 there is attributed to bof
  and G. Royle (MathOverflow, 2016).
- Pomerance-Selfridge, Mathematika 27 (1980), as cited. The case {1..N} -> {N+1..2N} is due to Daykin and Baines.
- Alon (1990) n^(3/2); Friedgut-Kahn's improvement (for Hamiltonian cycles; it transfers to paths); Adler-Alon-Ross
  (e - o(1)). No O(n!/2^(n-1)) upper bound is known, so (M_eps) is open, as the session says.

## Findings

### MAJOR

**M1. "max H(12) = 545,835 (OEIS A038375)" is unsourced. R_H(13) is typed FINITE-EXACT on that basis, and max H(13) is
misreported.**

The facts:
- OEIS A038375 lists only a(1..11).
- The value 545,835 entered the repo unsourced in S17's HYP-9532 (E), line 70.
- The repo's own THM-329 (opus-2026-05-27-S6) records a(12) >= 531,205 and a(13) >= 3,719,831, found by search and
  "not confirmed". I recomputed H = 3,719,831 for THM-329's 13-vertex tournament. It beats the Z_13 circulant
  (3,711,175, also recomputed).
- A rigorous bound is available: mapping a Hamiltonian path of T to (its last vertex v, the rest as a path of T - v) is
  injective. So max H(12) <= 12 · 95,095 = 1,141,140 < E[H(13)] = 1,520,268.75.
- R_H(14) < 1 is not proved.

Locations and exact fixes:
- **PROOF.md §5**, the two "Further" bullets and the next sentence ("max H(12) = 545,835 (OEIS A038375, ...) gives
  R_H(13) = 0.359; - the circulant lower bound max H(13) >= 3,711,175 (...) gives R_H(14) >= 0.349. So the largest
  value of order n - 1 lies below the mean of order n, and ...").

  Replace with:

  > For n = 13, 14 the maxima are not known: OEIS A038375 lists a(1..11) only. THM-329 (opus-2026-05-27-S6) found
  > tournaments with 531,205 paths (order 12) and 3,719,831 paths (order 13) by search; neither is proved maximal.
  > Mapping a Hamiltonian path to its last vertex v and a Hamiltonian path of T - v is injective, so
  > max H(m) <= m·max H(m - 1) and max H(12) <= 1,141,140. Hence 0.349 <= R_H(13) <= 0.751 (PROVED) and
  > R_H(14) >= 0.3495 (no upper bound below 1 is known). So for every n <= 13 the largest value of order n - 1 lies
  > below the mean of order n, and every n-tournament with H above max H(n - 1) is a first creator.
- **PROOF.md §3**, last paragraph: "which exceeds the circulant value 3,711,175 that is conjectured to be max H(13)".
  Replace with: "which exceeds 3,719,831, the largest known value of order 13 (THM-329; the Z_13 circulant gives
  3,711,175)".
- **THM-4986.**
  - Title (5): "max H(n - 1)/E[H(n)] < 1 at every computed n". Add "(n <= 13)".
  - Status (5): replace "max H(12) = 545,835 as cited in HYP-9532 (E); a circulant lower bound for max H(13)" with
    "for n = 13 only 0.349 <= R_H(13) <= 0.751 (THM-329's search value; max H(12) <= 12·max H(11)); R_H(14) >= 0.3495
    from THM-329's 3,719,831".
  - Body (5): replace "R_H(13) = 0.359; and R_H(14) >= 0.349 from a lower bound on max H(13)" with
    "0.349 <= R_H(13) <= 0.751; R_H(14) >= 0.3495".
- **README** line 50. Replace "0.40 at n = 11 and 0.36 at n = 13" with "0.40 at n = 11 and 12 and at most 0.75 at
  n = 13 (max H(12) is not known exactly)".
- **README** line 148, **HYP-9565 (C)** lines 55-56 and 59, and **PREREGISTRATION.md** line 52 (appended outcomes).
  Replace "max H(12) = 545,835 (OEIS A038375 as cited in HYP-9532 (E))" and "max H(12) = 545,835" with "531,205, the
  largest known value of order 12 (THM-329; not proved maximal)".
- **PREREGISTRATION.md** line 23 is pre-registered protocol and must not be edited. Append a dated erratum to the
  Outcomes instead:

  > Erratum (audit G1): the energy proxy used 545,835 as max H(12); the value is unsourced. It affects only search
  > efficiency, because every hit had to be certified.
- **HYP-9532 (E)** line 70 (S17's file; please route to that lane). Replace "B(h) >= 25 once h > 545,835 = max H(12)"
  with "B(h) >= 25 once h >= f(25) = 390,625, in particular once h > max H(12) (>= 531,205)".

**M2. The band tournaments are THM-4880's constant-letter ribbon tournaments (priority, uncredited).**

THM-4880 (mac-mini-2026-11-04-ribbon) defines R(a) on 1..n, with back arcs q -> q - s for s in a_q and every other arc
forward. With every letter equal to {1, ..., s} this is band_s(n). THM-4880 §F already lists:
- H for the constant letter `12` (= band_2(n)), with its order-12 recurrence;
- the series `123`: 1, 1, 1, 1, 9, 41, 175, 629, 2121, 6889 (= band_3(n)).

I recomputed both series term for term: band_2(n) for n <= 15 and band_3(n) for n <= 10. The latter continues 22,143,
72,073, 237,817. Lemma B3's key step (a run head is entered by a back arc, so it lies in the window) is THM-4880 Theorem
A's fragment argument. The session credits THM-4880 only for "finite-state transfer" and introduces "band tournaments"
as a new family.

Fix:
- In PROOF §2 after the Definition, and in THM-4985's Definitions, add:

  > band_s(n) is THM-4880's ribbon tournament R(a^n) with the constant letter a = {1, ..., s}. THM-4880 §F gives
  > H(band_2(n)) (letter 12, with its recurrence) and H(band_3(n)) for n <= 10 (letter 123). Lemma B3's proof is THM-4880
  > Theorem A's window argument made explicit. New here: Lemmas B1 and B2, the escaper certificates, and the (S-inf)
  > theorem.
- Also name the family in README §2 ("Band tournaments ...") and in HYP-9565 (A).

**M3. HYP-9565 (A)'s bounded-defect question was posed by Session230 (priority, uncredited).**

Session230 committed before this session's pre-registration: 2faff63934, 2026-11-08 20:36, against fb31d2c38a, 2026-11-09
01:06. Its text:
- `session230_creator_grammar_20261109/primes/REPORT.md` §5: "OPEN: Whether all prime noncreators admit a bounded-cut
  or bounded-defect derivation from creator parents at arbitrary order; whether the finite five-operation grammar
  persists".
- SYNTHESIS §6: "a bounded-defect operation that reaches them".

The README ("Is there a bounded 'defect' operation at all ...? That is HYP-9565 (A)") and HYP-9565 (A) do not credit
it.

Fix: in HYP-9565 (A), in the README §2 last paragraph, and in THM-4985's KNOWN list, add:

> The bounded-defect question is Session230's (primes/REPORT.md §5, SYNTHESIS §6). (A) makes it precise: escapers are
> the prime non-creators that THM-4942's grammar misses, and dec measures the defect.

### MINOR

- **m1 (PROOF §0, `defect.py` docstring).** dec(T) = min(lambda, module defect) is correct but the min is redundant.
  For n >= 3 every non-strong tournament is non-prime (the side of size >= 2 of a dominating cut is a module), so the
  module defect is always <= lambda and dec = module defect. All 193 classes and every random test satisfy this. Add:
  "(for n >= 3 the module defect never exceeds lambda, so dec equals the module defect)".
- **m2 (PROOF §1.2, last bullet; `twinquot.py`).** "One tested description never applies ... (0 classes)" is forced
  by the definition of an escaper. If reversing the arc between the unique distinguisher y and u makes {u, w} a module,
  the quotient is T - u (or T - w), which is never a first creator for an escaper. Replace with: "This description
  cannot apply to an escaper: its quotient is a vertex deletion of T, which is not a first creator."
- **m3 (PROOF §1.2 table; README).** The 106 re-connected sums are exactly the 106 classes with scores 3^6 8^6
  (verified). Their halves have minimum out-degree >= 2 (resp. minimum in-degree >= 2). Every 6-tournament with minimum
  out-degree >= 2 has H in {33, 37, 41, 43, 45}, except the source plus R5 with H = 15 (`aud_halves6.out`, which also checks the 106 = 3^6 8^6 identity). So "both
  sides are first creators" and the value set {33, ..., 45} are nearly automatic from the score sequence. Say so, so
  that the census is not read as a first-creator mechanism.
- **m4 (README §2).** "So at order 12 the creative operations are 'ordinal sum, then re-connect' and 'substitute, then
  add a vertex' ..." omits the 53 classes covered by neither. Add: "... and, for the remaining 53, 'a first creator of
  order 10 plus two vertices'".
- **m5 (HYP-9565 (C)).**
  - Its first evidence bullet is about order-13 extensions. Their deletions need not be 4-arc-strong, so the bullet
    says nothing about k*(12). Move it under Q1.
  - The status line's "annealing upper bound f_4(12) <= 229,725" mislabels the source. 229,725 comes from the explicit
    C3[R5, R5, TT_2]; annealing found only 232,001. Write "f_4(12) <= 229,725 (C3[R5, R5, TT_2]); annealing missed it".
  - "This is the lambda = 4 analogue of S17's lambda = 3 floor" should not suggest that f_4(12) = 229,725.
- **m6 (HYP-9565 (A), heuristic).** "So under (S-inf) these are non-creators" is false for bare (S-inf):
  H(A)·H(B) ~ (m!/2^m)^2 is superexponential, and (S-inf) only asks g(m)^(1/m) -> inf. Write "under (S-inf) with
  g(m) = max H(m - 1), e.g. HYP-9481 (a)".
- **m7 (HYP-9565 (B)).** The class-share version follows from (M_eps) only together with the classical fact that almost
  all tournaments are asymmetric (then the class and labelled shares have the same limit). Add this.
- **m8 (PROOF §5).** "(M_eps) asks for an upper bound o(m·m!/2^(m - 1))" is imprecise. (M_eps) is weaker: it allows
  c·m·m!/2^(m - 1) with c < 1/2, and is implied by o(·). Write "(M_eps) follows from max H(m) = o(m·m!/2^(m - 1))".
- **m9.** R_H(14) >= 0.349 from 3,711,175 is false as a lower bound (0.34873). It holds with THM-329's 3,719,831
  (0.34955); absorbed by M1.
- **m10 (THM-4986 title (2); PROOF §6 Theorem T2).** "The spanning-tree counts of simple graphs are exactly
  {1} ∪ {N >= 3}" should read "of connected simple graphs" (disconnected graphs give 0). The body of THM-4986 (2) is
  already correct.
- **m11 (Lemma T3).**
  - Neither conclusion uses the nondecreasing-floor hypothesis or "1 exactly on trivial blocks". "Counts = monoid
    generated by Irr" needs only multiplicativity and realizability. "v in Irr iff v in Irr(m) for some m with
    F(m) <= v" is the definition of F.
  - The working content is: if F is nondecreasing and unbounded, that set of m is the finite initial segment
    m <= B(v). Restate the lemma that way.
  - In the PROOF §6 instance table add "21 not in Irr(<= 6)" for 21, as THM-4986 (3) already does. Being in
    (f(6), f(7)) alone does not exclude 21.
- **m12 (README, "Two differences").** "the strong bands [f(m), max H(m)] leave crevasses" is wrong for 21. The strong
  bands overlap from order 5 on, and 21 lies inside the order-6 band [15, 45]. It is a hole of the strong order-6
  spectrum. Write: "the floors leave crevasses (f(m), f(m + 1)) that only the strong spectra of orders <= m can fill; 7
  and 21 are missed by them".
- **m13 (README dictionary, insertion row).** Singleton strong components can be internal, so "equality iff v is a
  singleton component (padding by a source or sink)" should read "(e.g. padding by a source or sink)".
- **m14 (new dictionary row; THM-4986 KNOWN list).** The fixed points 3, 4, 5, 6, 7, 10, 13, 22 of alpha are Nebeský's
  theorem (L. Nebeský, Časopis pro pěstování matematiky 98 (1973) 95-97). Azarija-Škrekovski completed the case 22.
  The data are in OEIS A182258. Replace "(Azarija-Škrekovski)" with "(Nebeský 1973; case 22 completed by
  Azarija-Škrekovski 2013; OEIS A182258)". My enumeration of every graph on <= 8 vertices gives alpha(10), alpha(13),
  alpha(22) >= 9, which strengthens the stated ">= 8".
- **m15 (Lemma B2; THM-4985 status (4)).** The case used, a coprime bijection {1..s} -> {s+1..2s}, is Daykin and
  Baines, Mathematika 10 (1963) 132-136. Cite it with Pomerance-Selfridge. Note too that the (S-inf) theorem needs only
  large n, so it uses only the elementary sigma(r) = rs + 1 range and no citation at all.
- **m16 (dictionary "Chan-Pak" rows).** Chan and Pak state the 7/21 completeness conjecture as Conjecture 6.2 and
  attribute it to bof and G. Royle (MathOverflow, 2016). Write "the bof-Royle conjecture (Chan-Pak, Conjecture 6.2)",
  in line with THM-1370's priority note.
- **m17 (THM-4986 KNOWN list).** Add THM-4805 (opus-2026-10-26-S60), whose "crevasse dictionary between untouchables
  and forbidden path counts" (ANALOGY) is an earlier dictionary of the same kind. Also note that "crevasse" is
  THM-1370's term (also used in S15's camion_busch note).
- **m18 (PROOF §0).** Escapers are exactly the lambda >= 3 part of S17's REST category (THM-4942 PROOF §3.1: strong,
  prime, not a first creator, no first-creator deletion). Say so; it ties the order-12 data to S17's census.
- **m19 (THM-4986 title and README).**
  - "labelled shares 0.419, 0.515, 0.572 ... against 0.589 ... 0.920 for tournaments" compares labelled shares with
    class shares. Label the measures.
  - "so a typical tournament exceeds the old spectrum" holds at the computed orders (n <= 13 with M1's bound). At
    every order it is conditional on (M_eps).
- **m20 (PROOF §1.3, last bullet).** "delete an in-neighbour of a vertex of in-degree 3" should read "an in-neighbour of a
  vertex of in-degree 3, or an out-neighbour of a vertex of out-degree 3. At n = 12, lambda = 3 means minimum semidegree 3
  (THM-4811 (B))".
- **m21 (PROOF §1.1-1.2 counts).** "all 4,095 cuts" should be 4,094 nonempty proper subsets. "module defect by
  enumeration of all 4,094 vertex sets" should be the 4,082 sets with 2 <= |M| <= 11.
- **m22 (HYP-9532 status edit).** The session appended a second "(D)" clause. Merge it into the existing (D) and route
  the edit to S17's lane. The M1 fix to line 70 of the same file is also needed.
- **m23 (PROOF §2, B3).** No numerical check of the bound (s + 3)^(2n) is possible at computable sizes: it exceeds
  max H(n) there, so any check is vacuous. The proof carries B3 (my run-count and injectivity test checks its core).
  Do not present a numerical check of the bound as evidence.

### Side note (outside this session's files)

`00-navigation/TANGENTS.md` T064 (kind-pasteur S14b) and `00-navigation/CONCEPT-MAP.md` line 567 assert
max H(T) = Θ(n!/2^(n-1)). That is not known: Alon's upper bound carries n^(3/2). If it were known, (M_eps) and the
conditional theorems of PROOF §5 would be unconditional. The session's "(M_eps) is OPEN" is correct. Those two lines
should be corrected.

## Priority hits

I grepped 01-canon, 05-knowledge and 00-navigation. A delegated sweep was re-checked by hand for every hit below.

| claim | earlier occurrence | status in the session |
| --- | --- | --- |
| band tournaments | THM-4880 ribbon words with constant letter {1..s}; H-series for s = 2, 3 (2026-11-04) | **not credited (M2)** |
| bounded-defect question | Session230 primes/REPORT.md §5, SYNTHESIS §6 (2026-11-08) | **not credited (M3)** |
| escaper notion, grammar persistence | HYP-9532 (F) (S17); S17's REST category | (F) credited; REST not mentioned (m18) |
| first-creator shares at orders 7-10 | HYP-9532 (C) (S17) | credited; order 11 is new |
| k*(n) | HYP-9532 (D) (S17); Session230 audit REVIEW: "if f_k(n) > max S(n - 1) ..." | credited |
| crevasse | THM-1370; S15 camion_busch note; THM-4866; THM-4805's crevasse dictionary | THM-1370 cited; THM-4805 not cited (m17) |
| dictionaries of the 7/21 objects | THM-4942 §5 (Moon <-> non-cut vertices and ears); THM-4940 note §5 | credited |
| spanning-tree twin, Sedláček, Azarija-Škrekovski, Stong, alpha | no earlier repo occurrence | new to the repo; literature correctly cited, except the fixed points (m14) |
| module defect, twin distance, decomposability defect | no earlier repo occurrence. The distinguisher-closure module test is an older repo method (session24 arithmetic, session189 audits) | new; the method is standard |
| Pomerance-Selfridge, arc-disjoint spanning strong subdigraphs | no earlier repo occurrence | new to the repo; add Daykin-Baines (m15) |
| A038375 values | THM-329 (a(12) >= 531,205; a(13) >= 3,719,831); THM-128 and THM-338 (circulant 3,711,175) | THM-329 not used; 545,835 unsourced (M1) |
| attachment form, witness search | HYP-2879 / THM-4097 / THM-4114 (via THM-4942 (1)) | credited |
| canonical forms (`canon.c`) | standard individualization-refinement | not an issue. I did not rely on `canon.c` |

## Typing review

- Correctly typed:
  - FINITE-EXACT: THM-4985 (1), (2) (examples; the search is EXPLORATORY), (3) and the structural census.
  - PROVED: (4).
  - PROVED-CONDITIONAL: (5) and (6).
  - THM-4986 (1) PROVED (with possibly-folklore disclaimed); (2) KNOWN; (3) PROVED; (4) ANALOGY.
  - Every typing statement in HYP-9565.
- Overclaims:
  - R_H(13) as FINITE-EXACT (M1).
  - "So the largest value of order n - 1 lies below the mean of order n" for n = 14 (M1).
  - The circulant 3,711,175 "conjectured to be max H(13)" (M1).
  - The "0 classes" doubled-vertex result presented as a finding (m2).
  - "a typical tournament exceeds the old spectrum" without an order range (m19).
- Pre-registered Q1-Q3: correctly reported as REFUTED, with the predictions unedited. Q1 is honestly qualified as a
  failed search, not a non-existence proof.
