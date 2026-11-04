# Independent adversarial audit: THM-4870 (the Thomassen merge), HYP-9462 and the merge note

**Auditor:** independent audit subagent of kind-pasteur-2026-11-04-S10, 2026-11-04.

**Audited:**
- `01-canon/theorems/THM-4870-thomassen-merge-reducible-or-obstruction-dichotomy-asymmetric-path-spectrum-3-7-21-square-sum-vertex-deletions.md`;
- `05-knowledge/results/thomassen_merge_hypohamiltonian_groups_20261104.md` (the note);
- `05-knowledge/hypotheses/HYP-9462-asymmetric-path-spectrum-3-7-21-and-noncritical-vertices-of-uniquely-hamiltonian-graphs.md`.

**Method.** All code in this folder was written for the audit. The author's programs were read only to learn what
they check; none was compiled or reused, and the author's outputs are compared against, never used as input.
- The tournament classes come from my own canonical labelling by individualisation-refinement. Every |Aut| is
  computed twice: as the IR leaf count and by a separate backtracking enumerator.
- The uniquely Hamiltonian (UH) graph census uses two independent methods (brute force and down-set).
- Single objects and certificates are re-checked in Python, using DFS path counts and networkx VF2, which share no code
  with the C programs.
- Compute: about 9 minutes of wall time in total, the largest job being the 5 min 23 s order-10 census below.

## Verdict: SOUND WITH CORRECTIONS

There is no blocking finding.
- Every finite claim in (A)-(D) reproduces exactly with independent code.
- The headline result (B), "asymmetric path spectrum = odd numbers minus {3, 7, 21} for H <= 1125", is confirmed by a
  stronger construction than the author's: explicit asymmetric tournaments of order <= 9 for all 560 values.
- The author's order-<= 10 census (107,773 strong classes with H <= 1125) is reproduced exactly.
- One stated formula is wrong: the error term in (A) is Theta(n·φ^(-4n)), not O(φ^(-4n)).
- The note miscounts the non-planar Thomassen graphs (six, not five) and gets the coincidence-ledger arithmetic wrong.
- HYP-9462 misstates the first order at which a chord joins two non-critical vertices (8, not 9).
- The words "exact analogue" / "exact counterpart" over-type an interpretive analogy inside a PROVED statement.

## Findings

| # | severity | location | problem | fix |
|---:|---|---|---|---|
| 1 | correction | THM-4870 title and (A); note §1 (iii) and proof | The error term is not O(φ^(-4n)). Exactly, E[r] − ((5+√5)n + 12 − 8√5)/10 = (72 + 161√5/5)(n − 8) φ^(-4n) + O(n φ^(-8n)), with 72 + 161√5/5 ≈ 144.0014 (sympy, `gf_mean.out`; the constant term is exactly −8 times the slope). The error is exactly 0 at n = 8, where both sides equal 26/5. The scaled error err·φ^(4n) grows linearly: 144.0, 288.0, …, 10368.1 at n = 9, 10, …, 80 (`gf_mean.out`). The cause is that the double pole at φ^2 also contributes a term linear in n. Leading term, slope (5+√5)/10, constant (6 − 4√5)/5 and the 10^(-13) check from n = 20 are all correct. | Replace O(φ^(-4n)) by O(n φ^(-4n)), or state the first-order term. |
| 2 | correction | THM-4870 title ("an exact analogue of the … theorem for planar triangulations"); THM-4870 (A) ("the exact counterpart"); note §6 ("the exact tournament analogue"); commit `aaaf77d896` ("tournament form of the 2026 … four-colour theorem") | The status line types the 4CT analogy as INTERPRETATION, but the title and body of a PROVED theorem call it exact. Concretely: (i) the 4CT theorem's reducible configurations are *pairwise non-touching*, which is what allows simultaneous reduction and the O(n log n) algorithm, whereas here the r reductions form ONE run of consecutive, mutually adjacent vertices and no simultaneous-reduction statement is proved; (ii) the 4CT statement is a lower bound with an unspecified constant, while r + k = n + 1 is an identity (the note itself says so); (iii) obstructing cycles of length <= 5 are not one-arc cuts. Only the "nested / non-crossing" part matches. | Say "analogue (INTERPRETATION)" in the title and body, and drop "exact". |
| 3 | correction | note §0 bullet "{6,4,8,9,12}" and §6 first line | "the five non-planar Thomassen graphs … 20 has 6, 24 has 4, 32 has 8, 41 has 9, 60 has 12" is wrong. MathWorld (fetched 2026-11-04) says "The 32- and 34-vertex Thomassen graphs both have graph crossing number 8 (E. Weisstein, Oct. 5, 2026)", and only the 94- and 105-vertex graphs are planar. So SIX Thomassen graphs are non-planar, with crossing numbers 6, 4, 8, 8, 9, 12 (on 20, 24, 32, 34, 41, 60 vertices), and {6,4,8,9,12} is their set of distinct values in vertex order. The identification survives; the count and the omission of 34 do not. | "the distinct crossing numbers of the six non-planar Thomassen graphs (20:6, 24:4, 32:8, 34:8, 41:9, 60:12)". |
| 4 | correction | note §7 (coincidence ledger) | "The landmark set {15, …, 42} has 15 of the 31 integers in [20, 50]" is false: 15, 16 and 17 lie below 20, so only 12 do. The expected number of hits among the 5 Thomassen orders in [20, 50] is 5·12/31 = 1.94, not 2.4, and P(>= 4 hits) (hypergeometric) is 0.060, not about 0.146 (`coincidence_ledger.out`). "Unremarkable" is defensible only because both the landmark list and the window were chosen after the fact (look-elsewhere effect). | Correct the numbers and give p ≈ 0.06 with the post-hoc caveat. |
| 5 | minor | HYP-9462, exploratory evidence (B), third bullet | "Two non-critical vertices can be adjacent through a chord (first at n = 9)" is wrong. Exhaustively there are 0 such graphs for n <= 7, **16 at n = 8**, 180 at n = 9 and 7,480 at n = 10 (both of my methods; the author's own confirm.out B1 line also reports `adjacent_NC_pairs=16` at n = 8). The Results section never corrects this. | "first at n = 8 (16 labelled graphs)". |
| 6 | minor | THM-4870 (C1); note §2 (C1) | "For N <= 30 the minimum degree is 1 (THM-4505)" misquotes THM-4505, which says "minimum degree >= 2 iff n >= 31". For N <= 6 the minimum degree is 0 (vertex 2 is isolated until N = 7). Degenerate case: Q_2 = 2K_1 is non-traceable while both of its vertex-deleted subgraphs (K_1) are traceable, so it is hypotraceable under the bare definition, and the leaf argument does not reach it. Q_1 = K_1 is degenerate for both notions. | Write "at most 1", and add the usual convention (graphs on >= 3 vertices; also handle isolated vertices) or restrict to N >= 3. The conclusion then holds: verified directly for 3 <= N <= 40; for N > 40 it inherits THM-4851 (a). |
| 7 | minor | note §5, first bullet | "there are F(2n − 6) of them, from n = 3 on", with "rigid": at n = 3 the only UH tournament is C_3, with \|Aut\| = 3, and F(0) = 0 ≠ 1. Also, "no asymmetric graphs on <= 5 vertices" ignores K_1. | "from n = 4 on"; "on 2..5 vertices". |
| 8 | minor | note §3, "Converse" bullet | "\|Aut±\| ≡ 2 mod 4, with the primitive exceptions C_3, Z_5{1,2} and P_7 (THM-4743)" reads as if these three violated the congruence. They satisfy it: \|Aut±\| = 6, 10, 42 (`aut_pm_check.out`). In THM-4743 they are the primitive *Gluck-exceptional* tournaments, and the congruence is an ingredient of that classification. | Reword. |
| 9 | minor | note §0, "Coincidence in the repository" | The handoff file is described as listing "the speed prefix [4, 6, 8, 9, 12] with sixth-speed lifts x1 = 14 + 13j". In fact it lists every one of the C(12,5) = 792 five-subsets of {1..12}. [4,6,8,9,12] has 95 lifts <= 188 in several residue classes mod 13 (14, 15, 16, 18, 20, …), and 105 is a lift of 462 of the 792 prefixes (`lrc_handoff_check.out`). The verdict "coincidence" is right and understated: the co-occurrence is guaranteed by the file's exhaustive design. | Say so. |
| 10 | minor | THM-4870 (B), second bullet; HYP-9462 (A); explore.out | "Strong tournaments of order >= 11 have H >= 225" is used to make the range below 225 complete, with no proof or citation. It is true. Busch, Electron. J. Combin. 13 (2006) #N3 (Moon's conjecture) gives the exact minima 3, 5, 9, 15, 25, 45, 75, 125, 225, … My census also gives an independent proof: no strong class of order 11 has H <= 224, and Moon 1966 together with H(T) >= H(T − v) carries this to every order >= 11. The same gap affects "UH order 15 has H > 1125" (Busch's minimum at 15 is 1875). | Cite Busch 2006, or add the two-line argument. |
| 11 | remark | THM-4870 (B), third bullet | "The only strong value <= 1125 among orders <= 10 that needs a symmetry is 3" holds only cumulatively. Restricted to orders <= 7, the values 171, 175 and 189 occur only in symmetric strong classes, and restricted to orders <= 8, the values 621 and 633 do. Asymmetric strong representatives appear one order later. The cumulative statement is verified at orders <= 9 (all H) and <= 10 (H <= 1125), see below. | None needed; optionally say "cumulatively". |
| 12 | remark | THM-4870 (B), first bullet; HYP-9462 "Proved" | The step "a strong component with H = 3 is C_3" needs one line: a strong tournament of order m >= 4 has a Hamiltonian cycle (Camion), hence at least m >= 4 Hamiltonian paths (one per deleted cycle arc), hence H >= 5 since H is odd (Rédei). | Add the line. |
| 13 | remark | note §2, "Planar and hypohamiltonian" | "The planar square-sum graphs (N <= 24) are exactly those of the leaf-ridden, non-Hamiltonian regime". The planar range N <= 24 (confirmed with networkx for N <= 40) is a proper subset of the leaf range N <= 30. | "lie inside". |
| 14 | remark | note §4, "The fives" | ADE is listed among the 2 + 3 splits that FIVE_LADDER types as cardinality coincidences. The rotation-group/ADE match is not a coincidence (McKay correspondence via the binary polyhedral groups). | Exclude that pair from the coincidence typing. |
| 15 | remark | note §0, triangulation-theorem bullet | "More than 8,200 D-reducible configurations" is not in the arXiv abstract or comments, which are the stated sources. It does appear in the paper's text (web search). Title, authors, v3 date (1 Oct 2026) and main theorem match arXiv:2603.24880. | None. |

## What was verified (all with the audit's own code)

### (A) Reducible-or-obstruction; exact mean (`gf_mean.py`, `tourn_audit.c` modes full / uh)

**From the generating function** t^2 x^4 (1 + (t^2 − t)x)/(1 − (1 + 2t)x + t x^2):
- it reproduces the published N(n, r) rows for n = 5..13;
- Σ_r N(n, r) = F(2n − 6), N(n, 2) = 1 and N(n, n − 1) = 2^(n−5), with support in [2, n − 1], for n <= 80;
- ∂_t G(x, 1) = x^4 (x^3 − 2x^2 − 3x + 2)/(1 − 3x + x^2)^2 exactly;
- the leading term of E[r] equals the claimed ((5+√5)n + 12 − 8√5)/10 symbolically;
- the error term is in finding 1.

**Own UH class census.** UH classes were taken from the complete class census for n <= 9, and by all-pattern
extension of UH classes for n = 10, 11. The extension is valid for the following reason.
- Moon 1966 gives a vertex v with T − v strong.
- Inserting v into a Hamiltonian cycle of T − v is injective, so hc(T − v) <= hc(T) = 1.
- Camion gives hc(T − v) >= 1.
- So T − v is UH, and every UH class of order n extends one of order n − 1.
- Counts 1, 1, 3, 8, 21, 55, 144, 377, 987 for n = 3..11 (Garey for n >= 4).
- Every UH class with n >= 4 is asymmetric.
- r + k = n + 1, with k the number of cycle arcs whose deletion leaves T non-strong.
- The non-critical set is one run with 2 <= r <= n − 1.
- The rows N(n, r) equal the published table.
- E[r] = 2, 3, 3.75, 94/21, 5.2, 853/144, 2506/377, 7275/987 for n = 4..11, equal to the generating-function values.

### (B) The asymmetric spectrum (`tourn_audit.c`, `strong_par.c`, `cert_search.c`, `verify_certs.py`, `analyze_B*.py`)

**Complete class census to order 9.** The class counts 1, 1, 2, 4, 12, 56, 456, 6880, 191536 match OEIS A000568, and
the strong counts 1, 0, 1, 1, 6, 35, 353, 6008, 178133 match A051337.
- \|Aut\| by IR and by backtracking agree on every class.
- \|Aut\| divides H throughout.

**Strong classes with H <= 224**, generated as strong extensions of strong classes (valid by Moon 1966 + monotonicity):
- 1, 1, 1, 6, 35, 353, 1233, 789, 225 at orders 1, 3, …, 10, and **none at order 11**;
- 2,643 excluding order 1, which matches the author's "strong classes read: 2643";
- the minima 3, 5, 9, 15, 25, 45, 75, 125 (orders 3..10) match Busch's formula.

**The spectrum below 225.**
- Holes of the full spectrum: {7, 21}. Holes of the asymmetric spectrum: **{3, 7, 21}**.
- Values needing symmetry: [3]. No strong class of order <= 10 has H = 7 or 21.
- \|Aut\| values: 1, 3, 5, 7, 9, 21.

**Certificates.** All 560 odd H <= 1125 other than 3, 7, 21 have an explicit asymmetric tournament of order **<= 9**
(`certs.txt`). Python re-counted H by DFS enumeration and confirmed a trivial automorphism group with networkx VF2:
0 failures.
- This is a stronger sufficient construction than the author's, which used strong classes of order <= 10 and UH
  tournaments up to order 14.
- With 3 forced (finding 12) and 7, 21 impossible (THM-1370; also re-derived from the census plus the order-11 bound),
  **(B) holds for H <= 1125.**

**Product logic.**
- H is multiplicative over strong components: exhaustive over all 33,867 labelled tournaments of order <= 6.
- \|Aut\| is the product of the components' \|Aut\|: exhaustive at order <= 5. The proof is that the condensation is
  transitive and therefore rigid.
- 30 random compositions of certificates (with and without C_3) behave as predicted.

**Small groups.**
- Strong order-6 classes have \|Aut\| in {1, 3} (29 and 6 classes).
- The first strong non-cyclic groups occur at order 7, exactly two:
  - \|Aut\| = 9 (Z_3 × Z_3: there is no element of order 9), H = 99, hc = 9, scores 2,2,2,3,4,4,4, self-converse;
  - \|Aut\| = 21 (F_21), H = 189, hc = 24 (P_7).
- C_3 ⇒ C_3 has H = 9.
- Aut±(C_3) ≅ D_3 and Aut±(Z_5{1,2}) ≅ D_5. Aut±(P_7) has order 42 and no element of order 21, so it is neither D_21
  nor cyclic, hence not a Maschke group.
- The note's group theory is correct: tournament groups have odd order, and the odd-order finite subgroups of O(3)
  are cyclic.

**"Needing symmetry", orders <= 9 (all H).** The values are 3 together with 2931, …, 3357 (all > 1125); within
H <= 1125 only 3 (see finding 11 for orders <= 8).
- \|Aut\| values with H <= 1125 at orders <= 9: 1, 3, 5, 7, 9, 15, 21.

**UH-product closure (secondary).** Using UH classes of orders 4..11 only, the holes are exactly the claimed 24-value
list up to 542; this is rigorous below 375, the minimum H at order 12. The claimed list is consistent: no claimed hole
is covered.

**Order-10 census, H <= 1125** (`strong_par.c`, parallel strong extension; 5 min 23 s). The counts are 1, 1, 6, 35,
353, 6008, 53275, **48094** at orders 3..10, so **107,773** in total excluding order 1. This equals the author's
"strong classes read: 107773" exactly.
- \|Aut\| values: **1, 3, 5, 7, 9, 15, 21**; \|Aut\| divides H throughout.
- Strong values attained only with symmetry, cumulatively:
  - orders <= 7: [3, 171, 175, 189];
  - orders <= 8: [3, 621, 633];
  - orders <= 9 and <= 10: **[3]**.
- The product closure of the asymmetric strong classes of order <= 10 misses exactly **{3, 7, 21}** in odd [1, 1125].
  This reproduces A1 without the UH input.

### (C) Square-sum graphs (`square_sum_audit.c`; DFS with sound pruning, cross-checked by a subset DP for <= 25 vertices)

**Degrees and traceability.**
- Minimum degree: 0 for N <= 6, 1 for 7 <= N <= 30, 2 for 31 <= N <= 40.
- Q_N is traceable exactly for N in {1, 15, 16, 17, 23} ∪ [25, 40], and Hamiltonian exactly for 32 <= N <= 40.

**(C2).** For N = 18, 19, 20, 21, 22, 24 the number of v with Q_N − v traceable is **1, 1, 1, 1, 2, 4**. The sets are
{18}, {18}, {18}, {18}, {2, 18} and {11, 18, 22, 24}, and the DP agrees.
- Q_N − v is Hamiltonian only for v = 18 at N = 23, 25 and 26, and Q_N − 18 is Hamiltonian there.
- At N = 27..30 no Q_N − v is Hamiltonian, which the note does not state.
- Q_31 − v is Hamiltonian exactly for v in **{2, 4, 13, 15, 26, 28, 29}**.

**(C3).** Q_32 and Q_33 have exactly one Hamiltonian cycle each, with non-critical sets **{2, 13, 15, 26, 29, 30}**
and **{2, 13, 15, 18, 26, 29, 30, 32, 33}**. Both sets are independent in Q_N, and neither has two cycle-consecutive
vertices.

**Planarity** (networkx): Q_N is planar exactly for N <= 24 within N <= 40, and hence for all N.

**(C1) argument.** A hypohamiltonian graph has minimum degree >= 3, a hypotraceable graph has no leaf, and
Hamiltonian implies traceable. The argument is correct apart from finding 6.

### (D) Uniquely Hamiltonian graphs on the fixed cycle (`uh_graph_audit.c`, `verify_witnesses.py`)

| n | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| UH graphs (brute force over all chord sets) | 3 | 11 | 69 | 561 | 6,903 | 113,107 | n/a |
| UH graphs (down-set enumeration) | 3 | 11 | 69 | 561 | 6,903 | 113,107 | **2,689,657** |
| two cycle-consecutive non-critical vertices | 0 | 0 | 0 | 0 | 0 | 0 | **1,300** |
| two non-critical vertices joined by a chord | 0 | 0 | 0 | 0 | **16** | 180 | 7,480 |

**n = 10.**
- The \|NC\| distributions agree with the author's.
- All 1,300 witnesses have \|NC\| = 6, and all 1,300 graphs with \|NC\| = 6 are witnesses.
- The fewest edges is 14. My lexicographically first minimum witness, {0,5}, {1,7}, {2,6}, {3,8}, is the author's
  witness rotated by one.

**Direct checks.**
- The quoted witness C_10 + {1,6}, {2,8}, {3,7}, {4,9} has exactly 1 Hamiltonian cycle, NC = {0, 2, 3, 5, 7, 8},
  and consecutive pairs (2,3) and (7,8).
- The n = 13 counterexample of HYP-9462 has 19 edges and a unique cycle 0..12, with NC = {2, 4, 5, 7}. The two
  quoted Hamiltonian cycles of G − 4 and G − 5 are valid.
- The "proved partial case" of HYP-9462 (B) (two shortcut forms give the second cycle c_(n−1) c_1 c_0 c_2 … c_(n−1))
  is correct.

### (E) Citations and typing

- **MathWorld:** see finding 3.
- **arXiv:2603.24880:** title, six authors, v1 25 Mar 2026, v2 7 May 2026, v3 1 Oct 2026. The abstract's theorem is
  "linearly many pairwise non-touching reducible configurations or pairwise non-crossing obstructing cycles" of length
  <= 5, as quoted.
- **arXiv:2204.01947** (Royle, Praeger, Glasby, Freedman, Devillers, J. Algebraic Combin. 57 (2023) 515-524): correct.
- **Maschke:** the planar Cayley groups are the finite subgroups of O(3). "Planar iff cyclic" is PROVED and correct
  (odd order excludes reflections and D_n, A_4, S_4 and A_5).
- **Arithmetic in §0:** {4,6,8,9,12} is the set of composite {2,3}-smooth numbers <= 12, and
  {ab mod 13 : a, b ∈ {2,3,4}} = {3,4,6,8,9,12}. Both are correct.
- **Coincidence typing:** correct (see findings 4 and 9).

## Reproduction (from this folder; gcc 15.2, Python 3 with sympy and networkx)

```
GCC=C:/Users/Eliott/scoop/apps/gcc/current/bin/gcc
python3 gf_mean.py > gf_mean.out                                      # (A) symbolic, 50 s
$GCC -O2 -o ta tourn_audit.c
./ta full 9 > full9.out                                               # all classes <= 9 (writes strong_classes.txt), 9 s
./ta uh 11 > uh11.out                                                 # UH classes <= 11 (writes uh_classes.txt), 5 s
./ta strongB 11 224 > strongB224.out && cp strongB_classes.txt strongB224_classes.txt   # 29 s
python3 analyze_B.py > analyze_B.out
$GCC -O2 -fopenmp -o sp strong_par.c && ./sp strongBpar 10 1125 > strongB1125.out       # order-10 census
python3 analyze_B1125.py > analyze_B1125.out
$GCC -O2 -o cs cert_search.c && ./cs 1125 200000 20261104 > certs.txt 2> certs.log      # 49 s
python3 verify_certs.py certs.txt > verify_certs.out
python3 multiplicativity_check.py > multiplicativity_check.out        # 33 s
$GCC -O2 -o sq square_sum_audit.c && ./sq 1 40 > square_sum.out       # 1.3 s
$GCC -O2 -fopenmp -o ug uh_graph_audit.c
./ug brute 4 9 > uh_graph_brute.out ; ./ug down 4 9 > uh_graph_down.out ; ./ug down 10 10 > uh_graph_down10.out
python3 verify_witnesses.py > verify_witnesses.out
python3 aut_pm_check.py > aut_pm_check.out
```

```
python3 coincidence_ledger.py > coincidence_ledger.out                # finding 4
# from the repository root:
python3 04-computation/experiments/thomassen_merge_20261104_audit_S10/lrc_handoff_check.py > .../lrc_handoff_check.out   # finding 9
```

`python3 uh_closure.py > uh_closure.out` runs the secondary UH-closure check (needs `uh_classes.txt`).

Some files were moved out of this folder after the runs. They were not deleted: they are in the session scratchpad.
- The compiled binaries (`*.exe`).
- The scratch copy `strongB_classes.txt`, identical to `strongB224_classes.txt`, which is kept.
- Two bulky regenerable class lists:
  - `strong_classes.txt` (4.6 MB), written by `./ta full 9` in 9 s and read by `analyze_B.py`;
  - `strongB1125_classes.txt` (2.4 MB), written by `./sp strongBpar 10 1125` in 5 min and read by
    `analyze_B1125.py`.

  Regenerate them with the commands above before re-running the two analysis scripts.
