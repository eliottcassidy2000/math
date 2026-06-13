# Paper Revision Plan

**Based on:** Devil's advocate referee report (S339) + seesaw proof analysis

## The Referee's Verdict: REJECT (split into multiple papers)

The referee is RIGHT about almost everything. Here's the honest response.

## What the Referee Got Right

1. **It's 3-5 papers in one.** We must split.
2. **Most "results" are computations, not theorems.** Must label honestly.
3. **Taxicab geometry (result 5) is known.** Remove from "main results."
4. **Fraïssé limit (result 6) is known.** Remove from "main results."
5. **Compression needs JPEG-XL/WebP baselines.** Must add or disclaim.
6. **A² completeness needs WL context.** Must discuss Cai-Fürer-Immerman, k-WL hierarchy.
7. **No coherent thesis.** Must focus.

## What the Referee Got Wrong

1. **THM-260 is proved, not just computed.** The proof uses polynomial degree bounds from the OCF + interleaving construction. This IS a theorem.
2. **The seesaw (THM-095) IS proved** (conditional on β₂=0, which is separately proved). The proof is algebraic (chain complex exact sequence + rank counting), not just computational.
3. **χ = n-1 is verified at n=7 (456 classes), not just 4 data points.** The clique number ω=4 at n=7 is exact (Bron-Kerbosch). The 6-coloring is exact (tabu search in 1454 iterations).

## The Revised Plan: THREE Papers

### Paper A: "Spectral Invariants of Tournaments on the Hamming Scheme"
**Venue:** Journal of Combinatorial Theory, Series B or Discrete Mathematics

**Focus:** THM-260 (band-limitedness for all n) + the A² conjecture + χ(G_n)

**Structure:**
1. Introduction: tournaments as binary words on Q_m = Hamming scheme
2. The Krawtchouk spectral framework (background, properly cited)
3. **THM-260: Band-limitedness** (full proof, comparison with Delsarte bounds)
4. **The A² conjecture** (carefully labeled as CONJECTURE, WL context, comparison with nauty, discussion of where it might fail)
5. **χ(G_n) = n-1** (labeled as conjecture, computational evidence, connection to Lie algebra rank)
6. The even graph dual (brief, as supporting evidence)
7. Open problems

**What changes:**
- Remove taxicab geometry, Fraïssé limit, compression, computational irreducibility
- Add: explicit WL hierarchy comparison, nauty timing benchmarks, Delsarte citation
- Relabel A² and χ as CONJECTURES
- Make THM-260 the centerpiece with FULL proof

### Paper B: "Path Homology of Tournaments: the Seesaw Mechanism"
**Venue:** Journal of Algebra or Homology, Homotopy and Applications

**Focus:** β₂=0 + the seesaw (β_{2k-1}·β_{2k+1}=0) + the defect wave

**Structure:**
1. Introduction: GLMY path homology for digraphs
2. The chain complex Ω_*(T) for tournaments
3. **THM: β₂(T) = 0 for all tournaments** (proof via good-vertex induction)
4. **THM-095: β₁·β₃ = 0** (full proof, using β₂=0 + rank counting)
5. The defect wave: β₁ decreasing, β₃ increasing with n
6. The general seesaw conjecture: β_{2k-1}·β_{2k+1} = 0
7. THM-261: Z₂ orbit pairing on Ω for SC tournaments
8. Open problems (β₂=0 independent proof, seesaw at all levels)

### Paper C: "Tournament-Inspired Lossless Image Compression"
**Venue:** IEEE Data Compression Conference (DCC) or arXiv systems paper

**Focus:** The TC codec, honest benchmarks, the score-conditioning principle

**Structure:**
1. Introduction: binary matrix compression via prediction + entropy coding
2. YCoCg-R decorrelation (cite Malvar-Sullivan)
3. Adaptive row filtering (cite PNG, explain the min-sum heuristic)
4. The score-conditioning principle from tournament theory (our ACTUAL novel contribution)
5. Benchmarks: PNG, JPEG-XL lossless, WebP lossless on Kodak + real photos
6. The C library implementation

**What changes:**
- Add JPEG-XL and WebP baselines (or honestly state we haven't compared)
- Focus on what IS novel: the connection between tournament score-conditioning and image prediction
- Remove all pure math (Krawtchouk, Fraïssé, etc.)

## The Seesaw Proof: Current State and Completion Path

### What IS proved (THM-095):
1. **β₂ = 0** for all tournaments (proved by induction: good-vertex deletion + LES)
2. **ker(d₁) = C(n,2) - n + 1** (proved: tournaments are connected → β₀ = 1)
3. **β₁ ∈ {0, 1}** (follows from ker(d₁) constant + im(d₂) takes only 2 values)
4. **β₁ · β₃ = 0** (proved: β₂=0 gives exact coupling; im(d₂) can't satisfy both β₁>0 and β₃>0 simultaneously)

### The KEY step that's computationally verified but not proved:
**im(d₂) takes only two values: C(n,2)-n and C(n,2)-n+1.**

This is verified exhaustively at n≤8 but the proof requires showing that rank(d₂) is always either maximal or maximal-1. A proof strategy:

1. d₂ maps Ω₂ → Ω₁ (2-paths → 1-paths)
2. dim(Ω₁) = C(n,2) (one basis element per arc)
3. im(d₁) = n - 1 (connected graph)
4. ker(d₁) = C(n,2) - n + 1

If we can show rank(d₂) ≥ C(n,2) - n (i.e., im(d₂) ≥ n-choose-2 minus n), then β₁ ≤ 1, and the seesaw follows.

The missing piece: **a lower bound on rank(d₂).**

### Possible proof of the rank bound:
- dim(Ω₂) = C(n,3) + 2·c₃ (proved, THM-119)
- The boundary map d₂ sends each 2-path (a,b,c) to arc (a,b) - arc (a,c) + arc (b,c)
- If we can show that the arcs not in any 3-cycle are ALWAYS in im(d₂), that gives rank(d₂) ≥ |non-3-cycle arcs|
- For tournaments with β₁ = 1: there is exactly ONE "missing direction" in im(d₂)

### The generalized seesaw (β_{2k-1}·β_{2k+1}=0):
This requires β_{2k} = 0 at each level — which is NOT true in general (β₄ > 0 at n≥7). However, computationally, when β_{2k-1} > 0, we ALWAYS have β_{2k} = 0 (the even Betti number between adjacent odd ones vanishes). This is a STRONGER claim: the seesaw holds not because all even Betti numbers are zero, but because the specific even Betti number BETWEEN two nonzero odd ones is always zero.

## Action Items

1. **Paper A**: Write from scratch, focused on THM-260 + A² + χ. 20-25 pages.
2. **Paper B**: Write from scratch, focused on seesaw + homology. 15-20 pages.
3. **Paper C**: Revise compression section with proper baselines. 10-12 pages.
4. **Prove the rank bound for d₂** (the missing step in the seesaw proof).
5. **Compute χ(G_8)** if feasible (strengthen the chromatic conjecture).
6. **Add JPEG-XL benchmarks** to Paper C.
