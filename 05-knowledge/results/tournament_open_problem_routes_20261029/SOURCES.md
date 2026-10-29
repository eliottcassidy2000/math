# Primary-source and status ledger

Checked2026-10-29. A publisher/author theorem is a CITED input; a preprint
claim is labelled as such. Search-engine dates were not used as publication
dates. Search snippets or problem aggregators supplied leads only, except
where the complete displayed primary theorem was available in indexed text
and the access limitation is explicitly noted.

| Source | What was checked and how it is used |
|---|---|
| D.R. Heath-Brown, *A mean value estimate for real character sums*, Acta Arith.72(1995),235–275, [DOI](https://doi.org/10.4064/aa-72-3-235-275), [Oxford manuscript](https://ora.ox.ac.uk/objects/uuid%3Ab188b365-d8d1-4267-8548-49f01e6cb6a7/files/m6af40476f1b860392271a2d8ecae593b) | Theorem1, printed p.237 / manuscript p.3, including both squarefree restrictions and arbitrary complex coefficients. Full statement available in indexed primary text; direct PDF fetch403. The cited theorem is the analytic input to P1, not reproved here. |
| Soundararajan–Young, *The second moment of quadratic twists of modular L-functions*, JEMS12(2010),1097–1116, [primary PDF](https://ems.press/content/serial-article-files/31734) | §2.3 independently quotes the same quadratic large-sieve inequality. No GRH-dependent L-function result is imported. |
| Tao, [prime number theorem in arithmetic progressions](https://terrytao.wordpress.com/2009/09/24/the-prime-number-theorem-in-arithmetic-progressions-and-dueling-conspiracies/) | Author exposition of the classical fixed-modulus theorem; used only at modulus4 to get N(X)~X/(2logX). |
| Nguyen–Scott–Seymour, *Induced subgraph density. VI. Bounded VC-dimension*, revised July18,2025, [author PDF](https://people.maths.ox.ac.uk/scott/Papers/ehvc.pdf) | §6, especially Lemma6.2/Theorem6.3: every two-colourable forbidden tournament has EH; all six-vertex tournaments follow. Earlier references to the last six-vertex case as open are stale. Their neighbourhood convention is stated explicitly; our P2 uses open OUT-neighbourhoods. |
| Chan–Pak, *Computational complexity of counting coincidences*, TCS1015(2024),114776, [paper §6.5](https://arxiv.org/html/2308.10214v2#S6.SS5) | Conjecture6.2 is exactly the Hamilton-count spectrum question; attribution to bof/Royle and finite-prefix evidence are already in the paper. No literature novelty for posing the problem. |
| Grinberg–Stanley, *The Rédei–Berge symmetric function of a directed graph*, [2307.05569](https://arxiv.org/abs/2307.05569) | Positive odd-power-sum formula and parity refinement are established inputs, not fresh open targets. The repository's odd-cycle specialization supplies P3's elementary H>=1+2c3 inequality. |
| Kühn–Mycroft–Osthus, *A proof of Sumner's universal tournament conjecture for large tournaments*, [author PDF](https://web.mat.bham.ac.uk/D.Osthus/exacttrees3.pdf) | The qualifier “large” matters: do not infer the full finite statement. |
| Benford–Montgomery, *Trees with many leaves in tournaments*, [2207.06384](https://arxiv.org/abs/2207.06384), JCTB170(2025),260–334 | For fixed alpha>0, large n, the leaf-sensitive host bound is(1+alpha)n+k; a separate result handles maximum degree<=cn in(1+alpha)n hosts. These are asymptotic results, not the sharp all-size tree statement. |
| Bang-Jensen–Bessy–Thomassé, *Disjoint3-cycles in tournaments*, [author PDF](https://perso.ens-lyon.fr/stephan.thomasse/liste/BTtournament.pdf), JGT75(2014) | Minimum outdegree2k-1 forces k disjoint triangles. The tournament Bermond–Thomassen case is solved. |
| Ma–West–Yan, *Lichiardopol's conjecture on disjoint cycles in tournaments*, EJC27(2)(2020),P2.52, [published PDF](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v27i2p52/pdf) | The minimum-outdegree(q-1)k-1 equal-length cycle conjecture is solved for all q>=3. It is not a new target for our ports. The competing [1706.01699](https://arxiv.org/abs/1706.01699) is explicitly withdrawn for a proof bug; search abstracts alone would mislead. |
| Cooper–Grzesik–Kabela–Král', *Packing and covering directed triangles asymptotically*, [1909.07120v2](https://arxiv.org/abs/1909.07120) | The directed factor2 analogue is already proved; their theorem gives1.8t+o(n²), with constructions near1.5t. Sharp variants require their own current status and graph-class specification. This is not undirected Tuza. |
| Chan–Grzesik–Král'–Noel, *Cycles of length three and four in tournaments*, [1902.00572](https://arxiv.org/abs/1902.00572); Ma–Tang, *Minimizing cycles in tournaments and normalized q-norms*, [2011.14142](https://arxiv.org/abs/2011.14142), Combinatorial Theory2(3)(2022),6 | Partial cycle-density minimization and the distinct2 mod4 behavior. Do not mix their cycle-density normalizations with simple-cycle counts in a fixed finite graph. No later full resolution found in this pass. |
| Fallon, *Equiangular tight frames in real symplectic space: Zauner's conjecture and the skew Hadamard conjecture*, [2509.14463](https://arxiv.org/abs/2509.14463) | PREPRINT: an explicit equivalence, involving tournament diamonds, provides a literal cross-domain route. We do not independently audit the whole paper or identify its conjecture with ordinary complex Zauner. |
| Cati–Pasechnik, *A database of constructions of Hadamard matrices*, [2411.18897](https://arxiv.org/abs/2411.18897) | Constructive database and skew-Hadamard existence scope; a specific claimed order292 needed repair. Paley constructions alone can miss whole arithmetic order families. |
| Chen–Lin–Ma–Wei, *Undecidability of polynomial inequalities in tournaments*, IMRN2025(11),rnaf131, [author/institution record](https://math.duke.edu/publications/undecidability-polynomial-inequalities-tournaments), [2412.04972](https://arxiv.org/abs/2412.04972) | Published undecidability for arbitrary polynomial homomorphism-density inequalities. It does not obstruct a fixed decidable fragment or finite-size exact enumeration. |
| Shapira–Yuster, *On Ramsey Properties of k-Majority Tournaments*, [2603.04174](https://arxiv.org/abs/2603.04174) | PREPRINT: improved n^{Omega(1/k)} transitive-subset bound for2k-1 majority orders. A specialized alternative to generic tournament Ramsey enumeration. |
| Aubian–Coulomb, *Clique Number of Tournaments II*, [2609.07481](https://arxiv.org/abs/2609.07481) | PREPRINT: multiple directed clique/chromatic conjectures proved or refuted. “Gyárfás–Sumner for tournaments” is not one stable interchangeable statement; do not import old open lists unchanged. |
| Letzter–Pokrovskiy–Williams, *On the gracesize of trees*, [2511.11331](https://arxiv.org/abs/2511.11331) | Recent primary work on an asymptotic version, not a full graceful-label existence proof. Repository shape-changing exchange theorems have an additional limitation. |

Additional status boundaries:

- Classical tournament deck reconstruction is false. Session39's
  [reconciliation](../session39_signed_phases_and_prime_cylinders_20261023/RECONSTRUCTION_PRECEDENTS.md)
  pins Stockmeyer/Kocay/McKay and distinguishes the observer languages. This
  pass inherits that audit rather than rereading all those papers.
- A September23,2026 [Seymour proof manuscript](https://github.com/openai/math/blob/main/preprints/A-proof-of-Seymours-second-neighborhood-conjecture-September-23-2026/paper.pdf)
  exists in the primary repository. The landing page was inspected; the proof
  was not audited here. Thus Seymour is a **status-review task**, not presented
  here as either a newly solved problem or an unqualified open target.
- The classical tournament case of Seymour was already a theorem. The
  Caccetta–Häggkvist Cayley/vertex-transitive lane in the May reflection also
  cites an already proved result. Partial digraphs need a missing-edge carrier;
  completion to a tournament can invent the cycle one is trying to prove.
- Equal cardinalities15,21,81 or matching tree pictures are not entries in
  this ledger unless a source, map and preserved predicate have been supplied.

Literature absence is not established by search. “No resolution found” is
a dated research assessment, not a theorem that no resolution exists.
