# PROBLEM LEDGER — the repo's novel results across problems, and the frontier for future work

> **PORTFOLIO INDEX, NOT HEADLINE STATUS (refresh warning 2026-07-21):** Use
> [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md) for current truth. Major later
> changes include THM-2022 proving NC2/GMC(2), explicit GMC(3) counterexamples,
> THM-2044/2045 separating rank-two Poisson from DC(2)/planar JC, partial NC2
> Lean coverage, and corrections MISTAKE-211–234. Historical priority/provenance
> claims below require primary-source and canon checks.

## Current portfolio at a glance

This table is the maintained entry point. Every row names a proved foothold and
the obstruction that still matters; the long inventory below is provenance.

| Area | Current truth | Highest-value next move |
|---|---|---|
| LRC(14) | Open. THM-2051/2052 give a bounded relation and finite atlas; THM-2053/2055/2056 give the deck, disk, owner fan, and Farey certificate. THM-2057 closes two AP one-tail planes; THM-2058 gives primitive owner intervals; THM-2059 joins clock packets by CRT. | Discharge every labelled interval by clocks, pair sums, Fejer, Euler, or rank gain. THM-2060/2061 remain reserved audits. |
| NC2 / GMC | NC2/GMC(2) are proved by THM-2022. The root-imported Lean spine is sorry-free at proved nodes; separate `GMC2NC2Capstone` typechecks the conditional skeleton with one `sorry`. GMC is false from dimension three. | Use HYP-8878's unique-minimum-channel lemma where applicable; wire the concrete capstone and formalize or explicitly cite DvdK on the residual. |
| Jacobian / Dixmier / Poisson | The two-pair Poisson conjecture is false. DC(2) and planar JC remain open; filtered first-order descent is obstructed, while THM-2049 proves the associated-graded correction complex acyclic. | Decide finite polynomial termination and the coupled `D` relations, or prove boundary poles prevent a polynomial lift. |
| Tournaments | Join laws and strong-core localization are proved; the invariant lattice is exact only through `n<=6`. Positive game equilibria need skew symmetry plus mod-2 support parity. GLMY path homology in degree at least three is under an active convention dispute. | Attack `H>=disc`, classify operation laws, resolve the path-homology convention case, and track positive-kernel support rather than conflating games, flows, and topology. |
| Reciprocal sequences | Support and indexed-multiplicity Dirichlet profiles are separated by a collision tax; Abel--Stieltjes/Dini/Bertrand boundary tests are available. | Build an operation calculus for support unions/products and connect collision structure to LRC denominator spectra. |
| Arrangement carrier | THM-2047 proves the labelled phase-height carrier, top wedge, Euler detector, Fejer formula, and paired deletion; THM-2050 proves period-14 local germs are globally blind. | Find a deletion/localization invariant preserving owner, side, height, and global termination that can force the AP-core branch. |
| Cross-domain wildcard | Many analogies survive only after naming the lost coordinate. MISTAKE-226 repairs diagonal energy; MISTAKE-230--234 quarantine universal entropy, composite-Paley clocks, modular-cusp LRC claims, and support-semigroup noncancellation. | Generate alternate objects/quotients, then demand a map, preserved predicate, hostile control, and explicit sidecar before promotion. |

## Legacy portfolio snapshot

The material below is a dated research inventory. It is useful for names,
failed routes, and provenance, but its status labels do not override the table
above, `CURRENT-FRONTIER.md`, canonical theorems, or recent corrections.

**Canonical consolidated ledger.** Created death-star-2026-07-20-S59u (HYP-8185),
owner-directed: given the Jacobian Conjecture was disproved externally, inventory
our NOVEL, non-obvious results across the surrounding problems, expand the list
with under-the-radar threads, and set priorities for future creative work.

**This is a CONSOLIDATION.** The owner prompt went fleet-wide and produced several
overlapping drafts — **klein-S332** (`PROBLEM-LEDGER-2026-07-20-klein-S332.md`, a
full peer ledger, merged here with credit), **kind-pasteur-S128c104** ("PROBLEM-
ATLAS", stub), **mac-mini-S140** ("PROBLEM-PORTFOLIO", stub). This file merges them
plus four death-star-S59u audit sweeps into one living index; the scattered drafts
should be retired in favor of this one. Grades: **PROVED** (in-repo proof) ·
**PARTIAL** (verified sub-case / bounded-exact with named residual) · **REFRAMED**
(new formulation relocating the difficulty) · **CONJECTURED** (evidence, no proof) ·
**UNTOUCHED** (named, not worked = a target). LRC(≤13) settled by citation.
**Future sessions: update grades, add problems, keep it honest.**

---

## A. The Jacobian ecosystem (van den Essen's world)

### A1. Jacobian Conjecture (JC_n) — DISPROVED (n≥3), verified in-repo.
THM-1300: the owner map F:ℂ³→ℂ³ has det JF≡−2 + triple collision (six independent
in-repo verifications). JC_n false ∀n≥3 by stabilization. JC_2 survives (A9).

### A2. Dixmier Conjecture (DC_n) — DC_{n≥3} constructively FALSE (PROVED); DC_1, DC_2 OPEN.
- **THM-1300 §1** — explicit A₃ Weyl endomorphism φ(X_i)=F_i, φ(D_j)=Σ B_jk D_k,
  B=(JF^T)⁻¹ over ℤ[1/2]; 18 Weyl/flatness identities verified; non-surjectivity by
  the module one-liner ⟹ A₃ has an explicit proper self-embedding. **The strongest
  Dixmier result: a held object, not just the classical contrapositive.** Double-
  constructed (death-star-S59m + mac-mini-S130); klein's symplectic route (A3).
- **DC_1** — REFRAMED (S59r, HYP-8160): A₁'s weight triple (N,q,p) IS the oriented
  3-cycle (observer=ℏ=the conserved +1); two-lens strategy (Rédei parity × 2D
  leading-form bound). OPEN (JC_2 ⟹ DC_1 via BKK); the surviving floor of the tower.
- **DC_2** — OPEN, but no longer computationally empty. HYP-8802 proves that
  naive Weyl ordering of the THM-2044 Poisson witness fails in five of six
  relations and gives finite exact corrections of the entire `R`-column:
  `M(Dq,R)=1`, `M(T,R)=M(Sq,R)=0`. No simultaneous `A_2` endomorphism is yet
  constructed; the coupled `T`-column/flatness termination gate remains.

### A3. Poisson Conjecture — FALSE already for two canonical pairs (PROVED).
**THM-2044** gives explicit `R,T,D,S in Q[x,q,p,z]` with two canonical bracket
pairs, every cross-bracket zero, and an exact three-point fibre. It is a
polynomial symplectic suspension of THM-1300, so PC(2) is false and padding gives
PC(n) false for every `n>=2`. This improves the older six-dimensional cotangent
lift (PC(3)) and is written as a canonical theorem with an exact verifier.
**Scope:** HYP-8802 finds nonzero Weyl-ordering anomalies, so PC(2) false does not
yet give an explicit `A_2` endomorphism; the safe Hamiltonian-dual construction
lands in `A_4`. THM-2045 proves the displayed first coordinate has no planar
Jacobian mate.

### A4. Shao's Vanishing Conjecture — **UNTOUCHED (target).** Not named anywhere (the
only "Shao" is a different tournament author). No content, no witness.

### A5. Generalized Differential-Operator Vanishing Conjecture — **UNTOUCHED (target).**
Not named ("differential operator" hits are a tournament-matrix atlas).

### A6. Zhao's Image Conjecture — **UNTOUCHED as independent work; COROLLARY-false, no witness.**
klein-S332's cascade: Zhao's Image Conjecture ⟹ JC, so JC_{n≥3} false ⟹ IC false
by contrapositive — **but this rests on the cited external implication (unverified
in-repo) and produces NO witness** ("nobody has seen these objects"). Grade honestly
as UNTOUCHED with a logical corollary. ★ Next: the witness-extraction pipeline
(F → Yagzhev via THM-1325 → de Bondt–van den Essen symmetrization → Hessian-vanishing
witness) — the first explicit IC-violating object.

### A7. Mathieu / Mathieu–Zhao Subspace Conjectures — **UNTOUCHED; same corollary, no witness.**
Same status as A6 (the relevant Mathieu-subspace statements imply JC ⟹ false by
contrapositive if the implication holds; no independent work, no witness). ★ Same
pipeline; the radical/Abel-Ruffini structure (A8) is the Mathieu–Zhao-adjacent tool.

### A8. Exact classification of Jacobian counterexamples — DEEP (the repo's strongest JC contribution).
Verdict (4 convergent authors): **sporadic at the core, familial by propagation —
the counterexample set is the monoid generated by rigid seeds.**
- **THM-1305** [PROVED/PARTIAL] — equivariant normal form; c₂=0 FORCES the cube
  A=v^{k+1} at every k; the rationals all DERIVED; infinitesimal RIGIDITY (moduli
  tangent 0); k=3 rung EMPTY (mod-p + Newton, both k=2-validated) ⟹ ISOLATED SPECIES.
- **THM-1320/THM-1325** [PROVED] — the +1 IS the Yagzhev X; the reduction hides the
  TORUS (dim≤3 equivariant cubic-linear Keller maps all injective); the SURGERY WALL
  (the +1's zero-hyperplane kills the cofactor field). [THM-1320 P3 amended trivial.]
- **Fine structure** [PROVED]: full ℂ*-torus (collision = λ=1 slice of an orbit,
  λ↦λ² double-cover; det's "2" = torus weight); S₃ monodromy / Chebotarev (½,⅙,⅓),
  4 independent proofs; fiber = 1+2 Rédei-shape; the **master quartic L=D=K** triple-
  derived (x/y/u coords), simultaneously fiber-cubic leading coeff + resolvent
  conductor + **Jelonek asymptotic variety** (first explicitly computed for any JC
  counterexample) — THM-1310/1315; resolvent = orientation double-cover, Rédei-sign =
  discriminant character (opus-S418); surjective + everywhere-étale, ramification
  entirely at infinity (THM-1315); the **Chebyshev trisection** reading (THM-1335:
  "JC₃ fails the way angle-trisection fails, cos(θ/3)∉ℚ(cosθ)"); the **cuspidal-cubic
  engine trichotomy** (mac-mini THM-1340: deg-1 unit / deg-2 impossible / deg-3
  minimal; all known counterexamples ride ONE projectively-rigid cusp).
- **Keller monoid & realization** (THM-1330): 𝒳_n=deg⁻¹({≥3}), no degree 2 (Campbell),
  finite factorization, classifying = the inverse-Jelonek problem; the **Smith
  selection rule** (klein-S325: deg-2 impossible over every ℂⁿ, deg-3 forces S₃).
- **Radical inverse & Abel–Ruffini dichotomy** (kind-pasteur, THM-1345-radical):
  every known counterexample is radical-invertible (Cardano, solvable S₃); an A₅
  quintic rung would be the FIRST non-radical one — the polynomial incarnation of
  quintic unsolvability. [REFRAMED/CONJECTURED — the A₅ rung is named, not built.]
- **Big open questions**: Mount-Everest seed-uniqueness at (3,3) [CONJECTURED];
  realizability of degrees 4–7 (the n=4 (z,w)-affine hunt has its first obstruction,
  THM-1340 §5) [PARTIAL negatives]; the order-{1,3} conjecture (z-affine ℂ³ Keller
  field degree ∈{1,3}) [CONJECTURED]; the resolvent-conductor = Jelonek conjecture.

### A9. Two-variable Jacobian Conjecture (JC_2) — PROVED (equivariant); full JC_2 OPEN, difficulty LOCATED.
- **THM-1345** [PROVED category-restricted + REFRAMED] — det JF={P,Q}; JC_2 = "a
  canonical pair is a coordinate system" = classical shadow of DC_1. Equivariant
  Keller ⟹ invertible for EVERY ℂ*-action (hyperbolic→linear [boxeph-S144],
  elliptic→triangular). Difficulty located: leading forms Poisson-commute
  {P_A,Q_B}=0, propagation down the weight filtration = the AMS-hard open content.
- boxeph-S144: the dim-2 no-go + the **mod-p decision theorem** (Keller automorphy ⟺
  bijective mod one large prime). mac-mini-S137: the **golden/worst-approximable
  degree corner** (Lamé-for-polygons) + the engine-dimension lemma (n=2 engine-
  starved — why n=3 fell first). klein-S329: the Euler–Zariski bootstrap
  (JC_2@3 ⟺ one ramification parabola cannot be pushed to infinity).

---

## B. The Lonely Runner cluster — LRC(14) NOT closed; deep exact structure around the wall.

Top novel results (full detail: the LRC14-* frontier docs + THM files):
1. **Covering-min deep well 14/183 = n/Φ₆(n)** [PROVED single-killer closed-form;
   PARTIAL multi-killer] — Eisenstein norm, Heegner −3; the global bound IS LRC(n)
   itself (THM-724/726, the Eisenstein reflection).
2. **The D-graded gate tower = primorial cascade** [PROVED 3 seals, Lean-checked;
   CONJECTURED general-N] — F_D(N) attains D/((N+1)D−1) iff N≡1 mod L_D, N≢1 mod 2D−1
   (2D−1 prime), L_next=L_D·(2D−1); a Proth-adjacent gate law. 8 out-of-sample tower
   confirmations (THM-1285/1286/1271).
3. **Kernel-exact Lean certificate spectrum** [PROVED] — 3/23, 4/127, 4/247, 4/367,
   6/1271 machine-checked kernel-pure (3/23 axiom-free).
4. **The K-ladder M(K_c(N))=c/(cN+1) + the 12m/13s ladders** [PROVED floor + seals;
   one shared residual with the tower] — THM-1295; the bottom spectrum = two ladders
   rooted at the AP.
5. **The cross-N first-gap band** [PROVED single-far N=6..13] — THM-1284; the band-law
   conjecture refuted by its own author's follow-up (N=31 = 4/127, opening the tower).
6. **The sharp measure-horn 1/(7L)** [PROVED local; CONJECTURED uniform] — THM-1132/1123.
7. **Rung theory + the θ-flow + the observer principle** [REFRAMED] — opus-S409/S410;
   Rédei = LRC + 1 (opus-2026-06-30); the OCF vacuum digit = the +1 = ℏ.
8. **Historical Heegner/Eisenstein analogy** [NOT AN LRC BRIDGE; MISTAKE-229/230/233]
   — the separate class-number-one, Eisenstein-norm, Moser-spindle, and
   three-distance facts survive, but no map shows that a Heegner `-7` field or
   class-number rigidity governs the LRC floor. Do not use this analogy as a
   proof route without an explicit loneliness-preserving construction.
**Frontier — three named walls**: Wall A = the inverse/rigidity Freiman/near-AP core
(HYP-7310 is a strong sufficient AP-extraction supplier, **not an equivalent
restatement of LRC(14)**; the covering-core gap is verified for 95% of spread cores,
5% rational-time-evasive residual); Wall B = the six-comb phase-transport wall; Wall C
= the bound-D shell. Bounded height: exhaustive to 55 (THM-1290), floor isolated at
all heights (THM-1289, published), δ ineffective.

---

## C. Geometry / discrete geometry

- **Unit distance / Hadwiger–Nelson / Moser spindle** [PARTIAL/DEEP] — the largest
  secondary vein (~70 files). THM-412 density quantization (w | r_Q(D), triangular
  lattice skips densities 4,5, first density 6 at the split prime D=7); the HN field
  tower (HYP-2276/2277): the most rigid χ=4 junctions land on the Heegner numbers
  {7,11,19,43,67,163}; 3n floor bounds (THM-421); u21=57, u22 (THM-431/440); the
  J₀ spectral floor χ(ℝ²)≥3.48. The self-audit `missed-important-problem-frontier-s657`
  ranks HN #1 for stride potential. ★ Promote to a first-class thread.
- **Three-gap / three-distance (Steinhaus)** [REFRAMED, clean] — the Eisenstein/cusp
  dichotomy IS the three-distance theorem; the LRC-AP tight locus ⟺ ≤3-gap config
  (HYP-2913, verified n=4..7).
- **Kakeya / Falconer (finite field)** [MODERATE] — parabola carrier size exact minima
  7,17,31; Kakeya-as-adaptive-graphic-rank.
- **Sphere packing dims 8 & 24 (E8/Leech via {7,21})** [MODERATE, cross-problem] —
  Fano→octonion→[8,4,4]→E8 and T₂₃→Golay→Leech chains; Cohn–Elkies magic function =
  Chebyshev equioscillation extremal. The "why 7 and 21" narrative.
- **Zaremba's conjecture** [MENTIONED] — the named target for JC_2's Lamé-for-polygons.

---

## D. Under-the-radar problems (PROMOTE these — genuine results, hidden by LRC file-naming)

- **Collatz conjecture** [DEEP — top promote] — the rapidity conservation law
  ln n = K ln2 − L ln3 − D(n) with a bounded sign-definite harmonic defect; the exact
  identity n·3^L·Π(3a_i+1)/(3a_i)=2^K (Fraction-verified); the θ=2 Cramér–Lundberg
  excursion exponent. Self-contained novel result set (collatz-rapidity-defect,
  collatz-iterated-log-tower; HYP-2147-2149). **The clearest hidden gem.**
- **Erdős Problem 592** [DEEP — $1000 problem] — the p=2 Schur-seam theorem (THM-469),
  R(n,2)=2n+1 linear conjecture, Chang-tower m=3 open; 21 scripts, a survey draft.
- **{7,21} forbidden H-spectrum** [DEEP, repo-signature] — the only unachievable OCF
  values; H=7 (THM-029) + H=21 (THM-079) proved; parallels the open LRC {12,24};
  E8/Golay/octonion home. (OPEN-Q-028.)
- **GLMY path homology of tournaments** [DEEP / HIGHER DEGREES DISPUTED] —
  `β₁∈{0,1}` and `β₂=0` are the convention-safe core. Claims in degree at
  least three, including `β₄(T₇)=6` and the n=8 seesaw, are not current canon;
  see the [active convention case](../02-court/active/CASE-path-homology-regularity-convention.md).
- **Tournament reconstruction conjecture** [MODERATE→DEEP] — the (OCF, det) key; first
  OCF-cospectral twin at n=6; degeneracy metrics. Ranked #2 in the repo self-audit.
- **Hadwiger's conjecture on the metagraph** [RESOLVED at n=7] — `G_7/Z_2`
  has a certified `K_6` minor and Hadwiger number at least `12`; the live
  question is growth of `h(G_n/Z_2)`, with the n=8 lower bound already `22`.
- **Caccetta–Häggkvist** [MODERATE] — the return-residue reframe.
- **Real-rootedness of tournament independence polynomials** [CLASSIFIED AT
  THE FIRST FAILURE] — proved for `n<=8`, refuted from `n=9` by THM-025; the
  open problem is characterization of the real-rooted subclass.
- **Erdős–Selfridge / odd covering systems** [MODERATE] — covering-systems ↔ danger-arc
  correspondence (macmini-S97).
- **Bernoulli 1806 fixed point / von Staudt–Clausen** [MODERATE] — 6→42→1806→1806.
- **Pisano periods** [MODERATE] — π_F(10)=24 (the JC tripling clock); LRC fibered band.
- **Doubly-regular / skew-Hadamard tournaments ↔ Golay/Leech codes** [MODERATE].
- **Lee–Yang zeros / chromatic roots (HN gadget)** [MODERATE].
- Also present [TANGENTIAL]: Erdős-Straus, Fermat-Catalan/Beal, Erdős-Moser, Erdős-870,
  Erdős-64 power-of-two cycles, Goldbach/twin-prime lenses, Zeckendorf, the figurate
  zoo (Moser/Faulhaber, HYP-8165/8170/8175), the "cancellation family" meta-result
  (THM-406, with an MRDP proof they are NOT one problem).
**NICHE additions (death-star-S59w miner, beyond the S59u DEEP top-promotes) — ranked:**
- **★ Alcuin number vs Robertson–Seymour** [MODERATE] — Alcuin = τ+1 BREAKS minor-closure (monotone under deletion, INCREASES under contraction; K_{2,4} witness); so {G: Alcuin ≤ k} has no finite forbidden-minor set. Self-contained, OPEN-Q-107 (the-plus-one-that-escapes-through-contraction).
- **Hadwiger's conjecture on the metagraph** [RESOLVED n=7, death-star-S59x] — G_7/Z_2 (272 vtx, ω=4<χ=6) has a certified K_6 minor (economical: {2},{3},{10},{16},{0,7},{19,28}, 8 vtx) AND Hadwiger number ≥ 12 ≫ χ=6 — Hadwiger holds with large margin; the metagraph's density is in its minors, not its cliques. NEW open: how does h(G_n/Z_2) grow (compute at n=8)? (klein-S198).
- **★ Real-rootedness of I(Ω,x)** [DEEP for tournaments] — PROVED real-rooted n≤8 (claw-free THM-020), DISPROVED n≥9 (THM-025, explicit Newton-violation); OPEN = characterize the real-rooted subclass (OPEN-Q-047/051).
- **EGZ zero-sum ↔ H-spectrum residue gaps** [PARTIAL] — H ≢ 2 (mod 5) at n=5, predicted for all prime n (a crisp OCF-machinery target).
- **Bernoulli 1806 / von Staudt–Clausen {2,3,7}** [MODERATE] — B_6 = 1/42 denominator primes {2,3,7} = the same {2,3,7} as the OCF/H-spectrum semigroup? niche, unworked.
- **Erdős–Selfridge / odd covering systems** [MODERATE] — the LRC danger-arc cover IS a covering system; import Hough's minimum-modulus bound (self-flagged most-promising bridge).
- **Zarankiewicz (cyclic book)** [MODERATE] — Z(m,m) achieved m≤7; general-m bipartite proof open (THM-922).
- **Caccetta–Häggkvist** [MODERATE] — the return-residue reframe (sumset-growth ↔ short-return).
- **Erdős–Moser two-step tower** [MODERATE] — T127=15, T255=23 exact; why does the recurrence persist?
- **Seymour 2nd-neighborhood + Sumner's conjecture** [GREENFIELD] — named tournament fronts, ZERO work, OCF/score machinery directly applies.
- **BRAND-NEW fronts (0 hits anywhere = pure greenfield on tournaments):** Manickam–Miklós–Singhi, the 1-2-3 conjecture, antimagic labeling, Tuza's conjecture, the toughness conjecture, Hadwiger–Debrunner (p,q), the Davenport constant, the cap set problem, Alon–Tarsi (currently only a tool). Each is a clean opening the repo has never touched.
- **★ LEDGER-MISS (S59z miner):** a paper-ready **Cayley/Delannoy identity cluster** is absent from ALL ledgers — `drafts/candidates-for-sharing.md`: the master GF Q(x)^m with duality k·g_k(m)=m·g_m(k); the CV²(H) variance with provable 1/n² cancellation (g₂(m)=g₁(m)²); OEIS-new W(n); THM-224 golden exceptional points (char poly λ³−λ²−xλ−x, EP eigenvalues 1/φ,−φ); the bilinear(Tustin)=Delannoy DSP bridge. These are the MOST Lean-able (pure identities) and conventionally publishable items in the repo (JCTA/EJC-shaped). FOLD IN.
- **TOP LEAN TARGETS (S59z triage, pure identity / finite decide):** THM-1300 §0 (JC counterexample kernel, det≡−2 + triple collision by ring/decide — the highest-leverage formalization, zero Mathlib AG); THM-1300 §1 (18 A₃ Dixmier identities by ring); THM-025 (real-rootedness disproof n=9); the {7,21} spectrum note; the Cayley/Delannoy identities. Already done+submittable: RedeiFromOCF.lean (axiom-free H-odd), the 4/367 & 4/247 kernel-pure certs.
- **HADWIGER on the metagraph — n=8 done (S59y):** h(G₈/Z₂) ≥ 22 (V=3528, ω=5, χ=7); h/χ steepens 2.0→3.1; a minor-dense/low-clique/low-chromatic family (publishable).
- **HYGIENE — SYSTEMIC THM-id collisions:** ~70 THM numbers doubled (THM-1345/1370/1375/1380/1385 all doubles; also 201/260/262/290/338/868/869/922). The MISTAKE-199 pattern at canon scale — a real citation obstacle needing a de-collision pass.
- HYGIENE: a parallel `PROBLEM-PORTFOLIO-2026-07-20.md` (mac-mini-S140) also exists; this file (PROBLEM-LEDGER.md) remains the canonical one — merge, don't fork.

**Repo self-audit to start from**: `07-reflections/missed-important-problem-frontier-s657.md`
(ranks HN > reconstruction > Kakeya/Falconer > sunflower > Caccetta–Häggkvist > … ).

---

## E. Hygiene — collisions to reconcile (this prompt + prior)

- **HYP-8185 is a 4-way collision** (death-star-S59u, klein-S332, kind-pasteur-S128c104,
  mac-mini-S140) — this file is the consolidation; the others' HYPs should point here.
- **THM-1345 is a DOUBLE** — `…poisson-reframing-dc1-shadow` (death-star-S59q) and
  `…plane-family-…-radical-inverse` (kind-pasteur-S128c101), both legitimate, both
  filed 2026-07-20. **Needs de-collision** (renumber the second, or merge).
- MISTAKE-199 (fleet-wide-prompt duplication) recurred here; caught mid-session this
  time. Lesson standing: grep concurrent same-prompt claims BEFORE deep work.

---

## F. Future priorities (for current and future sessions)

1. **★ THE WITNESS EXTRACTION (assessed S59z: novel, meaningful, the top publishable + Lean-worthy move — DO IT).** Transport F through Yagzhev (cubic-homog) + de Bondt-van den Essen (symmetric quartic, Hessian-nilpotent) to produce the FIRST EXPLICIT counterexample to Zhao's Vanishing Conjecture, the Image Conjecture, and a failing Mathieu subspace. Fleet-unanimous #1, genuinely UNTOUCHED ("no explicit witness exists anywhere", opus-S421). Meaningful because equivalences preserve TRUTH not WITNESSES — the objects are unbuilt. ROBUSTNESS RULE: verify the witness DIRECTLY (finite Δ^m(P^m) certificate — machinery built S59z, Lean-able), NOT via the VC⟺JC equivalence (unverified in-repo). Publishable standalone (conditional on THM-1300 + the JC disproof's acceptance). Feasibility gate = the dimension the F-reduction lands in — **PARTLY EXECUTED, gate NOT cleared by stacking (death-star-S61, exact ℚ(i) computation).** VERIFIED exactly: F Keller + triple collision; Yagzhev G (det JG≡1) with collision transported; **but JH(G) is NOT nilpotent** (nilpotency must be created). VALIDATED exactly: the cotangent-lift + de Bondt rotation machinery (T′=[[I,iI],[iI,I]], √2 cancels → exact over ℚ(i)) turns ANY nilpotent-Jacobian cubic map into a Hessian-nilpotent homogeneous quartic (∂ⱼP=(∇P)ⱼ and Hess(P) nilpotent both verified on 2 test maps). **BLOCKER (real math, not engineering):** producing a cubic-**homogeneous KELLER** reduction of F. Companion/stacking moves reduce degree and transport the collision but are Keller only ON the section {W=X^β} — off it det varies (first move breaks det≡1); naive homogenization also breaks Keller (det(I+x₀JH₂+JH₃) not const). So the earlier "6 helpers→N≈10→M≈20, feasible, no new math" (relayed from the reduction agent) **CONFLATED the easy degree-reduction with the hard Keller-preserving homogeneous reduction** — the witness EXISTS (BCW guarantees it) but its construction needs global-determinant control (genuine BCW/Drużkowski content). ≈20 was contingent on the correct reduction's size, not established. [≈76 was a separate per-monomial-heuristic error — MISTAKE-201; the feasibility overstatement is its addendum.] What remains: implement BCW homogeneous reduction (or Drużkowski cubic-linear) with det held globally constant, then feed the validated lift+rotation machinery. Refs: de Bondt *Symmetric Jacobians* arXiv:1206.2865; Zhao arXiv:0704.1689; BCW Bull.AMS 7 (1982); Drużkowski.
  **UPDATE death-star-S61b — the S61 global-det bug is FIXED at the framework level.** The correct BCW move is a SHEAR COMPOSITION E₁∘(F⊕id)∘E₂ with E₁,E₂ elementary shears (det J≡1 identically): concretely (vars,u)↦(F(vars)+h(u+g(vars)), u+g(vars)) = [shear vars by h(u)]∘(F⊕id)∘[shear u by g(vars)]. DEMONSTRATED exactly (`bcw_shear_globaldet_deathstar_S61b.py`): (A) det stays GLOBALLY constant =det JF=−2 at all random points for arbitrary g,h — unlike S61's companion move which was Keller only on-section; (B) with a perfect-power term, one shear drops the top degree (deg-6 (xy)³ → deg 5) while det stays globally =1. So the framework is right; the remaining specific work is decomposing each monomial into perfect-power-reducible pieces and iterating to cubic-homogeneous (research agent retrieving the exact BCW/Drużkowski shears). NOTE the whole downstream is DONE (klein THM-1435: de Bondt conjugation as a matrix identity, cotangent-lift symmetry one-line, positive+negative controls; kp THM-1430: symmetric ℂ⁶ object) — hand any cubic-homog nilpotent non-injective map to klein's pipeline and P falls out. This session did NOT complete the homogeneous reduction; it fixed the framework and demonstrated global-det preservation.
1. **The other untouched exotic conjectures — first contact.**
   Shao's vanishing (A4), the generalized differential-operator vanishing conjecture
   (A5), Zhao's image (A6), Mathieu/Mathieu–Zhao subspaces (A7): NAMED, unworked, and
   the repo owns the exact machinery — the witness-extraction pipeline (THM-1325
   Yagzhev + de Bondt–van den Essen), the equivariant/Poisson calculus (THM-1345), the
   radical/Abel-Ruffini structure, the constructive Weyl-endomorphism (THM-1300 §1).
   **The single highest-leverage move: produce the FIRST EXPLICIT WITNESS** for the
   Image / Mathieu-subspace / stable-Poisson conjectures from F — objects nobody has
   ever seen. One focused session each.
2. **The deepest live frontiers.** JC_2/DC_1 (the weight-filtration induction; the
   two-lens tournament strategy; is propagation length = the CF chain length?); the
   JC-classification realization program (degrees 4–7; the A₅ Abel-Ruffini rung; the
   resolvent-conductor=Jelonek conjecture); LRC(14)'s Wall A (the Freiman core).
3. **Promote the hidden gems** to first-class threads: Collatz, Erdős-592, {7,21},
   GLMY, reconstruction — each has real results buried under LRC file-naming.
4. **The connective tissue — integrate, don't silo.** The observer +1 (OCF vacuum
   digit = ℏ = Yagzhev X), the doubling/n-vs-2n ladder, parity-protection (Rédei-odd
   ↔ odd-degree), the Poisson/symplectic frame, and the Heegner/Eisenstein arithmetic
   run ACROSS these problems. Progress on one transfers; the ledger exists to make the
   transfers visible.
