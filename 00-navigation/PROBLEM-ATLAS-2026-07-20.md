# THE PROBLEM ATLAS — post-JC frontier ledger and priority map

> **SUPERSEDED DRAFT:** consolidated by [`PROBLEM-LEDGER.md`](PROBLEM-LEDGER.md).
> Use [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md) for post-2026-07-20 status;
> retain this file for credited provenance and discarded perspectives.

**kind-pasteur-2026-07-20-S128c104 (HYP-8185).** Owner directive: given the
Jacobian Conjecture is disproved externally by a singular counterexample, compile
the repo's NOVEL, NON-OBVIOUS results per problem; expand the list from
under-the-radar threads; prioritize creative, assumption-challenging progress.
Format per problem: STATUS → NOVEL REPO RESULTS (cites) → SHARPEST OPEN →
NEXT-SESSION SHAPE → PRIORITY. Maintain this file: update the per-problem blocks
as results land; date-stamp edits.

---

## A. Classifying Jacobian counterexamples EXACTLY
- **STATUS:** the fleet's 3-day swarm produced a near-complete structure theory
  of the known counterexample and strong rigidity toward "unique degree-3 seed."
- **NOVEL:** conic-pair fibers, depressed cubic (trace law Σxᵢ = 0), S₃ via
  resolvent √−L with pointwise χ(−L) law 0/8316, explicit Jelonek quartic
  {L=0} = cusp {P³+E²=0} — A₂ at both ends (THM-1310, THM-1335); master
  identities (map ≡ two lines); Chebyshev modulus m = essential-class invariant;
  trisection reading; radical inverse for the sporadic tower (THM-1345), now
  superseded globally by THM-3438's explicit non-radical `S_5` weighted atom;
  trace-polynomiality/centroid + propagation to F∘F;
  Smith selection rule (klein-S325: deg 2 impossible; d=3 ⟹ S₃; d=4 ⟹ A₄/S₄
  only; dihedral iff odd); surjectivity (mac-mini THM-1315); engine trichotomy
  (THM-1340); TRAP/Druzkowski strata, +1 = Yagzhev X (THM-1320/1325); ghost
  theorem (boxeph-S145, ANY dimension); rigidity from four independent probes
  (opus deformation 1+4-obstructed; death-star chart+W3-refutation; boxeph
  lifting uniqueness; my design boxes + equivariance-uniqueness "why cubic");
  towers/conjugates as the propagation monoid (klein S326/S327); resolvent =
  orientation cover of the fiber 3-tournament, Rédei sign = disc character
  (opus-S418).
- **SHARPEST OPEN:** the SEED CONJECTURE (unique deg-3 seed mod Aut×Aut) and the
  REALIZATION PROGRAM (T1549): THM-3438 now realizes `S_d` in every degree
  `d>=3`, including z-quadratic `S_4` and non-radical `S_5`; classification
  of other groups such as `A_4/A_5`, equivalence classes, and constrained
  z-affine/planar strata remains open.
- **NEXT:** the 2-jet design equations at d=4 (one hard session; solver
  patterns ready); Lean certification package v3 (masters + 108a²L + section +
  traces + det/collisions — plausibly the first formal analysis of any JC
  counterexample anywhere).
- **PRIORITY: HIGH** (headline-grade formalization + the d=4 door).

## B. Dixmier conjecture (DC_n)
- **STATUS:** false for n ≥ 3 (transfer, THM-1300); survivors = the linearly
  ordered ladder **DC₂ ⟹ JC₂ ⟹ DC₁** (c97) — refute one, everything above
  falls; prove one, everything below follows.
- **NOVEL:** the explicit A₃ endomorphism program (mac-mini/death-star
  HYP-8075/8080 — Weyl-relation refereeing; check current status before
  building on it); tournament Dixmier property holds n = 4, 5 with tame group
  {id, op} (mac-mini-S132); DC₁ VIA TOURNAMENTS (death-star-S59r HYP-8160): the
  A₁ weight-triple = an oriented 3-cycle, observer = ℏ = conserved +1, two-lens
  strategy (Rédei parity × 2D leading-form bound) — a genuinely new attack
  frame for the 1968 problem.
- **SHARPEST OPEN:** DC₁ itself; assumption-challenging: the fleet should hold
  BOTH directions live (the tournament atom could power a counterexample hunt
  for DC₂ as well as a proof route for DC₁).
- **NEXT:** finish the explicit A₃ witness (constructive ¬DC₃); push the
  S59r two-lens program to a first theorem about A₁ endomorphism normal forms.
- **PRIORITY: HIGH** (DC₁ is the bottom-of-the-tower crown).

## C. Stable Poisson conjecture
- **STATUS:** FALSE already at rank two. THM-2044 gives an explicit
  nonautomorphic endomorphism of `C[x,q,p,z]` preserving the standard Poisson
  bracket, with an exact three-point fibre. Hence PC(n) is false for every
  `n>=2` by identity padding; PC(1) is the surviving Poisson floor.
- **NOVEL:** THM-2044 reconstructs the rank-two witness as a polynomial
  symplectic suspension of THM-1300; det JF = {P, Q} Poisson/Hamiltonian reframing (death-star-S59q);
  UNDER-RADAR PRIOR WORK REDISCOVERED THIS SESSION: the anti-Poisson coimage
  atlas + coimage-Yoneda 2n−1 resonance (reflections s604/s605, old era) —
  unintegrated with the new JC corpus.
- **NEXT:** HYP-8802 direct `A_2` quantization: five cubic Moyal anomalies are
  nonzero, while finite corrections repair the entire `R`-column exactly.
  Solve the coupled `T`-column and flatness relations simultaneously or prove
  the correction tower cannot terminate;
  integrate s604/605.
- **PRIORITY: MEDIUM-HIGH.**

## D. Zhao's vanishing conjecture / generalized differential-operator vanishing
- **STATUS:** NEW target (no prior repo work). The JC-equivalent instances
  (Zhao 2007: VC for Laplace-type operators ⟺ JC) are now FALSE in the
  corresponding parameters; the broader operator families are NOT settled.
- **THE UNIQUE ASSET:** the repo holds the explicit F *and its anatomy* —
  constructing the EXPLICIT vanishing-conjecture witness polynomial (through
  the equivalence, from F's components/master identities) is a deliverable
  nobody outside can produce as concretely. Which operator classes survive =
  the new frontier map.
- **NEXT:** one literature+construction session: state the equivalence chain
  precisely, push F through it, exhibit the witness, and mark the survival
  boundary (which Λ, which degrees).
- **PRIORITY: HIGH (fresh, concrete, assumption-challenging).**

## E. Zhao's image conjecture and Mathieu subspaces
- **STATUS:** NEW target. The Mathieu-subspace FRAMEWORK is a definition and
  survives; conjecture instances implying JC fall; instances tied to JC₂/DC₁
  remain open and inherit the survivor-ladder structure.
- **NEXT:** map the fall-line through the image-conjecture family (which
  measures/subspaces); identify the minimal surviving Mathieu statements —
  they are then equivalent currency for the B/F survivors.
- **PRIORITY: MEDIUM.**

## F. Two-variable Jacobian Conjecture (JC₂) — THE main open positive target
- **STATUS:** open; now carries the whole classical weight (with DC₂ above it
  and DC₁ below it).
- **NOVEL — THE ASSEMBLED OBSTRUCTION CAGE (fleet synthesis):** any JC₂
  counterexample must have: field degree ≥ 3 with NON-GALOIS closure (Campbell
  + klein-S325 Smith: deg-2 étale self-covers impossible in ANY dimension);
  CUSPED Jelonek curve that cannot be pushed to infinity (klein-S329
  Euler–Zariski bootstrap; JC(2)@deg3 ⟺ one ramification parabola); it cannot
  be equivariant under ANY ℂ*-action (death-star-S59q theorem, elliptic +
  boxeph hyperbolic cases); it cannot have all-ghost asymptotic set
  (boxeph-S145 ghost theorem); Euler-ledger χ-profile constraints + d=2
  all-smooth impossible + the odd-degree conjecture (boxeph-S146, which also
  ties Keller degree structure to the H-spectrum holes {7,21}!); N₂ ≡ 0 as
  étaleness detector and Moh's degree distinction (klein-S329).
- **SHARPEST OPEN / ASSUMPTION-CHALLENGING:** the fleet should run JC₂ from
  BOTH ends: (i) close the cage (each obstruction is a rung; the cage is
  finite-looking at low degree), and (ii) HUNT inside the allowed region with
  the dim-3 architecture lessons (the killer existed in dim 3 because a plane
  engine had a compensating line; in dim 2 there is no spare line — that's an
  intuition, not a theorem — challenge it with 2-jet-style ansatz hunts).
- **PRIORITY: HIGHEST pure-math target.**

## G. Lonely Runner Conjecture LRC(14) — the flagship pre-JC corpus
- **STATUS:** first open case (13 runners); LRC(≤13) settled by owner directive.
- **NOVEL (highlights of ~100 sessions):** THM-651 shifted-tent k=8 leg; the
  finite-check FEASIBILITY LEDGER (~50 core-years; the two structural walls:
  I(13,p,1) sieve + the ×7 composite lift) (c87); tightness cages (HYP-7920
  12-speed; c89 k=13 J-separator, height 281,577); THM-1288 (S-T Conjecture 7.1
  refuted) + THM-1289 (G-K floor isolation); the eternal-inhabitant theorem
  (no extinction prime — MISTAKE-192's rule); shell collapse + covers census;
  the escape-atlas + ABSOLUTE-WINDOW law (klein: escape moduli live in [43,48]);
  cross-N band law (death-star THM-1284/1295) + kernel-exact Lean members
  (3/23, 4/127, 4/247, 4/367); uniform max|S| ≤ C·diam (klein THM-887);
  resonant-mode repair (THM-883); the mod-p descent pipeline (HYP-8020 — the
  same philosophy as BKK's Azumaya descent, the JC bridge).
- **SHARPEST OPEN:** Wall A (the n=12 AP-uniqueness inverse complex), Wall B
  (six-comb phase transport), the covering-core gap (HYP-7748); the small-witness
  CRT law with klein's fixed escape window.
- **PRIORITY: HIGH** (it remains the deepest sustained thread; the JC episode
  sharpened its instruments — certificate rungs, instrument gates, descent).

## H. Unit distance problem — UNDER-RADAR, prior work FOUND
- **STATUS:** dormant thread with forgotten canon: THM-408 (Moser-layered slabs
  have unit spines), THM-412 (density quantization of the unit-distance layer).
- **NEXT:** one revival session: re-read THM-408/412 + their era's reflections,
  restate in current vocabulary, connect to the Moser-circle/figurate thread
  (both Mosers!) and the incidence machinery; identify whether the layer
  quantization gives anything for the Erdős unit-distance exponent or
  Hadwiger–Nelson-adjacent statements.
- **PRIORITY: MEDIUM (revive + integrate).**

## I–L. The under-the-radar expansion (mined this session)
- **I. Erdős 592 (sum-free, $1000):** THM-469 dichotomy, t=7 wall, R(n,2) = 2n+1
  candidate — NOTE the cross-thread: 2n+1 is the Rosetta/Proth rim (HYP-8165/70);
  the owner's "2n+1 is key" lands here too. One session: re-run the R(n,2)
  candidate against the new shear/Proth laws. PRIORITY: MEDIUM.
- **J. Pentagonal hub:** γ_pent ≈ 0.206 Lyapunov, η²⁴ = code discriminant,
  [72,36,16] obstruction (THM-487/488/489) — curio: π_F(10) = 24 (HYP-8145).
  PRIORITY: MEDIUM-LOW (deep well, no fresh lever yet).
- **K. Kakeya K(A₅) = 15 / K(2I) = 30 (klein THM-870):** icosahedral needle
  directions — SAME GROUP as the realization program's quintic rung and the
  Pisano-60 decode: an A₅ TRIANGLE of threads (Kakeya × Keller × Fibonacci)
  that no one has braided. One creative session. PRIORITY: MEDIUM-HIGH (novel
  braid).
- **L. The rest of the sleeper list:** path homology β₂ = 0 (n=8 threshold);
  real-roots conjecture (dies at n=9 — WHY?); DRT n=11; cospectral invisible
  pairs → n=8 double-cones (T1546); toothpick A139250 substrate; Feit–Thompson
  solvability quantum (+8 = one 3-cycle, T1533); the (r,g) shadow lattice
  Moser-region law + 8 OEIS-absent sequences (T1533); Zeckendorf/F3 exchange
  walk; additive-basis ladder (HYP-2998–3000); the OEIS SUBMISSION PIPELINE
  (klein T1532 batch + the c103/S59t/S148 shear harvest — a dozen new
  sequences now waiting; an afternoon of engineering); mod_rank/circulant
  homology PyPI deliverables (the engineering mandate's standing debt);
  H-spectrum {odd}∖{7,21} — now a THREE-thread crossroads (Ham-paths ×
  boxeph-S146 Keller degrees × Rosetta T(7,2)/row-3-sum): the {7,21} holes
  deserve their own session. PRIORITY: the {7,21} crossroads and the OEIS
  pipeline first.

---

## Priority queue for the next five sessions (recommendation)
1. **JC₂ both-ends session** (F): cage-closing + in-cage hunting (highest).
2. **The 2-jet d=4 architecture** (A) — also feeds the Abel–Ruffini and A₅
   braids (K).
3. **Zhao fall-line with explicit witnesses** (D) — unique-asset play.
4. **DC₁ tournament-atom program** (B) — push S59r to a first theorem.
5. **The {7,21} crossroads + OEIS pipeline + unit-distance revival** (L/H) —
   one integrative housekeeping-and-discovery session.

Assumption-challenging norms (owner's standing instruction): hold both truth
values live for JC₂, DC₁, DC₂; treat "the counterexample is unique" as a
hypothesis with a hunt attached, not a belief; every negative search publishes
its detection floor (the instrument-gate rule).
