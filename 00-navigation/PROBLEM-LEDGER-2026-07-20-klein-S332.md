# THE PROBLEM LEDGER — what we have actually proven, and where assumption-challenging progress is next
(klein-2026-07-20-S332; owner directive: catalog novel non-obvious results per problem,
expand the list from under-the-radar repo threads, and set priorities for future sessions.)

Legend: ★ = owner-priority next attack. Every claim below is repo-verifiable (file refs).

## 1. Dixmier conjecture DC(n)
- **PROVED-EXPLICIT (novel):** DC(n) FALSE for n ≥ 3 with the explicit Weyl witness
  φ_F : A₃ → A₃ (x ↦ F, ∂ᵢ ↦ Σ(J⁻¹)ₖᵢ∂ₖ) — all relations verified symbolically;
  injective (A₃ simple), non-surjective via the classical DC⟹JC bridge (S323).
  Also: EVERY Keller map quantizes (étaleness forces [Dᵢ,Dⱼ] = 0) — the bridge's
  content made computational.
- **Open:** DC(1), DC(2). death-star HYP-8160 opened DC(1) VIA TOURNAMENTS
  (A1 weight-triple = oriented 3-cycle — the odd sector speaking Weyl).
- ★ Next: exhibit a CONCRETE element of A₃ ∖ im(φ_F) with certificate; push the
  tournament dictionary at DC(1).

## 2. Poisson conjecture / Kontsevich–BKK circle
- **NOVEL, NON-OBVIOUS (stated this session):** the cotangent lift Φ(x,ξ) =
  (F(x), J⁻ᵀξ) is an EXPLICIT symplectic Keller counterexample in ℂ⁶ (det ≡ 1,
  Φ*ω = ω exactly, 3-point collision — S323-verified). Hence **Φ* is an explicit
  Poisson-algebra endomorphism of ℂ[x₁..₃, ξ₁..₃] that is injective and NOT
  surjective — the Poisson Conjecture (endos of the symplectic Poisson algebra are
  autos) is FALSE at n = 3, explicitly.** The BKK equivalences (stable JC ⟺ DC ⟺ PC)
  are now all vacuously-consistent falsities with ONE explicit generator.
- ★ Next: write the one-page THM file with the Φ*-bracket verification cited; the
  Aut-group (not End) conjectures SURVIVE — mark the live boundary.

## 3. Zhao's web: Vanishing Conjecture, Image Conjecture, Mathieu(-Zhao) subspaces
- **THE CASCADE (novel, stated this session):** Zhao proved VC ⟺ JC; the Image
  Conjecture and the relevant Mathieu-subspace conjectures IMPLY JC. JC(n ≥ 3)
  false ⟹ **VC is FALSE, the IC is FALSE, and the corresponding Mathieu-subspace
  statements are FALSE** — by contrapositive, previously unstated in the repo
  (1 Zhao file repo-wide).
- ★ Next (assumption-challenging, concrete): EXTRACT THE EXPLICIT WITNESSES — run F
  through the Yagzhev/Drużkowski lift (death-star THM-1325 started: "+1 = the
  Yagzhev X") and Zhao's transfer to produce the first-ever explicit VC-violating
  pair (Λ, P) and an explicit non-Mathieu image. Nobody has ever seen these objects.

## 4. Classifying Jacobian counterexamples exactly
- **The frame is ours (novel):** the master quartic anatomy (D, disc = −4S²D, the
  universal fiber cubic; S324); the SMITH SELECTION RULE (degree 2 impossible over
  every ℂⁿ; monodromy needs self-normalizing stabilizers; deg-3 forces S₃; the
  allowed lattice through degree 7; S325/T1549); the Keller monoid picture
  (mac-mini THM-1330: units = deg 1, finite factorization, inverse-Jelonek);
  uniqueness in class (boxeph S142) + chart rigidity (klein shells; death-star
  THM-1305/1320); the essential-class invariant (kps THM-1335 trisection modulus);
  fiber spectrum {3,1}-never-2 = the census structural zero.
- ★ Next: the REALIZATION PROGRAM — does any A₄/S₄ quartic étale self-cover exist?
  (First genuinely new irreducible would live there.) Identify F∘F's wreath
  subgroup; the D₅ quintic question.

## 5. Two-variable Jacobian conjecture JC(2) — the last stand
- **Novel:** the EULER–ZARISKI BOOTSTRAP (S329): cover-degree-3 counterexamples
  need a CUSPED Jelonek curve (Fulton–Deligne + Smith), exact ledger
  χ(H⁻¹J) = 3χ(J) − 2; minimal model = cuspidal cubic + universal root cover;
  **JC(2)@3 ⟺ one ramification parabola cannot be pushed to infinity.** Plus
  death-star THM-1345 (JC₂ is a THEOREM in the ℂ*-equivariant category, every
  action) and boxeph's Euler ledger (d=2 all-smooth impossible; χ-profiles) ⟹ any
  counterexample is torus-asymmetric AND cusped-Jelonek.
- ★ Next: the candidate-Jelonek atlas (χ = 1 cusped curves, B₃-type π₁ — expected
  very short); Newton-polygon tracking of the parabola under Jung–van der Kulk.

## 6. Lonely Runner Conjecture, n = 14
- **Novel plank set:** THM-1290 (sub-gap (1/14, 3/41) EXHAUSTIVELY empty to height
  64; LRC(14) verified to height 55 — both re-certified after MISTAKE-194); the
  SANDWICH (cage 281,577 conditional / census unconditional / THM-1289 G-K floor
  isolation / Conj-1.5 roof); the LADDER-LOCKING LAW (band {6}∪{N≡1 mod 6}
  mechanism = opus-S409's gate; N = 19 CONFIRMED out-of-sample; tight locus
  {AP, GW} = the m ≤ 2 plateau of ONE ladder); the two-ladder bottom spectrum
  (all 8 sub-1/13 families on-ladder to h = 45); the escape atlas (absolute window
  q ∈ [43,48]); slack = D − k unification with Fan–Sun.
- ★ Next: Wall A (HYP-7310 AP-rigidity) unchanged as the core; gate arithmetic on
  the five k ≥ 3 rungs; effectivize the G-K δ; the small-witness CRT law (kps).

## 7. UNDER-THE-RADAR THREADS MINED (repo file counts; each = a standing problem)
- **Unit distance problem (Erdős)** — 70 files. ★ Audit thread: what is proven vs
  tangent; connect to the circle/chord machinery (Moser circle S330) and the
  escape-atlas distance methods.
- **Sidon sets / B₂ sequences** — 64 files (LRC witnesses ↔ Sidon structure).
- **Collatz** — 39 files (!). Audit for citable results.
- **Sylvester(–Gallai / Frobenius?)** — 36 files.
- **Kakeya (finite-field flavors)** — 35 files; the covering-spectrum "K(A₅) twin"
  memory (klein) — a REAL bridge to the LRC covering machinery.
- **Hadwiger(–Nelson, chromatic number of the plane)** — 28 files + Moser spindle
  12 files: unit-distance chromatic — natural home for tournament/parity methods.
- **Sum-product / BGK** — 19 files; kps HYP-8020: the mod-p DESCENT OF WALL A into
  sum-product territory — a live cross-problem pipeline.
- **Sierpiński/Proth** — 10 files + THM-1355 (the n·2^x+1 unification) + the
  S331 Mersenne antidiagonal law.
- **Prouhet–Tarry–Escott** — 4 files; kps HYP-7955: PTE is load-bearing for LRC
  moment blindness. ★ PTE size-13 ⟺ the cage's Newton depth — name it a problem.
- **Markov/Lagrange spectrum** — 1 file (opus S408 M_θ deformation: LRC ⟷ Markov
  spectrum; Fibonacci maximizers) — embryonic, high ceiling.
- **Goddyn–Wong / n = 12 uniqueness ("Tao's optimistic conjecture")** — 111 files;
  = Wall A's shadow; the repo's single deepest open engagement.
- **The H-spectrum problem (tournaments: which Hamiltonian-path counts occur)** —
  the repo OWNS this: Rédei parity, boxeph's {7,21}-impossibility H-template +
  Keller-degree-monoid dictionary (S146) — an under-named open problem the fleet
  is actively solving; promote it to first-class status.

## 8. Priorities for the next sessions (owner-aligned)
1. ★★ Zhao-witness extraction (the first explicit VC/IC/Mathieu counterexamples).
2. ★★ Poisson THM file (one page; already-verified ingredients).
3. ★ Realization program (A₄/S₄ existence; F∘F monodromy ID).
4. ★ JC(2): Jelonek atlas + Newton-polygon parabola tracking.
5. ★ LRC(14): the five-rung gate arithmetic + Wall A.
6. ★ Under-the-radar audits (unit distance, Collatz, Kakeya-K(A₅), Sidon) — one
   session each: harvest citable results, promote or retire.
