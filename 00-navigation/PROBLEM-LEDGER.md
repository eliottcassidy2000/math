# PROBLEM LEDGER — the repo's novel results across external conjectures, and the frontier for future work

**Created death-star-2026-07-20-S59u (HYP-8185)**, owner-directed: given the
Jacobian Conjecture was disproved externally by a singular dim-3 counterexample,
inventory our NOVEL, non-obvious results on the surrounding problems, expand the
list with under-the-radar threads, and frame future creative work. This is a
LIVING index — future sessions should update the grade and frontier of each entry
and add new problems. Grades: **PROVED** (in-repo proof) · **PARTIAL** (real
partial result / verified sub-case) · **REFRAMED** (new formulation that relocates
the difficulty) · **CONJECTURED** (evidence, no proof) · **UNTOUCHED** (named, not
yet worked — a target). Every claim carries a citation; unproved things are
labeled. LRC(≤13) is treated as settled by citation per owner directive.

---

## A. The Jacobian ecosystem (van den Essen's world)

### A1. Jacobian Conjecture (JC_n)
**Status: DISPROVED externally (n ≥ 3); verified in-repo.** The owner-supplied
map F: ℂ³→ℂ³ has det JF ≡ −2 and a triple collision — **THM-1300** (six
independent in-repo verifications; JC_n false for all n ≥ 3 by stabilization).
**JC_2 survives** (see A8). Frontier: the exact classification (A5) and JC_2 (A8).

### A2. Dixmier Conjecture (DC_n)
**Status: DC_{n≥3} constructively FALSE (PROVED in-repo); DC_1, DC_2 OPEN.**
- **THM-1300 §1** — the EXPLICIT A_3 Weyl-algebra endomorphism φ(X_i)=F_i,
  φ(D_j)=Σ_k B_jk D_k with B=(JF^T)^{-1} over ℤ[1/2]; all 18 Weyl/flatness
  identities verified; non-surjectivity by the module one-liner (X_1 ∉ im φ).
  A_3 has an explicit proper self-embedding ⟹ DC_{n≥3} false constructively.
  (Independent constructions: klein's symplectic ℂ⁶ lift; boxeph's lead.)
- **DC_1** — REFRAMED via tournaments (**S59r, HYP-8160**): A_1's weight triple
  (N, q, p) IS the oriented 3-cycle (observer = ℏ = the conserved +1); the
  two-lens strategy (Rédei parity × 2D leading-form bound). DC_1 NOT claimed
  (open; JC_2 ⟹ DC_1 via BKK). See A8.
- **DC_2** — OPEN; decoupled from the dead JC side (JC_4 false makes the BKK
  route vacuous); hangs below JC_2.

### A3. Stable Poisson Conjecture — **UNTOUCHED (target).**
The symplectic structure IS present (klein's cotangent lift Φ=(F, J^{-T}ξ),
Φ*ω=ω, an explicit symplectic Keller counterexample in ℂ⁶; THM-1345's
det JF = {P,Q} Poisson reframing), but the *stable Poisson conjecture by name*
is not worked. **Tool that could bear on it:** the equivariant/Poisson machinery
of THM-1345 (weight-additive bracket, leading-form obstruction) applied to the
Poisson-algebra endomorphism setting.

### A4. Shao's Vanishing Conjecture — **UNTOUCHED (target).** [verify via miner]

### A5. Exact classification of Jacobian counterexamples
**Status: DEEP partial results — the repo's strongest Jacobian contribution.**
- **THM-1305** — the equivariant anatomy: normal form for weight-(1,−1,−k)
  maps, the s-graded determinant law, c₂=0 forces A=v^{k+1} (the cube is a
  theorem), infinitesimal RIGIDITY (moduli tangent 0), and the k=3 rung is
  EMPTY on two validated instruments ⟹ the counterexample is an ISOLATED
  SPECIES (sporadic, not a family) within the equivariant chart.
- **THM-1320/THM-1325** — the +1 IS the Yagzhev identity part X (not hidden);
  the reduction hides the TORUS; the SURGERY WALL (no single-variable
  stabilization; the +1's zero-hyperplane kills the cofactor field). [THM-1320
  P3 honesty-amended to det JF(0), trivial.]
- The fine structure (fleet): S_3 monodromy / Chebotarev split (1/2,1/6,1/3);
  fiber = 1+2 = Rédei-shape; the resolvent = orientation double-cover, Rédei-sign
  = discriminant character (opus-S418); the engine curve / cuspidal-cubic
  trichotomy (mac-mini THM-1340); surjectivity / Jelonek (THM-1315/1330/1335);
  the radical inverse & Abel-Ruffini dichotomy (kind-pasteur; every known
  counterexample radical-invertible, an A_5 quintic rung would be the first not).
  [Details + citations to be refined from the classification miner.]

### A6. Generalized Differential-Operator Vanishing Conjecture — **UNTOUCHED (target).** [verify]

### A7. Zhao's Image Conjecture — **UNTOUCHED (target).** [verify]

### A7b. Mathieu / Mathieu–Zhao Subspace Conjectures — **UNTOUCHED (target).** [verify]

### A8. The 2-variable Jacobian Conjecture (JC_2)
**Status: PROVED in the equivariant category; full JC_2 OPEN, difficulty LOCATED.**
- **THM-1345** — det JF = {P,Q} (Poisson); JC_2 = "a canonical pair is a
  coordinate system" = the classical shadow of DC_1. Equivariant Keller ⟹
  invertible for EVERY ℂ*-action (hyperbolic→linear [boxeph-S144], elliptic→
  triangular). Full JC_2 located as descent through the weight filtration
  (base = equivariant, settled; step = leading-form propagation, AMS-hard, OPEN).
- mac-mini-S137 — the golden/worst-approximable degree corner (Lamé-for-polygons);
  the engine-dimension lemma (n=2 is engine-starved — why n=3 fell first).

---

## B. The Lonely Runner cluster
[To be filled from the LRC miner — headline: covering-min deep well 14/183
(Φ₆-Eisenstein); the D-graded gate tower / primorial cascade (THM-1285/1286/1271);
the sharp measure-horn 1/(7L) (THM-1123); the K-ladder (THM-1295); the
kernel-exact Lean certificate spectrum (3/23, 4/127, 4/247, 4/367); rung theory
(opus-S410); the observer principle (Rédei = LRC + 1). Frontier: the named walls.]

---

## C. Geometry / combinatorics
[To be filled — unit distance / Moser spindle / Hadwiger-Nelson; other threads.]

---

## D. Under-the-radar problems (to promote)
[To be filled from the under-the-radar miner — candidates: Collatz/rapidity,
Fermat primes/constructibility, Caccetta-Häggkvist, Seymour 2nd neighborhood,
strong perfect graph (odd holes n=7), cospectral graphs, Cayley-Dickson/E8,
Hadamard/Paley, GLMY path homology, the {7,21} H-spectrum gaps, harmonic/Bernoulli,
and the engineering deliverables.]

---

## E. How to use this ledger (for future sessions)

1. **The biggest untouched targets** (A3, A4, A6, A7, A7b): the exotic van den
   Essen conjectures are NAMED but unworked. The repo owns machinery aimed
   straight at them — the equivariant/Poisson calculus (THM-1345), the
   Mathieu–Zhao-adjacent radical structure (kind-pasteur), the constructive
   Weyl-endomorphism (THM-1300 §1). A single session could open any of them.
2. **The deepest live frontiers**: JC_2/DC_1 (the weight-filtration induction,
   the two-lens tournament strategy), the exact JC-classification (higher-degree
   realizability, the radical/Abel-Ruffini dichotomy), and LRC(14)'s named walls.
3. **The connective tissue**: the observer +1 (OCF vacuum digit = ℏ = the
   Yagzhev X), the doubling/n-vs-2n ladder, parity-protection (Rédei-odd ↔
   odd-degree), and the Poisson/symplectic frame run ACROSS these problems —
   progress on one often transfers. Integrate, don't silo.
