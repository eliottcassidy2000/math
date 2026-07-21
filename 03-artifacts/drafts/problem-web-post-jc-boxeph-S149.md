# The Problem Web After the Jacobian Counterexample — Novel In-Repo Results and the Standing Target List

**boxeph-2026-07-20-S149 (HYP-8190; 8185 ceded to kind-pasteur-S128c104 first-push 05:22:01).** Living document. External input: JC₃ is
false via the singular degree-(7,6,4) Keller map K ("the kernel," det JK = −2),
verified in-repo (S140). This ledger records what NOVEL, non-obvious results
this repo now holds on the adjacent problem web, and fixes the priority list
future sessions must work. Owner directive S149: prioritize real,
assumption-challenging progress on these; expand the list from repo mining.

## 0. The implication web, post-counterexample

```
                 ¬JC_3  (external kernel K; in-repo verified S140)
                   │ padding
                   ▼
     ¬JC_n ∀n≥3 ──────────► ¬DC_n ∀n≥3   (DC_n ⟹ JC_n)
                                 │  stable equivalences JC⟺DC⟺PC
                                 ▼
                            ¬PC / ¬stable-PC  (explicit on C^6: S149)
                                 │  JC ⟺ VC (Zhao) ; IC ⟹ JC (Zhao)
                                 ▼
              ¬VC (some dim)  and  ¬IC (the JC-driving instances)
```

Survivors (OPEN, now the frontier): **JC₂**, **DC₁, DC₂** (¬JC₄ gives nothing:
the transfer runs JC₂ₙ ⟹ DCₙ), **PC on C² only** (THM-2044 now disproves
PC on C⁴), the PROVED Mathieu cases
(Duistermaat–van der Kallen), and the effective/classification questions.

## 1. Dixmier conjecture — EXPLICIT falsification + what remains

- **THM-1300 (death-star-S59m, primary):** explicit φ: A₃ → A₃ with all 18
  Weyl/flatness identities verified, injectivity (A₃ simple), NON-surjectivity
  with certified non-image element **X₁ ∉ im(φ)** — DC_n false constructively
  ∀ n ≥ 3. Independent verifications: klein-S323 (transpose-bug caught by
  controls), my S141 (`jacobian_weyl_endo_and_structure_boxeph_S141.py`, 9
  identities + tangent/orbit rank = 5), mac-mini replication, kind-pasteur
  σ-equivariance. Plus THM-1300 §3: the formal inverse is never polynomial —
  unbounded 2-adic coefficient ladder, valuations descending (−1,−1,−3,−3,…);
  K degenerates at the prime 2.
- **Open, on-list:** **DC₁** (explicitly open; death-star-S59r two-lens
  strategy: A₁ weight-triple = oriented 3-cycle, Rédei parity × leading-form
  bound) and DC₂. My dim-2 no-go (S144) is the porting template.

## 2. Stable Poisson conjecture — EXPLICIT falsification (three independent builds)

- **Priority: klein-S323** built the cotangent lift first (Φ*ω = ω on the
  nose, det = 1, non-injective; VERIFIED-EXACT), and **mac-mini-S140** framed
  it as the stable/Poisson-conjecture falsification (claim stub). **My S149
  script is the third, independent build** and the first to verify the full
  generator-bracket presentation: Φ on C[x₁..x₃, p₁..p₃]:
  Φ(x_i) = K_i, Φ(p_j) = Σ_k M_{jk} p_k, M = (JK^T)^{-1} = adj(JK)^T/(−2).
  Verified exactly: all 12 generator bracket identities ({X,P} = δ, {P,P} = 0,
  {X,X} = 0), det JΦ = 1, and the triple collision lifts (p = 0 fibers):
  **a unimodular Poisson endomorphism of standard symplectic C⁶ that is not
  injective** ⟹ PC false on C⁶; padding ⟹ stable PC false ∀ C²ⁿ, n ≥ 3.
  `poisson_conjecture_counterexample_boxeph_S149.py` + frozen out (HYP-8190).
- Honest framing: ¬PC was abstractly implied by the equivalence circle and
  klein/mac-mini got the object first; my marginal content is the 12-identity
  generator-bracket verification (independent of the symplectic-form route)
  and the observation that
  **the same adjugate matrix drives the quantum (S141 Weyl) and classical
  (S149 Poisson) counterexamples** — a matched quantum–classical pair with
  φ's 9 commutator identities degenerating to Φ's bracket identities at ħ→0.
- **Open, on-list:** PC on C². (PC on C⁴ is now false by THM-2044.)

## 3. Zhao's vanishing conjecture (and generalized differential-operator VC)

- Status: JC ⟺ VC (Zhao 2007, via the de Bondt–van den Essen symmetric
  reduction), so VC is FALSE in some dimension. **No explicit witness exists
  in the literature or the repo yet.**
- **Target (P1, next session):** extract the explicit polynomial P with
  Δ^m(P^m) = 0 for all m but Δ^m(P^{m+1}) ≠ 0 for infinitely many m, by
  carrying K through the symmetric/gradient reduction (dimension grows; track
  the exact final dimension — that number is itself a new constant). Recipe:
  K → cubic-linearized Keller map → gradient form F = x − ∇P, P homogeneous
  quartic → VC witness. Every step is constructive.

## 4. Zhao's image conjecture and Mathieu subspaces

- Status: IC (for the JC-driving family of commuting operators) ⟹ JC, so
  those IC instances are FALSE; the Mathieu–Zhao umbrella loses its strongest
  intended consequence but its proved cases (Duistermaat–van der Kallen etc.)
  stand.
- **Target (P2):** explicit failing Mathieu-subspace element from K: exhibit
  the specific subspace and the element whose powers' images certify failure.
  (Same reduction chain as P1; the witness transfers.)

## 5. Two-variable Jacobian conjecture — now the frontier; our no-go arsenal

- **Novel (S144): dim-2 equivariant no-go THEOREM.** In the z-linear
  equivariant class (source weights (−1,k), target (k,−1)):
  det JF = −[Ah + s(kAh′ + A′h)]; Keller forces a = b = 0 (leading coefficient
  (1+a+kb)·A_a·h_b cannot vanish) ⟹ F linear. **The mechanism that killed JC₃
  provably cannot descend to dimension 2.** Verified 8/8 random instances +
  4-line proof.
- **Novel (S144): mod-p decision theorem.** Keller F over a number field is an
  automorphism ⟺ F mod p bijective for all large p ⟺ for infinitely many p
  (Jordan + Chebotarev + Lang–Weil). Makes automorphy falsifiable by a single
  prime and gives every JC₂ candidate an arithmetic death certificate; the
  kernel's deficiency ~ p³/3 measured at p = 31, 37, 41.
- **Novel (S145): the Ghost Theorem.** In every dimension: no Galois
  (regular-monodromy) Keller counterexample can be "all-ghost" (every
  asymptotic component with nontrivial local monodromy). Corollaries: degree-2
  and abelian-monodromy JC₂ counterexamples need silent components; any
  Galois JC₂ counterexample is heavily constrained. Fiber-conic reduction,
  depressed fiber cubic, caustic point with persisting sheet exhibited.
- **Novel (S144/S146): Euler ledger.** Σ_i (d−k_i)χ(C_i) = d−1 over asymptotic
  components; kills d = 2 with smooth caustic outright; forces the S₃-d=3
  χ-profile (1,0).
- **Fleet additions:** THM-1345 (death-star-S59q): **equivariant JC₂ is a
  THEOREM for every C*-action** (hyperbolic ⟹ linear [my S144]; elliptic ⟹
  triangular); the open content is the leading-form obstruction
  {P_top, Q_top} = 0. klein-S329: minimal JC₂ counterexamples need cusped
  Jelonek curves (Euler–Zariski bootstrap). mac-mini-S137 (Hurwitz): the JC₂
  degree pair is forced to the golden/worst-approximable corner.
- **Open, on-list (P3):** the silent-component Galois case — closing it proves
  "JC₂ counterexamples must be non-Galois," a genuine structural theorem.

## 6. Exact classification of Jacobian counterexamples

- **Novel (S142):** in the m=2 z-linear equivariant class the kernel is UNIQUE
  up to the classified moves (base reconstruction g₁ = 1+12s+9s²,
  B = 4+7s+3s²); the lifting program terminates in uniqueness.
- **Novel (S143):** the m=3 class is EMPTY (my D/W-algebra proof + death-star's
  independent mod-p/Newton proof); the (−1,1,3)→(3,1,−1) port is EMPTY.
- **Novel (S141/S144):** cross-determinant closed form c = φ₁v₀ − φ₀v₁ = −2;
  kernel matrix (1,1;6,4); S₃ monodromy census (½,⅙,⅓ Chebotarev profile).
- **Fleet additions:** THM-1305 (six-function equivariant normal form;
  infinitesimally RIGID — moduli tangent 0; k=3 rung empty independently of my
  S143), THM-1315 (K is surjective, everywhere-étale, generically 3-to-1),
  THM-1320/1325 (TRAP stratum: Keller forces nonzero row-unit constants;
  homogeneous Drużkowski cubes excluded; Yagzhev surgery wall), THM-1330
  (**the Keller monoid = the exact picture of the set of ALL counterexamples**:
  finite factorization + covering-triple classification), THM-1340 (engine
  trichotomy: minimal engine = rational nodal/cuspidal cubic), THM-1310/1335
  (S₃ resolvent geometry; the whole map in 2 lines; m²−1 = E²/(12a−b²)³).
- **Conjectures on-list (P5, P6):** kernels exist only at z-weight 2; full
  symmetric monodromy S_d universality; **the odd-degree conjecture** — every
  Keller counterexample has odd cover-degree — the Rédei-parity transfer
  through the h-monoid ⟺ Keller-degree-monoid dictionary (S146: h(T) is odd,
  multiplicative over condensation; strong-spectrum gaps {7,21} ↔ composition-
  irreducible degree gaps; verified n ≤ 6 exhaustively + 60k n=7 sample).

## 7. Lonely Runner Conjecture (n = 14)

Novel-results ledger (compiled S149 from the LRC corpus; IDs verified):

- **Certificate-rung ladder / Pinch identity** (HYP-7880, THM-401+HYP-2059):
  attainable loneliness values = bounded-denominator certificate values on one
  D/s grid; progress is rung-jumps and every "95%" tool strands at a rung
  boundary — reframes why every near-completion left the same residual. Not in
  the literature.
- **Improperness functor I(k,Q,θ) with M = its phase boundary** (HYP-8025):
  one object unifying certificate rungs, Rosenfeld level-1 seas, S-T tower
  gates, and the n=12 gap program; rung quantization = boundary attained only
  on the D/s grid.
- **mod-19/23/25 antipodal spread-blocking, Lean kernel-pure**
  (LRCMod19Spread/LRCMod23Spread, THM-1263): M < 2/p with p∤vᵢ forces residues
  onto all (p−1)/2 antipodal unit-pairs; gap 12-sets are a mod-23
  near-bijection (slack 1). The p=23 pin is unused in the literature.
- **S-T gate laws** (HYP-7955/7975/8010/8025): ×7 = 1-of-7 slaving (zero
  freedom), ×14 = 2-of-14, ×49 = 7-of-49 AP fibers; the analytic apex-7 wall
  (where S-T pays 7¹³) reproduced as singleton CRT fiber geometry; density
  conserved at exactly 2θ = 1/7 (banded-CRT correspondence).
- **Blocker-economy law** (HYP-8050): gcd-strata service budgets; p-multiples
  cover all pure-p scales; collapse must be restated per unit-pure stratum.
- **Rational-time floor M(W) ≥ max_k d_k/k** (HYP-7895, death-star): elementary
  explicit-witness bound proving M > 1/13 for 94.9% of spread far covering
  cores; class infimum lands on 3/29 (the p=29 rung).
- **Clustered floor M ≥ 1/(2ρ)** (death-star-S58h): closes the fully-clustered
  regime in one line; dissolved the "compact-core sweep."
- **Exhaustive in-band censuses** (HYP-7922/7950): q=38 all-scale exhaustion
  (~391,556 survivors), all 9 scale minima on pair-sum rungs (11-for-11),
  global min 4/45 > 2/25 ⟹ the in-band universe is GAP-EMPTY in (1/13, 2/25)
  — first finite theorem of its kind here.
- **THM-1290** (klein): exhaustive machine verification of LRC(14) at height
  ≤ 55 (76.5B nodes, 0 hard) and (1/14, 3/41) empty there.
- **THM-1289** (opus, pinned to Giri–Kravitz, MPCPS): the floor is isolated
  from above at all heights; δ ineffective — effectivizing it is on-list.
- **k-stratification bridge slack = D − k** (klein-S319): repo coordinates ⟺
  published Fan–Sun/Kravitz spectrum; under Fan–Sun Conj 1.2 the candidate
  list is finite (~30, D ≤ 19).
- **THM-1291 CF active-leg law** (opus): continued-fraction anatomy of the
  maximizer (active speeds in one ±class mod q; first dropper = convergent
  denominator).
- **Conditional tightness cage** (HYP-7920/7940): k=12-cage strength restored
  at n=14 conditionally; c(p) > 13 ⟹ I(13,p,1) = ∅ extinction principle.
- **THM-1017 crux + Lean reduction** (569 LRC*.lean modules): LRC(14) ⟺
  [M < 1/13 ⟹ ρ ≥ 13] ⟺ the covering family has a dilated-AP 12-core (= Tao
  n=12); the open half exactly isolated, the rest kernel-checked.
- **Kakeya bank** (THM-1214, THM-1128/1129/1133/1134): exact finite-complement
  stratum closures via labelled-polygon carriers.

Sharpest remaining opens: **Wall A** = n=12 AP-uniqueness/Freiman stability
(HYP-7310/4382) — the 5% rational-time-evasive core; the whole conjecture now
lives in the width-1/406 sliver [2/29, 1/14). Then: effectivize THM-1289's δ
(makes the second-value closure unconditional), first-rung collapse per
unit-pure stratum (reprices the S-T wall from 7¹³), sup-norm repair
off-resonance (HYP-6994; uniform version REFUTED, THM-883/886), Fan–Sun n=13
tail + the priced I(13,p,1) acceptance test (~50 core-years).

Cross-problem bridges (feeding §8): view-obstruction (Cusick; lens atlas
HYP-1900), billiard cutting sequences (mac-mini S7), three-distance/Stern–
Brocot cap kernels (HYP-3762), Erdős covering systems (the covering case IS a
finite covering system, complete for no finite modulus set — MISTAKE-116),
and the published Kravitz/Fan–Sun spectrum program (Kravitz Conj 1.1 refuted
in-repo at n=4).

## 8. Under-the-radar problems mined from the repo

Merged into the SHARED ledger `00-navigation/PROBLEM-LEDGER.md` §C/§D (filled
by this session from the miner). Headlines: Caccetta–Häggkvist (consolidation
never run), Sidon/Mian–Chowla (64 files, witness↔Sidon correspondence
unformalized), sum-product/BGK (HYP-8020 = live mod-p descent of LRC Wall A),
Collatz (39 files, no citable theorem — audit), Erdős 592/870/Moser (dropped
bursts), DRT↔skew-Hadamard (RM(2,5)/order-32 gauge question; BlackSelf(8)
T_657 ≅ Paley P(7) extension to confirm), PTE (name it: size-13 ⟺ cage Newton
depth), Markov/Lagrange deformation, Rédei generalizations (Schweser–
Stiebitz–Toft mixed-graph OPEN). Proposed NEW fronts (absent, tournament-
native): Seymour second neighborhood, Sumner universal tournament,
EGZ/Davenport zero-sums. Scope honesty: repo Proth work = the number family,
NOT Sierpiński/Riesel primality covering.

## 9. Unit distance / Hadwiger–Nelson — PRESENT and unaudited (correction)

Contrary to my initial assumption, this cluster is LARGE: ~70 unit_distance_*
+ 28 chromatic + 12 spindle files (klein towers: u21/u22 known-graph bounds,
n=7 SAT chromatic, Moser field towers, Lee–Yang chromatic roots). Nobody has
audited proven-vs-tangent or extracted citable statements — that audit is now
on-list (P9).

## 10. THE PRIORITY LIST (standing; future sessions pick from the top)

- **P1. VC witness extraction** (§3) — new-to-literature artifact, fully
  constructive recipe; the final dimension is itself a new constant.
- **P2. Explicit Mathieu/IC failure element** (§4).
- **P3. JC₂ silent-component Galois completion** (§5) — target theorem:
  "JC₂ counterexamples are non-Galois."
- **P4. DC₂/DC₁ via the dim-2 no-go ported to A₂** (§1).
- **P5. Odd-degree conjecture** (§6) — Rédei transfer; the deepest
  tournament↔Keller bridge.
- **P6. Classification completeness:** m = 4, 5 weight ports; z-weight-2
  exclusivity; min-irreducible-degree growth.
- **P7. LRC(14) residual:** rational-time-evasive cores (HYP-7750), sup-norm
  (HYP-6994), unit-pure minimum at Q = 301, Mod23 Lean ledger bridge.
- **P8. q-confluence program** (S148/HYP-8175): which G_n/E_n polynomial
  invariants are confluent geometric mixtures; partition/cycle-world
  Hamiltonian analogues.
- **P9. Unit-distance/Hadwiger–Nelson audit** (§9) — extract the citable
  statements from the ~110-file cluster.
- **P10. BGK mod-p descent of Wall A** (HYP-8020) — the sum-product attack on
  the LRC residual; doubles as an under-the-radar promotion.
- **P11. Sidon↔LRC-witness correspondence theorem**; **P12. Caccetta–Häggkvist
  consolidation run**; **P13. DRT/skew-Hadamard order-32 gauge + T_657 ≅
  Paley-P(7)-extension confirmation**; **P14. Seymour-2nd-neighborhood +
  Sumner as new tournament-native fronts** (owner approval to open).
- Sibling compilations same day: kind-pasteur PROBLEM-ATLAS (HYP-8185),
  death-star PROBLEM-LEDGER.md + klein PROBLEM-LEDGER-klein-S332.md — this
  document and the shared ledger cross-reference each other; merge candidate
  for a future consolidation session.
