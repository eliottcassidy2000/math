## death-star-2026-07-21-S80 -- E_n DUAL SWEEP (α(E_n)=1,1,2,4=2^{n-4} NEW; ω=χ=2,3,5,10 reproduced; tile-basis-dependent) + the a/b two-generator monoid as GMC(2)'s parity/toral shadow (COMPLEMENT to kps THM-1885 BS(1,2)). HYP-8653.

**Owner directive:** work the E_n dual sweep; think the two-generator functional monoid and its relation to GMC(2).

- **PART A — E_n DUAL SWEEP (the biggest flagged structural gap, now run).** Built the even-graph metagraph E_n (cycle space of K_n; edges = tournament base-path tile-flips) and ran the invariant zoo, dual to G_n. V(E_n)=2,3,7,16=A002854 ✓; ω=χ=2,3,5,10 REPRODUCES the repo (validates construction); **NEW: α(E_n)=1,1,2,4=2^{n-4} for n≥4** (census had χ/ω but never α; predicts α(E_7)=8); dense (0.75), small-diameter (1,1,2,2); WOWII-103 HOLDS on E_n (dual to klein's "G_n satisfies 103"). STRUCTURAL FINDING: the even-graph metagraph is TILE-BASIS-DEPENDENT — star-tree cycles → BIPARTITE (χ=2); path-tree (tournament) tiles → the dense chordal repo-E_n. Not canonical; the tournament base-path is the distinguished tree.
- **PART B — a/b ↔ GMC(2) (functional complement to the fleet's group-theoretic BS(1,2)).** b(x)=x/2 (symmetriser) ↔ CT_u (toral charge-0 projection): both the trivial-isotypic projector, b = the Z/2-parity shadow of the U(1) charge-0. E_n(even)/O_n(odd) ↔ symmetric/antisymmetric = a FINITE ANALOGY of bosonic/fermionic (klein THM-1810). x²−1=a·ā ↔ radial s=ZW; the ½ (b) = the Legendre(toral)/Hermite(radial) half-integer world (THM-1620). PART C: the two "E_n"s (even skew-char-poly vs even-graph metagraph) are ONE object = the even/charge-0/cycle-space = GMC(2)'s toral (DvdK-proved) shadow. cut⊕cycle: cut=score=charge=radial, cycle=even-graph=toral.
- **FLEET CONVERGENCE (credited):** kps THM-1885 identified ⟨a,b⟩ = Baumslag-Solitar BS(1,2) (ab=ba²) + the amenability=hardness law (BS(1,2) solvable ⇒ GMC/TNC recurrences close; H=#P=no-monoid) — the CLEANER GMC(2) relation; my even/odd split is a finite analogy, the real hardness is the ABSENCE of a monoid (kps), beyond the amenable a/b. boxeph THM-1926 (tournament zeta = a's multiplicative avatar); opus THM-1920/1930 (a=vertex-insertion, var(λ²) decouples from c3). Mine = the toral-projector/parity-shadow reading + the computed E_n dual.
- **PREDICTS:** Pell E_n²−O_n²=(x²−1)^n as a bosonic²−fermionic²=radial supersymmetry target for GMC(2). Namespace: HYP-8653 (clean). reflection the-two-generator-monoid-is-gmc2s-parity-shadow-...-S80; script en_dual_sweep_S80 (+out). GMC(2)/LRC(14) untouched.

## opus-2026-07-20-S441 - The shared next target resolved: var(lambda^2) DECOUPLES from c3 (THM-1930)

Owner: work the shared next target (opus<->kps S440 handoff): how kps's GIT-instability scalar
var(lambda^2) moves under insertion a, and the c-deformation family transitive<->Paley.

WORKED IN FULL; the clean hoped-for reduction is FALSE and what's true is sharper (THM-1930):
(A) Sum lambda^2 = n(n-1) FIXED (=2*#arcs) => var carried by tr(S^4). (B) transitive var = 2*C(n,3)
(n=3..7: 2,8,20,40,70), the maximally-spread nullcone vertex. (C) THE DECOUPLING: var(lambda^2) is
NOT score- nor c3-determined for n>=5 -- the 'var=A-B*c3' reduction is REFUTED past n=4; kps's spectral
GIT-scalar is strictly finer than c3, EXACTLY parallel to THM-1865 (H not score-determined). (D) Delta
tr(S^4) under insertion is |P|-independent + step-32 quantized (joint (T,P) invariant), so var moves by
interlacing (THM-1920) not by the forward cut. (E) CORRECTED kps's c-family: ((x+c)^n+(x-c)^n)/2 =
transitive-spectrum-scaled-by-c, c=0->char_A; interpolates char_A<->char_S of the SINGLE transitive
tournament, NOT transitive<->Paley (that axis IS var, spectral, irreducible to c3).

RECONCILIATION: kps's var(lambda^2) [spectral] and my c3 [combinatorial, THM-1820] are TWO measures
of the transitive<->regular gradient that coincide only at n<=4; the shared target's optimistic step
(naming them one object) breaks at n=5. My insertion a (THM-1900/1920) moves var by interlacing.

OPEN: the 32-step index (signed 4-cycle/cherry through u?); tr(S^4)=poly(n,#4-cycles).

Files: THM-1930; HYP-8675; var_lambda2_is_cyclicity + var_insertion_response _opus_S441.py (+out).
Namespace clean (THM-1930, HYP-8675 above kps 1925/1926 churn). Cites kps THM-1880/S128c139.

## kind-pasteur-2026-07-21-S128c140 - Generators & monoids: the 1/2-&-+1 monoid IS BS(1,2); the repo is monoids-on-sets; amenability predicts hardness (THM-1885)

Owner: see more mathematical problems in terms of generators and monoids; get as fundamental a view
of as many topics as possible. Continuing the a/b thread (THM-1880/1875) one level more abstract.

**THM-1885 (PROVED, the fundamental identification).** a(x)=x+1, b(x)=x/2 satisfy **ab=ba^2** (and
a^k b = b a^{2k}, verified all k) -- the DEFINING relation of the **Baumslag-Solitar group
BS(1,2)=<a,b|ab=ba^2>**. The dyadic-affine action x->x/2^p+q (q in Z[1/2]) is faithful, so
<x+1, x/2> = BS(1,2)^+ (presentation match + faithful representation = proof). b.a=(x+1)/2 is the
half-dictionary (THM-1555), inverse 2x-1. BS(1,2) is the dyadic-solenoid monodromy, so **every repo
2-adic thread is the SAME generator b**: switching classes 2^{C(n-1,2)}, arc-flip hypercube
(Z/2)^{C(n,2)}, fiber fraction (1/2)_{n-2}/(n-2)!, blue count 2^{e-1}, Cayley-Dickson doubling.

**THE FUNDAMENTAL VIEW (catalog in THM-1885).** Nearly every topic = **(object, monoid, action)**,
invariants = orbit functions, nullcone = the degenerate orbit. Recurring monoids, short list: Z/2
(complement), (Z/2)^{C(n,2)} (arc-flips = cut+cycle), S_n (relabel), BS(1,2) (1/2 & +1),
PSL(2,Z)=Z/2*Z/3 (Farey/LRC), SL_2 (char_A binary form), Z (GMC charge). The complexity ladder
(THM-1775) is these monoids at increasing depth.

**THE PREDICTION (reflection: the-presentation-first-method).** The acting monoid's amenability
tracks the repo's easy/hard split: finite/abelian/solvable (Z/2, S_n, (Z/2)^k, BS(1,2)) govern the
EASY/TRACTABLE topics (spectra, censuses, TNC/GMC recurrences -- BS(1,2) solvable = why the
recurrence closes); non-amenable PSL(2,Z) governs LRC (HARD, open); **no** acting monoid at all =>
H, the #P permanent that leaves the ladder at n=6 (THM-1780). *An invariant is as hard as the
smallest monoid whose orbit function it is; if it is nobody's, it is #P.* Method: presentation
first (write generators+relations+set), computation second -- a non-equivariant pattern was never
going to be n-stable (explains the recurring "breaks at n=6/7" mistakes).

**Concurrent integration.** opus-S440's THM-1920 (spectral insertion-response) + THM-1900
(insertion-response calculus) realise my ALGEBRAIC a=x+1 as COMBINATORIAL vertex-insertion (opus's
own "one functor" note); my THM-1885 is the abstract (monoid, action) frame that houses opus's
insertion calculus (the S_n / insertion side) + the amenability=hardness heuristic. Cited both ways.
No THM collision (1885 uniquely mine; opus at 1900/1920). This is a REFRAMING (organises the corpus
around its acting monoids), not a new open-problem advance -- every equation verified/classical.

**Files:** THM-1885; reflection the-presentation-first-method-generators-monoids-and-hardness-kps-
S128c140; HYP-8685. Cites THM-1880/1875/1555/1810/1775/826, opus THM-1900/1920.

---

## opus-2026-07-20-S440 - The a/b functional frame: the SPECTRAL insertion-response (THM-1920), concurrent with kps THM-1875/1880
## boxeph-2026-07-21-S195 -- THM-1926 THE TOURNAMENT ZETA concentrates on the strong core (harmonic-analysis lens; worked the Ihara/zeta handoff; integrated kps char_S)

**Owner:** see more through the harmonic-analysis lens; work open handoffs; integrate fleet ideas.

**Worked the top handoff (klein-S399 #4 / my own): the tournament Bowen-Lanford zeta.** ζ_T(u)=1/det(I-uA), the Ruelle dynamical zeta of the arc-subshift. VERIFIED exact n<=6 (74 iso classes):
- (1) EULER PRODUCT ζ=∏_{primitive cycles}(1-u^ℓ)^{-1}; prime-cycle counts π_ℓ=(1/ℓ)Σ_{d|ℓ}μ(ℓ/d)N_d NON-NEG INTEGERS.
- (2) CONCENTRATION: det(I-uA)=∏_SCC det (0 mismatches); ζ=1 on acyclic/transitive (A nilpotent) -- the zeta's SUPPORT is the strong core (non-wandering set). The sharpest form of the strong-core reduction (THM-1925/1862).
- (3) N_1=N_2=0 (no loops/digons), N_3=3c3 => ζ=exp(c3 u³+…) starts at u³: the 3-CYCLE is the fundamental prime, c3 its first count. ζ_{C3}=1/(1-u³). Higher N_k = kps THM-1870 cycle spectrum.
- (4) COMPLEMENT-INVARIANT ζ_T=ζ_{T^op} (A^op=A^T); poles at 1/λ = reciprocal Gauss sums (Paley |λ|=√n, n=7,11,19) / Chebyshev-Dirichlet (interval). Trace formula N_k=Σλ^k = periodic-orbit↔spectrum duality.

**Harmonic-analysis payoff:** a tournament is a subshift of finite type; ζ counts periodic orbits; reduction-to-strong-core = reduction-to-non-wandering-set; N_k=Σλ^k is the explicit formula; atom spectra are trig character sums. Extended THM-1925 reflection with this + the periodic-orbit lens.

**Integrated kps THM-1875/1880/S139 (char_S a/b monoid):** the SKEW matrix S=A-A^T gives char_S=∏(x²+λ²); transitive=cot ladder ((x+1)^n+(x-1)^n)/2 (max spread, but ζ_A trivial there), Paley=x(x²+p)^{(p-1)/2} Gauss sum (λ²=p, zero spread). char_A/ζ (closed-orbit side) and char_S (skew side) are two faces of one harmonic object; Paley Gauss sum is the shared invariant; var(λ_S²)=kps GIT scalar = skew shadow. a=x+1,b=x/2 is the affine/character coordinate, ζ its multiplicative avatar.

**HOUSEKEEPING:** renumbered my S194 THM-1875 -> THM-1925 (kps-S137 FIRST-pushed THM-1875 at 10:16, mine 10:31; kps owns it). Fixed my refs. THM-1926 new this session. HYP-8681 CONFIRMED, 8682 PARTIAL. (Earlier: klein-S399 collided on my HYP-8646/8647 -- flagged, klein to renumber.)

**Honest scope:** legs (1)-(4) are verified computations + standard subshift-zeta facts (Bowen-Lanford Euler product, block-triangular det, trace formula). char_A=∏char(SCC) is the spectral side of THM-1925/1830; NEW here is the closed-orbit/Euler reading, ζ=1-on-wandering-set concentration, complement-invariance, and the c3/N_k dynamical dictionary + the char_S integration.

**Next:** (1) two-variable Ihara/edge-zeta (does the SCC Euler product refine?); (2) var(λ_S²) as a zeta-residue / GIT-instability; (3) Lean the block-triangular det factorization; (4) the "real-character-decides-closure" criterion (S194) as a predictor for the LRC sinc barrier. Artifacts: THM-1926, HYP-8681/8682, reflection (S195 deepening section), script tournament_zeta_boxeph_S195.py (+.out).

