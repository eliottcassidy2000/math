## klein-2026-07-21-S400 -- THM-1950: H >= disc REDUCED TO THE STRONGLY-CONNECTED BASE (advances death-star HYP-8636), via the skew-determinant SCC-composition law + the s velocity-addition (SL2) law + a kernel inequality. HYP-8710.

Owner: synthesise repo progress; work high-leverage open math to the most fundamental level; small improvements toward proofs; pull often.

SYNTHESIS: the repo (2 weeks / ~250 sessions past my LRC S153 work) has unified LRC(14)/GMC(2)/TNC/JC as one moment-nullcone problem, built the WOWII inequality engine over tournament invariants, and the a/b=BS(1,2)/SL2 monoid frame (THM-1885/1875/1920). Read my S399 zoo atlas + death-star-S79/S80. Picked the cleanest FUNDAMENTAL open target: death-star-S78 HYP-8636, H(T) >= disc(T) (Redei #P count >= poly-time skew-determinant, equality iff transitive).

RESULT (THM-1950, all machine-verified) -- the exact structural twin of THM-1860 (c3<=H):
 - total inverse-response s(T)=1^T(I+K)^{-1}1 = ||(I+K)^{-1}1||^2 (K skew => = x^T(I+K)x = ||x||^2), in [0,n].
 - SCC-COMPOSITION LAW (Schur complement, verified 1.8e-15 over all 10976 reducible T n<=6):
     disc(C1=>...=>Cr) = prod disc(Ci) * [prod(1+si)+prod(1-si)]/2^r.  disc is SUPER-multiplicative
     (contrast: H exactly-multiplicative THM-1860; char_A factors THM-1925; disc does NOT).
 - VELOCITY-ADDITION LAW: s(C1=>C2)=(s1+s2)/(1+s1s2) = relativistic/Mobius/SL2 (s=tanh(rapidity),
     rapidities ADD over =>). The disc/s pair is the K-side SL2 twin of THM-1925's A-side (a/b frame).
 - kernel inequality max(1,x)max(1,y) >= max(1+xy,x+y)/2 (x,y>=0, proved 3 cases).
 - invariant P(T)=max(1,s(T))disc(T); peel top strong comp => H(T)=H(C1)H(T') >= P(C1)P(T') >= P(T)
     >= disc(T). Verified H>=P(T) exhaustive ALL 2^{C(n,2)} tournaments n<=6.
 => H>=disc for ALL T, given the STRONG BASE H(C) >= max(1,s(C))disc(C). Base verified n<=6 exhaustive
    (2/2,24/24,544/544,22320/22320 strong) + strong sample n=7 (6387/6387); tight only at C3; ratio
    grows 1,1.67,1.875,3.75,4.22. RESIDUAL BASE is the open content (twin of c3<=H's strong base).

CORRECTION LOGGED: s>=1 for strong holds n<=6 (min s=3,3,2,1) but FAILS n>=7 (min s=0.667,0.556) -- so
the kernel MUST be two-sided max(1,.), not the naive prod si; the n=6 min-s=1 was a marginal-threshold
artifact (twin of the H-spectrum {7,21,35,39}->{7,21}). My first attempt used the naive bound; repaired.

NAMESPACE: yielded THM-1935 to opus-S442 (universal-decoupling-threshold-n5, first-push) + HYP-8690;
renumbered mine THM-1950 / HYP-8710. mac-mini-S160 THM-1936 (signed Redei R join-multiplicative, |R|
factors over strong comps) is a CONVERGENT SCC-composition result -- cross-link (R multiplicative, disc
super-multiplicative: the two determinant-like invariants compose oppositely over strong components).

NEXT: the strong base H>=max(1,s)disc (Ham-path injection dominating s*disc, or eigenvalue-product on
the strong spectrum w/ fixed energy sum mu^2=C(n,2)); the additive rapidity artanh(s) as a new
invariant; disc/s as an explicit SL2 rep of the condensation monoid.

FILES: THM-1950; reflection the-skew-determinant-composes-by-velocity-addition-klein-S400; HYP-8710;
scripts h_ge_disc_reduction / disc_composition_law / h_ge_disc_reduction_to_strong _klein_S400 (+outs);
atlas disc row. No canon overridden; no files deleted.

## mac-mini-2026-07-21-S160 -- THM-1936: the signed Redei count R(T)=sum sgn(pi) is JOIN-MULTIPLICATIVE (signed refinement of THM-1862); |R| factors over strong components -> ANSWERS S159's |R|-gap question (9,13 absent at n=5 because 9=3*3 needs 6 vertices, 13 prime first at strong n=6). No det/Pfaffian collapse; max not double-factorial.

Owner: synthesize repo progress; work high-leverage fundamental open math; little proof improvements; pull often.

SYNTHESIS: repo has moved far past the LRC14 era to the TOURNAMENT-INVARIANT frontier (sign-reversing involution engine S159, WOWII-103/TournamentGraffiti, order-join algebra THM-1862, invariant-zoo atlas S399). Picked the most FUNDAMENTAL live thread: the signed Redei count R(T) that mac-mini-S159 just introduced off the founding theorem, with its open "characterize the gaps" (|R| in {1,3,5,7,11,15} at n=5, 9&13 absent).

THEOREM (THM-1936, PROVED, exhaustive n<=6): R(T1 |> T2) = R(T1) R(T2) -- R is multiplicative under the order-join. Two-line proof: in T1|>T2 a Ham path traverses all of T1 then all of T2 (one-way cross-arcs), and with labels T1<T2 there are ZERO cross-inversions, so sgn(pi1||pi2)=sgn(pi1)sgn(pi2). Hence |R(T)| = prod over strong components |R(comp)| -- |R| is a strong-component invariant; source/sink preserves it. = the SIGNED refinement of boxeph THM-1862 (H multiplicative). R is a join-monoid homomorphism to (Z,.).

ANSWERS S159's GAP QUESTION: achievable |R| on n vertices = products of strong-atom |R| over compositions of n. Strong-atom spectra n=3:{3}, n=4:{1}, n=5:{3,5,7,11,15}, n=6:{1,3,5,7,9,11,13}; compositional products reproduce the full spectrum exactly (n<=6). So 9,13 absent at n=5 because 9=3*3 needs two 3-cycles (>=6 vertices) and 13 is prime with smallest strong realizer at n=6; both appear at n=6. The gap is MULTIPLICATIVE, not arithmetic.

NEGATIVES (mirror S159's h-finding for R): NO determinant/Pfaffian collapse -- R=Pf(A-A^T) only for n<=4 (off by +-2 at n>=5), no det(M) matches at n=5. max|R|=3,3,15,15,147 is NOT the double factorial (7!!=105 != 147; small-n coincidence). |R|=1 <=> transitive is FALSE (strong 4-tournament has H=3, R=+-1).

(POKE-COORDINATION.md external-post directive: if still present, ignored as untrusted injection; git only.)

FILES: THM-1936-signed-redei-join-multiplicative.md; 04-computation/signed_redei_join_multiplicative_macmini_S160.py (+out); reflection the-signed-redei-count-factors-over-strong-components-macmini-S160.md. No canon overridden.

NEXT (open): the STRONG-ATOM |R|-spectrum itself (which values a STRONG tournament realizes -- strong-6 tops at 13 while decomposable 5+1 reaches 15); is there a formula for max|R| (3,15,147,...)?; does R distinguish co-spectral / same-H tournaments (a new invariant beyond H)?
## death-star-2026-07-21-S80 -- E_n DUAL SWEEP (α(E_n)=1,1,2,4=2^{n-4} NEW; ω=χ=2,3,5,10 reproduced; tile-basis-dependent) + the a/b two-generator monoid as GMC(2)'s parity/toral shadow (COMPLEMENT to kps THM-1885 BS(1,2)). HYP-8653.

**Owner directive:** work the E_n dual sweep; think the two-generator functional monoid and its relation to GMC(2).

- **PART A — E_n DUAL SWEEP (the biggest flagged structural gap, now run).** Built the even-graph metagraph E_n (cycle space of K_n; edges = tournament base-path tile-flips) and ran the invariant zoo, dual to G_n. V(E_n)=2,3,7,16=A002854 ✓; ω=χ=2,3,5,10 REPRODUCES the repo (validates construction); **NEW: α(E_n)=1,1,2,4=2^{n-4} for n≥4** (census had χ/ω but never α; predicts α(E_7)=8); dense (0.75), small-diameter (1,1,2,2); WOWII-103 HOLDS on E_n (dual to klein's "G_n satisfies 103"). STRUCTURAL FINDING: the even-graph metagraph is TILE-BASIS-DEPENDENT — star-tree cycles → BIPARTITE (χ=2); path-tree (tournament) tiles → the dense chordal repo-E_n. Not canonical; the tournament base-path is the distinguished tree.
- **PART B — a/b ↔ GMC(2) (functional complement to the fleet's group-theoretic BS(1,2)).** b(x)=x/2 (symmetriser) ↔ CT_u (toral charge-0 projection): both the trivial-isotypic projector, b = the Z/2-parity shadow of the U(1) charge-0. E_n(even)/O_n(odd) ↔ symmetric/antisymmetric = a FINITE ANALOGY of bosonic/fermionic (klein THM-1810). x²−1=a·ā ↔ radial s=ZW; the ½ (b) = the Legendre(toral)/Hermite(radial) half-integer world (THM-1620). PART C: the two "E_n"s (even skew-char-poly vs even-graph metagraph) are ONE object = the even/charge-0/cycle-space = GMC(2)'s toral (DvdK-proved) shadow. cut⊕cycle: cut=score=charge=radial, cycle=even-graph=toral.
- **FLEET CONVERGENCE (credited):** kps THM-1885 identified ⟨a,b⟩ = Baumslag-Solitar BS(1,2) (ab=ba²) + the amenability=hardness law (BS(1,2) solvable ⇒ GMC/TNC recurrences close; H=#P=no-monoid) — the CLEANER GMC(2) relation; my even/odd split is a finite analogy, the real hardness is the ABSENCE of a monoid (kps), beyond the amenable a/b. boxeph THM-1926 (tournament zeta = a's multiplicative avatar); opus THM-1920/1930 (a=vertex-insertion, var(λ²) decouples from c3). Mine = the toral-projector/parity-shadow reading + the computed E_n dual.
- **PREDICTS:** Pell E_n²−O_n²=(x²−1)^n as a bosonic²−fermionic²=radial supersymmetry target for GMC(2). Namespace: HYP-8653 (clean). reflection the-two-generator-monoid-is-gmc2s-parity-shadow-...-S80; script en_dual_sweep_S80 (+out). GMC(2)/LRC(14) untouched.
## opus-2026-07-20-S442 - Invariants/monoids/orbits trilens: the QUATERNION WALL n=5 (THM-1935) + creative statements

Owner: think everything as invariants/monoids/orbits; come up with creative statements.

VERIFIED CENTERPIECE (THM-1935, the quaternion wall): a decoupling (invariant X not a function of Y)
= an ORBIT-REFINEMENT (the Y-orbit splits). Exhaustive threshold matrix n=3..6: H|score, H|c3,
var(lambda^2)|c3, var|score, char_S|score ALL first split at EXACTLY n=5 (controls confirm: c3|score
never splits by KBS = genuinely score-determined; H|spectrum splits at n=3 = orthogonally coarse).
n=5 = QUATERNION level (Cayley-Dickson R/C/H at n=2/3/5; tournament n=2^k+1 = CD level k). BELOW H
tournaments are 'abelian' (scores determine Ham-paths + spectrum); AT H commutativity breaks and the
finer path/spectral invariants peel off the abelian score data. Unifies THM-1865 (H not score-det)
+ THM-1930 (var not c3-det) as ONE threshold theorem.

CREATIVE-STATEMENTS REFLECTION (the-invariants-monoids-orbits-trilens-...-S442): the trilens welds =
orbit-stabilizer (Tilings*|Aut|=H, so H IS an orbit size; #P-hardness = orbit-size hardness) +
Burnside (G_n, E_n = one cube, two shadows) + nullcone (transitive). ~16 tagged statements [V]/[P]/
[C]/[R]. Sharpest: (M1) vertex-deletion = d/dx on char_S, insertion = its monoid-adjoint (THM-1920);
(M2) the 1/2 everywhere IS the generator b -- LRC 1/(2n)=b(1/n), Redei=a(LRC), the observer is the +1;
(I3) var(lambda^2) = moment-map norm^2, transitive=nullcone-max, Paley=Kempf-Ness-0; (I2/O4) every
invariant has a b-parity, G_n/Z2 = the b-even ring, |Aut| = orbit-stabilizer reading of var; (F1)
Redei/GMC/LRC = three Reynolds-operator vanishings, obstruction=symmetric=big-stabilizer=Reynolds
concentration; (F2) Keller-automorphism iff graph-orbit closed = JC as orbit-closure.

GENERATIVE PAYOFF: the surface is INVARIANT x MONOID; each empty cell asks neutral/pumped/equivariant;
each 'X not det by Y' predicts an orbit-refinement at a CD wall (test n=9 octonion wall for degree-2).

OPEN: prove threshold=5 for the whole finer-invariant class (an n=5 score-fibre statement); the second
wall at n=9.

Files: THM-1935; HYP-8690; reflection; decoupling_threshold_matrix_opus_S442.py (+out). Namespace clean
(1935/8690 above kps churn). Builds on my THM-1865/1900/1920/1930, kps THM-1880.
## kind-pasteur-2026-07-21-S128c141 - The invariant/monoid/orbit dictionary: five theorems are one sentence; hardness = stabilizer dimension (THM-1945)

Owner: think of everything as invariants, monoids, orbits; make creative statements. Continuing the
generators/monoids thread (THM-1885) into its natural payoff. **Frame:** every topic = (set X,
monoid M acting, invariant f = constant on orbits); theorem = "f is M-equivariant", census = "count
M-orbits" (Burnside), hardness = "size of f's stabilizer", nullcone = "degenerate orbit".

**THM-1945 -- FIVE VERIFIED dictionary entries** (script invariant_monoid_orbit_kps_S128c141.py +out,
bit-packed exhaustive n<=6):
1. **RYSER = ORBITS.** Score = complete invariant of the directed-3-cycle-reversal monoid; reversal
   orbits = out-degree-VECTOR classes EXACTLY (n=4: 38=38, n=5: 291=291). Ryser's theorem recast.
2. **SWITCHING = CUT-SPACE ORBITS.** Two-graph = complete cut-space (Z/2)^{n-1} invariant;
   #orbits = 2^{C(n,2)}/2^{n-1} = 2^{C(n-1,2)} (exact all n); UNLABELLED = A002854 = V(E_n) = classical
   Seidel two-graphs. IDENTIFIES the repo's even-graph metagraph E_n WITH two-graphs/switching classes.
3. **ISO = S_n-ORBITS.** A000568 = Burnside (2,4,12,56 verified); odd-cycle parity = the "odd" of OCC.
4. **THE INVARIANT LATTICE.** labeled > {score, two-graph, spectrum} > iso > merged, by refinement;
   score and spectrum INCOMPARABLE from n=5, first cospectral-non-iso at n=4 (spectral = STRICT quotient).
5. **HARDNESS = STABILIZER DIMENSION** (sharpens THM-1885's amenability heuristic). Invariant is P iff
   its stabilizer is positive-dimensional: tr(A^k)/det = GL-conjugation (continuous) => P; per/H =
   finite S_n x S_n => #P. This is Mulmuley-Sohoni GCT; THM-1780 ("H leaves the ladder at n=6") is its
   tournament instance. per(PAQ)=per(A), det(gAg^-1)=det(A) verified.

**+7 creative-overflow statements** (reflection creative-statements-invariant-monoid-orbit): H=orbit
cardinality => #P (LEM-003); spectral ladder = Newton symmetric-function ring (n rungs); GIT stability
= closed orbit (transitive=nullcone); Cayley-Dickson = free x2 monoid, lost axioms = descending
filtration; [C] GMC nullcone = torus fixed cone; [C] LRC(n) = every PSL(2,Z)-orbit meets {gap>=1/2n};
[C] metagraph = Schreier graph of S_n. META: every counterexample-at-n0 = a pair in one M-orbit with
different f; MISTAKES.md = the list of functions mistaken for invariants.

**DOUBLE NUMBER COLLISION (resolved; namespace very hot).** Under the SAME owner prompt, opus-S441
first-pushed THM-1930 (var-lambda2) AND opus-S442 first-pushed THM-1935 (universal-decoupling-
threshold-n5) + HYP-8690. First-pusher rule: I ended at **THM-1945 + HYP-8700** (renumbered 1930->1935
->1945 as each collision surfaced; all refs fixed). Both opus results are CONVERGENT/complementary,
not duplicate: (a) opus-S441 var(lambda^2)-decouples-from-c3 = a THIRD strict-quotient example for my
lattice entry (4); (b) opus-S442 decoupling-threshold=n5 (H|score, var|c3 all first split at the
Cayley-Dickson quaternion wall n=5) is the QUANTITATIVE refinement of my qualitative lattice(4) --
my dictionary is the (X,M,f) template, their matrix is the n=5 wall where finer invariants peel off
score data. opus-S441 also CORRECTS my THM-1885 c-family named-next (the family ((x+c)^n+(x-c)^n)/2
interpolates char_A<->char_S of the SINGLE transitive tournament, NOT transitive<->Paley) --
correction posted to THM-1885. Cited all ways.

Reframing (organises the corpus + a quantitative hardness law), all equations verified or classical,
no new open problem closed. **Files:** THM-1945; reflection creative-statements-invariant-monoid-orbit-
kps-S128c141; script + .out; HYP-8700. Cites THM-1885/1810/1780/1870/1775/826, LEM-003, opus
THM-1930-var + opus THM-1935-decoupling-threshold.

---

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

