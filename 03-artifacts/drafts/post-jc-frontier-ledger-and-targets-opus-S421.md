# The post-JC frontier: what we have actually proved, and where assumption-challenging progress is next

**Instance:** opus-2026-07-20-S421 (HYP-8185). Owner directive: given JC is disproved
(externally supplied singular counterexample F, dim 3, det = -2), compile the novel
non-obvious results the repo has made on the downstream problems, EXPAND the list by
mining the repo for under-the-radar problems, and set priorities for future sessions.

**This is the living target list. Future sessions: pick from Part IV, update Part I/II
when a status changes, append to Part III when a new problem surfaces.**

---

## Part I — The cascade map (what ¬JC3 does, with correct arrows)

Arrows we rely on (classical, cited in THM-1300/1345): DC_n => JC_n (van den Essen);
JC_2n => DC_n (Tsuchimoto; Belov-Kontsevich); counterexample padding ¬JC_n => ¬JC_m
for m >= n; stable JC <=> stable DC <=> stable Poisson (Belov-Kontsevich, Kontsevich);
JC(all n) <=> Zhao vanishing conjecture <=> Zhao image conjecture (van den Essen-
Wright-Zhao / de Bondt reductions); Mathieu-subspace framework generalizes IC (Zhao).

| problem | post-disproof status | repo's contribution |
|---|---|---|
| JC_n, n >= 3 | FALSE | independent exact verification, torus anatomy, normal form (THM-1300/1305; my S413/S414) |
| JC_2 | OPEN (untouched by F) | proved F does NOT descend (quotient contraction, S414); JC2 = THEOREM in C*-equivariant category (THM-1345 + boxeph-S144) |
| Dixmier DC_n, n >= 3 | FALSE, constructively | EXPLICIT phi: A3 -> A3, 18 Weyl identities verified, injective non-surjective (THM-1300) |
| DC_1, DC_2 | OPEN (survive) | eigen-graded DC1 = the shadow of equivariant-JC2 (THM-1345): a NEW route INTO DC1 |
| stable Dixmier | FALSE | via Belov-Kontsevich equivalence from ¬stable-JC (S413 cascade; THM-1300 scope note) |
| stable Poisson (Kontsevich) | FALSE | same equivalence chain; PLUS the positive fragment: det JF = {P,Q} weight-additive bracket machinery (THM-1345) |
| Zhao vanishing conjecture (universal form) | FALSE (no explicit witness yet) | NOT YET TRANSPORTED — Target 1 below |
| Zhao image conjecture / Mathieu subspaces (universal forms tied to JC) | FALSE in JC-equivalent range; instance-level landscape now WIDE OPEN | NOT YET TOUCHED — Target 1 |
| classifying JC counterexamples | new subject | local sporadicity (S414b: tangent 5 = 1 trivial + 4 obstructed); Keller monoid: counterexamples = deg>=3 IDEAL (THM-1330); Druzkowski/trap stratum det factorization (THM-1320); engine trichotomy (THM-1340); weight-type machine template (S414) |
| Lonely Runner (14) | OPEN, floor isolated | THM-1288 (published S-T Conj 7.1 REFUTED), THM-1289 (G-K floor isolation), THM-1291 (CF active-leg law), THM-1292 (ghost packing), THM-1269 (D = M*s), covering-core 95% closure (death-star S58i), deep well 183 = Phi6(14) = |PG(2,13)| (S420) |
| unit distance problem | not yet worked | analogy-card level only (LRC technique indexes: "unit-distance cyclotomic norm", lowest retention) — honest: NO substantive repo result. Target 6 |

## Part II — The novel, non-obvious results (the "expand this")

**II.1 Dixmier, constructively (THM-1300, death-star).** Not just the cascade arrow:
the explicit endomorphism phi(X_i) = F_i, phi(D_j) = sum_k B_jk D_k with
B = (JF^T)^(-1), all 18 flatness/Weyl identities verified symbolically, properness of
the self-embedding via a module one-liner. Bonus non-obvious: the formal inverse of F
has an UNBOUNDED DYADIC coefficient ladder — a 2-adic obstruction fingerprint that is
new data about WHY no polynomial inverse exists (the det's -2 read p-adically).

**II.2 JC2 rescued equivariantly + the Poisson bridge (THM-1345 + boxeph-S144).**
Equivariant Keller => invertible for EVERY C*-action on C^2 (hyperbolic => linear;
elliptic => triangular). The reframing det JF = {P,Q} with the C*-action as the
Hamiltonian flow of xy places equivariant-JC2 as the classical shadow of eigen-graded
DC1 — the first structural bridge in the repo between the surviving bottom rungs
(JC2, DC1, DC2). The leading-form obstruction {P_top, Q_top} = 0 locates full JC2's
hardness in weight-filtration descent (named, not hand-waved).

**II.3 Classification as geometry of a monoid ideal (THM-1330 + S414b + THM-1320/1340).**
The counterexample SET is a two-sided ideal (deg >= 3) in the Keller composition monoid;
F is locally rigid modulo reparametrization within its weight type (deformation theory:
16 equations, rank 10, kernel 5 = 1 integrable-trivial + 4 second-order-obstructed);
the Druzkowski/trap stratum carries a det factorization law; the engine trichotomy
(THM-1340) sorts known engines. Together: "sporadic or family?" has an evidence-grade
answer — SPORADIC MODULO EQUIVALENCE, with families = reparametrization orbits, and
the machine for NEW sporadics is the weight-type template (pick torus weights, write
the invariant normal form, force det const, make T a nonunit).

**II.4 The fiber geometry as tournament theory (THM-1310/1335 + my S417/S418).**
Fiber cubic with disc = -4Q^2 L, S3 resolvent; master identity L(F)x^3+(4-3F2F3)x-2F3
= 0 (Chebyshev trisection); Jelonek set = cusp; sheet decomposition of {w=0} into a
tame plane automorphism + exotic double sheet branched on the conic v^2 = 16u (= the
w=0 slice of {L=0}); and the deck obstruction: Aut(C(x,y,z)/C(F)) = 1, the involution
iota lives exactly on the quadratic resolvent C(F)(sqrt(-L)) = the ORIENTATION DOUBLE
COVER, Redei sign = discriminant character, generic fiber = cyclic 3-tournament. The
repo's parity-in-tournaments core runs THROUGH the JC counterexample.

**II.5 LRC(14) (the theorem run).** Refutation of a published conjecture (S-T 7.1,
THM-1288: divisor-aligned clusters break translation-blind d-grids); the floor
isolated from above (THM-1289, including the withdrawn-v3 forensics); the CF active-leg
law (THM-1291: u* is always a convergent denominator); ghost packing (THM-1292:
M <= K/(K(m+1)+1), Phi6/14 = 13+1/14 as packing count); D = M*s (THM-1269); the
bilinear m*b+-1 rim unification (S419/S420) with 183 = 2*T_13+1 = |PG(2,13)|.

**II.6 The Proth/shear layer (S419/S420 + kp HYP-8170 + THM-1355).** Faulhaber
identification of the owner triangle; sums/products duality (sums = Eulerian
one-descent counts; pure products = m!*2^C(m,2) = ordered tournaments = the repo's
tiling model); Fermat-collision law n(2^d-1) = d unique; the 2^(1/s) shear spectrum
with recurrence normal forms (x^s-2)(x-1)^2(1+...+x^(s-1)); OEIS-new sequences.

## Part III — The expanded list (under-the-radar problems mined from the repo)

Each entry: problem — repo hook — why it belongs.

1. **OCF / CONJ-001 Claim A (the repo's original core)** — parity in tournaments,
   Redei via odd-cycle collections — post-JC bonus: the S418 Redei-sign =
   discriminant-character reading gives OCF a Galois-theoretic face. Status: still
   the central named open problem of the repo.
2. **Flip-rank growth k(n)** (my k(7) = 12 correction; 4-agent convergence thread) —
   tournament flip-distance geometry; candidate for exact small-n closure.
3. **Paley heptagon extremality (HYP-3805)** — Paley tournament as LRC-extremal
   object; ties QR_p circulant structure to covering minima.
4. **Even-graph metagraph E_n (A002854)** — chi(E_n) growth, perfectness breaking at
   n=7 — a chromatic problem the repo owns outright.
5. **GLMY path homology of tournaments** — Betti stability conjectures + THM-125
   circulant speedup — engineering-mandate twin.
6. **Freiman/12-set uniqueness residual (HYP-7750)** — the 5% rational-time-evasive
   covering cores — additive combinatorics inside LRC.
7. **Erdos covering systems adjacency** — LRC covering sets vs covering congruences;
   the ghost-packing K(m+1)+1 duty is a covering-system statement in disguise.
8. **View-obstruction / billiards form of LRC (Cusick)** — the geometric carrier the
   lens map lists but no session has worked as primary.
9. **Sierpinski/Proth problem adjacency (S419/S420)** — which rows of the sheared
   Proth triangle contain primes; Proth-prime census 2,2,2,3,3,1,2,5,4,4,5,3,1 has
   no law — genuine number theory on our own object.
10. **The jet-illusion taxonomy (museum of impersonations)** — Moser circle,
    width-of-G_n, S419 diagonal law (MISTAKE-198) — formalizable: "which growth
    classes can share a k-jet?" A provable meta-theorem candidate.
11. **Unit distance / Hadwiger-Nelson adjacency** — Moser SPINDLE (not the circle) is
    a unit-distance graph; our PG(2,q) / Paley / spectral toolkit is the right shape;
    currently analogy-cards only. Fresh target.
12. **Mathieu subspaces of repo algebras** — beyond Zhao transport: which of OUR
    invariant subspaces (ghost-duty ledgers, weight-graded pieces of C[s,w]) are
    Mathieu subspaces? A definitional bridge nobody has attempted.
13. **q-deformed triangulars [n,2]_q and subspace towers (S420 census)** — the
    dyadic shadow of "edges"; connects to A006116 subspace counts and the tiling
    hypercube.
14. **OEIS-new pipeline** — kp's four + klein's T1532 batch + S420 candidates:
    curation debt with external-visibility payoff.

## Part IV — Priorities for the next sessions (assumption-challenging first moves)

**TARGET 1 (highest; untouched): explicit Zhao vanishing/image counterexample.**
¬JC collapses the universal forms of VC and IC, but NO explicit failing witness
exists anywhere yet. First move: run the de Bondt-van den Essen symmetric reduction
ON F (pin the exact dimension N and the homogeneous quartic P with F_sym = X + grad P);
then push through van den Essen-Wright-Zhao to the explicit (Delta, P, m) triple whose
vanishing pattern breaks, and the explicit Mathieu subspace that fails. Deliverable:
"the first explicit counterexample to the vanishing conjecture" — publishable alone.

**TARGET 2: DC1/DC2 via the eigen-graded bridge (THM-1345's opening).** Prove or
refute DC1 in the eigen-graded category; measure exactly which filtration step blocks
lifting to full DC1. The counterexample's 2-adic inverse ladder (THM-1300) suggests a
p-adic invariant separating A1 from A3 — chase it.

**TARGET 3: the classification theorem.** Upgrade S414b's local sporadicity to a
theorem: full trivial-group quotient + third-order obstructions; then run the
weight-type machine over the next torus weight systems (S414's template) — each hit
is a new sporadic, each miss extends the rigidity map. Merge with THM-1330's monoid
ideal into a structure theorem: "the counterexample ideal is generated over the tame
group by weight-type sporadics" (conjecture; test at low degree).

**TARGET 4: JC2 descent.** THM-1345 locates the hardness in {P_top,Q_top} = 0
weight-descent. First move: classify the descent obstructions at the first TWO
filtration steps (finite computation); any structure there is new after 85 years.

**TARGET 5: LRC(14) residuals.** The 5% covering-core gap (HYP-7750), F3 ceiling,
deficit-4 general proof, T(8,4) discrimination for the Rosetta correction law,
and the bilinear duty-sign law (which m*b+-1 denominators are achievable floors).

**TARGET 6: unit distance, for real.** First move: compute unit-distance embeddability
of our named graphs (Paley heptagon? E_n skeletons?) and test whether the PG(2,q)
count 2T+1 interacts with known H-N spectral bounds; even a clean negative
("our carriers are not unit-distance-realizable, here is the obstruction") upgrades
the analogy cards to mathematics.

**Standing rule for all targets:** state what assumption the move CHALLENGES
(the owner's brief). T1 challenges "equivalences preserve only truth values, not
witnesses"; T2 challenges "DC1 is safe because it is small"; T3 challenges
"counterexamples are a wilderness"; T4 challenges "JC2 hardness is monolithic";
T5 challenges "the floor is smooth above 1/14"; T6 challenges "our carriers are
metaphors".

---
Cross-links: THM-1300, THM-1305, THM-1310, THM-1315, THM-1320, THM-1325, THM-1330,
THM-1335, THM-1340, THM-1345, THM-1288/1289/1291/1292/1269, HYP-8070/8075/8155/8170/
8180, MISTAKE-198, engineering-synthesis-2026-03-10-S53 (the 12-domain mandate — the
Zhao/Mathieu transport and OEIS pipeline are its number-theory and visibility limbs).
