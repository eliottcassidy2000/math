## death-star-2026-07-21-S82 -- H≥disc (HYP-8636): disc is a MEAN OF PFAFFIAN SQUARES + the strong base's crux is the REGULAR tournaments (toward klein THM-1950's open base). HYP-8697.

**Owner directive:** work HYP-8636 (H≥disc) and related ideas.

- **STATE:** klein THM-1950 reduced H≥disc to the strong base H(C)≥max(1,s(C))·disc(C) (s=1ᵀ(I+K)⁻¹1), all SCC-composition machinery proved; the base is the open residual.
- **(1) disc = mean_{S even} Pf(K[S])²** (VERIFIED n≤5). The minor expansion det(I+K)=Σ det(K[S])=Σ Pf(K[S])² (skew ⇒ each principal minor is a Pfaffian-square) makes disc a normalized SUM OF SQUARES over the 2^{n-1} even subsets: det(I+K)=1+Σ_{|S|≥2}Pf²≥1, disc=1 ⟺ transitive. Recasts the base as 2^{n-1}·max(1,s)·H ≥ Σ_S Pf(K[S])² — a target for a COMBINATORIAL INJECTION (the disc-side combinatorics the eigenvalue route hides).
- **(2) THE CRUX = REGULAR TOURNAMENTS.** Base ratio H/(max(1,s)disc) tight (=1) only at C3; tightest non-C3 strong are the REGULAR ones. Paley-7: 189/(7·8)=3.375, BELOW n=6's 3.75 and klein's stated n=7 min 4.22 (sampling artifact; my n≤6 mins 1,1.67,1.875,3.75 match klein exactly). So the min ratio is NON-MONOTONE and the base reduces morally to H(regular)≥n·disc(regular), disc(doubly-regular)=(n+1)^{(n-1)/2}/2^{n-1} — the Paley-is-the-wall / big-stabilizer pattern (S75/S76).
- **(3) NOT a literal per≥det.** per(I+K)=−2<4=|det(I+K)| at C3, so H≥disc is not per(I+K)≥|det|; the per≥det (bosonic≥fermionic, THM-1810) lives at the GAUSSIAN-MOMENT level, and the Pfaffian-mean is the fermionic side's finite residue.
- Structural progress toward klein's base + the regular crux + a ratio-non-monotone correction. No new theorem. Credits klein THM-1950. reflection h-ge-disc-the-pfaffian-mean-...-S82; script h_ge_disc_pfaffian_mean_S82 (+out). GMC(2)/LRC(14) untouched.

## opus-2026-07-20-S445 - H at the FORMULA/#P edge: the harmonic boundary (THM-1970 + reflection)

Owner: H not poly-determined is an edge case; maybe a more refined invariant is the real answer;
tournaments sit at the edge between what a formula expresses and what provably cannot (harmonic series).

WORKED IT (THM-1970 + reflection H-at-the-formula-sharp-P-edge). Tournament invariants = a DEGREE-
GRADED poly tower (score deg1 -> c3 deg3 -> var=SC4 deg4 -> tr(S^{2j}) deg2j -> char_S all-moments),
each a poly degree-k census. H is captured by NONE: H|char_S splits at n=3 (THM-1935, full spectrum
misses H); the H-defect within a degree-k census GROWS with n at fixed k (k3: 4->14, k4: 2->12 for
n=5,6), vanishing only at the deck k=n-1 (reconstruction). So H needs FULL-SUPPORT = the PERMANENT
(#P) to char_S's DETERMINANT (poly); the gap = the permanent/determinant boundary. Complete-signed-
relation => spectrum poly, path-count #P -- WHY tournaments sit on the edge.

REFINED OBJECT: scalar H is not even COMPOSITIONAL (H(C3[S1,S2,S3]) != f(H(Si)); block-H (1,1,1) ->
composites {3,...,2721}); the compositional refinement = the path-SYSTEM (linear-forest) polynomial
(categorifies H), functorial NOT a poly formula (none exists unless P=#P). The 'more refined answer'
buys composition, not complexity.

HARMONIC ANALOGY made exact: moment tower = partial sums; char_S = zeta(s>1) (formula); H = the pole
at s=1 (the edge); the char_S->H defect = the tournament gamma (the anomaly after all poly data),
physically = THM-805 resistance=harmonic-number, CLAUDE.md's gamma. OPEN: relative defect (small-n
0.53->0.62 => edge is REAL not measure-zero); the path-system transfer (resolves THM-1960 cyclic-H).

Honest correction mid-session: size-controlled test (n=5 blocks, same H, diff PH) gave EQUAL symmetric
composites, so the PH-non-composition is subtler than first claimed; the robust results are the
harmonic-edge defect table + H|char_S + H(C3[.]) not scalar-H-determined across sizes.

Files: THM-1970; HYP-8715; reflection; refined_H_and_harmonic_edge + PH_composes _opus_S445.py (+out).
Namespace clean (1970/8715). Builds on THM-1935/1940/1945/1960/1930, THM-805.

## opus-2026-07-20-S444 - Tournaments compose from REGULAR SEEDS: the spectral substitution law + octonion object C3[C3] (THM-1960)

Owner: consider the three recursion modes (A+B+C-D-E-F+G / A+B-C+D-E-F+G / A+B-C) + tournaments as
recursively composed of subtournament seeds; what iso-class seeds correspond to larger tournaments?
FRAME = modular/substitution (Gallai) decomposition; seeds = PRIME (no nontrivial module) tournaments.

THM-1960: (1) SEED CENSUS (modular-prime) = 1,1,1,0,3,15 for n=1..6; C3 = smallest nontrivial seed;
transitive = fully-linear; n=4 has ZERO seeds. CRUCIAL DISTINCTION (from the prior-work map): the
repo's order-join/SCC 'atoms' are STRONGLY-CONNECTED (THM-1862, 1,1,6,35,353) but modular-primes are
STRICTLY stricter -- the SC 4-tournament is an order-join atom yet has a module (0 modular seeds at
n=4). Substitution carves INSIDE strong components. (2) SPECTRAL SUBSTITUTION LAW: for T[S^m], skew =
S_T (x) J_m + I_k (x) S_S; when the SEED is REGULAR (row-sums 0 => block-mean in ker), nz-spec(T[S^m])
= [nz-spec(S) x k] U [m*nz-spec(T)]; verified regular C3,C5, FAILS irregular T3. => all even moments of
regular-seed substitution objects are closed forms in seed moments. (3) OCTONION OBJECT C3[C3] (n=9):
regular, char_S=x(x^2+3)^3(x^2+27), var=104, SC4=81=3^4, H=3159=3^5*13 -> exact degree-8 test object
for the octonion wall (THM-1935/1940). (4) the 3 modes = Mobius/Legendre/Eisenstein CHARACTERS
(+++---+ / ++-+--+ / ++-); ++- = the C3 base seed.

CONCURRENT/PRIOR-WORK: coordinates with boxeph stub THM-1955 ('which iso-classes come from smaller',
CLAIMED, circulant-character reduction DAG) -- I fill the modular-prime + spectral-substitution axis,
complementary not overriding. char_A=prod char_A(SCC) (THM-1830/1925) = the linear-quotient special
case; my skew regular-seed law is the CYCLIC-quotient generalization. H=prod_modules H (S531 apex-
recursion) = transitive-quotient case; my cyclic H(C3[C3])=3159 is new. Fixed a stale THM-1855 alias
-> THM-1862.

OPEN: H under cyclic substitution (the 13 in 3159); tr(S^8) octonion-wall test on C3[C_{2j+1}]; seed
census to n=7.

Files: THM-1960; HYP-8705; seed_tournaments_and_substitution + spectral_substitution_law _opus_S444.py
(+out). Namespace clean (1960/8705). Builds on THM-1920/1935/1940, cites THM-1862/1830/1955/442/830.

## death-star-2026-07-21-S81 -- CHASED THE PREDICTIONS: Pell supersymmetry is a real GMC(2) moment identity (E[sym²]−E[alt²]=E[(ZW)^n]=n!); α(E_n)=2^{n−4} confirmed to n=7 (α(E_7)=8); recursion modes resolve to the order-join. HYP-8696.

**Owner directive:** chase the predictions; think the A+B+C−D−E−F=G / A+B−C recursion modes as literally subtournaments — which iso classes come from smaller subtournament classes.

- **PELL SUPERSYMMETRY — CONFIRMED (the S80 prediction).** The a/b Pell identity E_n²−O_n²=(x²−1)^n (kps THM-1880) is the polynomial lift of an exact GMC(2) MOMENT identity: with sym_n=(Z^n+W^n)/2 (bosonic), alt_n=(Z^n−W^n)/2 (fermionic), exact Wick gives **E[sym_n²]−E[alt_n²] = E[(ZW)^n] = n!** (n=1..7), and generally E[sym(P)²−alt(P)²]=E[P·P̃] (charge-conjugate norm). Bosonic²−fermionic²=radial-norm, localizing on n! for P=Z^n — a positivity-flavored (radial mass n!>0), provable toral/parity companion to GMC(2)'s open radial gap (not the gap itself).
- **α(E_n)=2^{n−4} — CONFIRMED to n=7.** Built E_7 (V=54=A002854 exact, hash-canonical), α(E_7)=8=2³. So α(E_n)=1,2,4,8 at n=4..7; predicts α(E_8)=16. The G_n-dual α(G_n)=2,5,18 has no closed form — a metagraph/even-dual asymmetry.
- **RECURSION MODES — resolved to the order-join.** Signed-additive deck reading (±patterns on vertex-deleted subtournament invariants) is INCONCLUSIVE (recorded dead end). The clean "smaller classes → larger class" mechanism is the MULTIPLICATIVE order-join (condensation into strong atoms): H,|R| multiply (mac-mini THM-1936), disc super-multiplies (klein THM-1950 SL2 velocity-addition s(C1⇒C2)=(s1+s2)/(1+s1s2)); reducible classes are built from strong atoms; the ± are the R=Σsgn(π) permutation signs. So "which iso classes come from smaller subtournament classes" = the REDUCIBLE ones; strong tournaments are the atoms.
- **FLEET (credited):** klein THM-1950 reduced my HYP-8636 (H≥disc) to the strong base via this exact composition algebra — one strong-base inequality from proof. mac-mini THM-1936 (signed R multiplicative), kps THM-1885/1880 (a/b=BS(1,2)), boxeph THM-1926 (zeta). reflection chasing-the-predictions-...-S81; script pell_supersymmetry_and_deck_recursions_S81 (+out). GMC(2)/LRC(14) open; no LRC re-audit.
## kind-pasteur-2026-07-21-S128c142 - The tournament-invariant lattice DEFINITIVELY mapped; 3 conjectures resolved (THM-1965)

Owner: work reframings & conjectures, investigate more of them, chase definitive results. Took the
[C]/[R] statements from my THM-1945 dictionary and drove the tractable ones to definitive verdicts.

**METHOD:** reframe 9 tournament invariants as Sn-orbit functions ordered by "f refines g iff same-f
=> same-g"; compute the FULL refinement Hasse diagram + first-separation n, exhaustively over EVERY
iso class n=3..6 (bit-packed canon; script invariant_lattice_definitive_kps_S128c142.py +out).
Invariants: score, specA, specS, cyc=(c3..cn), H=#Ham-paths, R=signed-Redei, disc=|det(I+K)|/2^{n-1},
arb=arborescences, aut=|Aut|. VALIDATED: H odd (Redei); |R| spectrum n=5={1,3,5,7,11,15}, max 3,3,15,15
(= mac-mini THM-1936); disc = klein THM-1950.

**SIX DEFINITIVE FINDINGS (THM-1965):**
1. **CUT/CYCLE INCOMPARABILITY (headline).** score (cut-space/hierarchy) and cyc (cycle-space) are
   INCOMPARABLE from n=5 -- the lattice shadow of the GF(2) cut+cycle direct sum. arb (Laplacian/cut
   side) also incomparable to cyc. The project's core duality as an order-theoretic fact.
2. **cyc = the cycle-side MASTER invariant** (refines specA,specS,H,R,disc,aut; not score/arb). specA
   determines specS but not conversely (n<=6).
3. **No poly-time invariant refines H from n=5** (disc/arb/specA/specS/score all fail) = lattice
   restatement of THM-1780/1865 (H #P leaves the ladder). disc INCOMPARABLE to H (klein's H>=disc is a
   numeric bound, NOT a refinement) = the tournament permanent-vs-determinant pair.
4. **mac-mini's OPEN QUESTION ANSWERED:** R is INCOMPARABLE to both H and specA (R separates same-H
   from n=5, cospectral from n=6; H/specA separate same-R from n=4) => R a genuinely new coordinate.
5. **ISO RECONSTRUCTION:** {score,specA} pins iso n<=5, misses 10 at n=6; no proper subset of the 9
   determines iso from n=6; full 9-tuple does (n<=6).
6. **[C12] REFUTED:** metagraph G_n NOT regular (deg 3..14 n=6), NOT vertex-transitive, NOT Cayley --
   the transitive corner is a distinguished vertex, so the spine/ribs/sea + principal-line geometry IS
   the absence of symmetry. 1-WL canonical 34-class coloring at n=6.

**FLAGGED n<=6-only (conjecture n>=7):** cyc->H (cycle vector determines the Redei count -- SHARPEST
open lead, test n=7) and specA->specS. Both possibly fineness artifacts (cyc=32 of 56 values n=6).

**INTEGRATION:** answers mac-mini-S160's closing R-question; places klein-S400's H>=disc in the
lattice; reinforces THM-1780/1865; maps THM-1945's lattice entry (4) in full; refutes my own [C12].
Reflection creative-statements-...-S128c141 updated with resolution verdicts. Reframing + definitive
resolutions, all VERIFIED-EXACT n<=6, no new open problem opened beyond the flagged cyc->H lead.

**Namespace:** THM-1965 + HYP-8730 (hot: klein 1950, opus 1940/1955-region; took margin). No collision,
no canon overridden. **Files:** THM-1965; reflection update; script + .out; HYP-8730.

## opus-2026-07-20-S443 - Concrete next steps: var(lambda^2) is a 4-subtournament-census invariant (THM-1940) -- resolves THM-1930, pins the quaternion-wall mechanism
## boxeph-2026-07-21-S196 -- THM-1955 WHICH ISO CLASSES COME FROM SMALLER (reduction DAG); the owner's recursion modes ARE circulant characters

**Owner:** work handoffs; think about A+B+C−D−E−F=G, A+B−C+D−E−F+G, A+B−C as literal subtournaments; which iso classes come from smaller subtournament classes.

**The modes are circulant CHARACTER sums (verified).** Σ_kε_kω^{jk}=2λ_j+1 (j≠0) = 'A+B+C−D−E−F=G':
- +++−−− = sign(sin2πk/7) = interval {1,2,3} (Dirichlet character)
- ++−+−− = Legendre (k/7) = QR {1,2,4} = Paley-7 (Gauss character)
- ++− = ℤ/3 base = the 3-cycle (=Paley-3)
FINDING (REFINES my THM-1925): Re(λ_j)=−1/2 for ALL ℤ/p circulants (the pair {k,p−k} shares a cosine). The modes differ ONLY in the imaginary/sine part: Paley FLAT |2λ+1|=√7 (Gauss sum), interval SPREAD (Dirichlet, which vanishes at the roots of unity so its Re is −1/2 too). = kps's char_S spread-duality (THM-1880) on the adjacency eigenvalues; the −1/2 line is kps's b, the Gauss concentration is the GIT-polystable pole.

**Reduction DAG census (answers 'which iso classes come from smaller'), n=3..7:**
- REDUCIBLE (order-join composites; char_A/H/signed-R/ζ factor over strong components): 1,3,6,21,103
- STRONG-PRIME: 1,1,6,35,353;  of which CIRCULANT (mode/character-generated): 1,0,1,0,2
So the owner's recursions cover the reducible bulk (order-join = literal subtournaments) + a THIN circulant thread (character modes); 351/456 at n=7 are PRIME atoms from nothing smaller — which is why the general theory is hard.

**Deck recursion (literal subtournaments A..G = vertex-deleted):** char_A'(x)=Σ_i char_{T−v_i}(x) and Σ_i c3(T−v_i)=(n−3)c3(T), verified n≤6. The unsigned deck sum = the char-poly derivative; signed deck sums are twisted functionals (degenerate on vertex-transitive). Character modes live on ℤ/n (build the circulant atoms); the deck recursion lives on the vertices (relate a class to its subtournament deck) — the two faces of the reduction DAG.

**Integrated:** mac-mini THM-1936 (signed Rédei R join-multiplicative — reducible R = product over strong comps, the signed companion of my ζ/THM-1862); kps THM-1880 char_S spread-duality (the skew mirror of my adjacency flat/spread). Worked the reduction handoffs (char_A/H/R/ζ all factor; the DAG unifies them).

**Housekeeping:** THM numbers were flying — claimed THM-1955 (above the 1950 max), HYP-8720/8721 CONFIRMED. Added a refinement note to my THM-1925 leg (c). No collisions this session.

**Honest scope:** census + deck recursion + character-sum identity are verified exact n≤7 (n≤6 deck). char_A=∏char(SCC) is THM-1926; the NEW content is the reducible/circulant/prime CENSUS answering 'which come from smaller', the modes=characters identification, and the Re=−1/2-universal refinement of THM-1925.

**Next:** (1) the prime-atom interior (351/456 at n=7) — the genuinely-hard part no recursion reaches; can a finer (modular-decomposition / lexicographic) recursion carve it? (2) Lean the block-triangular char_A=∏char(SCC). (3) the two-variable Ihara zeta (THM-1926 handoff). Artifacts: THM-1955, HYP-8720/8721, reflection the-recursion-modes-are-characters-and-the-reduction-dag-boxeph-S196.md, script reduction_dag_recursion_modes_boxeph_S196.py (+.out).

