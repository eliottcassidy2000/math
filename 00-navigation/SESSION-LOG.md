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

