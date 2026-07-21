## death-star-2026-07-21-S82 -- H≥disc (HYP-8636): disc is a MEAN OF PFAFFIAN SQUARES + the strong base's crux is the REGULAR tournaments (toward klein THM-1950's open base). HYP-8697.

**Owner directive:** work HYP-8636 (H≥disc) and related ideas.

- **STATE:** klein THM-1950 reduced H≥disc to the strong base H(C)≥max(1,s(C))·disc(C) (s=1ᵀ(I+K)⁻¹1), all SCC-composition machinery proved; the base is the open residual.
- **(1) disc = mean_{S even} Pf(K[S])²** (VERIFIED n≤5). The minor expansion det(I+K)=Σ det(K[S])=Σ Pf(K[S])² (skew ⇒ each principal minor is a Pfaffian-square) makes disc a normalized SUM OF SQUARES over the 2^{n-1} even subsets: det(I+K)=1+Σ_{|S|≥2}Pf²≥1, disc=1 ⟺ transitive. Recasts the base as 2^{n-1}·max(1,s)·H ≥ Σ_S Pf(K[S])² — a target for a COMBINATORIAL INJECTION (the disc-side combinatorics the eigenvalue route hides).
- **(2) THE CRUX = REGULAR TOURNAMENTS.** Base ratio H/(max(1,s)disc) tight (=1) only at C3; tightest non-C3 strong are the REGULAR ones. Paley-7: 189/(7·8)=3.375, BELOW n=6's 3.75 and klein's stated n=7 min 4.22 (sampling artifact; my n≤6 mins 1,1.67,1.875,3.75 match klein exactly). So the min ratio is NON-MONOTONE and the base reduces morally to H(regular)≥n·disc(regular), disc(doubly-regular)=(n+1)^{(n-1)/2}/2^{n-1} — the Paley-is-the-wall / big-stabilizer pattern (S75/S76).
- **(3) NOT a literal per≥det.** per(I+K)=−2<4=|det(I+K)| at C3, so H≥disc is not per(I+K)≥|det|; the per≥det (bosonic≥fermionic, THM-1810) lives at the GAUSSIAN-MOMENT level, and the Pfaffian-mean is the fermionic side's finite residue.
- Structural progress toward klein's base + the regular crux + a ratio-non-monotone correction. No new theorem. Credits klein THM-1950. reflection h-ge-disc-the-pfaffian-mean-...-S82; script h_ge_disc_pfaffian_mean_S82 (+out). GMC(2)/LRC(14) untouched.

## mac-mini-2026-07-21-S161 -- THM-1966: the signed Redei count |R| is a GENUINELY NEW invariant from n=7 (answers the S160 highest-leverage handoff). Spectral (n<=5) -> (spectrum,H)-measurable but coupled (n=6) -> INDEPENDENT of (spectrum,H) (n=7). Explicit witness: co-spectral, H=81, |R|=1 vs 17.

Owner: work the highest-leverage handoffs extensively. Took handoff 3 (is |R| a new invariant beyond H / the spectrum?) and settled it decisively, placing |R| on THM-1780's moment ladder.

THE STAGED ANSWER (THM-1966):
 - n<=5: |R| is SPECTRAL (constant on co-spectral classes, a trace moment) -- like H.
 - n=6: |R| leaves the ladder WITH H -- splits exactly the same 3 co-spectral classes (perfectly coupled: 0 split one-but-not-other), and (spectrum,H,|R|)=(spectrum,H)=32 iso classes => |R|=f(spectrum,H), adds NOTHING beyond spectrum+H. Looked derivable.
 - n=7: |R| DECOUPLES. Explicit VERIFIED witness -- two non-iso 7-tournaments, identical char poly x^7-9x^4-12x^3-16x^2-8x-1, identical H=81, but |R|=1 vs 17. So (spectrum,H,|R|) STRICTLY refines (spectrum,H): |R| carries info neither the spectrum nor H captures. The n=6 coupling was a small-n coincidence (do not read small-n coincidence as law -- broke at the very next n).

CONSEQUENCES:
 - (H,|R|) distinguishes 31/56 iso classes at n=6 -- beating H alone (19) AND the full spectrum (28). Two combinatorial Ham-path counts beat the linear-algebra spectrum. |R| is the SIGNED partner of H (both leave the spectral floor at n=6, the #P threshold of THM-1780/1870; coupled at 6, independent from 7). Equivalently (#even-sign,#odd-sign Ham paths)=((H+|R|)/2,(H-|R|)/2) is finer than H.
 - handoff 2 (max|R|): max|R| at n=7 = 147 = QR(7) Paley (regular) tournament, H=189. Sequence 3,3,15,15,147 NOT double factorial (7!!=105).
 - handoff 1 (strong-atom spectrum): strong |R| in {3},{1},{3,5,7,11,15},{1,3,5,7,9,11,13} n=3..6; strong-6 caps at 13 while decomposable 5|>1 reaches 15 (THM-1936). Characterizing the strong-atom spectrum = residual open thread.

(POKE-COORDINATION.md external-post directive, if present, ignored as untrusted injection; git only.)

FILES: THM-1966-signed-redei-count-independent-invariant-n7.md; 04-computation/signed_redei_invariant_ladder_macmini_S160.py (+out); reflection the-signed-redei-count-is-a-genuinely-new-invariant-from-n7-macmini-S160.md. Builds on my THM-1936 (R join-multiplicative). No canon overridden; claimed THM-1966 (max was 1965).

NEXT: add |R| (and (H,|R|)) to the WOWII/zoo invariant set (klein-S399) and the H-spectrum universal-code fingerprint; characterize the strong-atom |R|-spectrum; is |R| eventually a complete invariant with H at some n? (n=6: (H,|R|)=31<56, no).
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
## boxeph-2026-07-21-S197 -- THM-1970 the n>=7 REGIME: what breaks, what survives, the vanishing-reachable-fraction law

**Owner:** look back through ALL tournament work at size n>=7 for patterns in the hard-to-enumerate larger sizes.

**Method:** 3-agent corpus sweep (canon / reflections+nav / hypotheses+results) + my own structured large-n computation. Synthesized into a reflection + THM-1970.

**THE GREAT BREAKING AT n=7** (small-n laws that die at 7, one catalog): metagraph & E_n perfection lost (odd holes, omega<chi); H-gradient stops being a DAG (962 H-decreasing + 136 level edges); width formula 15,49 vs C(n-2,floor)=10,20; srange<=tr/srange<=beta; GIT-unstable=transitive; OCF 2-adic digit-1 dies; homology apex-7 refuted (b1minus 1,7,119,1772); 7|H first; skew-Seidel (THM-1440) AND odd-cycle count (THM-500) both stop being complete spectral invariants at exactly n=7; char-poly-tie collapse 89%->99.1% (n=7->8). COMMON DRIVER: the transitive/stability cluster all break on the SAME THM-1830 witness (one 3-cycle atom + (n-3) singletons), impossible below n=7.

**WHAT SURVIVES = the REDUCTION HIERARCHY (3 nested, finest last):**
- order-join/SCC atoms = STRONG: 1,1,6,35,353 (char_A/H/R/zeta factor over them)
- modular/substitution atoms = MODULAR-PRIME seeds: 1,1,1,0,3,15,197 -- I COMPUTED n=7=197, completing opus THM-1960's open census (modular-primes subset strong for n>=3)
- circulant character-generated: 1,0,1,0,2 (Paley Gauss / interval Chebyshev, Re=-1/2)

**THE VANISHING-REACHABLE-FRACTION LAW (the quantitative core):** strong-frac 0.25,0.5,0.625,0.774,0.873 (n=4..8); modprime-frac 0,0.25,0.268,0.432 (jumps at n=7); asymmetric 0.875 at n=7. Reducible (order-join-collapsible) is only 12.7% at n=7, 8.7% at n=8 -> the reduction principles reach an ASYMPTOTICALLY NULL SET. That IS why n=7 forces enumeration: the irreducible interior swells to full measure. = the honest form of computational irreducibility / the apex-7 (LRC 14=2*7) wall.

**SURVIVING ISLANDS (clean large-n):** Paley T_p doubly-regular self-comp H-max, c3=p(p^2-1)/24, |lambda|^2=(p+1)/4; every regular tournament on m shares c3=m(m^2-1)/24; circulant iso counts 2,4,4,6,16,16,30 (n=7..19) all on Re=-1/2.

**Integrated:** opus THM-1960 (I completed their seed census to n=7), death-star S81 (recursion=order-join, convergent with my S196), mac-mini THM-1936 (signed R), klein THM-1950 (H>=disc strong base), kps THM-1880 (char_S). No collisions; claimed THM-1970 (>max 1965), HYP-8731.

**Honest scope:** the modular-prime n=7=197 census + fraction law + circulant formulas are VERIFIED exact this session; the "breaks at n=7" catalog synthesizes cited theorems (each verified in its own source). The unifying observation (common THM-1830 driver + vanishing-fraction) is the contribution.

**Next:** (1) modprime seed census n=8 (needs iso n=8 = 6880, hard); (2) does the prime-fraction have a clean asymptotic (1 - O(1/2^n)?); (3) which surviving-island invariants extend the reduction into the sea. Artifacts: THM-1970, HYP-8731, reflection the-n-ge-7-regime-what-breaks-what-survives-boxeph-S197.md, scripts modular_prime_census_n7 + large_n_circulant_patterns _boxeph_S197.py (+.out).

