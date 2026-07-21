## kind-pasteur-2026-07-21-S128c144 - The figurate reciprocal ladder and the harmonic edge: our sequences as sub-series of the harmonic series (THM-1990)

Owner: reciprocal of an integer sequence = a subset of the harmonic numbers; 1+1/2+..+1/5 > 2 while
sum 1/T_n = 2; study our sequences' reciprocal sums extensively and extend. Looked back at every
integer sequence in the corpus through the reciprocal-sum lens.

**PROVEN CORE (THM-1990) -- the figurate ladder.** Telescoping (exact rational, verified):
sum_{n=k}^N 1/C(n,k) = k/(k-1)*(1 - 1/C(N,k-1)) => sum_{n>=k} 1/C(n,k) = k/(k-1) for k>=2, from
1/C(n,k) = k/(k-1)[1/C(n-1,k-1) - 1/C(n,k-1)]. THE LADDER by dimension:
  k=1 vertices n = HARMONIC = DIVERGES (the edge); k=2 arcs/triangular = 2 (the owner's value);
  k=3 tetrahedral = 3/2; k=4 = 4/3; ... -> 1.
Arc count of K_n = C(n,2), so the corpus' central staircase sequence has reciprocal signature EXACTLY
2, and its telescoping 1/C(n,2)=2(1/(n-1)-1/n) IS the arc->vertex (dim2->dim1) reduction.

**THE HARMONIC EDGE = p=1 = dimension 1, UNIFYING THREE LENSES.** The ladder diverges only at k=1
(the ground set). This is the SAME p=1 boundary as THM-1980 (H's formula content collapses to one
bit = 'H at p=1') and THM-1870 (cycle counts turn #P at the Hamiltonian length). Three independent
lenses -- reciprocal convergence, 2-adic depth, cycle length -- place the marginal case at the same
p=1 corner. The project's dimensional ladder n->C(n,2)->C(n,3) IS the figurate reciprocal ladder
crossing that edge.

**RECIPROCAL SIGNATURE table (VERIFIED to named constants):** arcs 2, tetrahedral 3/2,
var(lambda^2)=2C(n,3) -> 3/4, squares pi^2/6, factorial e, central binomial 4/3+2pi sqrt3/27 (EXACT
18 digits), Catalan 2+4pi sqrt3/27, Fibonacci=reciprocal-Fibonacci const, Mersenne 2^n-1=Erdos-Borwein,
2^n=1. 2-ADIC/THETA: labeled tournaments sum 1/2^{C(n,2)}=1.6416325607... (partial theta q=1/2);
switching classes 2^{C(n-1,2)} = 1 + that (the +1 = extra n=1 term, PROVEN exact); SIGNED pentagonal
= Euler prod(1-2^-n)=(1/2;1/2)_inf (pentagonal number thm, THM-488 hub). Census fingerprints:
A000568=3.8535, A002854=3.0618, score=3.9325, tangent=1.5663, secant=2.2171.

**CONVERGENCE DICHOTOMY:** sum 1/a diverges iff a grows at most linearly; EVERY combinatorial repo
sequence (degree>=2) converges, only the linear ground-set ones (vertices, odds, H-spectrum) sit on
the harmonic edge -- a clean 'ground set vs structure' separation.

**EXTEND:** sigma(a)=sum 1/a_n as a sequence fingerprint invariant; telescoping = the a-monoid
(THM-1885, a:n->n+1) Mode-A action on figurate dimensions; a transcendence gradient (poly->rational,
exp->irrational analytic, lacunary->theta). NEXT: identify the theta constant 1.6416...; signed
reciprocal sums Sum(-1)^n/a_n for tournaments/even-graphs; Sum 1/H(T) as bridge to THM-1980.

Reframing + a proven ladder + a verified signature table; classical constants (figurate, Basel,
Erdos-Borwein, pentagonal thm) unified as one harmonic-subset classification of the corpus.
**Files:** THM-1990; reflection the-harmonic-boundary-p-equals-1-recurs-kps-S128c144; HYP-8750;
script reciprocal_sums_of_our_sequences_kps_S128c144 (+out). Namespace: THM-1990/HYP-8750 (clean).
## death-star-2026-07-21-S83 -- H≥disc: the REGULAR SUB-BASE reduced to ONE average (H(reg)≥Szele avg), proved n≥7 modulo that crux; the Pfaffian injection IS the even/odd (Ω/E_n) duality. HYP-8698.

**Owner directive:** work the Pfaffian injection and the regular sub-base (S82 handles toward klein THM-1950's open base H(C)≥max(1,s)disc(C) for HYP-8636).

- **REGULAR SUB-BASE H(reg)≥n·disc(reg): PROVED for n≥7 modulo one crux.** Chain H(reg) ≥(i) n!/2^{n-1} ≥(ii) n(n+1)^{(n-1)/2}/2^{n-1} ≥(iii) n·disc(reg).
  - **(iii) PROVED (AM-GM):** disc(reg)=∏(1+μ_j²)/2^{n-1}, Σμ_j²=C(n,2) fixed ⇒ ∏ maximized at equal μ_j²=n (doubly-regular=Paley) ⇒ disc(reg)≤(n+1)^{(n-1)/2}/2^{n-1}. Tight Paley-3,7,11.
  - **(ii) PROVED (elementary):** (n-1)!≥(n+1)^{(n-1)/2} — fails n=3,5, holds n=7 (720≥512), ratio increasing ⇒ all n≥7.
  - **(i) THE CRUX (conjecture, strongly evidenced):** every regular tournament has ≥ the Szele average n!/2^{n-1} Hamiltonian paths. Exhaustive n=3 (H=3≥1.5), n=5 (all 24 regular H=15≥7.5); samples n=7 (min 171≥79), n=9 (min 3243≥1418), huge margins. n=3,5 direct.
  - So the "regular is the wall" crux (S75/S76) is now a SINGLE tractable Ham-path statement (plausibly Moon/Alon/Busch), far easier than the eigenvalue-product original. Doubly-regular (Paley) tightest.
- **PFAFFIAN INJECTION:** aggregate 2^{n-1}H≥Σ_{S even}Pf(K[S])²=det(I+K) confirmed with room (slack 112,416 at n=5,6). STRUCTURAL READING: disc=Σ Pf² counts EVEN cycle-covers (cycle-space); H=I(Ω,2) counts via ODD cycles (OCF); so H≥disc = "the ODD (OCF) count dominates the EVEN (Pfaffian) count" — the even/odd, cut/cycle, E_n/Ω duality (S80). Per-subset injection open (Pf(K[S])²≤H(T[S])H(T\S) not clean; the right compatibility is subtler).
- Reduces HYP-8636's open crux to (i); no new theorem. Credits klein THM-1950. reflection the-regular-sub-base-...-S83; script h_ge_disc_regular_subbase_S83 (+out). GMC(2)/LRC(14) untouched.

## opus-2026-07-20-S446 - The path-cover polynomial is the refined compositional invariant; the formula/#P edge is REAL (THM-1975)

Owner: work the cleanest next computations (the refined-object + relative defect from THM-1970/1960).
BOTH resolved through one object -- the path-cover polynomial pc.

(A) pc(S,c) = # partitions of V(S) into c vertex-disjoint directed paths; pc(S,1)=H(S) is its TOP
coefficient. H(C3[S1,S2,S3]) is a FUNCTION of (pc(S1),pc(S2),pc(S3)) -- 0/20 pc-triples ambiguous over
all size-<=3 blocks -- but NOT of scalar H (3/4 ambiguous). So pc COMPOSES where scalar H does not:
H(C3[.]) = sum K(c1,c2,c3) prod pc(Si,ci), block-independent kernel K. RESOLVES THM-1970's 'more
refined than H = the real answer' (= pc) and THM-1960's cyclic-H (the 13 in H(C3[C3])=3159 = cyclic
interleaving of path-systems). pc is still #P (top coeff = H) -> the refinement is FUNCTORIAL, not a
complexity reduction, exactly as THM-1970 predicted.

(B) THE EDGE IS REAL: the relative H-defect defect_3(n)/Hbar(n) GROWS: 0.533 (n=5 exact), 0.622 (n=6
exact), >=0.915 (n=7, 400k-sample). H lives on the FAR side of the formula/#P boundary -- the part
invisible to a bounded poly census has positive+growing weight, NOT a measure-zero (gamma-like)
correction. Tournaments occupy the harmonic edge and CROSS it; the poly tower is the largest formula-
expressible shadow of an object past expression.

THM-1960->1970->1975 assemble one statement: tournaments built by substitution over seeds; their
formula-expressible invariants are a degree-graded poly ladder (score/c3/var/char_S); H is the first
invariant past the top, refined not by a formula but by the functorial path-cover polynomial that
composes over the seeds.

OPEN: exact kernel K; general prime-quotient pc-law; does defect_3/Hbar -> 1?

Files: THM-1975; HYP-8725; resolution banners on THM-1970/1960; path_cover_transfer + relative_defect_n7
_opus_S446.py (+out). Namespace clean (1975/8725).

## kind-pasteur-2026-07-21-S128c143 - THE 2-ADIC EDGE of H: the spectrum resolves H to a depth that DECAYS to one bit; Redei's parity is the LAST formula (THM-1980)

Owner: H being undetermined by poly-time invariants is an EDGE case; tournaments sit at the boundary
between formula-expressible and provably-not, like the harmonic series at p=1; maybe a more refined
invariant is the real answer. CHASED to a definitive 2-adic statement.

**THE INSIGHT MADE PRECISE.** Define d(n) = largest k with "H mod 2^k" constant on every cospectral
class (= a function of char_A; the spectrum is the poly-time ladder THM-1775/1780). Redei (H odd)
=> d(n) >= 1. **THM-1980, DEFINITIVE (exact integer arithmetic, no floats):**
- d(4)=d(5)= infinity  (the spectrum DETERMINES H completely at n<=5)
- d(6)=2  (spectrum pins H mod 4 but NOT H mod 8; the mod-8 bit turns #P at EXACTLY the cospectral
  class {13,17} where H leaves the ladder, THM-1780 -- both ==1 mod 4, differ by 4)
- d(7)=1  (26 sampled cospectral 7-buckets carry BOTH residues of H mod 4 => even the mod-4 bit is
  #P; d(7)=1 is DEFINITIVE by witness existence, immune to sampling incompleteness)
=> depth = inf, inf, 2, 1 for n=4..7, hitting the Redei floor at n=7: **ASYMPTOTICALLY THE SPECTRUM
PINS EXACTLY ONE BIT OF H (its parity); every higher bit is #P. Redei's theorem is the LAST formula.**

**H mod 4 is a real but FLEETING spectral invariant:** both residues occur ~equally (n=6: 32:24),
constant on cospectral classes for n<=6, NOT score/c3-determined from n=5 ((H-1)/2==c3 mod2 only at
n=4), then STOPS being spectral at n=7.

**TWO ORTHOGONAL EDGES MEET AT THE HAMILTONIAN OBJECT:** the LENGTH edge (THM-1870: cycle counts c_k
poly for k<=n-1, #P at the Hamiltonian length k=n) + this 2-ADIC edge (H's bits poly up to parity,
#P above). H is one length past the spectral cycle counts and one bit past a spectral formula -- the
marginal object on both axes (the harmonic-series-at-p=1 analogue). Mechanism: per=det mod 2 is
exactly one bit deep.

**OPEN (the owner's own next question):** does ANY poly-time invariant (not just the spectrum) beat
the parity bit asymptotically? The full poly-tuple (score,specA,specS,disc,arb_inv) determines H for
n<=6; its n>=7 depth is untested. A NO would PROVE Redei's bit is H's entire formula-expressible
content -- the exact statement of "tournaments at the edge."

**CORRECTION (MISTAKE-208):** the `arb` in the S128c142 lattice was arborescences ROOTED AT VERTEX 0
-- NOT iso-invariant (root depends on labeling; a random-label sampler exposed it via spurious
collisions). Proper arb_inv = sorted per-root tuple: |arb_inv|=55 at n=6 (nearly complete), refines
score, incomparable to specA/cyc/H -- STRENGTHENS THM-1965's cut/cycle story. Headline (score ⟂ cyc,
THM-1980) uses only exact invariants, unaffected. Banner added to THM-1965.

**Files:** THM-1980; MISTAKE-208; reflection update; HYP-8740; scripts H_two_adic_edge_v2 /
H_mod4_formula_and_n7 / exact_lattice_and_edge _kps_S128c143 (+outs). Cites THM-1780/1870/1965/1945/1775.
Namespace: THM-1980, HYP-8740 (hot: death-star 8697, THM-1970-region; took margin).

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
## boxeph-2026-07-21-S197 -- THM-1978 the n>=7 REGIME: what breaks, what survives, the vanishing-reachable-fraction law

**Owner:** look back through ALL tournament work at size n>=7 for patterns in the hard-to-enumerate larger sizes.

**Method:** 3-agent corpus sweep (canon / reflections+nav / hypotheses+results) + my own structured large-n computation. Synthesized into a reflection + THM-1978.

**THE GREAT BREAKING AT n=7** (small-n laws that die at 7, one catalog): metagraph & E_n perfection lost (odd holes, omega<chi); H-gradient stops being a DAG (962 H-decreasing + 136 level edges); width formula 15,49 vs C(n-2,floor)=10,20; srange<=tr/srange<=beta; GIT-unstable=transitive; OCF 2-adic digit-1 dies; homology apex-7 refuted (b1minus 1,7,119,1772); 7|H first; skew-Seidel (THM-1440) AND odd-cycle count (THM-500) both stop being complete spectral invariants at exactly n=7; char-poly-tie collapse 89%->99.1% (n=7->8). COMMON DRIVER: the transitive/stability cluster all break on the SAME THM-1830 witness (one 3-cycle atom + (n-3) singletons), impossible below n=7.

**WHAT SURVIVES = the REDUCTION HIERARCHY (3 nested, finest last):**
- order-join/SCC atoms = STRONG: 1,1,6,35,353 (char_A/H/R/zeta factor over them)
- modular/substitution atoms = MODULAR-PRIME seeds: 1,1,1,0,3,15,197 -- I COMPUTED n=7=197, completing opus THM-1960's open census (modular-primes subset strong for n>=3)
- circulant character-generated: 1,0,1,0,2 (Paley Gauss / interval Chebyshev, Re=-1/2)

**THE VANISHING-REACHABLE-FRACTION LAW (the quantitative core):** strong-frac 0.25,0.5,0.625,0.774,0.873 (n=4..8); modprime-frac 0,0.25,0.268,0.432 (jumps at n=7); asymmetric 0.875 at n=7. Reducible (order-join-collapsible) is only 12.7% at n=7, 8.7% at n=8 -> the reduction principles reach an ASYMPTOTICALLY NULL SET. That IS why n=7 forces enumeration: the irreducible interior swells to full measure. = the honest form of computational irreducibility / the apex-7 (LRC 14=2*7) wall.

**SURVIVING ISLANDS (clean large-n):** Paley T_p doubly-regular self-comp H-max, c3=p(p^2-1)/24, |lambda|^2=(p+1)/4; every regular tournament on m shares c3=m(m^2-1)/24; circulant iso counts 2,4,4,6,16,16,30 (n=7..19) all on Re=-1/2.

**Integrated:** opus THM-1960 (I completed their seed census to n=7), death-star S81 (recursion=order-join, convergent with my S196), mac-mini THM-1936 (signed R), klein THM-1950 (H>=disc strong base), kps THM-1880 (char_S). No collisions; claimed THM-1978 (>max 1965), HYP-8731.

**Honest scope:** the modular-prime n=7=197 census + fraction law + circulant formulas are VERIFIED exact this session; the "breaks at n=7" catalog synthesizes cited theorems (each verified in its own source). The unifying observation (common THM-1830 driver + vanishing-fraction) is the contribution.

**Next:** (1) modprime seed census n=8 (needs iso n=8 = 6880, hard); (2) does the prime-fraction have a clean asymptotic (1 - O(1/2^n)?); (3) which surviving-island invariants extend the reduction into the sea. Artifacts: THM-1978, HYP-8731, reflection the-n-ge-7-regime-what-breaks-what-survives-boxeph-S197.md, scripts modular_prime_census_n7 + large_n_circulant_patterns _boxeph_S197.py (+.out).

