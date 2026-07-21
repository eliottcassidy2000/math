## codex-2026-07-21 -- GMC(2) span-6 residual made exact: E6 = 466560*(ad)^3 on the B2 branch (HYP-8580 addendum)

Owner directive: synthesize past and incoming work, pull often, and work creatively toward the two-dimensional Gaussian Moment conjecture.  Pulled/rebased over incoming HYP-8760/THM-2010 sequence-valued tournament invariant work; it did not directly touch GMC(2), but its value-count/metagraph viewpoint informed the required Tournament Analysis as a proof-obligation tournament rather than a runner/charge tournament.

- **Concrete finite lemma progress.** Strengthened death-star-S73's span-6 `{+3,+1,-1,-3}` constant-coefficient stratum from a numerical B2 residual probe to an exact symbolic certificate.  For `P=aZ^3+bZ+cZbar+dZbar^3`, `E2=0` gives `bc=-6ad`.  In the all-nonzero branch set `u=ad`, `x=ac^3`, `y=b^3d`; exact reduction gives `E4=24*(x+y+54u^2)` and `E6=38880u(x+y)+2566080u^3`, hence under `E4=0`, `E6=466560u^3`, nonzero if `a*d!=0`.  The boundary cases `a=0` and `d=0` are already killed by `E2,E4`.  Therefore moments `2,4,6` rule out every two-sided nullcone member in this stratum.
- **Status discipline.** This upgrades HYP-8580's local proof quality but does **not** revive the false span-only global cutoff: THM-1770/1790 still force radial-degree detection depth to grow.  The usable synthesis is a finite lemma bank for constant/radial-bounded strata plus the Hermite/Sheffer/no-common-root route for the analytic bridge.
- **Artifacts.** `04-computation/gmc2_span6_symbolic_residual_codex_20260721.py` and `05-knowledge/results/gmc2_span6_symbolic_residual_codex_20260721.out` (source/output SHA-256 `b77300fe...65a15` / `3a677001...840e0`).  Results index updated.

## kind-pasteur-2026-07-21-S128c145 - More sequence-like invariants: ~10 candidate-NEW tournament-invariant sequences + OEIS ID's (THM-2010)

Owner: look for more sequence-like invariants. Hunted NEW integer sequences from our objects,
differentiated from the census/reciprocal sweep (opus THM-1985/1990) by targeting three
under-explored families computed exhaustively over iso classes n=3..6 (exact arithmetic): value-counts
|X|(n)=#distinct values of invariant X, extremal maxX(n), structural counts. Direct OEIS lookup (curl;
WebFetch was 403-blocked, curl works).

**KNOWN/identified:** iso=A000568; |score|=A000571; #strong=A051337; #regular=A007079; max H (Ham
paths)=3,5,15,45=A038375 (Moon 1968, 'max spanning paths', hard); #self-converse=2,2,8,12=A002785;
|specS| skew/Seidel char-poly count=1,2,2,6=Breen-Stover-Yates (arXiv 2406.09697, extends 11,50).

**CANDIDATE-NEW (no OEIS tournament match on 4 terms; [0]=ZERO hits at all):**
- |specA| adjacency char-poly count = 2,3,9,28  (the A-SIDE TWIN of the KNOWN skew |specS|; strictly
  finer, 28 vs 6 at n=6 -- adjacency keeps the score/hierarchy data the skew form discards; the
  obvious missing companion sequence)
- |H| Redei-spectrum size = 2,3,7,19; |cyc| = 2,3,9,32; |R| = 2,2,6,8; |disc|=|var(lam2)| = 1,2,2,5
  (equal, both skew-spectral)
- |arb_inv| = 2,4,12,55 [0] (near-complete invariant, 55 of 56 iso classes at n=6)
- max-total-arborescences = 3,10,55,333 [0]; max|R| = 3,3,15,15,147 [0]
- metagraph edges |E(G_n)| = 1,5,30,290 [0]; metagraph 1-WL colors = 1,2,10,34

**Observations:** (1) skew/adjacency ASYMMETRY -- the skew (Seidel) char-poly count is published, the
adjacency one is not; |specA| is the missing A-side sequence. (2) |disc|=|var(lam2)| take equally
many values (both skew-spectral). CAVEAT: 4 terms is thin; the four [0] are most confidently new.

Reference catalog (verified data + external ID's), not a proof; value = ~10 uncataloged
tournament-invariant sequences surfaced + the skew/adjacency gap, homed for OEIS submission +
reciprocal-signature study. NEXT: extend the [0] sequences to n=7 and submit to OEIS; closed form for
metagraph edges 1,5,30,290; reciprocal signatures (THM-1990). Files THM-2010; HYP-8760; script
new_sequence_invariants_kps_S128c145 (+out). Namespace THM-2010/HYP-8760 (max was 2000/8751; margin).

## death-star-2026-07-21-S84 -- H≥disc crux REDUCED to quasirandomness (binding case = Paley = quasirandom) + a rigorous bosonic≥fermionic positivity E[sym²]≥E[alt²]=E[|P|²]≥0 for GMC(2). HYP-8699.

**Owner directive:** keep reducing the H≥disc crux; pull often; also work GMC(2).

- **PART 1 — CRUX further reduced (S83's (i) H(reg)≥Szele avg → quasirandomness).** Binding case = DOUBLY-REGULAR (Paley, max disc): smallest H/(n·disc) among regulars (Paley-7: 3.38 vs rot 25; Paley-11: 35.6 vs 8457). Paley is QUASIRANDOM (K-spectrum ±i√n, ratio →0) ⟹ H=(1+o(1))·avg (quasirandom Ham-path counting lemma); measured H/avg≈2.0-2.4 (bounded). n·disc/avg=n(n+1)^{(n-1)/2}/n!→0 super-exp (0.71 at n=7, 0.069 at n=11), so large n is loose (needs only a tiny fraction of avg); small n direct. Min-strong (Busch) route EXCLUDED — doubly-regular disc ~(√n/2)^n is too big — so the crux genuinely needs regular=quasirandom=near-average. The crux is now standard pseudorandomness, not an eigenvalue-product mystery.
- **PART 2 — GMC(2) positivity.** The S81 Pell identity sharpens: for REAL-coeff P, E[sym(P)²]−E[alt(P)²]=E[P·P̃]=E[|P|²]≥0 (Bargmann norm) ⟹ E[sym²]≥E[alt²] — a RIGOROUS proof of klein THM-1810's bosonic≥fermionic at the squared-moment level. HONEST: orthogonal to the nullcone (one-sided P=Z has E[sym(P^m)²]=m!/2>0 despite E[P^m]=0); a Bargmann-PD handle on the TORAL side (S67/S77 toolkit), the open RADIAL gap unaffected — does NOT close GMC(2).
- No open problem closed; a genuine crux reduction + a small rigorous positivity. Credits klein THM-1950/1810. reflection crux-reduced-to-quasirandomness-...-S84; script crux_reduction_and_gmc2_positivity_S84 (+out). LRC(≤13) not re-audited.
## opus-2026-07-20-S448 - Folded the sequence-sweep findings into THM-1985: the LRC deep-well = harmonic/triangular (THM-805 bridge)

Folded the S447 background sequence-sweep (28-sequence master table + the repo's existing harmonic/
gamma/zeta appearances) into THM-1985 s4. THE KEY SYNTHESIS (verified exact):
- THM-805: the LRC deep-well base measure = a HARMONIC number over a TRIANGULAR number:
  m({1..k}; 1/(k+2)) = H_k / C(k+2,2). At k=12: H_12/C(14,2) = 6617/194040 (exact, = mac-mini
  THM-736's |G'({1..12})|). So the divergent-edge object (harmonic H_k) and the convergent figurate
  size (triangular C(k+2,2), sum 1/·=2 of THM-1985 s1) combine into the single most important LRC
  extremal measure -- the reciprocal-sum program meets the LRC flagship exactly.
- The LONELINESS FLOOR spectrum IS the harmonic series: M_floor(n)=1/n => sum = sum 1/n (diverges).
- THM-1926 tournament zeta = 1/det(I-uA) = Euler product over primitive cycles (cycles = primes).
- Staircase 3-cycles k(k-1) telescope to 1. Deep-well denominators Phi_6(n)=n^2-n+1 (cyclotomic),
  sum 1/Phi6=1.798 converges.
- REFINEMENT (agent): only pure binomials C(n,k) give rational k/(k-1); the 1+2^d / 2^{n-2}+1 /
  Cayley-Dickson 2^k+1 families converge to IRRATIONAL Erdos-Borwein-type constants. A002088, A001764
  absent from repo.

THE LOOP CLOSES: the figurate 2 (triangular=tournament size), the harmonic H_k (resistance/divergent
edge), and the zeta_T Euler product are three faces of one structure; the LRC deep-well = H_k/(triangular).

Files: THM-1985 s4 + depends_on (added THM-805/1926); harmonic_triangular_bridge_opus_S448.py(out).
Cites THM-805/1926/1370/1920/1930/1970/1975, kps THM-1980.

## opus-2026-07-20-S447 - Reciprocal sums = the harmonic-scale face of the poly/#P tower (THM-1985)

Owner: the reciprocal of an integer sequence is a subset of the harmonic numbers; study sum 1/a_n
for as many repo sequences as possible; figurate reciprocals (triangular=2), Abel-Dini, Bertrand;
1+1/2+..+1/5 > 2 already while sum 1/triangular = 2.

THM-1985: a sequence's GROWTH (its poly-tower position) = its reciprocal sum's CONVERGENCE. THREE
STRATA. (1) FIGURATE invariant-SIZES = char_S coefficients (THM-1920): sum_{n>=k} 1/C(n,k)=k/(k-1)
(exact k=2..6). Tournament sizes -> RATIONALS: arc=C(n,2)=triangular => sum 1/arc = 2 (the Downey-
Ong-Sellers triangular identity realized on the tournament -- the char_S subleading series sums to
exactly 2); c3max=C(n,3)->3/2; var-max=2C(n,3)->3/4. Degree-k invariant -> reciprocal-sum k/(k-1).
(2) COUNTING seqs (super-exp) -> fast transcendentals: sum1/A000568=2.8535, sum1/A038375=2.6293,
sum1/A051337=2.198, sum1/A002854=1.062; Cayley-Dickson sum1/(2^k+1)=0.7645 (Erdos-Borwein cousin),
H=1+2^(n-2) SC-neighbor -> 1.2645. (3) H-VALUE spectrum (odds minus {7,21}, THM-1370) ~linear =>
sum 1/H-value DIVERGES (harmonic-slow) = H's VALUE SET sits at the convergence/divergence boundary =
the reciprocal-sum face of THM-1970's formula/#P edge. ABEL-DINI: no series at the exact boundary =
kps THM-1980's 'Redei parity is the last formula'. BERTRAND boundary = sum 1/(n ln n).

THE PICTURE: rational k/(k-1) [figurate invariant sizes, poly, deep convergence] | transcendental
[counting sequences, the census] | DIVERGES [H-value spectrum, #P, the edge]. The reciprocal sum is
the harmonic-scale coordinate that recovers the poly/#P tower.

Concurrent: a background agent is sweeping the full repo sequence list (30+ OEIS A-numbers, growth
rates, existing harmonic/gamma appearances like THM-805 resistance=harmonic-number) -- findings fold
into THM-1985 next session.

OPEN: identify the counting constants (2.85, 2.63 -- e/pi/new?); the H-value density c in c*ln x;
Bertrand-scale repo sequences.

Files: THM-1985; HYP-8745; reciprocal_sums_of_repo_sequences_opus_S447.py (+out). Namespace clean
(1985/8745). Cites THM-1920/1930/1970/1975/1370, kps THM-1980, Downey-Ong-Sellers.

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
## boxeph-2026-07-21-S198 -- THM-1979 TOURNAMENT SPACE IS A SPECTRUM (single point -> continuum)

**Owner:** understand tournament space on n vertices as a spectrum from a single point (transitive) to a continuum; the maximally-spread score classes house the different structure.

**THE FRAME (verified n<=7).** Tournament space fibers over the SCORE SEQUENCE (Landau polytope; counts 2,4,9,22,59 for n=3..7). The spectral coordinate is score spread sigma^2=Var(scores) in [0,(n^2-1)/12]. Cyclicity is its EXACT affine image: c3 = n(n^2-1)/24 - (n/2)sigma^2 (= the classical c3=C(n,3)-sum C(s_i,2) restated). So score-spread and cyclicity are ONE axis, opposite directions.

- SINGLE POINT = TRANSITIVE (sigma^2 max, scores 0..n-1): fiber=1, c3=0, char_A=x^n (nullcone vertex), zeta=1 (my THM-1926, no periodic structure), reducible. The ordered structureless pole.
- CONTINUUM = REGULAR/near-regular (sigma^2->0): fibers SWELL -- max fiber 1,1,3,12,47 (n=3..7), -> inf with n; at low sigma^2 every class is STRONG and mostly MODULAR-PRIME; c3 max = n(n^2-1)/24. This is where the different structure lives (circulant/Paley thread THM-1955, |R|-independence mac-mini THM-1966 first at n=7, the whole irreducible interior THM-1978).
- MONOTONE LAW: fiber size, strong-frac, modprime-frac all run OPPOSITE to sigma^2. High-spread fibers are singleton reducible chains (strong=0); the regular center is all-strong. The structurally-richest score class is NEAR (not exactly at) the center: n=7 peak fiber 47 at sigma^2=4/7, score (2,2,3,3,3,4,4), all 47 strong, 29 modular-prime.
- LIMIT: the tournamenton spectrum -- transitive W=1_{x>y} (single ordered point) to quasirandom W=1/2 (positive-entropy continuum); degree function d(x): x -> 1/2; sigma^2 = integral (d-1/2)^2.

**UNIFIES the recent arc:** reduction principles = statements about the TOP (transitive/high-sigma^2, reducible); hardness = the BOTTOM (regular/low-sigma^2, irreducible). n=7 wall = first crack of the point opening into the continuum (THM-1830 atom). The clean theorems live at the point; the mathematics lives in the continuum.

**Housekeeping:** opus-S445 first-pushed THM-1970 (and THM-1975) so I renumbered my S197 n>=7-regime theorem THM-1970 -> THM-1978 (off the round-number grid to stop the collision churn). This session: THM-1979, HYP-8732. Integrated mac-mini THM-1966 (|R| new invariant from n=7). No overrides.

**Honest scope:** the fibration, the affine c3-sigma^2 law, the fiber/strong/modprime monotonicity, and the max-fiber sequence 1,1,3,12,47 are VERIFIED exact n<=7. The affine law is a restatement of a classical identity; the contribution is the SPECTRUM framing (single point <-> continuum via score spread) that unifies reducibility, cyclicity, structural diversity, and the n=7 breaks on one axis.

**Next:** (1) the max-fiber growth 1,1,3,12,47 -- asymptotics (does it track the # of tournaments with the balanced score seq, ~ 2^{C(n,2)}/(Landau count)?); (2) the exact location of the structural-diversity peak vs the regular center (off-center at n=7); (3) the quasirandom-neighborhood entropy as the continuum's "dimension". Artifacts: THM-1979, HYP-8732, reflection tournament-space-as-a-spectrum-single-point-to-continuum-boxeph-S198.md, script tournament_space_spectrum_boxeph_S198.py (+.out).

