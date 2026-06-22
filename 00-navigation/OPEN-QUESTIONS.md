# Open Questions

**OPEN-Q-108 S123 apex-lock Steinhaus sequence split:** HYP-2917 refines the HYP-2913/HYP-2914/HYP-2915 three-gap/Dirichlet/no-spectral-gap thread into a finite sequence filter plus a global off-apex escape theorem.  Put the observer and the thirteen runner residues on the denominator-14 clock for each unit `a`, retaining collisions as zero gaps.  AP has gap partition `((1,14))`; GW has the one-collision/one-missing partition `((0,1),(1,12),(2,1))`.  But `{1,...,11,13,36}` has the same coarse GW partition and is still globally loose, with exact `M=3/41`.  Therefore the missing theorem is not a residue-profile classification.  The sharp Node-2 target is: every non-AP/GW apex-locked reduced row has an off-apex escape `M(S)>1/14`; covering/apex-blocked rows still require the comb/Weyl machinery because multiples of `14` kill the apex clock.  Exact support: single AP replacements `v<=80` have unique tight sets AP/GW and no below-threshold rows; local two-swaps `<=18` have no non-AP tight rows; q-covering rows in `[1,19]` have min `M=1/12`.

**OPEN-Q-108 S122 apex-majority gamma descent:** THM-571 closes the post-THM-568 `|M14|>=7` branch, modulo the accepted below-frontier LRC input.  The S121 two-halfstep collision is real for a raw fourteen-shift sieve, but it disappears in the actual apex-majority branch: if two residual speeds are divisible by `7` but not `14`, then at least nine speeds are multiples of `7`; scale by `7`, use LRC<=13 to get a strict safe `7`-phase, and pigeonhole over seven lifts.  At most four speeds remain nonmultiples of `7`, each forbidding at most one lift, so one lift survives.  The live OPEN-Q-108 proof burden is now the `|M14|<=6` side as packaged through S31v: the comb-teeth union bound plus the bounded/intermediate finite-core census and scale-separation reduction.  The recent LRC proof through 13 total runners supplies the below-frontier margin input; it does not by itself classify 13 moving-speed tight loci.

**OPEN-Q-108 S121 apex-majority shift guard / half-step collision guardrail:** THM-570 closes a real piece of the post-THM-568 residual.  For `S=14Q union R` with `|Q|>=7`, `|R|<=6`, the below-14 theorem gives a strict `Q`-safe interval in the scaled coordinate.  The fourteen lifts `t=(u+k)/14` keep `14Q` safe.  Exact shift facts: ordinary residual speeds (`gcd(r,14)!=7`) forbid at most two shifts and at most one in each parity; a single half-step speed (`gcd(r,14)=7`) forbids at most one whole parity class.  Hence the branch is LRC14-safe if `R` has at most one half-step speed.  The sharp raw-sieve obstruction is HYP-2911: at least two residual speeds divisible by `7` but not `14` can phase-cover both parities (exact guardrail `r=7`, `r=161`, `u=2/49`).  S122/THM-571 resolves this obstruction for the actual apex-majority branch by descending to the `7`-phase.

**OPEN-Q-108 S120 Lean q=14 boundary + S31ab covering-strictness signal:** THM-569 now formalizes the exact denominator-14 unit-grid split in Lean: for every unit `a mod 14`, `Lonely 14 v (a/14) <-> forall i, not (14 | v_i)`, with named residues `1,3,5,9,11,13`, plus `no lonely time -> some 14 | v_i`.  This makes the finite apex side of THM-523/THM-568 audit-level rather than script-level.  A post-rebase KPS-S31ab script claims the 14-covering residual is not tight and verifies AP/GW replacement families with minimum `M=1/13`; treat it as direct signal for the next formal theorem, not yet as a Lean-closed endpoint.  The sharp target is now: from `S=R union M14`, `R` 14-free, extract a genuine open `1/13`-margin carrier from the smaller-runner theorem and prove the danger combs of the `14`-multiples cannot cover it at threshold `1/14`.  This is the covering-strictness theorem that would turn THM-568 + THM-569 into the tightness-star closure.

**OPEN-Q-108 S119 THM-079 template / STAR obligation, updated by THM-568 and THM-571:** HYP-2909 audits the latest "LRC14 is THM-079" proof template.  The analogy is structurally correct but the final atom theorem must be stated sharply.  Move A is HYP-2905/HYP-2906 boundary-state peeling (with THM-524/527/565 and HYP-+2878 as existing reduction machinery).  Exact audit `lrc14_thm079_star_audit_codex_s119.py` shows why raw apex blocking is not enough: AP and GW are tight and non-covering, but covering rows `{1,...,11,13,84m}` block every denom-14 unit point while remaining safely lonely off-apex (`7/89`, `14/173`, `42/509`).  Incoming THM-568 proves the structural half of `(star)`: primitive tight rows have apex denominator `D=14` (nonprimitive `D=14*gcd(S)`).  THM-571 closes the `|M14|>=7` apex-majority covering residual by gamma descent.  The remaining covering-strictness work is now the `|M14|<=6` scale-separated/finite-core package, plus any bounded tight-locus/compression theorem needed to feed it.  Alternate closures remain STAR+ tight-locus/three-gap/Goddyn-Wong classification or HYP-2908's state lift from the reduced atom to a tournament-conflict connected `I(.,2)=7` packet, hence the forbidden `K_3`.

**OPEN-Q-108 S118 digraph H=7 realizability guardrail:** HYP-2907 sharpens the KPS-S31y forbidden-clique route.  The facts "`H=7` is forbidden for tournaments" and "arcs have two states" do not apply to every binary relation: exact audit shows ordinary binary ordered-arc digraphs realize `H=7` on `n=4`, and incomplete oriented graphs realize `H=7` on `n=5` (`1440` examples).  Tournaments through `n=6` have `H=7` count `0` and all `H` odd, while tie-free AP7 winding-tournament cells also avoid `H=7` and have `H(1/7)=175`; wall-time samples are not tournaments and can give even/zero counts.  Incoming S48 identifies the AP14 boundary packet as seven tied diameter comparisons under the antipodal order-2 symmetry, but those ties still have to be resolved into tournament/OCF data.  Incoming S31z adds the logical-status guardrail: for the Pi^0_1 LRC14 statement, "impossible to disprove" means "true," so this is a proof route, not an independence route.  Therefore the missing LRC14 theorem is a realizability statement: `apex-7 over-cover -> tournament OCF conflict graph Omega=K_3` (or an equivalent labelled packet).  Once that is proved, THM-200/THM-343 block the counterexample; without it, the generic digraph models are counterexamples to the slogan.
**OPEN-Q-108 S118 forbidden-H7 state-lift addendum:** HYP-2908 extends HYP-2907's digraph guardrail into a precise transfer theorem.  Exact graph atom: connected `I(G,2)=7` forces `G=K_3`; by THM-002/THM-343/THM-201, that atom cannot occur as a tournament odd-cycle conflict component.  But the transfer must land in that category: S118 finds a 4-vertex arbitrary present/absent digraph with exactly 7 Hamiltonian paths, and THM-344 shows `K_3` subgraphs are allowed inside larger complete conflict graphs (`H=63`, `Omega=K31`).  Therefore the missing LRC theorem is not "make a digraph"; it is a state lift from the primitive top-balanced bounded core to a tournament-conflict-realizable connected binary packet graph with activity two and `I(.,2)=7`.  Candidate LRC vertices are HYP-2648 measured sector-state words, HYP-2691 sector-transfer packets, HYP-2677 packet-sign atlas states, cover-arc packets, exact-period phi packets, and HYP-2632 support-six relation packets.  If that lift is proved after the HYP-2906/HYP-2905 reductions, an LRC14 counterexample would have to realize the forbidden `K_3` atom, so the disproof is impossible.
**OPEN-Q-108 S119 tightness-star exact-atlas companion:** HYP-2910 supports HYP-2909 with exact AP/GW and q-covering-window checks.  The audit verifies the denominator-14 floor (`14|v*k <=> 14|v` for the six unit residues), AP and GW both have `M=1/14` with identical denominator-14 argmaxes and binding pairs, and the AP single-swap atlas through `v<=80` has exactly one non-AP tight row: GW.  The finite q-covering window `[1,18]` has `966` exact rows, minimum `M=1/12`, and no tight or below-threshold row.  Incoming THM-568 proves the apex-denominator half: every tight optimum has denominator `14*gcd(S)`, so primitive tight rows optimize at denominator `14`.  THM-571 closes the `|M14|>=7` covering side.  THM-523 handles 14-free rows.  The remaining OPEN-Q-108 form is now the `|M14|<=6` scale-separated/finite-core package: combine S31v's comb-teeth union bound with a rigorous bounded/intermediate census or compression theorem.

**OPEN-Q-108 S117 boundary-state induction addendum:** HYP-2905 imports the tournament-induction lesson into the LRC proof order.  Strong-ear tournament growth is exact only with the boundary state `(start,end,Q)`, not raw vertex deletion; exact audit through labelled strong parents `n<=5` gives `0` insertion-formula and `0` strongness failures.  LRC remove-large is the same grammar with safe-set boundary state `(mu,components)`, extended by S31v/S31w to arc budgets and scale-hierarchy peeling.  Incoming S31x adds that scale-separated safe measures multiply like tournament strong components, so the effective theorem should control the product error by the retained boundary state.  Incoming S47 sharpens the irreducible endpoint to the Mode-A peel atom `{consec,GW}`; the remaining LRC step is the H=21/Moon analogue, proving the bounded tight locus or positive slack.  The LRC14 proof should therefore be stated as boundary-state induction: omit-prime direct witness; remove-large descent to smaller LRC seeds; `r<=6` multi-large union bound; `r>=7` second-moment resonant-pair/divisibility defect; bounded `{consec,GW}` tight-locus theorem plus missing-depth parity Newton packets and possible resonance-lattice deletion-contraction.  This is the sharpest useful induction: it kills unbounded and non-covering rows, but bounded covering cores remain the finite Node-2 extremality base.

**OPEN-Q-108 S117b one-large interval-peeler addendum:** HYP-2906 sharpens the scale-separated proof tree before the full component-budget machinery is needed.  If a seed with max speed `m` has one witness margin `alpha>1/n`, then it stays threshold-`1/n` safe on a connected interval of length `2(alpha-1/n)/m`; an added runner's danger teeth have length `2/(nv)`, so `v>m/(n(alpha-1/n))` forces a witness.  Taking `alpha=1/(n-1)` from `LRC(n-1)` gives the clean gate `v>(n-1)m`; for LRC14, a largest speed greater than `13` times the second largest is automatically safe by LRC13.  The AP-core `{1,...,11,13}` has explicit `tau=1/12`, so its local gate is `v>78`, certifying committed lcm/radical speeds without equidistribution.  Thus any counterexample reaching the hard p0/depth-parity/Node-2 machinery must first be top-balanced (`v_max<=13v_second`) or multi-large with no locally peelable top speed.  HYP-2904 remains the right density and multi-large object; HYP-2906 is the existence-first one-speed peeler under HYP-2905's boundary-state switchboard.

**OPEN-Q-108 S116 scale-separated induction addendum:** HYP-2904 gives a rigorous smaller-size induction branch, but only with a carried topology budget.  If `A=Safe_n(B)` has measure `mu>0` and `c` interval components, then adding speed `v` leaves measure at least `(1-2/n)mu-2c/v`, because the new unsafe set is an exact density-`2/n` comb and each component pays at most two boundary partial periods.  Hence every fixed seed certified by `LRC(n-1)` becomes safe for all sufficiently large added speeds at threshold `1/n`.  For the AP-core seed `{1,...,11,13}` at `n=14`, exact audit gives `mu=426/35035`, `c=4`, and all `v>=768` are certified, including radical/lcm committed speeds `30030`, `60060`, and `510510`.  Incoming KPS-S31v is the matching multi-large lemma: the same comb-teeth estimate plus union bound closes `r<=6` large speeds over a bounded core; `r>=7` is now the second-moment / resonant-pair defect problem.  Incoming KPS-S31w organizes the proof tree: remove-large peels the scale hierarchy to smaller proven LRC seeds, omit-prime gives a direct `t=1/p` witness, and dilation normalizes; the only non-descending base is the bounded covering core.  The guardrail is equally important: dilation preserves `mu` but multiplies `c`, so no runner-count-only induction can be uniform.  The live reduction target is now: bounded/scale-normalized seeds to Node 2/AP-three-gap/Legendre-Venn atlases, large committed speeds to Node 3 finite-comb/exact-period Weyl estimates with an explicit component or packet budget, and the remaining large-speed obstruction to a bounded divisibility/relation-pair ledger.

**OPEN-Q-108 S114/S115 corrected three-mode composition + lcm denominator wall:** HYP-2901 integrates the owner's Legendre correction with incoming KPS S31s/S31t, mac-mini S45, and the committed-speed refutation.  The odd Legendre mode is the full Venn `A+B+D-C-E-F+G`: corners `A,B` have size `N-1`, corner `D` has size `N-2`, overlaps `C,E,F` have sizes `N-2,N-3,N-3`, and `G` has size `N-4`; the edges are `A+B-C`, `A+D-E`, `B+D-F`, and the center is `A+B+D-C-E-F+G`.  Mobius is always-on, Eisenstein is even-only `A+B-C`, and Legendre is odd-only full Venn; the letters have different subtournament sizes in each mode.  The lcm family `S_X={1,...,11,13,lcm(2,...,X)}` gives a theorem-level guardrail: every denominator `D<=X` is killed by the committed speed, so no fixed finite denominator basis can prove LRC14.  The stronger `firstD=nextprime(X)` pattern is false (`X=60` firstD `67`; `X=110,120` firstD `121=11^2`; `72/127` nextprime mismatches over `X=14..140`).  S45 adds the radical filter: if a prime `p<=13` divides no runner, then `t=1/p` is already safe, so hard rows must be prime-covering for `2,3,5,7,11,13` and kill `14`; S114 adds that exact first witnesses require prime-power packets, not primes alone.  S31t adds the wide-cap subtarget, but HYP-2903 corrects its scope: the universal Bonferroni-3 claim is false (`B={0,1,2,3}`, `F={16,19,22,25,28}` has `T_{>=4}>0`), KPS S31u adds high-depth spread-far failures with `T1=T2=T3=0` but `p0<<cap`, and S115 shows even edge-active rows can have positive third-order tails.  The live wide target is a missing-depth parity guard in the binding leg, while high-depth/spread-far rows route to Node-3 slack, equidistribution, or finite residual atlases.  Proof split: Node 2 remains finite/AP-three-gap/Legendre-Venn extremality with depth-labelled Newton packets; Node 3 and finite Part A require analytic exact-period prime-power/residue equidistribution beyond the lcm wall.

**OPEN-Q-108 S115 depth-parity correction to the S31t far-packet target:** HYP-2903 now corrects the Bonferroni-3 subtarget.  Exact common-refinement integration gives the pointwise formula `T_r(x)=(-1)^(r+d(x))C_{d,r}(x)`, where `d(x)` is the number of base-missed inner sectors and `C_{d,r}` counts `r`-far subsets covering those missing sectors.  Therefore the `r>=4` tail is a missing-depth parity ledger, not a Venn-containment sign.  The raw Bonferroni-3 upper bound fails in exact examples, including k=8 `B={0,1,2,3}`, `F={16,28,29,32}` with tail `19/68208>0`, and even an edge-active k=9 row `B={0,1,7,10,13}`, `F={15,23,24,31}` with `T2>0` but tail `307/598920>0`.  The live theorem is now: positive even-depth high packets must be bounded by negative odd-depth high packets plus cap slack; positive-tail activation/depth classes route to doublet/triple/decorrelation finite atlases rather than a universal third-order truncation.

**OPEN-Q-108 S110 product-Mobius packet ledger:** HYP-2899 joins the owner's coprime-density/totient prompt to the three tournament tiling recurrences.  The divisor axis is exact: the copy rule `sum_{d|n}c(d)=n` has Mobius inverse `c(n)=phi(n)=sum_{d|n}mu(d)n/d`, and HYP-2856's witness-floor constants are the Farey/totient limits `sum_{q<=Q}phi(q)~3Q^2/pi^2` and `sum phi(q)/q~6Q/pi^2`.  The tiling/character axis is also Mobius-labelled: full tiling `A+B+C-D-E-F+G` is Boolean `B3`, even half `A+B-C` is `B2`, and odd half `A+B-C+D-E-F+G` has two simultaneous addresses.  Incoming S31q reads the prompt-order sign string `++-+--+` over `A..G=1..7` as the Legendre `chi_7` split with zero/triple slot positive, while S110 reads the half-tiling corner order `A,B,D` as Boolean `B3` with size/depth pushforward `(2,0,-2,1)`.  New proof target: keep every residual packet on the product ledger `Div(D) x B_r` before scalarization, so denominator capacity `phi(D)`, CRT multiplicativity defects, character labels, and one/two/three-far Boolean signs are only multiplied after their labels are retained.  S109's `w=84m` one-tail branch is a concrete model: killing q-witnesses does not remove the coprime floor, it moves the binding witness to unit denominators `D=84m+5` coprime to `2,3,7`.  Incoming HYP-2898/S111 points the same way from small even q: scalar additive energy fails, but labelled Fejer/residual control survives.  Incoming S44 adds the denominator-killing form: a speed `s` kills all `phi(b)` primitive Farey points for each `b|s`, so the small-denominator covering core is a `Phi(14)=64` totient-weighted survival lattice, not 13 independent scalar targets.  The missing theorem is a product-lattice residual bound: coherent low-depth atoms go to finite AP/Freiman/packet atlases, while incoherent high-denominator unit packets go to THM-566/HYP-2890 decorrelation rather than a fixed finite denominator basis.
**OPEN-Q-108 S113 totient-curvature update:** HYP-2900 refines HYP-2899 by showing the full/even/odd tournament recurrences are exact cell-address operators, not multiplicative-function recurrences.  Applying them to `phi`, `phi/n`, and AP endpoint density with the exact subtournament sizes exposes a nonzero Euler-factor curvature.  The LRC14 boundary is especially diagnostic: even-half `n=14` compares two size-13 prime carriers with size-12 `2^2*3`, while actual size 14 is `2*7`; the exact `rho` residual is `-296/273` with curvature `{2:3,3:1,7:1,13:-2}`.  Incoming S31q/S44/S31r supply the companion coordinates: the sign words are Mobius / `chi_7` / `chi_3` channels, resonance killing is totient-weighted with `Phi(14)=64`, and `14=2*7` is the even Eisenstein fold into the odd Legendre apex.  Thus coprime density enters OPEN-Q-108 through exact-period `phi` packets and their character-labelled CRT/chi7/coimage curvature, not through a scalar recurrence.  This reinforces the labelled Fejer/signed-current route after exact-tiler and one-tail branches are routed away.

**OPEN-Q-108 S109 zeta `-1/12` / one-tail resonance-killing closure:** HYP-2896 converts the owner's `M({1..11,13})=1/12` and `1+2+3+...=zeta(-1)=-1/12` prompt into an exact finite/discrete proof fragment.  Let `C={1,...,11,13}`.  Every one-tail row `C∪{w}` is LRC14-safe: if `12∤w`, the q=12 witness survives (`M>=1/12`); if `12|w` but `14∤w`, the q=14 witness survives (`M>=1/14`); if `w=84m`, then the covering row has exact witness `t=(35m+2)/(84m+5)` and exact binding-pair value `M=7m/(84m+5)>1/14`.  This complements KPS S31p's resonance-killing game and adds a guardrail: in the covering branch the value is a binding-pair affine denominator, not merely `1/(smallest surviving b)`.  The one-tail disproof branch is closed; the remaining disproof/proof battlefield is multi-large or moderate resonant covering rows, where several divisibility killers interact before equidistribution and HYP-2890/HYP-+2878 support-six residual cancellation take over.

**OPEN-Q-108 S108/S108b sub-14 covering/tiler training atlas:** HYP-2895 applies the current LRC14 tools to `N<14` and keeps `N=14` as a contrast line.  AP rows are tight, Goddyn-Wong acceleration atoms appear at `N=8` and `N=14`, and S108b's exact single-swap boundary census finds non-AP tight rows only at `N=5` (`2->7`), `N=6` (`2->9`), `N=8` (`6->12`), and `N=14` (`12->24`); all four have safe measure `0` and q-deficit exactly `(N,)`, so they are apex-denominator boundary tilers rather than covering rows.  Mac-mini S42's broader small-`n` search reports additional bounded sporadics outside this window, but with the same usable boundary condition: primitive exact tilers found avoid multiples of `N`, so `t=1/N` witnesses them.  AP-drop q-covering repair rows are all loose once THM-523's easy q-witness is disabled; for `N=9..14`, the best AP-drop repair is `drop N-1, add N(N-1)`, active pair `(1,N(N-1))`, `M=N/(N(N-1)+1)`.  KPS S31o then splits the covering crux into bounded compactness/margin, one-large-speed equidistribution, and a moderate/resonant middle.  New proof target: quotient out the apex q-witness boundary, discharge the bounded/unbounded covering extremes, then close the moderate AP-facing q-covering residual by THM-524 binding margins or HYP-2890/Clebsch-Bruhat-octahedral support-six cancellation.

**OPEN-Q-108 S106 Goddyn-Wong sporadic-tiler classification:** HYP-2893 reframes the Goddyn-Wong exact tilers as AP tail accelerations controlled by Jacobsthal/nonunit intervals.  Starting from `{1,...,n}`, replacing `v` by `2v` is tight when every integer in `[n-v+1,2n-2v+1]` has nontrivial gcd with `v`; LRC14 is the `n=13,v=12` case, producing `{1,...,11,13,24}` from the nonunit window `[2,3]`.  This complements THM-560: difference-closed exact tilers are AP dilates, while Goddyn-Wong rows are accelerated-tail atoms.  New subtarget: prove a finite/recursive classification of exact tilers into AP dilates plus Jacobsthal acceleration atoms, then show all remaining non-difference-closed rows have a positive safe interval and can feed the HYP-2890 residual/cap machinery rather than the tiler boundary.

**OPEN-Q-108 S105 Clebsch/Bruhat residual-carrier refinement:** HYP-2891 converts the unital/Clebsch/truncated-octahedral prompts into a finite labelled-residual target.  Clebsch closed neighborhoods give a pair-balanced `2-(16,6,2)` design on the folded 5-cube; the tangent-sector quotient maps 64 residual masks to 16 classes but every class mixes missed depths, so it is a signed covariance carrier, not a scalar `q_t` quotient.  The truncated-octahedral graph is the Bruhat `S4` adjacent-swap carrier with 6 commutation squares, 8 braid hexagons, and two edge orbits (outer swaps vs middle swaps), matching the HYP-2889 local-compression failures.  New subtarget for the HYP-2890 residual leak: decompose `R_sf` over tangent Clebsch classes and Bruhat square/hexagon faces; prove square/commuting components nonpositive by design balance and isolate braid hexagons as the finite AP/Freiman low-depth atlas.
**OPEN-Q-108 S105 design/Hodge carrier split:** HYP-2892 refines HYP-2891 by turning the unital/Clebsch/truncated-octahedral prompts into a concrete carrier program for the HYP-2890 residual leak.  The `q=3` Hermitian unital is a `2-(28,4,1)` design, exactly matching the `C(8,2)=28` k=8 AP pair slots and giving `N^T N=8I+J`; Clebsch closed neighborhoods give a folded-parity `2-(16,6,2)` frame `N^T N=4I+2J`.  The truncated-octahedral graph is the `S4` adjacent-transposition/Bruhat graph with `6` square and `8` hex Coxeter faces, so failed one-step compression should be replaced by curl bounds on square/hex faces.  KPS S31m/S31n guardrail: this is not a sparse-design exact-tiler route; structured exact tilers are AP/dilates, and Goddyn-Wong `{1,...,11,13,24}` is a sporadic tight tiler.  New subtarget: attach `R_sf(E)-R_sf(AP)-Gamma_sf(A*(AP)-A*(E))` to Bruhat edges/faces for low-depth near-AP rows, then use unital/Clebsch block averages and HYP-2887/HYP-2636 curl/tail cancellation for the remaining labelled residual.
**OPEN-Q-108 S107 unital pair-slot realizability guardrail:** HYP-2894 tests the tempting literal map `q=3 unital points = C(8,2)` by enumerating all `S8`-invariant four-edge block systems on `K8`; none realize a `2-(28,4,1)` design.  AP8 sum/difference packets are visibly nonuniform, and THM-560/HYP-2893 now say the proof-critical tight locus is category-1: AP-dilates plus Goddyn-Wong acceleration/gap-collision atoms.  Therefore the remaining unital task is not a canonical `S8` design identification but a labelled or weighted map from AP/GW tight-locus packets into pair-frame coordinates for the HYP-2890 residual leak.
**OPEN-Q-108 S31m/S31n/S40 correction to S105:** Do not promote cut-side carriers into the final coverage invariant.  KPS S31m refutes the score-variance/Jensen coverage route and the sparse `PG(2,3)` design analogy; KPS S31n proves the structured diff-closed tilers are `Z_14 \ {0}` dilates but also verifies the sporadic Goddyn-Wong tight row.  mac-mini S40/S41 places Clebsch as folded-cube/cut-space, the truncated octahedron as the `S4` permutohedron/order side, and coverage as the finer observer-relative category.  Updated use of HYP-2891: a finite labelled residual atlas for covariance and compression faces that feeds tight-locus rigidity, not a standalone LRC invariant.

**OPEN-Q-108 S104 additive-energy tail refinement:** HYP-2890 turns KPS S31k's positive leading coefficient into a positive same-frequency packet with an explicit `1/m^4` tail: `Gamma_k(m)=C_{k,r}/m^4` on `r mod 7`, all constants positive for `k=8..13`, and tail after `H=12` `<=1.084e-6` at k=8/9.  This does **not** close by scalar additive energy: the full packet overpredicts AP by roughly 2x.  The sharpened analytic target is the residual-leak inequality `R_sf(E)-R_sf(AP)<=Gamma_sf(A*(AP)-A*(E))`, where `R_sf=p0-p0_decorr-Gamma_sf A*`.  Exact anchored scans show 0 violations for k=8 (`3432` rows, worst ratio `0.469`) and k=9 (`3003` rows, worst ratio `0.933` at `(0,2,4,6,7,8,10,12,14)`).  This integrates HYP-2889's AP-facing warning and HYP-+2888's scaling-invariant tiling signal: low-depth near-scale-reducible rows need a finite AP/Freiman residual-leak atlas, while high-depth labelled packets should route through HYP-2636/HYP-2887 cancellation.

**OPEN-Q-108 S31l/S104 signed-tail synthesis:** KPS S31l shows why the S104 residual cannot be bounded by termwise positivity: higher additive moment coefficients are mixed-sign (`k=9` has `s=3,s=4<0`; k=12 has negative `s>=4`).  Thus the positive same-frequency `s=2` packet is only the first convex carrier.  The closing theorem should be Jensen/Schur-convex over the AP-facing labelled difference profile, with the residual-leak inequality as the finite quantitative target.

**OPEN-Q-108 LEAN BOUNDARY (mac-mini-S27 retraction + codex S86g positive-p0 theorem): both NU and p0 routes are viable once the proof consumes only `0 < witnessG2`.** The witness route's event side is now fully concrete + sorry-free: `coverSet` (p0 event), `safeSet` (G_P event), `denseSet`, all measurable (`LRCDenseCovers`); `slowμ`=volume.restrict[0,1) is an `IsProbabilityMeasure`; the general Bonferroni `mu(A∩B)>=muA+muB-1` (`LRCBonferroniMeasure`); and the concrete floor `witness_floor_concrete : meas(G_P)-p0 <= meas(coverSetᶜ ∩ safeSet)` (`LRCWitnessFloorConcrete`; via Bonferroni(coverSetᶜ,safeSet)+complement identity). `LRCWitnessBonferroni` now has the corrected positive-p0 assembly `lrc14_from_p0_positive_wide_bound_split_nodes`: large clusters need only `0 < delta` with `p0<=cap-delta`, not the false k=8 comparison `witnessMP<=delta`; small clusters can still use the existing `m_P` floor. This is `Verify`-audited alongside `witness_margin_from_p0_wide_bound_shapes` and `large_witness_pos_from_p0_wide_bound_shapes`.

**OPEN-Q-108 KPS S31b concrete-p0 addendum:** `LRCFourteenSkeleton.p0` is now definitionally `DenseCovers.p0 E = (slowμ (coverSet E)).toReal`, not an opaque carrier.  `LRCP0Concrete` proves `0≤p0`, `p0≤1`, `p0=0` for fewer than six distinct speeds, monotonicity, and `wideBoundConcrete_of_decorrelation`, which packages hp0cap as the concrete cover-event inequality modulo the resonance/decorrelation residual.  This strengthens the p0 route's interface: the remaining `p0≤cap−delta` hypothesis now talks about the actual `coverSet` event consumed by the goodSet/Part-A bridges.

**OPEN-Q-108 S86g strict-cover addendum:** `LRCWitnessFloorConcrete.witness_pos_from_strict_cover_bound` now matches the exact output form of `LRCCoverBound.slowμ_coverSet_lt_cap`: `p0(E)<cap_k` and non-strict `cap_k<=meas(G_P)` imply `0<slowμ((coverSet E)^c ∩ safeSet P)`.  The p0 positivity route therefore does not need to move strictness onto `hmeasGP` or manufacture a separate delta at the concrete floor layer; the older margin lemma remains the finite-ruler error-budget version.

**OPEN-Q-108 S86g dense-complement addendum:** `LRCDenseCovers.coverSet_compl_subset_denseSet_compl` and `LRCWitnessFloorConcrete.dense_compl_witness_pos_from_strict_cover_bound` now transfer the strict hp0cap floor to `0<slowμ((denseSet E)^c ∩ safeSet P)` for anchored `0∈E`.  This is a formal max-gap proxy: it proves the p0 lower carrier lies in the complement of the 1/7-dense event.  The remaining readout is not another Bonferroni inequality; it is the cyclic sorted-gap bridge identifying `(denseSet E)^c` with the concrete `goodSet E` carrier used by `witnessG2`.

**OPEN-Q-108 S86g dense-complement margin addendum:** the same proxy now carries the quantitative p0 margin: `dense_compl_witness_margin_from_wide_bound` proves `delta≤slowμ((denseSet E)^c ∩ safeSet P)` from `p0(E)≤cap_k−delta`, `cap_k≤meas(G_P)`, and `0∈E`.  This is the interface needed by the finite-ruler error-budget side; the remaining open step is still the `denseSetᶜ` to `goodSet`/`witnessG2` readout.

**OPEN-Q-108 S86g phase-gap addendum:** `LRCDenseCovers.exists_phase_arc_empty_of_not_dense` now proves the finite cyclic-gap part of the readout: `¬Dense17` gives a phase with empty right arc `(0,1/7]` in `Int.fract(c-a)` coordinates.  The new `phaseGapSet` satisfies `(denseSet E)^c⊆phaseGapSet E`, and `LRCWitnessFloorConcrete` transfers both strict positivity and the p0 `delta` margin to `slowμ(phaseGapSet E ∩ safeSet P)`.  The remaining formal quotient to `GoodSet.goodSet E` is now specifically the speed-level identity from phase differences to `frac((b-a)x)` plus the finite-list witness packaging.

**OPEN-Q-108 S86g goodSet readout addendum:** the previous quotient is now closed.  `LRCGoodSet.phaseGapSet_subset_goodSet` and `denseSet_compl_subset_goodSet` prove the event readout into concrete `GOOD`; `LRCWitnessFloorConcrete.goodSet_witness_margin_from_wide_bound` gives `delta≤slowμ(goodSet E ∩ safeSet P)` from the p0 margin, and `goodSet_witness_pos_from_strict_cover_bound` gives positivity from strict hp0cap plus non-strict `hmeasGP`.  The witness-floor readout side is now concrete over `GOOD∩G_P`; the remaining hard nodes are the analytic p0/cap/Part-A inputs rather than this event quotient.

**OPEN-Q-108 S86g2 shape-readout addendum:** the concrete `GOOD∩G_P` carrier now plugs into the shape-indexed proof DAG.  `LRCEventMeasureBridge.shape_goodSet_witness_margin_from_wide_bound` proves `delta s ≤ witnessG2 s` from a readout equality `witnessG2 s = slowμ(goodSet(Eof s) ∩ safeSet(Pof s))`, anchored `0∈Eof s`, `p0(Eof s)≤cap s−delta s`, and `cap s≤slowμ(safeSet(Pof s))`; the strict-cover positivity analogue is also audited.  Thus the p0 route's remaining large-branch interface can be stated directly on `witnessG2`, leaving hp0cap/hmeasGP and finite-ruler Part A as the hard nodes.

**OPEN-Q-108 S86g2 Part-A goodSet-margin addendum:** `LRCWitnessPartA` now composes the concrete goodSet/safeSet readout with the finite-ruler Part-A budget.  `finite_witness_pos_from_goodSet_margin_shapes` and `finite_witness_pos_from_goodSet_margin_uniform_arc_bound_shapes` turn anchored `Eof`, `p0(Eof)≤cap−delta`, `cap≤slowμ(safeSet Pof)`, `witnessG2=slowμ(goodSet Eof∩safeSet Pof)`, and the finite `rhoK/arcCount/Vmax` error budget into positive finite witness density.  `lrc14_from_finite_partA_goodSet_margin_shapes` packages the same bridge into the conditional LRC14 assembly.  Remaining work is now the concrete node content: hp0cap, hmeasGP, concrete `arcCount/rhoK`, and the finite-ruler approximation inequality.

**OPEN-Q-108 S86g2 concrete-p0 Part-A addendum:** the same Part-A bridge now has wrappers whose margin hypothesis is the named concrete atom `DenseCovers.p0 (Eof s)≤cap s−delta s`, where `DenseCovers.p0 E=(slowμ(coverSet E)).toReal`.  `finite_witness_pos_from_goodSet_p0_margin_shapes`, its uniform-arc-bound variant, and `lrc14_from_finite_partA_goodSet_p0_margin_shapes` let the KPS S31b concrete p0 output feed the goodSet/safeSet finite-ruler route without unfolding p0 at each use.  This keeps the skeleton-facing hp0cap node aligned with the actual cover event.  Post-pull HYP-2840 reframes hp0cap through `p0≤L_y` and scalar `L_y` extremality; these wrappers are route-neutral and will consume either that Ly margin or the older decorrelation-style margin once formalized.

**OPEN-Q-108 current hard nodes:** the p0 route avoids the NU spreading lemma `hA` but still needs the concrete p0 margin `hp0cap`, the cap floor **hmeasGP** (`cap<=measGP`), and **hpartA** (slow-fast finite-ruler Part A). The NU route remains useful and stronger: `lrc14_from_bonferroni_split_nodes` uses Bonferroni + actual `nu` + cap floor and needs **hA** (`nu(E)>=nuConsec(k)`), **hmeasGP**, and **hpartA**. `LRCGapReach` closes the elementary geometric core of Part A (`>1/7` gap gives `>1/14` `nearInt` margin); the hard Part-A node is now the concrete rhoK/ruler approximation.

**OPEN-Q-108 STATUS: LEG C (GENUINE-WIDE) EXHAUSTIVELY VERIFIED (claude-opus-2026-06-22-S4, HYP-2825 corrects HYP-2817). ~1.59M genuine-wide doublet configs checked (CORRECTED from 1.16M — k=9 was not exhaustive in S3), 0 violations, ALL k=9..12 ALL gaps g=1..4, ALL bounded bases (exhaustive C(14,k-2) enumeration). Three-piece structure CLOSED: (I) frozen room Phi<cap + (II) Tornheim R-tail [M*_rig<=22] + (III) finite window [15,50] [exhaustive, all pass]. THM-527 rho*-CRUX VERIFIED FOR GENUINE-WIDE DOUBLETS (HYP-2826): rho*(P,E_co)>0 AND witness(P,E_co)>0 for ALL tested (k,P,B,g,M); global min rho*=2/147~0.0136>0; global min witness~0.483>0. LRC(14) reduces to: BOUNDED + SINGLE-FAR (THM-563, closed) + THIS LEG-C + L0 glue + Lean. See HYP-2817, HYP-2825, HYP-2826, reflection 07-reflections/lrc14-legC-closed-proof-synthesis-claudeopus-0622-S3.md.**

**OPEN-Q-108 — THE WIDE BOUND REDUCES TO CONCENTRATION EXTREMALITY OF L_yK8 (claude-opus-2026-06-22, HYP-2812). [LEAD]**
The cleanest closure of the whole wide region: **`max_E L_yK8` over ALL k-speed configs = `max_BOUNDED L_yK8` = MB < 10cap** (gK8=(10,0,0,1,0,0,10) on the miss-distribution, `L_yK8=10q0+q3+10q6`, `q0=p0`). VERIFIED EXACT over ~100k wide configs (incl. ALL binding families + small-M R-tail bumps): k=10 MB=5.286 vs MW=4.813; k=11 6.032 vs 5.632; k=12 6.641 vs MW=6.286=**E*** (the dichotomy-breaker — BELOW MB). **NO wide config exceeds the bounded max.** So gK8's BOUNDED cert is GLOBAL ⟹ `p0 <= cap` for all E with NO dichotomy, NO doublet, NO R-tail, NO frozen room — the wide leg collapses into the BOUNDED leg. MECHANISM/proof-route: gK8 charges the EXTREME miss-counts q0,q6; both are individually maximized by CONCENTRATED (bounded/slowest) configs (q6=all-missed maximized by slowest speeds; q0=coverage by the tight instance consec), and decorrelation smooths the miss-distribution to the MIDDLE (= the survival-middle-mass currency HYP-2701, now MONOTONE). REMAINING: prove concentration extremality (global extremality of consec under gK8 = a smoothing/majorization lemma on the 7-simplex). The explicit FALLBACK (if concentration extremality resists proof) is the generalized-doublet + Tornheim-R-tail route (HYP-2807/2808: max gw config is a generalized doublet {M,M+g}, R-tail = Mordell-Tornheim double sum |R_g|<~2.9 uniform). → HYP-2812, HYP-2811, HYP-2809, HYP-2807, HYP-2808, HYP-2701, THM-534, THM-538, mac-mini gK8/HYP-2810.

S79 correction: the small-`f` q6 contraction sublemma must use a boundary envelope, not a uniform `1/7` ratio.  Exhaustive single-far scans for k=10,11,12 through `15<=f<=60` give exact worst ratio `14/15`, while adjacent two-far frontier packets reach `7/8`; all such rows remain gK8-safe.  Add HYP-2822 to the proof DAG before applying asymptotic equidistribution.

S80 relation-depth refinement: the old one-peel dichotomy should be replaced by a two-peel relation-depth separator.  Exact `span<=18` genuine-wide scans have k=10 `over_Q=0`, k=11 `over_Q=0`, and k=12 `over_Q=4`, with every positive row at peel depth `2` and two-peel span `<=14`; the seven positive S78 `span<=20` witnesses also have depth histogram `d2:7`.  New subtarget: prove depth>=3 rows satisfy `p0<=Q(k-1)` or directly `p0<cap_k`, while depth 1 routes to single-far endpoint-period bounds and depth 2 routes to generalized-doublet/R-tail finite atlases.  Add HYP-2828 to the proof DAG as the covolume/Freiman separator target.  HYP-2828 is complementary to HYP-2823's exact degree-4 gK8 moment inequality `10-10S1+10S2-9S3+6S4<=10cap`, and should be treated as a finite resonant taxonomy for the sector-cap branch.  Do not confuse it with KPS HYP-2825/HYP-2826/HYP-2827, which attack the parallel 1/7 witness-floor route.


**OPEN-Q-108 — WIDE REGION UNIFIED BY gK8 + the GENERALIZED-DOUBLET / TORNHEIM-R-TAIL frame (claude-opus-2026-06-22, HYP-2807/2808/2809).**
Two convergent closures for the whole WIDE bound `p0(E)<cap_k`, k=8..12:
- **(CLEANEST) gK8 unifies all wide families.** The Delsarte dual `gK8=(10,0,0,1,0,0,10)` on the MISS-DISTRIBUTION `q_t=meas{exactly t of 6 sectors missed}` (`q0=p0`) gives `10*p0 <= L_yK8 = 10q0+q3+10q6` (trivial) with content `max_E L_yK8 <= 10*cap`. VERIFIED EXACT (HYP-2809) on ALL binding wide families with margin >= 0.138 (p0-units): genuine-wide maximizers, mac-mini's k=12 breaker `E*`, single-far near-cap plateau, AND dilated even-AP. **ONE moment cert bounds single-far + genuine-wide + dilated together -- superseding the binding/genuine-wide DICHOTOMY.** Remaining: `max_E L_yK8<=10cap` over ALL wide E (Delsarte LP feasibility for wide q-moments).
- **(EXPLICIT) generalized doublet + Tornheim R-tail.** The genuine-wide maximizer is a GENERALIZED doublet `{M,M+g}` (any base, any gap g -- HYP-2807); mac-mini's k=12 `E*` is the g=2 slice, NOT a new regime. THM-564's P/R split extends to gap g; the R-tail `R_g=M*(d2_g-d_inf)` is a convergent MORDELL-TORNHEIM double sum `<= (1/pi^3)*(#sector-pairs)*S`, `S~5.95`, giving `|R_g|<~2.9` ABSOLUTE & UNIFORM over (base,gap) (HYP-2808; empirical sup 2.24). => uniform `G~4.4`, cutoff `M*~28`.
- **SYNTHESIS (gK8 + R-tail):** applying the P/R framework to `L_yK8` (10x margin) absorbs the R-tail trivially; moment frozen room `Phi_Ly(B,g)<10cap` HOLDS (margin >=1.81), `M*_Ly~28`.
- **DEFINITIONAL FIX (reconciles kps HYP-2805 vs mac-mini-S7):** "genuine-wide" = IRREDUCIBLE (remove-any-one stays wide), NOT just `primitive(FULL E)+span>14`. kps's k=10 `{0,2,...,14,15,16}` (265/588, margin 0.1537) is BINDING -- remove 15 -> all-even -> `2*consec_9` -> bounded -- so THM-563's job. The true IRREDUCIBLE genuine-wide max at k=10 is 0.4423 (margin **0.162 >= 0.16**). gK8 makes the split moot. -> HYP-2807, HYP-2808, HYP-2809, THM-564, THM-563, THM-534, HYP-2805, mac-mini-S7.


**OPEN-Q-108 — CORRECTION TO THE GENUINE-WIDE MARGIN CLAIM (kind-pasteur-2026-06-21-kpswf9, HYP-2805): the consec doublet is NOT the genuine-wide maximizer; robust margin 0.16 FAILS at k=10 (true max 265/588, margin 0.1537), though `p0 < cap` holds everywhere.**
The THM-564 / HYP-2804 closure analyzes the CONSEC doublet `consec_{K-2}∪{M,M+1}` (k=10 cap−sup=0.16188 ≥ 0.16). But an EXHAUSTIVE genuine-wide search (all (k−2)-subset bounded bases × adjacent far pairs, filtering on `primitive(FULL E)` NOT `primitive(base)`) finds a HIGHER config: **k=10 `{0,2,4,6,8,10,12,14,15,16}` = 265/588 = 0.45068, margin 0.15372 < 0.16** (and k=9 `{0,4,6,8,10,12,14,15,16}`=321/980, margin 0.1667). The maximizing BASE is DILATED (gcd 2, = 2·consec_8) while the full set is primitive (15 odd) — so HYP-2804's base-primitivity sweep MISSED it. NET: (i) genuine-wide `p0 < cap_k` STILL HOLDS at every k (the actual LRC requirement; all margins > 0) — the leg closes; (ii) the ROBUST 0.16 reframe is UNATTAINABLE at k=10 (margin 0.1537). The frozen-law + correction machinery (HYP-2806) applies to these dilated-base doublets too (base=dilated even, far doublet {15,16}); the closure should target `< cap` at k=10 or fold the dilated bases into the doublet family. ACTION for any all-bounded-bases doublet sweep: filter on `primitive(FULL E)`. Scripts `04-computation/lrc14_genuine_wide_true_maximizer_kpswf9.py`. -> HYP-2805, HYP-2806, THM-564, HYP-2804, HYP-2795.
S77 addendum: `TournamentH7.LRCGenuineWideCorrection` now formalizes the exact arithmetic of this correction table (all rows below cap; k=10 smallest; `4/25` robust margin false; non-primitive base guardrail). This does not close HYP-2807's generalized-doublet finite window, whose current naive exact runner still needs optimization before it can be used as a certificate.

**OPEN-Q-108 — GENUINE-WIDE BINDING LEG (THE DOUBLET) VERIFIED-CLOSED via the almost-periodic P/R split (kind-pasteur-2026-06-21-Swf9, THM-564 / HYP-2804; CONVERGES with concurrent kps-S27 HYP-2799/2803 frozen-phase route).**
The genuine-wide leg's BINDING sub-case — the doublet maximizer `E_M=consec_{K-2}∪{M,M+1}`
(opus HYP-2797) — is now VERIFIED-closed with the LRC margin 0.16. The mechanism is the
**doublet analogue of THM-563**: center the deviation at the EXACT frozen plateau `Φ_frozen`
(not `bvd(base,2)`, which leaves a drift `d_inf` making `M·error→∞`, the HYP-2798 puzzle), then
`g(M)=M·(p0(E_M)−Φ)=P(M)+R(M)` exactly (inclusion-exclusion), with `P` exactly period-`7·lcm(base)`
(EXACT finite period-max, THM-563 sawtooth) and `R=M·(d2(M)−d_inf)=O(1/M)` (VERIFIED, Koksma on the
single residual far phase). Closure: `p0≤cap−0.16` for `M≥M*=⌈(period-max(P)+sup|R|)/H_K⌉` (TINY:
M*=13,24,55 for K=8,9,10) + an EXACT finite window `[15,M*]`. CLOSED K=8,9 (full period 420) and
K=10 (full period 2940, the BINDING case: `cap−sup_p0=0.16188`); K=11,12 window-verified (non-binding).
The crude `sup_g/15<H_K` FAILS at K=9 (razor-thin) — the careful `M*`-cutoff is what closes it.
REMAINING for a full PROVED status: the general bounded-base R-tail bound (a finite check, the
doublet analogue of THM-563's completed 12805-base periodmax check) + Lean formalization.
Scripts `04-computation/lrc14_doublet_{almostperiodic_PR,PR_closure}_kpswf9.py`,
output `05-knowledge/results/lrc14_doublet_PR_closure_kpswf9.out`. -> THM-564, THM-563, HYP-2804,
HYP-2799/2803 (kps-S27 convergent), HYP-2797, HYP-2798, HYP-2796, HYP-2795, HYP-2793.

**OPEN-Q-108 — LRC14 FRONTIER COMPRESSED AFTER BOUNDED AND SINGLE-FAR CLOSURE (codex-2026-06-21, HYP-2793).**
The current endgame should be handled as a proof DAG, not as one scalar search.
The bounded span<=14 leg is computationally closed for k=8..12 with split
completeness at reduced span 14; formalize the route/cap ledger.  The
single-perturbation / single-far leg is also computationally closed by the
complete `THM-563` bounded-base finite check:
`lrc_periodmax_THM563_general_check_COMPLETE_macmini_0621s7.out` checks all
`12805` primitive bounded bases `B subset [0,14]`, k=8..13, with `0` fails and
`0` skipped, proving `periodmax(B)<15*(cap_k-Plat(B))` everywhere.  The global
binding row is the k=9 even AP `(0,2,4,6,8,10,12,14)`, ratio `13.2805<15`.
The genuine-wide leg remains the live mathematical target and must use
relation-lattice / survival-middle-mass currency, because independent
`decorr_sup+err_sup` is false: room and resonance error anti-correlate.  Live
targets: (A) Lean/formal split + cap constants, (B) formal import/proof
compression for the completed periodmax certificate, (C) genuine-wide
pointwise room-vs-error or survival-currency signed-deviation theorem.  S77
now supplies the first periodmax formal kernel:
`LRCPeriodmaxCertificate.lean` proves the six per-k worst-row headrooms, the
k=9 worst-row comparison, the `12805`-base count checksum, and the k=8
`periodmax=2` normalization guardrail; full row enumeration remains the
Python/mac-mini audit. -> HYP-2793, THM-563, HYP-2788, HYP-2790, HYP-2792,
THM-557, THM-548, THM-546, HYP-2701, HYP-2684.

**OPEN-Q-108 — SINGLE-FAR CLOSED AS A FINITE PERIODIC MAX (mac-mini-2026-06-21-S6, THM-563).** The
signed-cancellation wall (HYP-2784: absolute bound 125× lossy) is COMBINATORIAL, not analytic.
`w·Δ_w = Σ_j Σ_{endpoints t of A_j} ±S_j(frac(w·t))` (exact Dedekind/sawtooth identity); the arcs `A_j`
depend ONLY on the base `B`, so `w·Δ_w` is EXACTLY PERIODIC in `w` with period `7·lcm(B)`, and
`sup_w Δ_w·w` = a finite exact period-max. For the binding consec bases: period-max = `1, 43/49, 1007/980`
(k=8,9,10), all `< 15·margin` ⟹ `Δ_w < margin` for ALL `w≥15` — no `w`-window finite check, no Koksma,
no reciprocity. The DILATED case (kps's single-perturbation reduction) closes via the CONTINUOUS period-max
(`contmax < 14·margin`: `1.0, 0.895, 1.028`). **COMBINED with kps HYP-2788** (near-cap ⟹ single-perturbation
⟹ single-far; genuine-wide ⟹ slack floor): the wide region closes via THM-563 + the slack floor, AVOIDING
the joint 2D ET-Koksma. Single-block domination (mine) + kps THM-557 confirm multi-block ≤ single-block.
Remaining: period-max ≤ 15·margin exhaustive over all bounded bases (dangerous k=8,9 verified, worst ratio
~10.8 < 15); kps's dichotomy rigor; HYP-2603 (consec maximizes Plat). → THM-563, THM-546, HYP-2788, HYP-2787.

**OPEN-Q-108 — PHASE-CARRIER / MAGIC-RANK FILTER (codex-2026-06-20-S67, HYP-2711/T942).**
The latest analogy batch does not add a proof by itself, but it clarifies which structures can be used safely.
Exact carriers are the mod-7 additive-character path integral for sector surjection, the death-chain Gibbs
transfer matrix, signed-incidence Hopfield energy, HYP-2707's Clifford/Gauss-sum tournament layer, and finite
even-page crossing atlases.  Literal Arnold cat maps, strict Cerny synchronization, beta-convexity,
low-order log-linear free energy, and raw path coherence are guardrails rather than proof routes.  The resulting
OPEN-Q-108 refinement is: keep the generated residual/phase profile until the final cap comparison, define its
LRC phase degree or magic rank, prove a Fubini-Study/projective-angle bound from the decorrelated death-chain
profile outside finite low-rank AP/cube-root/Freiman/squarefree atlases, and make the signed deviation smaller
than the HYP-2701 two-far boundary margin.  This is a filter on the HYP-2705 architecture, not an independent
claim that LRC14 is proved.
-> HYP-2711, HYP-2705, HYP-2707, HYP-2706, HYP-2704, HYP-2701, HYP-2702, HYP-2698, HYP-2684, HYP-2675, OPEN-Q-108.

**OPEN-Q-108 — TRUE-WIDE SURVIVAL MIDDLE-MASS GATE (codex-2026-06-20-S64, HYP-2701/T936).**
HYP-2695's true-wide cap-floor gate has an exact survival coordinate.  Since THM-556 gives
`U4=p0+p5+5p6`, the floor comparison `U4<=floor_k=(k-6)/7` is equivalent to
`p1+p2+p3+p4-4p6 >= (13-k)/7`; the cap comparison uses the same left side with right side
`1-cap_k`.  Exact S64 scan found no true-wide cap failures and no `k>=9` floor failures in the
audited boxes (k=8,9 B18; k=10,11,12 B17).  The only floor failures are three `k=8` rows and all
are cap-safe; worst `(0,3,6,9,12,14,15,18)` has floor debt `107/8820`, cap slack `295/3528`,
`p6=1/126`, and `far_count=2`.  The far-count ledger is the new proof split: every tight audited
leader is barely true-wide (`far_count=2`), while `far_count>=3` has larger margin in each layer.
Two-far addendum: the fully decorrelated death-chain boundary currency is already floor-safe for
every bounded core scanned, with margins k=9 `569/3430`, k=10 `5717/36015`, k=11 `5317/24010`,
k=12 `35543/123480`.  Actual two-far rows spend this by a negative signed deviation; the k=9
leader has boundary margin `18119/72030`, deviation `-6395/28812`, and remaining slack `29/980`.
Sharp target: prove `boundary margin >= -signed two-far deviation` for `k>=9`, using
Freiman/scale finite atlases for low relation-distance far pairs and signed Abel/BV bounds off
resonance; then prove a separate easier `r>=3` far-count margin lemma and reserve finite THM-535
dividend templates for true-wide `k=8`.  Near-equality rows must be lifted back to state-word,
transfer-tax, residual-profile, and relation-lattice coordinates before scalarizing.
-> HYP-2701, HYP-2695, HYP-2699, HYP-2700, HYP-2702, THM-548, THM-556, THM-535, HYP-2680,
HYP-2679, HYP-2696, HYP-2698, HYP-2684, HYP-2675, OPEN-Q-108.

**OPEN-Q-108 — RESIDUAL-PROFILE AUTOMATON CONE (codex-2026-06-20-S63, HYP-2698/T934).**
HYP-2697's generated residual-profile cone now has an exact coordinate.  Actual decorrelated contexts are
words `x -> w_x(R)` over the 64 residual masks, updated at fixed slow coordinate `x` by OR-convolution /
residual deletion.  In miss-zeta coordinates `Z_x(A)=Pr(A⊆residual)`, independent context clusters satisfy
the pointwise product law `Z_context,x(A)=prod_i Z_i,x(A)` before averaging over `x`; this is the cone
structure that arbitrary nonnegative residual weights destroy.  Exact S63 scout over coherent contexts from
integer partitions found that all S62 coordinatewise counterexamples lose against every generated context
with total nonzero size `m=7..11`, with worst margins `20/16807`, `37/16807`, and `199/24010`; a
near-consecutive frontier scan through size 6 found zero failures, worst `12027/2352980`.  New sharp target:
prove compression for miss-zeta product words.  The all-singleton context base case now reduces exactly to
the hit-count kernel `g_r(t)=7^-r sum_j (-1)^j binom(t,j)(7-j)^r`; prove that hit-count majorization first,
then prove coherent context merging does not reduce the margin, and keep sector labels available for THM-558
transfer-tax near-equality routing.
-> HYP-2698, HYP-2697, HYP-2696, THM-558, THM-557, HYP-2694, HYP-2684, HYP-2675, OPEN-Q-108.

**OPEN-Q-108 — ARBITRARY CLUSTER COMPRESSION CONE (codex-2026-06-20-S62, HYP-2697/T933).**
The arbitrary-shape part of HYP-2694 cannot be solved by coordinatewise stochastic dominance.  For a cluster
shape `C`, let `q_C(R)=Pr_{x,phi}(C covers residual sectors R)`.  Exact counterexamples, after S63's
pairwise-difference wall refinement: `(0,1,3)` beats `(0,1,2)` by `5/294` on several 3-sector residuals;
`(0,1,2,4)` beats `(0,1,2,3)` by `3/196` on several 4-sector residuals; `(0,1,2,3,5)` beats
`(0,1,2,3,4)` by `37/2940`, while full-cover difference is zero.
Thus arbitrary positive residual weights are too strong.  The new target is a generated-cone theorem:
characterize residual profiles `w_R` arising from actual decorrelated LRC contexts and prove
`Σ_R w_R q_C(R) ≤ Σ_R w_R q_K(R)` on that cone, where `K` is the coherent consecutive block.  Grid scouts
still put consecutive first for full-cover scalar in bounded sizes `6..9`, and split-context beams stay below
the one-block branch; exact singleton+size10 checks also favor the consecutive large block.  After this cone
compression, THM-557 split gaps and HYP-2684 joint carrier error become the remaining spendable margin ledger.
Incoming HYP-2696/THM-558 supplies the complementary transfer-tax account for unpaid one-missed closures after
sector-state insertion.
-> HYP-2697, HYP-2696, HYP-2694, THM-558, THM-557, HYP-2684, HYP-2675, OPEN-Q-108.

**OPEN-Q-108 — SINGLE-BLOCK EXTREMALITY: JOINT CARRIER GAP REMAINS (codex-2026-06-20-S61, HYP-2694/T931, THM-557).**
The HYP-2694 single-block wide-cover route now has a closed coherent-block core.  THM-557 proves, by exact
Fraction integration, that if anchor `0` is fixed and the `m=k-1` nonzero runners are partitioned into far
coherent consecutive blocks, then the one-part block `[m]` maximizes the shared-`x` decorrelated cover for
`m=7..11`.  Exact values/margins are `D_7=283/1470` with cap margin `1111/5880`, `D_8=629/2058` with
`111019/588588`, `D_9=16969/41160` with `102803/535080`, `D_10=30551/61740` with `184957/802620`, and
`D_11=71111/123480` with `34729/123480`.  Closest split is always `[m-1,1]`, with explicit split gaps
`1111/10290`, `374/5145`, `6561/96040`, `42661/864360`, `9047/172872`.  Single shifted blocks also have the
proved diagonal-freeze error `|p0({0}∪{M..M+m-1})-D_m|≤7*C(m,2)/M`, giving conservative large-M cutoffs
`779/1040/1312/1367/1369`; exact `M=19` rows are already below cap.  The remaining open target is now sharply:
prove arbitrary bounded cluster shapes compress to the coherent-block quotient, and prove the joint multi-carrier
decorrelation error is bounded by the available `cap margin + split gap`; then finite-check the small carrier gaps.
-> THM-557, HYP-2694, HYP-2684, HYP-2675, HYP-2695, OPEN-Q-108.

**OPEN-Q-108 — TRUE-WIDE CAP-FLOOR GATE (codex-2026-06-20-S60, HYP-2695/T930).**
HYP-2693's true-wide Bonferroni4 cap gate now has a sharper currency split.  THM-535 proves
`cap_k>=floor_k=(k-6)/7`; exact S60 audit shows the true-wide rows with `k>=9` appear to satisfy
the stronger floor gate `U4=p0+p5+5p6<=floor_k`, so exact cap minimizers should only be needed for
the `k=8` dividend row.  Exact boxes: k=8 B18 has `3` true-wide floor failures but `0` cap failures;
k=9 B18 (`27020` true-wide), k=10 B16 (`3432`), k=11 B16 (`3003`), and k=12 B16 (`2002`) have
`0` true-wide floor/cap failures.  Tight slacks: k=9 `U4=391/980`, floor slack `29/980`; k=10
`U4=307/588`, floor slack `29/588`; worst k=8 floor debt `107/8820` is covered by cap slack
`295/3528` and dividend `563/5880`.  New subtarget: prove the state-mass/decorrelation high-tail
floor lemma for true-wide `k>=9`, and route true-wide `k=8` through finite cap-dividend templates;
boundary/AP stays on HYP-2691.  Incoming HYP-2694's single-block decorrelated wide-cover extremizer is
the complementary partition-function route; HYP-2695 is the final-row Bonferroni currency split. ->
HYP-2695, HYP-2694, HYP-2693, THM-535, THM-556, HYP-2683, HYP-2684, HYP-2676, HYP-2677,
OPEN-Q-108.

**OPEN-Q-108 — INCLUSION-EXCLUSION-OVER-N COMPREHENSIVE VIEW + REDIRECT (mac-mini-2026-06-20-S5, HYP-2692).**
The LRC's inclusion-exclusions are one arithmetic skeleton indexed by `6=2·3`: N=7 sectors (=THM-534
moment-LP, optimal closes k=8-10, plain Bonferroni fails); N=2 quadratic χ (QR/NQR, Gauss sum √−7,
reality, Chebyshev bias ~70% non-universal); N=3 cube root C_3 (Eisenstein); keystone — C_3 orbit-sum
of 7th roots = Gaussian period `(−1+χ√−7)/2`, the correction's C_3 trace ∈ Q(√−7); N=runners danger
sieve DEAD (L=0 at tight {1..13}). **VERDICT:** the multiplicative arithmetic WASHES OUT on `p0`
(characters vanish, bias is archimedean); incl-excl over arithmetic N does not bound `p0`. **REDIRECT
(verified):** the lever is the summed **leading-order** residual `R_{s_0}`, `s_0=max(1,7−|B|)` — for the
dangerous true-wide leader (|B|=7) that is `R_1` (one-far, barely-far = THM-546/547 collar machinery,
the best-understood piece, not the d≥2 lattice); for sparse cores it is `R_3+`. `R_tot=p0−P_r` stays
within margin in all tested rows. The Q(√−7) C_3-orbits INDEX the resonance classes the height-weighting
sums over; smallness stays signed archimedean. → THM-534, THM-548, THM-551, HYP-2657, HYP-2684, HYP-2692.

**OPEN-Q-109 — The base-HP / grid-symmetric maximizer lemma (the one gap in HYP-2688).**
HYP-2688 (VERIFIED exhaustive n=3..9): the global H-maximizer is attained inside the
`2^{half(n)}` grid-symmetric (phi-self-converse) subcube of the tiling cube, giving a
`2^{(m-d)/2}` search reduction. The WEAK form ("the maximizer set CONTAINS a grid-sym point")
is what holds and what is needed; the STRONG form is REFUTED (non-grid-sym maximizers exist,
outnumbering grid-sym ones at n=5,6,7). `H=H^op` makes the maximizer set `rho`-invariant but a
size-2 `rho`-orbit has no fixed point, so invariance alone is insufficient. **The lemma to
prove:** *every regular self-converse global H-maximizer admits a base-Hamiltonian-path choice
`P0=n->...->1` under which its tiling is `rho`-fixed (grid-symmetric)*. The canon SC/regular
maximizer theorem gives the abstract self-converse symmetry (a `phi`-self-converse relabeling
exists); the gap is making that relabeling compatible with the FIXED base path. Proving it
upgrades HYP-2688 to a theorem and certifies the search speedup. → HYP-2688, THM-280, THM-552,
SC-maximizer canon.

**OPEN-Q-108 — TRUE-WIDE via DICHOTOMY-FINITE-REDUCTION (kps-2026-06-20-Sx-wf, INDEPENDENT route).**
A second, cluster-decorrelation attack on Region III (complementary to the THM-548 boundary-value picture).
[PROVED] `p0(λE)=p0(E)` for `gcd(λ,7)=1` (THM-531-B, re-verified for p0) + pigeonhole (a cluster of size
`≤6` has `p0=0` alone) ⟹ a true-wide set's `≥2`-cluster SHAPE-multiset family is FINITE. [VERIFIED/EXACT]
the correct SHARED-`x` decorrelation engine (RIGOR FIX: independent-`x` convolution under-counts, `17/343` vs
true `23/196`; only shared-`x` matches the `M→∞` limit) gives a worst decorrelated `p0_inf` over the WHOLE
finite family `= [k-1 consec]+[singleton] =` the THM-547 plateau `Qb(k-1)` `= 0.19660/0.36210/0.44789/
0.53125/0.60224 < cap_k` (margin `≥0.132`). So TRUE-WIDE decorrelates to NO WORSE than the (closed) boundary
collar. [REDUCTION] explicit gap cutoff `B=(6/49)·f·Vmax/margin = 682/1453/1774/1988/2034` (k=8..12,
conservative; signed Abel shrinks `5-76×`): gaps `>B` CLOSED, gaps in `(14,B]` a finite check gluing to span-14.
[VERIFIED] 0 violations of `p0≤cap_k` on >750 exact true-wide rows (margin `≥0.184`); `max(p0−p0_inf)~0.02-0.05`.
[SOLE GAP] the multi-cluster ERROR AGGREGATION `p0(E)≤p0_inf+Σ_far (6/49)V_i/g_i` (iterate of the PROVED
single-element THM-546) is verified-numerically but not yet closed-form — a pure PRODUCT/decorrelation bound,
not signed cancellation. Closing it finishes LRC(14). Files: `04-computation/lrc14_h2675_dichotomy-finite-
reduction_kps-Sx-wf.py`, `05-knowledge/results/...out`, reflection `07-reflections/lrc14-true-wide-decorrelates-
to-the-boundary-collar-plateau.md`. → HYP-2675, THM-547, THM-546, THM-531, OPEN-Q-108.

**OPEN-Q-108 — REGION III (true-wide) BOUNDARY-VALUE ARCHITECTURE (mac-mini-2026-06-20-S3, THM-548).**
The far-element process is a BOUNDARY-VALUE problem (the user's lead). Dictionary: bounded core `B`
= boundary point; far runner `w→∞` = radial approach; plateau `Φ(B)` = boundary function (Fatou,
the one-far limit PROVED to exist with rate `6/49`); two-far curvature `I_B(u,v)` = Bagemihl ambiguity
defect; resonance `mu+nv≈0` = ambiguous point = Freiman small relation. **Exact finite decomposition:**
`p0(B∪F) = Σ_{S⊆F,|S|≤6} Δ_S(B) = P_r(B) + R(B,F)`, where `P_r(B)=Σ_t prof_t(B)c_t(r)` is the
fully-decorrelated **Fatou boundary value** and `R` the resonance corrections. **VERIFIED this session:**
(a) `Φ₂(B)=(2p₂−p₁)/49`; (b) `P_r(B)≤cap_k` with margin GROWING in `r` (0.13→0.48 — boundary value is
safe); (c) two-far constant `C₂=13/1372=13/(2²·7³)` (parabolic analogue of one-far `3/49=3/7²`) — each
Abel order adds one power of the apex prime 7; (d) QR-reality of the PRODUCT kernel (licenses the signed
two-far bound); (e) resonance-gated `|I_B−Φ₂|≤C·V/resdist`, worst case `~0.013 ≪` margin — two-far
curvature is NEVER cap-threatening; (f) consecutive-pair curvature SATURATES (bounded). **CORRECTION**
(re-verified): the k=9 leader `(0,4,6,8,10,12,14,15,16)` has NEGATIVE curvature `−13/1470` and dilated
core `2·(0,2,3,4,5,6,7)` — it is the SCALE-INVARIANT branch, not a positive-synergy exception (HYP-2679
literal premise refuted). **REMAINING (honest):** the ONE genuine analytic gap is the SIGNED magnitude
bound for the `d≥2` relation lattice (absolute bound proven 5× lossy; needs Poisson/theta keeping
`(−1)^|T|` + the 7-vanishing); plus the divergent-resonance/stacking risk (sup `w|Δ_w|` grows with #scales,
closure relies on the offsetting tiny plateau — computational not yet analytic); plus the unrun finite
checks; plus the upstream glue (HYP-2603). LRC(14) NOT proved. → THM-548, THM-546, THM-547, HYP-2679,
HYP-2678, HYP-2637, HYP-2606.

**OPEN-Q-108 PROGRESS (mac-mini-2026-06-20-S2): 2 of 3 sector-crux regions now CLOSED.** The crux
`p0(E)≤cap_k` (k=8,9,10) splits by spread into three regions: **(I) finite half** `max(E)≤14` — PROVED
(kps-S19, 0 violations); **(II) boundary collar** `2nd-largest≤14, max=w>14` — CLOSED (THM-547) via
THM-546 sharpened `|Δ_w|≤(6/49)V(E')/w` for `w>w*=54/90/103` + a feasible finite check `14<w≤w*`
(k=8 verified: 19100 configs, 0 violations, worst margin 0.155); **(III) true-wide** `2nd-largest>14` —
OPEN, the Ruzsa/Plünnecke+Freiman program (HYP-2678: `d=1`→scale-invariance, `d≥2`→signed dimension
penalty). The signed (6/49) bound (THM-546 S2, from the QR reality HYP-2657, `6=−1`∈NQR mod 7) is the
shared engine codex's HYP-2676/2677 packet route also adopted. Only region III remains. → THM-546, THM-547,
HYP-2678, HYP-2676, HYP-2675.

**SIGNED MULTI-FAR BOUNDARY ADDENDUM (codex-2026-06-20-S51, HYP-2680/T919):** THM-548's two-far limit is the `s=2` member of an exact Newton/Stirling hierarchy.  If `p_t(B)` is the measure where bounded core `B` misses exactly `t` inner sectors, then the fully decorrelated `s`-far mixed term is `Phi_s(B)=7^-s sum_{t=1}^s (-1)^(s-t)t!S(s,t)p_t(B)`.  In particular `Phi_1=p1/7`, `Phi_2=(2p2-p1)/49`, and the signed three-far target is `Phi_3=(p1-6p2+6p3)/343`.  S51 corrects the truncation language: the Newton identity is over all far subsets, while the six-sector structure only limits the profile variables `p_0..p_6`.  Exact scout verifies the `P_r` boundary identity through `r=6` and shows `Delta_3-Phi_3` is governed by low-height triple forms.  For `(15,16,17)`, all `3003` bounded cores have exact relation `-u+2v-w=0` and no exact pair relation; deviations split `1999` positive, `1004` negative, with top abs deviation `40633081/445721640`, but direct cap margins remain large.  Incoming THM-548 simultaneous peel makes the `r=3` target precise: route `P_3`, one-far residuals, and two-far residuals by existing bounds, then prove the signed three-far relation-lattice bound for `Delta_{uvw}-Phi_3`.  Tail-rank addendum extends this to `r>3`: bound the signed order sums `R_s=sum_{|S|=s}(Delta_S-Phi_s)`.  The four-far bank `(15,16,17,18)` has `R2/R3` opposite signs in `1644/3003` cores and `R3/R4` opposite in `2053/3003`, so cancellation across Newton orders is a proof resource. -> THM-548, HYP-2680, HYP-2679, HYP-2678, HYP-2639, OPEN-Q-108.

**CUBE-ROOT ORDER-FILTER ADDENDUM (codex-2026-06-20-S52, HYP-2681/T920):** For a far triple, the seven packets `A,B,C;D,E,F;G` have actual correction `A+B+C+D+E+F+G`; the recursion `A+B+C-D-E-F+G` is the pair-tax shadow `H(1)-2(D+E+F)`.  Exact Eisenstein modes `S_omega=A+omega B+omega^2 C` and `P_omega=D+omega E+omega^2 F` preserve cyclic phase.  In the `(15,16,17)` all-core bank, actual residual signs are `+2821/-182`, while pair-tax shadow signs are `+1250/-1753`; actual residual, pair-tax shadow, Eisenstein imbalance, and direct `p0` choose different leaders.  Use the cube-root packet as a finite-atlas phase coordinate, not as a direct cap-risk scalar. -> HYP-2681, HYP-2680, THM-548, HYP-2677, HYP-2639, OPEN-Q-108.

**AP-TRIPLE PHASE-ATLAS ADDENDUM (codex-2026-06-20-S53, HYP-2682/T921; integrated with incoming KPS S19):** KPS S19 refutes the scalar `C(k)=sup w|Delta_w|` route and sharpens the surviving OPEN-Q-108 target to HYP-2675's Weyl/decorrelation route with exact plateau target `sup p0_decorr=Q(k-1)<cap_k`.  S53 tests the rank-one AP-triple branch `(m,m+1,m+2)`: the exact relation `u-2v+w=0` is fixed, but packet phase varies with core/offset and is not determined by `m mod 7`.  Named-core scan through `m=120` gives transitive core-family tournament `consec8 > direct_p0_leader_core > dilated > s51_top_dev_core > third_pocket_mixed_core`.  Selected all-core AP banks `m=15,16,22,28,42` all keep large direct-p0 margins; tightest direct row is `(0,9,10,11,12,13,14,15,16,17)`, `p0=2290763/5717712`, margin `1164997/5717712`.  Use AP packets as finite resonant phase/support labels inside the decorrelation/glue proof, not as a scalar rank/discrepancy bound. -> HYP-2682, HYP-2681, HYP-2680, HYP-2675, OPEN-Q-108.

**DECORRELATED PLATEAU-BOUND AUDIT (codex-2026-06-20-S53, HYP-2675):** Independent exact audit of the KPS S19 decorrelation foundation in THM-548/S51 language.  For bounded primitive bases `B subset {0,...,14}`, total `k=8..12`, base size `b`, and independent far count `r=k-b`, scan `P_r(B)=sum_t prof_t(B)c_t(r)`.  In every case the maximum occurs at `b=k-1`, `r=1`, `B={0,...,k-2}`, i.e. `Q(k-1)`: `Q(7)=289/1470`, `Q(8)=621/1715`, `Q(9)=1229/2744`, `Q(10)=65599/123480`, `Q(11)=14873/24696`, all below `cap_k`.  Remaining gap is explicit Weyl/decorrelation error plus finite bounded-gap glue, not the finite decorrelated comparison. -> HYP-2675, THM-548, HYP-2680, OPEN-Q-108.

**BV-FOURIER DECORRELATION ADDENDUM (codex-2026-06-20-S55, HYP-2684/T923):** Web/repo search on Weyl/decorrelation identifies a concrete analytic target for the remaining HYP-2675 error.  For a two-scale cluster coverage function `H(x,phi)`, the actual row and independent-anchor model differ by the exact resonance sum `int H(x,Mx)dx-int H(x,phi)dxdphi=sum_{s!=0}Hhat(-Ms,s)`.  If the LRC coverage function has mixed BV Fourier decay `|Hhat(r,s)|<=V_mix/(4*pi^2|rs|)`, the nonresonant error is `<=V_mix/(12M)`.  Thus the proof route is now: define the exact cluster coverage functions, prove an explicit mixed-variation budget, choose `G` with `V_mix/(12G)<cap_k-Q(k-1)`, and route low-height resonances to HYP-2682/HYP-2676 finite atlases.  Vertices for the tournament/quotient should be scale clusters or resonance equations, not runners. -> HYP-2684, HYP-2675, THM-546, THM-548, HYP-2682, HYP-2676, OPEN-Q-108.

**CUBE-ROOT PHASE/SUPPORT ADDENDUM (codex-2026-06-20-S54, HYP-2682/T921):** Holding relation rank fixed is still too coarse.  Exact scout scans consecutive rank-one far triples `(m,m+1,m+2)`, `m=15..35`, over all `3003` primitive bounded cores per triple (`63063` rows).  Every triple has `-u+2v-w=0`, but actual signs, pair-tax signs, cube-root A2 chambers, and direct-risk leaders vary with mod-7 phase/support.  Top-12 overlap: direct risk shares `6/12` rows with actual residual leaders, only `2/12` with pair-tax leaders, and only `2/12` with Eisenstein-norm leaders; pair-tax and Eisenstein share `5/12`.  Therefore the low-height resonant branch should first route finite keys `(far residues mod 7, actual/pair-tax/pair/triple signs, A2 chamber of S_omega-P_omega, bounded-core support/state data)`, and only then apply signed Abel/Koksma estimates.  KPS S19 now makes the global wide route decorrelation/coverage to the plateau `Q(k-1)`; HYP-2682 is the finite resonant router that keeps cube-root phase visible inside that route, not a replacement scalar constant. -> HYP-2682, HYP-2681, HYP-2680, HYP-2675, THM-548, HYP-2639, OPEN-Q-108.

**WIDE-ADDRESS REPAIR ADDENDUM (codex-2026-06-20-S55, HYP-2683/T922):** Exact B20 true-wide scout tests the old owner-private/compatibility-wall hidden gem in HYP-2675's direct-`p0` branch.  The scan keeps the true-wide leader `(0,4,6,8,10,12,14,15,16)`, `p0=321/980`, margin `11681/70070`, and audits `513` top/sample rows.  Scalar and additive summaries do mix proof states: scalar has `3` high/low mixed buckets, additive has `1`, and private mass alone has `3`.  The non-overfit repair is coarse missed-state mass: `state_mass=(support bucket, entropy bucket, binned p1,p2,p3)` has `286/513` buckets, `0` high/low mixed buckets, max width `52229/1113840`; residue-private is sharper (`480/513`, `0` mixed, max width `15027/340340`) but close to overfitting, while fine address is a row hash.  Pressure direction: high-risk rows average larger private mass (`~0.459` vs `~0.415`) and lower state entropy (`4.445` vs `4.731`) than low-risk rows.  New subtarget: prove a state-mass deficit lemma tying missed-state entropy/support to Freiman dimension or low-height compatibility packets; this finite resonant/address ledger then feeds HYP-2684's BV-Fourier nonresonant error before the final comparison to the `Q(k-1)` plateau. -> HYP-2683, HYP-2684, HYP-2675, HYP-2682, HYP-2648, HYP-2530, OPEN-Q-108.

**TRUE-WIDE BOUNDARY-CURVATURE ADDENDUM (codex-2026-06-20-S50, HYP-2679/T918):** The two-far boundary-function reframe has been tested exactly.  For a core `B` and far speeds `u<v`, define `I_B(u,v)=p0(B∪{u,v})-p0(B∪{u})-p0(B∪{v})+p0(B)`.  The k=9 scan over `B=(0)+6-subsets of [1,14]` and far pairs `15..24` checked `135065` primitive rows.  The direct-risk leader remains the HYP-2675 true-wide row `(0,4,6,8,10,12,14,15,16)`, with `p0=321/980` and margin `11681/70070`, but its two-far curvature is negative `-13/1470`; its core is `2*(0,2,3,4,5,6,7)`.  Largest positive curvature occurs at `(0,1,4,8,10,12,14,16,20)`, `I=307/1960`, but that row is safer (`p0=89/336`).  Reading: the live k=9 leader is a `d=1` dilated-core overlap row, not a high-dimensional positive two-far synergy.  Incoming THM-548 supplies the analytic companion target: `I_B(u,v)` should decorrelate to `Phi_2(B)=(2p2(B)-p1(B))/49`, with deviations governed by resonance frequencies `m*u+n*v`; nonresonant pairs should decay by signed BV/Abel bounds, and resonant pairs should Freiman/scale-reduce to finite atlas rows. -> THM-548, HYP-2679, HYP-2678, HYP-2675, THM-547, THM-531, OPEN-Q-108.

**PACKET-TOURNAMENT ADDENDUM (codex-2026-06-20-S49, HYP-2677/T917; integrated with incoming THM-546 S2 and THM-547):** HYP-2676's six packet signs have now been tested as K4 edge orientations and as six missed-sector vertices.  Exact atlas shows the top twelve B14 near-speed leaders all have `++++++`, identical transitive K4 type `scores=(3,2,1,0), c3=0, HP=1`, and identical cyclic-pair sector type `scores=(3,3,3,2,2,2), c3=8, scc=(6,), HP=41`; K4 signs locate the finite same-sign pocket but do not separate it.  The KPS third pocket has `++-+--`, `Delta=1171/452760`, Abel pressure `1171/133056` under the incoming rational bound `|Delta_w| <= (6/49)V/w`, a negative opposite-pair balance `3+4=-23/6468`, and cyclic-pair topology `scores=(4,3,3,3,2,0), c3=4, scc=(5,1), HP=15`.  THM-547 demotes the boundary collar to a closed/modulo-finite-check calibration branch, so the live OPEN-Q-108 subtarget is now sharper: classify finite `++++++` Ruzsa/Freiman packet models with opposite-pair balances and true-wide rows; prove that complements have a cyclic-pair arc flip, small exact pair mass, or enough THM-546 S2 signed Abel cancellation, then reattach HYP-2648 state-word and HYP-2639 relation-shell addresses. -> HYP-2677, HYP-2676, HYP-2674, THM-546, THM-547, HYP-2639, HYP-2648, OPEN-Q-108.

**SIGNED-PACKET/RUZSA ADDENDUM (codex-2026-06-20-S48, HYP-2676/T916):** HYP-2674's same-sign packet pocket has been refined into a finite-model versus signed-cancellation split.  Exact scout keeps the one-missed-sector packet telescope in Fractions and adds sumset excess, `K2`, `K3`, additive energy, squarefree profiles, THM-546 BV-budget share, and run-level cancellation.  The large positive near-speed rows in B14 are all `++++++` and low/small excess; the top is `(0,2,4,6,7,8,9,10), w=12`, `Delta=5347/30870`, excess `3`, `K2=9/4`.  Named finite models remain same-sign: B13 `997/5880`, HYP-2671 dyadic block `457/3920`, HYP-2672 doubled-odd `483281/5761028`.  The contrasting third-pocket row has sign `++-+--` and run cancellation `1171/15473`, suggesting the high-growth branch should be a signed packet estimate rather than an absolute one.  New OPEN-Q-108 subtarget: classify finite `++++++` packet models through Ruzsa/Freiman normalization, then prove signed Erdos-Turan packet cancellation or small mass on the high-growth complement, feeding the result back into HYP-2675's direct-`p0` wide/collar split. -> HYP-2676, HYP-2675, HYP-2674, THM-546, HYP-2638, HYP-2639, OPEN-Q-108.

**WIDE-RIDGE ADDENDUM (codex-2026-06-20-S47, HYP-2675/T915; integrated with THM-546):** exact direct-`p0`
scout confirms KPS's warning that `span>14` must be split before proving the
comfortable branch.  In the k=9 B20 scan (`125970` raw rows, `122922` primitive
`span>14` rows), the all-span leader is boundary
`E=(0,2,4,6,8,10,12,14,15)`, `p0=437/1176`, margin `20627/168168`, with
second-largest `14`; the true-wide (`second>14`) leader is
`E=(0,4,6,8,10,12,14,15,16)`, `p0=321/980`, margin `11681/70070`.  HYP-2671's
dyadic row has direct `p0=29/112`, margin `3769/16016`, so it is a
decoupled-Delta danger but direct-p0 comfortable.  New OPEN-Q-108 subtarget:
prove a boundary collar compression lemma for `second<=14`, then a true-wide
Freiman/GAP/state-word sector-cover deficit lemma for `second>14`, before using
the KPS post-25 packet tail.  Incoming THM-546 supplies the rigorous gapped
one-far decorrelation bound `|Delta_w|<=kappa V(E')/(pi^2 w)`; HYP-2675
identifies the ungapped finite ledgers where scale distance is absent. ->
THM-546, HYP-2675, HYP-2674, HYP-2653d, OPEN-Q-108.

**FINITE-HALF / PER-SECTOR ADDENDUM (codex-2026-06-20-S46 integrating incoming KPS 62fc2a58d):** KPS landed a stronger split after HYP-2674.  The finite sector half is now computationally certified for `max(E)<=14`, `k=8..12`: zero violations of `p0(E)<=cap_k`, with margins `cap_k-Q(k-1)` equal to about `0.185, 0.132, 0.157, 0.194, 0.255`.  The per-sector script proves/verifies the exact telescope over exact singleton-missed runs and the rigorous bound `w|Delta_w| <= (6/49)sum_s|R_s| <= (6/7)sigma(E')`; it also refutes any standalone bounded `w|Delta_w|` constant via `E'_M={0,1,2,3} union {M..M+3}, w=22M`, where `w|Delta_w|~0.08M`.  New synthesis: the remaining analytic input should be joint, not scalar: bounded bases close by `sigma(E')/w` decay after finite checking, and wide bases should have small plateau/p0 directly.  HYP-2674's `++++++` packet-alignment pocket is the bounded-near-plateau obstruction to classify inside that joint route. -> HYP-2674, HYP-2673, HYP-2653d, HYP-2671, OPEN-Q-108.

**PACKET-ALIGNMENT ADDENDUM (codex-2026-06-20-S46, HYP-2674/T913):** the corrected uniform `Delta_w` tail can be attacked through one-missed-sector packet sign words.  Exact scout decomposes `Delta_w=(1/w)sum_c[G0(w*b_c-s_c/7)-G0(w*a_c-s_c/7)]` by missed sector `s=1..6`.  The known risky rows are all same-sign packet alignments (`++++++`): the finite B13 pocket has `Delta_w=997/5880`, the HYP-2671 dyadic block has `Delta_w=457/3920` and margin gap `12223/784784`, and even the non-shell warning row has `++++++` but small absolute `Delta_w=11/315`.  In the dyadic family `E_s={0,1,2,4,8,3s,4s,5s}`, `s=3..120`, the `w=6s` spike is still finite at `s=4`; after `s>20`, the best sampled `Delta_w` is only `2539/64680`, leaving about `0.0929` k=9 margin.  New OPEN-Q-108 subtarget: classify finite `++++++` packet alignments, then prove the post-cutoff tail has packet sign cancellation or small same-sign mass. -> HYP-2674, HYP-2673, HYP-2671, HYP-2653d, OPEN-Q-108.

**CORRECTION ADDENDUM (codex-2026-06-20-S46, HYP-2673/HYP-2653d):** the global far-peel currency is now corrected: do **not** close OPEN-Q-108 by proving a bounded `C(k)=sup w|Delta_w|`.  HYP-2653d shows `Delta_w` has a small nonzero resonant floor along dyadic families, so `w*Delta_w` grows with scale.  The correct far-tail object is the uniform cap `sup_{max(E')>B} Delta_w(E',w) <= cap_k-Q(k-1)`, with `max(E')<=B` checked finitely.  KPS reports `B=14` already below margin but tight at k=9, `B=20` about `2.3x` safe, and the tight B14 row is exactly the codex HYP-2671 dyadic block `(0,1,2,4,8,12,16,20), w=24`.  The old finite-span numbers (`15`, `23`, `30`) remain useful diagnostics for why the stale normalization looked plausible, not proof targets. -> HYP-2673, HYP-2653d, HYP-2671, HYP-2653, OPEN-Q-108.

**UPDATE (codex-2026-06-20-S46, HYP-2673/T912):** OPEN-Q-108's "one open constant" now has a cleaner split, complementing HYP-2672's shell-full tail stratification.  The local shell-full route uses boundary-tax currency `Delta_w^+/p1(E')`: shell damage threshold `426/35035`, finite packet tax `2/5` (B13 leader gap `139/12810`), new-speed tax `1/3` (dyadic-block gap `206/12957`), and the corrected HYP-2672 tail split, not the raw B30 `p1/4` scout.  HYP-2672's B36 row above `1/4` but below `3/10` forces a finite intermediate ledger plus doubled-odd exception ledger before broad post-dyadic decay.  The global far-tail route uses the corrected HYP-2653d uniform `Delta_w` cap after a finite `max(E')` cutoff.  New subtarget: prove the local packet-tax stack after HYP-2661, prove the HYP-2672 exception/tail split, and prove `sup_{max(E')>B} Delta_w(E',w) <= cap_k-Q(k-1)` past a cutoff likely near `B=20`. -> HYP-2673, HYP-2672, HYP-2671, HYP-2670, HYP-2661, HYP-2655, HYP-2653d, HYP-2653, T912.

**BRIDGE ADDENDUM (codex-2026-06-20-S46, HYP-2673/KPS):** incoming KPS work identifies the codex p1-tax quantity and KPS far-plateau deviation exactly.  S46 verifier confirms `raw_wdelta(E',w)/w = p0(E' union {w})-Phi(E')`: at the dyadic-block extremizer both sides are `457/3920`, while the non-shell-full warning row `(0,2,3,5,6,15), w=18` has `Delta/p1=22/63>1/3`.  Therefore `Delta<=p1/3` is truly a shell-full theorem; non-shell rows must use the HYP-2661 shell-damage gate or the corrected uniform `Delta_w` far-tail route. -> HYP-2673, HYP-2671, HYP-2653d, HYP-2653, HYP-2661, OPEN-Q-108.

**UPDATE (codex-2026-06-20-S45, HYP-2671/T910):** the post-shell-gate "one open constant" is now localized as the shell-full new-speed `1/3` barrier.  In the B30 quotient, the new-speed maximum is `1371/4319` at `E'=(0,1,2,4,8,12,16,20), w=24`, with exact gap `206/12957` below `1/3`; this is the `m=4` spike of the family `E_m={0,1,2,4,8,3m,4m,5m}, w=6m`, while `m=3,5,6,...,24` are far lower.  Fold mass is a locator but not a certificate.  New subtarget: prove the dyadic block resonance is the only new-speed row near `1/3`, then prove all other shell-full new-speed rows have phase-packet cancellation below `p1/3`.  HYP-2672 supersedes the provisional far-tail `p1/4` guess. -> HYP-2672, HYP-2671, HYP-2670, HYP-2669, HYP-2668, HYP-2661, HYP-2643, T910.

**UPDATE (codex-2026-06-20-S45, HYP-2672/T911):** the shell-full far-tail constant needed correction.  Extending the exact quotient to B36 (`39680` rows) refutes the naive HYP-2670 target `max(E')>24 => Delta^+ <= p1/4`: the row `(0,1,2,4,8,14,26,34), w=38` has `Delta^+/p1=966562/3357319 > 1/4`.  It is still below `3/10` by `406337/33573190`, and all B36 rows with `max(E')>20` are below `3/10` once HYP-2671's dyadic block is treated separately.  The only post-dyadic rows above `1/4` are this doubled-odd tail packet plus two intermediate finite rows `(0,1,2,4,8,12,13,21), w=24` and `(0,1,2,4,8,14,20,24), w=26`.  A focused doubled-odd scout over `2912` rows `E'={0,1,2,4,8,2a,2b,2c}`, odd `a<b<c<=29`, found the same `(7,13,17;19)` packet as the unique row above `1/4` and no rows above `3/10`.  New subtarget: finite high pocket + HYP-2671 dyadic block + finite intermediate ledger + doubled-odd exception ledger + broad post-dyadic `3/10` decay. -> HYP-2672, HYP-2671, HYP-2670, HYP-2669, HYP-2661, T911.

**UPDATE (codex-2026-06-20-S44, HYP-2670/T909):** OPEN-Q-108's shell-full p1-tax half now has a sharper packet-gap formulation.  HYP-2670 scans the exact shell-1-full quotient `E'={0}+{1,2,4,8}+3` extras from `[1,30]`, `w=max(E')+1..max(E')+8`, for `20800` primitive rows.  The only row above `1/3` remains the B13 leader `(0,1,2,4,6,7,8,10), w=12`, `Delta^+/p1=997/2562`; every row with `max(E')>14` is below `1/3` with max `1371/4319` and gap `206/12957`, and B30 saw every row with `max(E')>24` below `1/4`.  HYP-2672 later refutes that last target at B36 and replaces it with an exception-ledger plus broad `<3/10` tail route.  New subtarget: after the HYP-2661/HYP-2666 shell-damage gate, prove a finite shell-full packet ledger for `max(E')<=14`, a new-speed `p1/3` decay lemma, and the HYP-2672 corrected tail split.  Concurrent KPS S17/THM-545 work strengthens the shell-damaged side, with k=1 tower deletion proved and k=2 wide scans reporting zero sub-threshold rows, so the two-gate route is now visibly splitting into a nearly closed local gate plus this post-gate packet tax. -> HYP-2670, HYP-2672, HYP-2669, HYP-2668, HYP-2667, HYP-2666, HYP-2661, HYP-2643, T909.

**UPDATE (codex-2026-06-20-S43, HYP-2669):** OPEN-Q-108 now has a sharper two-gate boundary-currency formulation.  HYP-2665 refutes raw `p1/3`, HYP-2667 refutes raw `3p1/8`, and HYP-2668 refutes global pre-gate `2p1/5`; the surviving target is ordered.  First apply the HYP-2661 shell-1/tower-deletion gate.  Then charge far endpoint imbalance to `p1(E')`, the single-missed-sector boundary mass, on the shell-1-full quotient.  HYP-2669 scans `E'={0}+{1,2,4,8}+3` extras from `[1,24]`, `9120` primitive rows, and finds `0` shell-full violations of `2p1/5`; the stable leader is `E'=(0,1,2,4,6,7,8,10), w=12`, `Delta^+/p1=997/2562`, exact gap `139/2450`.  New subtarget: close shell-1-damaged rows by HYP-2661 mouth/tower deletion, then prove a finite dyadic-even packet lemma around the B13 leader plus new-speed decay on the shell-1-full/nonlocal quotient. -> HYP-2669, HYP-2668, HYP-2667, HYP-2666, HYP-2665, HYP-2661, HYP-2662, HYP-2663, HYP-2664, HYP-2655, HYP-2658, HYP-2648, T908.

**UPDATE (codex-2026-06-19-S39, HYP-2664):** three-tail AP-tail proof order now explicitly uses the new HYP-2661 shell-1 gate before root-packet comb enumeration. A naive unbounded nested comb proof of HYP-2663 leaves a large exact residue region, but the first-tail cutoff frontier shows most of that burden is avoidable: among `715` four-hole AP bases in `({1..13}\H) union {r,s,t}`, `589` damage shell 1 and only `126` preserve `{1,2,4,8}`; `37/40` top crude comb bases and `87/100` top bases are shell-1 damaged. Global max cutoff is `R=308` at holes `(4,5,6,13)`, missing tower bit `4`; after applying HYP-2661, the shell-1-full max is `R=239` at holes `(3,5,6,13)`. New OPEN-Q-108 subtarget: prove the HYP-2661 shell-1 deletion lemma first, then run the remaining three-tail proof on shell-1-full root packets with mouth-owner pruning. -> HYP-2664, HYP-2661, HYP-2663, HYP-2654, HYP-2659, HYP-2660, T907.

**UPDATE (codex-2026-06-19-S38, HYP-2663):** an old hidden gem from HYP-2537 now gives a concrete AP-tail coordinate system for OPEN-Q-108: near-collar perturbations should be root packets, not uniform replacement shells.  The new exact scout scans `({1..13}\H) union {r,s,t}`, `|H|=4`, `14<=r<s<t<=35`, retaining AP holes, tail insertions, Glaisher odd-shell carry, and drop-6 mouth survival/damage.  It checks `1,076,482` exact primitive rows after `24,618` comb prunes and finds `0` below the AP one-hole second threshold `426/35035`; the best row `(3,4,10,12)->(17,19,20)` has `meas(G)=4309/255255`, margin `59/12495`, no old drop-6 mouth survivor, and shell-1 damage `{1:-4}`.  This independently supports KPS HYP-2661's carry-conservation law: sub-second rows should preserve shell-1 carry `{1:15}`.  New OPEN-Q-108 subtarget: prove a shell-1/mouth-damage packet theorem saying that damaging the shell-1 tower or the four drop-6 mouths, or opening genuinely new odd-shell carry, already pays at least `426/35035`; the only possible below-second rows must be finite full-mouth templates feeding HYP-2654/HYP-2659/HYP-2660 before HYP-2658 far-survival recursion. -> HYP-2663, HYP-2661, HYP-2654, HYP-2659, HYP-2660, HYP-2658, HYP-2655, HYP-2537, THM-541, THM-543, THM-544, T906.

**UPDATE (codex-2026-06-19-S36, HYP-2658):** fixed-observer core-gap survival bridge sharpens OPEN-Q-108 after the THM-541/542/543/544 near-collar layer.  For `G_C={t:||ct||>1/14 for all c in C}`, a genuinely far speed should give the survival limit `meas(G_{B union {w}}) -> (6/7)meas(G_B)`, the fixed-observer sibling of HYP-2644's sector-cover plateau and KPS HYP-2653's decorrelation route.  Tested bounded base ledgers have large collar margin: closest one-far limit is `313/11319`, with margin `5737/294294 ~= 0.01949` over `7/858`.  The local component ledger anatomizes the THM-543 exceptional row `(1,2,3,4,5,7,8,9,11,12,13,20)`: it is the `10->20` collar graft, adding two symmetric `1/1960` endpoint-owner bubbles `7->20` and `20->7` while preserving the old collar intervals; THM-544 proves the two-replacement AP-tail layer has no below-second rows.  New OPEN-Q-108 subtarget: after the proved collar/mouth/replacement layers, prove the fixed-observer `6/7` survival recursion for genuinely far rows, retaining HYP-2648/HYP-2652 addresses and using KPS HYP-2655's multiscale warning against a naive small uniform constant. -> HYP-2658, HYP-2655, HYP-2654, HYP-2653, HYP-2651, HYP-2644, HYP-2648, HYP-2652, THM-541, THM-542, THM-543, THM-544, THM-523, T903.
**UPDATE (codex-2026-06-19-S37, HYP-2660):** the odd/distinct partition identity, Witt ghost coordinates, and tournaments/even-graphs now point to the same OPEN-Q-108 proof quotient: close free data by a gauge before scalarizing.  Glaisher binary expansion says the LRC14 single-deletion layer should be addressed by dyadic tower words; `D=x d/dx` turns the Euler product into ghost sums `m[x^m]log prod(1+x^n)=sigma_odd(m)`; tournament bits become even graphs only after adding the forced parity root.  Exact one-hole rows confirm the low collar block is tower-addressed (`drop6=7/858`, `drop12=426/35035`, `drop10=1520/63063`, `drop4=97/4004`, `drop2=11/364`) but also warn that drop `8=2^3` is a high-level even outlier (`950/21021`), so the live object is `Glaisher tower word + CRT cell + endpoint-owner ledger`, not parity alone.  New subtarget: extend this tower-defect vocabulary to THM-543/544 AP-tail rows and prove every row below `426/35035` is one of finitely many tower/bubble templates before handing off to HYP-2658 fixed-observer far survival and HYP-2659 odd-shell carry. -> HYP-2660, HYP-2659, HYP-2658, HYP-2656, HYP-2655, HYP-2651, HYP-2648, HYP-2652, THM-541, THM-542, THM-543, THM-544, T905.

**UPDATE (codex-2026-06-19-S36, HYP-2658):** fixed-observer core-gap survival bridge sharpens OPEN-Q-108 after the THM-541/542/543/544 near-collar layer.  For `G_C={t:||ct||>1/14 for all c in C}`, a genuinely far speed should give the survival limit `meas(G_{B union {w}}) -> (6/7)meas(G_B)`, the fixed-observer sibling of HYP-2644's sector-cover plateau and KPS HYP-2653's decorrelation route.  Tested bounded base ledgers have large collar margin: closest one-far limit is `313/11319`, with margin `5737/294294 ~= 0.01949` over `7/858`.  The local component ledger anatomizes the THM-543 exceptional row `(1,2,3,4,5,7,8,9,11,12,13,20)`: it is the `10->20` collar graft, adding two symmetric `1/1960` endpoint-owner bubbles `7->20` and `20->7` while preserving the old collar intervals; THM-544 proves the two-replacement AP-tail layer has no below-second rows.  KPS HYP-2656 supplies the compatible CRT/Glaisher dyadic-tower explanation for the single-deletion layer.  New OPEN-Q-108 subtarget: after the proved collar/mouth/replacement layers, prove the fixed-observer `6/7` survival recursion for genuinely far rows, retaining HYP-2648/HYP-2652 addresses and using KPS HYP-2655's multiscale warning against a naive small uniform constant. -> HYP-2658, HYP-2656, HYP-2655, HYP-2654, HYP-2653, HYP-2651, HYP-2644, HYP-2648, HYP-2652, THM-541, THM-542, THM-543, THM-544, THM-523, T903.

**UPDATE (codex-2026-06-19-S33, HYP-2652):** speed-set invariant addendum to HYP-2650 and the HYP-2651 core-gap atlas.  Exact atlas over `640` primitive 13-speed rows at the LRC14 threshold supports the stack `CRT obligation -> additive anti-coset/relation density -> safe-component boundary-owner geometry -> binding denominator readout`.  Relation-density shadows are the strongest pre-address scalar signals (`three_term_count` corr with `M` `-0.919`, longest run `-0.896`, difference energy `-0.890`, active support6+ equal-subset relations `-0.722`), while residues/q-coverage are gates, not determinants.  Caveat from CASE-thm538: this is an active additive-relation proxy, not reliance on the disputed full zero-padded `K` support-six floor.  New near-tight tail laws to prove from THM-524: `{1..12,13a}` has `M=a/(13a+1)` for tested `a=2..9`; `{1..11,13,12a}` has `M=a/(12a+5)` for tested `a>=3`, with Goddyn-Wong at `a=2` the exceptional tight seed.  This reinforces HYP-2650's answer: scalar invariants work only with their address; here the address is endpoint owner geometry plus binding denominator. -> HYP-2652, HYP-2651, HYP-2650, HYP-2646, THM-524, T899.

**UPDATE (codex-2026-06-19-S32, HYP-2650):** invariant separation now gives a sharper answer to "what determines the LRC structure": scalar plus address.  Exact max-min scout over `1743` primitive rows shows coarse summaries (`sumset_excess`, `fold_count`, `fold_profile`, `gap_pattern`, and even bare optimal denominator `q`) mix exact `M` fibers, while addressed optimal-time words such as `(q,j,active runner values)` and `(q, folded residues)` separate in the tested bank.  LRC14 sector scout over all `1287` primitive `k=9` rows shows additive/fold summaries and transition signatures mix `L_y`, while missed-count histogram/state mass/full measured state word do not.  Histogram separates only because `L_y` is a direct valuation of it; HYP-2648's measured state word is still the richer carrier for transition complexity, signed transport, fold targets, and coimage phase.  New OPEN-Q-108 subtarget: define a canonical addressed wall/crossing sheaf unifying THM-524 clearance crossings with HYP-2648 sector state words, then prove a high-`L_y` template theorem before taking scalar valuations. -> HYP-2650, HYP-2648, HYP-2647, HYP-2646, THM-524, T897.
**UPDATE (codex-2026-06-19-S32, HYP-2651; sharpened by THM-544/HYP-2658):** OPEN-Q-108 now has an exact fixed-observer bounded atlas for the THM-523 core-gap crux.  The scout `lrc14_core_gap_atlas_codex_s32.py` scans every primitive positive 12-core `C subset [1,19]` at target `1/14` (`50,388` rows) without translation normalization and finds the unique minimum `meas(G_C)=7/858` at the drop-6 core `(1,2,3,4,5,7,8,9,10,11,12,13)`.  The second distinct `B<=19` value is `426/35035` at `(1,2,3,4,5,6,7,8,9,10,11,13)`, separated by `841/210210`.  THM-541 proves the single-hole collar, THM-542 proves one-tail mouth retention, THM-543 proves the one-replacement AP-tail layer with unique sub-`426/35035` exception `10->20`, and THM-544 proves the two-replacement AP-tail layer has no below-second rows; HYP-2658 records the component-owner bubble ledger and points the remaining genuinely-far rows to fixed-observer `6/7` survival plus HYP-2655 joint plateau/Delta recursion.  Guardrail: fixed-observer `G_C` is not freely translation-invariant, so Freiman normal forms are classifiers only until the predicate is preserved. -> HYP-2651, HYP-2658, HYP-2655, HYP-2654, HYP-2653, THM-541, THM-542, THM-543, THM-544, THM-523, HYP-2649, HYP-2648, HYP-2644, HYP-2638, T898, T903.
**UPDATE (codex-2026-06-19-S32, HYP-2651; sharpened by THM-544/HYP-2658 and KPS HYP-2656):** OPEN-Q-108 now has an exact fixed-observer bounded atlas for the THM-523 core-gap crux.  The scout `lrc14_core_gap_atlas_codex_s32.py` scans every primitive positive 12-core `C subset [1,19]` at target `1/14` (`50,388` rows) without translation normalization and finds the unique minimum `meas(G_C)=7/858` at the drop-6 core `(1,2,3,4,5,7,8,9,10,11,12,13)`.  The second distinct `B<=19` value is `426/35035` at `(1,2,3,4,5,6,7,8,9,10,11,13)`, separated by `841/210210`.  THM-541 proves the single-hole collar, THM-542 proves one-tail mouth retention, THM-543 proves the one-replacement AP-tail layer with unique sub-`426/35035` exception `10->20`, and THM-544 proves the two-replacement AP-tail layer has no below-second rows; HYP-2658 records the component-owner bubble ledger and points the remaining genuinely-far rows to fixed-observer `6/7` survival plus HYP-2655 joint plateau/Delta recursion, while KPS HYP-2656 explains the same single-deletion layer as an odd-base/dyadic-refinement tower.  Guardrail: fixed-observer `G_C` is not freely translation-invariant, so Freiman normal forms are classifiers only until the predicate is preserved. -> HYP-2651, HYP-2658, HYP-2656, HYP-2655, HYP-2654, HYP-2653, THM-541, THM-542, THM-543, THM-544, THM-523, HYP-2649, HYP-2648, HYP-2644, HYP-2638, T898, T903.
**UPDATE (codex-2026-06-19-S31, HYP-2648):** the next invariant layer is the measured cyclic state word `W(E)={(I,|I|,M_E(I))}` of missed inner seventh-sectors on wall atoms.  This contains the current scalar shadows as quotients: `p_t`, `L_y`, and `p0+p1/7` are histograms; HYP-2647 signed transport is a common-refinement coupling of two state words; HYP-2643 fold targets and HYP-2646 mod-7 reciprocal phases are addresses retained before valuation.  The S31 scout reproduces the AP9 -> `(0,1,2,3,4,5,6,7,9)` wall-transfer shadow exactly (`positive=9749/158760`, `negative=2659/39690`, `signed D-AP=-887/158760`) while showing the signed drop lives inside `76` state-pairs with `4393/5880` neutral mass.  New OPEN-Q-108 subtarget: prove a high-`L_y` state-word template theorem.  Rows matching the near-AP template should reduce to the HYP-2647/HYP-2643 transport lemmas; rows with higher state entropy/transition complexity should be forced into HYP-2638 Freiman small-excess, HYP-2639 relation-covered GAP slack, HYP-2646 signed coimage cancellation, or HYP-2644 far-element plateau contraction.  Tournament vertices should be state addresses/proof obligations, not raw runners or arcs. -> HYP-2648, HYP-2647, HYP-2646, HYP-2644, HYP-2643, HYP-2639, HYP-2638, T895.
**UPDATE (codex-2026-06-19-S31, HYP-2649):** the below-14 reproof target now has a modern training ladder.  Exact scout `lrc_below14_modern_reproof_codex_s31.py` verifies that AP rows `(1,...,N-1)` are tight for `N=3..13`, and in the one-step AP frontier (primitive size `N-1` subsets of `[1,N]`) AP is the unique tight/min row for every `N=4..13`; relaxing to `1/(N+1)` gives positive safe measure, with top-edge value `7/858`.  All `91` primitive 12-subsets of `[1,14]` satisfy `M>=1/13`, with unique tight row `(1..12)` and minimum strict safe measure at target `1/14` equal to `7/858` at `(1,2,3,4,5,7,8,9,10,11,12,13)`.  New OPEN-Q-108 subtarget: turn this into an AP-frontier fattening lemma, then extend from `[1,N]` to all AP-rich normal forms via Freiman/signed-wall transport.  The support-floor ladder explains why `N=14` is qualitatively new: it is the first even row with support floor `6`, where HYP-2646/HYP-2647 signed coset machinery becomes necessary. -> HYP-2649, HYP-2647, HYP-2646, THM-523, THM-538, T896.

**UPDATE (codex-2026-06-19-S30, HYP-2647):** the HYP-2642 weighted positive/negative wall-transfer ledger and HYP-2643 fold-target transport are now unified as a signed transport matrix on the common wall tiling, aligned with KPS HYP-2646's exact signed coset/reciprocal factorization.  For AP9 versus `(0,1,2,3,4,5,6,7,9)`, the scalar shadow is `positive=9749/158760`, `negative=2659/39690`, surplus `887/158760`; the hidden fold coordinate is the address shift `3/8 -> 3/9`.  The S30 scout keeps the moving endpoint sector address and verifies exact transport balance: old-sector row masses and new-sector column masses are all `1/7`, while pre-weight atom mass is `274/2205` positive, `2269/17640` negative, and `4393/5880` neutral.  New OPEN-Q-108 subtarget: define the addressed wall-transport matrix with buckets `(missed-sector set, N, fold target, sign type, residue/coimage shell)`, prove its row-balance checksum, and then prove that inside the k=9 near-AP transport class the endpoint defect `s=2` maximizes `L_y`; rows outside the class should route to HYP-2638 Freiman small-excess or HYP-2639 signed-shell cancellation.  This imports the tournament arc-flip lesson: do not scalarize positive/negative signs before retaining contraction/fold/residue address. -> HYP-2647, HYP-2646, HYP-2643, HYP-2642, HYP-2639, HYP-2638, T894.

**UPDATE (codex-2026-06-19-S29, HYP-2642):** the corrected KPS S12 binding bounded case has an exact wall-transfer certificate.  For `A=(0,1,2,3,4,5,6,7,8)` and `D=(0,1,2,3,4,5,6,7,9)`, common-wall refinement gives `L_y(A)-L_y(D)=887/158760=0.005587050`, while `cap_9-L_y(D)=39541/5675670=0.006966755`.  The AP-to-defect loss is a surplus of weighted negative wall transfers over positive ones: `2659/39690 - 9749/158760 = 887/158760`.  A one-gap scan through `s<=30` finds the endpoint row `F_s=(0,1,2,3,4,5,6,7,7+s)` best for every tested gap and the minimal gap `s=2` as the global top.  Endpoint rows have independent-extra-runner `L_y` limit `20698/46305`, with `F_2` higher by `22391/555660`; a proof of `|L_y(F_s)-20698/46305| <= 1/(7+s)` would reduce endpoint gaps to the finite check `s<=17`.  New OPEN-Q-108 subtarget: replace the tight finite near-AP check by three structural lemmas -- endpoint dominance, a discrepancy/residue envelope proving `F_s <= F_2`, and an AP-to-defect transfer pairing retaining at least the AP-to-cap slack `10441/7567560`.  This is complementary to KPS HYP-2641's far-element p0 plateau recursion. -> HYP-2642, HYP-2641, HYP-2638, HYP-2640, T890.
**UPDATE (codex-2026-06-19-S29, HYP-2641):** the corrected KPS S12 binding case has an exact wall-transfer certificate.  For `A=(0,1,2,3,4,5,6,7,8)` and `D=(0,1,2,3,4,5,6,7,9)`, common-wall refinement gives `L_y(A)-L_y(D)=887/158760=0.005587050`, while `cap_9-L_y(D)=39541/5675670=0.006966755`.  The AP-to-defect loss is a surplus of weighted negative wall transfers over positive ones: `2659/39690 - 9749/158760 = 887/158760`.  A one-gap scan through `s<=30` finds the endpoint row `F_s=(0,1,2,3,4,5,6,7,7+s)` best for every tested gap and the minimal gap `s=2` as the global top.  Endpoint rows have independent-extra-runner limit `20698/46305`, with `F_2` higher by `22391/555660`; a proof of `|L_y(F_s)-20698/46305| <= 1/(7+s)` would reduce endpoint gaps to the finite check `s<=17`.  New OPEN-Q-108 subtarget: replace the tight finite near-AP check by three structural lemmas -- endpoint dominance, a discrepancy/residue envelope proving `F_s <= F_2`, and an AP-to-defect transfer pairing retaining at least the AP-to-cap slack `10441/7567560`. -> HYP-2641, HYP-2638, HYP-2640, T889.
**UPDATE (codex-2026-06-19-S29, HYP-2643):** the fold-multiplicity route now has a sharper invariant, complementary to HYP-2642's wall-transfer certificate for the same k=9 endpoint defect. Total visible fold count is too coarse: AP9 and the binding near-AP row `(0,1,2,3,4,5,6,7,9)` both have `12` nontrivial folds. The discriminant is the fold target profile `F_E(c)=#{0<a<b in E:a+b=c in E}`: near-AP moves three folds from target `8` to target `9`, giving exact reciprocal transport `3/8-3/9=1/24`. In the bounded k=9 bank `max(E)<=13`, that row is the unique top non-AP and the unique tiny positive transport bucket. New OPEN-Q-108 subtarget: prove a clipped-AP fold-transport lemma, then convert this target-profile loss into a sector-distribution loss for `L_y`/`p0`; route larger transports through the existing Freiman small-excess/dimension and signed-tail envelopes. -> HYP-2643, HYP-2642, T891.

**UPDATE (codex-2026-06-19-S28, HYP-2640):** the correction-vs-relation-rank lens now has a concrete atlas. Exact height-2 rank is a phase switch, not a scalar ruler. Ternary dissociated rows have exact rank `0` and sit near the independent baseline, while AP, nearAP, d2 GAP, and third-pocket rows already saturate exact rank. At `k=8`, AP and third-pocket A both have `rE=6`, but `L_y-L_y^inf` is `0.308965` versus `0.013547`; the visible data drop from AP fold count `12` / exact relation count `1786` to third-pocket A fold count `3` / exact relation count `326`. New OPEN-Q-108 subtarget: prove a two-stage correction lemma -- low rank or uncovered weighted fibre gives the independent peel; saturated rank invokes inverse structure, and every non-AP saturated row must lose observer-fold multiplicity or signed mod-27 coimage alignment before HYP-2636/HYP-2633 tail scalarization. -> HYP-2640, T888.

**UPDATE (codex-2026-06-19-S26, HYP-2639):** after HYP-2637's weighted relation-fiber/GAP split, HYP-2638's Freiman small-excess finite pocket, and KPS S12's sumset-excess/Freiman-dimension scan, the HYP-2635 additive-energy lead has a sharper LRC guardrail. Direct "every element in a small relation -> one-dimensional Freiman GAP" is too strong for the third pocket: KPS relation-covered examples have every nonzero element in a small motif but `|S+S|=31` for `k=8`, above `3k-4=20`. AP versus shifted AP shows why raw energy is a decoy: both have `|S+S|=25` and energy `1469`, but AP has `36` observer-coupled visible folds and `M*n=1`, while shifted AP has `0` visible folds and `M*n=4.789`. New OPEN-Q-108 subtarget: prove a labelled relation-hypergraph regularity lemma using summand nodes `C=a+b`, multiplicand tests `C|w`, and sign type (balanced/even scalar vs observer-coupled/odd marked), then show every non-AP relation-covered pocket has Freiman-dimension slack or signed shell cancellation before absolute values. -> HYP-2639, HYP-2638, HYP-2637, HYP-2636, HYP-2635, HYP-2634, T887.

**UPDATE (codex-2026-06-19-S25, HYP-2634):** the HYP-2633 opposite bounded signs now have a first structural explanation. In the seed family `S_a=(1,a,8,a+7,15,22)`, finite HYP-2632 packets give QR weight `-25U` for both `a=2` and `a=4`, but the actual bounded reciprocal lift is negative only at `a=4`. The reason is an exact low-height defect sieve: every `S_a` has a universal positive height-2 relation, while only `a=4` has larger negative resonances with defects `2a-8` and `7a-28`. New OPEN-Q-108 subtarget: before proving residue-lift equidistribution / Abel summation, isolate finite low-height relation-defect zeros as a wall ledger; equidistribution should be required only on the residual. -> HYP-2634, T882.

**UPDATE (codex-2026-06-22-S101, HYP-2883):** the HYP-2632 finite packet now has a stronger exact form: a locally balanced signed graph.  On residue vertices `{0,2,3,4,5,6}`, put negative `4+2` loop weights `-4,-25,-18,-25,-18,-18` and `4+1+1` edge weights in `{0,1,8}`.  The zero edges are exactly the affine matching `a+b=2 mod 7`, and off that lane high/low is controlled by `chi_7(Q(a,b))`, `Q=ab*(1+3(a+b))-1`.  Exact audit `lrc14_repeated_packet_graph_codex_s101.py` verifies `loop(a)+sum_b edge(a,b)=0` at every vertex: scalar counting leaves `-108U+54U=-54U`, but incidence counting is perfectly conserved.  New OPEN-Q-108 subtarget: lift this local signed-current identity, not just the scalar signed ledger, through the reciprocal hyperplane sums after finite wall deletion. -> HYP-2883, HYP-2882, HYP-2632, HYP-2636, HYP-2617, T999.

**UPDATE (codex-2026-06-22-S102, HYP-2884):** the first lift test turns the HYP-2883 local balance into a precise divergence-defect obligation.  Using core `(1,8,15,22)`, `lrc14_packet_balance_lift_probe_codex_s102.py` compares start-aligned and raised-pair integer lifts of the repeated packet graph.  Through `H=12`, finite balance remains exact, but reciprocal lifts have nonzero vertex divergence.  Start-aligned `H=12`: `max|div|=0.00512112`, `L1 div=0.0193444`; raised-pair `H=12`: `max|div|=0.00191161`, `L1 div=0.00610376`.  New OPEN-Q-108 subtarget: delete coefficient-height `<=2` wall directions, then prove the lifted divergence is Abel-summable inside HYP-2636 additive-frequency shells. -> HYP-2884, HYP-2883, HYP-2633, HYP-2634, HYP-2636, T1000.
**UPDATE (codex-2026-06-22-S102, HYP-2886):** the exact-period witness side now has a packet atlas that preserves the same mod-7/affine layer.  For denominator `D`, a unit `a mod D` is safe iff `14*min(sa mod D,D-sa mod D)>=D` for every speed; the capacity is `phi(D)`, but the safe packet distribution must be read before scalarizing to `N(S,D)`.  Fixed finite bases are confirmed as charts, not closures: `divload_B90={1,...,11,13,84*lcm(1,...,90)}` kills `D=21,41,53,83,89` and first opens at `D=97`.  In mixed cases, quotient signal ranks `mod14 > mod7 > chi_7 x affine_pair > affine_pair > chi_7 > parity`, showing that HYP-2632/HYP-2883/HYP-2884's signed-current layer remains visible in actual rational witness packets.  New OPEN-Q-108 subtarget: after removing divisibility-killed denominators, lift the local signed-current balance on exact-period residue fibers and retain CRT multiplicativity defects as labelled atoms before the spectrum/L2/Part-A floor. -> HYP-2886, HYP-2884, HYP-2883, HYP-2882, HYP-2876, HYP-2865, HYP-2632, T1001.

**UPDATE (codex-2026-06-22-S103, HYP-2887):** the repeated-packet local-current defect now has a finite realizability carrier, complementary to HYP-2885's additive-energy extremality and HYP-2886's exact-period atlas.  The HYP-2632 nonzero graph is `K6` minus the affine-zero matching, equivalently the octahedron `L(K4)`: residues are tetrahedron edges and affine-zero pairs are opposite edges.  The octahedron cycle rank is `7` (eight triangular face curls with one dependence), giving a concrete apex-7 face-curl module.  Exact layer-cochain scan `lrc14_octahedral_current_realizability_codex_s103.py` over all `3^6=729` gauges at `H=10` finds best gauge `210210` with `L1 div=0.00225361` versus constant `000000` `0.0219283` and `111111` `0.00754806`; wall max correlates with divergence and mixed local signs correlate with cancellation.  New OPEN-Q-108 subtarget: prove an octahedral Hodge lemma for the lifted packet current after height-2 wall deletion: coherent triangular face curl is finite wall structure, incoherent spread current is HYP-2636 Abel-summable, and no harmonic one-current remains on the spherical carrier. -> HYP-2887, HYP-2885, HYP-2886, HYP-2884, HYP-2883, HYP-2636, HYP-2633, T1002.
**UPDATE (codex-2026-06-22-S103, HYP-2889):** the additive-energy majorization branch is now corrected.  Exact scout `lrc_additive_energy_majorization_codex_s103.py` refutes scalar monotonicity (`3137` p0 inversions), pairwise difference-profile monotonicity (`556752` violations), and naive one-step compression (`1177` profile-up moves with p0 down).  What survives is anchored AP-facing: the interval difference profile `(k-1,...,1)` majorizes every tested row in banks k=8 (`1716` rows), k=9 (`1287`), k=10 (`220`), and AP has `0` p0/L_y over-beaters.  Concurrent HYP-+2888/S39 pins the strict-threshold endpoint and refines it: exact tiling is scaling-invariant `d*{1,...,13}`, not arbitrary AP translates with the same additive energy.  KPS S31l also shows the higher additive-moment tail is mixed-sign, so the cap proof must prevent non-AP over-covering by signed cancellation rather than scalar energy monotonicity.  New OPEN-Q-108 subtarget: prove anchored interval difference-profile majorization, then split THM-534's `L_y` into an AP-facing Fejer part plus a labelled signed sector/Fourier remainder. -> HYP-2889, HYP-2885, HYP-+2888, HYP-2873, THM-534, T1003.

**UPDATE (codex-2026-06-19-S24, HYP-2633):** the HYP-2632 finite packet now has a reciprocal-lift guardrail. Exact cumulative support-six reciprocal sums through `H=16` on representative two-large supports show that finite `chi_7`/affine/Q packet signs do not yet control bounded analytic tail signs. Two QR packets with the same finite weight `-25U` lift to opposite signs (`42_QR_a2` has `+0.002676143676`, `42_QR_a4` has `-0.000130513735`); `411_low_26` has finite `+1U` but negative lift; and affine-zero `411_zero_02` has finite `0U` but nonzero bounded lift `+0.000411593461`. New OPEN-Q-108 subtarget: prove residue-lift equidistribution / Abel summation inside additive-frequency shells after low-height wall deletion, so the finite `-108U+54U` signed ledger can be used against the actual reciprocal hyperplane tail. -> HYP-2633, T881.
**UPDATE (codex-2026-06-19-S24, HYP-2636):** after HYP-2633's reciprocal-lift guardrail, HYP-2634's lift-opposition atlas, and HYP-2635's HYP-2607 frontier consolidation, the HYP-2632 two-large character packet now has an exact block-transfer skeleton for the reciprocal tail. For model faces `c1*n1+...+c4*n4+A*x+B*y=0`, group the support-six correction as `sum_s <Core_s(u,v), Pair_s^{A,B}(u,v)>` over exact additive channels `s` and the `6 x 6` residue matrix before applying absolute values. At `H=24`, this sharply reduces the visible envelope: QR/NQR `4+2` rows have block `L1/signed` near `1.05-1.11` versus raw/signed near `21`, and the affine zero-lane `4+1+1` row drops from raw/signed about `1420` to block `L1/signed` about `18.9`. Same-residue spread-core probes `(1,8,15,22)` preserve a weaker but real collapse: block `L1/signed` around `14-17` versus raw/signed `185-302`. New OPEN-Q-108 subtarget: prove the pair-line Dedekind/cotangent bound for `A*x+B*y=-s` channels and combine it with channelwise Cauchy/Schur over `sum_s ||Core_s||_2||Pair_s||_2` after finite wall deletion. This is the tail-side analogue of KPS S11's lesson to keep the full empty-sector distribution before scalarizing. -> HYP-2636, HYP-2635, HYP-2634, HYP-2633, T884.

**UPDATE (codex-2026-06-19-S23, HYP-2632):** the HYP-2630 two-large repeated tail now has a finite signed character kernel. Exact computation verifies `S_d(a)=sum_{a.r=0}C_d(r)=(1/7)sum_t C_hat(ta)` over all `159` projective support-six coimage classes with max error `1.10e-14`. In packet units `U=0.00955648353590534`, the `4+2` row is exactly Legendre: `S/U=-25` for QR `a=2,4`, `-18` for NQR `a=3,5,6`, equivalently `2S/U=-43-7chi_7(a)` for `a=2..6`. The `4+1+1` row has six `+8U`, six `+U`, and three zeros; the zero lane is affine, `a+b=2 mod 7`, and off that lane high/low is controlled by `Q(a,b)=ab*(1+3(a+b))-1`, with `S/U=8` iff `chi_7(Q)=+1` and `S/U=1` iff `chi_7(Q)=-1`. Companion Dedekind-shell computation expands the same kernel as explicit `D_T(ell)` factors and shows the blind two-large residue matrix misses exact zero rows. New OPEN-Q-108 subtarget: prove the two-large reciprocal hyperplane tail by exposing additive frequency/conjugate shells before absolute values, using the signed `chi_7`/affine/Q table and the repeated-kernel ledger `-108U+54U=-54U` instead of the `162U` absolute mass. -> HYP-2632, T880.

**UPDATE (codex-2026-06-19-S22, HYP-2631):** the HYP-2628 `Q=210 -> Q=1260` AP one-drop repair is now explained at the reduced-denominator level. Exact computation shows the two `Q=210`-blind AP drops have raw `Q=1260` strict-safe residues only in exact-period packets whose reduced denominators do not divide `210`: drop `6` uses `63,420,630`; drop `12` uses `12,315,630,1260`. Each component hit has the omitted AP speed as the only full-AP danger, so these are genuine transversality samples rather than accidental grid points. Caveat: drop `6` has a `q=98` witness outside the raw Hill carrier; the theorem-shaped claim is about the packet ledger retained by `1260=2^2*3^2*5*7`. New OPEN-Q-108 subtarget: extend the reduced-denominator mouth ledger from AP drops to the HYP-2626 repeated coimage tail, and only then project to squarefree masks / mod-7 coimage. -> HYP-2631, T879.
**UPDATE (codex-2026-06-19-S22, HYP-2630):** the HYP-2629 Euler-copy tail test is now a residue/phase split. Exact-period copy mass is uniform over the LRC14 unit seam: for raw `q=1260`, exact top-period packets give `48` copies per unit residue and the full `{2,3,5,7}` mask gives `96` per unit residue. Thus copy capacity thickens repeated-residue multiplicity patterns but cannot separate QR/NQR. Height `3` one-large wall enumeration still hits the same `85/116` nonzero `k=10` coimage classes as height `2`, so the `31` tail-only classes are not a missed one-large wall-height artifact. They are structurally multi-large: `94.382022%` of tail mass needs at least two large residue coordinates after the core `1..13`, and the `4+2` packet has identical copy capacities but different quadratic-character masses (`QR mean |S_9|=0.23891209`, `NQR mean |S_9|=0.17201670`). New OPEN-Q-108 subtarget: prove a two-large repeated-residue cotangent/Dedekind bound retaining the `chi_7` phase channel; stop raising one-large wall height as the main route. -> HYP-2630, T878.
**UPDATE (codex-2026-06-19-S21, HYP-2628):** the squarefree-profile route now has a canonical exact-period measure. The user's copy rule `sum_{d|n} c(d)=n` is `1*c=id`, hence `c=mu*id=phi`; on a `q`-grid it is exactly the reduced-denominator census, with `phi(d)` residues of exact denominator `d|q`. For LRC14, `phi(14)=6` is the HYP-2626 unit seam, and the raw Hill denominator `1260=2^2*3^2*5*7` should be projected through exact-period packets before squarefree compression. The full `{2,3,5,7}` mask mass of `1260` is `576`, equal to `phi(2520)`, the THM-523 half-clash denominator packet count before symmetry doubles `1/2520` to `1/1260`. New safe-center transfer row: `Q=210` catches `11/13` AP one-drop cores and misses exactly drops `6,12`; raw `Q=1260` catches all `13/13`, while full AP13 still has no strict-safe residue. New OPEN-Q-108 subtarget: explain this `210 -> 1260` repair and rewrite the HYP-2625/HYP-2626 transfer as an exact-period phi-packet matrix before mod-7 coimage projection. -> HYP-2628, T876.
**UPDATE (codex-2026-06-19-S21, HYP-2629):** after HYP-2628's exact-period packet law, the Hill-row scan isolates the crossing quotient loss. The squarefree copy profile `copy_mass_N(M)=sum_{d|N, mask(d)=M} phi(d)` has a prime-extension recurrence: adding prime `p` appends a shifted copy layer multiplied by `p-1`, so `{2,3}->{2,3,5}->{2,3,5,7}` is a genuine Euler-copy recurrence. At `K_14`, raw `P_14=1260` has full `{2,3,5,7}` copy mass `576`, while `cr(K_14)=315` has full copy mass `0`; division by `4` deletes the dyadic gate. New OPEN-Q-108 subtarget: re-index the HYP-2626 k=10 repeated-residue tail by Euler-copy mask mass and test whether it separates the quadratic-character packets beyond the raw `{7}` mask. -> HYP-2629, T877.

**UPDATE (codex-2026-06-19-S20, HYP-2627):** the squarefree divisor-profile route now has a complete-graph crossing source. For Hill's crossing product `P_n=floor(n/2)floor((n-1)/2)floor((n-2)/2)floor((n-3)/2)=4cr(K_n)`, `n=14` is exactly the first row with `rad(P_n)` divisible by `210`, and `q_14=(5,6,6,7)` has `P_14=1260=2^2*3^2*5*7`, `rad(P_14)=210`, `cr(K_14)=315`. This packages the LRC14 hierarchy: repeated `6` = mod-6 skeleton, `5` = mod-30 address, `7` = HYP-2626 coimage seam. Markov-Hurwitz `w^2+x^2+y^2+z^2=wxyz` is the recurrence archetype but not the direct carrier: `q_14` has pressure `73/630`, and generated positive Markov-Hurwitz solutions through max coordinate `10^8` have no coordinate divisible by `5`. New OPEN-Q-108 subtarget: re-index the repeated-residue HYP-2626 tail by the four-block squarefree profile `(5,6,6,7)` and test whether its character split becomes a crossing/Hurwitz pressure inequality. -> HYP-2627, T875.

**UPDATE (codex-2026-06-19-S20, HYP-2627):** the squarefree-profile clue now has a raw four-factor bridge. The direct Markov-Hurwitz/crossing identity is false: Harary-Hill tuples do not solve the normalized equation `a^2+b^2+c^2+d^2=4abcd` (`n=14` tuple `(5,6,6,7)` has defect `4894`). But the quotient object is live: the raw Harary-Hill product for `K_14` is `7*6*6*5=1260`, with squarefree core `210={2,3,5,7}`, exactly the HYP-2625/HYP-2626 mod-210 profile and the THM-522/HYP-2561 lonely-measure scale `1/1260`. The divided crossing value `315` loses the dyadic gate. New OPEN-Q-108 subtarget: derive the `1/1260` two-speed-clash denominator from the raw four-factor row `7,6,6,5`, or prove this is only a denominator coincidence. -> HYP-2627, T875.
**S20 addendum:** the denominator has now been derived algebraically: THM-523's half-clash is `15/36-2/5-1/70-1/504=1/(2*1260)`, doubled by symmetry to `1/1260`. The remaining question is geometric/proof-theoretic: why does the `12->36` perturbation select exactly the raw row `7,6,6,5`?
**UPDATE (codex-2026-06-19-S19, HYP-2626):** the support-six coimage target now has a prime-mask transfer coordinate. The exact seam law `(Z/14Z)^* -> F_7^*` shows that HYP-2617's projective mod-7 quotient is the LRC14 unit-action coimage. In the height `<=2` one-large wall census, k=10 mask `{}` already covers `73` classes / `72.120496%` signed mass, mask `{2,3,5}` adds nothing, and mask `{7}` reaches HYP-2624's `85` classes / `84.229179%`. Thus mod30 belongs to the moving-k spectral/primorial story, while the fixed LRC14 residual needs the prime-7 seam plus the signed mod-7 coimage tail. The remaining repeated packet splits by quadratic characters: `(1,1,1,1,a,a)` separates residues `a=2,4` from nonresidues `a=3,5,6`, and `(1,1,1,1,a,b)` reduces to a short list of character signatures. New OPEN-Q-108 subtarget: prove a repeated-root cotangent/Dedekind bound by these character cases after finite height-2/multi-large wall accounting. -> HYP-2626, T874.
**UPDATE (codex-2026-06-19-S18, HYP-2624):** the LRC(14) support-six coimage target has been narrowed after height-2 wall-addressing. Enumerating one-large support-six walls with coefficient height `<=2` and projecting to the HYP-2617 mod-7 coimage atlas hits every nonzero class for `k=8` and `k=9` (`46/46`, `79/79`). For `k=10`, height `<=2` walls hit `85/116` nonzero classes and `84.229179%` of signed coimage mass; the missed `31` classes are not arbitrary but are dominated by repeated-residue packets `(1,1,1,1,a,a)` and `(1,1,1,1,a,b)` plus a small zero-cusp halo. This is a routing lemma, not a proof: next steps are exact finite height-2 wall accounting and a repeated-residue cotangent/Dedekind reciprocal-tail theorem. -> HYP-2624, T872.

**UPDATE (codex-2026-06-19-S17, HYP-2622):** the LRC-spectrum lower-bound question is now an excess/height problem. For `M(S)=p/q`, set `e=p(k+1)-q`; then `M(S)-1/(k+1)=e/(q(k+1))`, so the rows that threaten the lower bound are small-excess, high-denominator rows, not simply rows below the doubled-top mediant. The S17 AP-defect audit finds only fixed `r=4,3,2` unit-excess branches at the top of normalized depth, and no `r>=5` growth in `k<=36,r<=12` or high-`r` probes at `k=31,61`. Integrating KPS S9's denominator lemma `q<=2 max(S)` gives `g(k)>=1/(2 max(S)(k+1))`; hence any true `o(1/k^2)` dip in `sigma_2(k)` must have `max(S)/k -> infinity`. Next route: search by `(excess, height ratio)` and prove the modular upper cover for the `r=3` residue-class witnesses. -> HYP-2622, T870.

**UPDATE (codex-2026-06-19-S16, HYP-2621):** the doubled-top LRC-spectrum family is now an exact computational packet: `D_k={1,...,k-1,2k}` has `M(D_k)=2/(2k+1)` for `k=2..30`, hence the gap depth is `(k+1)(2k+1)`. The lower-bound side immediately grew AP-defect constant branches: `A_{k,3}={1,...,k}\{k-1}∪{3(k-1)}` has `M=3/(3k+2)` for tested `k≡7,13,19,25 mod30`, but is AP-tight for tested `k≡1 mod30`; on that tight class, `A_{k,4}` takes over with `M=4/(4k+3)` for stored `k=31,61,91,121,151,181`. No `o(1/k^2)` dip appears; the sharp next question is whether the AP-defect ladder admits `r=r(k)->infty` with formulas `M≈r/(rk+c)`, or whether all realized constants stay bounded. -> HYP-2621, T869.

**UPDATE (codex-2026-06-19-S15, HYP-2618):** the OCF/noise-stability analogy resolves into a packet-address rule. `H(T)=I(Omega,2)` is hard-core activity `2`, equivalently `3^m*mu_{2/3}` for the independent odd-cycle-packet indicator; it is not a nontrivial `rho<1` noise-stability functional determined by `H` alone (same `H=23` at `n=6` can have different pair/noise spectra). Exact 3-voter majority profiles on `m=3,4,5` alternatives realize every tournament, so the forbidden `{7,21}` values are forbidden Condorcet-cyclicity packet spectra. LRC(14) lesson: keep the finite packet address, then bound the signed compatible sum. Highest-value next computation after HYP-2617: classify which of the `159` support-six coimage classes remain after height-1/height-2 wall deletion for `k=8,9,10`, then prove the signed reciprocal-tail estimate on that reduced non-null class list. → HYP-2618, T866.

**UPDATE (codex-2026-06-19-S12, HYP-2616):** the first low-height ledger piece is now exact: all one-large height-1 type-II support-six walls with bounded core `C⊆{1..B(k)}` are harmless. Exhaustive primitive rows: k=8 `226046`, k=9 `250264`, k=10 `54173`; `0` cap exceedances. Worst values stay well below cap (`0.220<0.381`, `0.372<0.494`, `0.480<0.604`). The visible span-only counterexamples (`21=sum(1..6)`, k10 `22`) move from analytic tail into a finite safe ledger. The residual is now height>=2 one-large walls, multi-large low-height walls, no-scale collapse, and the HYP-2613/HYP-2614 relative signed theta tail organized by the HYP-2615 signed-mass sequence spine. Files: `04-computation/lrc14_height1_typeII_wall_ledger_codex_s12.py`, `05-knowledge/results/lrc14_height1_typeII_wall_ledger_codex_s12.out`. → HYP-2616, T864.

**LRC(14) GAP-FREE-REDUCED to ONE Minkowski lemma (kind-pasteur-2026-06-19-S9/S10, THM-538/HYP-2608a/2611):** the first open Lonely Runner case is now reduced — GAP-FREE — to a SINGLE open analytic lemma, with everything else PROVED. PROVED chain: k≤7 pigeonhole; scale-invariance; the slow/fast reduction LRC(14)⟺meas(S7(E))≤cap_k (k=8..12; glue G1 global-witness soundness sound); THM-534 per-E moment dual meas(S7)≤L_y; the caps (cap_8=2243/5880,…,cap_12=6/7, each ≥(k−6)/7); **THM-538 the support-6 floor** (K(n)=0 unless the relation has ≥6 nonzero non-7 coords — explains the HYP-2606 5×-lossiness, the signed sum annihilates all support≤5); and the **bounded-spread finite certificate** (span≤B=16/15/13, consec the unique argmax, 11432/6435/715 sets, 0 exceedances, EXHAUSTIVE). THE SINGLE RESIDUAL = **HYP-2608(a) the wide-spread bound**: span(E)>B(k) ⟹ meas(S7(E))≤cap_k, k=8,9,10. It reduces (THM-538) to bounding the support-6 correction tail eps(B)<0.17 — but the free envelope Σ c1/|m| DIVERGES harmonically, so it needs a SUCCESSIVE-MINIMA / MINKOWSKI-SECOND-THEOREM count: |K(n)|≤c1⁶/(λ₁···λ₆) over the support-6 relation lattice Λ°(E), the lattice coupling making the harmonic sum converge. THE SINGLE HIGHEST-VALUE NEXT STEP: execute that Minkowski count (it converts 0-exceedances-over-40k into a gap-free proof of LRC(14)). 0 counterexamples anywhere; LRC(14) is almost certainly TRUE but NOT proved. Stranger-contraction (HYP-2610, = kps-S9, peel ×(1−r/7)) is the moment-side decoupling. → THM-538, HYP-2608/2610/2611, THM-534/535, MISTAKE-078.

**UPDATE (codex-2026-06-19, HYP-2612):** the first finite support-6 Minkowski count has been executed as a deleted anti-coset shell census (`Lambda(E)={n:sum n_i e_i=0}`, coordinate hyperplanes and nonzero 7-cosets deleted). It DOES NOT prove LRC(14), but it refines the residual: naive span-only decay is false because wide verifier rows can have height-1 large-involving support-six identities (`21=1+2+3+4+5+6`, and the k10 `22` row has a height-1 signed identity). Dissociated strangers behave well (`{0..6,97}` first type-II shell h=5; `{0..7,68}` h=3), while no-scale clusters have height-1 relations but small exact measure. The new highest-value sub-split: (A) finite low-height anti-coset/additive-energy ledger; (B) true deleted-coset theta/successive-minima tail after those resonances are removed; (C) cluster-collapse quotient for no-scale rows. Files: `04-computation/lrc14_support6_minkowski_count_codex_20260619.py`, `05-knowledge/results/lrc14_support6_minkowski_count_codex_20260619.out`, reflection `lrc14-support-six-anti-cosets-codex-20260619.md`. → HYP-2612, T861.

**UPDATE (codex-2026-06-19-S12, HYP-2614):** the "Minkowski count" target has been sharpened again. For exact six-support terms, `K(n_1,...,n_6,0,...)=C_d(n mod 7)/(n_1...n_6)` (verified on `320` random vectors, worst error `2.948e-23`), so the residual is a finite family of residue-addressed reciprocal sums on relation hyperplanes. S12 boundary/cusp ledgers show why the absolute count is too pessimistic: hard supports have huge abs/signed separation (`(1..6)` `0.920` vs `0.0317`, `(1,2,3,4,5,21)` `0.508` vs `-0.00234`, `(2,3,4,5,6,68)` `0.100` vs `8.94e-5`). Guardrail: simple residue-coordinate marginals are nonzero, so the proof must use relation-hyperplane summation by parts plus finite wall deletion. Highest-value next step: prove the residue-cusp tail theorem, a cotangent/Dedekind-style signed theta bound, then splice it with HYP-2612's finite low-height wall ledger and cluster quotient. → HYP-2614, T862.

**UPDATE (codex-2026-06-19-S14, HYP-2617):** the coimage target is now finite. For speed residues `a_i=e_i mod 7` and Fourier residues `r_i in F_7^*`, the support-6 coimage fiber is `sum a_i r_i=0 mod 7`; quotienting by scalar and permutation leaves `159` projective speed-residue classes, with zero-speed-residue histogram `{0:80,1:42,2:22,3:10,4:4,5:1}`. Zeros are not optional: a speed divisible by `7` is Fourier-live but relation-invisible. The named support atlas shows the k=10 height-one wall class `(0,1,1,1,2,4)` is coimage-null at `d=9` even though its absolute mass is large, matching HYP-2616's finding that the visible height-one wall is finite-ledger safe. Highest-value next step: delete/account for remaining low-height wall classes, then prove a signed reciprocal-tail estimate over the non-null projective coimage fibers. → HYP-2617, T865.

**UPDATE (codex-2026-06-19-S15, HYP-2619):** the "large absolute mass but tiny signed mass" clue now has an explicit alternating sequence atlas. Residue signs only exist after conjugation pairing `r <-> -r mod 7`; paired residue totals cancel through `d=10` and then have abs/net ratios `727,174.8,71.5,38.0,24.1,16.9` for `d=11..16`. Named supports split `raw/net` into cusp-collapse and shell-alternation factors; k=9 wide and k=10 wall rows require strong second-stage shell cancellation (`variation/net=20.76,9.99`). Extending coimage classes to `d=16` shows null counts stabilize at `3`, but max coimage fiber rebounds after `d=13`, so monotone max-fiber decay is false. Highest-value next step: class-by-class signed Dedekind/cotangent tail over non-null coimage fibers after finite wall deletion. -> HYP-2619, T867.


**THE 1/7 PIVOT — LRC(14) reduced to ONE lemma (kind-pasteur-2026-06-18-S5, THM-530/HYP-2602):** the residual's correct object is NOT the via-max density ρ*_{2/7} (REFUTED: μ_{2/7} has no floor — exact k=13 sets with μ_{2/7}<1/14; ρ*_{2/7}=0 on admissible (P,E) that are nonetheless lonely) but the GLOBAL-WITNESS density ρ*_{1/7}(P,E)=meas(G_P∩{maxgap{frac(e_i x)}>1/7}) (a free fast-phase exists ⟺ the cluster phases leave a 1/7-gap ⟺ M(S)≥1/14). Two-branch floor: **k≤7 PROVED unconditional** (pigeonhole μ_{1/7}=1 ⟹ ρ*=meas(G_P)≥m_P=14249/252252); **k≥8 union bound** ρ*_{1/7}≥meas(G_P)+μ_{1/7}(E)−1>0 contingent ONLY on **HYP-2602 (the 1/7-spread bound): μ_{1/7}(E)≥thr_k** (consecutive minimizes; ≥0.32 slack; binding k=8; VERIFIED, survived the descent that killed 2/7). With the upstream finite-Vmax/integer glue (THM-527-A), HYP-2602 ⟹ LRC(14). This is the single remaining analytic gap. → THM-530, HYP-2602, HYP-2592.  **[CONVERGENCE — mac-mini-2026-06-19-S1, 8-angle workflow]:** six independent angles now reduce LRC(14)-S3 to ONE scalar lemma HYP-2607 = THM-534's `consec maximizes the empty-sector moment functional L_y(E)=E[g(N)]`. THM-534 (LP dual, PROVED per-E meas(S7)<=L_y) CLOSES the dangerous rows k=8,9,10 (L_y(consec_8)=2633/7350<cap_8); THM-537 (Beurling/U4 moment LP) converges EXACTLY on it; THM-535 collapses the finite check to those same 3 rows; HYP-2606 proves the ABSOLUTE bound is 5x lossy so the closer must be SIGNED (L_y is). HYP-2607 is NON-separable (component moments fail) — in distributional form (k=8): consec maximizes P(N=0)+(1/10)P(N=3)+P(N=6), a convex-order/coupling on the empty-sector count N. → THM-534/535/537, HYP-2606/2607.


**Residual-case S3 OPEN-Q-108 update (kind-pasteur-2026-06-18-S3 / THM-526, HYP-2581, MISTAKE-076):** the LRC(14) covering reduction's last gap (case S3 = covering 13-sets, k=|{v>13}|≥2, spread≥13×) advanced. PROVED: the **k=2 slice** (exactly 2 large speeds ⟹ M≥1/14, via drop-max scaling V₂≥51 + bounded core Vmax≥63 + finite check ≤62 = 4865 sets, worst M=1/12) and **cluster-collapse Lemma A** (window W_K safe for the whole cluster, nonempty iff 13Vmin−Vmax>14Ks; closes single-gap clusters). CORRECTION (MISTAKE-076): the criterion C(S)=∃v:W(S\{v})>1/(7v) is SUFFICIENT but **NOT universal** — S\*={1,2,3,5,7,8,9,10,11,12,13,38,42} (covering S3, k=2) has C fail for all v yet M=2/23; so "prove C for all covering S" cannot close LRC(14), and a bounded-speed reduction is REFUTED (S3 infinite: AP family {t,2t,…,12t,V}). The residual is now UNIFIED and asymptotically TIGHT: the criterion margin's V₀→∞ limit-infimum is EXACTLY 1, lifted to a realized floor >1 by covering+primitivity discreteness — no compactness argument; closing it = a uniform positive density floor (Weyl/three-distance, route a) or multi-band CRT placement (route b). NEXT: the rigorous ρ*(Δ,P)≥c₀>0 equidistribution lemma. → THM-526, HYP-2581d. **ADVANCE (mac-mini-2026-06-18-S1, THM-527, HYP-2584..2586):** the slow-fast change of variables φ=frac(Vmax·τ) reformulates route (a) as a clean single-variable **lonely density** ρ*(P,E)=meas{x∈G_P : the k cluster phases {frac(e_i x)} (e_i=Vmax−u_i) have circular max-gap >2/7 ⟺ fit in a 5/7-arc}, with ρ*>0 ⟹ M(S)≥1/14 (PROVED in the w0→∞ limit; ρ_K→ρ* VALIDATED). **The "no-compactness/asymptotically-tight" obstruction DISSOLVES**: it is about the gap-WIDTH (margin→1), but the proof needs the gap-MEASURE ρ*, which IS compact — the extremal cluster has BOUNDED SPREAD (huge spread RAISES μ; verified), so the shape space is finite-dimensional/compact and inf ρ* is a positive minimum (no ρ*=0 found; consecutive-cluster exact floor 1/84). PROVED: k=3 (margin 4/3). CORRECTION: consecutive is NOT the globally extremal shape (HYP-2585). The isolated remaining crux = the rigorous uniform floor c₀>0 on the compact (P, bounded-spread shape) space. Files: 04-computation/lrc14_{rho_star_limit,threedistance_floor,exact_floor,spread_floor,broadfloor}_macmini_0618s1.py.
**Cluster-size split OPEN-Q-108/109 update (kind-pasteur-2026-06-18-S4 / THM-527, HYP-2581e/f):** the S3 residual sharpened. CONVERGENCE: mac-mini (THM-527 reservation) and kind-pasteur (slow-fast/offset-fit) independently reached the same reformulation [good x ⟺ x∈G_P AND the cluster phase-points {frac(e_i x)} leave a circular gap > 1/7 (global witness) / > 2/7 (via-max criterion)]. PIGEONHOLE: maxgap≥1/m ⟹ margin≥7/m−1; AUTOMATIC for |L|≤6 (global witness) / |L|≤4 (criterion); |L|≥7/|L|≥5 = the ρ* hard case. PROVED this session: THM-527 (fixed-small-part single-tight-cluster, explicit V0*, global-witness θ-sweep); AP-family {1..12,m} (M≥2/27, but k=1=S1); ALL-MULT7-LARGE window-collapse (conditional on w_max<14 w_min). The 7-adic angle is REFUTED as the uniform-floor mechanism (floor = small THM-524 binding pairs). REMAINING OPEN = the COORDINATED-GROWTH CORE (k≥3, no fixed bounded small part, exemplar {t,2t,…,12t,V}), asymptotically tight (M floors at 2/23 from above, criterion-margin limit-inf=1). Finish line = uniform ρ*(Δ,P)≥c0>0 (three-distance/Weyl) OR (THM-524) covering forbids binding denominator D=14q−r with small r. → THM-527, HYP-2581d.

**Residual-case S3 OPEN-Q-108 update (kind-pasteur-2026-06-18-S3 / THM-526, HYP-2581, MISTAKE-076):** the LRC(14) covering reduction's last gap (case S3 = covering 13-sets, k=|{v>13}|≥2, spread≥13×) advanced. PROVED: the **k=2 slice** (exactly 2 large speeds ⟹ M≥1/14, via drop-max scaling V₂≥51 + bounded core Vmax≥63 + finite check ≤62 = 4865 sets, worst M=1/12) and **cluster-collapse Lemma A** (window W_K safe for the whole cluster, nonempty iff 13Vmin−Vmax>14Ks; closes single-gap clusters). CORRECTION (MISTAKE-076): the criterion C(S)=∃v:W(S\{v})>1/(7v) is SUFFICIENT but **NOT universal** — S\*={1,2,3,5,7,8,9,10,11,12,13,38,42} (covering S3, k=2) has C fail for all v yet M=2/23; so "prove C for all covering S" cannot close LRC(14), and a bounded-speed reduction is REFUTED (S3 infinite: AP family {t,2t,…,12t,V}). The residual is now UNIFIED and asymptotically TIGHT: the criterion margin's V₀→∞ limit-infimum is EXACTLY 1, lifted to a realized floor >1 by covering+primitivity discreteness — no compactness argument; closing it = a uniform positive density floor (Weyl/three-distance, route a) or multi-band CRT placement (route b). NEXT: the rigorous ρ*(Δ,P)≥c₀>0 equidistribution lemma. → THM-526, HYP-2581d.

**Private-obligation OPEN-Q-108 update (codex-2026-06-17-S5 / HYP-2579):** after THM-526, exact classification suggests the residual is not general covering pressure but parked-runner private q-debt.  In a `103`-row primitive q-covering scout, `94` rows were arc-width certified, `9` remained, and every residual had a parked runner uniquely covering some q-obligation.  The run separates q-covering from unit-residue completeness: `{1..12,182}` is q-covering with `M=14/183`, binding `(1,182)`, and missing unit residue `13`, closer to `1/14` than the `7/89` unit-complete champion `{1..11,13,84}`.  New proof target: prove private q-debt forces the THM-524 crossing index `j>=D/14`, or recurse/delete when no parked runner has private debt.
**Easy-dominates-hard OPEN-Q-108 update (kind-pasteur-2026-06-17-S2 / THM-525, HYP-2573..2576):** the covering case of LRC(14) is LOCALIZED onto OPEN-Q-108 by a non-circular reduction Q⟹P (Q=uniform meas(G_C)≥c): a covering 13-set is an easy 12-core C (LRC(12)-lonely) plus a runner w≡0 mod 14 parked in section 0; STEP-1 closes the non-covering case unconditionally, STEP-2/3 rewrite the covering case via meas(G_C) and the decoupling floor. The reduction reaches NOTHING strictly weaker than Q (the k≥3 coordinated-growing-speed regime is Q unchanged). NEW SHARPENINGS: (i) the conjectured extremal min meas(G_C)=7/858 confirmed over ~135k cores (drop-6 AP core); (ii) a plateau datum — driving 2 coordinated parked speeds →∞ does NOT send L→0 (L plateaus ≈0.0238, coordinated cores keep meas≥~0.095≫7/858), so the feared "meas(G_C)→0" failure mode did not materialize; (iii) a SECOND named sub-gap **G2** (the transversality estimate: w's danger comb, measure 1/7, concentrated near {a/w}, cannot CONTAIN all of G_C) — distinct from the uniform-measure GAP A, nonempty in every computed case, no general argument; (iv) the "perfect middle of section 0" is a constructive certificate device (survivor-sufficiency + meas(T_w)=6/7 PROVED), NOT the optimizer (which edge-binds w at 7/89); (v) the naive LRC(12)+Lipschitz lever is REFUTED (safe-arc half-width ~1/v_max). ~105k covering 13-sets, zero counterexamples. THE NEXT STEP that would CLOSE it: a bounded-speed reduction (Tao Thm 1.3 / MSS) making the finite covering-check a certificate up to v_max≤V₀; else a direct attack on G2.

**Section-checkoff OPEN-Q-108 update (codex-2026-06-17-S3 / HYP-2570):** the user's region-first idea becomes a finite Hall problem.  In the slowest-runner gauge, connect each runner to the fixed loop sections where it has a lonely witness.  Compactified runner-section graphs have perfect matchings in all tested primitive rows (`n=4,5,6`) and in LRC14 AP/Goddyn-Wong rows, while strict-open matchings fail because endpoint sections carry wall debt.  The new local target is a wall-switch lemma: every open Hall packet should either gain a section by crossing a boundary or descend to the dihedral endpoint-mouth / observer-source machinery.

**Dihedral OPEN-Q-108 update (codex-2026-06-17-S2 / HYP-2569):** the drop-6 extremal has an exact endpoint-orbit explanation.  Danger endpoints are local dihedral-clock events `(14k+/-1)/(14v)`, and a safe component from `aRk` to `bLl` has mouth length `(a*(14l-1)-b*(14k+1))/(14ab)`.  The drop-6 safe set is two reflected mouth orbits inside the omitted speed-6 moat: `2*(1/728)+2*(5/1848)=7/858`.  In the tooth coordinate `x=84(t-1/6)`, speeds `13,12,11` have residue defects `-1,0,1` and clipped cover `[-1,-8/13] union [-1/2,1/2] union [8/11,1]`, leaving normalized length `49/143`.  In the two-delete/one-replacement scan through `w<=180`, every missing-6 row is at least `7/858+1/980`, and rows that damage old hexagon mouths force larger new mouth mass.  New proof target: a scale-invariant dihedral mouth-exchange inequality.

**Latest OPEN-Q-108 update (codex-2026-06-17-S1 / HYP-2568):** exact 12-core sweeps support the sharper subtarget `meas(G_C) >= 7/858`, with equality at the AP drop-6 core `{1,2,3,4,5,7,8,9,10,11,12,13}`. No tested coordinated family beats it (`13026` two-drop/one-replacement cores through `w<=180`, `3000` random primitive cores, greedy swaps from the sporadic core). The conditional speed-load tournament is transitive, so future attacks should move from runner vertices to safe components, endpoint events, q-grid obligations, or proof-obligation packets.

**Bonferroni transfer-tax OPEN-Q-108 update (codex-2026-06-20-S62 / THM-558, HYP-2696):** the sector route now has an exact local ledger connecting the insertion DP to the true-wide Bonferroni gate: `Delta U4=mass(1->0)-mass(5->4)-4*mass(6->5)` for `U4=p0+p5+5p6`. Positive cap pressure is exactly unpaid one-missed closure; five/six-missed transitions are the tax. Incoming THM-557/HYP-2694 supplies the complementary coherent-block compression route; this THM-558 update is the local final-row ledger after that route.  This refines the HYP-2675/HYP-2693 branch target: finite low-state unpaid closures should route to AP/dyadic/cube-root/Ruzsa templates, while true-wide high-state closures should be paid by high-tail tax or bounded by Weyl/BV decorrelation. It does not close OPEN-Q-108, but it supplies a signed local accounting law for the wide-sector proof route.

**Small-q proof-lab OPEN-Q-108 update (codex-2026-06-22-S111 / HYP-2898):** applying the current LRC14 machinery to smaller even denominators `q=8,10,12` and back to `q=14` selects the stable proof carriers. Exact bounded banks show Bonferroni floors and p0-cap margins are positive throughout, consecutive/AP is the `nu`/dense-set extremizer, and AP difference-profile majorization has zero failures. But scalar additive energy is already non-monotone in smaller q (`12706` p0 inversions, `12139` dense-set inversions), and p0 itself can have non-AP bounded leaders that are still cap-safe. This sharpens OPEN-Q-108: do not chase scalar additive-energy monotonicity or literal AP-maximizes-p0. The viable cap/floor branch is Bonferroni + cap-safe p0 + AP-facing Fejer/difference-profile + labelled residual leak.
**Three-mode address-sheaf OPEN-Q-108 update (codex-2026-06-22-S114 / HYP-2902):** the Legendre correction is now routed into the proof DAG, refining HYP-2901's lcm denominator wall. Odd half-tiling is a full 3-set Venn with slots `A,B:N-1`, `C,D:N-2`, `E,F:N-3`, `G:N-4`; the `N-2` terms cancel only in scalar cardinality, not geometry. LRC14's even Eisenstein chart samples child sizes `13,12` while the pronic fold exposes apex coordinate `7`, so the proof must retain local address labels plus exact-period packets before scalarizing. The raw lcm family `S_X={1..11,13,lcm(2..X)}` proves no fixed finite denominator basis can close the problem (`q_min>X`), but its first witness is not generally `nextprime(X)`. OPEN-Q-108 split sharpened: finite Node 2 = AP-hull/three-gap/Fejer with sector labels; analytic Node 3 = exact-period/Weyl/L2 floor after divisor-killed denominators are removed.

**Status codes:** 🔴 CRITICAL (blocks main proof) | 🟡 IMPORTANT (needed for paper) | 🟢 INTERESTING (worth exploring)

## OPEN-Q-108 🔴 The uniform fattening lemma — the ONE lemma that completes the singular-series proof of LRC(14)
**Status:** OPEN — the isolated crux (kind-pasteur-2026-06-17-S1, THM-523, from the prove/disprove dialectic). By THM-523 (decoupling floor + single-perturbation inf=1/1260 + quantization THM-522), `inf_S L(S)>0 ⟹ C'(14) ⟹ LRC(14)` reduces to: **∃ c>0 with meas(G_C) ≥ c for EVERY 12-subset C of distinct positive integers**, where `G_C = [0,1)∖∪_{v∈C}D_v` is the gap-1/14 lonely set of `C` (`D_v` = `v`'s danger arcs). Equivalently: **the primitive tight locus at n=13 is FINITE** (conjecturally exactly `{AP {1..13}, Goddyn–Wong T5 {1..11,13,24}}`). KNOWN: the decoupling bound `L(C∪{w})≥(6/7)meas(G_C)−r/(7w)` handles speeds growing ONE AT A TIME (floor 1/143) and iterates for `k` large entries while the residual `(13−k)`-core keeps positive measure; the ONLY uncontrolled regime is `k≥3` arithmetically-coordinated growing speeds (the drop-6 family minimizes at the large `w=69`). THE LEVER (re-verified): PROVEN LRC(12) gives exactly one 12-subset of `{1..14}` tight at gap 1/13, zero at gap 1/14 — converting the crux from EXISTENCE (`meas(G_C)>0`?) to TRANSVERSALITY (does the isolated gap-1/13 maximizer FATTEN to a uniformly-positive gap-1/14 measure?). LITERATURE: the fixed-n tight-locus classification is "widely open" (Perarnau–Serra arXiv:2409.20160); Goddyn–Wong's infinite tight family needs n→∞; NO published bound controls the safe-MEASURE (only the gap κ(n)), so this is original. The bounded-speed reduction (Tao Thm 1.3/MSS) makes the compactness scaffold rigorous in principle. Entry: THM-523/522, OPEN-Q-097 (the complementary analytic Abel/Bedert route), 04-computation/tight_locus_*_kps.py.


## OPEN-Q-107 🟢 The Alcuin "+1" as a non-minor-closed correction: forbidden-subgraph set + the cover-internal-edge mechanism (general graphs)
**Status:** OPEN (mac-mini-2026-06-15-S6, THM-520, HYP-2553..2555). Complementary to OPEN-Q-106 (kind-pasteur, Ω-specific). General conflict-graph→tournament map T_G (i<j: arc i→j iff edge else j→i; forward arcs = edges). VERIFIED exact n≤6: independent set ↔ reverse-transitive run (τ(G)=n−largest reverse-transitive run of T_G), #3-cycles=#ordered induced P₃, #HamPaths(T_G) odd (Rédei shadow). Csorba–Hurkens–Woeginger: τ ≤ Alcuin ≤ τ+1 (the +1 decided by CHW Lemma 4.3 / Thm 3.1). HEADLINE (THM-520): τ is minor-monotone (⟹ {τ≤k} minor-closed ⟹ finite Robertson-Seymour obstruction set, the Kuratowski {K₅,K₃,₃} analogue), but **Alcuin is subgraph-monotone yet NOT minor-monotone — it fails ONLY under edge contraction** (smallest witness: contract an edge of K₂,₄, Alcuin 2→3; mechanism = contraction creates an edge INSIDE a minimum vertex cover). So {Alcuin≤k} is NOT minor-closed. QUESTIONS: (a) prove the contraction mechanism in general (Alcuin jumps iff contraction over-commits a min cover); (b) since Alcuin is subgraph-monotone it HAS a finite forbidden-SUBGRAPH obstruction set — compute it for small k (k=1: edgeless n≥1, K₁,₃, …); (c) the +1 has NO clean T_G strong-connectivity signature (HYP-2555 refuted "+1 ⟺ non-strong" and "G-Ham-cycle ⟺ ∃-order-strong"; the only clean strong-order fact is ∃-order-T_G-strong ⟺ G neither empty nor complete, n≤7) — find the right tournament invariant for the +1 or prove it is genuinely order/cover-combinatorial (not a tournament-iso invariant); (d) n≥7 confirmation. Entry: THM-520, 04-computation/alcuin_tournament_macmini_0615s6.py, CHW SIAM JDM 24(3) (JSTOR 41642576), Robertson-Seymour, Moon 1966.

## OPEN-Q-106 🟢 The forbidden-sub-tournament characterization of "Ω(T) planar" (Kuratowski/Robertson–Seymour for conflict graphs)
**Status:** OPEN (kind-pasteur-2026-06-16-S1, THM-519, HYP-2551). "`Ω(T)` planar" is a HEREDITARY tournament property (tournament-vertex-deletion = `Ω`-induced-subgraph). In the `α(Ω)=1` regime `Ω=K_m`, planar ⟺ `m≤4`, so the minimal obstruction is the `n=5, H=11` tournament (`Ω=K₅` = 5 pairwise-overlapping odd cycles). (a) Enumerate the minimal "Ω-non-planar" tournaments — is the `K_{3,3}`-driven obstruction (needs `α(Ω)≥2`) also minimal at some `n`? (b) Is the forbidden set FINITE or infinite? Tournaments are not WQO under sub-tournament order (Chudnovsky–Seymour antichains), so possibly infinite — but is it finite within bounded-cutwidth/pathwidth tournaments (where Chudnovsky–Seymour DO get WQO)? (c) CHW prove small-vs-large-boat is poly-decidable for planar graphs; so on planar-`Ω` tournaments, `Alcuin(Ω)` is easy though `τ(Ω)=#oddcycles−ν_odd` is hard — exploit this. Entry: THM-519, THM-517, Csorba–Hurkens–Woeginger SIAM JDM 24(3) 2010, `04-computation/alcuin_conflict_graph_kps.py`.

## OPEN-Q-105 🟢 A closed form / θ-product for the Burnside core kernel B[m,t], and the unified metagraph-enumerator family
**Status:** OPEN (kind-pasteur-2026-06-15-S7, THM-516, HYP-2538/2539; renumbered from OPEN-Q-103 — collision with codex-S12's prime-tail-ladder OPEN-Q-103, a COMPLEMENTARY view of the same kernel). The 1-tail peeling isolates A000568's difficulty in the `n`-independent core kernel `B[m,t]=Σ_{μ:|μ|=m,ℓ=t,odd parts≥3}2^{e(μ)}/z_μ`, with `e(μ)=C(t,2)+½Σ_{d odd≥3}φ(d)M_d²` (positive-definite GCD/Euler-φ quadratic form → theta-function exponent; VERIFIED 1113 cores). (a) Does `B[m,t]` have a closed product/θ generating function over part-sizes, diagonalizing the φ(d)M_d² cross-coupling (HYP-2538)? The add-a-part recurrence `Δe=(p−1)/2+Σ_{d|p}φ(d)M'_d` closes only on the divisor-profile state — that IS the residual difficulty. (b) Do `G_n(x)` (graphs, A000088), `SC(n)` (base-4), `E_n` (A002854 even graphs) share the same core-kernel compression with their per-part rule (HYP-2539)? (c) Is the `P_n(x)=Σ_{odd-cycle σ}x^{#edge-orbits}` coefficient triangle (rows→A000246) a new OEIS sequence? Entry: THM-516, THM-505 (cores = A000009−3 OCF non-spectral family), codex-S12's THM-514 (the tail-ladder complement), `04-computation/burnside_core_kernel_phi_reframe_kps.py`.

---

## OPEN-Q-101 - Upgrade moment-shadow witnesses to genuine compatibility inequalities

**Status:** OPEN (monad-explorer-2026-06-15-S10, HYP-2530; extends OPEN-Q-099, HYP-2458, HYP-2457, HYP-2529). The new computation shows the flagship baby-Hodge hole `(c3,c5)=(8,10)` at `n=6` is already inside the simplest spectral/moment shadow: `(8,10)=(1/3)*(8,8)+(2/3)*(8,11)` in the same `c3=8` fiber. It also shows the odd Faulhaber moments are a positive Stieltjes sequence `S_{2r+1}(n)=sum_i i*(i^2)^r`, so their Hankel positivity does not explain why exact towers stop after `p=2`. The open problem is to write the missing retained variable explicitly as a compatibility inequality: on the tournament side, a flag-algebra / PSD / conflict-packet statement that cuts `c5=10`; on the Faulhaber side, a packet or integrality field playing the role of `D=alpha_2`; and on the repunit side, an atom-supply obstruction beyond scalar length. Files: HYP-2530, T822, `04-computation/baby_hodge_compatibility_wall_monad_s10.py`, reflection `the-moment-shadow-and-the-compatibility-wall-monad-s10.md`.

## OPEN-Q-103 - Does the A000568 prime-tail ladder close after finitely many divisibility statistics?

**Status:** OPEN (codex-2026-06-15-S12, THM-514). The `1`-tail is completely soluble:
`a(n)=sum_{m,t} B[m,t] 2^(C(n-m,2)+(n-m)t)/(n-m)!`, collapsing the `n=100` outer odd-part
sum from `444793` partitions to `834` active `(m,t)` states. The next rung is also exact:
peeling `3`-cycles only needs one extra statistic `c_3`, and the `3`-free kernel at mass
`100` uses `2049` active `(m,t,c_3)` states against `7551` residual odd partitions with
parts at least `5`. More generally, peeling a prime `p` from the current core needs only
`u=ell(nu)` and `c_p(nu)=# {parts divisible by p}` in the residual partition. The open
problem is whether this ladder stabilizes in a finite or controlled state family, or
whether iterating over odd primes inevitably reconstructs the full divisor-incidence DP.
Concrete targets: characterize the minimal sufficient statistic set after primes
`3,5,7`; prove or refute a finite-statistic closure theorem; and quantify the state growth
of the prime-tail ladder versus raw odd-partition growth.

## OPEN-Q-096 🟢 — The other faces of the master cycle-packing polynomial Φ

**Source:** monad-explorer-2026-06-15-S5, HYP-2514, reflection `the-master-cycle-packing-polynomial`, THM-505.

The spectrum (Sachs Coefficient Theorem, all-length **signed** `y_k=−x^{−k}`) and the OCF
`H` (odd-only **unsigned** fugacity-2 `y_k=2[k odd]`) are two evaluations of one master
disjoint-cycle-packing polynomial `Φ(T;{y_k}) = Σ_{packings}∏ y_{|C|}`. Open:

1. **The ALL face is the PERMANENTAL POLYNOMIAL — RESOLVED (monad-S6, THM-506, HYP-2515).**
   The all-length unsigned face of Φ, graded by **vertices**, is `per(xI+A) = Σ_m e_m^unsigned
   x^{n−m}` (`e_m^unsigned`=#packings on m vtx) — the permanental polynomial, the **unsigned
   twin of the char poly** `det(xI−A)`; the two differ ONLY by the cycle-parity sign
   `(−1)^#cyc→+1` (the det/per dichotomy; det spectral & poly-time, per non-spectral &
   #P-hard). Non-spectral first at n=6 (same wall), splits the same 47 cospectral classes as
   H at n=7 but strictly **finer** (recovers (c6,c7,…) vs H's one functional); **(char,perm)
   determines H iff n≤7**, breaking at n=8 via the D44↔D35 trade (within-class rank
   3→4). **EVEN FACE RESOLVED (monad-S7, HYP-2517):** the even face's SIGNED form IS a
   clean matrix function — the skew char poly `det(xI−S)`, `S=A−Aᵀ`, `=∏(x²+μ_j²)` with
   coeffs `Σ_W Pf(S[W])²` (Coates: odd cycles cancel under reversal) — but it is
   **SPECTRAL** (a function of charA; VERIFIED n≤6 exh, incl. the cospectral-different-H
   n=6 pairs; matches the known complement=converse spectral-DS equivalence). So the
   Pfaffian route recovers only the spectral shadow; the non-spectral even content
   `I(Ω_even,·)` (splits 3@n6,46@n7) has NO determinantal home, like H. **The det/per
   (Valiant) dichotomy = the spectral/non-spectral boundary, face by face**; the ODD face
   has no determinantal object at all (irreducibly non-spectral). **WALK-COUNT LINCHPIN
   NOW PROVED (monad-S8, THM-507/HYP-2518):** the clean general proof that walk counts
   `w_k=1ᵀAᵏ1` are spectral — exact closed form `1ᵀadj(xI−A)1=(−1)ⁿcharA(−x−1)−charA(x)`,
   equivalently `F(x)=∏(x+1+λᵢ)/∏(x−λᵢ)−1`, via the matrix-determinant lemma + the
   tournament identity `A−J=−(Aᵀ+I)` (the all-ones perturbation collapses to a
   transpose-shift with FORCED eigenvalues `{−1−λᵢ}`, no angle dependence; this is exactly
   why tournaments escape the cospectral walk obstruction `C₄⊔K₁` vs `K_{1,4}`). So the
   whole A-affine pencil determinant is now PROVABLY spectral, not just verified ≤n7;
   `w_k=C(n,k+1)+spectral cycle corrections` (`w_2=C(n,3)+2c₃`, `w_3=C(n,4)+(2n−3)c₃`);
   reciprocity `(1+F(x))(1+F(−x−1))=1` centred at `−1/2` (= fixed point of complementation
   on the spectral axis, same `−1/2` as THM-055/059/080). STILL OPEN: the
   permanental ROOTS as an invariant; the even/all dimension growth law; the general-n
   carrier deficit of (char,perm); the POINTED walk data `1ₐᵀ(xI−A)⁻¹1_b` / `M[a,b]` —
   where exactly does the spectral closed form break as we de-contract (handoff, this session).
2. **The signed odd face is "more spectral" — RESOLVED at n=7.** `sgn_odd =
   Σ_{odd packings}(−1)^{#cyc} = I(Ω,−1) = −χ̃(Ind Ω)`, the reduced **Euler characteristic of
   the odd-cycle packing complex**. At n=7 it equals `(1+e₃+e₅+e₆) + (c₆−c₇)` (verified
   3000/3000): non-spectral content is the **1-D** direction `c₆−c₇`, a projection of `H`'s
   2-D `(c₆,c₇)`. It splits 16 cospectral classes (iff `Δ(c₆−c₇)≠0`) vs `H`'s 47 (iff
   `Δ(2c₆+c₇)≠0`); `47 = 16 + 31`, the 31 gap = classes where `c₆~c₇` **covary** and the
   alternating sign `x=−1 ⊥ (1,1)` cancels them. **The fugacity `x` is a dial selecting a
   linear functional of the non-spectral level-sums `α_j`.** OPEN: the general-n analogue
   (the Euler-char direction vs the `H` direction in the `⌊n/3⌋`-D `α_j` space); does some
   `x` minimize the split count (the "most spectral" fugacity)?
3. **General-n Sachs-basis skeleton.** Prove `S = 1 − 2e₃ − 2e₅ + 4·Σ_{m even ≥6}e_m`
   (verified n≤9) for all `n ≤ 11`, and derive the `n=12` refinement when the `(3,3,3,3)`
   quadruple level switches on.

---

## OPEN-Q-093 🟡
**Can corrected trace vectors compute higher tournament cycle structure and compress H beyond n=6?**

HYP-2498 proves/validates the first trace correction boundary: `c_k=tr(A^k)/k` for `k=3,4,5`, while `tr(A^6)=6*c6+3*c3+6*p33_meet`, with `p33_meet` counting intersecting directed-triangle pairs. Exhaustive labelled `n=6` data shows the corrected low cycle vector `(c3,c4,c5,c6)` determines `H`, even though `score+c5+c6` still has a mixed `H` bucket. Build the support-type correction engine for `tr(A^7)` and beyond; test whether corrected trace vectors continue to determine or sharply compress `H` at `n=7`.

**Source:** HYP-2498, codex-2026-06-13.

**PARTIAL ANSWER (monad-explorer-2026-06-13, THM-500):** the `tr(A^7)` correction is
`tr(A^7) = 7*(c7 + TQ)`, i.e. `c7 = tr(A^7)/7 - TQ`, where `TQ` = #(directed-triangle,
directed-4-cycle) pairs with overlapping support (exact, 600/600). Odd analog of
`p33_meet`. BUT it does NOT compress `H` further: `TQ` (hence `c7`, `alpha_1`, `H`) is
the FIRST non-spectral quantity inside `alpha_1` — cospectral `n=7` tournaments realize
`c7 ∈ {4,5,10}` at identical `tr(A^k)` (THM-500). So at `n=7` the *uncorrected* trace
vector stops determining `H`; reconstructing `H` needs the overlap counts `TQ`/`DTP`
themselves, which are not power sums. The corrected-trace engine works, but the
corrections it requires are exactly the non-spectral support-geometry data.

**FULL ANSWER (monad-explorer-2026-06-14, THM-502 — the closed-walk census ladder):**
the correction engine has a *single generating principle*. `tr(A^k)` counts rooted
closed k-walks; loop-erasing gives a connected multiset of overlapping simple cycles
(parts ≥3) partitioning k, and each configuration `C` contributes `k/period(C)`. This
yields the **complete explicit ladder** through k=8 (the last k before triple configs):
`tr A^6=6c6+3c3+6 p33`, `tr A^7=7c7+7 TQ`, **`tr A^8=8c8+4c4+8 Q44+8 TF`** (Q44 =
overlapping 4-cycle pairs, TF = overlapping (triangle,5-cycle) pairs — NEW). The
distinct-pair coefficient is uniformly `k`; a doubled (k/2)-cycle contributes `k/2`.
**Conservation corollary** (the engine's structural content): within a cospectral
class, `c6+p33`, `c7+TQ`, `c8+Q44+TF` are spectral constants, so the simple top-cycle
count trades 1-for-1 against the overlap count — the exact reason corrected trace
vectors do NOT compress `H` past n=6: the corrections ARE the non-spectral
support-geometry. Confirms (via the exhaustive n≤6 spectral-horizon table) that `c6`
is non-spectral *from its onset* n=6, so `alpha_1` (THM-500) is the unique delayed
break. k=9 opens the first **triple** term (3+3+3); coefficient law verified, the
distinct-triple enumeration by overlap topology is the remaining open frontier.

**INVERTED + DIMENSION (monad-explorer-2026-06-15, THM-505):** the census ladder
*reconstructs `H` exactly*. Substituting the defects (`p33=W6−c6`, `TF=W8−c8−Q44`) into
the OCF `H=I(Ω,2)=Σ2^k α_k` gives a clean **spectral skeleton + carrier** split:
n=7 `H = [1+2c3+2c5+4C(c3,2)−4W6] + 4c6 + 2c7` (PROVED); n=8 adds `+4c8+4Q44`
(equiv. minimal-defect `+4c6+2c7−4TF`). The fugacity `x=2` sets the weights (`2^level`).
So the corrected-trace engine does NOT compress `H`, but it *coordinatizes* its
non-spectral content exactly. **Dimension answer:** the number of independent
non-spectral DOF of `H` is `n−5` (=0,1,2,3 for n=5..8); the carriers are the simple-cycle
counts `{c6,...,c_n}`, and every overlap defect (incl. the even `Q44`) is a spectral
function of them (n=8 probe: (c6,c7) insufficient/157 split buckets, (c6,c7,c8) determines
Q44/0 free).

**n=9 RESOLVED — DIMENSION BREAKS, LINEARITY NEGATIVE (monad-explorer-2026-06-15,
THM-505 extended):** n=9 closed form PROVED & verified 45000/45000:
`H = [skeleton] + 2c7 + 2c9 + 4c6 + 4c8 + 4Q44 + 8·T333` (the triple level α₃=T333 turns
on with weight 8=2³; full fugacity form `I(Ω,x)=SKEL(x)+(c7+c9)x+(c6+c8+Q44)x²+T333x³`).
(1) **dim ≠ n−5 at n=9 — it is EXACTLY 6:** nested chain (130000 samples)
`sig→+(c6,c7,c8)→+c9→+Q44→+T333` splits `14804→482→24→1→0`. So `(c6,c7,c8,c9)` does NOT
determine `H` (24 witnesses), and neither does `(c6,c7,c8,c9,Q44)` (1 residual split). The
minimal carrier set is the full `{c6,c7,c8,c9,Q44,T333}` (capped at 6 by the closed form),
so `dim_nonspec(H) = 6 > n−5 = 4`. BOTH overlap configs Q44 and T333 are INDEPENDENT
carriers, not spectral shadows — the dimension JUMPS 3→6 at n=9 (n=8 was the last size where
they coincided; chain8 at 60000 confirms (c6,c7,c8) determines H there). (2) **Linearity NEGATIVE:** `H` is universal-linear in the full
carrier set (incl. overlaps) but NOT a bounded-degree polynomial in the simple cycles
alone past n=7 (n=8 within-class linear & quadratic fits inexact). Three-stage degradation:
linear (n≤7) → non-polynomial-functional (n=8) → independent-correlations (n=9). The
non-spectral content of H is a correlation TOWER, not a flat vector. FILES: THM-505,
04-computation/ocf_nonspectral_n9_monad.py, reflection
the-overlaps-stop-being-shadows-the-correlation-tower.

**GROWTH LAW RESOLVED — A PARTITION FUNCTION (monad-explorer-2026-06-15-S3):** the dimension
growth law DOES have a closed form. In the basis-independent OCF *packing* basis (expand
`H = Σ_λ 2^{|λ|} N_λ` over length-multisets `λ` of disjoint odd-cycle packings, parts odd ≥3),
`dim_nonspec(H)(n) = #{partitions of s≤n into odd parts ≥3} − 3` = `Σ_{s≤n}[x^s]Π_{k odd≥3}1/(1−x^k) − 3`
= **1,2,3,5,7,9,12,15,19** for n=6..14 (increment = p_{odd≥3}(n)). VERIFIED by RANK of the
within-class carrier-delta matrix: dim = **3,5,7,9** at n=8,9,10,11 (every OCF carrier
independent, H in span, OCF holds 6000/6000 at n=10 and 5000/5000 at n=11, 704 split
cospectral classes). The new n=10 carriers are exactly the (3,7) and (5,5) disjoint pairs
`D37,D55`; n=11 adds `c11` and the new triple `T335` = {3,3,5}. CORRECTION: the n=9 dim
is intrinsically **5, not 6** — the trace-basis "6" over-counted because `c8` and `Q44` enter
`H` only via their sum `D35`. `dim ≤ #{λ}−3` PROVED; equality (no `N_λ` spectrally pinned)
VERIFIED n≤11, CONJECTURE general. FILES: 04-computation/ocf_nonspectral_n10_monad.py,
05-knowledge/results/{ocf_nonspectral_n10_n11_monad.out, ocf_nonspectral_n11_monad.out}, reflection
the-non-spectral-dimension-of-H-is-a-partition-function. NEW open: (1) PROVE no `N_λ` is
spectrally pinned (upgrades the law to a theorem); (2) where does the law first deviate, if
ever — does a tightness-pinning kill a carrier at large n (the rank could drop below the
partition count)? (3) the OEIS id of `1,2,3,5,7,9,12,15,19` / partial sums of partitions into
odd parts ≥3.

**OEIS RESOLVED + TWO-DIMENSIONS CORRECTION (monad-explorer-2026-06-15-S4):** (3) **the
sequence is A000009.** One-line GF identity: `Σ_{s≤n}[x^s]Π_{k odd≥3}1/(1−x^k) =
[x^n]Π_{k odd≥1}1/(1−x^k) = q(n)` (the cumulative `1/(1−x)` IS the missing odd part `k=1`),
and `q(n)` = #partitions of `n` into odd `=` distinct parts `=` **A000009**`(n)`. So
`dim(packing) = A000009(n)−3`, asymptotically `~ exp(π√(n/3))/(4·3^{1/4}n^{3/4})`
(super-polynomial). **CORRECTION to the headline:** `A000009(n)−3` is the non-spectral dim of
the **packing-count vector** `(N_λ)`, NOT of `H`. `H = I(Ω,2) = 1+Σ_{j≤⌊n/3⌋}2^j α_j` depends
only on the **level-sums** `α_j=Σ_{|λ|=j}N_λ` (it never sees the length-split of a level), and
`α_j=0` for `j>⌊n/3⌋`. So **`dim_func(H)(n) ≤ ⌊n/3⌋` (LINEAR), PROVED**, `= ⌊n/3⌋` for n≥7;
verified n=8 (level-sum rank 2 vs carrier rank 3, H in level-sum span). The fugacity-2
evaluation compresses `exp(√n)→n/3`. NEW open: (4) prove the levels `α_j` are non-spectrally
independent (⟹ `dim_func(H)=⌊n/3⌋` exactly). FILES: 04-computation/ocf_two_dimensions_monad.py,
05-knowledge/results/ocf_two_dimensions_n89_monad.out, reflection `H-reads-only-the-level-grading`.

## OPEN-Q-092 🟡
**Can Pollock's tetrahedral no-long-pair lemma be proved with dyadic carry-pair ledgers instead of single residues?**

HYP-2497 shows that the Sierpinski/Waring analogy is not a plain mod-2 obstruction: scanning `k < 4*2^e+16`, the single tetrahedral residues `{Te_k mod 2^e}` are all of `Z/2^e Z` for every `1<=e<=12`. But after HYP-2491's lift to defect pairs `r,r+tri(k) in D_4`, the tail pair-residue universe compresses sharply: for `k>=100`, observed pair classes stabilize at `168` by `2^8` while the possible pair classes grow as `4^e`, yielding a transitive dyadic compression tournament `12>11>...>3`. Prove the 2-adic surjectivity lemma, then use pair/carry/convolution constraints to rule out triangular defect self-correlations for `k>825`.

**Source:** HYP-2497, codex-2026-06-13.

## OPEN-Q-091 🟡
**Can Pollock's tetrahedral conjecture be proved by forbidding long triangular self-correlations in the four-defect set?**

HYP-2491 reframes Pollock's five-tetrahedral conjecture around `D_4`, the integers not representable by at most four tetrahedral numbers. For `n=Te_k+r`, a one-back descent works unless both `r` and `r+tri(k)` lie in `D_4`. The computation rediscovers the known `241` four-defects through `10^6`, largest `343867`, and the last triangular defect-pair separation among them is `3142 -> 343867 = 3142 + tri(825)`. Prove either the strong tail `D_4 subset [1,343867]` or the weaker no-pair lemma for all `k>825`; pair this with the width-3 finite shell certificate.

**Source:** HYP-2491, codex-2026-06-13.

## OPEN-Q-089 🔴
**Can LRC14 long blocking-height rows be split into peelable-carrier or balanced-cover-congruence cases?**

HYP-2481 shows that raw cumulative speed dominance grows with blocking height, but normalized dominance falls and the speed tournament becomes transitive in named hard packets. Prove a dichotomy: either some cumulative/private cover carrier can be peeled, transported, or converted into a Bprime/owner opening, or the lack of such a carrier forces balanced-cover congruences that enter the Q31/band-2 ramified portal of HYP-2471/HYP-2480. Immediate computational subtask: add leave-one-out support-criticality to `lrc14_blocking_height_dominance_codex.py` and test the five one-stranger evaders plus the two HYP-2470/HYP-2471 exception shapes.

**Source:** HYP-2481, codex-2026-06-13.

## OPEN-Q-090 🟡
**Can the source-deleted A000568 fingerprint be paired with Q31/band-cover ledgers to force the LRC14 portal?**

HYP-2486 isolates the clean A000568 layer in LRC: a threshold source-lift has the observer as source exactly at LRC-good states, and deleting that source leaves an ordinary moving-runner tournament class. Raw runner classes mix good and bad states, but the rooted source fiber is pure in the exact audits. For LRC14, attach this deleted-class fingerprint to each shell/blocking state together with Q27/Q31 obligations, divisor fiber, owner/Bprime debt, and HYP-2481 support loads. Prove that a long blocked-band walk either reaches a source-cone deleted class or its avoidance forces balanced-cover congruences, hence the HYP-2471/HYP-2480 Q31/7-ideal/13-clock portal.

**Source:** HYP-2486, codex-2026-06-13.

## OPEN-Q-001 -- RESOLVED
**The n=5 mystery: why does the per-path identity hold despite 5-cycles?**

**RESOLVED by THM-008:** The per-path identity holds trivially for n<=5 because mu(C) = 1 for ALL odd cycles C through v. For 3-cycles, the complement V\{v,a,b} has at most 2 vertices, which cannot form an odd cycle. For 5-cycles, C\{v} exhausts all of T-v, leaving 0 available vertices. The identity reduces to #TypeII = #TypeII. There is no "delicate balance" -- the identity is vacuous at n<=5.

**Additional detail (opus-S2):** More generally, mu(C) = 1 whenever cycle length L >= n-2 (THM-008 mu triviality bound). At n=6, mu(3-cycle) is in {1, 3}: mu=1 (76.7%) when 3 available vertices form transitive subtournament, mu=3 (23.3%) when cyclic.
## ~~OPEN-Q-001~~ RESOLVED
**The n=5 mystery: why does the per-path identity hold despite 5-cycles?**

**Resolved by:** opus-2026-03-05-S1 (THM-008)

**Answer:** At n=5, mu(C) = 1 for ALL cycles C through v (both 3-cycles and 5-cycles). This is because C\{v} leaves too few available vertices in T-v for any odd cycles to exist in the restricted conflict graph. Specifically, |Available| = n - L, and odd cycles need >= 3 vertices, so mu = 1 whenever L >= n-2. At n=5, both L=3 (available=2) and L=5 (available=0) satisfy this. The per-path identity holds not because of a deep structural coincidence but because the mu weights are trivially 1, reducing Claim A to a simple cycle-counting identity. See THM-008 for the full proof.

---

## OPEN-Q-002 -- RESOLVED
**Prove Claim A: H(T) − H(T−v) = 2Σ_{C∋v} μ(C)**

**RESOLVED by kind-pasteur-2026-03-05-S12:** Claim A is PROVED for all n.

**Proof:** OCF (H(T) = I(Omega(T), 2)) is proved by Grinberg & Stanley
(arXiv:2307.05569, 2023; arXiv:2412.10572, 2024, Corollary 20).
Their formula: ham(D̄) = Σ_{σ ∈ S(D), all cycles odd} 2^{ψ(σ)}.
For tournaments, D̄ = D^op (converse) and ham(D^op) = ham(D) by path reversal.
The RHS = I(Omega(D), 2) since independent sets in Omega(D) biject with
collections of vertex-disjoint odd directed cycles.
Therefore H(T) = I(Omega(T), 2). Combined with Claim B (THM-003, proved),
this gives Claim A. See CONJ-001, THM-002.

**Prior verification record:** n≤8 exhaustive (THM-015), n≤10 random sampling, all consistent.

---

## OPEN-Q-003 -- RESOLVED
**Characterize when the per-path identity holds at n=6**

**RESOLVED by THM-009:** The per-path identity fails for path P' iff some Type-II position (a,b) in P' has mu(v,a,b) > 1, which happens iff the 3 vertices V\{v,a,b} form a directed 3-cycle in T-v. This is a perfect binary separation: mu>1 at any TypeII position => always fails; all mu=1 => always holds.

---

## OPEN-Q-004 🟢
**Find a correct per-path formula for all n**

The 3-cycle-only formula (per-path identity) fails for n≥6. The natural generalization summing over all odd cycles overcounts. The maximal-embedding-only formula also fails. Is there a formula of the form Σ_{cycles C, (non-v consecutive in P')} f(C, P') = (inshat−1)/2 that works for all n?

**Note:** Since OCF/Claim A is now proved (Grinberg-Stanley), this is no longer blocking any main result. Downgraded from 🟡 to 🟢.

---

## OPEN-Q-005 -- RESOLVED
**Combinatorial proof of the C(L-2, 2k-1) distribution (THM-007)**

**RESOLVED (INV-029, opus-S5):** Bijective proof found. See INV-029 in INVESTIGATION-BACKLOG.md.

---

## OPEN-Q-006 🟢
**Asymptotic formula for Σ_C μ(C)**

The average Type-II contribution per L-cycle window is (L-4)/4, growing linearly with L. Does this yield an asymptotic formula for Σ_C μ(C) as a function of the cycle-length distribution of T? What happens for random tournaments as n→∞?

---

## OPEN-Q-007 🟡
**Full proof of Fix(σ) = 2^{m²} for self-evacuating SYT**

Verified for n=5 and n=7 (m=2 and m=3 respectively, giving 4 and 512 self-evacuating SYT). Full proof is conditional on a precise classical reference not yet pinned down. The identification with TSSCPPs may provide the reference.

---

## OPEN-Q-008 🟢 — PARTIALLY RESOLVED
**2-adic tower: what is the 2-adic valuation of H(T)?**

**PARTIALLY RESOLVED (opus-2026-03-05-S13):** v_2(H(T)) = 0 for ALL tournaments (this IS Redei's theorem — H(T) is always odd). Verified exhaustively at n≤6 and sampled at n=7 (5000 tournaments).

The mod-4 structure: H(T) ≡ 1 + 2·alpha_1 (mod 4) via OCF, where alpha_1 = #odd cycles in Omega(T). At n=3,4 this equals 1+2·c_3 (mod 4), but at n≥5 the relationship breaks because 5-cycles contribute.

**Reformulated question:** What is the distribution of H(T) mod 2^k for k≥2? Computations show it approaches uniform on odd residues as n grows. The OCF gives H mod 4 via alpha_1 parity, H mod 8 via alpha_1 and alpha_2, etc.

---

## OPEN-Q-009 -- RESOLVED
**Prove arc-flip identity: E(T) = H(T) - I(Omega(T), 2) is invariant under arc flips**

**RESOLVED by kind-pasteur-2026-03-05-S12:** E(T) = 0 for ALL tournaments (not just invariant).
OCF (H(T) = I(Omega(T), 2)) is proved by Grinberg-Stanley (arXiv:2412.10572, Corollary 20).
See THM-002, CONJ-001 for the complete proof chain.

**Historical work (preserved for reference):**

The project independently discovered and partially proved OCF via multiple routes:
- **THM-015**: Proved delta_H = delta_I as polynomial identity at n <= 8 (exhaustive)
- **THM-016/017**: Proved the even-odd split for all n (inductive proof via Claim B path identity)
- **THM-018**: Proved coefficient identity alpha_w^H = alpha_w^I symbolically at n <= 8
- **MISTAKE-008**: Correctly identified that even-odd split is necessary but NOT sufficient for OCF

The even-odd split (THM-016/017) was the strongest general-n result obtained internally.
The gap between even-odd split and full OCF is now bridged by the Grinberg-Stanley proof.

**Key structural facts discovered along the way:**
- All affected cycles contain {i,j} (complement unchanged by flip)
- At most one affected cycle in any independent set (A-clique)
- The swap involution (THM-014) gives adj(i,j)-adj'(j,i) = #U_T - #U_T'
- Even-odd split: delta decomposes equally between even-S and odd-S terms
- The s-coefficient identity (THM-018) reduces OCF to a per-vertex polynomial identity

See PROP-001, THM-013, THM-014, THM-015, THM-016, THM-017, THM-018.

---

## OPEN-Q-014 -- RESOLVED (DISPROVED)
**Prove Omega(T) is always perfect (and possibly claw-free)**

**DISPROVED by opus-2026-03-05-S7:**
- **Perfectness FAILS at n=8.** 53.8% of random n=8 tournaments have a C5 (5-hole) in the
  3-cycle conflict subgraph of Omega(T). Explicit counterexample constructed.
- **Claw-freeness TRIVIALLY holds at n<=8** (vertex counting: 3 pairwise disjoint odd cycles
  + 1 touching all three requires >= 9 vertices). FAILS at n=9 (90% of random tournaments).
- **Perfectness holds for n<=7** (0 failures in 1000 random trials).
- **OCF still holds** despite Omega(T) being imperfect (proved by Grinberg-Stanley).

The all-real-roots property of I(Omega(T), x) and log-concavity still hold empirically
at n<=6. Whether they hold at n>=8 (where Omega is imperfect) needs separate investigation.

See THM-019 (corrected), `04-computation/omega_c5_test.py`, `04-computation/omega_claw_fast.py`.

**Source:** opus-2026-03-05-S7 (disproof)

---

## OPEN-Q-010 -- RESOLVED (NEGATIVE)
**Per-path formula including 3-cycles AND 5-cycles at n=7**

At n=7, mu(5-cycle) = 1 always (V\{v + 4 cycle vertices} has 2 vertices, no odd cycles). So 5-cycle contributions are "trivially weighted" just like 3-cycles at n<=5. A per-path formula summing over both 3-cycle and 5-cycle embeddings (each with their mu weights) might work at n=7. Test computationally.

**Status (kind-pasteur-2026-03-05-S3):** NEGATIVE RESULT. The per-path formula does NOT simplify at n=7. The algebraic identity (inshat-1)/2 = #{TypeII} = #{3-cycle embeddings} is trivially satisfied, but this just restates THM-004+005 -- it does not encode 5-cycle information. Computing the actual A, B, D quantities (see test_n7_ABD.py) shows A=/=D in general. A=/=D means: total TypeII count (A) does not equal total odd-cycle mu sum (D). The 5-cycles contribute non-trivially even when mu=1. See T027 and OPEN-Q-011.

**Source:** FINAL_FINDINGS.md, Q3; kind-pasteur-2026-03-05-S3

---

## OPEN-Q-011 -- RESOLVED (statistical artifact, not structural)
**Near-cancellation of two error effects at n=6**

**Resolved by:** opus-2026-03-05-S2, confirmed by kind-pasteur at n=7

**Answer:** The near-cancellation is a statistical observation, NOT an exact identity. Computational verification (3000 pairs at n=6, opus-S2) shows:
- A = D exactly for only 836/3000 (28%) of pairs
- A - D ranges from -12 to +9 (mean ≈ 0)
- Mean(A-B) ≈ -Mean(B-D) is approximate, not exact

The decomposition A-B = -(B-D) does NOT hold pair-by-pair. The two effects cancel on average but not structurally. This is NOT a viable proof strategy for Claim A.

**Status (kind-pasteur-2026-03-05-S3):** PARTIAL ANSWER. At n=7, tested 1050 (T,v) pairs: mean A-D = 0.097 (near zero), but NOT zero in general (range -39 to 26). Mean A-B = -73.78, mean B-D = +73.88 (near-cancellation on average). The near-cancellation is STATISTICAL, not algebraic. The per-pair |A-D|<=1 holds only 13.1% of time. The decomposition Claim A = (A=B) + (B=D) does NOT yield two tractable sub-identities. The near-cancellation at n=6 was likely a low-n coincidence.

**Source:** FINAL_FINDINGS.md; kind-pasteur-2026-03-05-S3

---

## OPEN-Q-012 🟢
**Tower hypothesis: L-cycle corrections from (L+2)-cycles**

At n=2k, the first cycle whose mu can exceed 1 has length 2k-1. The "excess" mu from shorter cycles may be exactly compensated by contributions from cycles 2 vertices longer. Is there a recursive structure where L-cycle corrections are expressed in terms of (L+2)-cycle contributions, creating a tower that sums to Claim A?

**Source:** FINAL_FINDINGS.md, Q5

---

## OPEN-Q-013 🟡
**Correct formula for H(T_p) for Paley primes p ≡ 3 (mod 4)**

Both conjectures are FALSE for p=11:
- Original conjecture H(T_p) = p * 3^((p-1)(p-3)/8) gives 649,539 for p=11, not divisible by 55.
- Revised conjecture H(T_p) = |Aut(T_p)| * 3^((p-3)/2) gives 4455 for p=11 (off by factor 21.4).

**Known values (all confirmed by direct computation):**
- p=3: H=3, |Aut|=3, H/|Aut|=1
- p=7: H=189, |Aut|=21, H/|Aut|=9=3^2
- p=11: H=95095, |Aut|=55, H/|Aut|=1729=7*13*19 (Hardy-Ramanujan taxicab number)
- p=19: H=1,172,695,746,915, |Aut|=171, H/|Aut|=6,857,869,865 (computed opus-S5/S10)
- p=23: H=15,760,206,976,379,349, |Aut|=253, H/|Aut|=62,293,308,207,033=3*167*4567*27225299 (computed opus-S10, factored kind-pasteur-S1)

**Sequence H/|Aut|:** 1, 9, 1729, 6857869865, 62293308207033 — no obvious pattern. 3^k pattern breaks catastrophically at p=11. Factorizations are erratic. |Aut(T_p)| = p*(p-1)/2 confirmed for all p (affine QR group).

**ADDENDUM (monad-explorer-2026-06-07, HYP-2306): the "modular significance" of 1729 is REFUTED — and the erraticness is now EXPLAINED.** `the-tessellation.md` (Layer 6, opus-S131) read `1729 = r(11) = 7·13·19 = j(i)+1` as a *modular* fact (completely split in Q(√−3); appeared at the first genus-1 Paley prime). The sharpest test is **p=19, the next Paley prime, which is ALSO genus 1** (`X_0(11)=X_0(19)=`genus 1, `X_0(23)=`genus 2) — *not* the p=23 the reflection guessed. The structure does NOT persist: `r(19)=5·7·11·23·774463` has 5,11,23 INERT (≡2 mod3) and a large prime; `r(23)=3·167·4567·27225299` has 167,27225299 inert. **Mechanism:** 1729 is clean only because `H(T_11)=5·7·11·13·19` is an unusually smooth product of 5 small primes; `H(T_19),H(T_23)` carry large primes (774463; 27225299), so r(p) can never again be smooth / completely-split / near a j-value. The factorizations are erratic *because smoothness is a small-p artifact.* The genuine regularity of the sequence is ANALYTIC, not arithmetic: `H(T_p)·2^{p−1}/p! → e` (see ratio line below; the real law). Both H(T_19) and H(T_23) were INDEPENDENTLY re-verified here by a validated int64 Held-Karp counter (`04-computation/paley_H23_monad.py`). This severs the cross-lane "1729 spine" (tournament ratio ↔ S5 Moser-ladder record rung ↔ Klein's 1728): only the Moser-ladder 1729 is structural.

**Ratio H(P(p))/(p!/2^{p-1}):** 2.000, 2.400, 2.440, 2.527, 2.557 for p=3,7,11,19,23 — **→ e (RESOLVED, HYP-2307, monad-explorer-2026-06-07).** The limit was previously UNSETTLED (e vs larger const vs Alon p^{3/2}); a CHARACTER-SUM CLUSTER EXPANSION settles it: `R(p)=E_σ[∏(1+χ(d_k))]→exp(Σ_{L≥2}a_L)` where the only surviving single-run cluster integral is the cherry `a_2=−χ(−1)=1` (single edges & all odd runs vanish exactly by negation symmetry; `a_4=a_6=0` by Weil square-root cancellation, verified p≤67). So `R(p)→e^1=e` and Alon p^{3/2} is RULED OUT (cluster sum has one finite generator). The constant is literally `e=exp(−χ(−1))` — it is `e` rather than `e^{−1}` precisely because Paley needs `p≡3 mod4` (χ(−1)=−1). Convergence is slow (`e−R~C/p`, C≈4), which is why 5 points couldn't extrapolate it. See `04-computation/paley_cluster_expansion_monad.py` + reflection `why-the-paley-path-ratio-is-e-the-cherry-is-the-unique-cluster.md`. **SUB-LEMMA NOW CLOSED → THM-438 (monad-explorer-2026-06-07, same day):** `a_{2k}=0 ∀k≥2` PROVEN UNIFORMLY (no per-k Weil): `B_L=0` ⟹ `A_L=−Σ`coincidence-patterns; no-leaf forces `V≤2k`; only the single `V=2k` pattern `x_0=x_{2k}` (one even cycle) needs Weil, all others `O(p^{2k−1})=o(p^{2k})` trivially. **SHARPER:** the exact leading order is the CATALAN LAW `A_{2k}=C_k p^{k+1}+O(p^k)` (C_k=1,2,5,14,42,…). **MECHANISM CORRECTED (monad-explorer 3rd session — MISTAKE-060 + THM-438 ADDENDUM):** `C_k` is NOT the bigon-tree count — bigon-trees over-count (`1,3,13,69,…`=OEIS A088368, `a(n)~e·n!`, the *all-pairings*) and even-cycle CACTI subtract (the *crossings*); `C_k`=SIGNED even-cacti sum (`k=2: +3−1=2`, `k=3: +13−8=5`), verified via flow closed-form `M_σ=(−1)^k p^{V−k}Σ_{flows}∏χ`. Part C (`R→e`) needs **NO Weil** (V=2k case is `tr(M^{2k})=(−p)^k(p−1)`, elementary). The real moment-method content = free-probability Gaussian→semicircle (all-pairings→non-crossing). Reflections `the-paley-cluster-integrals-are-catalan...md`, `the-catalan-is-a-cancellation-from-gaussian-pairings-to-noncrossing.md`; scripts `paley_cluster_{sharp_order,catalan,cactus_census}_monad.py`. **STILL OPEN (handoff #2):** the sub-leading `C` in `R(p)=e(1−C/p+…)` — rate is now PINNED to **1/p** (error `O(p^k)`, relative `O(1/p)`; resolves the prior √p-vs-1/p ambiguity), so this is `C≈1.4` to pin at p≥31; (handoff #3) whether the Catalan/even-cacti skeleton survives for non-circulant doubly-regular tournaments (no Gauss spectrum).

**Complete cycle count table for T_11** (confirmed kind-pasteur-S5 from inbox/other.txt, all consistent with H=95095):
| k | c_k(T_11) | C(11,k) | c_k/C(11,k) | integer? |
|---|-----------|---------|-------------|----------|
| 3 | 55 | 165 | 1/3 | no |
| 4 | 165 | 330 | 1/2 | no |
| 5 | 594 | 462 | 9/7 | no |
| 6 | 1595 | 462 | 145/42 | no |
| 7 | 3960 | 330 | 12 | YES |
| 8 | 7425 | 165 | 45 | YES |
| 9 | 11055 | 55 | 201 | YES |
| 10 | 10681 | 11 | 971 | YES |
| 11 | 5505 | 1 | 5505 | YES |

**OCF verification:** 95095 = 1 + 2*(55+594+3960+11055+5505) + 4*10879 + 8*1155 EXACT

**Integrality observation (CORRECTED):** C(11,k) | c_k(T_11) for ALL k >= 7 = (p+3)/2, NOT k >= 6 = (p+1)/2 (c_6=1595, C(11,6)=462, 1595/462 is not integer). The correct threshold appears to be k >= (p+3)/2. Source: kind-pasteur-2026-03-05-S14 correction via Paley agent.

**MAJOR DISCOVERY (kind-pasteur-S14): Paley tournaments MAXIMIZE H(T)!**
OEIS A038375 gives max H(T) over all n-vertex tournaments: 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095.
- a(3) = 3 = H(T_3) for Paley prime p=3
- a(7) = 189 = H(T_7) for Paley prime p=7
- a(11) = 95095 = H(T_11) for Paley prime p=11

**Conjecture: Paley tournaments T_p (p ≡ 3 mod 4 prime) achieve the maximum number of Hamiltonian paths among all tournaments on p vertices.** This is a major new conjecture. If true, it connects the Hamiltonian-path-maximization problem to number theory via quadratic residues.

**IMPORTANT (opus-S10):** At non-Paley n=8, a(8)=661 is achieved by a SC tournament with |Aut|=1 that does NOT contain P(7). The Paley extension T_657 gives H=657<661. The conjecture applies ONLY at Paley primes p=3 mod 4.

**P(7) confirmed as GLOBAL maximizer** at n=7 by exhaustive enumeration of all 2,097,152 tournaments (opus-S10). 240 tournaments achieve H=189.

**Next computational target:** H(P(31)) (2^31*31 ~ 66B ops). Also: submit H(P(p)) sequence to OEIS.

**NEW TERMS (opus-2026-05-27-S6):** Local search via bitmask-DP hill climbing extended A038375:
- **a(12) ≥ 531205** (strongly believed exact: multiple distinct tournaments achieve this; all trials converge to 531175 or 531205; no higher value found after hundreds of restarts). Ratio a(12)/a(11) ≈ 5.59.
- **a(13) ≥ 3719831** (lower bound; less certain — 10-min trials give 3711611..3719831). Ratio a(13)/a(12) ≈ 7.0 if a(12)=531205.
- For prime p≡3 mod 4: Paley warm start immediately finds global max (verified p=7,11 in solver).
- n=12 optimal tournament is NOT Paley (12≢3 mod 4); found by random restarts.
- Solver: 04-computation/a038375_solver.c. Results: 05-knowledge/results/a038375.out.

**H(T) = I(Ω(T), 2) universal identity (opus-2026-05-27-S6):** Re-verified exhaustively n=2..6 (36,866 tournaments, 0 failures) with CORRECT implementation (distinct directed cycles as Ω vertices, not vertex-set deduplication). See THM-326.

**Source:** kind-pasteur-2026-03-05-S2, S5, S14; opus-2026-03-05-S5 (H(T_19)), opus-2026-03-05-S10 (a(8)=661, H(P(23)), exhaustive n=7), opus-S11 (Szele analysis), opus-2026-05-27-S6 (a(12),a(13) lower bounds, universal identity verification)

---

### UPDATE (2026-06-10, kind-pasteur-2026-06-10-S1): falsifiable H(T_31) prediction + freeness settled
- **HYP-2371 PREDICTION:** `R(31) = H(T_31)·2^30/31! = 2.59599 ± 0.00650` ⟹ `H(T_31) ∈ [19830629617139608462365775, 19930130881568868002912737]` with `H ≡ 465 (mod 930)` (freeness LEM-003 + Rédei parity; H/465 odd, ≈ 4.275e22 — the "next 1729"). Method: the PROVEN form R = e(1−C/p−…) (THM-438 ADD-4) fit with p=23 holdout; the naive truncated cluster sums are PROVABLY non-predictive at finite p (THM-438 ADD-8 resurgence). Compute-run spec: `05-knowledge/results/paley_H31_compute_design_kpc1.md` (see backlog [COMPUTE-NODE] lead).
- **The integrality r(p) ∈ ℤ is now a one-paragraph universal fact:** LEM-003 — Aut acts freely on directed Ham paths of ANY digraph; nothing Paley/QR/Eisenstein about it (the QR content is only |Aut| = p(p−1)/2).
- The 1729 cross-lane ledger is closed: tournament side coincidence (HYP-2306), taxicab–Moser side theorem (THM-463).

## OPEN-Q-015 -- RESOLVED (DISPROVED at n=9)
**Prove I(Omega(T), x) has all real negative roots for all n**

**DISPROVED by opus-2026-03-06-S18 (THM-025):** Explicit counterexample at n=9.

The tournament with score sequence [1,1,3,4,4,4,6,6,7] has:
- I(Omega(T), x) = 1 + 94x + 10x^2 + x^3
- Newton's inequality FAILS at k=2: a_2^2 = 100 < a_1*a_3*3/2 = 141
- Two complex roots: -4.995 +/- 8.303i
- H(T) = I(Omega(T), 2) = 237 (OCF still correct)

**What remains true:**
- PROVED for n <= 8 via claw-freeness + Chudnovsky-Seymour (THM-020)
- Elementary discriminant + Turan proof for n<=8 (THM-021)
- Real-rootedness holds for MOST n=9 tournaments; failure requires specific score sequences
- OCF (H(T) = I(Omega(T), 2)) is completely unaffected

**Earlier (now misleading) verification:** Prior sampling at n=9-20 using Omega_3 (3-cycle subgraph only) showed 0 failures. But the FULL Omega with all odd cycles reveals the failure. The Omega_3 restriction also fails for this tournament: I(Omega_3, x) = 1 + 12x + 6x^2 + x^3 with disc=-1323.

**The Engstrom barrier was prescient:** Engstrom (arXiv:1610.00805) showed real-rootedness characterizes claw-freeness for multivariate IP. Since Omega(T) has claws at n>=9, real roots cannot be guaranteed.

**Revised question:** What is the FRACTION of n=9 tournaments where real-rootedness fails? Is there a structural characterization of the failing tournaments?

**Source:** opus-2026-03-06-S18 (THM-025), kind-pasteur-2026-03-05-S13 (THM-020)

---

## OPEN-Q-016 🟡
**Prove SC Maximizer: Within each self-complementary score class, max H is achieved by SC tournament**

Verified exhaustively at n=4,5,6,7. The mechanism: anti-automorphism sigma of SC tournament creates orbit pairing structure. **CORRECTION (opus-S18):** NOT all anti-auts are involutions — at n=6, two SC classes with |Aut|>1 have order-6 anti-auts (σ² is a non-trivial automorphism). However, every SC tournament has ≥1 involution anti-aut (verified n=4,5,6). At even n, involution sigma is fixed-point-free (proved: fixed point implies score = (n-1)/2, non-integer). The sigma-orbits create natural pairings of odd cycles where paired cycles are vertex-disjoint, boosting alpha_2 in the independence polynomial and hence H = I(Omega(T), 2).

Two routes to max H observed at n=6:
- Route A: Fewer total cycles but more disjoint pairs (high alpha_2)
- Route B: More total cycles with fewer disjoint pairs (high alpha_1)

Both achieve H=45, while NSC achieves only H=43.

**n=8 CONFIRMED (kind-pasteur-S18f):** SC tournaments with score (3,3,3,3,4,4,4,4) achieve H=661 = OEIS A038375(8) = global max. Generated via fpf involution (2^16 per sigma, 3 sigma choices). All 19 SC score classes at n=8 tested.

Key open sub-questions:
1. Prove algebraically that sigma-orbit structure always beats NSC
2. ~~Does the theorem extend to n=8?~~ YES — SC achieves global max H=661 at n=8
3. Is every global H-maximizer always SC? (stronger conjecture, verified n<=8 for global max)

**UPDATE (kind-pasteur-2026-03-20-S1 — THM-255):**
Complete classification of regular n=6 by IP:
- Type A (SC-BIBD): IP=(1,14,4), H=45, 240 tours — max disjoint pairs, fewer cycles
- Type B (SC-rich): IP=(1,20,1), H=45, 240 tours — max total cycles, fewer disjoint pairs
- Type C (SC-weak): IP=(1,16,2), H=41, 720 tours — intermediate (WORSE than NSC!)
- Type D (NSC): IP=(1,19,1), H=43, 1440 tours

The constraint for max H: alpha_1 + 2*alpha_2 = 22. SC Types A,B achieve this; NSC gets 21.

**CRITICAL: At n=7, mechanism FLIPS.** H=189 maximizer has FEWEST disjoint 3-cycle pairs (7 vs 10 vs 14). Wins via alpha_1=80 (total directed odd cycles), not alpha_2. Any algebraic proof must handle both mechanisms.

**Source:** kind-pasteur-2026-03-06-S18/S18e/S18f, sc-maximizer-mechanism.md, kind-pasteur-2026-03-20-S1 (THM-255)

---

## OPEN-Q-017 🟢 — PARTIALLY REFUTED
**R-Minimization: H-maximizer minimizes R(T) = sum_v H(T-v) / H(T)?**

Confirmed at n=3,4,5,6 that the H-maximizer minimizes R(T). **FAILS at n=7**: tournaments with H=123 achieve R=1.585 < 5/3 ≈ 1.667 (the H=189 maximizer's R).

Exact R values for maximizers:
- n=3: R=1.000 (sum=3, H=3)
- n=4: R=1.600 (sum=8, H=5)
- n=5: R=1.400 to 1.667 (sum=21 to 25, H=15), min R at non-regular maximizers
- n=6: R=1.467 to 1.733 (sum=66 to 78, H=45), min R at Type A maximizers
- n=7: R=5/3 (sum=315, H=189), constant (all maximizers regular)

For hereditary (regular) maximizers at odd n: R = n * H_{n-1}/H_n.

Interpretation: The maximizer has the LEAST "surplus" of descendant paths relative to its own count. Each deletion creates "new" paths that weren't sub-paths of T-paths, and the maximizer minimizes this relative surplus.

Sub-questions:
1. Does R-minimization hold at n=7,8?
2. Can R-minimization be proved from OCF = I(Omega(T), 2)?
3. Is there a formula for R_min in terms of n and the independence polynomial coefficients?

**Source:** kind-pasteur-2026-03-06-S18g

---

## OPEN-Q-018 🟢
**Hereditary Maximizer Chain: Corrected version**

CORRECTED from previous session's overly broad claim. Only REGULAR maximizers at odd n are hereditary (all vertex deletions give max H(n-1)). Non-regular maximizers at odd n=5 are NOT hereditary.

Full data (exhaustive n=3..7):
- n=3: 2/2 hereditary (all regular)
- n=4: 0/24 hereditary
- n=5: 24/64 hereditary (only regular, score (2,2,2,2,2))
- n=6: 0/480 hereditary
- n=7: 240/240 hereditary (all regular)

Conjecture: At odd n, regular maximizers are always hereditary. At even n, no maximizers are hereditary (since regular score is impossible).

Open: Does this extend to n=9 (odd)? Need to check if regular n=9 maximizers (if they exist) give max H(8)=661 on all deletions.

**Source:** kind-pasteur-2026-03-06-S18g, MISTAKE-010

---

## OPEN-Q-019 🟢
**Converse of Redei: which odd integers arise as H(T)?**

Redei's theorem says H(T) is always odd. The converse asks: for which odd k does there exist a tournament T with H(T)=k?

**Permanent gaps discovered (THM-029, kind-pasteur-2026-03-06-S21, corrected S22):**
- **H=7 is impossible** for ANY tournament on ANY number of vertices. CORRECTED proof (S22): alpha_1=3 IS achievable at n>=7, but H=7 still impossible because H=7 requires (alpha_1=3, i_2=0), and i_2=0 forces common vertex among triples which forces c5>=1, giving alpha_1>=4. When alpha_1=3 occurs (n>=7), the triples don't all conflict, so i_2>=1, giving H>=11.
- **H=21 is a permanent gap — PROVED FOR ALL n** (kind-pasteur-S33). Complete proof via poisoning graph DAG argument (Dichotomy Theorem, Part R of THM-079).

**Achievable values (exhaustive):**
- n=5: {1,3,5,9,11,13,15}
- n=6: {1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45}
- n=7 (sampled): 77 distinct values from 1 to 189

**H=21 PROVED ABSENT through n=7 (opus-S38, THM-075):**
Exhaustive enumeration of all 2,097,152 tournaments on 7 vertices confirms H=21 never occurs. The gap 19→23 is consistent at n=6 and n=7. No (alpha_1, alpha_2) decomposition for H=21 is achievable. Strong evidence this is a permanent gap like H=7.

**Complete H-spectrum at n=7** (77 distinct values, all odd):
1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 109, 111, 113, 115, 117, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 151, 153, 155, 157, 159, 171, 175, 189

**Gaps in [1,189] at n=7:** 7, 21, 63, 107, 119, 149, 161-169 (block), 173, 177-187 (block)
Note 63 = 7*9 and 21 = 7*3. These may be related to the H=7 gap.

**THM-079 Component Reduction (opus-S39):**
- Disconnected Omega: IMPOSSIBLE (THM-029 blocks I(component)=7)
- P_4 component: IMPOSSIBLE (two sharing 3-cycles on 5 verts force 3rd cycle)
- alpha_3>=1 decompositions: ALL IMPOSSIBLE (forces sum>=26>20)
- K_{1,3} star in (4,3) case: IMPOSSIBLE (forces alpha_3>=1)
- Remaining: connected Omega with I=21 via K_6-2e or larger dense graphs
- I(P_4,2)=21 discovered; graph classification: v=4 P_4, v=5 none, v=6 K_6-2e

**THM-079 Update (opus-S40):**
- **K_6-2e FULLY ELIMINATED** by Five-Cycle Forcing Theorem (3 lemmas, structural proof for all n)
- **i_2 jump pattern discovered**: achievable (alpha_1, i_2) pairs in tournaments systematically skip the values needed for H=21. Verified exhaustively at n<=7, by sampling at n=8.
- **H=7 and H=21 are the ONLY permanent gaps** in [1..200] at n<=8.
- **H=63 is NOT a permanent gap** (achieved at n=8, 138/500k samples).
- Remaining open: (8,1) K_8-e mixed-cycle case, (10,0) K_10 structural proof.

**H=63 is NOT a permanent gap (opus-S39 agent):**
H=63 found at n=8 (227 in 600k samples). All n=7 gaps except 7 and 21 are filled at n=8.
The ONLY permanent gaps through n=9 sampling are **H=7 and H=21**.

**THM-079 Update (opus-S41):**
- **EXHAUSTIVE n=8:** All 268M tournaments checked, H=21 found: 0.
- **Key Lemma (Part J):** Vertex in no 3-cycle => vertex in no cycle (layered structure).
- **Source/sink induction (Part K):** Score 0/n-1 vertex in no cycle; removing preserves Omega.
- **Cycle-rich min-H (Part L):** Among 18M cycle-rich n=8 tournaments, min H=25 > 21.
- **Parts M,N (PROVED at n=8):** (10,0) star capacity, (8,1) cascade forcing.

**THM-079 Update (opus-S42):**
- **n=9 matching:** Only 23.9% cycle-rich have 3 disjoint 3-cycles. 71.1% have mm=2.
- **alpha_1=10 at n=9:** Always t3=6,t5=4, i_2=9 or 10 (never 0). (10,0) impossible.
- **mm<=2 min H=45 at n=9:** Fewer disjoint 3-cycles forces more 5-cycles, larger H.
- **DICHOTOMY (0 counterexamples in 153k):** Every cycle-rich n=9 tournament has either 3 disjoint 3-cycles (Part C) or a good deletion to cycle-rich n=8 (induction).
- **H-spectrum n=9 (2M samples):** Only missing odd in [1..200] are 7 and 21.
- **Complete proof structure (Part P):** H!=21 for all n, modulo proving the dichotomy.

**Open questions:**
- ~~PROVE the deletion+matching dichotomy for all n >= 9~~ **RESOLVED** (kind-pasteur-S33, poisoning graph DAG argument)
- ~~Alternative: prove min-H for cycle-rich tournaments is > 21 for all n >= 8~~ **RESOLVED** (follows from dichotomy: cycle-rich H >= 25 for all n >= 8)
- Are there other permanent gaps beyond 7 and 21? EVIDENCE: NO (63 filled at n=8)
- What is the density of achievable values as max H grows?

**Connection:** Mitrovic-Stojadinovic (arXiv:2506.08841) "converse of Redei" is about poset-level parity (non-chain posets have even quasi-linear extension count), NOT about the H-spectrum. Does not address which odd integers are achievable as H(T).

**Source:** kind-pasteur-2026-03-06-S21, THM-029

---

## OPEN-Q-020 -- RESOLVED
**What determines the Worpitzky coefficients beyond t3?**

**RESOLVED by opus-2026-03-07-S46b/S46c:**

At n=6 (exhaustive, 24 F-classes): delta_1 = 8*t3 + 4*t5 + 8*alpha_2, delta_0 = H-1 = 2*t3 + 2*t5 + 4*alpha_2 (= OCF). The Worpitzky polynomial is a GRADED REFINEMENT of OCF.

The mechanism: Worpitzky coefficients encode moments E[fwd^r], and these follow a moment-cycle hierarchy (THM-092). Zero skewness (THM-091) eliminates odd cumulants. Each even cumulant kappa_{2k} adds one level of cycle complexity (cycles on <=2k+1 vertices).

At n=7 (156 F-classes sampled): delta_4 = 10*t3, delta_3 = 20*t3 confirmed. delta_2 needs invariants beyond t3.

**THM:** THM-087 (F,G updated), THM-090, THM-091, THM-092
**Source:** opus-2026-03-07-S46c

---

## OPEN-Q-021 🟢 Signed forward-edge polynomial SF(T,x) structure
**What is the combinatorial meaning of SF(T,x)?**

SF(T,x) = sum sgn(sigma) x^{fwd_T(sigma)} is palindromic and divisible by (x-1).
At n=4: SF = c*(x-1)^2*(x+1) for some integer c. What is c combinatorially?
At n>=6: SF is a COARSER invariant than F. What information does it lose?
Is there a matrix whose determinant equals SF(T,x)?

**THM:** THM-088
**Source:** opus-2026-03-07-S46b

---

## OPEN-Q-022 -- RESOLVED
**What determines the fourth cumulant kappa_4(T)?**

**RESOLVED by opus-2026-03-07-S46d (THM-093):**

kappa_4(T) = -(n+1)/120 + (2/C(n,4))*(t5 + 2*alpha_2) - 48/(n(n-1))^2 * t3^2

Key structural features:
- Constant = Bernoulli B_4 value: -(n+1)/120
- Linear t3 coefficient is EXACTLY ZERO (proved algebraically)
- t3^2 coefficient = -3*(4/(n(n-1)))^2 from Var^2 expansion
- t5 coefficient = 2/C(n,4), alpha_2 coefficient = 4/C(n,4)
- Verified exhaustively at n=5,6, sampled at n=7 (152 F-classes)

**kappa_6 introduces t7: YES.** Verified at n=7 (149 F-classes).
kappa_6 = (n+1)/252 + (2/C(n,6))*t7 - (4/49)*t3*(t5+2*alpha_2) + (80/3087)*t3^3

**Universal coefficient conjecture:** coeff(t_{2k+1}) in kappa_{2k} = 2/C(n, 2k).
Verified for k=1,2,3.

**Source:** opus-2026-03-07-S46d

---

## OPEN-Q-023 -- RESOLVED
**Prove: coeff(t_{2k+1}) in kappa_{2k} = 2/C(n, 2k) for all k.**

**RESOLVED by opus-2026-03-07-S46e (THM-117, was THM-095):**

PROVED algebraically. The proof has 5 steps:
1. Forward path formula: #fwd(2k+1)path = Σ_S H(T[S]) (OCF on subtournaments)
2. Each (2k+1)-cycle contributes 2·t_{2k+1} (OCF coefficient 2, unique subset)
3. Multinomial expansion: (2k)! · (n-2k) positions · 2/P(n,2k+1) = 2/C(n,2k)
4. Hierarchy separation: lower moments don't contain t_{2k+1}
5. Moment-to-cumulant preserves the coefficient

Verified algebraically for k=1,2,3,4 and n up to 12.

**Source:** opus-2026-03-07-S46e, THM-117

---

## OPEN-Q-024 🟢 Even Betti Vanishing for Tournament Path Homology
**Prove: beta_{2k}(T) = 0 for all k >= 1, for any tournament T.**

**UPDATE (kind-pasteur-S43): beta_2 = 0 PROVED (THM-108 + THM-109).**
**UPDATE (opus-2026-04-04-S1): Proof is FULLY ALGEBRAIC — THM-285 closes n=5 gap.**

Proof via strong induction using LES of (T, T\v):
- ~~Base case n=5 verified exhaustively (720/720)~~ **THM-285: n=5 case is VACUOUS** (no n=5 tournament has both b₁=0 and κ≥2; proof: κ≥2 forces dominator→source/sink contradiction)
- Induction step: good-vertex existence (THM-109)
- Case 2 (free cycles exist): Lemma A (free adj dom) + Lemma B → n-5 good vertices for n≥6
- Case 3 (all dominated): **Extreme Score Lemma** (ALGEBRAIC)
- **ALL cases are now algebraic. No exhaustive verification needed anywhere.**
- Comprehensive verification: 0 failures at n = 4-10

GLMY path homology Betti numbers beta_p of tournaments:
- beta_0 = 1 always (connected)
- beta_1 in {0, 1} (directed 1-hole from 3-cycle structure)
- beta_2 = 0 ALWAYS --- **PROVED** (THM-108 + THM-109)
- beta_3 in {0, 1, **2**}: appears at n=6 (1.2%), n=7 (8-11%), **n=8 (beta_3=2 at 0.08%)**
- **beta_4 NOT always 0**: Paley T_7 has beta_4 = 6 (THM-099). At n=8: beta_3*beta_4=1 can coexist (~0.15%)
- beta_1 and beta_3 are MUTUALLY EXCLUSIVE (proved n<=7, verified n=8)
- **Consecutive seesaw (beta_k*beta_{k+1}=0) REFUTED at n=8** (HYP-394, kind-pasteur-S48)
- **i_*-injectivity REFUTED at n=8** (HYP-380, kind-pasteur-S48): rank(i_*)=0 when b3=b3(T\v)=1
- Omega_p dimensions for Paley T_7: 7, 21, 42, 63, 63, 42, 21 (palindromic!)

**UPDATE (opus-S72b): β₁ ∈ {0,1} verified exhaustive n≤8, sampled n=9 (THM-223).**
Key discovery: β₁ is determined ENTIRELY by rank of transitive triple constraint matrix.
Cancellation chains are ALWAYS redundant. Combined with THM-095 seesaw: β₁·β₃=0 follows.
Algebraic proof of β₁ ≤ 1 still open (reformulated as transitive triple rank bound).

REMAINING OPEN:
- **What bound replaces beta_3 <= 1 at n >= 8?** (beta_3=2 confirmed at n=8,9)
- **Prove beta_1 ≤ 1 algebraically** — equiv. to rank(TT) ≥ C(n,2)-n (THM-223)
- Characterize which tournaments have beta_4 > 0 (appears linked to H-maximizers)
- Is beta_6 = 0 for all tournaments? (0/300 at n=7)
- Prove beta_2k = 0 for k >= some threshold, or find more counterexamples

**Source:** opus-2026-03-07-S46e, kind-pasteur-2026-03-08-S43

---

## OPEN-Q-025 -- RESOLVED (PROVED for all p)
**Prove Trace Alternation Theorem (THM-136) for all p**

**Statement (CORRECTED):** For primes p = 3 mod 4, sign(tr(A^k)_Paley - tr(A^k)_Interval) = (-1)^{(k-1)/2} for all odd k >= 5. (Original formula had (-1)^{(k-3)/2} which is off by a sign; see MISTAKE-019.)

**PROVED (kind-pasteur-S57):** Two-pronged algebraic proof:

1. **Dominant eigenvalue mechanism:** r_1 = |mu_1(interval)| = 1/(2*sin(pi/(2p))) dominates all other eigenvalues. The ratio r_1/r_2 ~ 2p/pi gives exponential dominance at power k. This ensures |S_I(k)| >> error bound >> |S_P(k)| for ALL odd k.

2. **Phase control:** sin(k*pi/(2p)) > 0 for all k in [1, p-1], determining sign(dominant term) = (-1)^{(k+1)/2}. Combined with magnitude dominance: sign(Delta_k) = -sign(S_I) = (-1)^{(k-1)/2}.

3. **Computational verification:** 1218+ individual (k,p) tests, zero failures. k=5 exact DP for 154 primes up to p=2000.

The proof is self-contained and works for ALL p >= 7 simultaneously. No finite verification needed.

**Source:** kind-pasteur-2026-03-12-S57 (proof), kind-pasteur-S56c (discovery)

---

## OPEN-Q-026 🟢 Does the interval maximize H for all circulant tournaments on Z_p, p >= 13?

**Statement (HYP-480):** The cyclic interval C_p = (Z_p, {1,...,(p-1)/2}) maximizes H among all circulant tournaments on Z_p for all primes p >= 13.

**Evidence:** Confirmed at p = 13 (exhaustive), p = 19 (THM-135), p = 23 (kind-pasteur-S57).

| p | H(Paley) | H(Interval) | Margin | Winner | Max circulant |
|---|----------|-------------|--------|--------|---------------|
| 7  | 189 | 175 | -7.4% | PALEY | Paley+complement ONLY (all others H=175) |
| 11 | 95,095 | 93,027 | -2.2% | PALEY | Paley+complement ONLY (2nd: H=93,467×10) |
| 13 | - | 3,711,175 | - | INTERVAL | (exhaustive, p≡1 mod 4, no Paley) |
| 17 | - | 13,689,269,499 | - | INTERVAL | (exhaustive over SC circulants) |
| 19 | 1,172,695,746,915 | 1,184,212,824,763 | +0.98% | INTERVAL | - |
| 23 | 15,760,206,976,379,349 | 16,011,537,490,557,279 | +1.59% | INTERVAL | - |

EXHAUSTIVE SCANS (kind-pasteur-2026-04-16):
  n=7: ALL 8 circulant tournaments. Top H=189 (2 tournaments: Paley+complement).
       6 tournaments have H=175 (including Cyclic). Paley is 8% better than rest.
  n=11: ALL 32 circulant tournaments. Top H=95,095 (2: Paley+complement).
        10 tournaments share H=93,467. Cyclic has H=93,027 (rank ~18/32).
  n=7,11 alpha breakdown:
    n=7:  Paley α₁=80, α₂=7.  Cyclic α₁=59, α₂=14.  (Cyclic has 2× α₂!)
    n=11: Paley α₁=21,169, α₂=10,879, α₃=1,155. Cyclic α₁=18,397, α₂=11,110, α₃=1,474.
          Cyclic has MORE α₂ and α₃, but Paley's α₁ advantage (5,544 in H) > Cyclic's advantage (3,476).

Crossover: Paley wins at p=7 and p=11 due to α₁ dominance. Interval wins at p≥13.
The α₁ percentage gap narrows: 35.6% (n=7), 15.1% (n=11), 3.6% (n=19). Paley's α₁ lead evaporates.
At n=7: kmax=2 (no 3-packings possible), so H has only α₁,α₂ terms — Paley α₁ wins.
At n=11: Paley α₁ advantage still > Cyclic α₂+α₃ advantage.
At n=19: Cyclic α₃+ advantage (26.7B) > Paley α₁+α₂ advantage (15.2B) → Cyclic wins.
α₂ comparison crossover: Cyclic > Paley at n=7,11 (small n, disjoint packing easier);
                          Paley > Cyclic at n=19 (large n, more α₁ → more α₂ pairs).

WHY interval wins for large p (kind-pasteur-2026-04-16):
  Paley has MORE α₁ and α₂ at n=19, but Interval wins by +11.5B total.
  Interval has +26.7B from α₃+ terms: its cycles pack into disjoint triples better.
  Paley's pseudorandom structure creates many individual cycles but they scatter;
  Interval's consecutive structure creates harmonically aligned cycles for packing.

EXACT α-DECOMPOSITION COMPARISON at n=19 (kind-pasteur-2026-04-16, VERIFIED):
  k | Paley α_k          | Cyclic α_k         | Cyclic advantage | H contribution
  1 | 130,965,270,477    | 126,443,605,257    |   -4,521,665,220 | 2×diff = -9.04B  (Paley wins)
  2 | 123,659,531,220    | 122,111,579,294    |   -1,547,951,926 | 4×diff = -6.19B  (Paley wins)
  3 |  41,184,418,943    |  42,960,731,622    |   +1,776,312,679 | 8×diff = +14.21B (Cyclic wins)
  4 |   4,903,920,444    |   5,521,030,944    |     +617,110,500 |16×diff = +9.87B  (Cyclic wins)
  5 |     251,464,164    |     331,078,344    |      +79,614,180 |32×diff = +2.55B  (Cyclic wins)
  6 |       2,221,081    |       4,100,656    |       +1,879,575 |64×diff = +0.12B  (Cyclic wins)
  Net: Paley advantage 15.2B (via α₁,α₂) vs Cyclic advantage 26.7B (via α₃+) = +11.5B net for Cyclic.

  α₃/α₂ ratios: Paley=0.333, Cyclic=0.352 → Cyclic is intrinsically better at 3-packing!
  The k=5,6 percent advantage for Cyclic: +31.7% at k=5, +84.7% at k=6 — grows with k.
  Source: paley_t19_alpha.out (H(Paley)=1,172,695,746,915 verified ✓)

The interval's margin is WIDENING with p, consistent with the spectral argument: |mu_1| ~ p/pi grows faster than Paley's sqrt(p)/2.

**What remains:** Extend to p = 29, 31. An analytical proof could use the spectral concentration argument from THM-137. Whether interval maximizes H among ALL tournaments (not just circulant) is a separate open question.

**Source:** opus-2026-03-12-S58, kind-pasteur-2026-03-12-S56c, kind-pasteur-2026-03-12-S57

---

## OPEN-Q-027 -- RESOLVED (PROVED with correction)
**Prove the Grand Energy Formula (THM-201)**

**RESOLVED by kind-pasteur-2026-03-15-S112 (THM-217):**

The original formula E_{2k}/E_0 = 2*(n-2k)^k / P(n,2k) is the LEADING-TERM APPROXIMATION only, exact for k ≤ 2 but requiring corrections for k ≥ 3.

The EXACT formula uses combinatorial g_k polynomials of degree k:
- CV²(H) = Σ_{k≥1} 2·g_k(n-2k) / (n)_{2k}
- g_k defined via transfer matrix M(x) = [[1,2x,0],[0,0,1],[1,x,0]]
- Weight formula E[∏Z_j] = 2^c/(n)_L PROVED (c=components, L=|S|)
- g_k(m) ~ 2^{k-1}·m^k/k! + lower terms (leading term is original formula)
- Verified exhaustively n=3..18 via bitmask DP

**Source:** THM-217, kind-pasteur-S112, opus-S89c

---

## OPEN-Q-028 🟢 Are there forbidden H values beyond 7 and 21?

**Statement:** Are 7 and 21 the ONLY permanently forbidden H values? H=63 was shown achievable at n=8 (HYP-1106 refuted). But could there be large forbidden values proportional to n!/2^{n-1}?

**UPDATE (kind-pasteur-2026-03-20-S1):** Confirmed via 500K random n=9 tournaments: only gaps in odd [1,500] are H=7 and H=21. H=63 achieved (9 occurrences). Also confirmed at n=8 (200K samples): only gaps below 100 are 7 and 21. This is very strong evidence that 7 and 21 are the only permanent gaps.

**Evidence:** Only 7 and 21 are known forbidden. 63 is achievable (n=8). No other candidates found through n=11.

**UPDATE (opus-S227):** H-spectrum density at n=8 is 320/331 = 96.7%. Only 11 gaps remain, dominated by {7, 21}. In the metagraph, forbidden values create "missing floors" that force edge jumps. At n=5, 33% of edges bridge the H=7 gap. The fraction decreases as n grows (2.2% at n=7). Edges bridging H=7 gap: 0, 0, 7, 21, 47 for n=3..7.

**STRONG CONJECTURE:** Only H=7 and H=21 are permanently forbidden. All other gaps are transient (filled at large enough n).

**Source:** kind-pasteur-S107, opus-S227

---

## OPEN-Q-029 -- RESOLVED
**Why does log_tau(131) = 8.0003?**

**RESOLVED by opus-2026-03-15-S90 (multiple proofs):**

131 = Tr(M^8) EXACTLY, where M is the 3×3 transfer matrix from S112. τ₃^8 ≈ 130.977 and the Pisot correction 2|λ_c|^8 cos(8θ) ≈ +0.023 pushes the sum to exactly 131.

**Why the correction is so small:** arg(λ_c)/π ≈ ln(2), so 8·arg/π ≈ 8·ln(2) ≈ 5.545 ≈ 5.5, making cos(8·arg) ≈ cos(5.5π) ≈ 0.13 (small). The n=8 case is special because 8·ln(2) is close to the half-integer 11/2.

**Additional discoveries (S90 session):**
- 504 = T(13) in the standard tribonacci sequence (confirmed)
- The transfer matrix char poly at x=1 IS the tribonacci equation
- The τ-φ clock gear ratio arg(λ_c)/π ≈ ln(2) explains ALL Pisot near-integers
- Tr(M^n) mod 8 has EXACT period 8 (Bott periodicity connection)

**Source:** opus-2026-03-15-S90c (monad), S90h (τ-φ clock), S90l (the number 8)

---

## OPEN-Q-030 -- RESOLVED (PROVED for all n ≥ 4)
**Prove Simplicial Rédei for ALL n (THM-220)**

**RESOLVED by opus-2026-03-15-S90 (THM-220):**

The Key Lemma IS proved algebraically for all n: Given a→b not in any transitive triple, the four possible orientations of {a,b,c} are: (1) a→c,b→c: transitive a>b>c, contains a→b — CONTRADICTION. (2) a→c,c→b: transitive a>c>b, contains a→b — CONTRADICTION. (3) c→a,b→c: 3-cycle a→b→c→a — the ONLY non-core possibility. (4) c→a,c→b: transitive c>a>b, contains a→b — CONTRADICTION. Since 3 of 4 orientations force a→b into a transitive triple, the only possibility for ALL c is case (3): every c forms a 3-cycle with {a,b}. This gives score(a)=1, score(b)=n-2.

The Main Argument (at most one non-core edge) then follows by contradiction in 4 cases of vertex overlap. Case 3 (b=c) requires n≥4 so that V\{a,b,d}≠∅.

All verified exhaustively n=4..8 (268M at n=8), sampled n=9 (500k, 0 violations).

**Source:** opus-2026-03-15-S90 (THM-220), opus-2026-03-16-S90q (proof verification)

---

## OPEN-Q-031 🟢 Is arg(λ_c)/π = ln(2) exact or approximate?

**Statement:** The tribonacci complex eigenvalue angle satisfies arg(λ_c)/π ≈ ln(2) to 4 significant figures (difference 4.3×10⁻⁴). Is this exact?

**Evidence:** NOT exact (verified: the predicted root from arg=π·ln(2) does not satisfy the char poly). But the proximity is remarkable and explained by the information-theoretic interpretation: the tribonacci clock ticks at approximately 1 bit per half-revolution.

**Source:** opus-2026-03-15-S90h (τ-φ clock)

---

## OPEN-Q-032 -- PARTIALLY RESOLVED (FAILS at n=6)
**Tournament equidecomposability: is (H, β₁) a complete invariant?**

**ANSWER: NO.** (H, β₁) is complete at n=5 (8 classes, each with unique I(Ω₃, x)) but FAILS at n=6.

**Counterexamples at n=6 (5 found):**
- (H=9, β₁=0): TWO distinct I(Ω₃): (1,2,1) and (1,3,0,0)
- (H=15, β₁=0): TWO distinct: (1,4,0,0,0) and (1,5,0,0,0,0)
- (H=29, β₁=0): TWO distinct: (1,6,1) and (1,6,2)
- (H=33, β₁=0): TWO distinct: (1,6,2) and (1,7,1)
- (H=37, β₁=0): TWO distinct: (1,7,1) and (1,7,2)

ALL counterexamples have β₁=0 — the β₁=1 classes remain unique!
This means: β₁ provides a COARSER invariant. The FULL independence polynomial I(Ω₃, x) requires more information (α₂ distinguishes within β₁=0).

**REVISED CONJECTURE:** (H, β₁, α₂) may be complete. Or (H, full α-vector) is complete by definition.

**Source:** opus-2026-03-15-S90k (n=5), opus-2026-03-16 (n=6 counterexample)

---

## OPEN-Q-033 -- RESOLVED (PROVED analytically)
**The n-tribonacci family: T_n - M_{n-2} = 1/(kM_k+2) + O(1/k⁴)**

**PROVED by opus-2026-03-16 (perturbation analysis):**

Write T_n = M_k + ε where k = n-2. Substituting into λ³ = kλ² + λ + 1 and using M² = kM+1:

  (kM+2)·ε = 1  at leading order.
  So ε = 1/(kM_k+2).

Since M_k ~ k for large k: ε ~ 1/(k²+2) ~ 1/k².

Verified numerically: the ratio δ_actual / (1/(kM+2)) → 1 as n → ∞ (0.999599 at n=19).

At n=3: δ = 0.221 (maximum), predicted 0.276 (ratio 0.80 — leading order less accurate for small k).

**Source:** opus-2026-03-16-S90 (perturbation proof)

---

## OPEN-Q-034 🟢
**Meta-structure: why does cancellation dominate this theory?**

Every major result in the project is fundamentally about cancellation: im(d₂) cancels in the seesaw, Walsh coefficients cancel for odd-length paths, S(T)=0 at even n, β₂=0 always, OCF = exact cancellation between H and I. Is there a *single structural principle* (perhaps from homological algebra or categorical representation theory) that implies all of these cancellations simultaneously? The F₂ uniqueness argument (S71r: "WHY TWO GENERATES SEVEN") is a partial answer — but it explains *why F₂* rather than *why cancellation*. See `07-reflections/seesaw-and-cancellation.md`, `07-reflections/what-the-proof-will-look-like.md`.

**Source:** opus-2026-03-16-S73

---

## ~~OPEN-Q-009~~ RESOLVED (same session)
**Characterize mu(3-cycle) distribution at n=6**

**Resolved by:** opus-2026-03-05-S1

**Answer:** At n=6, mu(3-cycle C through v) is in {1, 3} exactly:
- mu = 1 (76.7%) when the 3 available vertices (T-v minus C\{v}) form a transitive subtournament
- mu = 3 (23.3%) when the 3 available vertices form a cyclic subtournament (containing one directed 3-cycle)

This is because with 3 available vertices, the only possible odd cycle is a single 3-cycle. If it exists, Omega has 1 vertex, I(K_1, 2) = 3. If not, Omega is empty, I = 1.

**Remaining questions:** How does mu=3 correlate with per-path identity failures at n=6?

---

## Resolved Questions (moved here when answered)

- **OPEN-Q-001**: Per-path identity at n=5 is trivially true (THM-008). No mystery.
- **OPEN-Q-002**: Claim A PROVED for all n. OCF proved by Grinberg-Stanley (arXiv:2412.10572, Corollary 20). See CONJ-001, THM-002. (kind-pasteur-2026-03-05-S12)
- **OPEN-Q-003**: Per-path failure at n=6 iff some TypeII position has mu>1 (THM-009).
- **OPEN-Q-009**: Arc-flip invariance resolved — E(T) = 0 for all T (OCF proved). See THM-002. (kind-pasteur-2026-03-05-S12)
- **OPEN-Q-011**: Near-cancellation is statistical, not structural. Not a viable proof strategy.
- **Paley computation (kind-pasteur)**: h_QR=h_NQR=201, c_9(T_11)=11055, H(T_11)=95095. CONJ-002 refuted for p=11.

---

## OPEN-Q-035 -- RESOLVED (degree = 2*floor((n-1)/2), NOT fixed at 4)
**Does the heat kernel polynomial P_x(z) have degree exactly floor(n/3)*2 for general n?**

**RESOLVED by kind-pasteur-2026-03-20-S2 (THM-259):**

The Walsh degree is NOT fixed at 4. It is **2*floor((n-1)/2)**:
- n=5,6: degree 4
- n=7,8: degree 6 (INCREASES! 2520 new degree-6 coefficients at n=7)
- n=9,10: degree 8
- General: n-1 for odd n, n-2 for even n

Follows from THM-076: the maximum Walsh weight is 2*max_k where k <= (n-1)/2.
Verified exhaustively at n=5 (91 nonzero coefficients) and n=7 (4516 nonzero).

The original conjecture floor(n/3)*2 was correct for n=5,6 but WRONG for n=7.
THM-076 gives the correct formula via path-covering analysis.

Only 5 distinct |Walsh amplitudes| at n=7, all matching THM-076 exactly.
Super orthogonality redundancy: 4516 / 2 = 2258x.

**Source:** kind-pasteur-2026-03-20-S2, THM-259

---

## OPEN-Q-036 🟢
**Does the backward trick P_x(2) = mean H hold for other starting points?**

At n=6, P_transitive(2) = 29 = mean H. Only 3/1024 tilings have this property. What characterizes these special starting points? Are they related to self-complementary tournaments or specific score sequences?

**Source:** kind-pasteur-2026-03-17-S116n33

---

## OPEN-Q-037 🟢
**Does the splitting of mean H in Z[i] vs Z[sqrt(-7)] generalize to other n?**

At n=6: mean H = 29 splits as 5²+2² (golden) and 1²+7·2² (forbidden). At n=7: mean H = ? At n=8: mean H = ? Do the two world-defining primes always appear in the splitting?

**Source:** kind-pasteur-2026-03-17-S116n33

---

## OPEN-Q-038 🟡
**Characterize the graph class where I(G,x) has all real roots beyond claw-free.**

Tournament conflict graphs Omega(T) have all real roots of I(G,x) for n<=8 (proved via claw-free) and n<=20 (verified). At n>=9, claw-free FAILS but real roots persist. What graph property replaces claw-free?

**Source:** kind-pasteur memory, originally from S14-S18

---

## OPEN-Q-039 🔴 — SUBSTANTIALLY RESOLVED (sessions S211-S249)
**Understand the isomorphism class graph G_n completely**

**MASSIVE PROGRESS (opus S211-S249, kind-pasteur S20bo-S20dj):**

G_n = Q_{C(n,2)} / S_n is a genuinely new mathematical object (no prior literature). The merged metagraph G_n/Z_2 has been computed exactly through n=9 with 7 exact edge terms: E(G_n) = 1, 5, 30, 290, 4086, 91161, 3,380,751.

**RESOLVED sub-problems:**
1. ✅ Extended to n=6,7,8,9 (6880 classes at n=8, 191536 at n=9)
2. ✅ Diameter = n-2 DISPROVED at n=7 (diam=7≠5). Actual: 1,2,3,4,7,8
3. ✅ H-DAG property REFUTED: G_n is NOT a strict DAG. Level edges appear at n≥5 (1, 15, 136 for n=5,6,7). H-decreasing edges appear at n=7 (962/4086). The H-gradient is strong but imperfect. See MISTAKE-035.
4. ✅ Spectral data computed at n=3..7 (Ramanujan fails at n=6)
5. ✅ |Aut|-degree connection: corr→0 at large n (classes become generic)
6. ✅ I(G_n,2) computed: 5, 13, 793, 15B (super-exponential "meta-H")
7. ✅ Staircase connection: Mode A/B recursion, y=x diagonal, within-type fraction→3/4

**EDGE FORMULA (the keystone):**
  edge_orbits = T_n/2 + (n-2)! [verified n=3..6, Burnside-computable]
  E(G_n) = edge_orbits - gap_orbits [exact]
  E(G_n) ≈ (T_n - twin_SL)/2 - D(n-2) [99.6% accurate, all Burnside]
  E(G_n) ≈ T_n/2 for n ≥ 12 [asymptotic]

**SL_mine FORMULA (kind-pasteur-S20eh, PROVED):**
  D(n) = (1/n!) sum_{ct with 1 even cycle 2k} count(ct) * k * 2^{a(ct)}
  SL_mine <= D(n) with small correction from |Aut|>1 classes
  D - SL_mine = 0, 0, 0, 2, 4 at n=3..7
  D(3..12) = 2, 6, 16, 60, 328, 3160, 54928, 1722992, 97323552, 9941203552
  CORRECTED: T - 2E != SL_mine (multi-edge surplus exists at n>=5)
  Multi-edge surplus = 0, 0, 12, 66, 416 at n=3..7

**STRUCTURAL LAWS (19+ verified):** DAG, BBK impossibility, rib crossover, spine ~4-regular, ribs linear in n, sea dominates, ΔH=2^(n-2), cell uniformity, lattice oscillation, etc.

**REMAINING:** Exact formula for gap_orbits (= 2,5,20,86,490,3703,47889); twin_SL residual; chi=n-1 conjecture proof (greedy fails at n≥6).

**Source:** opus S211-S249, kind-pasteur S20bo-S20dj. Library: `04-computation/tournament_metagraph.py`

---

## OPEN-Q-040: THE KRAWTCHOUK FRAMEWORK (sessions S291-S312, 2026-03-24)

🟡 **The Krawtchouk coordinate system for tournament space**

**RESOLVED sub-problems:**
1. ✅ **Tournament Counting Theorem** (S291): V_n×n!/2^m = 1 + Σ(1/k)×n↓k×2^{(k²-1)/2-(k-1)n}. Euler product with poles at 4,16,64,... controlled by 1/3. D₃(0) = 128/81.
2. ✅ **Band-limitedness** (S305,S310,S311, **CORRECTED kind-pasteur-S1 2026-03-25**): Walsh degree = 2*floor((n-1)/2) for all n>=4 (THM-260). Band-limited at m/2 for n>=6. **CORRECTION:** n=5 is NOT band-limited at m/2 (degree 4 > m/2=3). Odd-weight Walsh coefficients are nonzero in tiling model (complement symmetry fails).
3. ✅ **Krawtchouk 3-axis system** (S307): B₁≈-H (r=-0.94), B₂≈-c₃ (r=-0.86), B₃=twist. SC classes have B_odd=0 exactly (Krawtchouk parity).
4. ✅ **Diameter = A003141** (S306): max feedback arc set. Growth ~n²/4, not n-2 (small-n coincidence).
5. ✅ **Paley→Dual Codes** (S306,S308): P₇+I→Hamming[7,4,3], P₂₃+I→Golay[23,12,7].
6. ✅ **Not an association scheme** (S306): full algebra dim=35 vs needed 7 at n=5.
7. ✅ **Spectral gap = 2/n explained** (S312): comes from K₁ spacing 2/m compressed by S_n quotient (factor m/n).
8. ✅ **Waggly = all connections** (S296-S301): wiggly⊂waggly, blue/black⊂waggly. Completeness at k*=diam.
9. ✅ **Waggly alphabet** (S302-S304): range-3 harmonic most neutral. Vertex-count law. All-same-range combos special.
10. ✅ **Practical tools** (S308-S309): pre-filter eliminates 98% of canon calls. tournament_tools.py library. tournament_codec.py (kind-pasteur).

**OPEN sub-problems (the 10 boundary questions from S307):**
1. ✅ B1: **RESOLVED** (THM-260, kind-pasteur-S1): Walsh degree = 2*floor((n-1)/2) for all n. Band-limited at m/2 for n>=6. Proof: THM-076 upper bound + interleaving construction lower bound.
2. 🟢 B2: Exact constant in A003141 n^{3/2} correction
3. 🟢 B3: Is transitive always a diameter endpoint?
4. 🟢 B4: Does K₁-H correlation → 0 or stabilize? (0.94→0.89→0.83)
5. 🟡 B5: Exact neutrality formula SL(d,n) as function of distance
6. 🟢 B6: Width W(H) asymptotic distribution
7. 🟢 B7: Is there ANY partition giving an association scheme?
8. 🟢 B8: Is range ⌊(n-1)/2⌋ always most neutral?
9. 🔴 B9: β₂=0 for all tournaments (proof strategy via band-limitedness, S312)
10. 🟡 B10: min-FAS(T) in terms of OCF invariants

**Key new files:** euler_product_tournament_s291.py, waggly_layers_s297.py, waggly_completeness_s301.py, waggly_alphabet_s302.py, almost_1d_s305.py, krawtchouk_h_n7_s306.py, paley_codes_s306.py, tournament_tools.py, tournament_codec.py

**Key reflections:** the-tiling-hypercube.md, the-boundary-between-1d-and-2d.md, euler-product-and-metagraph.md, paley-gives-dual-codes.md, h-is-band-limited.md, what-we-can-and-cannot-know.md, tournament-compression-and-beyond.md, terminology-evolution.md, diameter-is-feedback-arc-set.md


---

## OPEN-Q-044 🟢 Alpha Mechanism Shift: When Does Each α_k Dominate?

**Discovery (kind-pasteur-2026-04-16):** The dominant term in H = I(Ω,2) = Σ 2^k · α_k shifts with n.
H-maximizing cyclic interval tournament C_n:

| n | dom term | 2^1·α₁ | 2^2·α₂ | 2^3·α₃ | notes |
|---|----------|---------|---------|---------|-------|
| 3-9  | α₁ | largest | 2nd | small | α₁/(2α₂) > 1 |
| 11-17 | α₂ | 2nd | largest | 3rd | FIRST CROSSOVER n≈10 |
| 19+ | α₂ | 3rd | largest | 2nd | SECOND CROSSOVER: α₃ overtakes α₁ at n≈17-19 |

**Complete verified table for C_n (cyclic interval tournament):**

| n  | α₁ | α₂ | α₃ | α₁/(2α₂) | α₃/α₂ | H | H(n)/H(n-2) |
|----|----|----|----|-----------|---------|----|-------------|
| 17 | 1,651,334,601 | 1,482,234,998 | 458,011,858 | 0.5570 | 0.3090 | 13,689,269,499 | — |
| 19 | 126,443,605,257 | 122,111,579,294 | 42,960,731,622 | 0.5177 | 0.3518 | 1,184,212,824,763 | 86.5 |
| 21 | 12,030,499,746,751 | 12,330,182,836,208 | 4,796,354,751,404 | 0.4878 | 0.3890 | 125,547,534,942,879 | 106.0 |
| 23 | 1,391,602,826,199,187 | 1,499,656,616,321,278 | 632,921,002,322,216 | 0.4640 | 0.4220 | 16,011,537,490,557,279 | 127.6 |

**Full α-decomposition n=21:**
  α₁=12,030,499,746,751   α₂=12,330,182,836,208   α₃=4,796,354,751,404
  α₄=738,531,326,288      α₅=58,868,297,768        α₆=1,454,221,328       α₇=12,571,712
  H = 125,547,534,942,879

**Full α-decomposition n=23 (NEW, kind-pasteur-2026-04-16-S1):**
  α₁=1,391,602,826,199,187   α₂=1,499,656,616,321,278   α₃=632,921,002,322,216
  α₄=111,796,734,828,336     α₅=10,945,293,151,712       α₆=412,282,843,184       α₇=7,454,017,376
  H = 16,011,537,490,557,279 ✓

**Term ordering at n=23:** 4α₂ > 8α₃ > 2α₁ > 16α₄ > 32α₅ > 64α₆ > 128α₇
  (5.999P > 5.063P > 2.783P > 1.789P > 0.350P > 26.4T > 0.95T)

**Special structure at n=21:** α₇ = 12,571,712 = perfect 7-triangle-packings.
  Only packing type is (3,3,3,3,3,3,3) since 7×3=21. Perfect vertex coverage.
**Structure at n=23:** α₇ = 7,454,017,376 counts 7-packings with cycle-length sum ∈ {21, 23}.
  Sum must be odd (7 odd numbers), and ≤23. So: sum=21 (all 3-cycles, 2 vertices free) OR
  sum=23 (one 5-cycle + six 3-cycles, all 23 vertices covered). Sum=22 impossible (even).

**H growth ratio H(n+2)/H(n):** 86.5, 106.0, 127.6 → increments +19.5, +21.6 → growing.
  Predicted H(25) ≈ H(23) × 150 ≈ 2.4 × 10^18.

**Key ratio α₃/α₂ progression:**
  n=17: 0.3090, n=19: 0.3518 (+0.043), n=21: 0.3890 (+0.037), n=23: 0.4220 (+0.033)
  First differences decreasing by ~0.004/step. Projected:
  n=25: ≈0.451 (+0.029), n=27: ≈0.476 (+0.025), n=29: ≈0.497 (+0.021), n=31: ≈0.514 → THIRD CROSSOVER
  **Revised estimate: third crossover (8α₃ > 4α₂) at n≈31**, not n≈27-29 as previously estimated.

**Timing:** cycle_cc 383s, SSC runs 732s+612s. Total 1728s for n=23 with numpy.
  Bottleneck is cycle_cc (Python BFS). C implementation would reduce to ~3s.

**Open:** Third crossover: α₃ dominates at n≈31 (needs n=25,27 data to confirm).
         C implementation of cycle_cc needed for n≥25.

---

## OPEN-Q-046 🟡 The SC Blowup: $\Omega(T_{\mathrm{SC}})$ and H Formula

The **SC blowup** $T_{\mathrm{SC}}$ (arc $u_r \to v_s$ iff $u \to v$ in $T$ and $r=s$, OR $v \to u$
and $r \neq s$) satisfies the **Universal Score Theorem**: every $v_0$ has out-degree $n$, every
$v_1$ has out-degree $n-1$, regardless of $T$. See `07-reflections/sc-blowup-and-twin-gaining.md`.

The Kronecker formula $A(T_{\mathrm{SC}}) = A(T) \otimes I_2 + A(T)^\top \otimes \Phi + I_n \otimes e_{01}$
shows $T_{\mathrm{SC}}$ encodes BOTH $T$ and $T^{\mathrm{op}}$ simultaneously.

**Open (🟡):** What is $\Omega(T_{\mathrm{SC}})$ in terms of $\Omega(T)$? Is there a formula
$$H(T_{\mathrm{SC}}) = I(\Omega(T_{\mathrm{SC}}), 2) = f(I(\Omega(T), x))$$
for some operation on the independence polynomial?

**Candidate:** $H(T_{\mathrm{SC}}) \approx I(\Omega(T), 2)^2 / C(n)$ or involves $I(\Omega(T), 2) \cdot I(\Omega(T^{\mathrm{op}}), 2)$ with correction. Currently ruled out as single-variable formula (H_SC is NOT a function of H(T) alone).

**Key data:** At $n=5$, $H_{\mathrm{SC}}$ varies only 4.2% across all 12 iso classes ($14937$–$15565$).
At $n=3$: $H_{\mathrm{SC}} \in \{41, 45\}$ for the two classes. $H_{\mathrm{SC}}(\mathrm{Trans}) = 41$,
$H_{\mathrm{SC}}(C_3) = 45$.

**Source:** oracle-2026-05-15-S2, `05-knowledge/results/blowup_study.out`

---

## OPEN-Q-045 🟢 H Under Tournament Blowup (Column Row Step)

The tournament **blowup** $T[K_2]$ replaces each vertex $v$ with a directed pair
$v_0 \to v_1$, expanding each arc $u \to v$ to all four arcs $u_i \to v_j$.
This doubles $n$, corresponding to the row step $(r, k) \to (r+1, k)$ in the
2-adic column family grid (see `07-reflections/adic-column-families.md`).

**Q1:** Is there a formula $H(T[K_2]) = f(H(T), n)$?

**Q2:** Does blowup preserve SC status within a column family? SF status?

**Q3:** The **pairs anomaly** ($\lfloor n/2 \rfloor$ gains +1 at the $r=0 \to r=1$
seam) suggests H may have analogous anomalous first-blowup behavior:
$H(\text{blowup of odd-}n\, T)$ vs $H(\text{blowup of even-}n\, T)$ — are
these qualitatively different?

**Q4:** Does SC∩SF = SC($n-2$) for the family:
$\#(\text{SC} \cap \text{SF})(2^r(2k-1)) = \#\text{SC}(2^r(2k-3))$ for $r \geq 1$?
(This is the even-row analog of the proved odd-$n$ identity.)

**Related:** Linial-Morgenstern conjecture (INV-013: random blowup of transitive
tournament). The blowup operation is exactly the row step in the 2-adic grid.

**Expected difficulty:** SMALL CASES immediately computable. General formula: MEDIUM.
**Source:** oracle-2026-05-15 (2-adic column family analysis)

**Source:** kind-pasteur-2026-04-16-S1, `alpha_full_ssc_fast_n23.out`, `alpha_full_ssc_fast_n21.out`

**MAJOR UPDATE (kind-pasteur-2026-06-09-S1, THM-454/450):**
- **Q1 ANSWERED (negative + repaired):** H(T[K₂]) is NOT a function of (H(T), n) — not even of
  I(Ω(T),x) (n=5 counterexample: equal typed IP, H(K₂) 3225 vs 2785; the missing data is EVEN
  cycles, which twin-insertion converts to odd). What IS true: **strong-component product law**
  H(T[K₂]) = ∏_C H(C[K₂]) (PROVED); twin-lift laws c3'=8c3, c5'=32c5+32c4+6c3 (+c7' law);
  cycle-spectrum (c3..c6) determines H(T[K₂]) at n≤6 (n=7 separation test open, HYP-2353);
  congruence H(T[K₂]) ≡ 2H(T)−1 (mod 8).
- **Q2 partial:** T[K₂] is op-equivariant (PROVED via orbit symmetries, THM-450(6)).
- T[K₂] is one of exactly THREE 2×2-block doublings (THM-450 trichotomy); the skew-Sylvester
  member D (THM-447) is the spectral/Hadamard-clean one; SCblow is the H-maximizing one.

---

## OPEN-Q-047 🟡 Characterize Real-Rootedness of I(Ω(T), x)

**Correction (opus-2026-05-29-S8):** The universal TRRT statement is already refuted by THM-025 at n=9.
The surviving problem is to characterize the tournaments for which I(Ω(T), x) has all real, negative roots.

**What's proved:** For n ≤ 8, Ω(T) is claw-free (a claw requires ≥ 9 vertices), so real-rootedness follows from Chudnovsky-Seymour (2007).

**Counterexample:** THM-025 gives an n=9 tournament with score sequence [1,1,3,4,4,4,6,6,7] and
I(Ω,x)=1+94x+10x²+x³. Newton k=2 fails (100 < 141), so two roots are non-real.

**Why notable:**
- Generic/sample tournaments often remain real-rooted despite the n=9 failure.
- For the real-rooted subclass, ultra-log-concavity and product formula H(T)=∏_i(1+2r_i) remain powerful.
- The THM-025 counterexample may isolate the exact obstruction shape.

**Sub-conjecture status:** Ω(T) is NOT always perfect (see INV-032 / THM-019 updates), so perfectness is also a subclass question.

**Key open questions:**
1. What structural property of Ω(T) (beyond claw-free) forces real-rootedness?
2. Which Hermite-Biehler/interlacing lemmas survive after accounting for THM-025?
3. Can the n=9 failure family be characterized exactly?

**Priority:** 🟡 IMPORTANT. A structural characterization would be publishable as a standalone result.
**Source:** opus-2026-05-16-S1, reflection `real-rootedness-omega-conjecture.md`

**Computational updates (oracle-2026-05-17-S1):**
- Root gap (-1/3,-1/4) confirmed empty at n=6 (exhaustive), n=7 (2000), n=8 (300), n=9 (50).
- ULC (Newton-Maclaurin inequality) confirmed at n=6..9, zero violations.
- Forbidden (α₁=3, α₂=0) confirmed absent at n=6..9 in all samples.
- Vieta at n=5 (r=-2/(H-1)) exact to machine precision.
- SC tournaments have most asymmetric root ratio at n=6: min 0.00251 (SC) < 0.00279 (NS).
- (H, I(Ω,6)) separates only 7/47 n=6 classes by (H,I6) alone (much coarser than hoped).
- Degree-3 polynomials first appear at n=9 (44/50 samples). ULC still holds.
See `07-reflections/root-spectrum-n6-computations.md`.

---

## OPEN-Q-048 🟢 Ultra-Log-Concavity for Tournament Independence Polynomials

**The theorem (proved):** If $I(\Omega(T),x)$ is real-rooted (proved universally only for $n \leq 8$; false universally from $n=9$ by THM-025), then $(\alpha_k/\binom{d}{k})_{k=0}^d$ is log-concave (ultra-log-concave), where $d = \alpha(\Omega(T))$.

**Proof:** Newton's inequalities for real-rooted polynomials with positive roots. Elementary symmetric polynomials $e_k(\rho_1,\ldots,\rho_d)$ satisfy Newton-Maclaurin: $(e_k/\binom{d}{k})^2 \geq (e_{k-1}/\binom{d}{k-1})(e_{k+1}/\binom{d}{k+1})$. Since $\alpha_k = \alpha_d \cdot e_{d-k}(\rho)$, ULC follows.

**Erdős context:** This is the tournament analog of the Heron-Rota-Welsh theorem (ULC for matroid Whitney numbers, proved by Adiprasito-Huh-Katz). Both prove ULC via underlying geometry (real-rootedness for tournaments, Hodge theory for matroids).

**Status:** PROVED conditional on real-rootedness. Computationally verified n=6..9.
**Priority:** 🟢 Interesting. Connects to the Huh-Katz matroid theory.
**Source:** oracle-2026-05-17-S1, computation `root_spectrum_fast.py`.

**NEW (oracle-2026-05-19-S1): UNCONDITIONAL proof of ULC at k=1 via Turán's theorem.**
For any tournament T: since bar_Omega(T) is K_{d+1}-free (max clique size = d = degree of I),
Turán's theorem gives alpha_2 <= (1-1/d)*alpha_1^2/2, which is exactly ULC at k=1:
   alpha_1^2 >= 2d/(d-1)*alpha_2.
No TRRT required. Equality iff I(Omega,x) = c*(x+rho)^d (all roots equal, Turán extremal).

**Also proved (conditionally on K4-free structure):** ULC at k=2 for complete tripartite
co-conflict graphs: (ab+bc+ca)^2 >= 3(a+b+c)*abc.
Proof: LHS-RHS = (1/2)[(ab-ac)^2+(ab-bc)^2+(ac-bc)^2] >= 0.
Verified: 0 violations in all n=9 samples (91/100 degree-3).
See `07-reflections/ulc-turan-unconditional-proof.md`.

---

## OPEN-Q-050 🟡 Unconditional ULC at k=2 via Kruskal-Katona

**Goal:** Prove $\alpha_2^2 \geq 3\alpha_1\alpha_3$ (ULC k=2, d=3) without assuming TRRT.

**Current status:**
- Proved for complete tripartite co-conflict graphs $K_{a,b,c}$ via the algebraic identity.
- Zero violations in n=9 random samples (91/100 degree-3 tournaments, 0 failures).
- Universal TRRT would have implied this via Newton's inequalities, but universal TRRT is refuted by THM-025.
- The "bad" counter-example ($K_4-e$ + isolated vertex, gives 25 < 30) does NOT occur in tournament conflict graphs.

**Approach:** Use the Kruskal-Katona shadow theorem for simplicial complexes, combined with the tournament-specific constraint that bar_Omega(T) arises from an actual tournament. The key step is showing that the $K_4$-free graphs that violate $\alpha_2^2 \geq 3\alpha_1\alpha_3$ cannot be co-conflict graphs of tournaments.

**Why hard:** The complement of a tournament conflict graph has special "tournament Ramsey" structure beyond just being $K_{d+1}$-free. Characterizing all graphs that arise as $\bar\Omega(T)$ is an open problem.

**Priority:** 🟡 Important. Would give the first unconditional ULC result beyond k=1.

---

## OPEN-Q-051 🟡 Interlacing Approach to Real-Rooted Subclasses

**Correction (opus-2026-05-29-S8):** Universal TRRT is false by THM-025, so this cannot prove a theorem for all tournaments as stated. The interlacing approach may still characterize or prove real-rooted subclasses.

**The proof strategy (computationally supported in tested subclasses):**
If for every cycle C* in Omega(T), I(Omega \ C*, x) interlaces I(Omega, x)
when deg(I_del) = deg(I_full) - 1, then real-rootedness follows by induction via Hermite-Biehler for the tournaments satisfying the hypotheses.

**The deletion-contraction:** I(Omega,x) = A(x) + x*B(x) where A = I(Omega\C*) and B = I(Omega-N[C*]).

**Computational evidence:** 444/444 verified at n=6 (stride 16 sampling), 0 failures.

**Why it's hard:** The proof needs to show B interlaces A for the specific structure of tournament conflict graphs. This is analogous to the Chudnovsky-Seymour claw-free proof but for non-claw-free graphs (n≥9).

**Connection:** For any subclass where Ω(T)'s independence complex is matroid/gammoid-like, Choe-Oxley-Sokal-Wagner stability may imply real-rootedness.

**Priority:** 🟡 IMPORTANT. Could characterize the real-rooted subclass or identify the THM-025 failure in the HB framework.
**Source:** oracle-2026-05-19-S1, `interlacing_investigation.py`.
See `07-reflections/interlacing-and-trrt-proof-strategy.md`.

**MAJOR UPDATE (oracle-2026-05-21-S1):** The Hermite-Biehler condition is MUCH more precisely established:
- Recursion I = A + xB VERIFIED: 5210 checks, 0 violations.
- B interlaces A when dA=dB+1: **3537/3537 = 100%, 0 failures at n=6,7.**
- No-HB-cycle cases: exactly d=2,alpha2=1 — proved real-rooted by Turán unconditionally.
- In the tested real-rooted regime, the HB route reduces to TWO lemmas: (A) existence of HB-cycle and (B) interlacing.
- Proof sketch for subclasses: induction on m, using Turán for base cases and HB for induction.
See `07-reflections/hermite-biehler-trrt-strategy.md`.

---

## OPEN-Q-052 🟡 Lemma A: Existence of HB-satisfying Cycle

For any tournament T with d≥2 and α₂≥2 (or d≥3), prove that there exists a cycle C* such that deg(I(Omega\\C*)) = deg(I(Omega-N[C*])) + 1.

Computationally: holds for ALL tested n=6,7 cases (except d=2,alpha2=1 which is handled by Turán).
Proof approach: if alpha2>=2 or d>=3, there are multiple maximum independent sets. A cycle C* NOT in all max sets satisfies the condition.

**Priority:** 🟡 IMPORTANT (one of two lemmas for the HB real-rootedness subclass program; universal TRRT is refuted by THM-025).

---

## OPEN-Q-053 🟡 Lemma B: B Interlaces A in Hermite-Biehler Recursion

Prove: for any tournament T and cycle C* with dA=deg(I(Omega\\C*)) = dB+1 = deg(I(Omega-N[C*]))+1, the polynomial I(Omega-N[C*],x) interlaces I(Omega\\C*,x).

Computationally: **3537/3537 = 100%, 0 failures at n=6,7.** Strongest computational evidence for any structural claim in this project.
Approach: multivariate stability, or direct interlacing via tournament Ramsey structure.

**Priority:** 🟡 IMPORTANT (other lemma for the HB real-rootedness subclass program; together with Lemma A it cannot imply universal TRRT because THM-025 refutes that statement).

**Update:** Extended to n=8 (107/107) and n=9 degree-3 (28/28). Cumulative: 3672 cases, 0 failures.
Key identity: B interlaces A iff A(-sigma)<=0 where sigma=root of B. This = I(Omega,-sigma)<=0
since B(-sigma)=0 and I=A+xB. So Lemma B is: independence polynomial of Omega is non-positive
at the root of the B-polynomial. This may be provable via Lee-Yang / Grace-Walsh-Szego theorem.

---

## OPEN-Q-049 🟢 Root Ratio as SC Detector

**Conjecture:** SC tournaments have the most asymmetric root ratio $\rho_2/\rho_1$ (minimum ratio) at each $n$.

**Evidence:** At n=6: SC min ratio = 0.00251 (H=45, α₁=20, α₂=1) < NS min = 0.00279 (H=43, α₁=19, α₂=1).

**Formula:** For $\alpha_2=1$ classes: ratio $= 1/\rho_1^2 \approx 4/\alpha_1^2$. SC tournaments maximize $\alpha_1$ (via SC Maximizer mechanism), hence minimize the ratio.

**Key insight:** SC asymmetry is measurable in the ROOT SPECTRUM. The SC blowup mechanism (anti-automorphism pairing of cycles) forces the polynomial toward the "maximally asymmetric" configuration.

**Status:** CONJECTURED, supported n=6 (exhaustive for SC, 2000 samples for n=7).
**Priority:** 🟢
**Source:** oracle-2026-05-17-S1.

## OPEN-Q-053 🔴 Prove HYP-1732: alpha2(Omega(T)) <= p*(m-p) for pair-partner C*

**Added:** opus-2026-05-22-S2

**Setup:** T tournament with d=alpha(Omega)=2, C* pair-partner from THM-311, p=#cycles disjoint from C*.

**Claim:** alpha2(Omega(T)) <= p*(m-p).

**Equivalences (all proved):**
- ⟺ B interlaces A in the Hermite-Biehler decomposition (Lemma B for d=2)
- ⟺ I(Omega, -1/p) <= 0 (via the identity A(-1/p)=I(-1/p), THM-313)
- ⟺ p lies between the two positive roots of I(Omega(T),x)

**Verified:** 1637 tests at n=7..11, 0 violations (opus-S2). **Strengthened (monad-compute-2026-06-06-S1):** 132,604,306 pair-partner tests over 291,788 distinct α(Ω)=2 tournaments at **n=7,8,9** (uniform random), 0 violations; both equivalent forms (combinatorial bound and quadratic I(Ω,−1/p)≤0) agree per test. **Min slack (bound−α₂)=0 ⟹ the bound is SHARP.** Caveat: the S1 run's n=8 layer ate the full time budget, so n≥10 was budget-skipped, not tested — uniform random at n≥10 almost never gives α=2, so n≥10 needs targeted low-cycle construction (prior n=10,11 coverage stands). Script `hyp1732_large_sample_monad_s1.py`.

**Proof status:** OPEN. Partial results:
- B-B pairs only occur between groups with disjoint portal sets (THM-315, proved).
- Key inequality: e_AB(b1)+e_AB(b2) <= p for each B-B pair (proved from K3-free).
- Full proof requires tournament-specific argument beyond K3-free structure.

**Note:** TRRT for d=2 follows from Turán-ULC WITHOUT this lemma. HYP-1732 would give an ADDITIONAL structural proof via HB induction.

## OPEN-Q-054 🟡 Lemma A for the UNIQUE max IS case (d>=3)

**Added:** opus-2026-05-22-S2

**Status:** THM-314 proves Lemma A for ALL non-unique max IS cases (all d>=2). Remaining gap: unique max IS at d>=3.

**Situation:** When S is the unique max IS of size d>=3: every C*∉S has d_A=d and d_B<=d-1 (Key Inequality). Whether d_B=d-1 depends on T[V\V(C*)] having enough disjoint cycles. Empirically: 0 failures at n=9..11.

**Proof approach:** Show that for SOME C*∉S, the sub-tournament T[V\V(C*)] supports an IS of size d-1 in Omega restricted to cycles disjoint from C*. For d=3 at n=9 (three disjoint triangles): equivalent to showing some 6-vertex sub-tournament has two vertex-disjoint odd cycles.


---

## OPEN-Q-055 🟡 Forbidden H-spectrum: Other universally forbidden H values beyond 7

**Added:** opus-2026-05-28-S5 (with THM-343 completion).

**Status:** THM-343 proves H(T) ≠ 7 for ALL tournaments. **H=21 — finite window now CLOSED (monad-compute-2026-06-04-S4, HYP-2200):** the HYP-2193 reduction (H=21 ⟹ a strong component with H=21 ⟹ c₃≤α₁≤10 ⟹ by Moon m≤12; THM-079 Part G killed m≤8) left only strong tournaments on m∈{9,10,11,12} with c₃≤10; these were **exhaustively enumerated (isomorph-free) and contain NO H=21** (min H = 75,125,225,375). So H(T)≠21 for all tournaments — {7,21} is the complete permanent H-gap set, modulo elevating the canon inputs to a formal THM-115 proof. (Even cleaner: the Busch lower bound p(7)=25>21, MISTAKE-053, gives H≥25 for every strong tournament on m≥7 directly.) H=63 is REFUTED as a universal gap: it is achieved at n=8.

**EXHAUSTIVE n=8 H-SPECTRUM (monad-compute-2026-06-04-S1, `h_spectrum_n8_exhaustive_monad.py`):** the complete census over all 2^28 = 268,435,456 labeled 8-vertex tournaments (census total verified = 2^28; all H odd). 320 distinct H values, range [1, 661]. **The only low odd gaps are {7, 21}** — every odd value in [23, 609] is achieved. H=35,39,49,63 all unlock at n=8 (counts 161280/188160/604800/80640). The remaining odd gaps {611,615,617,619,623,625,635,647,655} are high-end sparseness just below max H=661 (not permanent). This makes the n=8 forbidden set ∩[1,609] = {7,21} EXACT (previously only 100k sampling, HYP-1104), and exhaustively confirms H≠7, H≠21 at n=8 (upgrades the H=21 (8,1)/(6,2) cases from "strong n≥8 sampling" to exhaustion).

**H-UNLOCK TABLE n=3..9 — answers the "at what n does each value unlock?" sub-question (monad-compute-2026-06-04-S7, `h_unlock_table_monad_s7.py`):** for every odd H, `unlock(H)` = smallest n in {3..9} with some tournament achieving it, built from the EXHAUSTIVE per-level spectra (n=3..7 generated here, iso-class counts re-checked against A000568=2,4,12,56,456; n=8 from S1 census; n=9 from S6 iso-classes). **Unlock cascade** (distinct H / maxH / NEW-at-n): n=3 (2/3), n=4 (3/5, +1), n=5 (7/15, +4), n=6 (19/45, +12), n=7 (77/189, +58), n=8 (320/661, +243), n=9 (1520/3357, +1200). **27 transient gaps unlock** with explicit n: H=35,39 at **n=7**; H=63,107,119,149 and 161..187 (odd) at **n=8**; the nine n=8 high gaps {611,615,617,619,623,625,635,647,655} at **n=9**. *Precision fix to the S1 entry:* H=35,39,49 first appear at **n=7** (not n=8 — the S1 "unlock at n=8" referred to their n=8 census counts); only H=63 truly first unlocks at n=8. **Permanent-through-n=9 gaps**: 159 odd values ≤ maxH=3357, of which **LOW (≤609) = exactly {7,21}**; the other 157 are high-end sparseness ≥2883 just below the new max. Sampled n=10/n=11 cross-check: H=7,21 absent in both (consistent with permanent); 9/157 of the n=9 high gaps are already seen achieved by n≤11 sampling (transient sparseness, not permanent). Table saved at `05-knowledge/results/h_unlock_table_monad_s7.tsv` (one row per odd H ≤ 3357; blank = not achieved through n=9). No new HYP/THM minted (MISTAKE-053 discipline).

**ALL 157 n=9 HIGH GAPS UNLOCK AT n=10 (monad-compute-2026-06-04-S9, `h_high_gap_unlock_sampling_monad_s9.py`):** the 157 "permanent-through-n=9" HIGH gaps in [2883, 3355] (everything beyond {7,21}) were attacked with heavy bias-swept near-transitive sampling at n=10/11/12 (Held-Karp `H_count`; transitive base with each forward arc reversed w.p. p, p-grid calibrated so the achieved-H cloud sweeps the target window). **Result: all 157/157 are ACHIEVED at n=10** — every one is TRANSIENT, not permanent. The n=10 phase (167,600 samples, 9,365 in-window) hit all 157 by t=125s (~33k samples); a partial n=11 phase (20,800 samples) re-confirmed 157/157. H=7 and H=21 never appeared in any sample (consistent with permanent). This upgrades S7's "9/157 lower-bounded" to **157/157 transient**, so the n=9 high-end sparseness is a pure finite-level artifact and **{7,21} stand alone as the sole candidate permanent low gaps** (proved forbidden by THM-343/THM-079; {7,21} is the complete permanent H-gap set). Per-target table: `05-knowledge/results/h_high_gap_unlock_sampling_monad_s9.tsv` (all first_n_achieved=10). Sampling certifies achievability (concrete witnesses), never permanence. No new HYP/THM minted (MISTAKE-053 discipline).

**S652 speedup handoff (codex-2026-06-05-S652, HYP-2228):** before attempting a blind exhaustive `n=10` H-spectrum census, build a certified structured-witness menu.  THM-410 interval matchings give an additive low-`c3` ledger (`n=10`: `9496` matchings, `5538` with `c3<=10`), the general upset bitset identity handles near-transitive perturbations around a fixed order, and square/module substitutions give exact run-cover/macro-word H counts (`C3[C3]` has `H=3159` vs naive `81`).  This will not prove absence, but it can explain and certify large regions of the n=10 unlock cloud before a C/NumPy full A000568(10) node is spent.

**Evidence:**
- H=21: 0 occurrences at n≤7 (exhaustive as of S6). All four decompositions (10,0), (8,1,0), (6,2,0), (4,3,0) of α-vectors absent at n=6.
- H=63: absent at n≤7, but **achievable at n=8**. THM-344 (opus-S10) gives the exact n=8 census: exactly two n=8 isomorphism classes have H=63; both have |Aut|=1, score sequences (1,2,2,3,3,5,6,6) and (1,1,2,4,4,5,5,6), and Ω(T)=K31, hence H=I(K31,2)=63.
- S11 structural fingerprint: both H=63 classes are single-core. Every odd cycle contains one vertex; deleting it leaves the transitive tournament. The core signatures `1001100` and `1100110` have weighted count r=31. Complete-Ω class census n=3..8 has no r=3 or r=10; single-core signature search has no r=3 or r=10 through length 16.
- S12 projection-defect update: both H=63 classes are exact old-projection kills (delete the core vertex and Ω vanishes). A core-stratified complete-Ω census through n=8 still has no r=3 or r=10 in any core stratum, and the single-core target search now has no r=3 or r=10 through length 40.
- **SINGLE-CORE SIGNATURE GAP — RESOLVED (monad-compute-2026-06-04-S2, `single_core_signature_complete_monad_s2.py`):** the single-core odd-cycle count is `r(s)=Σ_{i<j, s_i=1, s_j=0} f(j-i-1)`, `f(0)=1, f(t)=2^{t-1}`, over bit strings `s` (core arc pattern relative to a transitive order). Stripping leading 0s / trailing 1s is `r`-invariant, and a canonical witness of length `L≥3` has `r≥2^{L-3}` (its first-1/last-0 pair). So every achievable `r∈(0,R]` has a witness of length `≤3+⌊log₂R⌋`; an exhaustive enumeration to that length therefore proves un/achievability for ALL lengths. Verdicts (complete to `R=2^17`): **r=3 (H=7) and r=10 (H=21) are PERMANENT single-core gaps** — unreachable at any length (witnesses would have length ≤6, all checked), upgrading S12's "absent through length 40" to a finite theorem. **r=31 (H=63) is reachable** (first length 7, matches THM-344's `1001100`/`1100110`). The single-core gap set is dense (~50%), so single-core complete-Ω is a strict sub-construction — it explains why H=63 unlocks this way while H=7/H=21 cannot, but does NOT by itself prove H=21 globally forbidden (that is HYP-1753/THM-079's job). NB also r=94 (H=189, THM-025's count) is a single-core gap.
- **SINGLE-CORE GAP-SET STRUCTURE (monad-compute-2026-06-04-S3, `single_core_gap_structure_monad_s3.py`, HYP-2199):** the single-core gap set `G={r : r≠r(s) for any string s}` — computed complete to R=2^20 — has **asymptotic density exactly 1/2** (dyadic-window densities converge monotonically to 50.0%; the gap set is PERSISTENT/INFINITE, NOT finite). Both `G`={3,6,10,14,17,20,21,24,27,28,29,33,…} and the achievable set {1,2,4,5,7,8,9,11,12,13,15,16,18,19,22,23,25,26,30,31,…} are **NOVEL to OEIS** (no match). **No simple closed form:** not a residue-class union (mod≤12), not Thue-Morse (gaps 50.1% odious — popcount-parity-independent), not a Beatty sequence (gap-differences span 1..12+), achievable-set not an additive semigroup (1+2=3∈G) nor doubling-closed; only the powers of two are guaranteed achievable (`1·0^k→2^{k-1}`). So single-core complete-Ω carries no arithmetic structure that would single out {7,21} — reinforcing that the GLOBAL {7,21} gap is HYP-1753/THM-079's job (all Ω shapes), not the single-core picture's.
- Pattern correction: the apparent sequence {7,21,63} = {7·3⁰,7·3¹,7·3²} is a finite-n mirage. The 7·3^k universal obstruction terminates at k=1.

**Sub-questions:**
- ~~Prove HYP-1753 (H≠21 for all n).~~ **FINITE WINDOW CLOSED computationally** (monad-compute-S4, HYP-2200): exhaustive strong c₃≤10 enumeration on m=9..12 finds no H=21 (min H 75/125/225/375); combined with THM-079 (m≤8), Moon, THM-029, H-multiplicativity this completes the H≠21 case analysis. Remaining: a theorist should confirm the reduction chain (and/or the Busch p(7)=25>21 bound) and elevate THM-115 from conjecture to theorem.
- Prove HYP-1755 (Strong Key Lemma: 3 pairwise-int 3-cycles force a 4th INSIDE their vertex union). [No longer needed for H≠21, but still of independent interest.]
- ~~Prove or refute the single-core signature gap: r_core(s) never equals 3 or 10.~~ **RESOLVED** (monad-compute-S2, above): proven for ALL lengths — r∈{3,10} unreachable, r=31 reachable.
- Explain structurally why the two THM-344 classes are the first complete-core unlocks for H=63 while H=7 (K3) and H=21 (K10) remain blocked.
- Decide whether projection-kill/near-kill defects are the right invariant for separating complete-Ω unlocks from non-real-root residues.
- Is the forbidden set finite? At what n does each forbidden value "unlock"?

**Tools:** SCC decomposition + Moon-Moser + Moon-Camion (as in THM-343 proof). Strong Key Lemma. Score sequence analysis. Independence-vector enumeration. THM-344 n=8 class census.

---

## OPEN-Q-056 🟡 Merged Bucket Transport Excess

**Added:** kind-pasteur-2026-05-29-S5

**Question:** After THM-345's forced parity constraints and THM-346's general bucket-balance law, what controls the excess transport above the parity lower bound?

For each Hamming layer `d`, THM-345 gives:

- bucket sizes `B_M`;
- row sums `B_M*C(m,d)`;
- symmetry of `W_d`;
- even diagonal;
- forced cross-outflow parity.

The actual cross-line mass is much larger than the parity minimum. Is that excess determined by spine/ribs/sea type, H-gradient position, bucket size, or a new invariant?

**Next steps:**
- Label `W_d(M,N)` entries by SC-SC / SC-NS / NS-NS.
- Compare excess over parity lower bound by H-gradient and principal-line distance.
- Test whether generic NS-NS sea entries are approximable from bucket sizes alone.
- Package normalized bucket transport as a Tournament TDA feature.

**Source:** THM-345, THM-346, INV-194, `04-computation/merged_bucket_constraints_s5.py`, `04-computation/tiling_quotient_bucket_balance_s5.py`.

**Files:** 04-computation/{thm343_complete_proof,h_spectrum_forbidden,forbidden_h_n7,h21_structure}_s5.py; `04-computation/h63_counterexample_audit_s8.py`; `04-computation/omega_extreme_fingerprints_s11.py`; `04-computation/projection_defect_bridge_s12.py`; `05-knowledge/results/omega_extreme_fingerprints_s11.out`; `05-knowledge/results/single_core_signature_targets_s11.out`; `05-knowledge/results/projection_defect_bridge_s12.out`

---

## OPEN-Q-057 🟢 Exact value of N* — the smallest N whose unit-distance maximum beats 3N

**Status:** N* ∈ [25, 28] (THM-431, sharpening THM-421's [17,32]). PROVEN floor N*≥25 (u(n)≤3n for all n≤24, via AMP arXiv:2412.11914 exact n≤21 + upper bounds u(22)≤61,u(23)≤66,u(24)≤72); PROVEN ceiling N*≤28 (realizable u(28)≥85>84). The dispatched n=21 campaign is itself SETTLED: **u(21)=57** (AMP, proven optimal; extremal graph = K₃□W₇, the unit-triangle × unit-wheel Cartesian product, 57=3·7+3·12).

**The sharp target is n=27 = 3³.** The best known construction *ties* exactly there: u(27) ≥ 81 = 3·27. The best-construction deficit u≥(n)−3n runs −6,−5,−4,−3,−2,**0**,+1 for n=22..28, closing to a clean tie at 27 before breaking through at 28.

**To settle N*, either:**
1. **Lower the ceiling:** find an exact-integer construction beating 3N at n∈{25,26,27} (is u(27)=81 or >81?). **It MUST be NON-PRODUCT** — see THM-433 below.
2. **Raise the floor:** prove an upper bound u(n)≤3n for n∈{25,26,27} (AMP's current upper bounds 78,84,90 exceed 3n=75,78,81 — they would need improvement).

**SHARPENED (THM-433, monad-explorer-2026-06-07-S1):** average degree is ADDITIVE under the Cartesian/Minkowski product, `avgdeg(G□H)=avgdeg(G)+avgdeg(H)`. Over the proven optima u(n) (n≤21, all factors of N≤42 are ≤N/2≤21, so this is EXACT) the product family caps at `P(N)≤3N` for **every N≤31**, ties only at **{27,30}**, and first strictly beats 3N only at **N=32** (W₁₆□K₂, 98>96). Since N*∈[25,28] sits strictly below 32, **the crossover graph is necessarily NON-PRODUCT (irreducible / Moser-lattice).** The tie at n=27=3³ is the **Cartesian cube K₃^□3** (avgdeg 2+2+2=6); `G₉□K₃`, `G₁₀□K₃` give the ties at 27,30. ⟹ **No product can give 82 at n=27** (additivity caps it at 81); the suggestion in the old handoff to seek "a product config with 82 edges" is impossible — only a non-product (Moser) config can decide u(27)>81. Bonus: u(32)≥98 (was 97). Files: THM-433, `04-computation/unit_distance_product_cap_s1.py` (+`.out`), reflection `average-degree-is-additive-and-the-crossover-is-irreducible-s1.md`.

**Note (THM-431-C):** the √7 Eisenstein family (THM-421's construction lane) is the WRONG family — it only beats 3N at n≈39 (disk) / 32 (anneal). The first crossing is boundary-dominated (THM-431-C) AND irreducible/non-product (THM-433) — two independent reasons it evades the "structured" families. Any attempt to lower the ceiling should use the Moser lattice, not √7 disks or products.

**UPDATE (THM-432, monad-explorer-S711) — the n=27 tie IS the Hamming graph H(3,3).** The Erdős product (Minkowski sum) `K₃□K₃□K₃ = K₃^□3 = H(3,3)`: 27 points, **6-REGULAR**, exactly 81=3·27 unit distances (verified exact in ℚ(√3)). The mysterious "3³ tie" is literally the 3-fold product of unit triangles; it ties (not beats) because a product of triangles is forced 6-regular and `6=κ`. Product criterion: `G□H` beats 3N ⟺ `ρ(G)+ρ(H)>3` (avg degrees sum > κ). Census (proven u(a)): smallest product TIE at N=27 & 30, smallest product BEAT at N=32 (K₂□G₁₆, U=98>96). **Since N*∈[25,28]<32, N* is NOT a product — it is an irregular rigid blob** (consistent with AMP's Moser-lattice extremals). The best product per n is *tight with the global optimum exactly at n=27* and loses by only 1–3 elsewhere. ⟹ strong structural evidence (not proof) that **u(27)=81, hence N*=28** (HYP-2299). Concrete next probe: is the u(28)≥85=81+4 crosser `H(3,3)+1` (a 28th point unit-distant from 4 of its vertices)? — pure products are futile below N=32. Also independently reproduced AMP's *proven* u(21)=57 extremal as exact W₇□K₃. See `04-computation/unit_distance_product_crossover_monad_s711.py`, reflection `07-reflections/symmetry-saturates-irregularity-violates-the-hamming-tie-s711.md`.

**CUBE ANGLE-RIGIDITY (THM-437, monad-explorer-2026-06-07-S6) — the cube cannot be tuned past 81.** The obvious route to `u(27)>81` is to *tune the three rotation angles* of `H(3,3)=K₃^□3` so it gains accidental (non-product) unit distances. **PROVED impossible:** the 81 product edges are angle-independent; any *extra* unit distance needs a sum of triangle-edge unit vectors (one per differing factor) of length 1, and the complete solution set of the 3-factor condition `cos u+cos w+cos(u−w)=−1` is exactly `{t₂≡0}∪{t₃≡0}∪{t₂−t₃≡0} (mod 60°)` — each a **collision locus** (two factors align in the Eisenstein lattice ⟹ two of the 27 points coincide). So for every angle choice: 27 distinct points ⟹ exactly 81 unit distances. The 3N tie at n=27 is *angle-rigid*, not a generic-angle artifact. This **closes the "just tune the cube" idea** and complements THM-432/433: even non-product perturbations of the cube are stuck at 81, so `N*`'s extremal graph (if ≤27) must be a genuinely irregular blob — neither product nor tuned cube. **Scope:** rules out the cube family only; does NOT prove `u(27)=81` (AMP upper bound still 90). Files: `01-canon/theorems/THM-437-cube-angle-rigidity-at-81.md`, `04-computation/unit_distance_cube_angle_rigidity_monad_s6.py` (LEM-A/B/C, 0 rogues), reflection `the-cube-tie-is-angle-rigid-accidental-edges-collide-s6.md`.

**HARBORTH CORRECTION (monad-explorer-S6).** The S4 entry's "a 27-cell triangular/penny blob gives ≈78 (deficit −3)" is **wrong by 15**: the exact max triangular-(Eisenstein-)lattice patch is `⌊3n−√(12n−3)⌋` = **63** at n=27 (deficit −18), confirmed by an exact greedy patch search matching Harborth at every n=22..28 (49,52,55,57,60,63,65). The flat triangular patch is far from optimal at these n; the route to 81 is the *3-layer cube* (H(3,3)), not a flat patch + O(1) surplus. So the S4 "concrete residue" (triangular patch + 4 off-lattice → 82) is mis-scaled. (Numbers in `05-knowledge/results/unit_distance_augment_cube_monad_s6_partAB.out`.)

**Source:** THM-431, THM-432, THM-437, HYP-2298, HYP-2299, monad-explorer-2026-06-07-S710/S711/S6; AMP arXiv:2412.11914; `04-computation/unit_distance_3n_floor_sharpen_s710.py`, `04-computation/unit_distance_product_crossover_monad_s711.py`, `04-computation/unit_distance_cube_angle_rigidity_monad_s6.py`.

**Sharpening note (S2, HYP-2300):** the PRODUCT family first beats 3N only at N=32 (S1 THM-432), while N*(true) ∈ [25,28] is irreducible. The gap `32 − N*(true) ≥ 4` is the "irreducibility premium" — the unit-distance face of the integrality gap χ>χ_f (opus-S699g Vitali wall). The Cartesian-product trichotomy (HYP-2300) proves products are structurally a UD-only lower-bound device (avgdeg AMPLIFIES under []; LRC's lonely density DEGRADES, HN's χ NEUTRAL), so the crossover graph at N* MUST be irreducible — no product search can find it. Pinning N*(true) makes the premium exact.

**RESONANT-PRODUCT UPDATE (THM-493, monad-explorer-2026-06-13) — the "non-product" crosser IS a product at the RESONANT angle; the bonus is the crossing.** The Moser lattice `L_t=ℤ[ζ₆]⊕ω_t·ℤ[ζ₆]` is literally the Minkowski product of the triangular lattice with a copy rotated by the **Moser angle** `ω_t`. At a generic angle the product is Cartesian (THM-433, Δ=0); at `ω_t` the **transverse unit vectors** of THM-434 appear as extra diagonal edges, giving the EXACT count `U(G⊞_t H)=e(G)|H|+|G|e(H)+Δ_t`, `Δ_t=½Σ_{N(α)=t}m_α(G)m_α(H)` (a correlation of the factors' `√t`-displacement spectra). **Constructive `u(28)≥85`:** `W₇⊞₃R` (Eisenstein rosette × unit rhombus, Moser angle t=3) has `48+35+Δ₃=48+35+2=85>84` on 28 points — the SAME product graph has only 83 (=P(28)) at a generic angle, and the `Δ₃=2` transverse edges ARE the entire crossing `83<84<85`. So THM-433's "non-product crossover" = "product + the non-additive transverse bonus." **Why 27 holds (sharper than THM-433/437):** `27=3³` forces a size-3 factor, and the densest 3-point UDG `K₃` is `√t`-FREE for every t ⟹ zero resonance bonus (this re-explains THM-437's cube angle-rigidity). `28=4·7` is the first composite whose factorization (rhombus×W₇) gives a `√3`-bearing *and* edge-dense factor pair. A curated exact two-factor resonant search finds NO beat at n=25,26,27 (best 72,61,75 < 3n; K₃^□3=81 ties with bonus 0) — evidence for `u(27)=81, N*=28`. To settle: an *upper* bound `u(27)≤81` — and THM-493 says the obstruction at `3³` is **arithmetic** (no edge-dense `√t`-factor at size 3), not merely geometric. Files: `01-canon/theorems/THM-493-resonant-product-decomposition-unifies-thm433-thm434.md`, `04-computation/resonant_product_{bonus,Nstar_search}_monad.py` (+`.out`), reflection `the-resonance-bonus-is-the-crossing-and-27-is-bonus-hostile.md`.

**FREE-PATCH CONFIRMATION + DECISIVE t=4 CONTROL (HYP-2461, monad-explorer-2026-06-13, agent-mathworker).** Independent of THM-493's curated 2-factor search, I ran the SAME exact-integer annealing densest-patch search (the one reproducing Engel's table on `L_3`, here re-deriving `u(28)=85>84`) over the WHOLE lattice (free, non-product) across the bridge family `L_t` (`t=2,3,4,5,13,21,31,49`; THM-434 unit counts 12..30), `n=21..30`, every count exact-recounted. Result: a clean **tie-vs-crossing dichotomy.** (i) The `81`-tie at 27 is reached by EVERY transverse-bearing lattice (`t=3,4,13,21,31,49`) and NEVER beaten; transverse-FREE `t=2,5` cap at `78` (can't even build the tie). So the tie is the carrier-robust 6-regular `H(3,3)` and **no free patch in any bridge carrier reaches 82** — strong non-product evidence for `u(27)=81, N*=28`. (ii) **The decisive control:** `t=4` (`√−15`, 18 units, rosette angle **29.0°, geometrically CLOSER to the 30° bisector than Moser's 33.6°**, same 6 transverse vectors) ties `81@27` but caps at `83<84@28` — it does NOT cross. So the `u(28)=85` crossing is NOT a "good-angle band"; it is the SPECIFIC ARITHMETIC of `√−11` (THM-493's `Δ₃`), singular to `t=3` among all tested carriers. (iii) unit-vector COUNT is a red herring (30-unit `t=49` ≈ 24-unit `t=13` ≈ 18-unit Moser; all within 1-3 of Engel). Bonus: `L_13`=the `√13` layer reaches `u(21)=57` (AMP's optimum lives there) ⟹ one `L_t` family holds BOTH the n=21 optimum (`t=13`) and the crossover engine (`t=3`). Open follow-up: does the exact-30° lattice `ℤ[ζ₁₂]` (NOT an `L_t`) cross at 28, or is `√−11` arithmetically singular? Files: `04-computation/unit_distance_bridge_lattice_family_monad.py`, `..._bridge_t4_probe_monad.py` (+`.out`s), reflection `the-unit-distance-tie-is-carrier-robust-the-crossing-is-resonant.md`. Complements/credits THM-493.

**BISECTOR HANDOFF RESOLVED + RATIONAL-COSINE CHARACTERIZATION (THM-494, monad-explorer-2026-06-13).** Answered HYP-2461's next-explorer question: does the exact-30° bisector `ℤ[ζ₁₂]` (the geometrically perfect interleaving the `L_t` family brackets — t=3→33.6°, t=4→29.0° — but can never hit) cross `3N` at `n=28`? **NO.** Exact-integer densest-patch search (same engine, calibrated to the exact triangular maximum): `ℤ[ζ₁₂]` caps at **78@27** (cannot even build the 81 tie) and **83@28** (does not cross) — *bit-for-bit the transverse-free `t=2,5` profile*. **Reframe (the durable content):** the bisector fails not because 30° is a "bad angle in a band" but because it is **OFF THE RESONANCE LADDER ENTIRELY.** THM-494 (PROVED): a glued pair of triangular lattices at rotation `ω=e^{iθ}` carries a *transverse* unit vector `α(1−ω)` iff `|1−ω|²=2−2cosθ=1/N(α)`, i.e. iff `cosθ=(2t−1)/2t` — **the Moser ladder is exactly the rational-cosine rotations.** `cos30°=√3/2` is irrational (`|1−ζ₁₂|²=2−√3≈0.268`, bracketed between t=3's 1/3 and t=4's 1/4 but off-ladder), so by **Kronecker** `ℤ[ζ₁₂]` has exactly 12 unit vectors (two 30°-hexagons, zero transverse). **Niven's theorem** makes it a clean dichotomy: rational *angle* (cyclotomic `ζ₁₂`@30°, `ζ₈`@45°) and rational *cosine* (the ladder) are disjoint except at the degenerate 60°. ⟹ the crossing lives on rational-cosine/irrational-angle rotations; the geometrically perfect rational-angle bisector is arithmetically barren. **Third independent confirmation of THM-493:** the bisector hits exactly `P(28)=83` (generic product cap) and falls short by exactly `Δ₃=2`. The crossing gates are now nested: (1) on-ladder = rational cosine [bisector fails here], (2) Loeschian `r_E>0` = transverse exists [`t=2,5` fail here], (3) crossing-resonant at n=28 [only `t=3` passes]. Open: characterize gate 3 — is `t=3`/`√−11` unique, or merely first? (Heegner class-number-1 texture unused.) Files: `01-canon/theorems/THM-494-transverse-resonance-is-rational-cosine-the-bisector-is-off-ladder.md`, `04-computation/unit_distance_zeta12_bisector_monad.py` (+`.out`), reflection `the-perfect-bisector-is-off-the-ladder-rational-cosine-not-rational-angle.md`.

**GATE-3 RESOLVED — THE CROSSING NORM IS THE SMALL FACTOR'S CHORD (THM-495, monad-explorer-2026-06-13).** Answers THM-494's open "is `t=3`/`√−11` unique or merely first?" The THM-493 bonus `Δ_t(G,H)=½Σ_{N(α)=t}m_α(G)m_α(H)` is nonzero **only if `t` is a shared chord norm of both factors** (THM-495(A), PROVED corollary) — so admissible `t` ⊆ ChordSpec(small factor). **(i) FORCED-UNIQUE at 28:** `28=4·7`, the dense 4-factor is the rhombus `R=K₄−e` with `ChordSpec(R)={1,3}`; its only non-unit chord is `√3`, so `Δ_t(R,W₇)=0` for ALL `t≠3` (exact scan t=2..59: lone survivor t=3, Δ₃=2). `t=3` is not "first" — it is the ONLY admissible norm. **(ii) DOMINANT everywhere:** across the whole 2-factor triangular family `n=24..49`, `Δ₃≥Δ₄≥Δ₇` in every case, because `√3` is the triangular lattice's second-nearest-neighbour (the most abundant non-unit chord). **(iii) UNIFIES THM-437 ⊕ THM-493 combinatorially:** `27=3³` routes every factorization through a size-3 factor = the **chord-free triangle** (`ChordSpec(K₃)={1}`), so it gets ZERO resonance bonus and can only tie 81 — re-deriving the cube angle-rigidity by counting chords, no `cos u+cos w+cos(u−w)=−1` calculus, and pinning the reason: 3 is prime and its optimal factor is chord-free. **The whole `27→28`/`N*` boundary = chord-free vs chord-bearing smallest factor.** Strong new combinatorial support for `N*=28` (HYP-2299). PROVED (A,B,C) + VERIFIED in-family (D); does NOT prove u(27)=81 (AMP bound still 90). Open (HYP-2466): prove `m_α`-domination for all dense patches; bridge to free-patch crossings. Files: `01-canon/theorems/THM-495-resonant-crossing-norm-is-the-small-factor-chord-spectrum.md`, `04-computation/resonant_crossing_chord_spectrum_monad.py` (+`_partB.py`, +`.out`s), reflection `the-crossing-norm-is-the-small-factors-chord-the-triangle-is-chord-free.md`. New: HYP-2466.
**CAPACITY ATLAS UPDATE (HYP-2467, codex-2026-06-13).** Exhausting every connected triangular-lattice factor patch through size `9` gives an exact small-factor resonance ledger for THM-493's bonus and complements THM-495's chord-spectrum theorem. In this carrier class `27=3*9` maxes at `75<81`: `K3` is edge-dense but resonance-free, while the resonance-bearing 3-point paths reach only `69/70` against all 9-patches. `28=4*7` reaches the known `85>84` crosser with generic `83` plus `Delta_3=2`; `30=5*6` ties and `32=4*8` crosses. This does not prove `u(27)<=81`, but it turns the next proof step into a compression theorem: any 82-edge 27-point Moser patch must evade connected small-factor resonance capacity. Files: `04-computation/unit_distance_resonance_capacity_atlas_codex.py`, result `.out`, HYP-2467, OPEN-Q-085, T807.

**LATTICE-PERFECTION GATE — WHY 27 NOT 28, AND THE FIRST IMPERFECT SIZE IS 9 (THM-496, HYP-2468, monad-explorer-2026-06-13).** Orthogonal axis to THM-495's chord-spectrum: define size `k` **lattice-perfect** iff `Harborth(k)=u(k)` (the triangular lattice attains the planar unit-distance max). **PROVED/VERIFIED (all connected patches k≤9, 77359 at k=9):** `Harborth(k)=u(k)` for `k≤8`, and **`k=9` is the FIRST imperfect size** (`u(9)=18>16=Harborth(9)`). Since resonant-product factors live in `ℤ[ζ₆]`, a resonant product matches the generic Cartesian cap only when every factor size is `≤8`. ⟹ the **resonant cap at `27=3·9` is `75`, NOT the `81` generic tie** — resonance strictly HURTS at 27 (the `81` is the generic/off-lattice cube; triple-confirmed with THM-495 & codex HYP-2467). The complete 2-factor gate is the conjunction: (i) lattice-perfect factorization (parts ≤8), (ii) chord-bearing factor (size ≥4), (iii) `Δ_t>gap(n)=3n−P_gen(n)`. `n=24,25` pass (i)+(ii) but fail (iii) (gap 6,5 ≫ exhaustive max bonus 2); `n=26,27` fail (i)+(ii) (13,9 imperfect); **`n=28=4·7` is the FIRST to clear all three** (LP, rhombus carries √3, gap=1<Δ₃=2 → 83+2=85). The exhaustive 2-factor resonant maxima `68,72,65,75` at `n=24..27` are now EXACT (upgrading THM-493's curated search). **Deep link:** `u(9)=18` needs `K₃□K₃` at a GENERIC angle (lattice collapses to 16) — the smallest "product needs an irrational angle" / integrality-premium instance — and the 27-optimum is `K₃□G₉` (generic cube), so the imperfection at 9 **propagates multiplicatively** to 27 (HYP-2468). Does NOT prove `u(27)=81` (AMP bound 90); it pins the `27→28` boundary to the first lattice-imperfect size. Files: `01-canon/theorems/THM-496-lattice-perfection-gate-resonant-crossover.md`, `04-computation/lattice_perfection_gate_monad.py` (+`.out`), reflection `the-lattice-perfection-gate-nine-is-the-first-imperfect-size.md`.

**LATTICE-LANE CONFIRMATION (HYP-2301, monad-explorer-2026-06-07-S4) — the [28,32] gap from a SECOND, independent family.** A systematic exact-integer densest-patch sweep over SIX single-norm lattice families {penny t=1, knight t=5, √7, √13, grid t=25, grid t=65} (anneal calibrated to the repo's √7=97@32, every patch exact-recounted) finds **NO single-norm lattice beats 3N at N≤28**; the earliest is **√7 at N=32** (exactly where products bottom out — the convergence), √13 at 33, while the *higher-degree* knight (deg8), t=25, t=65 cross *much* later (>60). Governing law: a **degree–radius tension** `N_cross ∝ ρ·t·(deg/(deg−6))²` (radius² × a degree-excess penalty that is ∞ at κ=6, 16 at deg8, 4 at deg12), minimized uniquely by √7 (deg12 at minimal norm 7) — so the "32" rung is the genuine min over ALL single-norm lattices, not a √7 artifact, and the "irreducibility premium" [28,32] equals the "cost of regularity" from the lattice side too. **Punchline (corrected):** the degree–radius tension IS the 2-D kissing bound — a rank-2 lattice cannot carry deg>6 at radius 1. Engel's u(28)=85 is NOT a 2-D lattice patch (triangular gives 65, best √7 gives 83); it lives in the **rank-4 Moser ring M_L=ℤ[ζ₆,ω₃]** whose non-torsion unit ω₃=(5+i√11)/6 (cos 5/6) packs **18 unit vectors at radius 1**, escaping the tension. So [28,32] = the cost of staying rank-2; the right hunt for u(27)>81 is a dense M_L patch (exact ℚ(√3,√11) machinery in `unit_distance_moser_lattice_u21_monad_s4.py`), NOT a denser 2-D lattice and NOT a product. **VERIFIED — ceiling now self-contained:** a densest-patch search run DIRECTLY in M_L (graph-BFS + anneal in ℤ⁴ with the 18 unit-vector offsets, exact |z|²=1 recount over ℚ(√3,√11)) reproduces Engel's ENTIRE deficit table from scratch — u(M_L)=60,64,68,72,76,**81 (tie 27)**,**85 (beats 3N at 28)**,89,93 for n=22..30 — so THM-431's previously CITED ceiling N*≤28 is now backed by explicit exact-integer coordinates found here. Files: `04-computation/unit_distance_3n_crossover_{families,focus,moser_crossover}_s4.py`, reflection `the-3N-crossover-is-won-by-the-densest-layer-plus-surplus-not-a-high-degree-layer-s4.md`.

**SUB-QUESTION (1) ANSWERED — NEGATIVELY (THM-437, monad-explorer-2026-06-07-S5).** "Is the u(28)≥85 crosser literally `H(3,3)+1`?" → **NO**, for the generic realization. Exact `ℚ(√3)` circumcircle enumeration over the faithful generic K₃^□3 (27 pts, 81 edges, 6-regular; triangles rotated by Pythagorean angles 3-4-5 & 5-12-13): the ONLY unit circles through ≥3 vertices are the 27 Eisenstein hexagons, **each centered ON an existing vertex** ⟹ no off-vertex point is unit-distant from ≥3 vertices ⟹ any added 28th point has degree ≤2 ⟹ `H(3,3)+1pt ≤ 83 < 85`. Not even a one-point perturbation of the product — genuinely irreducible. (Generic-realization caveat; special-angle is a separate finite check.) **Also new — the product-defect profile** δ(N)=u(N)−bestproduct: δ=0 (product-optimal) at {6,8,9,12,21}, δ>0 (irreducible) at {4,10,14,15,16,18,20}, all δ≤2 ⟹ irreducibility is the RULE below threshold but always by ≤2 edges; **N* = first N where this O(1) surplus also lifts α=2u/N past κ=6** (tangent at 27 ⟹ generic prediction N*=28). α superadditive over multiplication (=Erdős bound); principal line α(3^j)=2j tangent to κ=6 at 27=3³. Files: THM-437; HYP-2304; `04-computation/unit_distance_product_defect_monad_s5.py` (+`.out`); reflection `the-product-defect-profiles-irreducibility-s5.md`.

## OPEN-Q-058 🟡 The Tournament Barba Problem (n ≡ 1 mod 4): prove max det(I+S) = 2(n-1)^((n-1)/2)

**Status:** OPEN — but the LOWER (construction) half is now PROVED: THM-475 (claudebox-2026-06-11-S1), the DRT flag construction. For every n ≡ 1 mod 4 with a DRT on n−2 (all orders under the skew-Hadamard conjecture; unconditionally for n−2 ∈ Paley/doubling-tower/GF(27) orders), Flag(DRT(n−2)) = DRT + two stacked apexes attains 2(n−1)^((n−1)/2) with EXACTLY the conjectured spectrum x(x²+n−2)^((n−3)/2)(x²+2n−3) — verified exactly at n = 9, 13, 17, 25, 29; at n=9 the flag char poly equals the unique char poly of all 216 exhaustive maximizer classes. Remaining open: the UPPER bound (no tournament beats the flag). Strong evidence (mac-mini-2026-06-10-S2, HYP-2389). Exhaustively exact at n=5 (32) and n=9 (8192 = 2^13, 216 classes ALL sharing spectrum x(x²+7)³(x²+15)); hill-climb HIT the conjectured 5971968 = 2·12⁶ at n=13 in <1s with exactly the predicted spectrum x(x²+11)⁵(x²+23), >1M restarts found nothing higher. The conjectured extremal spectrum is two-level: flat base n−2 with multiplicity (n−3)/2 plus ONE excited pair at 2n−3. The n ≡ 2 mod 4 analog without skew-EW shows the same (n−3)-base + (2n−3)-excited shape (n=6: (y−3)²(y−9)). This is the missing congruence class of the maximal-determinant theory for skew-type matrices: n ≡ 3 mod 4 is Reid–Brown/DRT (THM-472), n ≡ 0 mod 4 is skew-Hadamard, n ≡ 2 mod 4 is the Armario/Greaves–Suda skew E-W theory (2n−3 square condition), and n ≡ 1 mod 4 appears genuinely untreated (literature + OEIS negative, 2026-06-11). Proof routes: integrality/Galois constraints on the char poly of S + the trace identity Σμ² = n(n−1)/2; or a Greaves–Suda-style spectral rigidity argument. A proof would be a publishable companion to Klanderman et al. LAA 707 (2025).

## OPEN-Q-059 🟡 Tournament Ky Fan: replace Fan's magnitude order by an arbitrary tournament

**Status:** OPEN, literature-confirmed empty (2026-06-11 search). Ky Fan's lemma counts ALTERNATING simplices — monotone label chains with alternating signs, i.e. antidirected paths in the TRANSITIVE tournament on label magnitudes — and guarantees an odd number of them. The tournament-side parity results that exist (Rédei = all-forward type; Forcade 1973: every orientation type has odd count when n = 2^k; El Sahili–Abi Aad 2020: antidirected Hamiltonian paths ≡ 2 mod 4 at even order, proving Grünbaum's conjecture) have no Fan-style topological formulation. QUESTION: is there a Z₂-equivariant/simplicial statement in which the linear order of Fan's labels is replaced by an arbitrary tournament T, with the alternating-simplex count controlled by an invariant of T (H(T)? the orientation-type parities?)? A positive answer would make Rédei/Forcade theorems shadows of a Borsuk–Ulam-type theorem. Entry points: Prescott–Su's constructive proof (path-following = the project's transfer-matrix style), the bistellar-move proof (arXiv:2308.07103), the s690 double-cover reading of tournaments (odd sections of the pair double cover), and THM-474 (tilings = switching classes — the gauge in which the base path P₀ IS Fan's linear order). Related new data: x ↦ Sx is an odd tangent field whose hairy-ball singularity is the Pfaffian vector w, kept off the sum-zero sphere by Rédei parity (HYP-2398).

## OPEN-Q-060 🟢 The odd Mallows–Sloane partner: what does A049313 count, the way A002854 counts Euler graphs?

**Status:** OPEN — sharpened by THM-479 (claudebox-2026-06-11-S2): the count splits as A049313(n) = N_odd(n) + N_lev(n) (odd-order branch + even-level branch, both separately integral for n ≥ 3, N_lev = 0 at odd n; values in switching_classes_level_burnside_cbx2.out; neither branch in OEIS). Any "second incarnation" must respect this 2-adic branch split — graphs:Euler graphs :: tournaments:(odd-branch object ⊔ even-level object)? Note Babai–Cameron Lemma 3.1: the even-level branch is symmetry WITHOUT fixed member tournaments, so the partner object cannot be "a distinguished member per class" at even n (Mallows–Sloane's even-n non-bijectivity, verified quote, is the same wall). (Originally flagged by the two-graphs literature sweep, 2026-06-11.) Mallows–Sloane: #two-graphs = #switching classes of graphs = #EULER GRAPHS (A002854 — which equals the project's even-graph metagraph node counts V(E_n)). The tournament analog A049313 (#switching classes of tournaments up to iso = #oriented two-graphs: 1,1,2,2,6,12,79 for n=2..8, Babai–Cameron Thm 7.2, summed over LEVEL permutations — constant 2-adic valuation across cycles) has NO known second combinatorial incarnation. Find the natural class of "odd directed objects" equinumerous with it. The project owns the natural toolkit: THM-474 (tilings = labeled switching classes), the even-graph metagraph E_n, and the level-permutation 2-adic seam. A bijective answer would complete the even/odd duality square: graphs:even-graphs :: tournaments:???.

## OPEN-Q-061 🟡 The extremal [72,36,16] code as a tournament-gauge problem

**Status:** OPEN (claudebox-2026-06-11-S4, HYP-2415). One of the most famous open problems in
coding theory — does an extremal Type II self-dual [72,36,16] binary code exist? (Sloane 1973;
$\$$-history; still open 2026.) THM-481's eQR tournament-gauge ladder C(I+S(Paley_q)) is
EXTREMAL Type II at q = 7, 23, 31, 47 (lengths 8, 24, 32, 48; minimum distances 4, 8, 8, 12 =
4⌊n/24⌋+4, all verified exactly) and FIRST FAILS at **q = 71**: eQR(72) has d = 12 < extremal
16. Since order 72 ≡ 8 (mod 16), the tournament gauge C(I+S(H)) of EVERY skew-Hadamard matrix
H of order 72 is a Type II [72,36] code (M. Hall §17.3). **Sufficient route (not an
equivalence):** if any order-72 skew-Hadamard / doubly-regular-tournament-switching has gauge
minimum distance 16, the famous code exists. Paley (the highest-symmetry tournament) gives only
12; there are very many other skew-Hadamard matrices of order 72 (Đoković–Kotsireas catalogues).
**Program:** compute (or bound below, via partial-distance / coset-leader methods) the gauge
minimum distance of known order-72 skew-Hadamard classes; characterize which tournament
spectral feature of H lifts d from 12 to 16. A sharp tournament-theoretic handle on a famous
coding open problem, and the natural continuation of the THM-480/481/482 gauge line. Repo
bridge: THM-484 (24 = involution modulus; the eQR ladder is extremal exactly while the Gleason
extremal d = 4⌊n/24⌋+4 stays at the Golay/√-ramped value, jumping to 16 at the 3rd multiple of
24 where Paley loses it). Task t-0120.

## OPEN-Q-062 🟢 A Bombieri–Vinogradov level-of-distribution for the LRC multiplier orbits

**Status:** OPEN (claudebox-2026-06-11-S5, HYP-2416 part B). The Elliott–Halberstam exponent θ
measures how deep into the modulus range q ≤ x^θ the primes equidistribute among residue classes
(Bombieri–Vinogradov: θ=1/2 unconditional; EH: θ→1). The LRC window lemma (S625) is the same shape
for the multiplier orbits (ℤ/m)*: a good multiplier (a residue avoiding every runner's width-2
danger band) survives once the shell m is deep enough. QUESTION: formulate the LRC analogue
precisely — for a random/typical speed set, control on AVERAGE over shells m ≤ M the discrepancy
between the danger-band-avoidance count and its expectation, and identify the "level" M = M(n) up to
which this holds. What is the LRC analogue of the √-barrier θ=1/2 (conjecturally the gap between the
easy M>1/(2n) and the optimal 2/(2n−1)), and is there a large-sieve/Bombieri-type proof at "θ=1/2"?
A positive answer would import the bounded-gaps technology (GPY/Maynard–Tao multidimensional sieve)
into the LRC/covering frontier. Repo bridge: HYP-2416 (the dictionary), THM-406/S561 (ρ = the sieve),
the S625 window lemma, THM-415 (the optimal 2/(2n−1)). Honest: this is the right QUESTION; a proof
needs analytic large-sieve input the repo does not yet have. Task t-0121.
**Status:** OPEN (flagged by the two-graphs literature sweep, 2026-06-11). Mallows–Sloane: #two-graphs = #switching classes of graphs = #EULER GRAPHS (A002854 — which equals the project's even-graph metagraph node counts V(E_n)). The tournament analog A049313 (#switching classes of tournaments up to iso = #oriented two-graphs: 1,1,2,2,6,12,79 for n=2..8, Babai–Cameron Thm 7.2, summed over LEVEL permutations — constant 2-adic valuation across cycles) has NO known second combinatorial incarnation. Find the natural class of "odd directed objects" equinumerous with it. The project owns the natural toolkit: THM-474 (tilings = labeled switching classes), the even-graph metagraph E_n, and the level-permutation 2-adic seam. A bijective answer would complete the even/odd duality square: graphs:even-graphs :: tournaments:???.

## OPEN-Q-064 🟡 Random pentagonal interior-zero theorem and zero-Lyapunov sign laws

**Status:** OPEN (codex-2026-06-11-P1). Let `D_eps(q)=1+sum eps_g q^g` over generalized pentagonal exponents `g=k(3k+-1)/2`. Euler's signs factor as `prod(1-q^n)`, so `1/D` is the partition generating function and has zero ordinary Lyapunov exponent. Random signs on the same support experimentally have positive finite-window reciprocal growth. Prove (or refute): a random pentagonal sign denominator almost surely has a zero in `|q|<1`, giving positive reciprocal Lyapunov exponent. Secondary classification problem: which deterministic sign laws on pentagonal support have zero reciprocal Lyapunov exponent? The all-plus control has low finite-window slope, so uniqueness of Euler is NOT safe. Entry points: Jensen formula for random analytic functions, Rouché/small-ball estimates on two radii, and finite truncation root certification. Files: HYP-2424, T783, `04-computation/pentagonal_lyapunov_code72_codex.py`.

## OPEN-Q-063 🟡 Tutte/matroid support gate for the extremal Type II [72,36,16] target

**Status:** OPEN (codex-2026-06-11-P2). The length-72 scalar Gleason enumerator is healthy (`A_16=249849`, `5-(72,16,78)` minimum design), and Type II formal scalar positivity persists through the stored `24..240` ladder. Use Greene's theorem to recast the code existence problem as binary self-dual matroid support realization at a Tutte specialization. Build a leakage diagnostic: first forbidden low dual weight, first design-incidence failure, first neighborhood obstruction, first automorphism-forced contradiction. The goal is a support-building Tournament Analysis whose vertices are construction moves and whose observable is `(low-weight suppression, design/neighborhood compatibility)`, expected to be nontransitive where scalar cancellation and realizability trade off. Files: HYP-2425, HYP-2429, HYP-2430, T781, `04-computation/cancellation_gate_atlas_codex.py`.

## OPEN-Q-065 🟢 Dirichlet-character version of the Euler-product ghost atlas

**Status:** OPEN (codex-2026-06-11-P3). The ordinary `q`-product atlas separates exponent schedules, Witt ghosts, and coefficients for eta/primes/Mobius/Liouville/random signs. Build the Dirichlet analogue `prod_p(1-chi(p)p^{-s})` for true characters and random completely multiplicative signs, then compare carriers: Dirichlet zero pressure, ordinary coefficient leakage, ghost irregularity, and partial-sum cancellation. The first target is a two-observable Tournament Analysis that is no longer transitive. Files: HYP-2431, HYP-2432, T782, `04-computation/euler_product_ghost_atlas_codex.py`.

## OPEN-Q-066 🟡 The 72 support bridge between Nebe lattices and binary Type II codes

**Status:** OPEN (codex-2026-06-11-P4). The scalar theta gate and scalar Gleason gate both pass at dimension/length 72: the lattice row kills `q^1,q^2,q^3` and starts with `6218175600 q^4`, while the code row kills weights `4,8,12` and starts with `249849 y^16`. Nebe's extremal 72-dimensional even-unimodular lattice exists; the binary `[72,36,16]` code remains open. Find the retained support bridge or obstruction: lattice polarizations, frame data, Z4/code lifts, binary matroids, skew-Hadamard gauges, or the `5-(72,16,78)` design incidence layer. Files: HYP-2433, HYP-2434, HYP-2435, T784, `04-computation/theta_code_lattice_gate_codex.py`.

## OPEN-Q-067 🟡 Complete or kill the order-5 branch of the extremal [72,36,16] code

**Status:** OPEN (codex-2026-06-11-P5). The order-5 fixed projection has been reduced to a tiny exact gate: for automorphism type `5-(14,2)`, the projected fixed code must be `e8+e8` with the two fixed coordinates split across the summands; the `d16+` branch is excluded because every marked pair lies in a tetrad and lifts to weight `12`. Thus the fourteen 5-cycles split into two heptads with Fano-plane tetrads, producing exactly `14` fixed minimum words and `49967` moving minimum-word orbits. The next problem is the nonfixed `F_16` component: enumerate or obstruct Hermitian self-dual `[14,7]` candidates compatible with the split-heptad fixed boundary, binary minimum distance `>=16`, and the residual `5-(72,16,78)` design ledger. Files: HYP-2439, HYP-2440, HYP-2441, T785, `04-computation/order5_fixed_projection_72_codex.py`.

## OPEN-Q-068 🟡 Prove the LRC14 Q27 resource bound beyond one stranger

**Status:** OPEN (codex-2026-06-11-P6, HYP-2444/HYP-2438). The one-stranger family `S(r)=7*{1..12} union {r}` is now closed by the fibered band-1 lattice `Q27={d*m:d|14,m<=27}`: all 936 primitive rows have a Q27 witness, and the two rows whose first plain witness is `q=41` are caught at `q=91`. The residue mechanism is explicit: the 7-core covers 8/9 classes of `(Z/27)^*/+-`, misses `+-10`, and every plain q<=27 shell blocker also has `r mod 13=0`. The open problem is to lift this from one stranger to arbitrary primitive multiple-of-14 configurations: prove that blocking Q27 consumes independent resources across shell-27 classes, low clocks such as 13, divisor fibers `d in {1,2,7,14}`, and B' safe-component gaps, so that 13 runners cannot block all Q27 and B'(any). First computational target: two-stranger rows with a resource-vector output, constrained to keep low divisor clocks covered; the naive pair of one-stranger blockers over `7*{1..11}` is too easy because all 28 such pairs have a q=12 witness.

## OPEN-Q-069 🟡 Transfer Church's diagonal Frobenius support gate to LRC14 and the [72,36,16] support problem

**Status:** OPEN (codex-2026-06-12, HYP-2445). Church's product-quotient counterexamples show the scalar/support split in geometric form: Shioda supersingularity is too coarse, while diagonal symmetric forms on every partial Frobenius twist force curve descent or finite exceptional types. Formalize the shared support-gate lemma: scalar quotient `Q`, retained channel `S`, and descent/exception rule `D`. Test two transfers: for LRC14, can Q27 blockers be forced to spend independent resources or descend to Bprime/owner-deletion exceptions; for `[72,36,16]`, can the minimum-word `5-(72,16,78)` support ledger use the `D7` index `78` and `D6/A4` index `91` as incidence-arithmetic probes? Files: HYP-2445, T789, `04-computation/product_quotient_support_gate_atlas_codex.py`, reflection `shioda-product-quotient-obstructions-and-support-gates.md`.

## OPEN-Q-070 🟡 Build the irreducible-prime certificate tournament

**Status:** OPEN (codex-2026-06-12, HYP-2448; extends HYP-2447). Formalize the finite/infinite tournament suggested by Bunyakovsky/Buniakowski plus the Singh/Cohn/Iravanian reverse certificates. Vertices should be certificate states, not just polynomials: fixed divisor, local residue obstructions, least Singh/Cohn value depth, trace-subset survivor profile, and Newton/non-Archimedean support data. Edges orient toward smaller unresolved factorization ambiguity after normalizing degree and fixed divisor, with richer retained support as tie channel. First tasks: replace the floating real-trace scout by exact algebraic trace lattices; build `C(f;X)` for a larger polynomial family and measure edge flips as `X` grows; translate the same carrier to LRC14 Q27 resource vectors and to `[72,36,16]` support/matroid/design construction moves. Files: HYP-2448, HYP-2447, T792, `04-computation/irreducible_prime_carrier_tournament_codex.py`, reflection `irreducible-prime-carriers-and-certificate-tournaments.md`.

## OPEN-Q-071 - Build the marked coefficient-row irreducibility tournament

**Status:** OPEN (codex-2026-06-12, HYP-2449; extends HYP-2447/HYP-2448). The coefficient-sign tiling is a genuine fixed-path tournament carrier, but the finite scout shows bare unmarked tournaments and sign vectors are too coarse for irreducibility. Formalize a marked coefficient-row state `R(f;P,X)` consisting of skip-row signs, coefficient magnitudes, local zero-prime residues, p-adic valuation/Newton-slope data for primes `P`, Cohn base/evaluation addresses, and Singh value-depth certificates up to `X`. First tasks: implement exact Newton-row tournaments for Eisenstein/Dumas/Perron criteria; measure edge flips as primes and evaluation depth are added; compare Cohn prime rows against low-Omega composite rows with identical sign tournaments; transfer the fixed-divisor row detector to LRC14 Q27 resource rows. Files: HYP-2449, T793, `04-computation/coefficient_tiling_prime_irreducible_codex.py`, reflection `coefficient-tiling-and-prime-irreducible-addresses.md`.
## OPEN-Q-072 🟡 Classify irreducible coefficient-magnitude slices in the tiling quotient

**Status:** OPEN (codex-2026-06-12, HYP-2450; extends HYP-2448). The coefficient-tiling quotient maps fixed-path tournaments to count profiles `c_d` and centered magnetizations `A_d=2c_d-(N-d)`. Cohn gives one rigorous lane: a positive-degree prime base-value of the diagonal-count profile certifies irreducibility of the digit polynomial. The open problem is the magnetization lane: characterize magnitude vectors `(|A_d|)` whose sign slices are forced irreducible, forced reducible, or have bounded factor patterns. The pilot at `N=6` finds the parity-minimum slice `(1,0,1,0,1)` has only 8 distinct polynomials and all are irreducible, while the full fixed-path quotient has 91/120 profiles with hidden `H` variation. Transfer target: attach the lost fiber data to LRC14 Q27 resource ledgers and to `[72,36,16]` support/matroid/design realization. Files: HYP-2450, HYP-2448, T794, `04-computation/coefficient_tiling_prime_bridge_codex.py`, reflection `coefficient-tilings-and-the-prime-irreducible-bridge.md`.

## OPEN-Q-073 - Build split-survivor ledgers for polynomial rows and LRC14 resources

**Status:** OPEN (codex-2026-06-12, HYP-2451; extends HYP-2449/HYP-2450). Reducibility is a hidden convolution lift of the coefficient row, so the live state should record which degree-split rectangles survive after each local gate. First tasks: add Newton/valuation certificates to the `38` degree-4 irreducibles with no small mod-p blocker through `31`; extend split-survivor signatures to degree `5` and `6` with cached finite-field factorizations; add Singh-depth/Cohn-depth only for rows that survive residue and valuation gates; transfer the same survivor ledger to LRC14 Q27 denominator/resource fibers, replacing scalar `q blocked` with surviving local lift obligations. Files: HYP-2451, HYP-2450, HYP-2449, T795, `04-computation/convolution_lift_irreducibility_carrier_codex.py`, reflection `convolution-lift-split-survivors-and-hidden-factor-grids.md`.
## OPEN-Q-074 🟡 Build bounded integer convolution-lift obstructions beyond degree 5

**Status:** OPEN (codex-2026-06-12, HYP-2452; extends HYP-2451/HYP-2450). Reducibility can be encoded as an integral hidden tiling problem: find nontrivial factor coefficient rows `b_i,c_j` whose multiplication grid has diagonal sums `a_k=sum_{i+j=k} b_i c_j`. The HYP-2452 pilot gives an exact integer oracle for primitive degree `<=5`, with zero mismatches against Sympy on `3856` degree-4 rows and `2016` degree-5 rows, complementing HYP-2451's residue/valuation split-survivor carrier. The open problem is to push this beyond the linear/quadratic-factor range without falling back to a black-box factorizer: encode bounded degree splits as SAT/ILP/SMT feasibility, add Newton-slope boundary constraints for sparse/multivariate rows, and use Singh-style low-`Omega(f(m))` factor-capture witnesses as pruning. Transfer target: treat LRC14 blocker ledgers and `[72,36,16]` weight-enumerator coefficients as boundary totals whose hidden support/incidence lifts must exist. Files: HYP-2452, HYP-2451, HYP-2450, T796, `04-computation/convolution_factor_capture_tiling_codex.py`, reflection `convolution-factor-capture-and-hidden-coefficient-tilings.md`.

## OPEN-Q-075 - Build moment-lift resource ledgers for LRC14 shells

**Status:** OPEN (codex-2026-06-12, HYP-2453; extends HYP-2443/HYP-2444/HYP-2452). The triangular-tower computation reframes the user's two towers as moment-balanced shell splits: `A_n` is the square-shell first-moment split and `B_n` is the triangular-shell second-moment split. The first two moments give exact integer starts `n^2` and `2n^2+n`; higher moments require a fractional address with leading term `(p-1)(p-2)/(12p)`. Addendum: A covers every positive integer, while B only covers alternating triangular shells; whole-equation side-aligned containment is the Pell family `T_n=2T_m`, and the exact whole-side equality `B_3.L=A_4.R=[21,24]` is unique. The open problem is to transfer this to LRC14 by enriching Q27 blocker rows with moment/resource data: blocked unit twists, owner supports, divisor fibers, raw moment defects, and the fractional or fiber address needed to lift a scalar blocked shell into an actual support proof. First tasks: prove the higher-moment expansion to more terms, extend the floor-sqrt/Beatty classifier into a useful Q27 address ledger, and compare AP/V*/2AP plus HYP-2444 one-stranger residuals under the new ledger. Files: HYP-2453, HYP-2444, HYP-2452, T797, `04-computation/triangular_tower_moment_bridge_codex.py`, `04-computation/triangular_tower_overlap_families_codex.py`, reflection `triangular-towers-moment-lifts-and-fractional-addresses.md`.
## OPEN-Q-076 🟡 -- PARTIALLY RESOLVED
**Prove the triangular power-center bracket and finish the 78/90 support transfer**

**Status:** PARTIALLY RESOLVED (codex-2026-06-12, HYP-2454; addendum to HYP-2453). The user's ordinary and square towers are exact interval power balances with centers `2T_n` and `4T_n`. The finite scout still verifies that for `3<=p<=8` and `n<=40`, the positive root of `D_p(C,n)=0` lies between `2pT_n` and `2pT_n+1`, but the new Faulhaber packet now proves the exact odd-moment reformulation
`c^p = 2 * sum_{r odd} binom(p,r) c^(p-r) S_r(n)` with `c=a+n` and `u=n(n+1)`, derives
`a_p(n)=p*n^2+(p-1)*n+(p-1)(p-2)/(12p)-(p-1)(p-2)(2p^2-4p-1)/(180 p^3 n(n+1))+O(n^-4)`,
and records the square-pyramidal cuboid identity `6*(1^2+...+n^2)=n(n+1)(2n+1)`. It also shows that the `n=1` live higher-power balance factors are irreducible through `p<=20`, so the tower split stops exactly where Bernoulli/Faulhaber corrections begin. The remaining open part is the global theorem: prove or refute the bracket for all `p>=3`; solve the Pell-style endpoint boundary families controlling overlaps between the first square-shell partition and the second square-balance tower; and turn the special row `Q_L(3)=[21,22,23,24]`, with ordinary shadows `90=S1(4)` and `78=C(13,2)`, into an actual support-ledger constraint for LRC14 and the `[72,36,16]` `5-(72,16,78)` minimum-design problem. Files: HYP-2454, HYP-2453, T798, `04-computation/triangular_power_balance_towers_codex.py`, `04-computation/triangular_power_faulhaber_asymptotic_codex.py`, reflections `triangular-power-balance-towers-and-additive-square-bridges.md` and `faulhaber-odd-moments-and-square-pyramidal-cuboids.md`.

## OPEN-Q-077 🟡 Build a common hidden-lift feasibility engine across irreducibility, LRC, unit distance, and code72

**Status:** OPEN (codex-2026-06-13, HYP-2455; extends HYP-2452/HYP-2444/OPEN-Q-057/HYP-2454). Recent work says the scalar boundary total is not the proof object: polynomial coefficients need convolution factor grids, LRC q-blocking needs runner/Pisano/divisor/owner support, unit-distance products are reducible baselines before Moser-irreducible fibers, and `[72,36,16]` weight enumerators need support-design incidence. Build a shared lift-feasibility data model with boundary totals, candidate hidden cells, local gates, surviving allocations, and proof owners. First tasks: degree-6 bounded ILP/SAT for HYP-2452, multi-stranger LRC allocation ledgers beyond one-stranger Q27, product-reducible versus Moser-irreducible `N=27/28` unit-distance fibers, and a binary incidence-lift encoding for the `[72,36,16]` `78/90` support address. Files: HYP-2455, T799, `04-computation/boundary_lift_analogy_atlas_codex.py`, reflection `boundary-lift-irreducibility-transfer.md`.

## OPEN-Q-078 - Build a Beatty-Pell style Q27 address ledger for LRC14

**Status:** OPEN (codex-2026-06-13, HYP-2456; concrete instance of HYP-2455; extends HYP-2241/HYP-2443/HYP-2444/HYP-2453). The triangular crossover word now has an exact hidden-address normal form: a Beatty shell address `d_m`, a Pell/carry remainder `r_m`, and state inequalities whose rare equality walls are Pell atoms. Build the analogous LRC14 ledger for Q27 blockers. For each candidate row and denominator, record `(q, shell class, unit quotient class, divisor fiber, owner support, carry residue, endpoint/boundary atom, opening or deletion target)` rather than only "q blocked." First tasks: run this on AP/Vstar/2AP and the HYP-2444 one-stranger family; measure whether visible strict/wall/open status becomes pure after adding owner/carry/private-deletion fields; compare the remaining boundary atoms to the triangular `LR/RL` zero-density wall grammar. Files: HYP-2456, HYP-2455, HYP-2453, HYP-2241, HYP-2443, HYP-2444, T800, `04-computation/triangular_tower_beatty_pell_decomposition_codex.py`, reflection `beatty-pell-crossover-word-and-lrc-address-ledgers.md`.

## OPEN-Q-079 - Prove the Faulhaber anchor expansion/bracket and port odd-wall ledgers to LRC14

**Status:** OPEN (codex-2026-06-13, HYP-2457; sharpens HYP-2454 and complements HYP-2456). The midpoint defect for the power-balance anchor is exactly `D_p(c,n)=c^p-2*sum_{r odd} binom(p,r)c^(p-r)S_r(n)`, so only odd Faulhaber moments survive. The stored computation verifies the formal expansion `c=p*n(n+1)+alpha_p+beta_p/(n(n+1))+gamma_p/(n(n+1))^2+...`, with all displayed corrections divisible by `(p-1)(p-2)` and hence exact recovery of the p=1/p=2 towers. First tasks: prove a uniform fixed-`p` remainder after `gamma_p`; use it, or a sharper direct inequality, to prove HYP-2454's bracket `D_p(p*n(n+1),n)<0<D_p(p*n(n+1)+1,n)` for every `p>=3`; compare the p=2 square-pyramidal cuboid packing against higher simplex/cuboid carriers; and build an LRC14 analogue where odd walls, owner support, shell-27 class, divisor fiber, carry residue, and endpoint atom replace scalar "q blocked" status. Files: HYP-2457, HYP-2454, HYP-2456, T801, `04-computation/triangular_faulhaber_anchor_expansion_codex.py`, reflection `faulhaber-anchors-square-pyramids-and-bernoulli-addresses.md`.

## OPEN-Q-080 - Build an odd-moment compatibility lift analogous to OCF alpha packets

**Status:** OPEN (codex-2026-06-13, HYP-2458; extends HYP-2457/HYP-2456 and OCF). HYP-2457 isolates the odd Faulhaber anchor expansion, but OCF warns that odd atom counts are not the full object: `H(T)=I(Omega(T),2)` needs compatible packets `alpha_k` of vertex-disjoint odd cycles. Build an explicit finite compatibility lift whose one-particle shadow is the odd moment list and whose packet terms record coexistence of shell, carry, endpoint, owner-support, and support-design atoms. First targets: add odd-atom compatibility fields to the HYP-2456 Beatty/Pell side states, run them against LRC14 Q27 AP/Vstar/2AP and HYP-2444 one-stranger rows, and test whether code72 `78/90` support packets behave more like OCF `alpha_k` than like scalar moments. Files: HYP-2458, HYP-2457, T802, `04-computation/faulhaber_odd_moment_ocf_bridge_codex.py`, reflection `faulhaber-odd-moments-and-ocf-cycle-packets.md`.

## OPEN-Q-081 - Build a parity-typed Q27 ledger for LRC14

**Status:** OPEN (codex-2026-06-13, HYP-2459; extends HYP-2458/HYP-2444/HYP-2443 and the complement-Walsh line). The projector rule is exact: midpoint anti-symmetrization keeps odd Faulhaber channels, while tournament converse keeps even Walsh channels for invariant scalars. The open LRC task is to type every Q27 ledger field as `even_scalar`, `odd_marked`, `transported`, or `compatibility_packet`. First targets: AP/Vstar/2AP and HYP-2444 one-stranger rows; split source/sink or start/end fields into sum and difference before quotienting; then test whether remaining primitive rows either get a strict witness, descend to the known wall atoms, or expose an odd owner/carry/deletion opening. Files: HYP-2459, HYP-2458, T803, `04-computation/parity_projector_channel_atlas_codex.py`, reflection `parity-projectors-and-even-odd-channel-gates.md`.

## OPEN-Q-082 - Prove LRC14 Q27 hard-resource independence

**Status:** OPEN (codex-2026-06-13, HYP-2463; sharpens OPEN-Q-081/HYP-2444/HYP-2443). The complete hard-replacement hull is exact: all `77520` rows formed by replacing `k-1` core speeds in `7*{1,...,12}` by `k` hard residues from `{260,351,442,611,702,793,962,1053}` have Q27 witnesses, with only ten plain-shell misses and only two non-original residuals caught at `q=30,34`. Prove the compression theorem suggested by this: any primitive LRC14 row with no Q27 witness can be parity-typed and compressed to this hard-replacement hull without losing blockedness, unless it opens a low clock, divisor-fiber witness, AP/Vstar/2AP descent, or odd owner/Bprime deletion. Files: HYP-2463, T804, `04-computation/lrc14_parity_typed_q27_ledger_codex.py`, result `lrc14_parity_typed_q27_ledger_codex.out`, reflection `lrc14-hard-resources-do-not-stack.md`.

## OPEN-Q-083 - Prove the LRC14 two-stranger compression-resource lemma

**Status:** OPEN (codex-2026-06-13, HYP-2464; refines OPEN-Q-082/HYP-2463). The bounded two-stranger stress test deletes one runner from `7*{1,...,12}` and adds two distinct non-core speeds up to `13*84`, scanning `6,868,368` primitive rows. Only `877` block every plain shell `q<=27`, and all `877` have Q27 witnesses. The residuals are broader than the hard-residue hull (`636/877` use no old hard residue), but every residual has at least one added `13`-clock speed, no residual deletes `7,21,49`, and the late cases are divisor fibers including `91=7*13` and `161=7*23`. Prove the unbounded version: any primitive row blocking the plain shell must either compress to these resource coordinates or open a low clock, divisor-fiber witness, AP/Vstar/2AP descent, or odd owner/Bprime channel. Files: HYP-2464, T805, `04-computation/lrc14_two_stranger_compression_stress_codex.py`, result `lrc14_two_stranger_compression_stress_codex.out`, reflection `lrc14-two-stranger-compression-stress.md`.
## OPEN-Q-084 - Force any LRC14 Q27 blocker below the nine-core threshold or into descent

**Status:** OPEN (codex-2026-06-13, HYP-2465; strengthened by HYP-2470). HYP-2465 proves an exact bounded near-core lemma: in the carry window `1..1092`, no primitive replacement row retaining at least nine speeds of `7*{1,...,12}` can block Q27. HYP-2470 decides the first below-nine-core boundary in corrected form: retaining exactly eight core speeds forces either Q27 or a plain witness `q<=41`. The next proof task is to show every possible LRC14 row either normalizes into this carry window and therefore must delete at least five core speeds if it has no Q27/no plain-41 witness, or else opens a known side channel. Then analyze the below-eight-core regime: prove it forces low clocks, divisor-fiber witnesses, AP/Vstar/2AP descent, owner-private/Bprime deletion, or a support-load contradiction. Files: HYP-2465, HYP-2470, T806, T809, `04-computation/lrc14_near_core_q27_setcover_codex.py`, `04-computation/lrc14_eight_core_q27_setcover_codex.py`, results `lrc14_near_core_q27_setcover_codex.out`, `lrc14_eight_core_q27_setcover_codex.out`.

## OPEN-Q-085 - Prove the small-factor resonance capacity gap for unit distances

**Status:** OPEN (codex-2026-06-13, HYP-2467; refines OPEN-Q-057/THM-493/THM-495/HYP-2462/HYP-2466). Exact connected-factor atlas through size `9` proves the resonant product carrier separates the `27 -> 28` gate: `27=3*9` maxes at `75<81`, while `28=4*7` reaches `85>84`. The size-3 obstruction is now explicit: `K3` is edge-dense but has no non-degenerate norm-`t` displacement, while resonance-bearing 3-point paths lose too much generic edge budget. Prove this capacity lemma analytically, then lift it from connected triangular factors to arbitrary dense rank-4 Moser patches; any 82-edge 27-point counterexample must evade this compression. Files: HYP-2467, T807, `04-computation/unit_distance_resonance_capacity_atlas_codex.py`, result `unit_distance_resonance_capacity_atlas_codex.out`, reflection `unit-distance-resonance-capacity-and-the-27-28-gate.md`.
## OPEN-Q-086 - Prove the Church-style LRC14 descent theorem

**Status:** OPEN (codex-2026-06-13, HYP-2469; upgraded by HYP-2470). Church's arXiv:2508.14876 gives the proof template: scalar quotient is not enough; a retained side channel on every twist/fiber plus finite exceptions or strict descent carries the obstruction. LRC14 transfer: raw plain-shell blocking is the scalar shadow, while Q27 obligations, 13-clock debt, deleted-core address, shell-27 class, divisor fiber, owner/Bprime support, and support-load data form the retained channel. Certified finite anchors now cover one-stranger, hard-stack, two-stranger residual, near-core `|D|<=3`, and the corrected eight-core shell-41 boundary. Prove the remaining descent theorem: any primitive row with no Q27 and no plain `q<=41` witness either normalizes into those finite atlases, opens a named exception, or strictly lowers a resource rank. First tasks: below-eight-core typed MILP/set-cover for `|D|>=5`, outside-window Bprime/divisor/carry normalizer, and formal exception catalogue. Files: HYP-2469, HYP-2470, T808, T809, `04-computation/lrc14_church_frobenius_descent_codex.py`, `04-computation/lrc14_eight_core_q27_setcover_codex.py`, results `lrc14_church_frobenius_descent_codex.out`, `lrc14_eight_core_q27_setcover_codex.out`.

## OPEN-Q-087 - Prove the LRC14 below-eight-core or outside-window descent

**Status:** OPEN (codex-2026-06-13, HYP-2470; strengthened by HYP-2471 and corrected by THM-497). HYP-2470 decides the first below-nine-core finite boundary in corrected form: in the carry window, retaining at least eight core speeds forces either a Q27 witness or a plain witness `q<=41`. HYP-2471 adds that the two Q27-only four-deletion exceptions both die over Q31, preserving the divisor/fiber-ladder explanation. THM-497 adds the crucial warning that this is a near-core theorem, not a global plain-shell ceiling: dilated-band cardinality permits covers, and kps1 scouts produce non-dominant rows blocking all plain shells through `41`. The remaining proof task is therefore typed below-eight-core/outside-window descent: any normalized row with no Q27/no useful fiber witness and no retained structural opening must delete at least five core speeds, leave `1..1092`, violate replacement normalization, or land in a named AP/Vstar/2AP/owner-Bprime/low-clock exception. Files: HYP-2470, HYP-2471, THM-497, T809, T812, T813, `04-computation/lrc14_eight_core_q27_setcover_codex.py`, `04-computation/lrc14_q31_exception_probe_codex.py`, `04-computation/lrc14_resource_climb_kps1.py`, result files, reflections `lrc14-eight-core-exceptions-open-at-shell41.md`, `lrc14-q31-fiber-repair-for-eight-core-exceptions.md`, and `lrc14-covering-cardinality-permits-structure-forbids-kps1.md`.

## OPEN-Q-088 - Prove the LRC14 ramified irreducibility-transfer portal

**Status:** OPEN (codex-2026-06-13, HYP-2480; sharpens HYP-2470 and imports HYP-2451/HYP-2452 tactics). HYP-2480 shows that the two Q27-feasible four-deletion packets in HYP-2470 share an Eisenstein/Newton-like valuation pattern: `12/13` speeds are divisible by `7`, there is exactly one primitive non-7 escape, and that escape is divisible by `13` (`936` or `1066`). Both packets open at the missing plain-shell layer (`q=33` or `q=31`) and also have Bprime/positive safe-measure certificates. Prove this as a lemma rather than a census: a Q27-feasible near-core packet with 7-ideal occupancy plus a single 13-clock primitive escape must open at `q in {31,33,41}` or Bprime. Parallel tasks: extract dual/Farkas set-cover certificates from Q27 infeasibility, formalize the factor-capture obligation budget, and use Cohn/Perron dominance to normalize outside-window speeds. Files: HYP-2480, T810, `04-computation/irreducibility_tricks_lrc_transfer_codex.py`, result `irreducibility_tricks_lrc_transfer_codex.out`, reflection `irreducibility-tricks-and-lrc14-ramified-local-gates.md`.

## OPEN-Q-094 ✅ WITHDRAWN (already closed by HYP-2271) — The Tournament Frobenius Number

**Status:** WITHDRAWN by mac-mini-2026-06-11-S2 on the same day it was posed. The
question "is the forbidden-H set finite?" is ALREADY ANSWERED in canon: HYP-2271 /
HYP-2180 prove the H-spectrum is a **co-finite multiplicative numerical semigroup with
exactly 2 gaps {7,21} (genus 2, Frobenius number 21)**. Mechanism: H multiplicative
over strong components (Moon) + strong-minimum minH_strong(m) = 3,5,9,15,25,45,75,…
= Busch (2006) (= Moon's upper bound), strictly increasing, so minH_strong(m) ≥ 25 > 21
for all m ≥ 7; combined with exhaustive m≤6 strong spectra, {7,21} are the only
permanent gaps. (The formula m²−5m+9 is a 4-point coincidence that fails at m=7 —
MISTAKE in 01-canon/MISTAKES.md; true value 25 not 23.) The lesson logged: SCOUR
CANON before posing — this is exactly the genus-2 numerical-semigroup result already
abstracted in polarized-delta-fields-band-gaps-and-numerical-semigroups-s699.md. The
genuinely-new lens (additive Goldbach vs multiplicative Euler, the s=2 segment bridge)
is HYP-2424, not an open question.

## OPEN-Q-095 🟢 Is there a tournament invariant that is a square exactly when an alternating pairing is present (the Ш-square analog beyond THM-174)?

**Status:** OPEN/exploratory (mac-mini-2026-06-11-S2, HYP-2420). THM-174 gives det(skew S)=Pf² (even n) — the literal "alternating ⟹ square" mechanism shared with |Ш(E)| being a perfect square. BSD's Poonen–Stoll shows the square can degrade to 2×□ when the pairing is only antisymmetric (a ±1 defect), mirrored heuristically by THM-442's H²−Pf²=8Q correction. QUESTION: is there a NATURAL finite abelian group / pairing attached to a tournament (a "Ш-analog") whose order is exactly det(I+2A) or related to Q, and whose alternating-vs-antisymmetric type is controlled by a tournament parity (SC/transpose-self/blue-black)? If so, the project would have a combinatorial model of the Cassels–Tate square phenomenon with an explicit, computable defect. Entry: the cokernel of I+2A or of S as a finite abelian group (Smith normal form), its induced pairing, and whether SNF squares track det=Pf². Cf. Klanderman et al. (Smith normal form of skew D-optimal designs).

## OPEN-Q-097 🔴 inf L(S) > 0 over the LRC(14) dilated-AP cores — the archimedean |T|≥3 sinc-lattice bound

**Status:** OPEN, the C'(14)/LRC(14) prize. **[UPDATE kind-pasteur-2026-06-16-S7, THM-522/MISTAKE-075/HYP-2561]: TWO new exact levers + an inf correction + a compactness reframe.** (i) `L` is SCALE-INVARIANT `L(cS)=L(S)` (τ↦cτ measure-preserving) and QUANTIZED `L∈(1/(14·lcm S))·ℤ` ⟹ `L>0 ⟹ L≥1/(14·lcm S)`. (ii) So inf `L>0` ⟺ the L-minimizers have BOUNDED lcm (compactness): quantization makes any bounded-lcm family automatic, scale-invariance kills dilation, THM-518 stranger-decoupling kills one-entry-→∞; the open piece is configs with lcm→∞ at bounded shape. (iii) The inf was OVERESTIMATED: `inf L ≤ 1/1260 ≈ 0.000794` (NOT 0.0052), at the minimal single-element perturbation `12→36` of the tight AP = `{1,2,…,11,13,36}` — the prior search restricted to multiple-of-14 strangers and missed the SPORADIC-tight perturbations (the tight locus includes `{1..11,13,24}`, not just the AP). NEW PROGRAM: classify the (conjecturally finite, bounded-lcm) tight locus; then quantization+compactness give `inf L>0` with constant `1/1260` (HYP-2561). The Abel/Bedert-level-bound program below is the COMPLEMENTARY analytic route. ──────── (original, mac-mini-2026-06-14-S1, THM-503, HYP-2521/2522): Reduction chain (THM-398/501): LRC(14) ⟸ C'(14) ⟸ inf_S L(S) > 0 over primitive multiple-of-14 S, where L(S) = (6/7)¹³ + Σ over 7-primitive exact additive relations of (6/7)^{13−|T|}(−1)^{|T|} Π s(t_v), s(t)=sin(πt/7)/(πt). NEW structure (THM-503): (a) only relations with all coefficients coprime to 7 contribute (s(7j)=0); (b) the |T|=2 corrections are absolutely convergent (|P|≤g²/3v_av_b) and the almost-Sidon class (pairwise mass < 36/49) is PROVED loose; (c) L is the ARCHIMEDEAN singular integral — β_p(S)=L(S) for every prime p, so positivity is NOT a product-of-local-densities statement (HYP-2503 corrected). The remaining open core: inf L ≈ 0.0237 is attained at the dilated-AP cores d·{1,…,12}∪{r} (HYP-2521), whose suppression is driven entirely by the |T|≥3 7-primitive relations (e.g. 1·7−2·14+1·21=0), a CONDITIONALLY convergent multidimensional sinc-lattice sum (this open core is now sharpened by the concurrent THM-504: the cancellation is the cross-level (−1)^{|T|} alternation, not within-level). ~~PROVE: |Σ_{|T|≥3} corrections| < (6/7)¹³ uniformly~~ — **this literal target is FALSE** (mac-mini-2026-06-15-S5, THM-518). The per-level masses Λ_k = (6/7)^{13−k}Σ_{|T|=k}∏s **GROW**: Λ₂≈+0.11, Λ₃≈−0.55, Λ₄≈+1.17, so Σ_{|T|≥3}≈+0.62 ≫ (6/7)¹³=0.135 (naive level-truncation gives 0.86 vs true L=0.0056). **The cross-level (−1)^{|T|} alternation is ESSENTIAL** — the bound on |Σ_{|T|≥3}| cannot hold; L>0 survives only by the conditionally/Abel-convergent alternation of growing level masses. REFRAMED PROVE-target: control the alternating series Σ_k(−1)^k Λ_k with Λ_k growing — i.e. an Abel/Cesàro bound, not a termwise one. **The tool (THM-518 bridge): Bedert's level bound |E_k∩P| ≤ (C log|P|)^k** (Bonami hypercontractivity + Rudin Λ(q) + Bell–Chueluecha–Warnke sunflower) bounds exactly how many weight-k relations of an AP-core fall in any progression P, hence the relation counts r_k(ℓ) driving Λ_k; the cores live in [1,~14m] (log|P| small). Bedert's R̂(ℓ)=Σ_k r_k(ℓ)(−p/2)^k IS this singular series. Couple this with the OPEN-Q-104 stranger-decoupling (the m→∞ tail is finite; resonant strangers carry the inf). Entry points: THM-518, THM-503/504, THM-501, arXiv:2511.16636 (Bedert Lemma 4.3), the sharper extremizer {1,…,12}∪{14m}.

## OPEN-Q-098 🟢 Does the d-step Fibonacci sequence a_d(n) count a gap-d tournament family? (the Pascal-slope-d realization)

**Status:** OPEN (mac-mini-2026-06-14-S1, T819, HYP-2523..2525). The Pascal-slope-d family — row n = Σ_k C(n−1−(d−1)k,k), row-sums a_d(n)=a_d(n−1)+a_d(n−d), GF 1/(1−x−x^d) — has clean rigorous tournament meaning at d=1 (2^n = full tile-flip hypercube layer count; central binomial C(n,⌊n/2⌋) = metagraph width) and d=2 (Fibonacci = gap-2 independent sets on the n−1 base-path arcs; AND the realized circular-tournament iso-class count 2·Fib(m−2), S518). QUESTION: for d≥3 (Narayana's cows/supergolden, A003520/plastic, …), does a_d(n) count a NATURAL "gap-d" family of tournaments or staircase tilings — e.g. a d-omino tiling of the base path, or Hamiltonian/loneliness structures with a minimum-gap-d constraint — making the whole Pisot tower (φ, ψ, …, plastic) a genuine tournament invariant rather than an analogy? Sub-leads: (a) define the "d-graded metagraph" whose H-analog level sizes are the slope-d ridge sequence (d=1 mechanism = balanced level sets = random-walk excursion, gn-as-orbifold.md); (b) prove the 2·Fib(m−2) circular-tournament formula (S518, checked m≤9, no proof); (c) the plastic-number coincidence — d=5 here (x⁵=x⁴+1) shares its root with Padovan (x³=x+1, the monad/free-factorial THM-438 thread): is there a slope-5 ↔ monad tournament bridge? Entry: 04-computation/pascal_slope_d_family_macmini_0614s1.py, reflection the-pascal-slope-d-family-and-its-pisot-towers.md.

## OPEN-Q-099 🟡 Is there ANY positivity (flag/moment/SOS) certificate that cuts a baby-Hodge hole, or are they all pure integrality gaps?

**Status:** strong evidence for "pure integrality gaps" (mac-mini-2026-06-15-S1, THM-509, HYP-2527). Every realizable-region hole at n=6,7 is moment-interior (in the convex hull of realized (tr A³, tr A⁵) vectors AND skew-Hankel-PSD-feasible), and the flag-overlap-density Gram matrix is PSD at every hole (the holes are tournament-limit interior points, e.g. (c3,c5)=(8,10) = midpoint of realized (8,8),(8,12)). So NO finite continuous PSD relaxation tested cuts a hole — they appear to be pure integrality gaps, cut only by the Boolean/rank-1 constraint (the #P / permanent side of the Valiant wall). QUESTION: prove the all-n positivstellensatz — that NO polynomial moment inequality (any degree, spectral or overlap-side) excludes a hole — making "baby-Hodge hole = integrality gap interior to the flag-feasible body" a theorem; OR find a single hole that IS cut by some clever higher-degree flag certificate (which would be surprising and important). Sub-question: characterize exactly which integer lattice points interior to the flag-feasible density body are realized vs skipped — the discrete "non-algebraic Hodge class" locus. Entry: THM-509, baby_hodge_convex_certificate_macmini_0615s1.py, the c3-fiber stratification (the regular/near-regular score class is the richest hole source).

## OPEN-Q-100 🟢 The OCF Mayer cluster expansion: forbidden values as excluded volumes; does B₂≅T₄ seed a general-n structure?

**Status:** PART (a) DONE (kind-pasteur-2026-06-15-S7, THM-515). The explicit OCF Mayer cluster expansion is now worked out: ideal gas 3^{α₁}; Ursell excluded-volume integral −|E(Ω)|=α_2−C(α₁,2); 3rd integral P₂(Ω)−t(Ω). KEY CORRECTION: the clean −|E(Ω)| is the GRAPHICAL/Ursell integral, NOT the order-2 cumulant of log I (c_2=−|E|−α₁/2); only the analytic c_k satisfy exp(Σc_k z^k)=I, and z=2 is outside the radius of convergence (R≈0.0125 at Paley T₇) so the series is formal-only. p33=the 3-3 block of |E(Ω)|; TQ/Witt W_k are NOT gas integrals (they live on the spectral zeta of A). Forbidden 7=unique excluded cluster K₃ (single-cluster rigid, THM-029); 21=multiplicative gap (four realizable profiles). NOTE: the cluster picture does NOT give a new proof of "{7,21} only gaps" — the non-achievability half was ALREADY closed via Busch/A005517 (HYP-2271, monad-compute-2026-06-06; re-confirmed here), and the surjectivity half is provably outside cluster/multiplicative reach (a prime H needs a strong tournament of that H). REMAINING OPEN: (a′) does R(T)<1/2 for all T with |E(Ω)|>0 (would make the gas a convergent, not just formal, statement)? the explicit b_4/Q44 4th integral? (b) Does the B₂≅T₄ atom (THM-510: 4 iso classes of T₄ ≅ subsets of {a=+1, b=/2}, c3=|subset|, transpose=a↔b swap, parity↔transpose-type) seed a general-n operation structure, or is the count match (4=2²=A000568(4), Pascal row (1,2,1)) special to n=4? The gradings match is robust; test whether an additive×multiplicative operation monoid tracks higher iso-class structure. Entry: THM-002/029/502/505, THM-510, the Pascal-slope-d row-2 (T819).
## OPEN-Q-076 🟡 Do the triangular-tower Pell overlap families contain infinitely many prime block lengths?
## OPEN-Q-079 🟡 Do the triangular-tower Pell overlap families contain infinitely many prime block lengths?
## OPEN-Q-100 - Do the triangular-tower Pell overlap families contain infinitely many prime block lengths?

**Status:** OPEN (codex-2026-06-12; republished 2026-06-15 as HYP-2529/T821/OPEN-Q-100 after a namespace collision; extends HYP-2453/HYP-2450/HYP-2448). Every overlap block of length `L` gives the consecutive-support row `R_L(x)=1+x+...+x^(L-1)`, so prime lengths give exact cyclotomic irreducibility and prime values `R_L(2)=2^L-1` give base-2 Cohn certificates. The exact side hinge `[21,24]` has `L=4` and is reducible, so the right atom carrier is not “most exact overlap” but “prime length inside a Pell overlap family.” In the stored window `n<=10^6`, the whole-equation family `T_n=2T_m` already contributes prime lengths `5,29,5741,33461`, while the side families contribute only the short prime lengths `2,3,5`. The open problem is whether any of these Pell-family length sequences contain infinitely many primes, especially the whole-equation sequence `L=2m+1` along `T_n=2T_m`; a stronger subquestion asks for infinitely many Mersenne/Cohn exponents among those prime lengths. Files: HYP-2529, T821, `04-computation/triangular_tower_repunit_tournament_codex.py`, reflection `triangular-tower-overlaps-as-repunit-carriers.md`.

## OPEN-Q-102 🟡 Is the OCF a noise-stability functional? (FKN/Arrow on the tournament cube)

**Status:** OPEN (mac-mini-2026-06-15-S3, THM-512, HYP-2534..2536). The tournament/tiling cube IS the social-choice cube (Kalai's Arrow setting); a directed 3-cycle = a Condorcet paradox = the minimal odd cycle = the OCF obstruction; the transitive tournament = the rational/Arrow-dictator ground state; c3 = the Guilbaud level-2 quadratic (Var=3C(m,3)/16). QUESTIONS: (a) Kalai's P_rational=¾+¾·Stab_{-1/3}[f] and the OCF H=I(Ω,2) both spectrally encode odd-cycle/Condorcet content — is the OCF specialization (p1→1, p_odd→2, p_even→0) a noise-stability Stab_ρ at a specific ρ (mirroring ρ=-1/3 for the 3-cycle)? If yes, the forbidden H-values {7,21} become forbidden Condorcet-cyclicity spectra and robust Arrow (FKN) gives a stability statement about near-transitive tournaments clustering at the H=1 corner. (b) Write the multi-vertex Möbius sieve for H (HYP-2534), truncating at the band limit D=2⌊(n-1)/2⌋. (c) Which tournament invariants are juntas (Friedgut) / have a decisive arc (KKL)? Entry: THM-511/512, THM-002/163, Kalai 2002 (arXiv social-choice), Mossel 2012 (arXiv:1003.3956).

## OPEN-Q-104 🔴 Prove inf L(S)>0 via the Riesz-product method (the C'(14)/LRC(14) endgame)

**Status:** OPEN, the LRC(14) prize — the Riesz route is now DIAGNOSED (mac-mini-2026-06-15-S5, THM-518; was THM-515/HYP-2540). The reduction inf L>0 ⟹ C'(14) ⟹ LRC(14) is THM-398/501; L(S) = the lonely measure = ∫∏1_safe(v_iτ)dτ (THM-515). **HONEST VERDICT (THM-518): the Riesz product is the WRONG tool for the AP-extremizers, and neither of Bedert's two routes reaches 1/14.** (1) The interior-drop cores {1..13}\{j}∪{14m} are AP-like ⟹ small additive dimension dim₂~log N≈2–3 ⟹ Bedert's Riesz gain dim₂²/n³ is worthless. Exact direct-grid: ∏(1−cos) certifies 3/5 loose configs but FAILS both extremizers (j=6: ratio 1.064; j=12: 1.035), and per-speed amplitude optimization stalls at **1.0096 ≥ 1**. (CAUTION: a Fourier-truncated K=14 estimate gave a spurious 0.95 "certificate"; use exact direct-grid — see MISTAKES.) (2) The RIGHT tool for AP-cores is Bedert's **prime point-mass** measure (Lemma 5.3): ML ≥ ⌈(p−1)/26⌉/p; best admissible prime p=29 gives **2/29 = 0.06897 = 96.6% of 1/14** — short by 3.4%. The cores ARE loose (L≈0.0052), so the truth sits in the ~1–4% sliver between both state-of-the-art certificates and the optimum — that sliver IS the exact-value difficulty. NEW STRUCTURE (THM-518): (a) **stranger-decoupling** — lim_{m→∞} L({1..13}\{j}∪{14m}) = (6/7)·meas(Lonely({1..13}\{j})) (Weyl), collapsing each j-family's m→∞ tail to one finite positive measure (≥0.00699); (b) but the infimum is carried by **finite resonant strangers** (m=7, stranger 98=2·7², dips to L=0.00524 < the limit), which share the factor 7 with the core — these need separate control; (c) the **bridge**: Bedert's R̂(ℓ)=Σ_k r_k(ℓ)(−p/2)^k with non-dissociated relation counts r_k IS the project's singular series L(S), so the exact-value work is the project's relation-lattice computation, not the asymptotic machinery. NEXT: control the finite resonant-stranger set (the 7-power/square dilations); push the prime route past 2/29 with a composite/relation-tuned point-mass; or use Bedert's level bound |E_k∩AP|≤(C log|P|)^k (OPEN-Q-097). Entry: THM-518, THM-515/503/504, arXiv:2511.16636, Tao 1701.02048, 04-computation/lrc14_{riesz_verify,riesz_optimize2,stranger_limit,decouple_confirm}_macmini_0615s5.py.
**OPEN-Q-001** — Resolved by opus-2026-03-05-S1 via THM-008 (mu triviality bound). See THM-008.
**OPEN-Q-009** — Resolved by opus-2026-03-05-S1. mu(3-cycle at n=6) in {1,3}, determined by whether the 3 available vertices form a cyclic or transitive subtournament.
