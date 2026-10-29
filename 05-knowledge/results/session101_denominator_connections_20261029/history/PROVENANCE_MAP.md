# Session101 inheritance map: 101, 35, 175

Pinned tree: `e871dd4513c7dcaa75cb29bb9847aafac35a795c` (2026-10-29).
This is a **provenance and scope map**, not a claim that similarly numbered objects are equivalent. Numbered theorem IDs, table row numbers, decimal fragments, sample bounds, and duplicated inbox broadcasts do not count as mathematical occurrences. Different notes from one dated research campaign are identified as such, not inflated into separate sessions.

## Source verification and status

`SOURCE_BLOBS.json` records 40 selected historical paths. 36 locally materialized files have byte-identical Git blob IDs to the pinned tree. The one changed locally read file, `collatz_catalan_ramsey_20260927.md`, was instead read directly from the pinned tree; its corrected opening retracts a former tournament interpretation (the orientation was only a value total order), replaces it by an inversion graph, repairs a propagating-carry claim and a factor-two asymptotic error, and credits classical Ramsey results. The two-edit note and both later aliquot notes were also read directly by `git show` at the pin. All40 final read-blob entries match. This is stronger than assuming the older checkout is current. Source files in `05-knowledge/results` often contain elementary proved statements without a theorem ID. Their own PROVED labels below are reported with actual scope; historical reflections are not promoted to canon. No claim here implies Collatz convergence.

## Fast route to the strongest original sources

- [101 prime-square resonance and sharp squarefree cutoff](https://github.com/eliottcassidy2000/math/blob/e871dd4513c7dcaa75cb29bb9847aafac35a795c/05-knowledge/results/collatz_guards_20260921_squarefree.md)
- [175-to35 partial5-adic retention under two word edits](https://github.com/eliottcassidy2000/math/blob/e871dd4513c7dcaa75cb29bb9847aafac35a795c/05-knowledge/results/collatz_two_swap_repetition_classification_20261020/TWO_EDIT_REPETITION_CLASSIFICATION.md)
- [Weighted-mediant103/175 integrality hostile](https://github.com/eliottcassidy2000/math/blob/e871dd4513c7dcaa75cb29bb9847aafac35a795c/05-knowledge/results/collatz_mediant_tree_K_20260930.md)
- [Order7 Paley/Fano/rotational classification](https://github.com/eliottcassidy2000/math/blob/e871dd4513c7dcaa75cb29bb9847aafac35a795c/05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md)
- [175/35 Hamiltonian-path multiplicities](https://github.com/eliottcassidy2000/math/blob/e871dd4513c7dcaa75cb29bb9847aafac35a795c/05-knowledge/results/h_spectrum_multiplicities_20261021.md)
- [91,175,325 prime-incidence triangle](https://github.com/eliottcassidy2000/math/blob/e871dd4513c7dcaa75cb29bb9847aafac35a795c/05-knowledge/results/crossroads10_20260926_arithmetic.md)
- [27/35 aliquot quotient and doubling hostile](https://github.com/eliottcassidy2000/math/blob/e871dd4513c7dcaa75cb29bb9847aafac35a795c/05-knowledge/results/grid_aliquot_overlap_20261028/README.md)

## A. The primary 101 lineage: the gap prime is an active obstruction

### A1. September21 arithmetic guards — strongest older mechanism

Source: `05-knowledge/results/collatz_guards_20260921_squarefree.md`, §§1,3–6; companion `collatz_guards_20260921_synthesis.md`.

**PROVED (elementary p-adic lemma + finite linear-form sieve), with FINITE-EXACT modular certificate.** The exact valuation word `(1,1,5)` has return
`R(x)=(27x+19)/128`, gap `128-27=101`, and fixed point `19/101`; its three phases are `19/101,79/101,169/101`. Because the gap is divisible by101 and carry19 is a101-adic unit, R acts as a single cycle modulo every101-power. In the specified row `n≡5 mod9`, r repeats can have all `3r+1` nodes squarefree iff `r≤4985`; at4986 some displayed node is divisible by101². Modulo10201, survivor count is
`h(r)=(3795-r)_+ +(1419-r)_+ +(4986-r)_+`.
All other primes admit a simultaneous squarefree avoiding class;13 needs the repaired residue `f+13 mod169`, not the fixed point f itself. Finite-word sources still exist at every r: the cutoff is squarefreeness, not realizability or convergence.

**General mechanism:** if `F(x)=αx+β`, `α≡1 modp`, β a unit, then `v_p(F^t(x)-x)=v_p(t)`; hence a single `p^a`-cycle. For `(ax+b)/d`, sufficient conditions are `p|(d-a)` and `p∤db`. The SAME missing divisibility blocks an integral fixed point and forces finite prime-square coverage. Prime-power multiplicities remain essential if no such unit-carry prime exists.

**Transfer to current task:** retain `(gap,carry,phase roots)`, not only reduced denominator. After appending a step, test valuations at5 and7 and the new affine action; the old101 action does not automatically transfer. This is a true operation-level inheritance, not a common-number association.

### A2. Session100, October29 — denominator persistence and parameter conjugacy

Source in pinned working tree: `05-knowledge/results/session100_selective_collatz_targets_20261029/bijections/VALUATION_CONJUGACY_AND_INTEGER_SIDECAR.md`, §§3–4; audited in its `audit/INDEPENDENT_AUDIT.md`.

**PROVED:** identical valuation word `(1,1,5)` gives the positive5x+1 integer cycle `13→33→83→13`, but the3x+1 rational cycle `19/101→79/101→169/101`. Valuation conjugacy preserves the word and division clock, not integrality. Every finite backward prefix of a reduced rational terminal p/q with odd q>1 retains q in its reduced denominator: `gcd(2^a p-q,q)=1`; dividing by odd b cannot remove q. This is **prepending inverse history**, not recomputing the periodic point of a newly appended word. The present101→35/175 example is therefore not a counterexample to denominator persistence: it changes the periodic closure problem.

### A3. October19 reset-return campaign — 101 changes type

Five distinct notes in the same campaign: `collatz_odd_adjacent_return_20261019.md`, `collatz_maximal_run_reset_20261019.md`, `collatz_bank_reset_20261019.md`, `collatz_return_rank_obstruction_20261019.md`, and the inherited switching-height work.

The branch `H(x)=(81x+101)/512` on `x≡43 mod512` has fixed point `101/431`. Here101 is the **carry/numerator**, while431 is the slope gap/denominator. In the odd-adjacent-return note the bounded-head reset lemma and two all-L actual descent-receipt families are PROVED; universal entry is OPEN. Exactly `floor(v_2(y-101/431)/9)` complete H blocks can be consumed in the displayed family. Maximal-run and return-rank notes develop analogous phase budgets. This is useful as an explicit warning that “the101 obstruction” must be typed by its coordinate; number identity alone loses its role.

### A4. October14 successor-holonomy session

Source: `collatz_successor_holonomy_20261014b.md`, explicit terminal family near §6:
`n_s=(32·4^(10+27s)-101)/81`.

**PROVED** exact native cells, modular-return families, and explicit ROOT-certified subfamily; universal source coverage OPEN.101 is the carry subtracted during terminal decoding, not a common periodic denominator. The finite-group return does not intrinsically determine the sign of real descent. Connection: integer seams and source guards must accompany any finite residue recurrence.

### A5. September25–26 gates/sporadic-cycle repair

Sources: `procgen_gates_20260925_gate_equidistribution.md`, resolved correction §final; `procgen_sporadic_20260926_free_and_sporadic_cycles.md`, opening and census.

**FINITE-EXACT:** a missing primitive3x+d cycle at d14303 has least element101, total period2155 with1092 odd steps. It lay beyond the earlier period2000 scan; adding it repaired the Belaga–Mignotte table comparison.101 is a cycle **minimum**, not denominator (after scaling by d the denominator role belongs to14303). Lesson: a bounded clock census is not an all-cycle statement. The same clock-lattice machinery can conceal long sporadic cycles beyond apparently empty search windows.

### A6. September27 Catalan/Ramsey session — legitimate symmetry, not a tournament at101

Source: [corrected pinned Catalan/Ramsey note](https://github.com/eliottcassidy2000/math/blob/e871dd4513c7dcaa75cb29bb9847aafac35a795c/05-knowledge/results/collatz_catalan_ramsey_20260927.md), §§Paley table.

**FINITE-EXACT / CITED classical Ramsey bound:** the Paley **graph** P101 is a self-complementary graph with clique number5, yielding `R(6,6)>101`.101≡1 mod4, so quadratic residues define an undirected graph, not the usual Paley tournament. The note also proves Catalan/Dyck critical-ladder identities. The shared “101” does not supply a map from the graph to a Collatz denominator; a character/spectral analogy is weaker than an intertwining map.

## B. 175: genuine denominators, repeated-word gaps, and tournament path counts

### B1. September30 weighted-mediant / fixed-point session

Sources: `collatz_mediant_tree_K_20260930.md`, §4 Theorem F; `collatz_pentagon_fixed_points_argument_styles_20260930.md`, blocked strategy P1.

**PROVED:** a reduced rational fixed point m/d with d odd and3∤d corresponds exactly to a primitive numerator cycle of `m↦oddpart(3m+d)`; the denominator is preserved along that orbit. Concatenated words combine fixed points by a clock-weighted mediant, not a unimodular Stern–Brocot move. Explicit hostile:
`-17=[(-1539)(-1)+1400(103/175)]/(-1539+1400)`.
This splits the negative integer cycle word as `(1,1,1)(2,1,1,4)`; the second block has fixed point103/175. Integrality of a whole word does NOT descend to its blocks. In the current append experiment, the175 denominator is consequently a familiar rational-cycle carrier, not an unexplained accident.

The same note's **FINITE-EXACT** bounded census finds every d≤100 coprime to6 as a rational-cycle denominator, using words with total exponent≤22 plus direct searches for the four missing denominators53,65,67,79. The general every-d assertion remains conjectural.

### B2. September21 Pillai/convergent cycle-gate session

Source: `collatz_mod6_20260921_pillai_convergents_cycle_gates.md`, nonprimitive clocks.

The repeated constant word `(2,2,2,2)` has gap `4^4-3^4=175`, carry175 and reduced fixed point1. More generally repeated all-two words have the gap sequence1,7,37,175,781,…; cancellation is complete for the constant word.175 is thus a **raw closure gap**, not necessarily a reduced denominator. This is the exact hostile to identifying `2^D-3^J` with the reduced denominator without taking `gcd(C,2^D-3^J)`.

### B3. September22 Paley/Fano/octonion session

Source: `collatz_mod6_20260922_paley_fano_octonion_design.md`, S17–S19.

**FINITE-EXACT full order7 enumeration:** the2640 labeled regular tournaments split into three classes:

| class | H | automorphism order | cyclic triangles | skew polynomial |
|---|---:|---:|---:|---|
| Paley |189|21|14|`x(x²+7)^3`|
| rotational R7 |175|7|14|`x(x^6+21x^4+35x²+7)`|
| third regular |171|3|14|`x(x²+7)(x^4+14x²+1)`|

R7 contains no Steiner triple system among its cyclic triples; Paley has two disjoint Fano systems; the third has two sharing three blocks. Equal order, regular score sequence, and triangle count do not determine H, spectrum, or Fano incidence. Here175 and35 occur together **inside one coherent tournament calculation** (path count and characteristic coefficient), but no relation to the Collatz reduced denominators has been established.

### B4. June1 oracle S518/S522/S523 — early polygon lineage, historical only

Sources: reflections `lrc-as-a-tournament-and-its-fibonacci-isoclass-menu-s518.md`, `tournaments-as-polygons-lrc-twin-bridge-s522.md`, `lrc-is-a-round-tournament-necklace-s523.md`.

The useful verified seed: orient an odd regular polygon by the forward semicircle; the heptagon gives R7 and H175, whereas even regular polygons have antipodal ties and are NOT tournaments. P7 has the same regular scores but is not the same circular realization. The reflection's extrapolated Fibonacci count and broad LRC equivalence are **not inherited as proved here**; these are historical synthesis/seeding surfaces. The fixed R7 count is independently supported by B3.

### B5. July1 kind-pasteur S9 and opus heptagon sessions — symmetry dictionary

Sources: `the-atoms-are-a-heptagon-the-census-is-a-polygon-tournament-and-the-extremizer-is-the-pentagon-kps.md`; `the-heptagon-tournament-6-units-plus-antipode-are-roots-of-z7-minus1-a-self-complementary-D7-tournament-Cayley-dual-to-the-runners-opus-20260701.md`.

**Historical verified constructions:** odd residues modulo14 form a heptagon; rotations preserve R7, reflection reverses every arrow, so oriented Aut(R7)=C7, not D7. The full geometric D7 includes anti-automorphisms. Multiplication by units modulo14 is a different action and does not automatically preserve heptagon edges. This is a direct useful model for distinguishing valuation-word cyclic rotation (an orbit symmetry) from word reversal (generally not a Collatz dynamical symmetry).

### B6. October22 Sierpinski/Steiner/Collatz session

Source: `sierpinski_steiner_collatz_tournaments_20261022.md`, §7.

Explicitly types `H(R7)=175=7·25`, `H(P7)=189=7·27`: because these tournaments are vertex-transitive,25 and27 are their per-start-vertex path counts. This is an exact orbit-count mechanism. The accompanying coincidences with27's Collatz orbit and70 halvings are explicitly NUMEROLOGY. A numeric factor5 between35 and175 should similarly not be mistaken for a symmetry quotient without a specified five-to-one map.

### B7. September22 order-laws / October5 receipts / October4 rank sessions

Sources: `collatz_procgen_20260922_order_laws.md`, `collatz_receipts_four_questions_20261005.md`, `collatz_branch_toll_rank_20261004.md`.

175 and35 are actual orbit hubs. The order-law scan's plus crossings all touch the27 highway; the receipts scan reports175 as2938 product-join hits. These are **FINITE-EXACT frequency effects**, not algebraic multiplicativity. The rank paper supplies the **PROVED guarded family** `223+20736t ←175+16384t` by cutting a child to its actual successor; this is source-preserving family transport, not universal coverage. It is another distinct role for175: an integer source/intermediate vertex, not denominator.

## C. 35: composite closure, higher faces, and failed coarse transports

### C1. October6 cycles/tubes audit: primitive3x+35 cycles

Source: `collatz_cycles_tubes_debt_walk_openai_20261006.md`, cycle-count correction near §related literature.

**Reported corrected finite cycle evidence:**3x+35 has at least nine positive cycles with minima7,13,17,25,35,133,161,1309,2429; two are primitive, minima13 and17. Scaling those primitive cycles by1/35 places them on the rational3x+1 denominator35 carrier (Theorem F above). “At least” is essential; no complete all-cycle classification is asserted. This supplies actual prior denominator35 content, not just the session number35.

### C2. October7 completion campaign: parent35 is the canonical hostile

Sources: `collatz_completion_paper24_23_20261007b.md`, §fixed-parameter correction; `collatz_completion_synthesis_20261007b.md`.

**PROVED explicit obstruction:** inverse-one operations give `35→23→15`, but the child15 has empty required request interval `(Dmin,Dmax)=(2,1)` in the named scheme. The authentic receipt `35 --(1,5)→5←(2,3)--17` cannot simply be transported to15. A supplied ROOT path `35 --(1,5,4)→1` does ground the relevant obligation. Missing data are source, word seam, and actual terminal receipt—not clock counts. This is the same type of quotient failure as equal `(J,D,q)` failing to determine the next denominator, but not the identical theorem.

### C3. August21 higher-face/Kneser session: the two35s really are dual

Source: reflection `tournament-subset-decks-higher-order-faces-and-kneser-middle-duality-codex-20260821.md`, §§1,4.

**PROVED additive-incidence identity:** for n7, both triangle and four-set layers have35 addressed cards. Complementing the four-sets gives the Kneser matrix `KG(7,3)`; its inverse has entries
`(-1)^|S∩R|/[4 binom(3,|S∩R|)]`.
Thus the two35-dimensional layers are losslessly dual for **addressed additive face fields over characteristic zero**. This does NOT say that unlabelled tournament-type histograms are equivalent or that independent higher orientations come from edge tournaments. Missing sidecar: overlap addresses. Strong higher-order connection to the current task: a quotient that forgets the carry or phase may be repairable by a small addressed coordinate, but equality of dimensions is not itself a reconstruction proof.

### C4. Spectral-gap historical session, repaired by October21 census

Sources: reflection `spectral-tournamentgraffiti-a-forbidden-eigenvalue-count-and-h-dominates-the-skew-determinant-deathstar-S78.md`; current result `tournament_pseudoprimes_camion_fermat_20261021.md`, §6.

**FINITE-EXACT:** H35 is missing at order6 but occurs at order7; holes at order6 are7,21,35,39, while at order7 they include7,21,63,107,119,149 (plus ceiling run). “35 forbidden” without specifying order is false. The same result derives that a composite H whose every nontrivial factorization uses forbidden7 or21, including35, must have only one non-singleton strong component. This uses the exact product law for H under ordered strong-component joins;35 is an arithmetic factorization obstruction for a tournament decomposition, not a cycle denominator mechanism.

### C5. October1 two-sheet clock tournament session

Source: `procgen_tcpc_20261001_tournament_clock_prime_collatz.md`, §3.

**FINITE-EXACT existence certificates:** the two-sheet clock D17 on34 vertices is all-odd (every arc lies on an odd number of Hamiltonian paths), although35=5·7 is not a prime power so the usual Paley construction is unavailable. D7 and D19 similarly cover orders14 and38. General all-order existence is not settled by these three cases. This is a positive example where a refined two-sheet carrier bypasses a prime-power obstruction without pretending that35 itself is a field order.

### C6. October21 scale-credit/tournament-compression session

Source: `scale_credit_and_tournament_compression_20261021a.md`, inheritance table and guarded inverse rays.

The native inverse rays `s→64s+35` and `s→64s+49` shift clock by6 and rotate phase. Three phases restore phase, but fixed-parent rays provide only boundedly many annulus addresses, not linear supply.35 is an affine intercept; the source/address obstruction survives phase closure. This is directly useful against treating a common appended clock or periodic phase as evidence of universal source coverage.

## C7. October20 two-edit repetition classification — exact175/35 ancestry

Source: `05-knowledge/results/collatz_two_swap_repetition_classification_20261020/TWO_EDIT_REPETITION_CLASSIFICATION.md`, Theorem R2 and §6.

**PROVED elementary repetition-factor classification; PROVED mechanical consequence using independently audited one-swap theorems; FINITE-EXACT hostile controls.** For arbitrary binary base v and repetition d≥2, write `D=(2^q0)^d-(3^p0)^d=H D0`. Exactly two cyclic adjacent unequal-bit swaps preserve full repetition factor H in the carry iff the second undoes the first, or d2 and the edited edges are one base period apart. Primitive two-edit frontier remains OPEN.

The printed hostile is `v=0011,d=2,H=25,D=175`, with
`00110011→00101011→00011011`, final carry680 and `gcd(H,680)=gcd(175,680)=5`. Hence its fixed-point fraction is680/175=136/35. This is a prior exact appearance of **raw175, surviving5, reduced35** under a different word operation. “Every prime of H is destroyed” would be false; only prime-power depth is lost. This is the closest older prime-power-cancellation analogue of the current same-clock append separation, without equating swaps and appends.

## C8. October21 opus S50 — actual35 and175 class multiplicities

Source: `05-knowledge/results/h_spectrum_multiplicities_20261021.md`, §3, independently reproduced census.

**FINITE-EXACT:** at order7 H35 has2 isomorphism classes; H175 has exactly1, the rotational regular tournament. All order7 tournaments whose H is divisible by7 are strong, while none at order≤6 has7|H. Other order7 multiples have multiplicities49×7,77×10,91×15,105×4,133×6,147×1,189×1. **General symmetry identity:** class count at fixed H is congruent modulo2 to its self-converse class count, because converse pairs the rest. Unique H175 therefore forces self-converseness at this order. The empirical explanation of certain spectral holes by product borrowing was retracted; that story is not inherited.

## C9. September26 crossroads ten-vertex arithmetic —175 in an odd prime-incidence cycle

Source: `05-knowledge/results/crossroads10_20260926_arithmetic.md`, actual-orbit construction near lines153–171; inherited `THM-4493-collatz-gcd-height-rank-certificate.md`.

**PROVED, independently audited by the root coordinator and geometry lane**, with FINITE-EXACT bounded controls. The actual odd27 orbit contains `91=7·13`, `175=5²·7`, `325=5²·13`. Each is the product of its pairwise gcds with the other two, so neither private primes nor strict gcd dominance certify them. The odd exclusive-prime-incidence triangle nevertheless forces all coefficients of a multiplicative relation to vanish. Seven further actual orbit nodes have private primes, giving the explicit ten-node independent family. The general proof uses endpoint prime exponents; it does not assert that arbitrary graph realizations are chronological Collatz paths. A symmetric gcd/incidence graph should not be cosmetically tournamentized. The new Session101 comparison of rational independence with integral Smith/squareclass saturation can cite this proved older triangle, rather than rebrand its independence as new.

## C10. October14 number/tournament trichotomies —35 in the aliquot geometry

Source: `05-knowledge/results/number_tournament_bridges_20261014.md` with its exact output `number_tournament_bridges_20261014.out` (output was used as a discovery pointer, not a new rerun).

For35 the nontrivial divisors are5,7,35; `σ(35)/35−2=−22/35`. The note identifies this deficient/abundant/perfect sign with the angle excess of the triangle with angles π/5,π/7,π/35 (the exact dictionary requires τ(n)=4). This is a role-preserving geometric encoding of the aliquot ratio, not a statement about Collatz35 denominator. A targeted search of aliquot-named materialized files found no substantive175 occurrence; no175 aliquot claim is manufactured merely because it is easy to compute an aliquot sum.

## C11. Session66 October26 and October28 aliquot fibre/dilation sessions

Sources: `05-knowledge/results/session66_root27_aliquot_graph_interfaces_20261026/session66-aliquot/REPORT.md`; `05-knowledge/results/grid_aliquot_overlap_20261028/README.md`, §retained carry.

**FINITE-EXACT complete inverse fibre:** Session66 certifies `s^{-1}(13)={27,35}`, not just two witnesses; its proof bounds make the finite enumeration exhaustive. **PROVED elementary transport hostile (Oct28):** doubling these sources gives54 and70, whose aliquot sums are66 and74. Equal current target13 does not determine next target under the same doubling operation. The missing source/carry coordinate repairs it. This is an unusually close cross-problem analogue of equal periodic denominator101 not determining the new denominator under the same append operation: specify quotient, lawful operation, separating witnesses, and retained coordinate. No theorem transfers between aliquot and Collatz from this analogy alone.

## Coverage count (conservative, by substantively distinct prior sessions/campaigns)

| number | documented sessions/campaigns | counting rule |
|---|---:|---|
|101|6|A1 Sep21; A2 Session100; A3 Oct19 (all reset notes together); A4 Oct14; A5 Sep25–26 corrected gates/sporadic campaign together; A6 Sep27|
|35|12|C1 Oct6; C2 Oct7; C3 Aug21; C4 spectral-gap lineage plus Oct21 current census counted as one; C5 Oct1; C6 Oct21 scale-credit; C7 Oct20; C8 Oct21 opusS50 distinct from pseudoprime campaign; C10 Oct14; B3 Sep22 spectral coefficient; C11 Session66 Oct26; C11 Oct28 quotient/dilation|
|175|12|B1 Sep30; B2 Sep21; B3 Sep22 Fano; B4 June1 oracle grouped despite3 sessionIDs; B5 July1 kind-pasteur/opus grouped; B6 Oct22; B7 has3 distinct campaigns (Sep22 order-laws, Oct5 receipts, Oct4 rank); C7 Oct20; C8 Oct21S50; C9 Sep26 crossroads|

These counts deliberately undercount multiple notes in one campaign and exclude the October29 triggering append example. They are coverage counts for this map, not claims that every use is independent or equally relevant. Cross-field101 Paley graph is typed as nonmechanistic; strongest101 mechanism is A1. Strongest35/175 arithmetic connection is C7, while B3/C8 are actual tournament counts. There are40 byte-verified source files, not40 independent sessions.

## D. Higher-order synthesis: what genuinely transfers

1. **Marked operation state is the common mechanism.** Collatz needs `(J,D,C,source/phase)` to transport a word; tournament higher decks need overlap addresses; round-tournament reflection needs orientation/anti-orientation type. They are examples of the same *information-loss diagnosis*, not a claimed common mathematical object.
2. **The101 squarefree theorem suggests a new exact experiment.** For each original and appended word, compute all gap-prime valuations and intermediate affine zero roots. Compare orbit decompositions modulo p² for p5,7,101. A reduced denominator alone cannot say whether a prime-square obstruction is one full affine orbit, multiple cycles, or cancellation.
3. **Closure vs transport must be separated.** Forward movement around one rational cycle preserves reduced denominator. Prepending a finite inverse word cannot remove its terminal denominator. Appending a letter and solving a NEW periodic fixed-point equation can change denominator. Word splitting can destroy integrality (the -17/103/175 hostile). These four operations must never be conflated.
4. **Useful symmetry comparison:** cyclic rotation of a periodic word changes its marked numerator/phase while retaining its rational orbit and reduced denominator. R7 rotations retain H; reflection is anti-isomorphism. Appending at a chosen seam is a marked operation and may break rotation equivariance unless the seam is transported too. This is the most promising precise tournament/Collatz analogy to test.
5. **35 and175 are not interchangeable shadows.** They can be raw gap, reduced denominator, integer vertex, path count, polynomial coefficient, or face-space dimension. Each role above is substantive; no unconstructed map identifies them. The factor175/35=5 in the target experiment should be explained by prime-adic carry cancellation first.

## E. Exclusions and research boundary

Ignored as noise: theorem/session IDs101/35/175, page numbers, rows, scan endpoints, decimals0.101/1.175, binary strings101, counts such as175 tests, and mirrored broadcast copies. Separate notes within the October19 and October7 campaigns are not called independent discoveries. Historical proof claims without a current correction check are marked as historical rather than canon. This map recovers over a dozen genuinely distinct dated/session lanes, with40 source files, but does not claim an exhaustive repository census or an all-cycle theorem.
