# LRC(14) FRONTIER — 2026-07-15 (corrected synthesis)

**Correction (2026-07-14).** The original S297 version of this file overlooked the scale-quotient
obstruction already present in its parent commit (HYP-6780 and the corrected THM-758). Its claims
that the `f>=4` branch was globally finite-decidable, that all such families have `M>=0.097`, and
that LRC(14) was assembled were false. This file now records the logically current frontier.
The standard 14-runner case remains open.

**External-status note (checked 2026-07-14).**  The April 2026 arXiv preprint
[Sungkawichai--Trakulthongchai, *Eleven, twelve, and thirteen lonely runners*](https://arxiv.org/abs/2604.23906)
claims a computer-assisted proof for `k=10,11,12` nonzero speeds.  In the
usual convention with one stationary runner, that makes the case studied here
(`k=13`, fourteen runners, threshold `1/14`) the first unresolved case if the
preprint is accepted; it is currently an arXiv v1, not a peer-reviewed theorem.
The finite-height input used below is independently aligned with the published
zonotopal finite-checking theorem of
[Malikiosis--Santos--Schymura](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/linearly-exponential-checking-is-enough-for-the-lonely-runner-conjecture-and-some-of-its-variants/A51A991DE89B8C9C2E2FF13FBD4501DA),
which bounds the relevant velocity sum by a linearly exponential quantity but
does not itself settle this case.

## 0. The honest reduction

LRC(14) = **NON-COVERING** (THM-366, via LRC(≤13), settled) + **COVERING**, and by
klein-S309's THM-758 far-count split:

> **COVERING = [≤3 elements > 14] ∪ [≥4 elements > 14]**
> - **≤3-far ⟹ ≥10 speeds in {1..14} ⟹ kps THM-738** (the 1001-body exact-ℚ tree, PROVED).
> - **≥4-far remains OPEN in general.** THM-755 gives a finite interval only after a core is
>   fixed. Under dilation `P -> cP`, the good-set measure is unchanged, its positive-length component count is
>   multiplied by `c`, and the cutoff obeys `v*(cP)=c v*(P)`. Therefore no uniform raw speed
>   cutoff follows.

HYP-6780 gives an explicit unbounded primitive covering ray

`V_c={c,2c,...,12c,13c+1}`

for infinitely many `c`, with `f=13`, lying below the scaled THM-755 edge and satisfying the exact
value `M(V_c)=1/13`. This refutes `f>=4 => M>=0.097` and the sampled cutoff near `500`; it does not
refute LRC. The correct target is a classification modulo scale, retaining normalized core shape
and the last runner's residue/offset.

THM-760 subsequently closes this displayed ray and, more generally, every primitive family in
which twelve speeds share a common divisor: lifting a core witness through the divisor produces a
sheet on which the one coprime exceptional runner is safe. The genuine scale residual therefore
has at least two exceptional residue classes, or has no twelve-speed common-factor core.

**THM-761 (opus-S299) now covers the multi-exception half of that residual:** for V = cP ⊔ W
with r = |W| ≤ 6 exceptions (any gcd strata), the union of bad sheets cannot cover the sheet
cycle once `Σ_a g_a(⌊c/(7g_a)⌋+1) ≤ c−1`; uniformly, coprime exceptions dodge at every c ≥ 43,
the exact failure sets for c ≤ 42 are computed, and r = 1 closes at EVERY scale and EVERY
stratum (extending THM-760 to non-coprime exceptions). Related: **codex-S3 REFUTED the uniform
q ≤ 25 good-period finish** (family 26·{1..12} ∪ {339}: first witness q = 27 — witness
denominators scale with c; bounded-q banks are scale-blind, sheets are scale-native; the same
family closes by THM-760/761 in one line). The scale residual after THM-761: (i) r ≥ 7 decks
(union wall structural; the zero-defect subcase is a cyclic tiling of Z_c, but maintained covers
may have overlap; the wall is realized at c = 7 by a family
that is still lonely via a non-sheet route), (ii) c ≤ 42 small-scale criterion failures (a named
dilation has bounded inflation, but `c*` does not bound the global normalized domain), (iii)
gcd-descent bookkeeping, and (iv) scale-free transverse fragmentation. HYP-6830's proposed
`r_P≤B(c*)` bridge is refuted generally and inside exactly `f=4` by
`{1,...,9,15,110,N,1092N}`; the open splice must be peel-relative
(MISTAKE-145 and the exact four-far falsifier below).

Status cautions: THM-724's addendum closes its genuine single-killer case, but THM-726 still relies
on an unproved global far-element monotonicity statement; THM-741 is explicitly `CLAIMED` with an
unfinished evidence checklist. Neither may be used as an unconditional assembly lemma.

## 1. Status by piece

| piece | status | artifact |
|---|---|---|
| non-covering | PROVED | THM-366 + LRC(≤13) citation |
| ≤3-far | FINITE-EXACT AS RECORDED; independent rerun/Lean transcription open | kps THM-738 + THM-758 split (klein-S309) |
| ≥4-far branch | OPEN modulo scale-normal classification | corrected THM-758 + HYP-6780 |
| sampled raw bands | FINITE-EXACT FOR THE STATED BANKS ONLY | mac-mini-S105; klein-S312 |
| proposed uniform `q<=25` good period | FALSE; exact coherent and gcd-incoherent residuals | THM-762/764; MISTAKE-143 |
| (H)-bands, bottom cores | CLOSED (complete sweep) | THM-756 (4,032 pairs; AP/GW corners) |
| safe-peel | Parts A/B PROVED; irreducible-tiled Part C empirical | mac-mini THM-753 |
| aligned monotonicity | PROVED | mac-mini THM-751 |
| shadow tiles | PROVED | THM-748 (klein), THM-749, THM-754 clean-slot |
| named coherent/cluster families | PROVED at their stated scopes | THM-668/737/739/740 |
| 12-speed common-factor core + one coprime exception | PROVED, all scales | THM-760 |
| (13−r)-speed common-factor core + r ≤ 6 exceptions | PROVED when the exact Σg_a budget fires; c ≥ 43 is uniform only for coprime exceptions | THM-761 (opus-S299) + battery |
| r = 7 unramified deck stratum (`7\mid c/g_a` for every owner) | PROVED above reduced-shape bound `max_a(w_a/g_a)>=7 max(P)`; the S299 wall is closed at a switching time | THM-771; corrected core of THM-767; MISTAKE-146 |
| prime lens `c=7`, any unramified owner count | Exact token polynomial: coverage iff `X^7-X` divides `product(X-k_a)`; seven-owner exact states map to all 25 masks at heptagon node `n7-a267`; any covered `r=8` wall is a simple event with a seven-owner heptagon stalk | THM-773 + exact 5,040-state/3,003-profile audits |
| prime-lens endpoint transport | Pairwise midpoint clocks are centered mechanical words with an Euclidean parity cocycle; centered Beatty ranks reconstruct every simultaneous wall and drive the exact `F_7` skew product.  The named r=8 row has 10 simple covered walls with palindromic owner word `162,108,108,206,197,197,206,108,108,162` | THM-778 + 6,400-pair/five-movie exact audit |
| merged-node tiling fibres | A fixed-path tiling maps to its canonical class and converse-merged node; the inverse node fibre is intrinsically `union HP(T)/Aut(T)`.  All n=3..7 atlases round-trip.  For `n7-a267`, `H=175` and `|Aut|=7` explain the exact 25 masks. | THM-781 + 33,866-tiling exact audit |
| prime-lens r=8 blocking chain | Full blocking is exactly piece surjectivity + wall rainbow + no simultaneous walls; the redundancy cocycle gives an exact event-word test and rooted states form an `A8` torsor. Raw wall length is unbounded by persistent-stalk refinement. Active-period count is also unbounded at `ceil(f/g)=2`: THM-794 repeats a full seven-visitor packet whose return is diagonal translation, hence zero in `F_7^8/Delta`. The live target is the stalk- and holonomy-quotiented packet path with metric lift. | THM-779/784/788/794; MISTAKE-147/149; exact refuter/A8/holonomy audits |
| r=8 visitor/extent laws | Phi recurrence, signed balance, and factor-two span laws are proved. A per-tuple finite bound follows if the balanced-cluster hypergraph has fractional speed-weight `<g`, or if the complete-`g`-period target lies outside its exact packet polytope. The dense `(65,64)` row closes with two periods; exact `(69,29)` data survive both marginals. THM-794 refutes the universal extent and active-period compactness claims while satisfying both positive tests: `1-W_*/g=42/g` and `dist(r,P)=49/(2g)` both vanish. The missing carrier is the common Beatty-ordered collision path, reduced deck holonomy, and metric/core incidence. | THM-783/786/788/794 + exact transversal/packet/holonomy certificates; MISTAKE-147/148/149 |
| eight-owner `c=7` buffer rigidity | The chamber/rainbow condition is subsumed by THM-779's integer token walk. The global `1/7+O(gcd/w)` partner-buffer law remains VERIFIED but must be restricted to the core-safe exit set; whole-core blocking and the general-density stalk/holonomy-quotiented switch bound remain OPEN. | HYP-6840 + THM-779/783/784/786/788/794 |
| seven-owner deck defect / ramified residue | Exact identity `F=Q+Omega-sigma`; exact tilings are chamber-locked, KCL necessity is WITHDRAWN, and mirror coincidence is diagnostic. Primitive `c=21` row realizes `(0,12,12,0)` | THM-771 + corrected THM-767 + exact audits |
| raw positive-length fragmentation bound `r_+(P)<=B(c*)` | REFUTED by an exact four-far family and an independent census; peel-relative `rho` is measured at most `9.335` only on the stated bank | THM-798; HYP-6830 correction; MISTAKE-145 |
| safe-measure floor / normalized band bridge | `rho<=12/(pi|G'_P|)` and `|G'_P|>=1/(91 maxP)` PROVED; phase pigeonhole gives the height-free floor `182^(-12)`; exact `maxP<=18` floor is `7/858`, unique at `{1,...,13}\{6}`, while that sharp global value remains CONJECTURAL | THM-777/780 |
| positive good-set state `(mu,r_top)` + frequency `N` + proportional peel `aN` | PROVED conservative transition `mu_N>=6mu/7-2r_top/(7N)`, `r_top,N<=N+r_top`, hence eventually capped; fully lacunary flags are terminal in every far-count stratum, with factors `412,405,394,27,17,14,13,13,13,13`; an exact 2,002-core floor gives a complementary factor-19 cone on `f=4` | THM-799; THM-798 is a full-family instance |
| primitive tight 12-speed locus | UNIFORMLY FINITE (`sum A<=78^11`), not classified | THM-763 |
| hereditary primitivity of tight 12-sets | PROVED; every leave-one-out core is primitive | THM-765 |
| unique-largest-13-multiple tight branch | IMPOSSIBLE by explicit prime-grid perturbation | THM-768 |
| arbitrary binding scale `p/(13s)` | PROVED packet split; shallow iff full nonzero residues; deep exceptions obey exact sheet capacity, with `s=2` parity and `s=3` colour criteria | THM-769 |
| two-/three-sheet equality quotients | PROVED primitive divisor transfer and speed bounds; two-sheet core contains multiples of every `2,...,12`, three-sheet core of every `2,...,11` | THM-772 |
| two-sheet metric residual | PROVED exact folded diamond `||(x+y)tau/2||+||(x-y)tau/2||>=11/13`; sharp measure cap `8/117`; all quotient cores in `[1,19]` closed against unbounded odd exceptions | THM-774 + exact certificate |
| odd-exception divisor-grid / signed-wall selector | PROVED for every odd divisor `q`: deep unit classes outside the explicit opposite-exception shell force global erosion escape. At mandatory `q=13`, one-sided wall leakage eliminates double-13 and full support; the sole survivor of this prime-grid signed-wall gate is the exact signed complement `U mod 13=(Z/13Z)^*\{+/-y}`, with `M(U)>=2/13`, `x<=2B-1`, and `y<=B-1`. A gate-sharp row silences every exception-divisor grid but escapes at denominator 17; THM-803's anti-grid catches it. | THM-797 + exact 10,971,770-grid / 352,716-signed-profile replay |
| q=13 anti-grid / all-component selector | PROVED: every odd exception divisor has a mandatory half-grid; at `q=13` this forces full parity-twisted support, and the complete universal even anti-grid ladder is exactly `d=2,4,6`. The full erosion predicate is equivalent to finitely many owner-labelled component endpoints and folded cusps. THM-817 sharpens the bound from `200B^2+22B` to the adaptive `2c_E N_R+2W-2g<=20B^2+22B-2g`, where `N_R` is the exact number of return cells. A sharp signed-complement row passes every divisor/anti-grid and global-maximizer test but escapes at the nonmaximal singleton `7/22`. Uniform failure of the selector remains open. | THM-803/817 + exact half-grid/support/cell replays |
| return-cell selector compression / atomic-stalk boundary | PROVED / VERIFIED: the mandatory central return component gives a selector of size at most `42B-2` when the closed return set is connected. In general every return component is one signed max-speed cell inside a tooth, so `N_R<=B` with an exact gcd label sieve and the adaptive selector above. An explicit primitive, divisor-complete, exact signed-complement family has `N_R=3+1440n=Theta(B)` while passing the present arithmetic and scalar gates; connected, bounded-component, and sublinear return geometry therefore cannot be inferred from those gates. Uniformly over positive odd exception pairs, THM-821 proves that each return-cell/deep-component atomic verdict factors through the pair-indexed circular sum arc and its endpoint/cusp margin. The fixed `(13,5)` 9,974-atom audit has 492 successes and 9,482 failures and explicit mixed fibres for return-only, deep-only, owner-only, width, and event-shape shadows. The linear-satellite family has a uniformly failing central stalk at that pair, but all-core and all-pair noncontainment remain open. | THM-807/817/821 + exact connected-return/grid/cell/stalk replays |
| fixed-ratio full-return radius budget | PROVED / VERIFIED ALL-SIZE: for `(x,y)=(13,5)`, the folded target is exactly the two radius-`2/169` cells centred at `5/13,8/13`. If `E` is nonempty compact and `R=-R` is compact with `0 in R`, then `E+R` lies in the target iff `rho_C(E)+rho_0(R)<=2/169`; the common-dilate intrinsic criterion is `rho_(C,d)(E)+rho_d(R)<=2/169`. Its phase rewrite `max_E||13dt||+13max_R||dr||<=2/13` is equivalent only together with `dE subset H`; alone it is merely a necessary thickness tax, with `E=R={0}` the guardrail counterexample. The no-switch lemma is what makes the intrinsic scalar lossless. It does not apply to one nonsymmetric satellite, force budget failure for every core, or classify arbitrary odd ratios; owner intervals remain the transport carrier. All 214 replayed cores fail the intrinsic global budget, but this is not a uniform arithmetic proof. | THM-824 + corrected exact 12,159-packet / 214-core replay |
| folded-target component and switch frontier | PROVED SHARP / VERIFIED ALL-SIZE: for `A>B`, reduce by `g=gcd(A,B)` to coprime opposite-parity half-frequencies `(alpha,beta)`. THM-831 gives every primitive target component by an exact Bezout-offset ball formula. The nonempty primitive symmetric-return predicate has a lossless no-switch radius factorization exactly for `4<=alpha<=9` (sixteen antipodal two-ball rows); `alpha<=3` is empty, while every `alpha>=10` has a strictly negative ordered switch triple, already at offsets `+/-1,+/-3`. The abstract factorization also requires the sharp winding guard `4h<1`, automatic for every folded component but false for arbitrary large balls. A common gcd creates raw deck switches while multiplication by `g` preserves the quotient-scaled criterion. The exact replay checks 518 primitive pairs, 6,160 compact packets, all endpoint and winding cases, and the four `alpha=10` witnesses. This classifies the legal scalar quotient; it does not force an LRC core to violate it or factor one nonsymmetric satellite. | THM-831 + exact formula/switch/packet replay |
| fixed-ratio thin owner shells | PROVED SHARP FOR THIS LOCAL STEP + FINITE-EXACT: for ten distinct speeds with THM-836's exact signed-complement residues and guarded `(13d,5d)` containment, the two directional exit coefficients are both permutations of `(9,20,31,42,53,64,75,86,97,108)`. In the band `B<13d`, the two coefficient-nine owners force `delta(d)<=11s(13d+s)/(2(117d+11s))`, `s=2B-13d`. Modular alignment empties both `s=1` and `s=3`, hence `13d<=2B-5`. Exact local `s=5` owner survivors show this two-sided packing argument stops there; no global core survivor is asserted. | THM-836 + 1,502,751-shell Fraction replay |
| shallow full-residue locus through lift height 12 | FINITE-EXACT over `13^12` conceptual packets; 13 dilates, unique primitive row `{1..12}` | THM-770 + exact owner-CSP |
| shallow AP Hamming-one star at arbitrary height | PROVED scale-free: every proper residue-preserving one-coordinate lift of `{1,...,12}`, and of every unit AP dilation, is loose. Every residual non-AP shallow packet differs from every AP dilation in at least two labelled coordinates. | THM-795 + exact core-threshold/atom/deck replay |
| shallow AP Hamming-two star at arbitrary height | PROVED scale-free: oriented half-open splice decks force both replacements back to the common AP scale at exact tightness; the normalized proper double lift then has the sharp floor `M>=2/25>1/13`. Thus every residual shallow packet has labelled Hamming radius at least three from every AP dilation. | THM-800 + exact 600,756-row replay |
| shallow AP Hamming-three star at arbitrary height | PROVED scale-free: THM-804 forces every hypothetically tight three-replacement lift to common AP scale. THM-806 then proves every proper normalized triple lift loose, so the entire residue-preserving Hamming-three star is closed at arbitrary height and scale. | THM-804/806 + exact deck/collar replays |
| scale-one Hamming-three collar closure | PROVED / FINITE-EXACT: a universal owner collar and half-open handoff digraph force one replacement into `[14,24]`; settled 10-/11-speed bounds plus exact two-comb geometry give `v<=262`, `w<=12v`; the legacy superset replay through `v<=381` rejects all 5,713,539 canonical rows. The closest certificate surplus is `2366/835198`, and the height-one minimum is `2/21`. | THM-806 + exact component replay |
| shallow Hamming-four closure | PROVED / FINITE-EXACT: THM-810's oriented deck dichotomy is exhaustive. THM-815 closes the common-scale branch twice: a sharp safe-component discrepancy ladder (`x<=105,v<=118,w<=83,z<=50`) with two 35,640-row terminal sweeps, and an independent collar/doubling reduction with 768,735 C++-certified rows. THM-816 closes every arbitrary lift on the order-three coset interface by a 7,909-state residual-component recursion plus an independent 132,510-row endpoint replay. The lift-invariant `q=39` clock lies on the boundary of the strict core-safe set; the proofs find surviving open components elsewhere. Hence the full AP-centred, residue-preserving Hamming-four star is loose at every scale and height. | THM-810/815/816 + exact deck/component/endpoint/collar replays |
| scale-one Hamming-five closure | PROVED / FINITE-EXACT, UNIFORM IN LIFT HEIGHT: THM-820 first puts every hypothetical tight row in a doubling box or the exceptional `a{1,2,4,8}` top-four branch. THM-845 then applies THM-815's longest-safe-component cap separately at every labelled prefix. The exact tree has 40,590 / 612,221 / 111,675 / 7,255 / 9 states at depths one through five; the exceptional branch has 415 / 178 / 1 / 0 states from depth two onward and dies before completion. All nine doubling-branch terminals retain an explicit strict-safe interval, with exact maxima between `1/10` and `2/15`; hence every proper scale-one residue-preserving five-coordinate lift of `[12]` is loose. The 772,543-row certificate is byte-stable at `-O3`/`-O0`, sanitizer-clean, and independently checked by closed-danger-union and Fraction endpoint replays. This does not descend arbitrary AP scale. | THM-815/820/822/845 + exact component/certificate/endpoint replays |
| Hamming-five codec / continuation boundary | FINITE-EXACT EMPTY + PROVED OPERATION LAW: all `C(12,5)2^5=25,344` height-at-most-two packets are loose. `H0=H1` has 3,810 exact-`M`-mixed three-face fibres and is not insertion-Markov: the same labelled insertion splits the canonical liar pair's handoff relation. Literal exact endpoint words make the bounded three-face join injective and update exactly under every monotone comb addition by intersection. They are not deletion-Markov and do not determine exact `M`; a fixed-cardinality common-deletion counterexample forces retention of the labelled tooth bank for reversible/transport operations. | THM-822/840 + exact C++/Fraction replays |
| Hamming-five arbitrary-scale deck boundary | PROVED STRUCTURAL / BOUNDED COMMON-SHEET BANK CLOSED: every scalar cover has `sum 1/D_i>3/13`, hence some `D_i<=21`, but the family `(1,2,3,5,10):(2,5,13q+1,2,13q+1)` refutes any all-order cutoff. THM-823 classifies every common-sheet survivor with all effective orders at most twelve: all order one, one all-order-three orbit, or one-plus-an-order-three quartet. THM-845 closes all order one after normalization, THM-844 closes all 96 all-order-three contexts, and THM-847 closes all 96 proper mixed contexts (the zero-height face is THM-816). Thus the complete bounded bank is empty; common-sheet languages with an order above twelve remain open. | THM-816/823/837/844/845/847 + exact C++/Fraction replays |
| all arbitrary-height order-three Hamming-five contexts | PROVED / FINITE-EXACT EMPTY: THM-844 applies THM-815's longest-component cap to every labelled prefix in all 96 contexts left by THM-823. The exact recursion visits 28,876 states, finds no covering prefix, and dies at depth four with no depth-five terminal; it strictly improves the earlier `K/L` cap at every visited state. THM-837's original 75,371-state single-context calculation becomes a 213-state subcase. Raw and longest-component-marginal comb tournaments are transitive in every context but flip 492 edges; the theorem-bearing carrier is the active component--comb incidence with exact endpoints, remaining CRT progressions, and last speed. | THM-815/823/837/844 + exact Fraction replay |
| all arbitrary-height mixed order-one/order-three Hamming-five contexts | PROVED / FINITE-EXACT EMPTY: after normalization the order-one speed is `3(b+13h)`, with proper lifts exactly `h>=1`, so all five obligations lie in labelled step-39 progressions. THM-847 exhausts the `3*8*4=96` contexts in 31,715 exact states, with no covering prefix or depth-five terminal; a standalone replay and a literal census of all 63,360 marked mixed sheet presentations agree. Both tournament gauges are transitive but flip 395 edges. The faithful carrier is again active component--remaining-comb incidence with endpoints, progression labels, and last speed. | THM-815/816/823/847 + two Fraction replays |
| scale-one Hamming-six finite recursion / seven-comb wall | PROVED finite-decidability, not emptiness: THM-815's discrepancy recursively bounds every next proper scale-one lift through radius six, with initial cap `x_1<=468` at radius six. Seven remaining combs are the precise density barrier for this method because `13-2m` changes sign. Arbitrary five-/six-coordinate deck ramification is not classified. | THM-815 Part C + exact root-component scans |
| dyadic deletion descent from the two-sheet packet | PROVED: every imprimitive-deletion branch is a factor-2 seam and a finite dyadic quotient chain with binary safe-child fibers, a unique first `Z/4` seam, primitive divisor-complete quotients, and a hereditarily primitive terminal base; terminal exclusion remains open | THM-775 |
| ten-even/two-odd locus through speed height 100 | FINITE-EXACT EMPTY; all 1,225 odd-pair bad-atom hypergraphs have transversal number 12, independently regenerated atlas | THM-776 |
| ten-core phase-cell / erosion packet | PROVED anchored/symmetric return packet, pointwise thickness tax, and global gap/Kneser budgets `mu(E)+sum min(g_i,4/(143B))<=mu(H)` and `mu(E)+mu(R)<=mu(H)`. Exact liars show that fixed anchors, raw component tournaments, exception-divisor grids, and signed residue support all lose the escape predicate. THM-803 now constructs the exact all-component selector on `K_U=E_U+closure(R_U)` and a sharp row for which every grid and global maximizer is silent but the nonmaximal component `7/22` escapes. The remaining theorem is uniform failure or incompatibility of those finite selector obligations, not construction of the selector. | THM-782/789/797/803 + exact trap/erosion-liar/anti-grid/component certificates |
| even-maximum two-sheet collar | PROVED rational blocker clock, top-tooth incidence, and a `Z/13` moving-edge carrier. Its exact quotient is an `A_12` root-current walk in the 50,388-state seven-chip simplex, with coverage iff all singleton cut capacities remain nonnegative. The tropical block transfer `T(W)=(c_W,b_W)` composes exactly and preserves survival with the actual initial allocation. FINITE-EXACT/UNIFORM-IN-MULTIPLIER: at quotient height 24, `c=1` tears by `3/8`, `c=3` by `1/7`, and every odd `c>=5` fails in the initial chamber. Uniform quotient-height exclusion remains open. | THM-792 + exact root-current/tropical `w=13c` certificates |
| n=12 sporadic branch | OPEN globally. Closed: the bounded shallow slice; the full AP-centred Hamming-one through Hamming-four stars at arbitrary height and scale; the proper scale-one Hamming-five chart; every effective-order-at-most-twelve common-sheet H5 survivor language (all-one, all-three, and mixed one-plus-three); every two-sheet core in `[1,19]` with unbounded odd exceptions; and the other named finite banks. THM-844/845/847 make the bounded H5 sheet bank empty. THM-840 proves exact endpoints are Markov only for monotone insertion; reversible/transport operations need the labelled tooth bank. Still open are common-sheet H5 languages with an order above twelve and other unbounded deck ramifications, non-AP H6, and the deeper branches. On the deep side THM-824/831 classify the radius scalar and ternary switch boundary, while THM-836 excludes the first two thin owner shells but stops sharply at local `s=5` survivors. Remaining: unbounded-order arbitrary-scale H5 descent; non-AP H6; uniform radius/sum-arc exclusion or transport; the dyadic/collar residual; and higher sheets. | THM-759/763/765/766/768/769/770/772/774/775/776/782/786/789/792/795/797/800/803/804/806/807/810/815/816/817/820/821/822/823/824/831/836/837/840/844/845/847; HYP-6820 |
| tail lanes | EXACT (identity + scans) | THM-750 closed budget; U1 discharged (S283) |

## 2. The Lean ledger

- **LRCClosedBudget.lean: 47 declarations, 0 sorries, all [propext, Classical.choice,
  Quot.sound]** — the (H)-edge chain machine-checked END TO END, geometric to spectral:
  pairOverlap → acorr_eq_model → geometric_disc_eq_discB → grid_deficit / raabe_B2 → discB;
  capped_envelope_kernel → Fourier envelopes → spectral_thm755. THM-731/732/755 have complete
  Lean faces; kps's exact-ℚ certificates have no prose links below the band edge.
- Fleet Lean assets: LRCDetunedDispatch (THM-668), LRCShadowGap (THM-748), the certificate
  supply (THM-693–698), SmallClusterFull, LRC13Citation, kps decide-bottoms.
- The formalized chain proves the stated geometric/Fourier/capped-envelope lemmas. It does **not**
  supply the missing uniform scale-normal classification, so there is no sound assembly theorem
  to wire yet.

## 3. Remaining work items

1. **Prove a scale-normal structure theorem.** Split normalized families into coherent dilation
   packs, additive/hierarchical clusters, and an incoherent residual, while retaining scale residue
   and killer offset. Raw far count and diameter are not invariants of this quotient.
   *Progress (S299--S303/S2 audit): THM-761 proves the exact-budget `r<=6`
   sheet regime (`c>=43` is uniform only for coprime exceptions), and THM-771
   closes the unramified `r=7` lane above its reduced-winding bound. The raw
   bound `r_+(P)<=B(c*)` is false even in exactly `f=4` (THM-798). THM-799
   proves that one high-frequency insertion over a fixed positive-mass base
   cannot cause safe-mass decay and gives a proportional-peel terminal state.
   THM-777 gives the rho bridge and bounded census; THM-780 now supplies the
   unconditional height-free floor `182^(-12)`. Thus the ratio coordinate is
   uniformly bounded, but the open composition must still retain endpoint
   owners, divisor/gcd sidecars, and deck-tiling residue. THM-799 additionally
   closes every fully lacunary far-count flag, so unresolved infinity lies on
   clustered comparable-scale faces. Only the sharp global floor value and
   a practical cutoff remain conjectural.*
2. **Attack persistent translate covers with their metric stalks retained.** At level one,
   `I(13,p,1)` is exactly a 13-translate cover of `F_p^x/{±1}` by the
   strict-danger set. In the tight `s=2` quotient, THM-772 now forces a
   primitive divisor-complete ten-core and THM-774 turns the two odd colours
   into a sharp folded diamond. THM-775 forces every imprimitive-deletion route
   into a finite dyadic quotient chain with binary safe-child fibers, and
   THM-776 excludes the full height-100 packet. THM-789 shows that neither scalar measure nor
   arbitrary fixed-anchor refinement can prove noncontainment: one admissible
   deep anchor traps the full Bohr return set even though another deep time
   escapes.  Its global gap/Kneser budgets tax every tight component layout,
   but an exact two-pair liar shows that raw component order, raw escape signs,
   and total eroded measure still lose the signed tooth slope. THM-797 now
   supplies a uniform selector on odd exception-divisor grids; its q=13
   signed-wall upgrade removes double-13 and full support and pins the sole
   survivor to the exact signed complement of `+/-y`, with `M(U)>=2/13`,
   `x<=2B-1`, and `y<=B-1`.  Its sharp post-wall row silences every
   exception-divisor grid but escapes at denominator 17, so selecting only
   global maximizers or exception divisors is still insufficient. THM-803 adds
   the mandatory half-grid ladder, full parity-twisted support at `q=13`, and
   an exact finite all-component endpoint/cusp selector. Its sharp `U_*` row
   passes all grids and both maximizers but escapes at the nonmaximal singleton
   `7/22`. THM-817 now decomposes `closure(R_U)` into signed maximum-speed
   cells and replaces the coarse selector size by
   `2c_E N_R+2W-2g<=20B^2+22B-2g`.  Its exact family with
   `N_R=3+1440n` proves that the present arithmetic/scalar gates allow linearly
   many satellites; the cell/deep-component incidence must therefore be used
   rather than assuming connectedness. Prove that this adaptive selector fails
   uniformly on every surviving packet,
   or derive a global incompatibility among its owner-labelled obligations;
   alternatively prove a uniform bad-atom transversal
   lower bound, quantitatively land the terminal dyadic base in a certified
   region, or force a transversal/base/effective-order descent. For higher
   sheets, classify which colour covers persist under lifts and evade the
   omit-one gcd reduction.
   *Prime-lens refinement:* THM-778 supplies the full Euclidean endpoint word;
   THM-779/784 rule out raw wall count, and MISTAKE-148 withdraws fixed-index
   de-phasing. After persistent stalks are contracted, THM-788 gives the
   decorated normal form `E_0,V_1,...,V_A,E_A`, where `E` blocks absorb
   free fastest-owner refinement and `V` packets are ordered zero-sum visitors.
   THM-794 proves this is only an intermediate normal form: a full visitor
   packet can repeat as a diagonal deck translation, making `A`, genuine
   switches, and the proposed universal extent unbounded/false at `R=2`.
   THM-802 upgrades that example to an affine pumping theorem: every
   fundamental unequal one-fast-owner count class `(1,...,1,k)`,
   `k=2,...,6`, occurs existentially (with non-effective onset) with
   `Theta(L)` prefix-legal diagonal-return copies inside any prescribed open
   core-safe interval.  First quotient
   packet return maps modulo common sheet translation; then retain the exact
   centered-mechanical owner word, repeated-owner insertion, redundancy roots,
   and metric/core phase cell.  Bare normalized collision SCCs are downstream
   telemetry, not the legality predicate. A continued-
   fraction digit is useful only through this labelled fibre action.
3. **Classify higher deck ramification, finish Hamming six,
   then cross the seven-comb wall.**
   THM-815/816 close every AP-centred Hamming-four packet.  THM-815's
   discrepancy recursion makes the scale-one radius-five and radius-six
   charts finite exact trees, with first caps `146` and `468`. THM-820 gives
   radius five two explicit uniform collar boxes and THM-845 empties both;
   THM-844/847 also empty the bounded all-three and mixed common-sheet
   languages. Classify effective orders above twelve, enumerate radius six
   after separating the genuine `2[12]` equality orbit, and build the
   radius-six analogue of THM-810's oriented deck classification. At
   scale-one radius seven the mean danger density is
   `14/13`, so the same potential ceases to decrease. Seek a replacement
   potential using overlap debt, owner diversity, or signed component/comb
   incidence. The parallel six-exception/seven-exception sheet
   wall is a structural clue, not yet an identification theorem.
4. **Reproduce and formalize finite tiles.** Independently rerun THM-738's complete bank; finish
   rather than promote THM-741; attach machine-verifiable certificates and exact scope metadata.
5. **Hygiene.** Resolve theorem-ID collisions and require the status vocabulary `proved`,
   `finite-exact for stated bank`, `verified sample`, `conditional`, or `open`.

Small-period correction: S312's `120/120` observation cannot be a terminal
uniform lemma.  THM-762 supplies exact residuals first witnessing at `q=26`
and `q=27`, while replay of S105 gives `91/8260` rows with `q_min>25` and a
bank maximum of `38`.  The exact replacement through `q=28` is a
zero-owner/signed-unit-pair blocker deck, not a denominator-only tournament.

Tight-locus correction: the incoming THM-763 uses the zonotopal finite-checking
argument with a strict-interior refinement to prove that every primitive tight
twelve-speed tuple has `sum A<=78^11`.  This makes the sporadic branch honestly
finite, but not computationally exhausted.  THM-765 removes every imprimitive
leave-one-out core uniformly.  THM-766 adds the uniform projective
bound `a_11/a_1>=72/7`; below `12` it places the top speed in one of eleven
danger-tooth cones and retains an exact core-maximizer residue band.  The
remaining equality problem is simultaneous component-tooth/splice-lattice
coherence, not raw height.

Binding-scale refinement: a tight rational maximizer is not described by its
prime residue alone.  THM-769 writes every reduced maximizer as `p/(13s)` and
splits the speeds into an on-sheet core `sU` and off-sheet exceptions.  The
familiar complete nonzero residues are exactly the shallow `s=1` fibre.  In a
deep fibre the exceptions must cover every lift of the entire loose set of
`U`; two exceptions force `s=2` and persistent opposite parity, while the
three-exception equality edge is an `s=3` persistent colour cover.  THM-768
separately rules out a largest speed that is the unique 13-multiple.  THM-770's
exact owner-CSP settles all shallow labelled packets through lift height twelve
(`13^12` conceptual rows), leaving only the permitted dilates and the unique
   primitive AP.  THM-795 proves the complete unbounded one-coordinate direction:
every proper residue-preserving Hamming-one lift of every unit AP dilation is
loose.  THM-800 closes the next stratum uniformly: at exact tightness, oriented
half-open splice-deck capacity forces two replacements to share the AP scale,
and every normalized proper double lift has `M>=2/25`. THM-804 forces exact
tightness at Hamming radius three to common scale, and THM-806 closes that
scale-one base: one replacement lies in `[14,24]`, the others satisfy
`v<=381,w<=12v`, and an exact `5,713,539`-row component sweep has no tight
packet. THM-810 identifies the first ramification at radius four: either
common scale, or the order-three quartic coset `<5>` which normalizes directly
into an `s=3` deep packet. THM-815 closes the common-scale chart twice: by a
recursive safe-component/danger-comb discrepancy ladder and by an independent
collar-cycle reduction with `768,735` exact component-containment rows.
THM-816 closes every lift of the structured three-sheet interface by a finite
dynamic residual-comb recursion. Thus AP-centred Hamming four is fully loose.
THM-820 makes the complete scale-one Hamming-five chart uniformly finite. A
minimal four- or five-cycle gives the doubling box through `x<=388`, while the
exceptional `a{1,2,4,8}` top SCC satisfies `x<=228`, `v<=1986`, and
`max_top<=7944`. THM-845 exhausts both branches by row-wise longest-component
caps: its exceptional branch dies early and all nine terminals are loose. This
is arbitrary-height scale-one closure, not an arbitrary-scale descent. The exact height-one minimum `2[11] union {11}` at
`1/12` reaches the doubled AP `2[12]` after the sixth odd lift, exposing the
adjacent scale-change face. THM-815 also makes scale-one Hamming six recursively
finite with initial cap `468`; its density coefficient changes sign at seven.
The non-AP H6 tree, all-scale deck interfaces, radius at least seven, the
two-sheet selector, and other deep colour covers remain open.

THM-823 now proves why the radius-five all-scale interface cannot repeat the
radius-four scalar argument.  Five-colour scalar capacity has robust infinite
cones; the explicit two-large-order family survives for every `q>=3`, although
its first `D=40` row has no common-sheet unit word and covers at most `29/40`
sheets.  The exact bounded no-order-one bank collapses instead to the forward
flags of the three-coset cycle
`H -> 2H -> 4H -> H`: all five orders are three, four labels form a THM-810
coset, and the fifth lies in the doubled coset.  Every least CRT lift is loose,
but the lift-invariant mod-39 clock supplies only the boundary `1/13` at
arbitrary height.  THM-837 first closed one directed-flag/parity word by a
global `K/L` recursion.  THM-844 replaces that bound by the longest active
component at every prefix and exhausts all 96 words: 28,876 exact states, no
covering prefix, and no depth-five terminal.  Thus this entire bounded
all-order-three survivor language is loose at arbitrary lift height.  The
remaining deck problem is common-sheet and metric, not scalar capacity.
THM-845 also closes the all-order-one/common-scale language after division;
THM-847 closes all 96 proper mixed one-plus-order-three contexts by the same
longest-component recursion, with `h=0` reducing to THM-816.  Therefore every
common-sheet survivor language with all effective orders at most twelve is
loose.  Only unbounded-order common-sheet languages remain on this deck edge.

The two-sheet edge is now considerably narrower.  THM-772 proves that its
ten-speed quotient is primitive, has a multiple of every modulus `2,...,12`,
has no 13-multiple, and bounds both odd exceptions by `11 max(U)`.  THM-774
then folds eligibility and opposite parity into the single exact inequality

```text
||(x+y)tau/2||+||(x-y)tau/2||>=11/13.
```

The permitted diamond has sharp normalized measure `8/117`, attained at
reduced ratio `x:y=9:1`.  This is a necessary filter, not a closure: some
explicit loose cores have smaller measure, so the next theorem must compare
the labelled component locations and widths with the diamond teeth.  That
comparison is now finite-exact for every ten-core `U subset [1,19]`: the widest
component intrinsically caps each odd exception, and all `767,700` permitted
core/odd incidences fail individual coverage.  This closes the whole low-core
slice with odd speeds unbounded, but supplies no global bound on `max(U)`.
THM-775 shows that every failure of hereditary primitivity is not arbitrary:
it is a dyadic seam, giving a literal `2+1+1` partition on the first four
sheets and then a canonical binary assigned-ownership tower with binary
safe-child fibers, ending at a hereditarily primitive quotient. THM-776 reverses the remaining finite
quantifiers: every odd pair through height 100 induces a bad-atom hypergraph
of transversal number 12, too large for the ten-speed quotient core.  The
uniform metric substrate is no longer missing: THM-789 strengthens THM-782 to
a symmetric simultaneous-return packet of measure `2*72^(-10)` and a component
of length at least `72^(-10)/(5 max(U))`, and derives the pointwise thickness
tax and Bohr erosion forced by tightness.  Globally it proves
`mu(E)+sum min(g_i,4/(143B))<=mu(H)` and
`mu(E)+mu(R)<=mu(H)`, with
`mu(R)>=max(4/(143B),2*72^(-10))`. The exact row
`U={1,2,3,5,7,8,9,10,11,12}`, `(x,y)=(13,9)` goes further: its entire return
set `(-1/858,1/858)` is trapped at `t_0=4/17`, while `14/19` is deep and
escapes. Thus even symmetric local incidence at an arbitrary anchor cannot
finish the job.  On the same core, `(13,9)` and `(17,13)` have identical raw
signed component-tournament fingerprints and identical scalar eroded-diamond
measure, but different return-thickened component incidence. The uniform
residual is now smaller: THM-797's odd-divisor shells discharge every row with
a deep grid class outside the opposite-exception shell, and its q=13 signed
walls leave only the exact ten-residue complement of `+/-y`. THM-803 then
forces full parity-twisted support and silence on the universal `26,52,78`
anti-grid ladder, and replaces informal all-component search by an exact
endpoint/cusp selector of quadratic size. Its sharp row passes every grid and
both maximizers but escapes on the nonmaximal singleton `7/22`. The remaining
theorem is therefore uniform failure or incompatibility of this signed-tooth,
owner-labelled **all-component selector**. THM-817 makes the signed-tooth
coordinate literal: each return component is one maximum-speed cell, and the
selector has adaptive size `2c_E N_R+2W-2g`. Its exact family with
`N_R=Theta(B)` proves that current arithmetic/scalar gates cannot force
connected or sublinear return geometry, while an explicit central escape
keeps the family out of the tight branch.  The target is the cell-by-deep-
component incidence (or
quantitative seam-guard bounds that put the reconstructed ten-core/full packet
inside a certified base), not an unstructured search over ten-even/two-odd
tuples.  THM-817 resolves the selector's disconnected-return input exactly:
every component of `closure(R_U)` is one signed cell inside a tooth of
`B=max(U)`, with endpoint owners reconstructed by interval intersection.  If
`N_R` cells survive, the exact selector needs at most
`2c_E N_R+2W-2g<=20B^2+22B-2g` tests.  This is an adaptive representation
theorem, not compactness: the primitive divisor-complete signed-complement
family `B_n=506+360360n` has `N_R=3+1440n`.  The remaining argument must use
the labelled cell/deep-component incidence, not hope that the satellites
become connected or sublinear.

THM-821 resolves that incidence at the atomic boundary uniformly over every
positive odd folded pair `x>y`.  With
`a=(x+y)/2`, `b=(x-y)/2`, each obligation factors through the pair-indexed
circular sum arc `((x,y),C+R_k)`: its minimum occurs at an arc endpoint or a
cusp `j/(2a),j/(2b)`.  The pair cannot be forgotten.  In the fixed
`(13,5)` diagnostic census, the deterministic THM-817 bank plus the complete
three-cell `U_0` row gives `9,974` atoms, with `492` successes and `9,482`
failures. Exact return data alone and exact deep data alone both mix verdicts;
the disconnected satellites make the same deep interval fail for `k=-1` and
succeed for `k=+1`. Exact input intervals or their exact sum arc have zero
mixed fibres. Owners are redundant for one numerical evaluation after the arc
is formed, but remain indispensable ancestry for moving an endpoint under
descent. Thus the next global theorem acts on pair-indexed signed sum arcs with
owner provenance; no component-width, event-type, or bare-cell ranking can
replace them.

THM-824 identifies a different quotient boundary for the **global symmetric
return union**.  The folded target for `(13,5)` is exactly the union of the two
closed balls of radius `2/169` about `5/13` and `8/13`.  A no-switch lemma
shows that, for every nonempty compact deep set `E` and compact symmetric
return set `R=-R` containing zero,

```text
E+R subset H_(13,5)
iff rho_{ {5/13,8/13} }(E)+rho_0(R)<=2/169.
```

For `(13d,5d)` the intrinsic criterion becomes
`rho_(C,d)(E)+rho_d(R)<=2/169`.  If the independent centre-cell condition
`dE subset H_(13,5)` is retained, it has the phase rewrite
`max_E ||13dt||+13 max_R ||dr||<=2/13`; without that condition the rewrite is
false, as `E=R={0}` shows.  The phase inequality by itself is the familiar
necessary thickness tax.  The intrinsic radius criterion evaluates from the
two exact extrema in time linear in the stored components and return cells.  This does
not contradict THM-821's liar fibres: a single signed satellite is not
symmetric, and the global conjunction can forget numerical owner labels only
after the full exact sets have been assembled.  Owners still govern their
transport.  All 214 reconstructed core packets fail both the direct global
test and the radius budget, with zero mismatches, so the result changes the
shape of the remaining proof rather than closing it.  The live target is now
either a uniform arithmetic violation of this budget or a genuinely
satellite-level obstruction where the symmetry collapse is unavailable.

THM-831 classifies that symmetry collapse for every reduced fold.  If
`alpha>beta` are the coprime opposite-parity half-frequencies, the primitive
target has an explicit disjoint Bezout-offset ball decomposition.  Its
symmetric no-switch radius factorization holds exactly for
`4<=alpha<=9`: sixteen two-ball types.  From `alpha=10` onward the component
offsets `q=1,3` carry the gauge-invariant ternary relation behind the switched
triple `(-1|-3,+1)`.  A common gcd adds a raw deck progression and therefore
raw switches, although multiplication by the gcd still gives an exact
quotient-scaled radius theorem.  The method boundary is now proved; what
remains is arithmetic exclusion on the sixteen viable quotient types and a
three-hyperedge-aware argument elsewhere.

THM-829 independently makes one owner-transport mechanism exact.  A primitive
centre column `v` and Bezout owner row `b` move under a unimodular action as
`(v,b)->(Av,bA^{-1})`; reflection conjugates the action, and only
`+/-I,+/-R` commute with the reflection branch.  This confirms that a bare
inverse residue cannot be the cross-denominator state.  It does not yet
identify a THM-817 cell endpoint with a primitive witness column or preserve
the LRC predicate, so the missing construction is a fibre product between
this arithmetic stalk and the metric signed-cell atlas.

On the prime-seven face, THM-778 reconstructs the complete endpoint event
schedule that THM-773's moments forget. THM-779 converts full eight-owner
blocking into an integer token walk; the exact criterion is proved, while its
uniform exit bound remains open. Truth compression, event transport, and
long-run exclusion are therefore three distinct obligations.

THM-792 now clocks the remaining even-maximum collar more sharply.  Its
eleven-speed maximizer has denominator at most `4R-2`; safe mass occupies many
top teeth with repeated disjoint flank owners; and a forced odd `13c` exception
produces a labelled moving edge cover on thirteen quotient sheets.  The entire
height-24 multiplier lane is exact: `c=1` has 101,850 initial covers and tears
by `3/8`; `c=3` has 2,528 and tears by `1/7`; every odd `c>=5` fails in the
initial chamber because `13c>2 max(U)`.  This is uniform in the multiplier, not
in quotient height.

With initial excess `e^0` and cumulative grouped root current `C`,
coverage is precisely `e^0+C>=0` coordinatewise, and
`K=K_0+<e^0,C>+||C||^2/2`.  Thus the finite state is the 50,388-point
seven-chip simplex plus a dead state.  Labels can be compiled into the
`A_12` root word, and its tropical transfer `T(W)=(c_W,b_W)` composes exactly,
but neither word nor actual initial allocation can be discarded.  The uniform
theorem is an intersection problem between a regular safe-current language and
the thin arithmetic language generated by divisor-complete rational clocks.
THM-802 supplies the precise warning from the prime-seven carrier: affine
phase cells can pump prefix-legal return loops.  The collar target must rule out
or classify arithmetic-realizable safe tropical returns up to `Z/13` sheet
rotation, not merely show that the abstract current graph is finite.

The collar overlap degrees carry the exact energy
`K=sum_j binom(d_j-1,2)`.  A simple event satisfies
`Delta(-K)=d_departure-d_entry-1`, the THM-785 endpoint-defect flux, and `8K`
has THM-787's step-eight increments.  This is a useful current coordinate but
not a quotient: two certified height-24 cores have the same indexed degree
vector and `8K=80` yet first tear at different events and sheets.

THM-789 changes the folded residual in a complementary direction.  Symmetric
return packets improve the safe-measure and component-width floors, and exact
gap/Kneser erosion budgets constrain the whole component layout.  Its `U_0`
example proves that the full local Bohr return set can remain trapped at one
deep time even though another escapes.  Its component-tournament liar then
proves that raw margin order/sign and total eroded measure still do not encode
which component escapes: the next folded theorem must transport signed affine
tooth addresses or eroded margins while selecting globally.

The exact max-peel tooth atlas rules out a tempting shortcut.  In the
exhaustive slice `A subset [1,20]`, `M(A\{w})>1/12`, `M(A)<=1/10`, all
`2,453` escaping rows still have winding one.  Pure endpoint ownership and
even a transitive safe-component phase tournament with one Hamiltonian path
have explicit liars.  Any anti-cover theorem must retain exact midpoint-width
slack and owner incidence, not just winding or phase order.

THM-768--770 sharpen the arithmetic split.  A tight twelve-set cannot have
its maximum as the unique multiple of `13`.  The no-multiple (shallow) branch
is exactly a full nonzero residue transversal; every multiple-owner (deep)
binding point has scale `s>=2` and at least two off-sheet tighteners satisfying
THM-769's exact sheet-capacity inequality.  THM-770 exhausts all `13^12`
shallow packets through lift height twelve by a lossless `24,008`-cell
endpoint-owner CSP: the only zero-defect rows are `c*{1,...,12}` for
`c=1,...,12,14`, hence the only primitive row is the AP.  This closes the
shallow sporadic slice through `max A<=168`. THM-795 additionally closes all
higher shallow lifts at labelled Hamming distance one from any AP dilation;
THM-800 closes every proper residue-preserving labelled Hamming-two lift and
gives the sharp normalized floor `2/25`. THM-804/806 close the complete
Hamming-three star at arbitrary height. THM-810 then splits the first live
Hamming-four chart into common-scale packets and one exact order-three quartic
coset interface with the `s=3` deep branch. THM-815 closes the common-scale
alternative by two independent exact reductions, and THM-816 closes the coset,
so every AP-centred residue-preserving Hamming-four packet is loose. Their
shared discrepancy state also makes scale-one radii five and six finite exact
decision trees. THM-820 sharpens radius five to a recursive doubling box and a
metrically capped top-four SCC box; THM-845 closes both boxes at arbitrary lift
height. THM-844 then transports the same longest-component principle through
the complete bounded all-order-three common-sheet orbit, closing all 96
arbitrary-height contexts. THM-847 closes the parallel 96-context mixed orbit,
so THM-823's complete effective-order-at-most-twelve bank is empty. Two
labelled five-cycles with identical live/tie
fingerprints have different integer centres and maxima, and all THM-844 binary
tournaments remain transitive despite 492 gauge flips. Thus the
residue-obligation, collar-cycle, and component-order tournaments are
telemetry: the proof-facing state is the operation-indexed residual interval
union with its remaining labelled tooth bank, endpoint owners, and exact
widths. Seven is the first scale-one radius where mean danger density exceeds
one. Unbounded-order H5 decks, non-AP H6, deep two-sheet, and general
higher-sheet branches remain.

For the historical viewpoint atlas behind this synthesis, use
`00-navigation/LRC-LENS-MAP.md`; for the current predicate-preserving carrier,
use HYP-6820 section F and LTT-434--436.  The lens map is navigation, not a
claim that every listed analogy has theorem status.

## 4. The one-line frontier

> **LRC(14) is open.  The AP-centred Hamming-four star, proper scale-one
> Hamming-five chart, and entire effective-order-at-most-twelve Hamming-five
> common-sheet bank are closed; scale-one Hamming six is finite but contains
> the doubled-AP equality face.  The live obstruction is unbounded-order
> arbitrary-scale deck ramification plus signed-cell/component compatibility
> in the deep and higher-sheet branches.**

*Controlling corrections: HYP-6780, MISTAKE-143, MISTAKE-149, THM-762/764,
THM-768--770, THM-794/795/797/800--804/806/807/810/815--817/820--824/831/836/837/840/844/845/847,
and the updated THM-758. Earlier
S297/S310/S312 closure language and the companion S297 reflection must be read through these corrections.*
