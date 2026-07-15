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
| q=13 anti-grid / all-component selector | PROVED: every odd exception divisor has a mandatory half-grid; at `q=13` this forces full parity-twisted support, and the complete universal even anti-grid ladder is exactly `d=2,4,6`. The full erosion predicate is equivalent to finitely many owner-labelled component endpoints and folded cusps, at most `200B^2+22B`. A sharp signed-complement row passes every divisor/anti-grid and global-maximizer test but escapes at the nonmaximal singleton `7/22`. Uniform failure of the selector remains open. | THM-803 + exact 71,673-half-grid / 4,096-support replay |
| shallow full-residue locus through lift height 12 | FINITE-EXACT over `13^12` conceptual packets; 13 dilates, unique primitive row `{1..12}` | THM-770 + exact owner-CSP |
| shallow AP Hamming-one star at arbitrary height | PROVED scale-free: every proper residue-preserving one-coordinate lift of `{1,...,12}`, and of every unit AP dilation, is loose. Every residual non-AP shallow packet differs from every AP dilation in at least two labelled coordinates. | THM-795 + exact core-threshold/atom/deck replay |
| shallow AP Hamming-two star at arbitrary height | PROVED scale-free: oriented half-open splice decks force both replacements back to the common AP scale at exact tightness; the normalized proper double lift then has the sharp floor `M>=2/25>1/13`. Thus every residual shallow packet has labelled Hamming radius at least three from every AP dilation. | THM-800 + exact 600,756-row replay |
| shallow AP Hamming-three scale descent | PROVED scale-free: if a residue-preserving three-replacement lift is tight, all three deck orders are one and all replacements share the AP scale. Thus every off-scale Hamming-three branch descends to a genuine scale-one triple lift; looseness of that normalized chart remains open. | THM-804 + exact 296,640-grid / 11,143,660-capacity replay |
| scale-one Hamming-three collar attack | CLAIMED / AUDIT IN PROGRESS: scratch work proposes one lift `<=24`, then bounds the others by `v<=381`, `w<=12v`; the component-containment proof and corrected full replay are unfinished. No closure may use this row yet. | THM-806 reservation |
| dyadic deletion descent from the two-sheet packet | PROVED: every imprimitive-deletion branch is a factor-2 seam and a finite dyadic quotient chain with binary safe-child fibers, a unique first `Z/4` seam, primitive divisor-complete quotients, and a hereditarily primitive terminal base; terminal exclusion remains open | THM-775 |
| ten-even/two-odd locus through speed height 100 | FINITE-EXACT EMPTY; all 1,225 odd-pair bad-atom hypergraphs have transversal number 12, independently regenerated atlas | THM-776 |
| ten-core phase-cell / erosion packet | PROVED anchored/symmetric return packet, pointwise thickness tax, and global gap/Kneser budgets `mu(E)+sum min(g_i,4/(143B))<=mu(H)` and `mu(E)+mu(R)<=mu(H)`. Exact liars show that fixed anchors, raw component tournaments, exception-divisor grids, and signed residue support all lose the escape predicate. THM-803 now constructs the exact all-component selector on `K_U=E_U+closure(R_U)` and a sharp row for which every grid and global maximizer is silent but the nonmaximal component `7/22` escapes. The remaining theorem is uniform failure or incompatibility of those finite selector obligations, not construction of the selector. | THM-782/789/797/803 + exact trap/erosion-liar/anti-grid/component certificates |
| even-maximum two-sheet collar | PROVED rational blocker clock, top-tooth incidence, and a `Z/13` moving-edge carrier. Its exact quotient is an `A_12` root-current walk in the 50,388-state seven-chip simplex, with coverage iff all singleton cut capacities remain nonnegative. The tropical block transfer `T(W)=(c_W,b_W)` composes exactly and preserves survival with the actual initial allocation. FINITE-EXACT/UNIFORM-IN-MULTIPLIER: at quotient height 24, `c=1` tears by `3/8`, `c=3` by `1/7`, and every odd `c>=5` fails in the initial chamber. Uniform quotient-height exclusion remains open. | THM-792 + exact root-current/tropical `w=13c` certificates |
| n=12 sporadic branch | OPEN globally. Closed: the bounded shallow slice; the full Hamming-one and Hamming-two stars at arbitrary height; every off-scale residue-preserving Hamming-three branch by descent; every two-sheet core in `[1,19]` with unbounded odd exceptions; the full two-sheet speed box through 100; every non-signed-complement `q=13` profile; and the forced-`w=13c` quotient box through height 24. Remaining: genuine scale-one Hamming-three lifts (THM-806 only claimed) and radius at least four; uniform failure of THM-803's exact signed-complement selector; an unbounded-height collar tear; and higher-sheet packets. | THM-759/763/765/766/768/769/770/772/774/775/776/782/786/789/792/795/797/800/803/804; HYP-6820 |
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
   `7/22`. Prove that this selector fails uniformly on every surviving packet,
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
3. **Reproduce and formalize finite tiles.** Independently rerun THM-738's complete bank; finish
   rather than promote THM-741; attach machine-verifiable certificates and exact scope metadata.
4. **Hygiene.** Resolve theorem-ID collisions and require the status vocabulary `proved`,
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
and every normalized proper double lift has `M>=2/25`. THM-804 now proves the
next descent: exact tightness at Hamming radius three forces
all three replacements to share the AP scale.  The genuine scale-one
Hamming-three chart remains open, as do radius at least four and the deep
colour-cover branch.

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
owner-labelled **all-component selector** (or
quantitative seam-guard bounds that put the reconstructed ten-core/full packet
inside a certified base), not an unstructured search over ten-even/two-odd
tuples.

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
gives the sharp normalized floor `2/25`. THM-804 proves common-scale descent
for three replacements, leaving genuine scale-one Hamming-three lifts, radius
at least four, and every deep branch. The residue-obligation and
sheet-margin tournaments are transitive telemetry; the proofs live in the
endpoint-owner hypergraph and the missing-splice sheet-danger deck with
oriented core-safe germs.

## 4. The one-line frontier

> **LRC(14) is open. The proved `f<=3` tile and the per-core capped envelope are substantial,
> but the `f>=4` branch requires a scale-normal structural classification, not raw enumeration.**

*Controlling corrections: HYP-6780, MISTAKE-143, MISTAKE-149, THM-762/764,
THM-768--770, THM-794/795/797/800--804,
and the updated THM-758. Earlier
S297/S310/S312 closure language and the companion S297 reflection must be read through these corrections.*
