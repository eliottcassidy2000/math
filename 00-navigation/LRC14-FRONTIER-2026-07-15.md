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
| prime-lens r=8 blocking chain | Full blocking is exactly piece surjectivity + wall rainbow + no simultaneous walls; the redundancy cocycle gives an exact event-word test and rooted states form an `A8` torsor. Raw wall length is unbounded: a divisor-complete persistent seven-stalk family has `2m` covered walls, and an independent simpler family has `21N`. The live target is the stalk-quotiented owner-switch/visitor skeleton. | THM-779/784; MISTAKE-147; exact refuter/A8 audit |
| r=8 visitor/extent laws | Phi recurrence, period-sum, single-visitor break, cluster balance, no-companion extent pierce, signed visitor-set difference, and the factor-two fixed-companion span are proved. If the six companion speeds sum to less than `g`, an explicit extent bound follows; THM-788 contracts empty fastest periods to active normal form. The old factor-one span and `sum c<f` threshold are false/unsupported; general alternating companions remain open. | THM-783/786/788; MISTAKE-147/148 |
| eight-owner `c=7` buffer rigidity | The chamber/rainbow condition is subsumed by THM-779's integer token walk. The global `1/7+O(gcd/w)` partner-buffer law remains VERIFIED but must be restricted to the core-safe exit set; whole-core blocking and the general-density stalk-quotiented switch bound remain OPEN. | HYP-6840 + THM-779/783/784/786/788 |
| seven-owner deck defect / ramified residue | Exact identity `F=Q+Omega-sigma`; exact tilings are chamber-locked, KCL necessity is WITHDRAWN, and mirror coincidence is diagnostic. Primitive `c=21` row realizes `(0,12,12,0)` | THM-771 + corrected THM-767 + exact audits |
| raw positive-length fragmentation bound `r_+(P)<=B(c*)` | REFUTED by an exact four-far family and an independent census; peel-relative `rho` is measured at most `9.335` only on the stated bank | THM-794; HYP-6830 correction; MISTAKE-145 |
| safe-measure floor / normalized band bridge | `rho<=12/(pi|G'_P|)` and `|G'_P|>=1/(91 maxP)` PROVED; phase pigeonhole gives the height-free floor `182^(-12)`; exact `maxP<=18` floor is `7/858`, unique at `{1,...,13}\{6}`, while that sharp global value remains CONJECTURAL | THM-777/780 |
| positive good-set state `(mu,r_top)` + frequency `N` + proportional peel `aN` | PROVED conservative transition `mu_N>=6mu/7-2r_top/(7N)`, `r_top,N<=N+r_top`, hence eventually capped; fully lacunary flags are terminal in every far-count stratum, with factors `412,405,394,27,17,14,13,13,13,13`; an exact 2,002-core floor gives a complementary factor-19 cone on `f=4` | THM-795; THM-794 is a full-family instance |
| primitive tight 12-speed locus | UNIFORMLY FINITE (`sum A<=78^11`), not classified | THM-763 |
| hereditary primitivity of tight 12-sets | PROVED; every leave-one-out core is primitive | THM-765 |
| unique-largest-13-multiple tight branch | IMPOSSIBLE by explicit prime-grid perturbation | THM-768 |
| arbitrary binding scale `p/(13s)` | PROVED packet split; shallow iff full nonzero residues; deep exceptions obey exact sheet capacity, with `s=2` parity and `s=3` colour criteria | THM-769 |
| two-/three-sheet equality quotients | PROVED primitive divisor transfer and speed bounds; two-sheet core contains multiples of every `2,...,12`, three-sheet core of every `2,...,11` | THM-772 |
| two-sheet metric residual | PROVED exact folded diamond `||(x+y)tau/2||+||(x-y)tau/2||>=11/13`; sharp measure cap `8/117`; all quotient cores in `[1,19]` closed against unbounded odd exceptions | THM-774 + exact certificate |
| shallow full-residue locus through lift height 12 | FINITE-EXACT over `13^12` conceptual packets; 13 dilates, unique primitive row `{1..12}` | THM-770 + exact owner-CSP |
| dyadic deletion descent from the two-sheet packet | PROVED: every imprimitive-deletion branch is a finite dyadic quotient chain with binary safe-child fibers, a unique first `Z/4` seam, primitive divisor-complete quotients, and a hereditarily primitive terminal base; terminal exclusion remains open | THM-775 |
| ten-even/two-odd locus through speed height 100 | FINITE-EXACT EMPTY; all 1,225 odd-pair bad-atom hypergraphs have transversal number 12, independently regenerated atlas | THM-776 |
| ten-core phase-cell / erosion packet | PROVED anchored return packet; symmetrization gives safe measure at least `2*72^(-10)` and a component of length at least `72^(-10)/(5 max(U))`; Bohr erosion and a pointwise odd-exception thickness tax hold.  An exact core traps every natural local refinement at `t=4/17` but escapes at `14/19`, so global deep-component selection—not fixed-anchor refinement—is the residual. | THM-782/789 + exact trapping certificate |
| even-maximum two-sheet collar | PROVED rational blocker clock `q<=4R-2`, center/effective-order dichotomy, uniform occupied top teeth with repeated disjoint flank types, and a `Z/13` moving-edge carrier.  FINITE-EXACT: forced `w=13` has no packet for `U subset [1,24]` (1,144,066 cores; zero full-word survivors).  Uniform tear remains open. | THM-792 + exact 117-event-group certificate |
| n=12 sporadic branch | OPEN globally; bounded shallow slice, every two-sheet core in `[1,19]` with unbounded odd exceptions, the full two-sheet speed box through 100, and the forced-`w=13` quotient box through 24 are empty; unbounded shallow descent, global eroded-packet noncontainment, uniform collar tear, and higher-sheet packets remain | THM-759/763/765/766/768/769/770/772/774/775/776/782/786/789/792; HYP-6820 |
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
   bound `r_+(P)<=B(c*)` is false even in exactly `f=4` (THM-794). THM-795
   proves that one high-frequency insertion over a fixed positive-mass base
   cannot cause safe-mass decay and gives a proportional-peel terminal state.
   THM-777 gives the rho bridge and bounded census; THM-780 now supplies the
   unconditional height-free floor `182^(-12)`. Thus the ratio coordinate is
   uniformly bounded, but the open composition must still retain endpoint
   owners, divisor/gcd sidecars, and deck-tiling residue. THM-795 additionally
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
   escapes. Prove the global erosion failure
   `E_U not subset H_(x,y) minus R_U`, using component/owner incidence to select
   the right deep component; alternatively prove a uniform bad-atom transversal
   lower bound or quantitatively land the terminal dyadic base in a certified
   region. For higher
   sheets, classify which colour covers persist under lifts and evade the
   omit-one gcd reduction.
   *Prime-lens refinement:* THM-778 supplies the full Euclidean endpoint word;
   THM-779/784 rule out raw wall count, and MISTAKE-148 withdraws fixed-index
   de-phasing. After persistent stalks are contracted, THM-788 replaces the
   word by decorated normal form `E_0,V_1,...,V_A,E_A`, where `E` blocks absorb
   free fastest-owner refinement and `V` packets are ordered zero-sum visitors.
   Bound `A`, retain the redundancy roots, then intersect the resulting metric
   interval with the core-safe base. A continued-fraction digit is useful only
   through this labelled fibre action.
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
primitive AP.  Neither result proves an unbounded shallow descent or eliminates
the deep colour-cover branch.

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
tax and Bohr erosion forced by tightness. The exact row
`U={1,2,3,5,7,8,9,10,11,12}`, `(x,y)=(13,9)` goes further: its entire return
set `(-1/858,1/858)` is trapped at `t_0=4/17`, while `14/19` is deep and
escapes. Thus even symmetric local incidence at an arbitrary anchor cannot
finish the job. The uniform residual is now a global deep-component
selection/noncontainment theorem (or
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
produces a labelled moving edge cover on thirteen quotient sheets.  For `c=1`
and `U subset [1,24]`, an exact sweep of all 1,144,066 cores and 117 grouped
event times leaves zero survivors.  This is a new finite slice, not the missing
uniform automaton tear.

The collar overlap degrees carry the exact energy
`K=sum_j binom(d_j-1,2)`.  A simple event satisfies
`Delta(-K)=d_departure-d_entry-1`, the THM-785 endpoint-defect flux, and `8K`
has THM-787's step-eight increments.  This is a useful current coordinate but
not a quotient: two certified height-24 cores have the same indexed degree
vector and `8K=80` yet first tear at different events and sheets.

THM-789 changes the folded residual in a complementary direction.  Symmetric
return packets improve the safe-measure and component-width floors, while its
exact `U_0` example proves that the full local Bohr return set and every
literal phase refinement can remain trapped at one deep time even though a
different deep time escapes.  The next folded theorem must compare or select
global deep components; further refinement of a fixed anchor is insufficient.

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
shallow sporadic slice through `max A<=168`, not the higher shallow lifts or
any deep branch.  Its residue-obligation tournament has all `66` burdens tied;
the theorem lives in the owner-incidence hypergraph, not the transitive tie
quotient.

## 4. The one-line frontier

> **LRC(14) is open. The proved `f<=3` tile and the per-core capped envelope are substantial,
> but the `f>=4` branch requires a scale-normal structural classification, not raw enumeration.**

*Controlling corrections: HYP-6780, MISTAKE-143, THM-762/764, THM-768--770,
and the updated THM-758. Earlier
S297/S310/S312 closure language and the companion S297 reflection must be read through these corrections.*
