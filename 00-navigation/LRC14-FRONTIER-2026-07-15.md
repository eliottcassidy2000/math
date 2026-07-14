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
>   fixed. Under dilation `P -> cP`, the good-set measure is unchanged, its component count is
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
(union wall structural; tight case = cyclic tilings of Z_c; wall realized at c = 7 by a family
that is still lonely via a non-sheet route), (ii) c ≤ 42 small-scale criterion failures (band
inflation bounded there — the finite regime), (iii) gcd-descent bookkeeping, (iv) families with
no scale structure at all — the capped envelope's natural domain, pending the
peel-relative safe-measure/endpoint replacement for HYP-6830's now-refuted
raw fragmentation⟺divisibility complementarity.

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
| (13−r)-speed common-factor core + r ≤ 6 exceptions, any gcds | PROVED above exact thresholds (c ≥ 43 uniform; exact per-(r,c) sets; Σg_a budget) | THM-761 (opus-S299) + battery |
| r = 7 unramified deck stratum (`7\mid c/g_a` for every owner) | PROVED above reduced-shape bound `max_a(w_a/g_a)>=7 max(P)`; the S299 wall is closed at a switching time | THM-771; corrected core of THM-767; MISTAKE-146 |
| prime lens `c=7`, any unramified owner count | Exact token polynomial: coverage iff `X^7-X` divides `product(X-k_a)`; seven-owner exact states map to all 25 masks at heptagon node `n7-a267`; any covered `r=8` wall is a simple event with a seven-owner heptagon stalk | THM-773 + exact 5,040-state/3,003-profile audits |
| prime-lens endpoint transport | Pairwise midpoint clocks are centered mechanical words with an Euclidean parity cocycle; centered Beatty ranks reconstruct every simultaneous wall and drive the exact `F_7` skew product.  The named r=8 row has 10 simple covered walls with palindromic owner word `162,108,108,206,197,197,206,108,108,162` | THM-778 + 6,400-pair/five-movie exact audit |
| prime-lens r=8 blocking chain | Full blocking is exactly piece surjectivity + wall rainbow + no simultaneous walls; consecutive wall owners must follow the collision-pair hop chain.  Adversarial census finds runs through 5, but a universal exit bound remains open.  THM-778's named ten wall hits are all isolated runs. | THM-779 proved criterion; `K0=5` verified, not universal |
| seven-owner deck defect / ramified residue | Exact identity `F=Q+Omega-sigma`; exact tilings are chamber-locked, KCL necessity is WITHDRAWN, and mirror coincidence is diagnostic. Primitive `c=21` row realizes `(0,12,12,0)` | THM-771 + corrected THM-767 + exact audits |
| raw fragmentation bound r_P ≤ B(c*) | REFUTED twice (exact falsifier + census); surviving peel-relative invariant ρ = v*/maxP measured ≤ 9.335, extremal at {1..12} | HYP-6830 correction; MISTAKE-145 |
| safe-measure floor / normalized band bridge | `rho<=12/(pi|G'_P|)` and `|G'_P|>=1/(91 maxP)` PROVED; exact `maxP<=18` floor is `7/858`, unique at `{1,...,13}\{6}`; the same global floor is CONJECTURAL | THM-777 |
| primitive tight 12-speed locus | UNIFORMLY FINITE (`sum A<=78^11`), not classified | THM-763 |
| hereditary primitivity of tight 12-sets | PROVED; every leave-one-out core is primitive | THM-765 |
| unique-largest-13-multiple tight branch | IMPOSSIBLE by explicit prime-grid perturbation | THM-768 |
| arbitrary binding scale `p/(13s)` | PROVED packet split; shallow iff full nonzero residues; deep exceptions obey exact sheet capacity, with `s=2` parity and `s=3` colour criteria | THM-769 |
| two-/three-sheet equality quotients | PROVED primitive divisor transfer and speed bounds; two-sheet core contains multiples of every `2,...,12`, three-sheet core of every `2,...,11` | THM-772 |
| two-sheet metric residual | PROVED exact folded diamond `||(x+y)tau/2||+||(x-y)tau/2||>=11/13`; sharp measure cap `8/117`; all quotient cores in `[1,19]` closed against unbounded odd exceptions | THM-774 + exact certificate |
| shallow full-residue locus through lift height 12 | FINITE-EXACT over `13^12` conceptual packets; 13 dilates, unique primitive row `{1..12}` | THM-770 + exact owner-CSP |
| two-sheet deletion recursion | PROVED every imprimitive deletion is a factor-2 seam; exact first `Z/4` tiling and finite dyadic descent to a hereditarily primitive divisor-complete core | THM-775 |
| ten-even/two-odd locus through speed height 100 | FINITE-EXACT EMPTY; all 1,225 odd-pair bad-atom hypergraphs have transversal number 12, independently regenerated atlas | THM-776 |
| n=12 sporadic branch | OPEN globally; bounded shallow slice, every two-sheet core in `[1,19]` with unbounded odd exceptions, and the full two-sheet speed box through 100 are empty; unbounded shallow descent, uniform folded noncoverage, and higher-sheet packets remain | THM-759/763/765/766/768--772/774--776; HYP-6820 |
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
   *Progress (S299/S300/S3 audit): THM-761 proves its exact sheet-budget regime and THM-767
   closes compatible r=7 decks above `w_max>7 sum(P)` by an event witness. The raw bound
   `r_P≤B(c*)` is false; `P_N={1,...,11,N}` has `c*=1` and `r_P≥N/77-O(1)`.
   The surviving endgame coordinate is peel-relative boundary intensity
   `r_P/(max(P)|G'_P|)`, with endpoint owners and divisor/gcd sidecars retained.
   THM-777 proves `rho<=12/(pi|G'_P|)`, the height-decaying Lipschitz floor,
   and the exact shape floor `7/858` through `maxP<=18`; its global asymptotic
   floor is still a named conjecture, not a compactness theorem.*
2. **Attack persistent translate covers with their metric stalks retained.** At level one,
   `I(13,p,1)` is exactly a 13-translate cover of `F_p^x/{±1}` by the
   strict-danger set.  In the tight `s=2` quotient, THM-772 now forces a
   primitive divisor-complete ten-core and THM-774 turns the two odd colours
   into a sharp folded diamond.  Prove that no such diamond contains the whole
   loose-component word, or force a dyadic/effective-order descent.  For higher
   sheets, classify which colour covers persist under lifts and evade the
   omit-one gcd reduction.
   *Prime-lens refinement:* THM-778 supplies the full Euclidean endpoint word.
   On the named eight-owner row, attach `n7-a267` masks and redundancy roots to
   the ten covered walls and seek a descent on the five-wall half word.  A
   continued-fraction digit is useful only through its labelled fibre action.
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
sheets and then a canonical binary assigned-ownership tower ending at a
hereditarily primitive quotient.  THM-776 reverses the remaining finite
quantifiers: every odd pair through height 100 induces a bad-atom hypergraph
of transversal number 12, too large for the ten-speed quotient core.  The
uniform residual is now a scale-free transversal/noncontainment theorem (or a
descent of the dyadic terminal core into that finite base), not an unstructured
search over ten-even/two-odd tuples.

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
