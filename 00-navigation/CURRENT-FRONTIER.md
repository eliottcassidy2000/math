# Current Mathematical Frontier

**Rolling state — refreshed 2026-07-21.** This file supersedes dated frontier
snapshots as a statement of present status. A linked theorem is the proof
source; this document records how the pieces compose and what remains.

Status vocabulary:

- **PROVED:** a complete in-repo proof or exact reduction is recorded.
- **CITED:** an external theorem is imported with its hypotheses and provenance.
- **FINITE-EXACT:** exhaustive over the stated finite universe.
- **VERIFIED:** reproducible evidence, not a general proof.
- **CONDITIONAL:** implication depends on a named open input.
- **OPEN / REFUTED / SUPERSEDED:** literal meanings.

## LRC(14)

### Headline

**OPEN.** In the standard reduction there are 13 distinct nonzero relative
speeds, corresponding to 14 total runners. External work settles at most 13
total runners. The residual is not a missing routine finite run and not a
uniform small-period lemma.

The two most useful decompositions of the remaining obstruction are:

```text
global LRC(14)
  -> finite-height / spectrum / pinning shell
  -> compact-covering structural shell
       -> extract a tight or AP-like 12-speed core
       -> use the proved post-extraction far-element exit
```

and

```text
safe certificate
  = seed + selector + preserver + pointwise exit.
```

We have strong selectors and preservers. The genuinely missing objects are the
structural supplier/seed and a lossless exit.

### What is closed

- **CITED:** LRC through 13 total runners. The newest step is
  [Sungkawichai–Trakulthongchai](../05-knowledge/reference/CORE-PAPERS.md#eleven-twelve-and-thirteen-lonely-runners),
  a computer-assisted proof for 10–12 nonzero speeds.
- **PROVED:** [THM-763](../01-canon/theorems/THM-763-strict-finite-height-for-tight-lrc-instances.md)
  gives every primitive tight 12-speed tuple the finite bound
  `sum(a_i) <= 78^11`, conditional only on the settled lower-dimensional LRC
  citation used in its descent.
- **PROVED in its locus:**
  [THM-1171](../01-canon/theorems/THM-1171-twelve-term-ap-tight-rigidity.md)
  handles twelve-term AP tightness. It is not a classification of arbitrary
  twelve-speed tight tuples.
- **PROVED / FINITE-EXACT in stated strata:**
  [THM-1284](../01-canon/theorems/THM-1284-crossN-firstgap-band-and-singlefar-classification.md)
  closes the first-gap, single-far `N=12` stratum; the AP-centered common-scale
  H6 program is closed through `c=36`, with `c=37` prime-excluded and `c=38`
  next.
- **PROVED-BY-EXHAUSTION:**
  [THM-1290](../01-canon/theorems/THM-1290-subgap-exhaustive-census-bounded-height.md),
  after the full MISTAKE-194 rerun, proves LRC(14) for maximum speed at most 55
  and emptiness of `(1/14,3/41)` through maximum speed 64.
- **CITED + proved translation:**
  [THM-1289](../01-canon/theorems/THM-1289-floor-isolated-from-above-by-GK.md)
  imports the Giri–Kravitz one-sided accumulation theorem: the floor `1/14` is
  isolated from above by some ineffective `delta>0`. Whole-window finiteness is
  not a consequence of the published theorem; it needs their Conjecture 1.5 or
  another input.
- **PROVED with scope:**
  [THM-1291](../01-canon/theorems/THM-1291-cf-active-leg-law-proved.md)
  makes the first beating integer a continued-fraction convergent denominator.
  Identifying it with an active speed needs its additional hypothesis H.
- **REFUTED:**
  [THM-1288](../01-canon/theorems/THM-1288-c71-refuted-divisor-aligned-clusters.md)
  refutes Sungkawichai–Trakulthongchai Conjecture 7.1 literally.
- **PROVED:** [THM-819](../01-canon/theorems/THM-819-primitive-harmonic-law-for-interval-cores.md)
  gives the primitive harmonic good-measure law.

### What is false or exhausted

- **REFUTED:** uniform `q <= 25`. [THM-762](../01-canon/theorems/THM-762-small-denominator-signed-pair-deck.md)
  and [THM-764](../01-canon/theorems/THM-764-covering-small-period-signed-pair-deck-and-q25-refutation.md)
  give `26*{1,...,12} union {339}` with first good period `q=27`.
- **REFUTED method:** averaging a good-period count. MISTAKE-127/129/130 show
  that existence is governed by the maximum over periods; the resonant tight AP
  is the mandatory hostile control.
- **CORRECTED search frame:** HYP-8815/MISTAKE-221 gives a proved necessary
  kernel: after gcd normalization, a counterexample is primitive Cover14 with
  `M<1/14`, and its maximum-deletion core is non-AP, hence has Schur-triple
  count at most `65` by THM-1017/730. Since its lonely-time measure is zero,
  THM-731 also forces `disc_v>=6|G'_{~v}|^2` for every peel; small discrepancy
  is a safety signal. THM-2048 sharpens this by a discrete fiber tax: with
  `mu=|G'_{~v}|`, `theta={7vmu}`, and `r_v` interval components, every
  zero-measure packet must satisfy
  `6(vmu)^2+theta(1-theta)/7<=r_v^2/3` for every peel.
  This is a strict improvement over the old uniform tail test: for the
  primitive Cover14 row `{1,8,11,12,14,17,22,26,35,40,54,90,93}` at `v=93`,
  the old inequality is inconclusive while the integer tax violates the new
  necessary inequality and forces a positive-measure lonely interval.
  “Near-AP,”
  “anti-golden,” Fibonacci-foil, CF-blocker, and
  full-autocorrelation iff claims remain heuristic. THM-1002 proves every
  maximizer lies on a pair-sum ruler `t=p/(v_i+v_j)`. Enumerating all such
  rulers makes the corrected S206 replay exact: all fifteen displayed rows are
  safe, while `{1,...,12,5460}` proves the old `q<=1200` cutoff was incomplete
  (`92/1197 < 420/5461`). This finite bank does not identify the global
  covering minimum.
- **REFUTED arrangement shortcut / exact replacement:** MISTAKE-223 separates
  the THM-1820 Fourier annihilator from a standard toric-complement/OS layer
  formula. Shi counts and bounded short-relation counts do not preserve LRC.
  THM-2047's signed phase-height complex is lossless; for `delta>0`,
  `chi(G_delta)=#components`, so it sees isolated tight witnesses that volume
  misses.
- **PROVED finite-circuit alternative:** THM-2051 closes the coarse analytic
  “middle” sought by THM-935/946. THM-965 gives the sharp universal pair-
  covariance floor `delta_(a,b)>=-6/637`, so pairs can be paid exactly. If a
  thirteen-speed row has no exact support-three-through-five relation with all
  coefficient magnitudes at most `2^20`, whole-product Fejer--BV approximation
  gives continuous `BONF5>0` and therefore a positive-measure strict lonely
  set. The proof avoids the unresolved termwise `T4/T5` strip/slab bounds. The
  LRC residual is now localized to finitely many genuine higher-relation
  templates, but descent or classification within those hyperplanes—including
  the relation-rich AP boundary—remains open.
- **PROVED rank descent and explicit residual atlas:** THM-2052 combines
  THM-763's finite primitive-height bound with a pigeonhole/code argument.
  Every hypothetical counterexample has eleven independent support-at-most-three
  relations of height `91^6`, hence is rank-twelve finite-box or lies in one of
  finitely many two-anchor stars `a_i v_p+b_i v_q+c_i v_i=0`, `c_i!=0`.
  Modulo dilation each star has one projective parameter; pair clusters and
  cross-cluster triples give the same split. THM-2053 proves
  `max_i|a z_i-b u_i|<=(a^2+b^2)/91`, removing `||(a,b)||>=91L`, while adjacent
  normalized columns raise the parameter torus elementarily to margin `1/13`.
  Every row obeys `M(v)>=1/13-R/(2N)`, `N|(v_i+v_j)`; the exact deck `D_N(m)`
  removes longitudinal fibers, with a divisibility down-set and AP cutoff `156`.
  THM-2055 restricts failure to hull-owner sectors while non-hull deck sidecars
  remain essential. THM-2056 Kelvin-inverts them to `(1/91)R^{-1}K^o` and proves
  `2u.v>=A_p(u)+A_p(v)` for whole acute unimodular cones. THM-2057 closes two
  entire one-tail planes: `{a,2a,...,11a,13a,w}` by a `12a` / `14a` / `84a`
  certificate, and `{a,2a,...,12a,w}` by a `13a` / `14a` / `182a` certificate.
  More generally, a missing clock `N<=14` forces `Na|w` over `aC`. Despite
  `640690` primitive determinant failures, its arithmetic carrier has three
  leaves. THM-2059 makes arbitrary-`N` core/tail composition an exact CRT
  histogram dot product; its zero-mode split isolates disjoint-support zeros.
  THM-2060 gives every tail bin the sharp floor `q-ceil(q/7)` unless `a|w`
  and proves the multi-tail capacity law. THM-2064 independently repackages
  that law on one common finite fiber: a safe core clock lifts whenever
  `sum_w ceil(t_w/7)/t_w<1`. Both isolate the sole primitive
  imprimitive-core two-tail capacity tie `2C union {x,y}`, with `x,y` odd.
  THM-2061 folds that exact seam, pins its core through `2..14`, bounds its
  tails, and closes `max(C)<=19`. THM-2066 attaches the exact binary
  owner word to every safe clock: a counterexample needs complementary tail
  words. The clock bank `15..34` empties that complement graph for all `59,880`
  primitive divisor-complete quotient cores in `{1,...,24}`, extending the
  finite closure to `max(C)<=24`. THM-2062 puts an exact hereditary CRT wheel
  on THM-2058's coprime intervals.
- **PROVED relative whole-product decorrelation:** [THM-2054](../01-canon/theorems/THM-2054-relative-fejer-whole-product-decorrelation.md)
  says bounded scalar relations must lift to vector-character relations. Two
  telescopes bound the whole product; rowwise `H=2^19` clears every recorded
  pinned-base margin. MISTAKE-080/082 still require shape-specific torus models.
- **Scope separation:** [THM-1149](../01-canon/theorems/THM-1149-compact-essential-crown-and-farey-blocker.md)
  proves that tight deletion and an all-loose essential crown are different
  branches. Equality classification after a tight deletion cannot be applied
  before extracting that deletion.
- **Local-comb ceiling:** THM-1252 through THM-1274 close or sharply saturate
  most purely local six-comb return arguments. The live residue is global
  endpoint/child transport or a phase-located turn tax, not another unlocated
  local-return charge.

### Exact live residuals

1. **AP-core / tight-deletion supplier (highest leverage).**
   [THM-1017](../01-canon/theorems/THM-1017-ap-core-bridge-reduction.md) proves
   `AP core -> far element -> LRC(14)`. The extraction of that AP core from the
   compact Cover14 residual is open (the current HYP-7310 form). Non-AP and
   deeper multi-defect twelve-speed branches are precisely what uniform
   “sporadic emptiness” has not removed.
2. **Six-comb global transport.** Complete endpoint/child transport, or prove a
   phase-aware turn tax. A local return lemma without a global location sidecar
   repeats an exhausted pattern.
3. **Effective height/spectrum bridge.** Extend exact pinning beyond the current
   maximum-height shell or make the Giri–Kravitz isolated gap effective. A raw
   larger census is valuable only if its new survivors/filters expose structure.
4. **Frobenius-safe certificate.**
   [THM-2041](../01-canon/theorems/THM-2041-frobenius-stability-of-exact-period-projectors.md)
   supplies whole-packet preservation. The ranked
   [HYP-8800](../05-knowledge/hypotheses/HYP-8800-lrc14-face-carry-frobenius-transfer.md)
   routes are: THM-671 `B5` supply; familywise Fejer/Toeplitz rigidity;
   characteristic-3 endpoint-labelled period-14 propagation; an adaptive
   resolved-phase sheaf; and conductor plus owner-current glue. Raw Ramanujan
   energy is diagnostic only.
5. **Signed wall-word / Euler deletion route.**
   [THM-2047](../01-canon/theorems/THM-2047-phase-height-toric-arrangement-for-lrc.md)
   consolidates the earlier threshold tope/Cech interfaces into a
   general-height carrier. Opposite signs recover pair-sum rulers; deletion is
   exact; and the top-vertex boundary-layer coefficient depends only on the
   extreme active slopes. The live target is a deletion--restriction or nerve
   argument forcing `chi(G_{1/14})>0`; arrangement cohomology by itself does
   not provide that positivity. [THM-2050](../01-canon/theorems/THM-2050-period14-top-germs-do-not-determine-global-loneliness.md)
   sharpens the no-go: AP13 and `12->26` have identical complete local germs at
   all six unit points but maxima `1/14` and `1/12`. The exact
   [HYP-8841](../05-knowledge/hypotheses/HYP-8841-lrc14-noetherian-first-exit-termination.md)
   audit shows THM-2048's peel tax fires on deep `12->96` and covering
   `12->84`, and independently reproduces the genuine Cover14 gain at peel
   `93`, but misses `12->26`, `12->36`, and P10+K33. It is a pruning rule, not
   the termination height. Once `M>1/14`, THM-2047 bounds a first strict
   exit by `q<=2 max(S)`; the sole Wall-A branch is owner-labelled Euler-point
   survival when no positive pair-sum margin occurs.
   [THM-2051](../01-canon/theorems/THM-2051-fejer-bv-small-relation-alternative-for-lrc14.md)
   now proves that every remaining hard row has a genuine support-`3..5`
   relation of coefficient height at most `2^20` after triangle-packing the
   exact pair covariances.
   [THM-2052](../01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md)
   goes much further: every hypothetical counterexample already has eleven
   independent bounded three-support relations. It either reaches the
   rank-twelve finite box or belongs to a finite two-anchor,
   one-projective-parameter star atlas. The exact HYP-8841 bank also shows why
   incidence is load-bearing: all eight controls, tight and strict alike,
   already have a height-one three-term circuit (the first seven share
   `1+2=3`).
   [THM-2053](../01-canon/theorems/THM-2053-rank-two-parameter-plane-geodesic-terminal.md)
   supplies two complementary filters before rank twelve: its exact deck
   discards whole fibers, while THM-2055 prunes determinant failure to hull-
   owner sectors. THM-2058 then enumerates their intersection as labelled
   primitive packets and coprime intervals. HYP-2108 supplies the
   endpoint gate inside those cells: for core-safe
   components with midpoint `m_i` and length `l_i`, prove for some peel `w`
   that `P_w=max_i(||w m_i||+(w/2)l_i-1/14)>=0`. Equality already preserves a
   weak endpoint; strict positivity opens an interval. HYP-2896 is the exact
   one-parameter model, and THM-2054 is the analytic complement on vector-
   character cells whose bounded scalar resonances all lift. The sharpened
   target is **Euler survival, rank gain, off-resonance decorrelation, or exact
   deck/fan discharge**.
   HYP-8845 halves the covering branch by mirror symmetry: one survivor on
   `[0,1/2]` automatically supplies its partner and `chi>=2`.

The characteristic-7 result
[THM-2043](../01-canon/theorems/THM-2043-period14-parity-hasse-jet-completeness.md),
is **PROVED, but is a carrier no-go rather than LRC(14)**. It identifies
`F_7[C_14]` with the two length-seven local rings at `+1` and `-1`. THM-2058's
GW-based family keeps that packet, blockedness through `q=13`, threshold `14`,
any prescribed finite lift depth, and the same one-tail hull-owner sector, yet
has the uniform strict phase `47/113` with integer slack `13`. Thus raw jets,
threshold, fan owner, and fixed height precision are jointly insufficient.
Exact labelled phase packets or an adaptive resolved certificate are the
honest positive carrier; the
smallest useful atlas now visibly includes the incompatible scales
`{14,27,41}`.

### Fresh perspective prompts

Do not default to runners or arcs as vertices. Test gaps, fixed sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier
packets, matroid circuits, proof obligations, and survivor states. Ask whether
the hard object is a set of speeds, a deletion/core extension, a phase-current
process, a spectrum stratum, or a finite-state certificate. Every quotient must
name the LRC predicate it preserves and the coordinate a sidecar must restore.

## NC2 and Gaussian moments

### Headline

- **PROVED:** [THM-2022](../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md)
  proves NC2 and hence unrestricted GMC(2).
- **REFUTED for higher dimensions:** explicit GMC(3) counterexamples now exist;
  GMC is false for every dimension at least 3. See the current external preprint
  entry in [`CORE-PAPERS.md`](../05-knowledge/reference/CORE-PAPERS.md).
- **FORMALIZATION PARTIAL:** the former descent and face-construction gaps are
  closed, alongside the Wick, normalized-moment, dilation, Kummer/Lucas,
  Frobenius, and abstract three-case residue layers. DvdK remains an explicit
  cited premise. The remaining internal step is the concrete normalized-channel
  instantiation of the residue assembly and the final conditional theorem
  `DvdK1 -> NC2`; the aggregator does not yet assert end-to-end `NC2`.

### Why THM-2022 works

For a fixed exact support, algebraic torus descent reduces a hypothetical
complex null point to a number field. The lowest balanced Wick face has a
nonzero complete constant term by the one-variable Duistermaat–van der Kallen
theorem. At a good prime, Kummer kills non-dilated channels, strict face height
kills dilated off-face channels, and Lucas plus Frobenius leaves the full face
residue `Q^p`. The success is **whole-layer preservation**, not atomwise
separation.

### Boundaries and repaired bridges

- [THM-2040](../01-canon/theorems/THM-2040-the-de-factorialization-principle.md)
  is now a retired pointer to THM-2022/2033, not an independent theorem. The
  exact symmetric-wall calculation survives in the historical proof record;
  the general principle is prime-local initial form, not a global
  common-factorial or Vandermonde factorization.
- [THM-2033](../01-canon/theorems/THM-2033-the-nc2-wall-is-the-confluent-transitivity-vandermonde.md)
  is a valid determinant/Vandermonde identity for a special moment matrix, not
  for the general scalar moment.
- [THM-2023](../01-canon/theorems/THM-2023-gmc2-hyperbessel-boundary-laguerre-polya.md)
  proves the boundary hyper-Bessel function is Laguerre–Pólya; it is now a
  refinement rather than the missing NC2 step.
- MISTAKE-211/212/214/215 bar four recurring overclaims: scalar moments do not
  separate atoms; a tournament source is sufficient, not iff; Vandermonde node
  repetition is not tournament score repetition; and no global
  de-factorialization exists.

### Live work

Package THM-2022 for publication and formalization; generalize the scalar-grade
prime-block criterion (positive rational Gamma radial shapes are already in
scope); identify where vector-valued Wick grades require a different exposed
face; and transfer the four-gate `seed/selector/preserver/exit` design without
pretending the preserver supplies the seed.

## Tournaments

### Reliable structural toolkit

- [THM-1805](../01-canon/theorems/THM-1805-the-vandermonde-is-a-signed-tournament-sum-intransitivity-cancels.md):
  Vandermonde is a signed tournament sum; transitive score permutations survive
  and directed 3-cycles are the fundamental cancellation atom.
- [THM-1862](../01-canon/theorems/THM-1862-order-join-reduction-principle.md):
  under order-join, Hamiltonian-path count `H` is multiplicative and directed
  triangle count `c3` is additive.
- [THM-1936](../01-canon/theorems/THM-1936-signed-redei-join-multiplicative.md):
  signed Rédei data is join-multiplicative.
- [THM-1880](../01-canon/theorems/THM-1880-the-a-b-functional-frame-chebyshev-pell-companions.md)
  and [THM-1885](../01-canon/theorems/THM-1885-the-half-and-plus-one-monoid-is-baumslag-solitar.md):
  the transitive skew-characteristic recurrence is Chebyshev–Pell-like under
  `x -> x+1` and `x -> x/2`, producing `BS(1,2)^+`.
- [THM-1926](../01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md):
  the Bowen–Lanford zeta factors over the strong core.
- [THM-1940](../01-canon/theorems/THM-1940-var-lambda2-is-a-4-subtournament-census-invariant.md):
  `var(lambda^2)` is a census of four-vertex subtournaments.
- [THM-1965](../01-canon/theorems/THM-1965-the-tournament-invariant-lattice-mapped.md):
  the mapped invariant lattice is exact through `n <= 6` only.
- [THM-1966](../01-canon/theorems/THM-1966-signed-redei-count-independent-invariant-n7.md):
  at `n=7`, signed Rédei magnitude is independent from spectrum plus `H`.
- [THM-2013](../01-canon/theorems/THM-2013-coordinates-for-the-continuum-cyclic-temperature.md)
  and [THM-2016](../01-canon/theorems/THM-2016-the-deep-continuum-and-the-reducibility-ceiling.md):
  cyclic temperature, the reducibility ceiling, and the invariant-resistant hot
  center are the best current decomposition. `H` is an empirical thermometer,
  not a proved coordinate of every shell.

### Live work and limits

- [THM-1950](../01-canon/theorems/THM-1950-h-ge-disc-reduced-to-strongly-connected.md)
  reduces `H >= disc` to strongly connected tournaments and verifies finite
  cases; the global inequality is open.
- Determine operation-response laws before proposing new scalar invariants.
  Work on substitution, join factorization, SCC cores, switching, duality, and
  minimal independent coordinates at `n>=7`.
- A literal Paley bridge exists only in THM-640's prime quadratic-residue scope.
  Composite modulus 14 and repeated score/node heuristics do not inherit it.
- For applications outside tournament theory, first prove that the pairwise
  relation is intrinsic and target-preserving. A forced total orientation can
  destroy exactly the ties or magnitudes the original problem needs.

## Integer sequences

### Headline

[THM-2000](../01-canon/theorems/THM-2000-support-harmonic-abel-dini-figurate-surface.md)
and [THM-2005](../01-canon/theorems/THM-2005-support-dirichlet-automatic-tournament-atlas.md)
replace “the reciprocal sum of a sequence” with two explicit objects:

```text
support profile       D_A(z) = sum_{a in support} a^{-z}
indexed multiplicity  = D_A(z) + collision tax.
```

At `z=1`, Abel–Stieltjes and logarithmic block occupancy give the exact
convergence criterion. Density zero is insufficient; the critical family
continues through every Bertrand/iterated-log scale. Egyptian splitting
conserves mass exactly at `z=1`.

### Established landmarks

- Ordinary and centered polygonal supports have digamma clocks; square and
  centered-octagonal confluences lead to trigamma forms.
- The maximum-cyclic-triangle support mass is exactly
  `75/4 - 24*log(2)`.
- The condensation-hazard support is `{5,6,...}` with profile
  `zeta(s)-H_4^(s)` and prefix-product parity shuffle tax
  `67/4 - 24*log(2)`.
- The A000568 support mass has a certified bracket narrower than
  `3.11e-44`.
- [THM-2010](../01-canon/theorems/THM-2010-new-tournament-invariant-sequences.md)
  records four-term candidates and a no-match search. That is not proof of
  sequence novelty.
- **EXACT identities / heuristic bridge:** corrected HYP-8820 places caterer
  and cake as Pascal-prefix sums, Moser and bagel as fixed binomial row
  functionals, and Fibonacci as a shallow diagonal. Algebraically,
  `bagel(n)=cake(n+1)-2`, hence `bagel(n)-cake(n)=T_n-1`. Klein-S313's
  full-rank gap diagonals exactly realize the g-bonacci kernels; finite-rank
  shadows eventually differ. MISTAKE-222 blocks the stronger claim that a
  shared array or matching `-1` identifies the torus and shadow boundaries or
  transfers an LRC/JC predicate. The live bridge test is an explicit common
  boundary-cell/Euler-characteristic valuation, not another prefix match.

### Live work

Classify profiles under support operations rather than compare only their value
at one; study analytic continuation/abscissae and automatic/Mahler structure;
track collision taxes for census sequences; and feed support identities back
into tournament operation laws and LRC residue packets.

## Other active portfolio

- **Jacobian/Dixmier:**
  [THM-1300](../01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md)
  exactly verifies an owner-supplied three-dimensional Keller map and its triple
  collision; [THM-1315](../01-canon/theorems/THM-1315-keller-counterexample-surjective-etale-3to1-caustic.md)
  analyzes its geometry. Public provenance remains unsettled in canon;
  MISTAKE-205 forbids the old Alpöge–Mathew attribution. Keep exact verification,
  public announcement, and discovery credit separate.
  The linked external certificate independently replays the resulting `A_3`
  cotangent pullback; it does not decide `DC(2)`.
  [THM-2044](../01-canon/theorems/THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension.md)
  disproves the two-pair Poisson conjecture, while
  [THM-2046](../01-canon/theorems/THM-2046-first-order-cotangent-pullbacks-cannot-cross-the-DC2-wall.md)
  proves that this witness cannot descend through multiplication positions and
  first-order dual momenta. [HYP-8803](../05-knowledge/hypotheses/HYP-8803-A3-pullback-versus-A2-quantum-descent.md)
  localizes the remaining nonfiltered gate at extension across `x=0` in an
  exact Ore-Weyl chart. [THM-2049](../01-canon/theorems/THM-2049-the-DC2-Ore-boundary-correction-complex-is-acyclic.md)
  proves its simultaneous `(S,T)` associated-graded correction map is onto in
  every relevant beta grade, kills the grade-six obstruction guess, and gives
  a formal beta-adic lift. Finite polynomial termination and coupling to `D`
  are now the decisive gates. THM-2063 proves a complete planar normal form:
  if any nonzero output-pencil member is affine along any linear source fiber,
  the Keller map is tame with an explicit inverse. Thus a hypothetical JC(2)
  counterexample has fiber degree at least two for every source direction and
  every output-pencil direction. MISTAKE-228 relabels THM-1330 correctly as a
  necessary-structure atlas, not a classification. MISTAKE-229 keeps the valid
  binary symmetric-Hessian subcase but blocks the unproved NC2/GMC-to-JC chain
  and the claim that VC(4), leading-form/Jelonek descent, and
  Lame-for-polygons are equivalent. `DC(2)` and planar `JC(2)` remain open.
- **Gaussian higher dimensions:** THM-1490 is one verified higher-dimensional
  construction; newer three-real-Gaussian examples supersede any claim that
  dimension four is sharp.
- **Full portfolio:** use [`PROBLEM-LEDGER.md`](PROBLEM-LEDGER.md). Every session
  should keep at least one niche or wildcard item outside its anchor, especially
  when an operation, obstruction, or sidecar transfers across domains.

## Cross-domain connection discipline

The most reusable current bridges are not literal object identifications:

| Mechanism | Proven source | Legitimate transfer question |
|---|---|---|
| Whole-layer Frobenius / orbit norm | NC2 balanced face and TNC monodromy orbit | Exact finite-abelian packet preservation and uniform-incidence norms transfer; LRC still lacks the seed and exit, and `-1` in the packet stabilizer forces a signed orientation sidecar. |
| Quotient + sidecar | spectrum/SCC/support quotients | Which lost coordinate restores the target predicate? |
| Bulk / boundary / null | THM-2058 phase height | Strict templates occur on all large prime grids; tight support lies on finitely many level-14 pair-sum clocks; subthreshold packets vanish. |
| Operation-response | tournament joins/support unions | Which observables add, multiply, localize, or collide? |
| One-sided accumulation | LRC spectrum | Can qualitative isolation be made effective by finite state/pinning data? |
| Support Dirichlet profile | integer sequences | Do tournament/LRC census supports have functional or automatic structure? |

HYP-8810's new “JC(2) and LRC(14) share AP-rigidity” reflection is a useful
**wildcard frame, not a proved common reduction**. The LRC side has the precise
one-way THM-1017 supplier. The claimed planar-JC continued-fraction/AP wall must
still be stated as an exact theorem with a map and preserved predicate before
the two residuals can be called identical. Its safe present use is as a prompt:
compare how AP/transitive/one-sided extremals sit in different nullcones, and
test where the functionals destroy the analogy.

A proposed bridge is mathematical only after it names its map, preserved
predicate, information loss, sidecar, and a falsifying control.
