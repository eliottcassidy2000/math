# Core papers: imported results, consumers, and guardrails

> **Freshness:** checked **2026-07-21** unless dated later; recheck **PREPRINT/RADAR/SEMINAR ONLY** before priority claims.

## Fast frontier snapshot

- **LRC:** fourteen runners is open; April 2026 reaches twelve nonzero speeds.
- **Moments:** [THM-2022](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md) gives GMC exactly through 2; [THM-2801](../../01-canon/theorems/THM-2801-sharp-special-image-boundary-and-beta-shift-witness.md) gives SIC only at 1; JC-bearing restrictions remain.
- **Tournament citations:** 2412.10572 Irving--Omar; 2307.05569 Grinberg--Stanley; 2406.09697 Klanderman et al.
- **Flows:** [THM-2177](../../01-canon/theorems/THM-2177-planar-counterexample-to-goemans-unsplittable-cost-flow-conjecture.md) refutes the exact-cost conjecture; cost-free `d_max` and planar `2d_max` survive.
- **Reciprocals:** support differs from indexed multiplicity; figurate formulas omit repo extensions.
- **JC/DC:** THM-1300 refutes `JC(n>=3)`; `JC(2),DC(1),DC(2)` are open. THM-2071 is local; THM-2044 refutes the two-pair Poisson claim.
- **Hopf preprints:** [Brendle--Hung `S2 x S2`](CORE-PAPERS-BRENDLE-HUNG-S2XS2-2026-08-24.md) and [Hopf/S6](CORE-PAPERS-HOPF-S6-2026-08-24.md) are distinct **CLAIMS / UNDER AUDIT**; THM-3990/3991/3993 prove only local obstructions.
- **Kakeya `R^4`:** **OPEN**; `3.0543` is PREPRINT v1. See the [source/status and THM-4035 boundary audit](CORE-PAPERS-KAKEYA-4D-2026-08-24.md).
- **Apex cubic:** [2608.22870v1](https://arxiv.org/abs/2608.22870) is a computer-assisted **PREPRINT UNDER AUDIT**; see [audit](CORE-PAPERS-APEX-CUBIC-2026-08-25.md).
- **Higher-order truth:** [Ye--Xu preprint audit](CORE-PAPERS-HEYTING-HIGHER-ORDER-TRUTH-2026-08-29.md).
## Rule 30

- **Sources/status:** [announcement](https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/); [active prizes](https://rule30prize.org/) checked 2026-08-15; openness is repo inference.
- **Imported:** [Rowland](https://doi.org/10.25088/ComplexSystems.16.3.239) gives `v_m`; [Christol et al.](https://www.numdam.org/articles/10.24033/bsmf.1926/) applies only to [THM-3471](../../01-canon/theorems/THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum.md) finite sums.
- **Repo frontier:** [THM-4204](../../01-canon/theorems/THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation.md) gives inverse fibres; [THM-4210](../../01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md) gives the forward current tree. Physical closure and prizes remain open.

## Lonely Runner Conjecture

### Sungkawichai--Trakulthongchai — *Eleven, twelve, and thirteen lonely runners*

- **Primary / freshness:** [arXiv:2604.23906v1](https://arxiv.org/abs/2604.23906), submitted 2026-04-26. **PREPRINT v1; not peer reviewed as of this check.**
- **Imported role:** a computer-assisted divisibility sieve proves LRC for `k=10,11,12` nonzero speeds; it also gives a finite-field polynomial shortcut when `k+1` and a sufficiently large auxiliary modulus are odd primes, and identifies the initial improper-tuple computation as the extension bottleneck.
- **Repo consumers:** [LRC14 frontier](../../00-navigation/LRC14-FRONTIER-2026-07-15.md), [finite-check feasibility ledger](../../00-navigation/LRC14-FINITE-CHECK-FEASIBILITY-LEDGER-2026-07-19.md), [THM-1288, literal Conjecture 7.1 refutation](../../01-canon/theorems/THM-1288-c71-refuted-divisor-aligned-clusters.md), and [THM-523](../../01-canon/theorems/THM-523-lrc14-singular-series-proof-skeleton.md).
- **Does not prove:** `k=13`, the standard fourteen-runner case. Conjecture 7.1 is false as literally written by THM-1288's exact family; do not use it as an inverse theorem. The prime shortcut also does not cross the composite `k+1=14` lift wall.

### Vaaler — *Some extremal functions in Fourier analysis*

- **Primary / freshness:** J. D. Vaaler, *Bulletin of the American Mathematical
  Society (N.S.)* **12** (1985), 183--216, Theorem 19. **PUBLISHED / stable.**
- **Imported role:** supplies the one-dimensional degree-`H` trigonometric
  majorant/minorant sandwich for an interval that THM-2085 tensors with signed
  coordinate defects.
- **Repo consumer:** [THM-2085](../../01-canon/theorems/THM-2085-explicit-height-57-rank-seven-selberg-gate.md).
- **Does not prove:** the relative-Hunter inequality, the signed tensor
  bookkeeping, `H=57`, optimality of that height, or LRC(14). Those are
  repo-derived arguments and constants.

### Ungar — *2N noncollinear points determine at least 2N directions*

- **Primary / freshness:** Peter Ungar, *Journal of Combinatorial Theory,
  Series A* **33** (1982), 343–347,
  [DOI 10.1016/0097-3165(82)90045-0](https://doi.org/10.1016/0097-3165(82)90045-0).
  **PUBLISHED / stable; bibliographic record checked 2026-07-21.**
- **Imported role:** historical alternate lens only.  Applying the
  even-cardinality direction bound to the symmetric signed-column
  configuration can supply a nonradial secant whose perpendicular projection
  has a repeated absolute speed.  The current proof of THM-2053 does **not**
  depend on Ungar: its adjacent-normalized-column construction produces the
  full-support repeat projection elementarily.  The standing lower-dimensional
  LRC citation, not this paper, is the remaining external input to the torus
  floor `M_T>=1/13`.
- **Repo consumers:**
  [THM-2053, rank-two geodesic terminal](../../01-canon/theorems/THM-2053-rank-two-parameter-plane-geodesic-terminal.md),
  [HYP-8846, finite tangent-disk completion](../hypotheses/HYP-8846-lrc14-pointed-plane-transport.md).
- **Does not prove:** LRC(14), the repeat-projection lemma now used in
  THM-2053, the determinant gate, necessity of that gate, or emptiness of any
  tangent disk. The anisotropic estimate and disk identity are separate
  in-repo arguments; disk membership remains only an uncertified case.

### Malikiosis--Santos--Schymura — *Linearly-exponential checking is enough...*

- **Primary / freshness:** [arXiv:2411.06903v2](https://arxiv.org/abs/2411.06903),
  published as *Forum of Mathematics, Sigma* **13** (2025), e164,
  [DOI 10.1017/fms.2025.10107](https://doi.org/10.1017/fms.2025.10107).
- **Imported role:** reduces LRC up to `n+1` runners to an explicit finite
  integer-velocity range of order
  `binom(n+1,2)^(n-1) <= n^(2n)`, dramatically improving Tao's earlier
  `n^(O(n^2))` range.  The zonotope formulation is the standard finite-check
  backend used in the repo.
- **Repo consumers:**
  [finite-check feasibility ledger](../../00-navigation/LRC14-FINITE-CHECK-FEASIBILITY-LEDGER-2026-07-19.md),
  [LRC14 frontier](../../00-navigation/LRC14-FRONTIER-2026-07-15.md),
  [THM-599 torus-band theorem](../../01-canon/theorems/THM-599-torus-band-theorem.md).
- **Does not prove:** LRC(14), a feasible small practical search, or the shifted
  theorem unconditionally in every dimension; the paper's shifted statement
  retains its stated Lonely Vector Problem dependency.
### 2026-08-23 primary pins: lonely-runner polyhedra and Khinchin flatness
- Beck--Hoşten--Schymura [arXiv:1606.01783](https://arxiv.org/abs/1606.01783) pins the line-plus-box formulation; Codenotti--Freyer [arXiv:2307.09429](https://arxiv.org/abs/2307.09429) pins integral-dual width; Averkov--Hofscheier--Nill [arXiv:1911.03511](https://arxiv.org/abs/1911.03511) records Barvinok's bound used by [THM-3743](../../01-canon/theorems/THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction.md). None proves LRC(14), sharpness/sparsity, or recursive slicing.
- [Banaszczyk (1993)](https://doi.org/10.1007/BF01445125) supplies THM-4009's input `2 mu(L) lambda_1(L*)<=d`; the inball/sidecars are local and LRC(14) is open.

### Giri--Kravitz — *The structure of Lonely Runner spectra*

- **Primary / freshness:** [arXiv:2304.01462v4](https://arxiv.org/abs/2304.01462),
  published in *Mathematical Proceedings of the Cambridge Philosophical
  Society*, [DOI 10.1017/S0305004125101497](https://doi.org/10.1017/S0305004125101497).
  Version 3 was withdrawn; version 4 says “Fixed error from previous version.”
- **Imported role:** the pinned use is published Theorem 1.4: the relevant
  subtorus spectra have only upper accumulation points.  Under the exact
  translation `D(T_v)=1/2-M(v)`, this proves that every attained maximum-
  loneliness value is isolated immediately above; in particular there is an
  ineffective `delta>0` with no thirteen-speed value in
  `(1/14,1/14+delta)`.
- **Repo consumers:**
  [THM-1289, the pinned citation and translation](../../01-canon/theorems/THM-1289-floor-isolated-from-above-by-GK.md),
  [THM-1290 bounded census](../../01-canon/theorems/THM-1290-subgap-exhaustive-census-bounded-height.md),
  [near-miss ledger](../../00-navigation/LRC14-NEAR-MISS-LEDGER-AND-SANDWICH-2026-07-19-klein-S319.md).
- **Does not prove:** an effective `delta`, finiteness of the full
  `(1/14,3/41]` window, or LRC(14).  The equality of the entire accumulation
  set with the preceding spectrum was withdrawn and is Conjecture 1.5 in the
  audited v4/published text.  **Do not rely on the shorter arXiv landing-page
  abstract for this distinction; use the PDF-level audit in THM-1289.**

### Fan--Sun — *Amending the Lonely Runner Spectrum Conjecture*

- **Primary / freshness:** [arXiv:2306.10417v2](https://arxiv.org/abs/2306.10417),
  revised 2026-01-30 and published in *Electronic Journal of Combinatorics*
  **33** (2026), P1.38, [DOI 10.37236/13840](https://doi.org/10.37236/13840).
- **Imported role:** gives an infinite `n=4` counterexample family to Kravitz's
  original spectrum conjecture and replaces the single rung `s/(ns+1)` by a
  wider `s/(ns+k)` stratification.  The repo identifies this `k` with its
  near-floor slack coordinate and uses the paper's gcd arguments as a model.
- **Repo consumers:**
  [near-miss ledger, Sections 4--5](../../00-navigation/LRC14-NEAR-MISS-LEDGER-AND-SANDWICH-2026-07-19-klein-S319.md),
  [LRC14 proof map](../../00-navigation/LRC14-PROOF-MAP.md),
  [formal spectrum-window arithmetic](../../04-computation/lean/TournamentH7/TournamentH7/LRCSpectrumWindow.lean).
- **Does not prove:** the amended conjecture for general `n`, bound the order
  `k` for thirteen speeds, or make the LRC(14) sub-gap finite.

### Kravitz — *Barely lonely runners and very lonely runners*

- **Primary / freshness:** [arXiv:1912.06034](https://arxiv.org/abs/1912.06034),
  published in *Combinatorial Theory* **1** (2021),
  [DOI 10.5070/C61055383](https://doi.org/10.5070/C61055383).
- **Imported role:** introduces the maximum-loneliness spectrum and the rigid
  rung family `s/(ns+1)`, proves the proposed sharpened statement for `n<=3`,
  and treats specified “one speed much faster” cases.  The rung/pair-sum ruler
  is a basic coordinate system for the repo's first-window work.
- **Repo consumers:**
  [investigation backlog, spectrum section](../../00-navigation/INVESTIGATION-BACKLOG.md),
  [near-miss ledger](../../00-navigation/LRC14-NEAR-MISS-LEDGER-AND-SANDWICH-2026-07-19-klein-S319.md),
  [formal spectrum-window arithmetic](../../04-computation/lean/TournamentH7/TournamentH7/LRCSpectrumWindow.lean).
- **Does not prove:** the proposed spectrum description for all `n`; Fan--Sun
  refute it already at `n=4`.  A Kravitz rung is a candidate value, not a
  classification of all near-floor LRC(14) families.

### Tao — *Some remarks on the lonely runner conjecture*

- **Primary / freshness:** [arXiv:1701.02048v4](https://arxiv.org/abs/1701.02048),
  published in *Contributions to Discrete Mathematics* **13** (2018).
- **Imported role:** supplies the first general finite-checking reduction
  (`n^(O(n^2))`-scale integer speeds), a large-`n` improvement over the union
  bound, and structural estimates for velocities of size `O(n)`.
- **Repo consumers:**
  [THM-523](../../01-canon/theorems/THM-523-lrc14-singular-series-proof-skeleton.md),
  [finite-check feasibility ledger](../../00-navigation/LRC14-FINITE-CHECK-FEASIBILITY-LEDGER-2026-07-19.md),
  [THM-610 covering/deep-hiding dichotomy](../../01-canon/theorems/THM-610-covering-deep-hiding-dichotomy.md).
- **Does not prove:** LRC at any new fixed frontier case, a small usable bound
  for LRC(14), or a uniform lower bound on the **measure** of lonely times.

### Perarnau--Serra — *The Lonely Runner Conjecture turns 60*

- **Primary / freshness:** [arXiv:2409.20160v3](https://arxiv.org/abs/2409.20160),
  revised 2025-08-12. **SURVEY.**
- **Imported role:** the default panoramic map for formulations, equivalences,
  historical results, tight instances, applications, and open problems through
  its revision date.
- **Repo consumers:**
  [THM-523](../../01-canon/theorems/THM-523-lrc14-singular-series-proof-skeleton.md),
  [THM-759 tight-ratio bound](../../01-canon/theorems/THM-759-tight-instance-ratio-bound.md),
  [LRC14 finish map](../../00-navigation/LRC14-FINISH-MAP-2026-07-11.md).
- **Does not prove:** a new runner case or tight-locus classification.  It
  predates the April and July 2026 papers in this map, so it is not by itself a
  current-frontier certificate.

### Bedert — *Riesz products and the Lonely Runner Conjecture: A wider gap of loneliness*

- **Primary / freshness:** [arXiv:2511.16636v1](https://arxiv.org/abs/2511.16636),
  submitted 2025-11-20. **PREPRINT v1.**
- **Imported role:** proves the asymptotic lower bound
  `1/(2n)+1/n^(5/3+o(1))` by Riesz-product methods.  The repo mines its
  Fourier/resonance mechanisms and uses it as the modern asymptotic comparison.
- **Repo consumers:**
  [THM-515 theta form](../../01-canon/theorems/THM-515-lrc-singular-series-is-the-lonely-measure-theta-form.md),
  [THM-518 route diagnosis](../../01-canon/theorems/THM-518-stranger-decoupling-and-the-bedert-two-route-diagnosis.md),
  [LRC14 proof map](../../00-navigation/LRC14-PROOF-MAP.md).
- **Does not prove:** the conjectural `1/n` threshold, LRC(14), or any lower
  bound for the **measure** of the lonely-time set.  A bound on maximum
  loneliness must not be silently converted into a measure bound.

### Goddyn--Wong — *Tight instances of the lonely runner*

- **Primary / freshness:** [author-hosted paper](https://www.sfu.ca/~goddyn/Papers/063tight_lonely_runner.pdf),
  *Integers* **6** (2006), A38. **PUBLISHED / stable.**
- **Imported role:** constructs an infinite family of tight vectors and
  completely characterizes the one-coordinate accelerations of the canonical
  arithmetic progression.  In repo notation this includes
  `{1,...,11,13,24}` and the `n == 1 (mod 6)` doubling family.
- **Repo consumers:**
  [THM-1065, exact doubling law](../../01-canon/theorems/THM-1065-doubling-family-mod-six-characterization.md),
  [THM-709, acceleration criterion](../../01-canon/theorems/THM-709-doubling-tight-locus-singleton.md),
  [THM-523](../../01-canon/theorems/THM-523-lrc14-singular-series-proof-skeleton.md).
- **Does not prove:** classification of every tight vector.  In particular it
  does not establish uniform emptiness of the twelve-speed sporadic branch or
  show that AP/Goddyn--Wong are the complete LRC(14) tight locus.

### Lee — *Lonely runners in real life: Sharp bounds for time-dependent velocities*

- **RADAR / PREPRINT v1:** [arXiv:2607.16082v1](https://arxiv.org/abs/2607.16082) proves the sharp `2^(-n+1)` extreme-runner guarantee for ordered time-dependent velocities with divergent relative displacement.  No theorem depends on it, and it proves nothing about constant-velocity LRC(14), `1/14`, integer spectra, or finite checking.

## Continued fractions and Khinchin content
- [Khinchin/DS pins](continued-fractions-khinchin-duffin-schaeffer-pins.md) and [Weil pin](stern-depth-kloosterman-weil-pin.md): no named-constant, composite, or LRC claim.

## Jacobian, Dixmier, and Poisson

- **Florit--Smith / Milne:** [Florit--Smith, §4.14](https://arxiv.org/abs/2101.00917) corroborates `Jac(y^2=x(x^4+1))~E_8000^2`; [Milne, CM 1.9/1.10/3.12--3.13](https://www.jmilne.org/math/CourseNotes/CM.pdf) pins primitive-type classification for 4012's wt `9--11` audit.
- **Dokchitser / guardrail:** [*Models of curves over DVRs*, arXiv:1807.00025v2](https://arxiv.org/abs/1807.00025), Definitions 3.7/3.9/3.12 and Theorem 3.14, supplies the general face/edge model used by [THM-4045](../../01-canon/theorems/THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go.md) before its seven exact nodes are resolved. It proves none of the in-repo face arithmetic, Keller map, or `JC(2)` claims; Florit--Smith/Milne likewise do not prove face-stability.

### Degtyarev — weight-eight/nine torus sextics and the stereographic trigonal model

- **Primary / freshness:** A. Degtyarev,
  [*Irreducible plane sextics with large fundamental groups*](https://arxiv.org/abs/0712.2290),
  especially Sections 1.1 and 3.2--3.6, and
  [*Oka's conjecture on irreducible plane sextics. II*](https://arxiv.org/abs/math/0702546).
  **PUBLISHED-era stable sources; exact cited sections pinned 2026-08-24.**
- **Imported role:** over the complex plane, the first source classifies the
  seven simple weight-eight singularity configurations and their four torus
  structures, as well as the simple weight-nine `9A_2` row and its twelve
  torus structures;
  the stereographic construction supplies the two non-simple
  `J_{2,0}+4A_2` and `J_{2,3}+3A_2` families from the four-cuspidal trigonal
  curve on `Sigma_2`.
- **Repo consumers:**
  [THM-3943](../../01-canon/theorems/THM-3943-rational-weight-eight-four-torus-sextics-have-no-one-place-line.md)
  and
  [THM-3945](../../01-canon/theorems/THM-3945-nonsimple-weight-eight-j-sextics-have-no-one-place-line.md).
  Their complete genus-and-line synthesis is
  [THM-3948](../../01-canon/theorems/THM-3948-classified-weight-eight-nine-sextics-have-no-polynomial-normalization.md).
- **Does not prove:** a characteristic-zero base-change theorem making that
  complex classification exhaustive, the repo's binary-sextic line-norm contradictions, the
  uniform no-one-place-line statements, a classification of non-line affine
  boundaries, a Keller obstruction for every torus sextic, or `JC(2)`.  Those
  exact deductions are separated from the cited classification in
  THM-3943/3945/3948 and their companions or cited genus ledger.

### Miyanishi -- affine pseudo-planes and the canonical-trivial exception

- **Primary/imported:** Miyanishi [arXiv:1504.07179](https://arxiv.org/abs/1504.07179), Lemma 2.5.2, constructs the affine-plane Galois pseudo-cover; Theorem 2.5.10 excludes type `(d,n,r)` with `r!=2` from mapping etale to `A2` via `Pic(X)=Z/d`, `K_X~(r-2)F_0`.
- **Modern character notation:** Dubouloz--Palka, [arXiv:1701.01425v2](https://arxiv.org/abs/1701.01425), equations (1.1), (2.5), and (5.1)--(5.2), distinguish `S(k,r,a)` by the deck character `a` and construct nonproper etale self-maps using Belyi--Shabat data.
- **Consumers:** [3785](../../01-canon/theorems/THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable.md) identifies `S(3,3,1)`; [3788](../../01-canon/theorems/THM-3788-dubouloz-palka-standard-chart-containment-obstruction.md) excludes standard plane charts for its `alpha=0` factors. Neither imports old-type `r=2`, proves `JC(2)`, or closes nonstandard charts.

### Gao / Shaska / Meng--Yang — the 2026 Jacobian boundary

- Gao [arXiv:2608.00222v1](https://arxiv.org/abs/2608.00222) gives tangent-sweep counterexamples in dimensions `>2`; [THM-1300](../../01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md)/3438 verify the degree-four map/lift.
- Shaska [arXiv:2607.20210v2](https://arxiv.org/abs/2607.20210) excludes planar equivariant signatures (THM-3686/3691 are analogues); Meng--Yang [arXiv:2607.22198v2](https://arxiv.org/abs/2607.22198) leaves only `JC(2)` and `HC(4)`, with `HC(4) => JC(2)`.

### Nagata / Shastri / Guccione--Guccione--Horruitiner--Valqui — planar degree and infinity gates

- **Primary/imported:** Nagata's corrected [Theorem 7.3](https://repository.kulib.kyoto-u.ac.jp/server/api/core/bitstreams/9ef8e868-5526-4830-b19f-543c0af09e7c/content) gives `Omega(deg R)>=3` for every counterexample-pencil member; Shastri's [criterion](https://doi.org/10.18910/6794) forces two top-form roots; [arXiv:2204.14178v1](https://arxiv.org/abs/2204.14178) gives reduced height `>=108`, only `(72,108)` below `125`.
- **Consumer/boundary:** [THM-3550](../../01-canon/theorems/THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor.md) supplies an independent elementary floor.  None proves `JC(2)`, existence of `(72,108)`, or height `125`.

### Planar Keller target geometry — Lang; Jelonek--Lasoń; Nguyen; Shioda; Miranda--Persson

- **CITED:** [Lang](https://doi.org/10.1016/0022-4049(91)90128-O) proves Newton similarity; [Jelonek--Lasoń](https://arxiv.org/abs/1411.5011), Theorem 3.2/Corollary 3.5, give `D-1`/chart-width bounds; [Nguyen](https://arxiv.org/abs/math/0305088), Theorem 1, gives component degree pair `(M d,M e)`. THM-4122 divides out cover inflation; THM-4124 turns integral similarity into a target shear.
- **CITED:** [Shioda](https://rikkyo.repo.nii.ac.jp/records/10027) supplies the Mordell--Weil framework; [Miranda--Persson](https://www.math.colostate.edu/~miranda/preprints/Miranda-Persson1986_Article_OnExtremalRationalEllipticSurf.pdf), Theorem 4.1/Table 5.2, classify `X_211` as `II*+I1+I1` with section group one. THM-4120 uses this only on one smooth theta seam.
- **Boundary:** raw `M` is not intrinsic; `rho<=2` at `(72,108)` needs width six; base change can enlarge Mordell--Weil; and `Q^2-cP^3` is not a target automorphism. First route any escaping divisor into the nonproperness set.

### Belov--Kanel--Kontsevich — *The Jacobian Conjecture is stably equivalent to the Dixmier Conjecture*

- **Primary/imported:** [arXiv:math/0512171](https://arxiv.org/abs/math/0512171) proves `JC(2n) -> DC(n)` and the stable JC/DC equivalence via finite-characteristic reduction.
- **Consumer/guardrail:** [THM-2044](../../01-canon/theorems/THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension.md). It gives neither `PC(2)->DC(2)`, direct quantization, nor a fixed-rank route from that four-variable Poisson map to `JC(2)` or `DC(2)`.

### Adjamagbo--van den Essen — *On the equivalence of the Jacobian, Dixmier and Poisson Conjectures in any characteristic*

- **Primary/imported:** [arXiv:math/0608009](https://arxiv.org/abs/math/0608009) supplies the characteristic-sensitive stable JC/DC/Poisson framework via Azumaya and reduction methods.
- **Consumer/guardrail:** THM-2044 and its rank-scope reflection. Symbol replacement need not preserve Weyl relations or rank; no implication here settles `DC(2)` or `JC(2)`.

### de Bondt — *Symmetric Jacobians*

- **Primary/imported:** [arXiv:1206.2865v3](https://arxiv.org/abs/1206.2865) gives stable reductions to specified Jacobian symmetries; dimension and the full symmetry class are load-bearing.
- **Consumer/guardrail:** [HYP-8905](../hypotheses/HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs.md) and MISTAKE-237. It does not reduce full `JC(2)` to a binary symmetric problem or give an NC2/GMC implication.

### Zhao — Hessian-nilpotent polynomials and the Vanishing Conjecture

- **Primary / freshness:** [*Some Properties of and Open Problems on Hessian
  Nilpotent Polynomials*, arXiv:0704.1689v2](https://arxiv.org/abs/0704.1689),
  published in *Annales Polonici Mathematici* **93** (2008), and [*A Vanishing
  Conjecture on Differential Operators with Constant Coefficients*,
  arXiv:0704.1691v2](https://arxiv.org/abs/0704.1691), published in *Acta
  Mathematica Vietnamica* **32** (2007).
- **Imported role:** develops Hessian-nilpotent polynomials, their inversion
  pairs, and the Laplacian/constant-coefficient Vanishing Conjecture framework.
  The relevant JC reduction uses homogeneous Hessian-nilpotent targets and the
  eventual vanishing of `Delta^m(P^(m+1))`; `Delta^m(P^m)=0` encodes Hessian
  nilpotence in the cited setting.
- **Repo consumers:** [HYP-8905](../hypotheses/HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs.md),
  [THM-2063](../../01-canon/theorems/THM-2063-one-fiber-linear-planar-keller-pairs.md), [THM-2801](../../01-canon/theorems/THM-2801-sharp-special-image-boundary-and-beta-shift-witness.md), and MISTAKE-237.
- **Does not prove:** equality with complex Gaussian moment functionals, an
  `NC2 -> GMC(2) -> JC(2)` chain, or full JC(2) from the binary homogeneous
  calculation. VC(4), planar Jelonek/leading-form descent, and Newton/Lame
  descent remain separate programs.
- **In-repo GVC(3) witness (2026-08-03, provenance UNKNOWN):** a supplied
  homogeneous `Lambda = Delta^6` witness (`Delta = 4 d_x d_y + d_t^2`,
  `P = A C^2` deg 12, `Q = x^2`) is FINITE-EXACT through `m = 4` with exact
  nonvanishing constants and tight exponent; absent from indexed literature
  after documented search, plausibly a homogeneous lift of Long's GMC(3)
  example via `E[f] = (e^{Delta/2}f)(0)`; dimension-minimal if the all-`m`
  induction holds (de Bondt proved homogeneous-operator GVC in 2 variables).
  See `05-knowledge/results/gvc3-delta6-counterexample-verification-boxeph.md`
  and `gvc-tv-provenance-hunt-boxeph.md`.

### Lee--Li — *On the two-dimensional Jacobian conjecture: Magnus' formula revisited, IV*

- **Primary / freshness:** [arXiv:2408.01279v1](https://arxiv.org/abs/2408.01279), submitted 2024-08-02. **PREPRINT.**
- **Imported role:** develops inner polynomials for a planar Jacobian pair and
  constrains the northeastern Newton-polygon vertex. The repo uses this as the
  current positive framework in which to test THM-2045's weighted-sector /
  exponent-semigroup obstruction on exposed inner edges.
- **Repo consumers:**
  [THM-2045](../../01-canon/theorems/THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate.md),
  [corrected JC/LRC reflection](../../07-reflections/jacobian-and-lonely-runner-two-nullcones-that-diverge-boxeph-S205.md).
- **Does not prove:** planar JC, an arithmetic-progression reduction, or that
  every inner edge has THM-2045's unique constant-producing sector.

### Cheng--McKay--Wang / Moskowicz — the planar Jacobian centralizer theorem

- **Primary:** Cheng--McKay--Wang, [*Younger mates and the Jacobian conjecture*](https://doi.org/10.1090/S0002-9939-1995-1257100-4), *Proc. AMS* **123** (1995), Theorem 1 over `C`; characteristic-zero field form in Moskowicz, [*The two-dimensional Centralizer Conjecture*](https://arxiv.org/abs/1802.04685v2), Theorem 2.3.
- **Imported role:** `Jac(P,B) in k*` and `Jac(P,w)=0` imply `w in k[P]`; [THM-2230](../../01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md) derives that every response fiber is exactly one target-shear orbit.
- **Does not prove:** planar JC or that `P` is a coordinate. The unit-mate and characteristic-zero hypotheses are essential.

### Serre--Tate — *Good reduction of abelian varieties*

- **Primary / freshness:** [Annals of Mathematics **88** (1968), 492--517](https://annals.math.princeton.edu/1968/88-3/p05), Section 1. **PUBLISHED / stable; checked 2026-08-24.**
- **Imported role:** Neron--Ogg--Shafarevich and the resulting invariance of
  potential good reduction under isogeny, used only in THM-3997's independent
  generic-fibre exclusion of the zero-residual reduced `2:3` cell.
- **Does not prove:** the two elliptic models or their `j`-invariants, the
  Keller-induced function-field map, any nonzero-residual exclusion, `JC(2)`,
  or the full Jacobian Conjecture; those remain in-repo arguments or open.

### Han--Pan--Chen — *Normal forms of elements in the Weyl algebra and Dixmier Conjecture*

- **Primary / freshness:** [arXiv:2407.11291v1](https://arxiv.org/abs/2407.11291),
  submitted 2024-07-16. **PREPRINT.**
- **Imported role:** organizes nilpotent or semisimple `A_1` elements by Joseph
  normal-form pairs `(k,n)`, reformulates the remaining Dixmier obstruction in
  those coordinates, and proves the needed nonexistence when `gcd(k,n)=1` or
  one entry is prime; it records analogous planar-Poisson results.
- **Repo consumer:** the `DC(1)` control discussion in the
  [rank-two boundary reflection](../../07-reflections/rank-two-poisson-suspension-leaves-two-rigidity-gates-dc2-and-planar-jc-codex-20260721.md).
- **Does not prove:** DC(1), DC(2), planar JC, or termination of THM-2044's
  observed six-weight quantization cascade.

## Gaussian Moments Conjecture / NC2

### [BGKP](https://arxiv.org/abs/2310.18020)
- Thms. 1.4/2.11: `m<=n` gives Schur-ratio monotonicity including zeros; THM-3110 slides only, not its signed bank.

### Derksen--van den Essen--Zhao — *The Gaussian Moments Conjecture and the Jacobian Conjecture*

- **Primary / freshness:** [arXiv:1506.05192](https://arxiv.org/abs/1506.05192),
  *Israel Journal of Mathematics* **219** (2017), 917--928.
- **Imported role:** distinguishes operator `E_n` from scalar `F_n` and proves
  `GMC(2n) => ker(F_n) MZ => SIC(n)`, plus the global GMC-to-JC implication.
- **Repo consumers:**
  [THM-2022, GMC(2)](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md),
  [THM-1490, earlier explicit GMC counterexample](../../01-canon/theorems/THM-1490-the-gaussian-moment-counterexample-verified-proved-shortened-and-obstructed.md),
  [THM-2801, sharp SIC boundary](../../01-canon/theorems/THM-2801-sharp-special-image-boundary-and-beta-shift-witness.md),
  [GMC2 finish map](../../00-navigation/GMC2-FINISH-MAP-2026-07-21.md).
- **Does not prove:** any converse: THM-2022 proves `GMC(2)` while THM-2801
  refutes `SIC(2)`.  For the true endpoint, van den Essen--Wright--Zhao
  [arXiv:1008.3962v2](https://arxiv.org/abs/1008.3962), Theorem 2.8, gives
  `SIC(1)` directly from the prime ideal `(xi)`; neither route proves `JC(2)`.

### Duistermaat--van der Kallen — *Constant terms in powers of a Laurent polynomial*

- **Primary / freshness:** [author PDF](https://webspace.science.uu.nl/~kalle101/powers.pdf),
  *Indagationes Mathematicae* **9** (1998), 221--231,
  [DOI 10.1016/S0019-3577(98)80020-7](https://doi.org/10.1016/S0019-3577(98)80020-7).
- **Imported role:** classifies Laurent polynomials whose every positive power
  has zero constant term.  Its one-variable theorem is a stronger historical
  route to the nonzero constant-term seed; the current internal dependency is
  the elementary Galois-orbit proof in THM-2067.
- **Repo consumers:**
  [THM-1630, exact citation identification](../../01-canon/theorems/THM-1630-tnc-is-duistermaat-van-der-kallen-theorem-2.md),
  [THM-1645, angular/radial split](../../01-canon/theorems/THM-1645-gmc2-angular-layer-is-dvdk-the-gap-is-purely-radial.md),
  [THM-2067, internal bare-existence replacement](../../01-canon/theorems/THM-2067-galois-orbit-product-closes-one-variable-dvdk.md),
  [THM-2022](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md),
  and [THM-2070, failed-bypass correction](../../01-canon/theorems/THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation.md).
- **Does not prove:** the Gaussian radial/factorial step, Frobenius survival of
  a complete Wick face, or an effective first return. Lean now root-imports
  the height supplier and positive/two-charge/fixed-support unique-channel
  leaves; general complex `DvdK1` remains the sole formal premise.
  HYP-8931's proposed class bypass is vacuous (MISTAKE-240). THM-2070 refutes
  the general semigroup/saddle bypass: support return is not noncancellation.

### Edo--van den Essen — *The Strong Factorial Conjecture*

- **Primary / freshness:** [arXiv:1304.3956v2](https://arxiv.org/abs/1304.3956);
  *J. Algebra* **397** (2014), 443--456, [DOI](https://doi.org/10.1016/j.jalgebra.2013.09.011).
  **PUBLISHED / stable; source and journal record checked 2026-07-24.**
- **Imported role / notation gate:** Definitions 2.1/2.7 and Conjectures 2.4/2.8 index `FC(m)`/`SFC(m)` by ambient variables; `N(f)` is support/window size. Thus `t` univariate terms mean `SlotSFC_1(t)`, not `SFC(t)` (MISTAKE-350).
  The Gaussian identity
  `E[(W+Zh)^(2j)]=binom(2j,j)L((sh)^j)` identifies its `n=1` two-charge case.
- **Repo consumers:** [THM-1790](../../01-canon/theorems/THM-1790-the-emp-floor-detection-depth-at-least-degree-plus-one.md)
  and [HYP-8765](../hypotheses/HYP-8765-gmc2-radial-channel-return-tower.md).
- **Does not prove:** SFC, the full HYP-8765 cutoff, or a uniform NC2 depth.
  THM-2022 instead uses a coefficient-dependent good prime; THM-1790 proves
  the complementary projective lower bound.

### Erman--Smith--Várilly-Alvarado — *Laurent polynomials and Eulerian numbers*

- **Primary / freshness:** [arXiv:0908.2609v2](https://arxiv.org/abs/0908.2609),
  *Journal of Combinatorial Theory A* **118** (2011), 396--402,
  [DOI 10.1016/j.jcta.2010.02.006](https://doi.org/10.1016/j.jcta.2010.02.006).
- **Imported role:** studies generic constant-term equations through toric
  geometry and Eulerian numbers and records the sharp Sturmfels effective
  question.  In the repo this is the target “first nonzero constant term by
  exponent-width `m+n`,” with Newton-polygon branch count as its geometry.
- **Repo consumers:**
  [THM-1630 effective redirect](../../01-canon/theorems/THM-1630-tnc-is-duistermaat-van-der-kallen-theorem-2.md),
  [THM-1650 Newton polygon](../../01-canon/theorems/THM-1650-newton-polygon-of-the-effective-dvdk-bound.md).
- **Does not prove:** the Sturmfels `m+n` effective bound, NC2, or GMC(2).
  Generic regular-sequence/degree information is not a uniform first-return
  theorem for every Laurent polynomial.

### Baricz--Singh — *Zeros of some special entire functions*

- **Primary / freshness:** [arXiv:1702.00626v2](https://arxiv.org/abs/1702.00626),
  *Proceedings of the AMS* **146** (2018), 2207--2216,
  [DOI 10.1090/proc/13927](https://doi.org/10.1090/proc/13927).
- **Imported role:** the positive-parameter hyper-Bessel / `0F_q` zero theorem
  places all zeros on the negative real axis.  Gauss multiplication makes the
  repo's `Phi_(p,q)(x)=sum x^k/((pk)!(qk)!)` an exact instance.
- **Repo consumer:**
  [THM-2023, Laguerre--Polya boundary theorem](../../01-canon/theorems/THM-2023-gmc2-hyperbessel-boundary-laguerre-polya.md).
- **Does not prove:** GMC or NC2.  It handles the `Phi` sharp boundary, not the
  opposite Mittag--Leffler/root-of-unity-filter `Psi` boundary, and it does not
  exclude cancellation at an allowed negative-real zero without further work.

### Mangoubi--Kadets--Weller Weiser — common roots of Legendre polynomials

- **Primary / freshness:** [IAS seminar announcement, 2026-02-18](https://www.ias.edu/math/events/special-year-learning-seminar-41).
  **SEMINAR ONLY:** no public paper or proof was located on 2026-07-21.
- **Imported role:** the announcement states that only finitely many Legendre
  polynomials vanish at a fixed nonzero point, with quantitative bounds to be
  discussed, and attributes the work to Dan Mangoubi, Borys Kadets, and Adi
  Weller Weiser.
- **Repo consumer:**
  [THM-2021, Legendre refinement](../../01-canon/theorems/THM-2021-gmc2-legendre-finite-recurrence-closure.md).
- **Does not prove:** a citable published finite-recurrence
  theorem until a paper/proof is public.  It is no longer a dependency of
  NC2: [THM-2018's recurrence mechanism](../../01-canon/theorems/THM-2018-gmc2-resonance-algebraic-egf-and-proportional-shadow.md)
  and [THM-2022's full proof](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md)
  supersede that need.  It also does not establish Stieltjes' stronger “at
  most one” claim.

### Long — *Small Counterexamples to the Gaussian Moments Conjecture*

- **Primary / freshness:** [arXiv:2607.18186v1](https://arxiv.org/abs/2607.18186),
  submitted 2026-07-20. **RADAR / PREPRINT v1.**
- **Imported role:** gives explicit `P,Q` in **three** independent standard real
  Gaussian variables with `E(P^m)=0` and `E(QP^m)=m!` for every `m>=1`; `P`
  has five terms and degree four.  A six-term cubic example in four variables
  is also supplied.  Hence GMC is false in every dimension at least three.
- **Repo consumers:** pair it with
  [THM-2022](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md)
  for the exact dimensional boundary, and with
  [THM-1490](../../01-canon/theorems/THM-1490-the-gaussian-moment-counterexample-verified-proved-shortened-and-obstructed.md),
  whose older outside example only reached four real variables and made no
  priority claim.
- **Does not prove:** failure of GMC(2).  The paper explicitly says its small
  Gaussian examples were not derived from the announced Jacobian map; do not
  conflate the direct counterexamples with the high-dimensional
  Jacobian-to-GMC route.

## Knots and compositional relations

Detailed source/import/guardrail records for Brittenham--Hermiller, Zakharov, Rybin--Zhang--Luo, Kominers--Mrazovic--Pomerance--Sole, Krenn--Gu--Soltesz, Schubert, and Owens--Strle are preserved in the
[compositional-relations source sidecar](CORE-PAPERS-COMPOSITIONAL-RELATIONS.md). THM-2176 derives the continuation profile, positive interaction coboundary, and pure-bypass split from cited crossing certificates.
THM-2191 uses prime-decomposition cancellation and the four-ball crossing metric to construct catalytic Gordian localization and close the `9_10` reverse direction.

## Tournaments ([Moon 1966 sidecar](CORE-PAPERS-TOURNAMENT-STRONG-EARS-2026-08-25.md))

**Schmerl--Trotter — PUBLISHED / stable; checked 2026-08-26.** [*Discrete Mathematics* **113** (1993), 191--205](https://doi.org/10.1016/0012-365X(93)90516-V) classifies the three critical prime tournament families at odd order `>=5`, and none at even order.
THM-4169 imports only prime deletion; all other claims are in-repo.

### Grinberg--Stanley — *The Rédei--Berge symmetric function of a directed graph*

- **Primary / freshness:** [arXiv:2307.05569v1](https://arxiv.org/abs/2307.05569),
  submitted 2023-07-10. **DRAFT:** the authors describe portions of Sections
  3--8 as an outline.
- **Imported role:** develops the Rédei--Berge symmetric function, proves
  power-sum positivity for tournaments, recovers Rédei/Berge Hamiltonian-path
  theorems, and yields the odd-cycle-collection formula
  `H(T)=I(Omega(T),2)` plus its mod-4 refinement.
- **Repo consumers:**
  [THM-002 OCF](../../01-canon/theorems/THM-002-ocf.md),
  [THM-466 higher Rédei tower](../../01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md),
  [CONJ-001/Claim A](../../01-canon/theorems/CONJ-001-claim-a.md).
- **Does not prove:** that OCF was first discovered in this repository, that
  Hamiltonian-path count classifies tournaments, or any LRC/Paley equivalence.
  Cite this paper as Grinberg--Stanley; **do not attach arXiv:2412.10572 to
  their names.**

### Irving--Omar — *Revisiting the Rédei--Berge Symmetric Functions via Matrix Algebra*

- **Primary / freshness:** [arXiv:2412.10572v3](https://arxiv.org/abs/2412.10572),
  revised 2025-10-23. **PREPRINT / edited manuscript.**
- **Imported role:** gives a matrix-algebra consolidation and extension of
  Chow, Grinberg--Stanley, and Lass, including Hamiltonian-path consequences.
  Its Corollary 20 is the repo's cleanest explicit citation of OCF, attributed
  there to Grinberg--Stanley.
- **Repo consumers:**
  [THM-002 OCF](../../01-canon/theorems/THM-002-ocf.md),
  [THM-466](../../01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md),
  [canonical definitions](../../01-canon/definitions.md).
- **Does not prove:** priority for the underlying Grinberg--Stanley formula or
  a complete tournament invariant.  **Correct attribution:** this arXiv ID is
  Irving--Omar, not Grinberg--Stanley.

### Klanderman--Montee--Piotrowski--Rice--Shader — *Determinants of Seidel Tournament Matrices*

- **Primary / freshness:** [arXiv:2406.09697v1](https://arxiv.org/abs/2406.09697),
  published in *Linear Algebra and its Applications* **707** (2025), 126--151,
  [DOI 10.1016/j.laa.2024.11.011](https://doi.org/10.1016/j.laa.2024.11.011).
- **Imported role:** studies the attainable square roots of Seidel
  determinants, their gaps and bounds, and exact random-tournament moments.
  Theorem 7.7 identifies the expected characteristic polynomial with the
  signless matching polynomial; the paper also supplies the known count of
  distinct Seidel characteristic polynomials used by the sequence atlas.
- **Repo consumers:**
  [THM-473 involution/Hermite corollary](../../01-canon/theorems/THM-473-involution-average-hermite.md),
  [THM-2010 sequence atlas](../../01-canon/theorems/THM-2010-new-tournament-invariant-sequences.md).
- **Does not prove:** the repo's `E det(I+S)=involutions` interpretation,
  adjacency-spectrum counts, or novelty of any four-term candidate sequence.
  THM-2010 now uses these five primary-record authors; older historical logs
  may retain the superseded attribution.

## Arrangements and tournament games

Detailed source/import/guardrail records for De Concini--Procesi, Moci,
Stanley (arrangement surveys), and Fisher--Ryan are preserved verbatim in the
[arrangements/games sidecar](CORE-PAPERS-ARRANGEMENTS-GAMES.md). All four are
consumed by corrected reflections (S209/S210, MISTAKE-223/224 lineage): toric
arrangements and Shi/braid counts are audit controls, not LRC carriers, and a
tournament-game pure optimum means a Condorcet winner, not transitivity.

## Analytic number theory

Williams' fixed-modulus Mertens theorem, THM-3793's exact import, and its
no-support-asymptotic guardrail are recorded in the [analytic-number-theory sidecar](CORE-PAPERS-ANALYTIC-NUMBER-THEORY.md).

## Reciprocal integer sequences

The [reciprocal-sequences sidecar](CORE-PAPERS-RECIPROCAL-SEQUENCES.md) records
the figurate/toothpick sources used by THM-2000/2005; neither proves the
support--multiplicity collision law or Abel--Stieltjes/Dini criterion.

## Pythagorean trees and square-pyramidal intersections

The [Pythagorean sidecar](CORE-PAPERS-PYTHAGOREAN.md) records Berggren,
fixed-hypotenuse, and Bennett imports used by THM-3334/3335; none supplies an
LRC, tournament, Jacobian, or skew-EW construction.

## Unstable homotopy

Ivanov--Mikhailov--Wu ([arXiv:1506.00952](https://arxiv.org/abs/1506.00952),
*HHA* **18** (2016) 337--344; `pi_n(S^2)` nontrivial for `n>=2`) is recorded in
the [homotopy sidecar](CORE-PAPERS-HOMOTOPY.md); repo consumers THM-3204/3205
(continuant gate, odd-primary lambda-algebra engine). It does not prove any
new homotopy group, `p=2` statement, or LRC/GMC consequence.

## Probability and extraction

### Kontorovich — *TV Homogenization Inequalities*

- **Primary / freshness:** [arXiv:2601.04079v3](https://arxiv.org/abs/2601.04079),
  submitted 2026-01, revised 2026-02. **PREPRINT.**
- **Imported role:** homogenization (each Bernoulli parameter mapped to the
  block mean) reduces TV distance up to a universal constant. Lemma 1.4
  PROVES `delta_N <= 2(delta_I + delta_J)` for a block partition; the
  product form `delta_N <= delta_I + delta_J - delta_I delta_J` is the
  paper's stated CONJECTURE ("no pathway" via its methods). The repo PROVES
  that conjecture's first case (`|I|=|J|=1`, i.e. `n=2`) by an exact
  positivity certificate, with the equality face characterized — see
  `05-knowledge/results/tv-fusion-homogenization-lemma-boxeph.md`.
- **Does not prove:** the product-form conjecture beyond `n=2` (OPEN); any
  AMM 12592 deadline bound (transfer refuted in the repo note); "rigid face"
  is repo/user terminology, not the paper's.
- **Repo consumer:** [THM-3291](../../01-canon/theorems/THM-3291-two-block-tv-homogenization-rigidity.md)
  proves the `n=2` case from a box constraint plus AM-GM and classifies the
  equality locus; its 62 nontrivial equality points are boxeph's 351 minus the
  degenerate faces, an exact two-implementation cross-confirmation.

### Zhao's GVC in three variables — the supplied object, and its dictionary

- **Provenance UNRESOLVED — do not cite.**  The three-variable object
  `rho=t^2+xy`, `A=rho+x^2`, `C=(rho^3-t^2A^2)/x`, `P=AC^2`, `Q=x^2`,
  `Delta=4d_x d_y+d_t^2` was supplied with the identifier `arXiv:2606.17854`.
  That identifier resolves to Ajwani--Gajjala--Raman--Ray, *Counterexamples to
  Wegner's Conjecture for Rectangles* (cs.CG), which contains none of it.  The
  mathematics is correct; the source is unknown.  Attach no priority claim.
- **Dictionary (reusable).**  For a nondegenerate quadratic form `rho` with
  Laplacian `Delta` and `L(f)=(exp(Delta/2)f)(0)`, a degree-`2k` form has
  `L(f)=Delta^k f/(2^k k!)`.  So a *Generalized Vanishing Conjecture* statement
  about `Delta^j` in `n` variables **is** a *Gaussian Moments Conjecture*
  statement in `n` variables; they are not separate lanes.  Zhao's equivalence
  with the Jacobian Conjecture is for `j=1` only, so no statement about powers
  of `Delta` touches `JC`, `JC(2)`, or THM-1435's VC-witness bracket.
- **Repo consumers:** boxeph's FINITE-EXACT instance verification
  (`05-knowledge/results/gvc3-delta6-counterexample-verification-boxeph.md`)
  and [THM-3290](../../01-canon/theorems/THM-3290-archimedes-flatness-and-the-gmc3-gvc3-counterexample-family.md),
  which proves both all-`m` statements that verification left open and
  generalizes them to an infinite family.

## August 2026 preprint intake
- [August 14 intake](CORE-PAPERS-INTAKE-2026-08-14.md): preprints/Rule 30.
- August 25--27: [p-adic/matching](CORE-PAPERS-INTAKE-2026-08-25-PADIC-ZETA-MATCHING-LOGIC.md), [superseded specialization](CORE-PAPERS-INTAKE-2026-08-26-PADIC-ZETA-DENSITY.md), [density/percolation](CORE-PAPERS-INTAKE-2026-08-26-PADIC-DENSITY-PERCOLATION.md). Long is **PREPRINT / UNDER AUDIT**; no `u-f` map. Chen--Rosu is unrefereed; Cerf unrelated.

## Maintenance rule

For updates, refresh the primary link/version and exact imported theorem; name
the first consumer; sharpen **does not prove** and attribution before routing.

## External intake

[2026-07-28 owner puzzle bundle — **CITED-ABSTRACT only**](CORE-PAPERS-INTAKE-2026-07-28.md).
