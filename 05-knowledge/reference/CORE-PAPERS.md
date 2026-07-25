# Core papers: imported results, consumers, and guardrails

> **Freshness:** primary records checked **2026-07-21**.  An arXiv version in
> this file is the version visible on that date; recheck entries marked
> **PREPRINT**, **RADAR**, or **SEMINAR ONLY** before making a priority or
> state-of-the-art claim.

This is a role map, not a general bibliography.  Read an entry to learn exactly
what this repository imports from a source, where that input is consumed, and
what the source does **not** establish.  Repository theorem numbers collide in
places, so the links below, not a bare `THM-N`, are the canonical addresses.

## Fast frontier snapshot

- **LRC:** fourteen runners (thirteen nonzero speeds) remains open; the April 2026 computation reaches twelve nonzero speeds.
- **Gaussian moments:** [THM-2022](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md) proves the two-real-Gaussian case; Long's July 2026 preprint gives three-dimensional counterexamples, so the checked boundary is `true through 2 / false from 3`.
- **Tournament attribution:** arXiv:2412.10572 is Irving--Omar, arXiv:2307.05569 Grinberg--Stanley, and arXiv:2406.09697 Klanderman--Montee--Piotrowski--Rice--Shader.
- **Unsplittable flows:** [THM-2177](../../01-canon/theorems/THM-2177-planar-counterexample-to-goemans-unsplittable-cost-flow-conjecture.md) refutes Goemans' exact `d_max` cost conjecture; the cost-free `d_max` and planar cost-preserving `2d_max` theorems remain intact.
- **Reciprocal sequences:** separate support from indexed multiplicity; external figurate formulas do not supply the repo's Abel--Dini/support-Dirichlet extensions.
- **Jacobian/Dixmier/Poisson:** THM-1300 refutes `JC(n)` for `n>=3`; `JC(2)`, `DC(1)`, and `DC(2)` remain open. THM-2071 closes quadratic-fiber pencil cells only; THM-2044 separately refutes the two-pair Poisson conjecture.

## Lonely Runner Conjecture

### Sungkawichai--Trakulthongchai — *Eleven, twelve, and thirteen lonely runners*

- **Primary / freshness:** [arXiv:2604.23906v1](https://arxiv.org/abs/2604.23906),
  submitted 2026-04-26. **PREPRINT v1; not peer reviewed as of this check.**
- **Imported role:** a computer-assisted divisibility sieve proves LRC for
  `k=10,11,12` nonzero speeds.  The paper also gives a finite-field polynomial
  shortcut when `k+1` and a sufficiently large auxiliary modulus are odd
  primes, and identifies the initial improper-tuple computation as the
  extension bottleneck.
- **Repo consumers:**
  [LRC14 frontier](../../00-navigation/LRC14-FRONTIER-2026-07-15.md),
  [finite-check feasibility ledger](../../00-navigation/LRC14-FINITE-CHECK-FEASIBILITY-LEDGER-2026-07-19.md),
  [THM-1288, literal Conjecture 7.1 refutation](../../01-canon/theorems/THM-1288-c71-refuted-divisor-aligned-clusters.md),
  [THM-523](../../01-canon/theorems/THM-523-lrc14-singular-series-proof-skeleton.md).
- **Does not prove:** `k=13`, i.e. the standard fourteen-runner case.  Its
  Conjecture 7.1 is false as literally written by the exact family in THM-1288;
  do not use it as an inverse theorem.  The prime shortcut also does not cross
  the composite `k+1=14` lift wall.

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

- **Primary / freshness:** [arXiv:2607.16082v1](https://arxiv.org/abs/2607.16082),
  submitted 2026-07-17. **RADAR / PREPRINT v1.**
- **Imported role:** for locally integrable, almost-everywhere strictly ordered
  time-dependent velocities with divergent pairwise relative displacement,
  proves the sharp `2^(-n+1)` guarantee for each extreme runner and constructs
  examples showing intermediate runners may remain arbitrarily close to
  someone forever.  This is a useful stress test for which LRC arguments use
  linear flow rather than order alone.
- **Repo landing point:** no theorem currently depends on it; compare the
  constant-speed assumptions in the
  [LRC14 frontier](../../00-navigation/LRC14-FRONTIER-2026-07-15.md) before
  importing any mechanism.
- **Does not prove:** anything about the standard constant-velocity LRC(14),
  its `1/14` threshold, integer spectra, or finite checking.

## Jacobian, Dixmier, and Poisson

### July 2026 explicit dimension-three Jacobian certificate — public announcement, no paper

- **Primary / freshness:** [original public announcement](https://x.com/__alpoge__/status/2079028340955197566),
  2026-07-19, and [exact displayed map and certificate](https://jacobianfun.org/jacobian-explained).
  **ANNOUNCEMENT ONLY:** no stable paper or peer review was located by 2026-07-21;
  do not infer discovery credit beyond the primary wording.
- **Imported role:** the displayed map has constant Jacobian determinant `-2`
  and an explicit triple collision, so it refutes `JC(3)` and, by padding,
  `JC(n)` for all `n>=3`.
- **Repo consumer:** [THM-1300](../../01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md)
  independently checks both identities exactly and records the downstream Weyl
  and equivariant anatomy.
- **Does not prove:** `JC(2)`, `DC(1)`, or `DC(2)`.  Verification, provenance,
  publication, and priority are separate claims; use THM-1300's byte-zero
  correction before repeating any historical attribution in its body.

### Belov--Kanel--Kontsevich — *The Jacobian Conjecture is stably equivalent to the Dixmier Conjecture*

- **Primary / freshness:** [arXiv:math/0512171](https://arxiv.org/abs/math/0512171),
  submitted 2005-12-08.
- **Imported role:** proves the finite-direction implication `JC(2n) -> DC(n)`
  in characteristic zero and, with the known implication in the other stable
  direction, identifies the stable Jacobian and Dixmier conjectures. Reduction
  to finite characteristic is the central mechanism.
- **Repo consumers:**
  [THM-2044 rank-two Poisson suspension](../../01-canon/theorems/THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension.md),
  [DC2/JC2 boundary reflection](../../07-reflections/rank-two-poisson-suspension-leaves-two-rigidity-gates-dc2-and-planar-jc-codex-20260721.md).
- **Does not prove:** `PC(2) -> DC(2)`, a direct quantization of a given
  symplectic polynomial map, or a fixed-dimension equivalence that turns
  THM-2044 into a planar Jacobian or `A_2` Dixmier counterexample.

### Adjamagbo--van den Essen — *On the equivalence of the Jacobian, Dixmier and Poisson Conjectures in any characteristic*

- **Primary / freshness:** [arXiv:math/0608009](https://arxiv.org/abs/math/0608009),
  submitted 2006-08-01.
- **Imported role:** supplies the characteristic-sensitive stable equivalence
  framework among canonical polynomial, Weyl, and Poisson endomorphism
  conjectures, using Azumaya and reduction techniques. It is the correct
  background for interpreting a Poisson counterexample's stable consequences.
- **Repo consumers:** THM-2044 and the
  [rank/scope separation reflection](../../07-reflections/rank-two-poisson-suspension-leaves-two-rigidity-gates-dc2-and-planar-jc-codex-20260721.md).
- **Does not prove:** that replacing classical brackets by Weyl commutators in
  the same displayed formulas preserves the relations, that stable transport
  preserves the same rank, or that THM-2044 settles DC(2) or JC(2).

### de Bondt — *Symmetric Jacobians*

- **Primary / freshness:** [arXiv:1206.2865v3](https://arxiv.org/abs/1206.2865),
  revised 2013-10-24; published in *Open Mathematics*.
- **Imported role:** proves reductions of the Jacobian Conjecture to polynomial
  maps with specified Jacobian symmetries. In particular, it is enough in the
  stated stable reduction to handle homogeneous maps whose Jacobian has all
  square symmetries. Dimension and the precise symmetry class are
  load-bearing; a planar input need not remain a binary symmetric problem.
- **Repo consumers:** [HYP-8905 corrected map](../hypotheses/HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs.md)
  and [MISTAKE-237](../../01-canon/MISTAKES.md).
- **Does not prove:** the binary symmetric case is equivalent to full JC(2),
  that a reduction into four variables can be solved in two variables, or any
  NC2/GMC-to-JC implication.

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
  [THM-2063](../../01-canon/theorems/THM-2063-one-fiber-linear-planar-keller-pairs.md),
  and MISTAKE-237.
- **Does not prove:** equality with complex Gaussian moment functionals, an
  `NC2 -> GMC(2) -> JC(2)` chain, or full JC(2) from the binary homogeneous
  calculation. VC(4), planar Jelonek/leading-form descent, and Newton/Lame
  descent remain separate programs.

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

### Derksen--van den Essen--Zhao — *The Gaussian Moments Conjecture and the Jacobian Conjecture*

- **Primary / freshness:** [arXiv:1506.05192](https://arxiv.org/abs/1506.05192),
  *Israel Journal of Mathematics* **219** (2017), 917--928.
- **Imported role:** defines GMC in its standard Mathieu-subspace form and
  proves that the Jacobian Conjecture would follow from GMC; it also records a
  counterexample to a broader Moments Vanishing statement.
- **Repo consumers:**
  [THM-2022, GMC(2)](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md),
  [THM-1490, earlier explicit GMC counterexample](../../01-canon/theorems/THM-1490-the-gaussian-moment-counterexample-verified-proved-shortened-and-obstructed.md),
  [GMC2 finish map](../../00-navigation/GMC2-FINISH-MAP-2026-07-21.md).
- **Does not prove:** GMC in any open dimension.  The conjecture is now false
  from three real Gaussian variables onward; the implication to the Jacobian
  Conjecture is not a converse and supplies no Gaussian witness by itself.

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
- **Imported role:** Definition 2.7/Conjecture 2.8 say a `t`-monomial has a
  nonzero moment in every `t` consecutive powers. The Gaussian identity
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

Detailed source/import/guardrail records for Brittenham--Hermiller, Zakharov, Rybin--Zhang--Luo, Krenn--Gu--Soltesz, Schubert, and Owens--Strle are preserved in the
[compositional-relations source sidecar](CORE-PAPERS-COMPOSITIONAL-RELATIONS.md). THM-2176 derives the continuation profile, positive interaction coboundary, and pure-bypass split from cited crossing certificates.
THM-2191 uses prime-decomposition cancellation and the four-ball crossing metric to construct catalytic Gordian localization and close the `9_10` reverse direction.

## Tournaments

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

### De Concini--Procesi — *On the geometry of toric arrangements*

- **Primary / freshness:** [arXiv:math/0505351v4](https://arxiv.org/abs/math/0505351),
  revised 2005-06-08.
- **Imported role:** defines a toric arrangement as a finite family of
  codimension-one subtori or their cosets and computes complement cohomology.
  It supplies the standard object against which HYP-8830's LRC terminology is
  audited.
- **Repo consumer:** the
  [corrected phase-height reflection](../../07-reflections/orlik-solomon-is-a-repo-wide-pattern-toric-arrangements-are-the-lrc-lens-boxeph-S209.md).
- **Does not prove:** that the infinite Fourier annihilator
  `{k:k.v=0}` is a layer poset, that a thickened safe set is an ordinary
  complement, or that complement cohomology retains LRC height, wall side,
  owner, sign, or deletion data.

### Moci — *A Tutte polynomial for toric arrangements*

- **Primary / freshness:** [arXiv:0911.4823v5](https://arxiv.org/abs/0911.4823),
  revised 2010-11-09; final journal version in *Transactions of the AMS*.
- **Imported role:** introduces the multiplicity Tutte polynomial, proves
  deletion--restriction and positivity, and obtains characteristic and
  Poincare polynomials of a toric arrangement as specializations.
- **Repo consumer:** corrected HYP-8830 and MISTAKE-224.
- **Does not prove:** that the cutoff statistic `N_R`,
  [THM-1820's sinc-weighted Fourier series](../../01-canon/theorems/THM-1820-lrc-is-a-moment-nullcone-problem-relation-lattice-pairing.md),
  or an LRCMod count is an arithmetic Tutte or
  Mobius specialization. A finite character list and its actual layers must be
  stated first.

### Stanley — *Hyperplane Arrangements, Interval Orders and Trees*

- **Primary / freshness:** [author-hosted survey PDF](https://math.mit.edu/~rstan/papers/nas.pdf),
  version dated 1995-12-01.
- **Imported role:** records the braid characteristic polynomial and its `n!`
  real chambers, and the Shi formulas
  `chi(q)=q(q-n)^(n-1)` and `r=(n+1)^(n-1)`. These are the exact classical
  controls replayed in the corrected S209 computation.
- **Repo consumer:** corrected HYP-8830, MISTAKE-224, and the
  [S209 arrangement reflection](../../07-reflections/orlik-solomon-is-a-repo-wide-pattern-toric-arrangements-are-the-lrc-lens-boxeph-S209.md).
- **Does not prove:** that Shi walls are LRC safety walls, that braid
  cohomology is a per-tournament invariant, or that finite-field arrangement
  counts solve the LRC(14) AP-core extraction problem.

### Fisher--Ryan — *Tournament games and positive tournaments*

- **Primary / freshness:** [DOI 10.1002/jgt.3190190208](https://doi.org/10.1002/jgt.3190190208),
  *Journal of Graph Theory* **19** (1995), 217--236.
- **Imported role:** the tournament game has a unique optimal strategy with
  odd-cardinality support. The repo's elementary parity reading restores the
  essential second coordinate: an even principal tournament payoff block is
  nonsingular modulo two.
- **Repo consumer:** the
  [corrected antisymmetry reflection](../../07-reflections/antisymmetry-is-the-hinge-tori-odd-functions-saddles-and-tournaments-boxeph-S210.md).
- **Does not prove:** that pure optimum is equivalent to transitivity, that
  every intransitive tournament has recurrent replicator dynamics, or that a
  game saddle is a Morse saddle. Pure optimum means a Condorcet winner.

## Reciprocal integer sequences

### Downey--Ong--Sellers — *Beyond the Basel Problem: Sums of Reciprocals of Figurate Numbers*

- **Primary / freshness:** [author-hosted preprint](https://www.d.umn.edu/~jsellers/downey_ong_sellers_cmj_preprint.pdf),
  [archived copy](https://web.archive.org/web/20130529032918/http://www.math.psu.edu/sellersj/downey_ong_sellers_cmj_preprint.pdf),
  *College Mathematics Journal* **39** (2008), 391--394,
  [JSTOR 27646686](https://www.jstor.org/stable/27646686).
- **Imported role:** supplies the classical telescoping program and exact
  reciprocal sums for figurate-number families, including the triangular
  identity `sum 1/T_n=2`.  These are the seed rows for the repo's polygonal,
  simplex, digamma, and trigamma mass surfaces.
- **Repo consumers:**
  [THM-2000 support-harmonic/figurate surface](../../01-canon/theorems/THM-2000-support-harmonic-abel-dini-figurate-surface.md),
  [THM-2005 support-Dirichlet atlas](../../01-canon/theorems/THM-2005-support-dirichlet-automatic-tournament-atlas.md).
- **Does not prove:** the repo's support-versus-multiplicity collision law,
  Abel--Stieltjes/Dini iff criterion, iterated Bertrand boundary, full
  support-Dirichlet profile, or tournament reciprocal atlas.

### Applegate--Pol--Sloane — *The Toothpick Sequence and Other Sequences from Cellular Automata*

- **Primary / freshness:** [arXiv:1004.3036v2](https://arxiv.org/abs/1004.3036),
  *Congressus Numerantium* **206** (2010), 157--191.
- **Imported role:** develops the toothpick cellular automaton and its dyadic
  recurrences and product-like generating functions.  The repo uses the exact
  A139250 formula as an automatic-sequence test case and as a model for
  operation-generated tournament-adjacent sequences.
- **Repo consumers:**
  [THM-2000](../../01-canon/theorems/THM-2000-support-harmonic-abel-dini-figurate-surface.md),
  [toothpick bridge tangent](../../00-navigation/TANGENTS.md),
  [Perron/toothpick experiment](../../04-computation/perron_toothpick_klein_S315.py).
- **Does not prove:** convergence or a closed form for the support reciprocal
  mass, a tournament interpretation, or novelty of any repo-generated
  sequence.  Cellular-automaton multiplicities must still be collapsed before
  applying the support-harmonic lens.

## Maintenance rule

When a source changes version or a new paper lands:

1. update the per-entry freshness line and primary link;
2. state the exact theorem imported, not just the topic;
3. add the first canonical consumer that actually uses it;
4. sharpen the **does not prove** line before changing any frontier claim; and
5. record attribution corrections here immediately, then repair stale consumer
   files in a separate coherent edit.
