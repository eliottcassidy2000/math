# Limits as observer reconstruction: six bridges and their missing sidecars

Session272, 2026-11-11. **PROVED** below means an elementary lemma proved in this note, not a new Yang–Mills existence theorem or canon promotion. **CITED/TRUSTED** means the supplied source is accepted without re-verification. **CONDITIONAL/OPEN** distinguishes the proposed applications. No additional continuum probability measure is claimed constructed.

## Source and inheritance

External source, fixed throughout: [OpenAI-YM-20261028, commit 3456e5e75cf8d4cd558955a5c7638fba2d68f464](https://github.com/yang-mills-group/OpenAI-YM-20261028/tree/3456e5e75cf8d4cd558955a5c7638fba2d68f464), especially `PublicWilson.lean:WilsonFamily.Converges` and `ChallengeSpec.lean:JointExistence`, `Hierarchy.Reconstructed`, `Hierarchy.Properties`. These files were read directly. The fixed Wilson family precedes every requested increasing volume sequence; the hierarchy may depend on its chosen refinement. Mesh convergence is full at each suitable fixed volume. The lattice limit is the **unselected** hierarchy, which is a character mixture; the analytic continuation and Wightman representation select a character. Tests in the Wilson assertion are compact and separated. These are interfaces to use, not defects to re-audit.

Inherited local results read:

* `05-knowledge/results/session267_gaussian_planar_20261111/gaussian/GAUSSIAN-NULL-PACKETS.md`, especially 267-G1/G2: finite-support charge and uniform marked-observer budgets. Conditional there on inherited NC2 where stated.
* `05-knowledge/results/session267_gaussian_planar_20261111/SYNTHESIS.md`, §§1–3: angular convolution, non-invariant source flags, finite rational-support certificates versus arbitrary finite prefixes.
* `05-knowledge/results/session267_gaussian_planar_20261111/meta/QUOTIENT-GUARDRAILS.md`: zero-character projection fails to commute with products.
* `05-knowledge/results/information_observers_collatz_20261014a.md`, §2: verifier-relative sufficiency is constancy on observer fibres; Fisher preservation for one tilt does not preserve every target.

Closest mechanism: observer-fibre sufficiency. Canonical hostile: scalar expectation forgets reciprocal charges. Corrected near miss: compactness gives limit points, not a unique limit. Least-used sidecar: contact-supported distributions and the difference between an observable algebra and its linear span. Concept board: compactness / observer fibres / products / contacts / phase / positivity.

## 1. Compactness is a uniqueness amplifier, not a uniqueness theorem

**PROVED (LM1).** Let K be a compact metric space and x_p in K. Let A be a family of continuous complex-valued observers. Assume A separates points of K and for every f in A the scalar f(x_p) has a limit. Then x_p converges to exactly one x in K.

Proof. Sequential compactness gives a cluster point x. Every cluster point y satisfies f(y)=lim f(x_p)=f(x), so separation gives y=x. If x_p fails to converge to x, a subsequence stays a fixed distance away; a convergent further subsequence yields another cluster point. Contradiction.

A more source-faithful formulation replaces scalar convergence by: every increasing subsequence has a further subsequence whose limit lies in a common subset K_0, and every two points of K_0 have identical observer values. Separation on K_0 suffices, provided all such further limits are in the same compact metric carrier. It is essential to compare hierarchies through a common carrier, not identify their separately presented label or Hilbert spaces by fiat.

**PROVED (LM1b, finite certificates at fixed resolution).** For every epsilon>0 there exist f_1,...,f_m in A and delta>0 such that

    max_j |f_j(x)-f_j(y)| < delta  implies  d(x,y)<epsilon   (x,y in K).

Proof. The compact set {(x,y):d(x,y)>=epsilon} is covered by the open sets {|f(x)-f(y)|>0}. Choose a finite subcover. The maximum of their absolute differences has positive minimum on that compact set; take delta to be that minimum. If the set is empty the assertion is trivial.

This is the continuum analogue of a verifier-sufficient finite statistic. A growing finite bank can suffice at each fixed resolution even when no one finite bank suffices exactly. No computable m or delta follows from qualitative compactness.

**CONDITIONAL program.** Put candidate YM unselected hierarchies in one metrizable locally convex distribution space; prove precompactness in a specified topology; choose observers that determine limits modulo agreed contact terms; prove candidate independence on those observers. LM1 then upgrades subsequential volume convergence to full convergence **in that quotient**. Recovering more requires a contact sidecar. The trusted existence statement alone does not supply this extra uniqueness assertion.

## 2. Ordered limits are a two-clock debt, not a free diagonal limit

**PROVED (LM2).** Write a_{N,p} in a metric space, lim_N a_{N,p}=b_p for each p, and b_p→b. For a chosen p(k)→infinity, the diagonal a_{N(k),p(k)}→b follows if d(a_{N(k),p(k)},b_{p(k)})→0. In particular, a uniform bound d(a_{N,p},b_p)<=e_N with e_N→0 suffices for every N(k)→infinity. Proof: triangle inequality.

**REFUTED shortcut, exact hostile.** Let a_{N,p}=1_{N<=p} for N,p>=1. Each fixed-p mesh limit is zero and the subsequent volume limit is zero. The diagonal N=p is identically one. This does not contradict either ordered limit.

For countably many scalar tests f_j, a diagonal can always be chosen *after* its target volumes: choose N(k)>=k such that the first k fixed-volume errors at p(k) are all <1/k. This proves convergence on that countable test family, not every test automatically. Extension to all tests needs an equicontinuity bound and a dense test family; an uncountable token set cannot silently be diagonalized.

**Bridge to S267.** A finite prefix does not bound a delayed pole; fixed-volume convergence does not bound a moving-volume error. In both cases the missing object is a uniform budget, not an additional finite census. A running coupling uniform in p is one source condition, but is not by itself a uniform error bound on every renormalized insertion.

## 3. Separated moments are a quotient by contact information

**PROVED (LM3).** Let Delta be the union of pairwise collision diagonals in (R^d)^n. The restriction map from distributions to the open complement of Delta has kernel exactly the distributions supported in Delta. This follows directly from the definition of distributional support: restriction is zero iff the distribution vanishes on every test supported away from Delta.

Thus separated observers distinguish a distribution only modulo contact-supported terms. When the target is already a smooth function on the separated locus, the observers do determine that function there: a smooth difference with zero pairing against every compact test is zero. These two statements are compatible.

**PROVED hostile across ALL orders.** Let W be centered real Gaussian white noise on Schwartz distributions over R^d, with covariance

    E[W(f)W(g)] = integral f(x)g(x) dx.

Compare its law with the deterministic zero field. Every odd moment distribution is zero. Every even moment distribution of W is the sum over pairings of products of delta(x_i-x_j); every summand is supported on a partial collision diagonal, hence vanishes upon restriction to the pairwise-separated configuration space. Therefore W and zero have identical separated n-point distributions for every n>=1, but different laws: E[W(f)^2]=||f||_2^2 versus zero. White-noise existence here is the standard Gaussian construction, not a YM model or a claim about its short-distance trace test.

**Surprise.** Even infinitely many perfect separated observers can miss an entire positive covariance, not merely an arbitrary signed counterterm. Positivity of probability laws does not remove the kernel. This is a continuum version of S267's missing source flag: completing an observer quotient does not complete the original object.

**CONDITIONAL repair.** Specify a contact extension class with a uniqueness condition, or supply all overlapping smeared moments/characteristic functionals. An OPE controlling separated approaches may constrain extensions but should not be declared to choose every delta-supported term without an explicit extension theorem.

## 4. Three positivity notions must not be interchanged

Classical positivity is L(P*P)>=0 for every polynomial P in commuting field variables, with * actual complex conjugation. OS reflection positivity is L((Theta P)*P)>=0 for positive-time P. Hilbert-space/operator positivity is omega(A* A)>=0 on a specified noncommutative star algebra, with domains needed for unbounded fields. These are different predicates and carriers.

**PROVED finite hostile (LM4).** Let Z be a symmetric random sign and (X_-,X_+)=(Z,-Z). Its probability law is invariant under swapping the two time sites. Classical positivity holds by expectation of absolute squares. For positive-time P=X_+, the OS form is E[X_-X_+]=-1, so OS positivity fails. Reflection invariance plus classical positivity is insufficient.

**PROVED Gaussian warning.** For circular Z with E|Z|²=1, E[Z^m]=0 for every m>=1 while E[Z*Z]=1. S267 null-power escape is about a complex bilinear moment functional, not zero norm, zero random variable, or zero operator. Its marked observer Q=conjugate(Z) restores the missing reciprocal charge immediately. This is why it is useful for detecting lost tests but cannot substitute for an OS argument.

**CONDITIONAL reconstruction routes, kept separate.** On a compact configuration space K, a normalized positive functional on C(K) gives a probability measure (Riesz representation). For a field on a nuclear test space, a normalized continuous positive-definite characteristic functional gives a measure on the dual under the usual Bochner–Minlos hypotheses. Neither conclusion follows merely from positive separated moments; nor does an operator/GNS state by itself identify a classical commuting field law. Conversely a classical law does not imply OS positivity, as LM4 shows.

## 5. Phase mixing is a barycentre, not phase identification

Let C be a compact metric phase space and h_j in C(C). The observed mixture map is nu↦(integral h_j dnu)_j.

**PROVED (LM5).** If the complex linear span of {1,h_j,conjugate(h_j)} is uniformly dense in C(C), this map is injective on Borel probability measures. Proof: equality on this span extends by uniform approximation to every continuous function; uniqueness of the representing measure gives equality.

It is not enough that the h_j separate phase points unless products are also available. On C={-1,0,1}, h(x)=x separates points, but

    nu=(-1:1/4, 0:1/2, 1:1/4),
    rho=(-1:1/3, 0:1/3, 1:1/3)

are distinct full-support probability measures with the same observed mean zero. Adding h² separates this pair. A unital self-adjoint **algebra** separating points is dense by Stone–Weierstrass, but an arbitrary separating observer family is not that algebra.

**Bridge to the exact Gaussian convolution.** In S267, expectations of U² and V² vanish separately while E[U²V²]=8 for U=H+iG,V=H-iG with independent unit-variance real Gaussians. The same lost-product mechanism explains why phase-point separation need not imply phase-measure separation. Full support supplies no uniqueness repair.

For the trusted hierarchy, characterContinuity and compactness are appropriate starting interfaces. To apply LM5 one must prove a density/product-generation assertion for its actual phase observables; a product of phase expectations is not automatically another correlation. Even uniquely identifying the unselected phase measure would not identify its distinguished selected phase or prove all phases equivalent.

## 6. A constructive bounded-observer measure theorem

**PROVED (LM6; compact product version).** Suppose for each finite set I of coordinate labels one has probability measures mu_{r,I} on [-1,1]^I, consistent under coordinate projection for each r. Assume every mixed monomial moment converges as r→infinity. Then there is a unique consistent limiting family mu_I, and mu_{r,I} converges weakly to mu_I for each I. For a countable coordinate set this determines a unique probability law on [-1,1]^N.

Proof. Compactness gives subsequential probability limits on every finite cube. Continuous monomials preserve their limiting values. Polynomials are uniformly dense in continuous functions on the cube, so any two cluster laws agree. LM1 gives full finite-dimensional convergence. Projection consistency passes to the limit by testing a continuous function of the projected coordinates. The standard countable product extension theorem supplies the final law and cylinder sets determine uniqueness.

This is a reusable version of 'moments construct a measure' with every sidecar visible: **bounded common carrier, all mixed moments including conjugates/repetitions, and projective consistency**. For unbounded renormalized fields, replace boundedness by tightness and enough uniform integrability to pass the required moments, plus a determinacy condition. Moment convergence alone on unbounded carriers does not justify those replacements.

**CONDITIONAL YM program.** Choose actual smeared lattice observables in a common coordinate system, prove joint tightness and consistency, control overlapping/repeated moments, and establish determinacy. Only then identify their law with the separated hierarchy and treat the contact extension. This is a constructive target parallel to S267's marked-observer compiler, rather than a claim that finite lattice probability measures automatically converge on an unspecified continuum space.

## Connection ledger and stopping boundary

| Source | Target/map | Preserved | Lost | Required sidecar / cheapest test |
|---|---|---|---|---|
| finite observer sufficiency | compact candidate→observer values | specified predicates | other fibres | separation; LM1b finite-resolution test |
| S267 uniform charge budget | two-clock mesh/volume approximation | fixed-index limit | moving-index rate | uniform error or adaptive diagonal; LM2 hostile |
| angular projection | separated distribution restriction | off-diagonal hierarchy | contact terms | contact extension; white noise versus zero |
| marked Gaussian moments | classical/OS/operator reconstruction | exact chosen functional | star/reflection/domain | LM4 two-site sign law |
| reciprocal-charge products | phase→barycentre observables | observed averages | phase law and selection | product/density test; three-point full-support laws |
| recursive observer compiler | consistent finite-dimensional laws | all retained joint tests | tightness/regularity | bounded cube theorem LM6, then genuine field sidecars |

All six lemmas are reusable elementary results or explicitly named standard representation steps. None proves volume uniqueness, a new continuum YM probability measure, universality, or JC(2). The most promising next task is **observer sufficiency on the actual candidate-limit class**, with separate contact and phase components, rather than another compactness-only extraction.
