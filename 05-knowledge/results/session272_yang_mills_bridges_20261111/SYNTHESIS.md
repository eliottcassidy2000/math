# Yang–Mills through observer completeness, phases, flux, and spectral visibility

**Session272, 11 November 2026. Research synthesis.** The supplied construction
is **CITED / TRUSTED**. New statements here are **PROVED abstract lemmas**,
**FINITE-EXACT toy computations**, or **CONDITIONAL programs**, as marked.
The five requested Yang–Mills extensions remain **OPEN** in this session.

## 1. The strongest organizing connection

The repository's Gaussian and tournament advances repeatedly show that the
right state is the one closed under the next operation. Session267 needs
marked Gaussian moments and source flags rather than a scalar expectation.
Session271 needs rooted resolvents or path-cover tensors rather than an
unrooted statistic. The Yang–Mills construction makes the same question much
more consequential: which extra coordinates make an existing hierarchy
sufficient for the next reconstruction step?

The answer is not one universal invariant. We identified distinct coordinates:

- **volume comparison:** a common carrier and convergence-determining tests;
- **probability:** contact terms, full mixed moments, positivity and tightness;
- **phase:** products of character observables, not merely their averages;
- **gap:** quantitative visibility of normalized low-energy states;
- **Wilson loops:** extended support, flux sectors, representations and scale;
- **particles:** the joint energy–momentum spectrum and isolated mass shells.

These are concrete mathematical objects with exact tests below. Their roles
cannot be exchanged just because the same word, such as spectrum or moment,
appears in both subjects.

## 2. What the trusted construction already supplies

The pinned source is
[3456e5e75cf8d4cd558955a5c7638fba2d68f464](https://github.com/yang-mills-group/OpenAI-YM-20261028/tree/3456e5e75cf8d4cd558955a5c7638fba2d68f464).
Twenty-three defining/interface files were read, with deeper reads of phase,
light-spectrum, heat, OPE, contacts, and fixed-volume profiles. No Lean rebuild
or source proof re-verification was undertaken.

The quantifier order is unusually useful: **one Wilson family is chosen before
all increasing volume sequences**, then an extracted sequence and hierarchy
are supplied for each. At suitable fixed volumes the mesh limit is full. The
remaining volume limit can refine. The unselected lattice-limit hierarchy is
a phase mixture; the selected Wightman theory is a character phase. These are
not the same limiting object by definition.

The deeper construction gives a reconstruction in every character, the same
graded null ideals, a unique zero-energy vacuum and scalar weak commutant in
each phase. For each fixed word combination A satisfying the source's
finite-support condition, constants M_A and M'_A work in every phase:

    ν_(i,A)(0<E≤u) ≤ M_A u²       (0<u≤1),
    0≤heat_(i,A)(t) ≤ M'_A t⁻²   (t≥1).

Here heat subtracts the full zero-energy projection.

The spectral interval is **energy**, not invariant mass. OPE coefficients are
shared across phases and truncation degrees; each finite order and admissible
finite collision tree has its estimate. Bounds survive bounded-total-variation
integration with varying geometry and phase. This is much stronger than a
pointwise OPE, but does not assert convergence of the infinite OPE at fixed
radius. Contact distributions extend separated correlations, with explicit
boundaries on positivity and uniqueness.

The detailed map, exact declarations and limitations are in
[SOURCE_READING.md](SOURCE_READING.md). Trusting these results leaves more
room for productive work, not less need to respect their quantifiers.

## 3. Six positive reconstruction bridges

[Limits and measures](limits_measure/README.md) gives proofs and exact controls.

**Compact observer theorem.** A sequence in a compact metric candidate space
converges if every observer in a continuous point-separating family converges.
At each prescribed accuracy a finite observer bank suffices. This supplies a
precise upgrade from subsequences once candidate-independence is proved, but
not the missing candidate-independence itself.

**Two-clock theorem.** Fixed-volume mesh convergence plus full volume
convergence permits a suitably chosen diagonal on countably many tests.
Arbitrary diagonal convergence needs a moving-volume error bound. The hostile
array 1_(N≤p) has mesh-first iterated limit zero but diagonal value one
(and the reverse iterated limit is one).

**Contact-kernel theorem.** Restriction to pairwise-separated configurations
forgets exactly distributions supported on collision diagonals. Remarkably,
zero field and Gaussian white noise have the same separated correlations at
*every order*, while their probability laws differ. This is positive,
probabilistic lost information, not just an arbitrary signed counterterm.

**Phase barycentre criterion.** A point-separating collection of character
observables need not determine a probability measure on phase space. Its
linear span must be dense, or a product-closed algebra must be retained. Two
different full-support laws on {-1,0,1} have the same mean of the separating
coordinate x; x² distinguishes them. This is the same lost-product mechanism
as Session267's E[U²]=E[V²]=0 but E[U²V²]=8.

**Bounded mixed-moment construction.** Consistent finite-dimensional laws on
[-1,1]^I with convergence of every mixed moment converge to unique laws and
construct a countable-coordinate measure. The unbounded field version needs
new tightness, uniform-integrability and determinacy inputs. This theorem
names an attainable intermediate construction instead of assuming a
continuum space exists already.

**Full-support phase rigidity.** A continuous nonnegative defect with zero
phase average vanishes in every phase. In particular zero variance of all
character-algebra generators forces a scalar algebra. A quantitative version
uses a lower ball-mass modulus. This converts the source's full support from
a descriptive feature into a rigidity tool; proving the required vanishing
is an additional obligation. See [phase rigidity](limits_measure/PHASE_RIGIDITY.md).

## 4. The mass-gap bridge is observability, not just irreducibility

[Spectral proofs](spectral_loops/PROOFS.md) give a sharp countermodel. In phase n,

    H_n=diag(0,ε_n), A_n=ε_n σ_x, B_n=σ_z, ε_n=2⁻ⁿ.

The vacuum is cyclic, its zero-energy space is one-dimensional, and the
observable algebra is M₂, so the commutant is scalar. For every *fixed word*
P, excited spectral weight is O_P(ε_n²); hence the quadratic low-energy and
inverse-square heat estimates hold uniformly in n. Nonetheless gap(H_n)→0.
Normalizing the visible excitation costs ε_n⁻¹. The example is finite
operator algebra, not a claimed Wightman model or a refutation of the source.

Two rigorous sufficient upgrades are available:

1. A **common positive exponential rate** for heat correlations on a total
   family gives a spectral gap. Constants and onset times may depend on the
   vector. Rates may not: a dense algebra of finitely supported vectors can
   have a different positive rate for each vector while the gap is zero.
2. A uniform visibility estimate

       ||1_(0,u](H)|| ≤ K Σ_j ||1_(0,u](H)x_j||²

   turns quadratic bounds M_j u² into a gap of at least
   min{1,(KΣM_j)^(-1/2)}. The real target is to control K using physically
   allowed probes, with scale and phase explicitly retained.

This is a literal logical bridge: positivity makes a projection's norm zero
or one, while quantitative observation prevents its nonzero part from hiding.
The source's localized commutator/detector argument is a concrete place to
look for additional coercivity. The present session does not prove it.

## 5. A complete loop laboratory: why a local area law has a global sidecar

In a finite two-dimensional Z₂ gauge disk, plaquette signs are independent
with mean q=tanh β. A loop enclosing a faces satisfies exactly

    E W(∂S)=q^a.

On a torus with F faces, the product of all plaquettes is constrained to +1:

    E W(∂S)=(q^a+q^(F−a))/(1+q^F).

The second term is the complementary spanning sheet. A noncontractible
fundamental loop has expectation zero when the two toron labels are freely
summed, and those labels are invisible to local plaquettes. The underlying
cochain map is literal: star flips/cuts in tournament switching become Z₂
vertex gauge changes, cycle pairings become Wilson signs, and H¹ stores what
curvature forgets. It is not a nonabelian Yang–Mills proof.

Exact enumeration covers all 4,096 link assignments on the disk and 256 on the
torus, every one of the 16 face subsets in each, plus transfer powers through
20. At q=1/3 the torus expectations for one and two faces are 15/41 and 9/41,
not 1/3 and 1/9. This is an inexpensive hostile test for any proposed local
factorization or area-law transfer.

The positive transfer object is the two-state parity kernel, with eigenvalues
1 and q. It composes; one scalar entry does not compose without its missing
channel. This recovers the path-cover/port-state mechanism in an exact gauge
toy. Moving to the source requires continuum loop approximation and
renormalization, nonabelian representation/path ordering, and physical
flux-sector coercivity uniform in its weak-coupling running regime.

## 6. An explicit ultraviolet coefficient with three different infrared worlds

[UV/IR separation](uv_ir/UV_IR_SEPARATION.md) constructs a positive continuous
mass-squared density

    ρ_UV(ds)=D/(16π²) · 1_(s≥e²) s²/(log s)² ds.

Integrating the four-dimensional free covariance against this measure gives
exactly the source's leading short-distance coefficient:

    lim_(r→0+) r⁸ log²(1/r) C(r)=3D/π⁴.

Any finite positive spectral measure supported below e² can be added without
changing that coefficient. Three choices give a gapped continuum without an
isolated mass, a gapless continuum, or an isolated mass shell below the UV
continuum. The constant follows from the exact Bessel integral
∫u⁶K₁(u)du=384; the asymptotic proof and an independent review are supplied.

These are positive Minkowski spectral covariance models with the displayed
separated-point Euclidean continuations. They do not provide a positive
Euclidean covariance on all tests or a continuum probability measure: the
full contact extension needs additional choices and positivity conditions.
Nor are they identifications with the source's physical stress tensor or full
hierarchy. They turn its explicit
continuous-mass boundary into a useful test family. A much sharper ultraviolet
coefficient cannot substitute for measuring infrared support. Source-specific
stress identities and higher correlations may distinguish them, and must be
retained if that is the intended next theorem.

## 7. Recovered cross-domain map, in sixteen connections

1. S267 Gaussian marked observers → full mixed/contact moment reconstruction.
2. S267 non-invariant source flags → admissible contact/positivity constraints.
3. S267 delayed poles → nonuniform mesh/volume diagonals.
4. S271 contextual states → convergence-determining observer banks.
5. S271 positive mixtures → phase-law versus phase-point distinction.
6. S271 rare-event/low-bit loss → UV accuracy versus IR support.
7. S271 rooted resolvents → stress/operator labels beyond a scalar correlator.
8. THM-2842 measured variance jets → quantitative spectral visibility.
9. THM-1975 path-cover channels → exact flux-parity transfer kernel.
10. THM-479 cut/holonomy coordinates → Z₂ gauge cochains and torons.
11. THM-174 Pfaffian square → covariance/higher-cumulant loss warning.
12. Compactness → finite-resolution observer certificates.
13. Distributional restriction → contact-term kernel and white-noise witness.
14. Full-support character measures → positive-defect rigidity.
15. Positive Laplace transforms → low-energy tails, with a gap upgrade only
    after a common rate or visibility condition.
16. Källén–Lehmann mass mixing → independently variable ultraviolet tails,
    mass gaps, and particle atoms.

Every connection has its carrier, preserved property, loss, extra condition,
and cheap decisive test in the lane ledgers. Historical results are used only
in their cited scope; in particular a VERIFIED finite tournament calculation
is not promoted to a universal compositional theorem.

## 8. Prioritized next research

The most promising next theorem is **source-compatible observability**:
construct controlled probes of every low-energy sector, retaining their
normalization as cutoff, volume and character vary. In parallel, isolate the
**candidate-limit observer algebra** and prove volume comparison estimates
for it, rather than adding another extraction argument.

For probability, test **overlapping-support positivity of the actual contact
prescription** before moment-determinacy work. For confinement, work in the
**charged flux sector**, not just the vacuum sector. For particles, study the
**joint translation PVM** and invariant-mass atoms after gap/sector control.

[OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) turns all five requested directions into
specified hypotheses, first tests, and stopping criteria. Nothing in this
packet claims those extra hypotheses for the four-dimensional construction.

**Incoming-research guardrail.** The concurrent Stein-module audit
(c440c4a4 / 182fd029, MISTAKE-764) repaired a claim that a conjectural ladder
was intermediate when it was equivalent to JC(2). No such claim is imported
here. A sufficient bridge in this packet may be as hard as its target; its
value is a precise interface and decisive tests, not an unsupported promise
of lower difficulty.
