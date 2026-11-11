# A logarithmic ultraviolet signature does not select an infrared spectrum

**Session272. PROVED IN NOTE (positive spectral covariance family); CONDITIONAL
as a diagnostic for the source's physical stress-trace hierarchy.** This is an
explicit construction of a useful countermodel, not a challenge to the trusted
source and not a new Yang–Mills solution. No external novelty claim is made.

## 1. Why this is a useful bridge

The trusted `TraceBridge.lean` / `YangMills.Hierarchy.TraceTest.shortDistance`
declaration in `ChallengeSpec.lean` specifies a
centered grade-four two-point function with

    lim_(r→0+) r^8 log²(1/r) C(r) = 3D/π⁴,

where D is the real Lie-algebra dimension. The README expressly does not
exclude all continuous-mass generalized free fields. We turn that boundary
into an explicit family: exactly this leading coefficient can coexist with a
positive gap and no isolated mass, with no gap, or with an extra isolated mass.
This makes the ultraviolet/infrared distinction testable rather than verbal.

The internal analogy is to Session271's very accurate large-scale census
forecast, which does not recover exact low bits or rare-event information, and
to Session267's Gaussian scalar observer, which does not determine the
operation-compatible state. Here the lost coordinate is spectral support at
bounded mass, not precision in the ultraviolet coefficient.

## 2. Definition and exact normalization

For s>0 and r>0 use the four-dimensional free Euclidean covariance

    Δ_s(r) = sqrt(s) K_1(r sqrt(s))/(4π² r).

Let D>0 and c=D/(16π²). Define the positive ultraviolet spectral density

    ρ_UV(ds) = c · 1_(s≥e²) s²/(log s)² ds.

For any finite positive measure σ supported in [0,e²), put

    C_σ(r) = ∫ Δ_s(r) [ρ_UV(ds)+σ(ds)],

where Δ_0(r)=1/(4π²r²). All pointwise integrals converge for r>0 because
K_1 decays exponentially at infinity. The spectral measure has polynomial
growth, so the positive mass-shell superposition defines a tempered
Minkowski two-point distribution and its generalized free Wightman field.
The displayed Euclidean function is its separated-point continuation. It
admits distributional extensions across coincidence, but those extensions
require contact choices: the naive Euclidean Fourier integral
∫ρ(ds)/(p²+s) diverges for this density. In particular, no ordinary positive
Euclidean covariance on all test functions, continuum probability measure,
or positive contact prescription is asserted here. The proof below uses
only the separated-point spectral integral, not a physical stress tensor.

**Proposition U1.** For every such σ,

    lim_(r→0+) r^8 log²(1/r) C_σ(r) = 3D/π⁴.

**Proof.** Set L=log(1/r), u=r sqrt(s). Then

    C_UV(r) = c/(8π² r^8) ∫_(er)^∞ u^6 K_1(u)/(L+log u)² du.

Multiplying by r^8 L² leaves the integrand factor
L²/(L+log u)². Take 0<r<e^(−2), so L>2 and er<sqrt(r).
On u≥sqrt(r), L+log u≥L/2, so this factor is at most 4.
For each fixed u>0 it tends to 1. Apply dominated convergence on (0,∞)
to the integrand with indicator 1_(u≥sqrt(r)), dominated by
4u^6 K_1(u). On er≤u≤sqrt(r), set u=r v and use
K_1(u)≤1/u (proved below). Its contribution to the rescaled integral is at
most

    L² r^6 ∫_e^(r^(-1/2)) v^5/(log v)² dv
       ≤ L² r^6 ∫_e^(r^(-1/2)) v^5 dv
       ≤ (L²/6) r³ → 0,

since log v≥1 throughout this interval.

The Bessel integral is exactly

    ∫_0^∞ u^6 K_1(u) du = 384.

Hence the limit is 48c/π²=3D/π⁴. For the added compact spectral measure,
uK_1(u)≤1 yields C_σ(r)≤σ([0,e²))/(4π²r²), so its rescaled
contribution is O(r^6 L²) and vanishes. QED.

For completeness the two Bessel facts can be obtained without asymptotic
fitting. The positive integral representation

    K_1(u)=(u/4)∫_0^∞ t^(-2) exp[-t-u²/(4t)] dt

gives the moment 384 by Tonelli and the substitutions v=u²/(4t): the
inner ∫u^7 exp[-u²/(4t)]du is 768t⁴, then the outer integral is
192Γ(3)=384. Substituting v=u²/(4t) in the representation also gives
uK_1(u)=∫_0^∞ exp[-v-u²/(4v)]dv≤1.

## 3. Three exactly different infrared profiles

The proposition holds for all of these choices:

1. **Gap and no isolated mass:** σ=0. The mass-squared support starts at e²
   and is continuous there. The generalized-free one-particle space has a
   continuous mass spectrum starting at e; it has no isolated mass shell.
2. **No gap and no isolated mass:** σ(ds)=1_(0<s<1)ds. Positive spectral
   weight occurs arbitrarily near zero and there is no atom at any mass.
3. **Isolated particle mass:** σ=wδ_(m²), 0<m<e, w>0. There is an isolated
   atom at m² below the continuous ultraviolet spectrum. In the generalized
   free Fock realization this is a stable one-particle sector. To call it
   isolated from *all* multiparticle energy-momentum spectrum one also keeps
   the usual multiparticle threshold distinction; that is separate from the
   scalar covariance atom statement.

They are not the same theory. They share one stringent, non-power-law
short-distance two-point signature. The source's stress-tensor identification,
all higher correlations, Ward/translation properties and the Wilson origin
may distinguish them. Nothing here asserts that any one is its hierarchy.

A Gaussian *field* can have a nontrivial continuous mass spectrum. The source's
finite specified Wick-profile exclusion is therefore not permission to label
all Gaussian-field countermodels impossible. Conversely, supplying a generic
positive spectral kernel is not supplying the source's physical stress tensor.
Both directions of that distinction are essential.

## 4. From a non-freeness test to a discriminating test design

The existing trace diagnostic asks for a logarithmic ultraviolet coefficient.
The complementary measurements must retain a different coordinate:

- for a mass gap, test support of every energy measure uniformly across a
  dense/cyclic observable class, as in the spectral lane;
- for an isolated particle, test an atom or isolated mass-shell projection,
  not merely positivity of a two-point function or an exponential envelope;
- for the physical stress trace, retain the translation/stress identities
  and its correlations with other operators, not just its scalar C(r);
- for non-Gaussianity beyond the finite Wick test, retain connected higher
  cumulants in the correct field and phase, together with the contact
  prescription used to define them.

The tensor/operator label and its lawful identities act like the rooted port
sidecar in Session271: the same scalar spectrum at a forgotten port can lead
to different compositions. This is an analogy with a specified lost coordinate,
not a map from tournaments to Yang–Mills fields.

## 5. Numerical probe, not numerical proof

`check_uv_ir.py` evaluates the dimensionless ratio

    R(L)=(1/384)∫_(exp(1−L))^∞ u^6 K_1(u)[L/(L+log u)]²du.

It approaches 1. For efficiency and numerical stability the integral is
truncated to [exp(−L/2),100] and the omitted interval has the analytic bound
L² exp(−3L)/2304; the upper-tail error is separately bounded by a positive
one-dimensional integral and reported. These are numerical diagnostics of an
analytically proved limit. The output reports working precision and quadrature
sensitivity and does not claim directed interval arithmetic.

## Sources and inheritance

- Trusted source pinned at `3456e5e75cf8d4cd558955a5c7638fba2d68f464`:
  [README](https://github.com/yang-mills-group/OpenAI-YM-20261028/blob/3456e5e75cf8d4cd558955a5c7638fba2d68f464/README.md),
  [TraceBridge](https://github.com/yang-mills-group/OpenAI-YM-20261028/blob/3456e5e75cf8d4cd558955a5c7638fba2d68f464/TraceBridge.lean),
  [ChallengeSpec](https://github.com/yang-mills-group/OpenAI-YM-20261028/blob/3456e5e75cf8d4cd558955a5c7638fba2d68f464/ChallengeSpec.lean).
- [Session271 ensemble proof](../../session271_predictive_tournaments_20261111/ensemble/PROOF.md): positive mixture transport and rare-event amplification.
- [Session271 synthesis](../../session271_predictive_tournaments_20261111/SYNTHESIS.md): rooted resolvent observer and recursion boundaries.
- [Session267 Gaussian/planar synthesis](../../session267_gaussian_planar_20261111/SYNTHESIS.md): scalar projection versus operation closure.

**Open transfer:** identify a source-compatible, renormalization-stable
observable that couples with a controlled lower bound to every low-energy
sector. The ultraviolet coefficient alone cannot perform that task.
