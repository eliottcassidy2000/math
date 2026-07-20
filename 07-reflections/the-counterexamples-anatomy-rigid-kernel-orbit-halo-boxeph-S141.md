# The counterexample's anatomy: rigid kernel, orbit halo (boxeph-2026-07-19-S141)

Owner brief: build the explicit A₃ endomorphism; understand sporadic-vs-families as
deeply as possible; archaeology of the rationals; wild hypotheses welcome (labeled).

## A. The explicit proper self-embedding of A₃ (DELIVERED, fully machine-verified)

φ: A₃ → A₃,  φ(x_i) = F_i,  φ(∂_j) = Σ_k G[k][j] ∂_k,  G = (JF)⁻¹ = −½·adj(JF)
(all 77 monomials of G dumped in the frozen .out). Verified exactly: J·G = I (gives
[φ∂_j, φx_i] = δ_ij) and all nine identities [D_i, D_j] = 0 (gives [φ∂_i, φ∂_j] = 0)
— φ preserves every Weyl relation. Injective since A₃ is simple; NOT surjective:
a preimage of x_i would yield P with P(F) = x, i.e. a polynomial left inverse,
contradicting the verified collision. **The first explicit Keller-type proper
self-embedding of the third Weyl algebra — the Dixmier conjecture's dim-3 death
made constructive.** (DC₃ was already false abstractly via ¬JC₃; this is the object.)

## B. The map's generative anatomy (all verified exactly)

Let u = 1+xy. Then 4+3xy = **1+3u**, and:
- **z-linearity**: F = a(x,y) + z·b(x,y) with b = (u³, 3xu², −x³) — Pascal row
  (1,3,3,1) with a sign-twisted tail; a = (y²u(1+3u), y + 3xy²(1+3u), 5x−3xu).
- **ℂ*-equivariance**: F(t⁻¹x, ty, t²z) = (t²F₁, tF₂, t⁻¹F₃) — source weights
  (−1,1,2), target (2,1,−1) (a weight REVERSAL); the S140 σ/τ pair is t = −1.
- **The collision curve**: for every t ≠ 0, {(0,0,−t²/4), (±t⁻¹, ∓3t/2, 13t²/2)}
  ↦ (−t²/4, 0, 0) — verified at t = 2, 3, ½. The famous three points are the t=1
  slice: the counterexample PHENOMENON is a 1-parameter family inside one map.
- **The descent**: on invariants (s,r) = (xy, x²z) → (vw, uw²), F induces the plane
  map (P,Q) (coefficients in the .out) whose Jacobian is **non-constant**
  (8 − 24s + 18s² − 8r + 12sr + 2r²): Keller-ness dies exactly at the torus-fixed
  defect — a precise mechanism for why this construction needs dimension 3.
- With S140: generic degree 3, **full S₃ monodromy** (Chebotarev ½/⅙/⅓), σ = a
  transposition; fibers only {1,3} at all tested p.

## C. Sporadic or family? The exact local answer

- **Trivial closure** (families for free): Aut∘F∘Aut orbits, stabilization F⊕id
  (⟹ cubic-linear counterexamples in higher dim via Bass–Connell–Wright normal
  form), products. These are halo, not new kernels.
- **The equivariant Keller tangent at F** (weight-respecting, degree-capped
  perturbations with d/dε det = 0; exact linear algebra): **dimension 5** of 20.
- **The orbit tangent** (equivariant source reparameterizations JF·V and target
  reparameterizations W∘F, projected to the same class): **dimension 5, and the
  spans coincide.** Every infinitesimal equivariant Keller deformation is a
  reparameterization: **F is locally rigid modulo equivalence in its class.**
- Cross-check: all 5 tangent representatives FAIL to integrate as raw monomial
  paths (det breaks at 2nd order, both signs) — consistent with pure-orbit tangents
  in a curved parameterization. VERDICT: **sporadic-flavored — a rigid kernel** at
  this support/weight class. Honest scope: other weight classes, higher degree
  caps, and non-equivariant deformations remain open (named leads).

## D. Archaeology of the rationals + wild hypotheses (WILD = labeled speculation)

- **13/2**: the collision z-coefficient AND death-star's clustered-floor threshold
  (M < 1/13 ⟹ ρ > 6.5, S58h); 13 = the tight AP length of the LRC thread. WILD:
  both are "the first obstruction above the doubling floor."
- **deg f₁ = 7**: the apex-7 of the entire project (14 = 2·7 walls; the ×7 gate's
  1-of-7 rigidity). WILD: degree-7 is the minimal ℂ*-equivariant realization for
  the same reason 7 > √14 blocks the polynomial shortcut in S-T.
- **det = −2**: the factor-2/doubling motif (THM-882 doubling law F_N = 2G_N; the
  ×2 gate; JC₂ₙ ⟹ DCₙ). The Keller-normalized form diag(−½,1,1)∘F has det 1.
- **(1,3,3,1)**: b = (u³, 3xu², −x³) is a Pascal row with one sign flipped and one
  slot suppressed — WILD: the staircase-triangle row of the tournament project,
  with the flip playing the complement involution.
- **Weight reversal (−1,1,2) ↔ (2,1,−1)**: WILD: the transpose/complement duality
  of the metagraph, appearing as source↔target duality of the counterexample.
- **S₃ on 3 sheets with a ℤ/2 inside**: the two tournaments on 3 vertices (3-cycle
  vs transitive) — WILD but crisp: the monodromy's cycle types are exactly the
  orbit census our SC/NS machinery counts.
- **THE LIFTING PROGRAM (serious conjecture):** ℂ*-equivariant z-linear Keller maps
  of reversed weights are classified by plane data (P,Q) with prescribed
  torus-fixed defect; the rigid kernel F is the minimal solution, and NEW kernels
  (if any) correspond to other defect classes — the right place to hunt further
  counterexamples and the right no-go frame for JC₂ (a free ℂ*-action cannot exist
  on ℂ², so THIS mechanism cannot descend — a provable-looking JC₂ obstruction).

Files: jacobian_weyl_endo_and_structure_boxeph_S141.py + .out (all claims frozen).
