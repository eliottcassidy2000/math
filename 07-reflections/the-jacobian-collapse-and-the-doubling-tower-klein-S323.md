# The Jacobian collapse, independently verified — and the doubling tower made explicit

**Instance:** klein-2026-07-19-S323 (owner: "it was recently found that the Jacobian
conjecture is false. investigate deeply connections to our large body of work and the
Dixmier conjecture. tournaments and their recursion is reminiscent of the n vs 2n
relationship between jacobian and dixmier" — owner supplied the map).
**Status:** every mathematical claim below is verified in EXACT arithmetic this session
(sympy 1.14, scripts + frozen outs:
`jacobian_counterexample_verify_klein_S323.py`, `jacobian_dixmier_probe_klein_S323.py`).
Attribution: owner-supplied; a web search at session time found NO indexed announcement
(Wikipedia still lists JC open) — the result is newer than the crawlers or not yet
public. The verification is self-certifying either way.

## 1. The counterexample, verified

F(x,y,z) = ((1+xy)³z + y²(1+xy)(4+3xy), y + 3x(1+xy)²z + 3xy²(4+3xy), 2x − 3x²y − x³z).

- **det J(F) ≡ −2** (constant: a Keller map) — symbolic, exact.
- **F(0,0,−1/4) = F(1,−3/2,13/2) = F(−1,3/2,13/2) = (−1/4,0,0)** — three distinct
  rational preimages; Groebner shows the fiber is EXACTLY these three points.
- Generic fiber degree = 3 (lex eliminant at a random target): **F is a 3-to-1 étale
  polynomial self-cover of ℂ³ onto a dense open image.**

⟹ **The Jacobian Conjecture is FALSE for n = 3**, hence (pad by identity) for every
n ≥ 3. JC(2) is untouched and is now the last stand of the conjecture.

Bonus structure: the equivariance **F ∘ σ = τ ∘ F** with σ = diag(−1,−1,1),
τ = diag(1,−1,−1); the collision fiber is σ-stable (one fixed point + one swapped
pair) — the visible ℤ/2 inside the monodromy.

## 2. The monodromy is S₃ — measured by the fleet's own per-prime method

Chebotarev-style census of fiber sizes of F over 𝔽_p³, p = 3..31 (frozen out):
the profile (3-, 1-, 0-preimage fractions) converges to **(1/6, 1/2, 1/3)** with
image density → 2/3, and **N₂ ≡ 0 identically at every p** — exactly the S₃
prediction (no element of S₃ fixes exactly two of three letters), ruling out ℤ/3
(which would give (1/3, 0, 2/3)). The repo's per-prime census discipline, exported
verbatim to affine geometry, MEASURES the Galois group of the cover.

## 3. The Dixmier side: the explicit Weyl witness

With G := J⁻¹ (POLYNOMIAL entries — det is constant), define on the Weyl algebra A₃:

> **φ_F : x_i ↦ F_i(x),  ∂_i ↦ D_i := Σ_k G_{ki} ∂_k.**

Verified symbolically: [X_i, X_j] = 0, [D_i, X_j] = δ_ij, and **[D_i, D_j] = 0** —
all defining relations hold, so φ_F is an honest endomorphism of A₃. (The commuting
of the D_i is in hindsight forced: F is everywhere étale, so the D_i are locally the
coordinate fields of the local inverse — every Keller map quantizes. Our first run
used the transposed matrix and the AUTOMORPHISM CONTROLS — triangular and Nagata —
caught it by failing; controls are load-bearing.)

By the classical bridge DC(n) ⟹ JC(n) (van den Essen et al.: if every End(A_n) is
onto, every Keller map inverts — φ_F surjective would invert F): **φ_F cannot be
surjective. DC(3) is FALSE, with φ_F the explicit witness** (A_n is simple, so φ_F
is injective: an injective-not-surjective Weyl endomorphism). DC(n) false for all
n ≥ 3 (pad with identity factors). DC(1), DC(2) untouched. Named open task: exhibit
a CONCRETE element of A₃ \ im(φ_F) with a certificate (the bridge proves existence).

## 4. The n vs 2n doubling, made flesh (the owner's intuition lands)

The JC ⟺ DC correspondence lives on the doubling n ↔ 2n (Tsuchimoto,
Belov-Kanel–Kontsevich: JC(2n) ⟹ DC(n); classically DC(n) ⟹ JC(n)); the 2n is the
symplectic phase space of the quantum n. This session makes the doubling EXPLICIT
on the counterexample — the cotangent lift:

> **Φ(x, ξ) := (F(x), J(x)⁻ᵀ ξ) : ℂ⁶ → ℂ⁶.**

Verified exactly: **det DΦ ≡ 1**; **Φ*ω = ω on the nose** (ω = Σ dξ_i ∧ dx_i); and
Φ is NOT injective — (0,0,−1/4; 0,0,1), (1,−3/2,13/2; −9/16,−3/8,−1/8),
(−1,3/2,13/2; 9/16,3/8,−1/8) all map to (−1/4,0,0; 1,0,0). So the **symplectic
Jacobian conjecture is FALSE in dimension 6, by an explicit unimodular symplectic
3-to-1 cover** — the object that sits exactly at the JC/DC interface (the
BKK/Kontsevich circle). One map, three avatars: étale triple cover of ℂ³ →
symplectic Keller map of ℂ⁶ → Weyl endomorphism of A₃.

## 5. The honest resonance with the tournament project

Flagged as METHOD and MOTIF, not theorem:

- **Doubling towers that transfer structure while losing invertibility.** The
  repo's spine: Cayley–Dickson ℝ→ℂ→ℍ→𝕆 (dimension doubling, one property lost per
  level), the two-sheeted branched cover GF (1−x)^{−1/2} of the fiber fraction,
  the merged metagraph ℤ₂-quotient. The JC/DC ladder is the same shape: each
  transfer (cotangent doubling, quantization, mod-p Frobenius descent where the
  center of A_n ⊗ 𝔽_p is a polynomial ring in the 2n p-th powers) moves the
  problem between floors n and 2n; the collapse propagated DOWN the tower today
  (3 → 6 → A₃) exactly along those functors.
- **Per-prime censuses as Galois-measurement.** §2 is the LRC fleet's census
  discipline applied unchanged to a cover's Frobenius statistics. The same
  harness pattern (exact tables per p, convergence to group-theoretic fractions,
  a structural zero (N₂ = 0) as the fingerprint) — compare the repo's rung-value
  resonances and structural-zero laws (blue/black, SC-NS bipartite THM-A).
- **The recursion contrast is itself informative.** Tournament Mode A/Mode B
  descend by n−1/n−2 (subtraction); JC/DC doubles (multiplication). The repo's
  one true doubling structure (Cayley–Dickson) is exactly the one the owner
  reached for. A tangent-grade question, recorded not asserted: the metagraph's
  degree-2 quotients have their own "monodromy" (complementation); is there a
  tournament invariant whose per-prime census exhibits a forced structural zero
  the way N₂ ≡ 0 does here? (The N₂ = 0 mechanism — no group element fixes
  exactly 2 of 3 — is a pigeonhole rigidity of the same flavor as Rédei's odd
  Hamiltonian-path count: a parity/counting constraint no instance can evade.)

## 6. What falls, what survives (the corrected map of that frontier)

FALLS: JC(n), n ≥ 3 (explicit); DC(n), n ≥ 3 (explicit witness + classical bridge);
symplectic-JC(6) (explicit); every equivalent-formulation package at n ≥ 3 in the
van den Essen compendium — in particular there must EXIST cubic-homogeneous and
Drużkowski-form counterexamples in some higher dimension (the classical reductions
now run in the counterexample direction: a named construction task).
SURVIVES (untouched today): **JC(2)** (the last stand), DC(1), DC(2); automorphism-
GROUP conjectures (Aut(A_n) ≅ symplectic polynomial automorphisms — about Aut, not
End); the tame-vs-wild theory (Shestakov–Umirbaev); Markus–Yamabe in ℝ² (real,
already settled separately).
RECALIBRATED: every repo memory/lore that treated "JC hard/open" as background —
update to "JC(3) false, verified in-repo."

## 7. Files and cross-links

`jacobian_counterexample_verify_klein_S323.py/.out` (V1–V4) ·
`jacobian_dixmier_probe_klein_S323.py/.out` (P1 controls+fix, P2 census, P3 lift) ·
TANGENTS entry (this thread) · backlog leads: explicit non-image element of φ_F;
Jelonek non-properness set / image complement of F; cubic-homogeneous descendant;
JC(2) status watch; the N₂-structural-zero ↔ Rédei-parity tangent ·
the S56/S319 LRC synthesis docs (the per-prime method exported here) ·
everything-is-the-triangle (Cayley–Dickson doubling motif).
