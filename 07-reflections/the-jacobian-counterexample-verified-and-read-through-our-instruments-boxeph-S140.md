# The Jacobian counterexample: verified exactly, then read through our instruments
**boxeph-2026-07-19-S140** · owner communication: JC recently falsified; investigate
connections to this project and to Dixmier; "tournaments and their recursion is
reminiscent of the n vs 2n relationship between Jacobian and Dixmier."

## 1. Verification (machine-checked, self-contained — provenance-independent)

F(x,y,z) = ((1+xy)³z + y²(1+xy)(4+3xy), y + 3x(1+xy)²z + 3xy²(4+3xy), 2x − 3x²y − x³z).

Exact rational polynomial arithmetic (script + frozen out, S140):
- **det JF ≡ −2** identically (symbolic 3×3 determinant of the 9 partials);
- **F(0,0,−1/4) = F(1,−3/2,13/2) = F(−1,3/2,13/2) = (−1/4,0,0)** — three distinct
  preimages. A Keller map that is not injective: **JC₃ is FALSE.** These two facts
  alone suffice; nothing else needs to be trusted.
- Stabilization (append identity coordinates): **JC_n false for every n ≥ 3.**
  JC₂ is untouched and is now the last Jacobian island.

## 2. The map's own structure (new measurements, this session)

- **ℤ/2-equivariance (exact):** F∘σ = τ∘F with σ = diag(−1,−1,1), τ = diag(1,−1,−1).
  The collision fiber over the τ-fixed point (−1/4,0,0) is σ-organized as
  {fixed point} ∪ {2-orbit}: **3 = 1 + 2.**
- **The mod-p spectrum (full enumeration, p = 3..29):** fiber sizes are ONLY 1 and 3
  at every prime; the proportions of (1-fibers, 3-fibers, missed points) converge to
  **(1/2, 1/6, 1/3) — the class equation of S₃ read by Chebotarev.** So F is
  generically **3-to-1 with full S₃ monodromy**: a non-Galois étale-away-from-nothing
  triple self-cover of affine 3-space (fibers drop only by sheets escaping to
  infinity — non-properness, as étaleness forces). Since étale fibers are bounded by
  the geometric degree, the exhibited fiber is exactly {3 points}. The σ-involution
  realizes a transposition sheet-swap over the τ-fixed locus.
- **Open exact law:** T(p) = #(3-fibers) = 3, 18, 51, 205, 342, 776, 1089, 1947,
  3934 for p = 3..29 — T(p)/p³ → 1/6 with a structured secondary term; the exact
  point-count law of the collision threefold {F(a) = F(b), a ≠ b} is a fresh,
  well-posed object (its class in the Grothendieck ring; is T quasi-polynomial?).

## 3. The Dixmier cascade (standard implications, now with a live input)

- **DC_n ⟹ JC_n** (classical; a Keller map F lifts to the Weyl endomorphism
  x_i ↦ F_i, ∂_j ↦ Σ_k ((JF)⁻¹)ᵀ-entries·∂_k, well-defined since det JF is a unit
  and the formal-inverse chain rule gives the commutation). Contrapositive with §1:
  **DC₃ is FALSE, hence DC_n false for all n ≥ 3** (stabilize the endomorphism by
  identity on the extra Weyl generators).
- **Explicit corollary worth constructing (lead):** the induced endomorphism of A₃
  from THIS F is injective (Weyl algebras are simple) but **not surjective** — the
  first explicit proper Weyl-subalgebra embedding of A₃ into itself of Keller type.
  Writing it down explicitly (the (JF)⁻¹ entries are computable — det = −2 makes
  them ½·(cofactors)) is a mechanical, publishable-grade exercise.
- **The surviving islands and their new logic:** JC₂ open; DC₁ and DC₂ open.
  Tsuchimoto/Belov-Kanel–Kontsevich give JC_{2n} ⟹ DC_n; with JC₄ now false that
  implication says nothing about DC₂ — **DC₂ is decoupled from the dead part of the
  Jacobian hierarchy and hangs only below JC₂** (¬JC₂ would kill it via DC₂ ⟹ JC₂).
  The BKK "stable equivalence" JC ⟺ DC is now the equivalence of two false
  statements — its live content migrates entirely to the n ≤ 2 islands.

## 4. The n ↔ 2n web, and where this project genuinely touches it

The BKK/Tsuchimoto bridge runs through **reduction mod p**: A_n over 𝔽_p is Azumaya
of rank p²ⁿ over its center 𝔽_p[x_i^p, ∂_i^p] — a polynomial ring in **2n**
variables whose generators are p-th powers. An endomorphism of A_n induces a
Poisson endomorphism of this 2n-dimensional center; that is the n → 2n doubling.
Structural rhymes with this project (labeled: analogies at the method level,
theorems only where cited):

- **The p-power center ↔ our degeneracy channels.** The center's coordinates exist
  only mod p (the inseparable/Frobenius direction). In our gate framework the same
  role is played by the all-divisible channels (dilated cores, S134/S139 blocker
  strata): the mod-p world acquires a "fold-down" direction with no continuum
  counterpart, and in both settings the interesting objects are the ones
  transverse to it (our primitive strata; their Azumaya locus).
- **The doubling recursions.** The Cayley–Dickson tower of the tournament project
  (ℝ→ℂ→ℍ→𝕆, tournaments at n = 2^k+1, each level trading a property) is our native
  n ↔ 2n ladder; the Weyl-center doubling trades noncommutativity for dimension the
  same way each CD step trades a property for dimension. The owner's instinct is
  right at the level of *mechanism*: both hierarchies convert structure into
  doubled coordinates under a limit/reduction.
- **The ℤ/2 bookkeeping.** The counterexample's collision census (3 = 1 fixed + one
  2-orbit under σ) is exactly the SC/NS decomposition discipline of the merged
  metagraph G_n/ℤ₂ — involution-fixed objects plus swapped pairs. Our census
  machinery applies verbatim to collision correspondences of polynomial maps.
- **The instrument transfer (the real one).** The mod-p fiber census that read off
  S₃ in §2 IS the improperness-mod-p ladder of the LRC thread pointed at a
  different continuum statement. One tool, two conjectures: reduce mod p, enumerate
  the finite shadow exactly, read the invariant (monodromy there, gate laws here),
  and treat divisible channels as a separate stratum. That instrument just
  extracted the counterexample's hidden Galois group in two seconds of arithmetic.

## 5. Leads filed

(i) exact T(p) law / the collision threefold's zeta data; (ii) the explicit
non-surjective A₃-endomorphism from F (mechanical construction); (iii) the JC₂
island: does the σ-equivariant *method* behind F admit a 2-variable analogue, or
is there an obstruction (this is where our involution-census machinery could do
real work); (iv) engineering: package the mod-p Chebotarev fingerprinter (map →
monodromy statistics) alongside the tournament H-spectrum fingerprint — same
product family. Honest boundary: beyond §4's method-level bridges, no
object-level theorem connecting tournaments to JC/DC is claimed anywhere above.
