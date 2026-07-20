# The salvage catalog: what IS true in the ruins of JC — and the trisection reading of the counterexample

**kind-pasteur-2026-07-20-S128c99** (HYP-8130, THM-1335) · owner: *"pull creatively
from other agents and work creatively and fully understand every aspect of each
counterexample of the jacobian conjecture, or to find sharper versions of it that
are true."*

## 1. The counterexample census, one day in

Explicit counterexamples now in the repo: **F** (the seed); **G1 = F∘(x, y+x², z+y)**
(radical 1+xy+x³) and **G2 = F∘Nagata** (klein-S326 — honest conjugates, same
essential class, same modulus m); **F∘F** (field degree 9, the first genuinely
inequivalent descendant); and the ladders F^{∘m}. Every rigidity probe agrees the
degree-3 seed is unique: opus's deformation theory (tangent 5 = 1 trivial + 4
obstructed), death-star's chart rigidity + W3 refutation, boxeph's lifting
uniqueness, my design-box uniqueness — now extended to the κ-free/S-free box where
the entire solution variety decodes into scalings, translations, and the z-shift.
"Every aspect of each counterexample" today means: every aspect of F, transported.

## 2. THE SALVAGE CATALOG — sharper versions of JC that are TRUE

(each entry: statement — status — owner)

1. **Injective ⟹ automorphism.** Classical (Białynicki-Birula–Rosenlicht /
   Ax–Grothendieck). The failure mode is exactly non-injectivity of degree ≥ 3.
2. **Proper ⟹ automorphism.** Classical: a proper étale map of simply-connected
   ℂⁿ is a trivial cover. Sharpened here: the minimal failure of properness is a
   HYPERSURFACE, and for F explicitly the irreducible quartic {L=0} — which is
   also a CUSP {P³ + E² = 0} (THM-1335). klein's H6: Jelonek degree as a
   monodromy invariant.
3. **Galois ⟹ automorphism.** Campbell 1973 — now a special case of:
4. **The Smith selection rule** (klein-S325): étale self-covers of ℂⁿ have
   trivial deck group; degree 2 is IMPOSSIBLE; monodromy needs self-normalizing
   point stabilizers. d=3 forces S₃ (four fleet proofs = one topological
   necessity); d=4 allows only A₄/S₄; dihedral iff odd degree; Frobenius always.
   THE map sits in the minimal allowed cell.
5. **Degree-3 arithmetic is abelian-by-quadratic and fully decidable** (THM-1310
   + 1335): fiber rationality is decided by ONE quadratic character χ(−L)
   (0/8316 violations), and the resolvent field ℚ(√−L) has conductor = the
   Jelonek set. Conjecture (sharp, falsifiable): *for every degree-3 Keller
   cover, the quadratic resolvent of the fiber cubic ramifies exactly on the
   asymptotic variety* — arithmetic and geometry of the failure are ONE object.
6. **Trace-polynomiality / the polynomial centroid** (THM-1335, NEW): all three
   coordinate fiber-traces of F are polynomial — (0, 3b/2, the explicit
   weight-(−2) quartic) — though the full trace form has L-poles from degree 2.
   Conjecture: every Keller map admits trace-polynomial coordinates. A
   "sharper JC that is true" in integrated form: *the fiber-AVERAGE of the
   identity map is polynomial* (rank 2 here, image in {x=0}).
7. **Equivariant JC is TRUE in the relevant weights** (my c98 §8 + death-star
   THM-1325 torus-blindness + boxeph weight-3 emptiness): for weights
   (1,−1,1−k), k≠3, the linear part is forced singular — no Keller maps at all
   (vacuous JC); at cubic-linear equivariant dim ≤ 3, all maps injective
   (death-star). The counterexample threads the UNIQUE equivariant window k=3.
8. **The order-{1,3} line-congruence conjecture** (mine, evidence accumulating):
   a z-AFFINE Keller map of ℂ³ has field degree 1 or 3. Evidence: every k≠3
   design box empty (c98 + today's refined κ-free/x²u-stratified boxes at k=4,5
   — verdict recorded in the .out), boxeph's "kernels only at z-weight 2",
   Smith's d=2 exclusion. If true: degree-4 seeds (A₄/S₄ per Smith) require
   z-degree ≥ 2 — "the 2-jet architecture" — a sharp structural dichotomy.

## 3. The trisection reading (the "fully understand" capstone)

The master identities (THM-1335 §2) compress F to two lines. Scaling the u-cubic
to 4T³ − 3T = m gives the CHEBYSHEV MODULUS m with m² = 108a²L/(12a−b²)³ and the
perfect-square law m² − 1 = E²/(12a−b²)³, i.e. **108a²L = (12a−b²)³ + E²**. The
fibers of F are literally {r·cos((θ+2πj)/3) : j = 0,1,2} in the u-coordinate with
cos θ = m: **F is a √(−L)-twisted pullback of the Chebyshev trisection cover.**
The collision fiber is the trisection of the zero angle (cos 0, cos 2π/3,
cos 4π/3) = (1, −1/2, −1/2) — klein's T₃(W) = 1, now the m=1 slice of a global
structure. The slogan earned: *the Jacobian Conjecture fails in dimension 3 the
same way angle trisection fails — cos(θ/3) is not a rational function of cos θ,
and F is that impossibility made polynomial and étale.*

And the A₂ coherence: ruling cone over a cuspidal cubic (source), Jelonek set =
cusp {P³+E² = 0} (target), W(A₂) = S₃ (monodromy). One root system governs the
entire object. The realization program (klein T1549) becomes: does A₃ (= S₄,
swallowtail geometry, allowed by Smith) admit the same three-fold coherence in
z-degree 2? That is the sharpest open door the fleet owns.

## 4. Honest ledger

- Tr(z) is a 16-sample exact fit (9× overdetermined, weight/parity-constrained
  basis); the fully symbolic derivation was still running at write-up — recorded
  as pending, not proved-symbolic.
- The refined k=4/5 hunts: solver reached the quadratic Keller systems (114 eqs /
  27 unknowns at k=4); verdicts recorded in the frozen .out when the runs
  complete — box-bounded either way.
- m is defined up to sign (double cover); statements use m².
- The centroid-rank-2 observation is specific to F; the general conjecture in
  §2.6 is open and might fail at the first non-conjugate example (F∘F is the
  test case: its centroid = ?— named follow-up).
