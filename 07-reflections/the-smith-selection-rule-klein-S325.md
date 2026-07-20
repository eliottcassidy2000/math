# The Smith selection rule: the topology above the fleet's Jacobian anatomy

**Instance:** klein-2026-07-19-S325 (owner: "explore these wild conjectures and pull
from other agents to formulate and explore more").
**Computation:** `jacobian_smith_selection_rule_klein_S325.py` + frozen out (degrees
2–7 tabulated exactly). **Topology input:** P.A. Smith's theorem (classical, 1930s–40s:
a ℤ/p action on a finite-dimensional mod-p acyclic space has nonempty fixed set) —
citation-solid; everything else below is elementary group theory + the fleet's
verified anatomy.

## 1. The rule

ℂⁿ is contractible; the deck group of an étale cover acts freely. Smith ⟹ no
finite-order homeomorphism of ℂⁿ is fixed-point-free. Hence, for ANY étale
polynomial self-cover F : ℂⁿ → ℂⁿ of generic degree d:

> **(R1) F has NO nontrivial deck transformation. Galois covers of degree ≥ 2 do
> not exist: generic degree 2 is IMPOSSIBLE, over every ℂⁿ.**
> **(R2) The monodromy group G ≤ S_d must be transitive with |Fix(G₀)| = 1 —
> the point stabilizer fixes only its own point (deck group N_G(G₀)/G₀ = 1).**

## 2. The table (computed, degrees 2–7)

| d | DEAD | ALLOWED |
|---|---|---|
| 2 | C₂ | — (impossible) |
| 3 | C₃ | **S₃** (unique!) |
| 4 | C₄, V₄, D₄ | A₄, S₄ |
| 5 | C₅ | D₅, F₂₀, A₅, S₅ |
| 6 | C₆, S₃-reg, D₆ | PSL(2,5), PGL(2,5), A₆, S₆ |
| 7 | C₇ | D₇, F₂₁, F₄₂, PSL(3,2), A₇, S₇ |

Emergent laws inside the table: **dihedral D_d allowed ⟺ d odd** (odd-gon
reflections fix exactly one vertex); **every Frobenius group in its standard
action is allowed** (stabilizers fixing only their point IS the Frobenius
property); **every 2-transitive action is allowed**; **regular (= Galois)
actions never**.

## 3. What this does to the fleet's converged anatomy

The four independent S₃ results of the last 24 hours — kps THM-1310 (resolvent
√(−L)), mac-mini THM-1315 (syzygy + van der Waerden), boxeph-S142 ("S₃
universal" in the lifting class), klein-S323 (mod-p Chebotarev census) — are
**instances of a single topological necessity**: at degree 3, S₃ is the only
group Smith allows. Likewise the k = 3 rigidity seen three ways (klein's k=2
ansatz death, kps's singular-linear-part for k ≠ 3, death-star's equivariant
"sporadic in chart") sits above the fact that **degree 2 was never available**:
the minimal possible counterexample degree is 3, and F occupies the unique
allowed monodromy slot there. The counterexample is not just verified — it is
the **first inhabitant of the minimal cell of a classification that topology
fixes in advance.**

## 4. The wild conjectures, updated by exploration

- **H1 (was: all minimal counterexamples share S₃/{3,1}/−4S²D)** → the S₃ part
  is now a THEOREM-modulo-citation (R2 at d = 3). The {3,1} fiber spectrum and
  the −4S²D discriminant shape remain conjectural for OTHER degree-3 examples —
  but boxeph's uniqueness ("THE counterexample of its class") and death-star's
  chart-sporadicity suggest there may be essentially no others at degree 3.
- **H2 (was: {2,3}-smoothness)** → TRANSFORMED by death-star's THM-1305 decode:
  the rationals are values of k-polynomial laws (Φ = 4t+6 = (k+1)t + 2k at
  k = 3, Ψ = (1+t)(2+t), 4 = k+1, 2 = k−1). The smoothness was the k = 3 shadow
  of degree-polynomial structure, not a universal prime law. H2 retired in its
  original form; its successor: **the constants of any degree-d example are
  {p ≤ d}-smooth** (still wild, still open).
- **H5 (the next-degree existence)** → sharpened into the **REALIZATION
  PROGRAM**: which ALLOWED monodromies are realized by étale self-covers of ℂ³?
  Known: S₃ at d = 3 (F). Open and now well-posed: A₄ or S₄ at d = 4 (kps's
  linear-part obstruction kills F's chart at k ≠ 3, so a d = 4 example needs a
  genuinely different design); D₅ at d = 5 would be spectacular (a dihedral
  étale quintic cover). Competing wild answers: (W-abundance) every allowed
  group is realized; (W-rigidity) only 3^m-degree compositions of F-conjugates
  exist — the composition monoid's monodromies (F∘F at 9:1, subgroups of
  S₃ ≀ S₃ obeying R2) are the immediate test bed, measurable by the S323
  Chebotarev harness applied to F∘F. Next session's first computation.
- **NEW (H6, from the table):** the Jelonek set {D = 0} (kps THM-1310 proved it
  exactly) is a QUARTIC for the minimal example. Conjecture: the Jelonek degree
  is a monodromy invariant — minimal allowed configurations have minimal
  asymptotic-variety degree.

## 5. Cross-links

THM-1305 (equivariant anatomy; the k-polynomial decode) · THM-1310 (fiber
geometry; Jelonek = {L = 0}; the (Q, L) = klein's (S, D) up to convention) ·
THM-1315 (surjectivity — klein-S324's conjecture PROVED — étale, 3:1, caustic)
· boxeph-S142 (lifting uniqueness; S₃ universal; −1/4 = c/(2φ₁)) · opus-S414
(torus normal form) · klein T1547/T1548 + S323/S324 reflections (verification,
master quartic, minimality mechanism — now upgraded by Smith) · the
structural-zero epistemology (N₂ ≡ 0 = the {3,1} spectrum = R2's shadow in the
census).
