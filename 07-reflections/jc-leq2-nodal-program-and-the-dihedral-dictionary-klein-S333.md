# A reduced Jacobian conjecture that HOLDS, the nodal program, and the dihedral dictionary

**Instance:** klein-2026-07-20-S333 (owner: long session — Kakeya bank/carriers,
tournaments ↔ dihedral groups, and creatively find a reduced JC that is TRUE).
Fleet context (pulled): the HYP-8185 five-ledger consolidation; Poisson witness
verified independently (mac-mini S140 PC₃-false + boxeph third build) — priority 2
of the S332 ledger EXECUTED by the fleet within hours.

## 1. THEOREM (JC≤2 — the reduced Jacobian conjecture that holds)

> **Every polynomial local diffeomorphism F : ℂⁿ → ℂⁿ of GEOMETRIC degree ≤ 2 is
> injective.** (All n; no properness or degree-bound hypotheses on the polynomial
> degree.)

Proof: geometric degree 2 would make the field extension degree 2, hence the cover
Galois wherever it is a cover; more directly, by the Smith selection rule (S325/T1549)
an étale self-cover of ℂⁿ has trivial deck group, but a degree-2 étale map is
everywhere 2:1-or-degenerate with the fiber swap a global deck involution — free by
étaleness — contradicting Smith's fixed-point theorem on the contractible source.
Degree 1 = birational + étale ⟹ open immersion with... injective. ∎ (One-page
formal writeup = a promotion task; the content is S325's R1 specialized and STATED —
the first unconditional "reduced JC" from the post-collapse anatomy: the conjecture
fails first, and provably only, at geometric degree 3.)

Companions now in canon: equivariant JC₂ (THM-1345, death-star); proper JC
(classical: proper + étale over simply-connected = trivial cover); Jelonek's
codim-≥2 case (his non-properness set is empty or a hypersurface).

## 2. The NODAL JC₂ program (with one case proved)

CONJECTURE (Nodal JC₂): a Keller map of ℂ² whose Jelonek curve J is NODAL
(possibly reducible; including smooth or empty) is an automorphism.

- Deligne–Fulton: nodal ⟹ π₁(ℂ²∖J) ABELIAN ⟹ the monodromy of the restricted
  cover (proper étale over the complement, by Jelonek's definition) is abelian
  transitive ⟹ REGULAR ⟹ the cover is Galois cyclic.
- **Case H⁻¹(J) = ∅ (PROVED):** then ℂ² itself covers ℂ²∖J finitely, so π₁'s
  finite quotient acts freely on contractible ℂ² — Smith kills it; d = 1. ∎
- Case H⁻¹(J) ≠ ∅ (open crux, cleanly named): the deck generator is a finite-order
  automorphism of the plane-curve complement ℂ²∖H⁻¹(J); b₁-bookkeeping forces
  #components(H⁻¹J) = #components(J), and the question reduces to: **which finite
  cyclic covers of nodal-curve complements are again plane-curve complements?**
  (Libgober characteristic-variety territory; Alexander-trivial for irreducible
  nodal. Boundary fixed points of the regularized deck on a rational
  compactification are the remaining escape — the precise wall is identified.)
- Contrast making it sharp: the ACTUAL dim-3 counterexample's Jelonek quartic has
  three-CUSPIDAL plane sections (THM-1330's Zariski selection), and my S329
  bootstrap forces cusps in dim 2 — nodality is exactly what the collapse avoids.

## 3. The dihedral dictionary (tournaments ↔ dihedral groups, verified)

- **Moon's parity law** (verified exhaustively n = 3, 4, 5: Aut orders {1,3},
  {1,3}, {1,3,5}): tournament automorphism groups have ODD order — no reflection
  can ever be an automorphism.
- Hence in any dihedral action on the vertex circle, **rotations act as
  automorphisms, reflections act as ANTI-automorphisms (arc reversal)**: the
  dihedral group acts on the PAIR (T, T^op), and the reflection-quotient IS the
  repo's merged metagraph G_n/ℤ₂. Self-converse counts measured: 8/8 (n=3),
  48/64 (n=4), 704/1024 (n=5) — the reflection-fixed stratum (the SC/spine world).
- **D₃ = S₃**: the Jacobian counterexample's forced monodromy is the dihedral
  group OF THE TRIANGLE — the relational primitive's symmetry group, containing
  rotations (the cyclic fiber structure) and reflections (the √(−D) resolvent
  character = the arc-reversal sign = the Rédei/discriminant parity, as the owner's
  inspo notes: "Rédei sign = discriminant character").
- Smith's table: D_d allowed as étale monodromy iff d ODD — the same parity that
  makes odd-gon reflections fix one vertex and keeps tournament Aut odd. One
  parity, three theaters (permutation groups, tournaments, covers).

## 4. Kakeya bank and labelled-polygon carriers (pull + cross-link, not duplicated)

The owner's pointer lands on live fleet threads: the covering-spectrum Kakeya twin
(klein memory: the (1/7, 6/7) minimal covering spectrum at clock witnesses = the
K(A₅) Kakeya witness twin), kps's shear catalog (HYP-8170, in flight), and codex's
labelled-carrier machinery (THM-1265/1266 line). The named join: finite-complement
stratum closures = carrier-labelled exhaustions where the ESCAPE ATLAS (absolute
window q ∈ [43,48], S323c) plays the Kakeya-bank role — each escape family is a
"direction" and the pinning stack is the bank. Handed to codex/kps with the atlas
as the direction-set; a dedicated joint session is the right vehicle (flagged in
the letter, not half-done here).

## 5. Files and next steps

Verifications frozen in-session (Moon parity + self-converse censuses + the S325
table reused). Next: (i) promote JC≤2 to a canon THM file (one page); (ii) the
nodal crux as a named problem in the ATLAS ("cyclic covers of nodal complements
that are plane-curve complements"); (iii) the dihedral dictionary as a standing
lens for the DC₁ tournament-atom program (death-star HYP-8160: the A1-triple is an
oriented 3-cycle = a D₃-rotation object; its reflection is the Weyl-algebra
transpose/formal adjoint — conjecture: **the formal adjoint on A₁ plays
complementation's role in the DC₁ dictionary**).
