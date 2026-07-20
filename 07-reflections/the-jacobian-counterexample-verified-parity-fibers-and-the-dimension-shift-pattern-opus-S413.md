# The Jacobian counterexample, verified in-repo: parity fibers, the Dixmier cascade, and the dimension-shift pattern

**Instance:** opus-2026-07-19-S413 (owner: the Jacobian conjecture was recently refuted;
investigate connections to our body of work and the Dixmier conjecture; "tournaments and
their recursion is reminiscent of the n vs 2n relationship"). **HYP-8065, new tangent
thread** (no prior JC/Dixmier/Weyl content in canon — grep-checked). Script + frozen out:
`jacobian_counterexample_verification_opus_S413.py`.

## 1. The verification (self-contained; no literature trust needed)

F(x,y,z) = ((1+xy)³z + y²(1+xy)(4+3xy), y + 3x(1+xy)²z + 3xy²(4+3xy), 2x − 3x²y − x³z).

- **det J(F) = −2 identically** (symbolic, sympy expand of the 3×3 determinant).
- **F(0,0,−1/4) = F(1,−3/2,13/2) = F(−1,3/2,13/2) = (−1/4, 0, 0)** — three distinct
  points, exact rational arithmetic. A polynomial self-map of ℂ³ with constant nonzero
  Jacobian that is not injective: **JC₃ is false.** The fact is self-verifying; whatever
  the provenance, the map settles the matter.

## 2. The structure we found in it (new, proved here)

- **ℤ/2-equivariance:** under σ(x,y,z) = (−x,−y,z), F₁ is even and F₂, F₃ odd (symbolic).
  So σ permutes every fiber over the fixed target axis {(c,0,0)}.
- **The fixed line is rigid:** F(0,0,z) = (z,0,0) — a bijection of the axis onto itself.
- **PARITY LEMMA (three lines):** every finite fiber over the axis has cardinality
  **1 + 2k (odd)** — one preimage on the fixed line, the rest in σ-pairs. Non-injectivity
  therefore requires fiber ≥ 3, and the discovery realizes the minimum: its collision
  fiber is {fixed point} ∪ {one mirror pair}. Verified at a second fiber: F⁻¹(1,0,0) =
  {(0,0,1)} ∪ {(±i/2, ±3i, −26)} — again 1 + 2. The counterexample appears to be a
  **degree-3 étale self-cover** (two independent fibers of size exactly 3), the smallest
  size the parity permits. This is the repo's mirror-involution parity lore (S407-D; the
  Rédei involution-pairing mechanism; "#components odd ⟺ fixed-point membership")
  appearing verbatim in affine algebraic geometry — the one connection in this session
  that is a theorem, not an analogy.

## 3. The Dixmier cascade (written precisely)

With ¬JC₃ in hand: (i) **stability** (F × id is again Keller and non-injective) gives
¬JC_n for every n ≥ 3; (ii) van den Essen's arrow D_n ⟹ JC_n gives, contrapositively,
**¬D_n for every n ≥ 3 — the Dixmier conjecture fails for all Weyl algebras A_n, n ≥ 3**;
(iii) the Tsuchimoto/Belov-Kanel–Kontsevich arrow JC_{2n} ⟹ D_n loses its forward use
but leaves the surviving frontier sharp: **JC₂, D₂, and Dixmier's original D₁ remain
open**, connected by the last bridge standing, JC₂ ⟹ D₁. The n-vs-2n web, post-refutation,
collapses onto dimensions one and two — the same "small dimensions are where rigidity
survives" shape as Jung–van der Kulk tameness in dim 2 vs Nagata wildness in dim 3; the
counterexample lives exactly at the tame/wild boundary.

## 4. The dimension-shift pattern (the owner's directive, graded RHYME with one hook)

The owner's instinct: JC_{2n} ⟹ D_n rhymes with tournament recursion. Graded honestly:

- The BKK mechanism behind the ×2 shift is passage to the **char-p center of A_n**, a
  symplectic polynomial ring in 2n variables — a doubling functor from the quantum object
  to a commutative shadow. The repo's own ×2 functor is the **SC blowup / voltage lift**
  (THM-378: double-round-robin doublings classified by triangle parities). RHYME: both
  attach a ℤ/2-layered commutative double to a noncommutative/oriented object, and in
  both the classification data is a PARITY (triangle parities there; the σ-fiber parity
  of §2 here). One precise hook worth a future session: whether THM-378's parity
  classification has an Azumaya-style reading.
- The repo's own dimension-SHIFT reductions: LRC(14) ⟺ the n−2-speed inverse statement
  (THM-1017); Mode B staircase reduction removes both legs (n → n−2); the Cayley–Dickson
  tower's 2n+... ladder. The shared pattern — "the conjecture at size n is governed by a
  rigidity statement at a shifted size, through a functor that adds/removes a layer" —
  is real as a pattern; no theorem is claimed across domains, and per the S405 grading
  discipline it is recorded as a rhyme with named hooks, not an identity.

## 5. What this changes for the repo

Nothing in the LRC program logically depends on JC/Dixmier — the audit found zero
canon dependencies. What the event contributes: (a) a live demonstration that
sixty-year conjectures fall to EXPLICIT SMALL WITNESSES verifiable in exact arithmetic
— the repo's entire referee discipline vindicated at field scale (the counterexample is
to JC what a slack-1 family would be to LRC(14): one certificate, checkable by anyone);
(b) the parity lemma as a genuinely portable tool — equivariant maps with rigid fixed
loci have odd fibers over fixed values, minimal exotic fiber 3 — the same 1 + 2k
signature as our witness-component parity; (c) a caution transferred from our own
MISTAKES ledger to the wider world: the JC survived every average-case/degree-reduction
attack and fell to a structured construction — the "sampler-blind-to-extremal" genus at
historic scale.

## Cross-links

S407-D + LEM-020 (the mirror/parity lore this realizes) · THM-378 (the repo's ×2 parity
functor) · THM-1017 (the repo's own n↔n−2 bridge) · S399 triage (why explicit witnesses
beat average-case attacks) · TANGENTS entry (this thread) · script + frozen out (S413).
