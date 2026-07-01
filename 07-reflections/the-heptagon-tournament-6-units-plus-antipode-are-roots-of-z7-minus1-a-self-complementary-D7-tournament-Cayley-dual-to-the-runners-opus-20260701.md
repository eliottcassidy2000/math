# The heptagon tournament: the 6 unit atoms + the antipode 7/14 ARE the roots of z⁷=−1 (a regular heptagon with dihedral D₇, order 14 = the LRC's 14); the natural "beat-the-next-3" tournament on them is SELF-COMPLEMENTARY (Aut=C₇, ι=converse-iso — it lives on the repo's SC spine), has exactly 14 cyclic triangles (=|D₇|), sits at tiling Hamming-distance 6 (=#units) from the transitive, and CAYLEY-TRANSFORMS EXACTLY (U⁷=I, U∈SO(7)) back to the 7th roots of unity — so the ±1 tournament and the runners-on-a-loop are Cayley-dual. And at the AP's lonely time t=1/14 the 13 runners FACTOR by parity (the CRT 14=2·7): ODD runners → the z⁷=−1 heptagon (the SC tournament, Verblunsky (0,…,0,1)), EVEN runners → z⁷=+1 minus the origin → the HARMONIC Verblunsky 1/(6−j) (a nested LRC(7))

*opus-2026-07-01-S14. Owner: extend the runners-on-a-loop / Verblunsky synthesis to tournaments and dihedral
groups — for the 6 unit atoms build a tournament on 7 vertices (vertex 0↔1/14, vertex 6↔13/14, added vertex at
7/14), apply the tiling model, make the finish concrete, formulate and refine invariants. This bridges the LRC
world to the repo's central tournament-tiling world through the self-complementary / dihedral structure.*

## The object: the heptagon H₁₄ (all verified)
The 7 vertices `v_k = (2k+1)/14`, k=0..6 (vertex 0=1/14 … vertex 6=13/14, added vertex 3=7/14=1/2) are **exactly
the roots of z⁷=−1** (`z_k^7 = e^{πi(2k+1)} = −1`), a regular heptagon. Its symmetry is the **dihedral D₇ of order
14 — the LRC modulus itself**:
- **rotation** ρ: k→k+1 (add 1/7 = 2/14) — the C₇ part;
- **reflection** ι: x→−x, i.e. k→6−k — the antipode. It **fixes vertex 3 (=7/14=1/2)** — *that is why the added
  vertex is special: it is the reflection center* — and pairs the 6 units into 3 ι-pairs {0,6},{1,5},{2,4}.
- The 6 units are exactly the non-fixed vertices: (2k+1) coprime to 14 ⇔ 2k+1≠7 ⇔ k≠3. Units = {1,3,5,9,11,13}.

## The tournament and its invariants (all verified, rigorous)
Order the vertices cyclically and take the natural geometric tournament **R₇ = "beat the next 3"**: i→j iff
(j−i) mod 7 ∈ {1,2,3}. Refined invariants:
- **SELF-COMPLEMENTARY (SC).** |Aut(R₇)| = 7 (the rotations C₇); there are also **7 converse-isomorphisms**
  (T≅Tᵒᵖ), the witness being ι:k→6−k. So the reflections of D₇ act as **arc-reversal** (anti-automorphisms):
  **D₇ acts on R₇ with rotations preserving and reflections reversing** — R₇ is self-complementary and **lives on
  the repo's SC SPINE** (the complement ℤ₂ = ι). |Aut|+|converse-iso| = 7+7 = 14 = |D₇|.
- **Triangle census:** score-regular (all out-degrees 3); **14 cyclic triangles = |D₇| = the LRC 14**; 21 transitive
  triples = C(7,2) = #edges; **175 Hamiltonian paths** (odd, Rédei ✓).
- **Skew-spectrum:** the ±1 skew-adjacency S=A−Aᵀ (circulant) has eigenvalues `2i·Σ_{k=1}^3 sin(2πmk/7)`,
  |λ|_max=4.381. (Its arithmetic sibling, the **Paley/QR tournament** {1,2,4}, has |Aut|=21 and the **√7 Gauss-sum**
  spectrum ±i√7 — the doubly-regular one; connects to the repo's QR_p engineering thread.)
- **CAYLEY DUALITY (exact).** U=(I−S)(I+S)⁻¹ ∈ SO(7), and **U⁷=I** with spectrum **exactly the 7th roots of unity**
  (angles k/7). So the Cayley transform maps R₇'s skew-adjacency to an order-7 rotation whose eigenvalues ARE the
  7th roots — the **±1 tournament and the 7 runners on the loop are Cayley-dual**. The uniform measure on that
  spectrum has Verblunsky (0,…,0,1) [Haar/symmetric].
- **TILING model.** With base Hamiltonian path 0→1→…→6 (the cyclic order) and 15 tiles = the chords (i,j), j≥i+2,
  R₇ = the transitive tournament with **exactly 6 tiles flipped** (Hamming distance 6 = #units), the flipped tiles
  being the **long chords {4,5,6}** (the "back half"). Paley sits at distance 7 (chords {3,5,6}).

## Making the finish concrete: the CRT/parity factorization (verified)
At the AP's lonely time t=1/14 the 13 runners `v/14` **factor by parity — the CRT `ℤ/14 ≅ ℤ/2 × ℤ/7`**:
| half | runners | land on | Verblunsky | meaning |
|---|---|---|---|---|
| **ODD** | {1,3,5,7,9,11,13} | roots of **z⁷=−1** (the full heptagon) | (0,…,0,1) [Haar] | the **SC D₇ tournament**, Cayley-dual to the 7th roots |
| **EVEN** | {2,4,6,8,10,12} | roots of **z⁷=+1** minus origin | **1/(6−j)** exact | the **harmonic law at n=7** = a nested **LRC(7)** |
- The two halves are the two SIGNS of `z⁷=±1`, and **the ± is the parity is ι is the complement ℤ₂**. The `7` is
  the heptagon/D₇; the `2` is the antipode/complement.
- **M=1/14 is set by runners 1 & 13 = vertices 0 & 6 = the extreme ι-pair** (the nearest odd/14 points to the
  origin). Vertex distances to the origin: {1/14,1/14, 3/14,3/14, 5/14,5/14, 1/2} — the ι-pairs at equal radii, the
  gap at the outer pair, the antipode at 1/2.
- So **LRC(14)'s lonely configuration is exactly (parity ℤ/2) × (heptagon D₇)** — a concrete multiplicative/CRT
  descent, with a self-complementary tournament in one factor and a nested harmonic-Verblunsky LRC(7) in the
  other. This ties to the covering-min's arithmetic (mac-mini-S77 / klein-S68's Φ₆ Eisenstein-CRT; kind-pasteur-S7's
  (ℤ/N)\* Ramanujan skeleton — the same units, here carrying a tournament).

## Refined invariant list (the deliverable)
A signature of the extremal, spanning both worlds:
1. **Vertex/symmetry:** V₇ = roots of z⁷=−1 = units(ℤ/14) ∪ {antipode}; symmetry D₇ (order 14), ι = complement.
2. **Tournament type:** self-complementary; |Aut|=C₇; 14 cyclic triangles; 175 Ham paths.
3. **Cayley invariant:** Cayley(R₇)=SO(7) order-7 element, spectrum = 7th roots (U⁷=I); Verblunsky (0,…,0,1).
4. **Tiling invariant:** Hamming distance 6 to transitive; flipped tiles = long chords {4,5,6}.
5. **CRT factorization:** parity split — odd↦SC-tournament heptagon, even↦harmonic-Verblunsky nested LRC(7).

## Honest scope — what's proved vs the open crux
- **Proved/verified (all above):** the geometry (z⁷=−1), the D₇/SC structure, the triangle counts, the exact Cayley
  duality (U⁷=I), the tiling distance, the CRT/parity Verblunsky factorization. These are rigorous.
- **The finish itself (open):** this **characterizes** the AP's lonely configuration as the maximally-symmetric
  (SC / D₇) object, but does not yet **prove** the covering-min lower bound. The natural conjecture it sharpens is
  a **SYMMETRY-EXTREMALITY principle**: *among covering sets, minimal M is achieved by the configuration of maximal
  symmetry* (the SC/D₇ heptagon) — a tournament-side restatement of "AP/GW are extremal." Proving symmetry ⇒ minimal
  gap is exactly the open crux (the for-all-covering-S band-dodging, OPEN-Q-108). So this is a new **invariant/
  characterization + a precise conjecture**, not a closure.
- **Value:** it makes the "finish" concrete as a *structural factorization* and connects the LRC endgame to the
  repo's central objects (SC spine, tiling model, Cayley-Dickson/CRT descent, QR_p) — giving the team a
  tournament-theoretic handle (SC classification, Rédei parity, the metagraph) on the covering-min extremal.

## Status
- **VERIFIED (opus):** 6 units + antipode = roots of z⁷=−1 (D₇, order 14); R₇ self-complementary (Aut C₇ + ι),
  14 cyclic triangles, Cayley→7th roots exact (U⁷=I∈SO(7)), tiling distance 6; t=1/14 CRT parity split =
  odd(SC tournament, Haar Verblunsky) × even(harmonic Verblunsky, nested LRC(7)).
- **CONJECTURE (opus, open):** symmetry-extremality — minimal covering M ↔ maximal (SC/D₇) symmetry; = the
  covering-min crux, not closed.
- **Bridges:** LRC ↔ tournament-tiling world via SC/dihedral + Cayley duality; the units carry a tournament
  (extends kind-pasteur-S7's (ℤ/N)\* skeleton and my HYP-3795 runner-cloud OPUC).

Related: HYP-3795/opus (runner-cloud harmonic Verblunsky — the even half here IS its n=7 case), HYP-3793/kp-S7
((ℤ/N)\* Ramanujan skeleton — the same units), HYP-3792/mac-mini-S77 + HYP-3800/klein-S68 (Φ₆ Eisenstein-CRT — the
covering-min arithmetic), the repo SC spine + merged metagraph G_n/ℤ₂ (ι = complement), the tiling model + Rédei
(01-canon), OPEN-Q-108 (the for-all-covering-S crux this sharpens). HYP-3802 (this). Scripts:
04-computation/lrc_heptagon_{tournament_setup,cayley_verblunsky,aut_sc,crt_parity_decomp}_opus_20260701.py (+.out).
