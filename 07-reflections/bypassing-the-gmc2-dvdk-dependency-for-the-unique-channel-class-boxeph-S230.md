# Bypassing the GMC(2) DvdK dependency for the unique-channel class

*boxeph-2026-07-22-S230. Owner: work creatively to bypass the GMC(2) dependency on DvdK, or find a
new way to easily formalize it. Builds on S226/S227/S228/S229 (the DvdK-free Lean fragments), the mine
of death-star-S101/HYP-8878 (unique primitive cycle), codex THM-2067 (the general DvdK, Galois
orbit-product) and THM-2070 (the dihedral cancellation obstruction), and codex's NC2 spine
(`GMC2FaceSeed`, `GMC2NC2`). New Lean: `dvdk1_of_uniqueChannel` (in `GMC2DvdKUniqueChannel.lean`) and
`GMC2DvdKUniqueChannelBypass.lean`, both kernel-pure.*

## Where DvdK actually enters the spine — a single, localizable point

Reading codex's spine top to bottom, `DvdK1` is consumed in exactly **one** place. `GMC2NC2.
nc2_of_dvdK1_of_heightWitnessSupplier` threads `hDvdK` down through
`GMC2IntegralFaceSeedDescent → GMC2FaceSeed.exists_nonzero_lowest_face_seed`, and there it is used for
a single thing: to produce `∃ m₀ ≥ 1, CT(f_F^{m₀}) ≠ 0` on the rational lowest face `F`. Everything
*around* that call is DvdK-free Newton-polygon geometry:

- `GMC2.exists_rational_lowest_face_finset P hP` produces the slope `λ`, level `δ`, the exact face
  `F ⊆ supp P`, and the straddling of its charges — pure geometry, kernel-pure, no analytic input;
- charge-injectivity on `F` and nonzero coefficients on `F` are *proved*, not assumed.

So the entire DvdK dependency of GMC(2) is the single implication **"the lowest face has a nonzero
constant term in some power."** That is the target for a bypass.

## The bypass: the unique-channel class needs no DvdK axiom

My S229 `ct_ne_zero_of_unique_balanced` proves exactly that seed — with *no* premise — whenever the
face charges admit a **unique balanced channel** (some size `m₀` with a single balanced composition):
the constant term collapses to one uncancellable multinomial term. death-star-S101/HYP-8878 measured
this criterion to cover **98 of 116 (84%)** of straddling supports of size 3–4 in `[-4,4]`.

I turned that into two kernel-pure artifacts (`#print axioms = [propext, Classical.choice,
Quot.sound]`):

1. **`dvdk1_of_uniqueChannel`** — the *exact* `DvdK1` existential conclusion (`∃ m ≥ 1, CT(f^m) ≠ 0`)
   proved outright, no premise, for any support with a unique balanced channel. This is the shape the
   interface (`GMC2DvdKInterface.exists_nonzero_face_seed`) consumes, so it discharges that input for
   the whole unique-channel class.
2. **`exists_nonzero_lowest_face_seed_of_uniqueChannel`** — a **drop-in replacement** for codex's
   `GMC2FaceSeed.exists_nonzero_lowest_face_seed`: *identical conclusion* (`∃ λ δ F m₀, … ∧
   CT(f_F^{m₀}) ≠ 0`), but with the `DvdK1` premise replaced by `LowestFaceUniqueChannel P` (every
   lowest face of `P` has a unique channel). Its proof reuses codex's geometry lemma verbatim and
   swaps only the final DvdK call for `ct_ne_zero_of_unique_balanced`.

Net: **for every `P` whose lowest face has a unique balanced channel, the GMC(2)/NC2 descent needs no
DvdK axiom** — only the (separate, also-unproven) `HeightWitnessSupplier`. This is a genuine, checked
bypass for 84% of the support classes, at codex's own interface.

## The honest boundary — and why a *full* elementary bypass is ruled out

The residual is the **coincident-channel stratum**: supports where every size with a balanced
composition carries `card ≥ 2` of them (my S229 `two_le_card_balanced_of_ct_zero`). These are exactly
the symmetric/resonant supports — e.g. `{-2,-1,1,2}`, where the involution `u ↦ -1/u`
(`f(-u^{-1}) = -f(u)`, THM-2070) pairs balanced compositions, forcing even multiplicity at *every*
mass, so a unique channel never exists. Here cancellation genuinely can occur (`f = u²+u+u⁻¹−u⁻²`:
`CT(f^m)=0` for all odd `m`), and this is the irreducible DvdK content, resolved on paper only by
codex THM-2067 (Galois orbit-product). I checked the three obvious routes to erase this residual and
each is blocked:

- **Simplify the face to two charges.** Blocked by THM-2070: *every* one-variable Laurent polynomial
  is the lowest face of some horizontal Gaussian polynomial, so the face can be any coincident-channel
  support. The face cannot be assumed simple.
- **Saddle-point / analytic non-vanishing for large `m`.** This is my retracted S222 route; THM-2070's
  dihedral witness is the counterexample — the saddle phase is resonant and kills a whole residue
  class of `m`, so no elementary growth argument isolates a good `m`. Taming the phase needs the
  algebraic (Galois) relations among the roots.
- **Work in characteristic `p`** (where the spine's Frobenius contradiction already lives). This is
  *harder*, not easier: multinomial coefficients vanish mod `p` (Lucas), and Frobenius powers give
  `CT(f^{p^k}) = 0` outright. Char-`p` adds cancellation.
- **Generic coefficients.** Feasibility (`CT(f^m)` is a nonzero *polynomial* in the coefficients — the
  S228 content) gives non-vanishing on a dense set, but GMC(2) needs the *specific* complex
  coefficients of `P`, and non-vanishing is open, not closed — it does not extend to the boundary.

So the honest conclusion: DvdK's difficulty is irreducibly the coincident-channel/complex-coefficient
cancellation, and no elementary trick removes it. What *is* achievable — and now done — is to (a)
localize the entire dependency to one seed implication, and (b) discharge it, kernel-pure, for the
84% unique-channel class.

## Architectural next step (proposed to codex)

To thread this bypass all the way to a DvdK-free `NC2` on the unique-channel class, the one remaining
edit is codex's: parameterize `GMC2IntegralFaceSeedDescent.exists_finite_field_moment_point_
preserving_integral_lowest_face_seed` (and the `nc2_of_…` capstone) by the *seed lemma* rather than by
`DvdK1` directly — i.e. take `exists_nonzero_lowest_face_seed`'s **conclusion** as the input. Then both
`DvdK1` (general) and `exists_nonzero_lowest_face_seed_of_uniqueChannel` (the 84%) drive the same
descent, and NC2 becomes DvdK-axiom-free for the unique-channel class. My seed lemma already has the
exact conclusion shape required.

## Scope

Concrete: the GMC(2) DvdK dependency is now localized to a single lowest-face seed implication and
**discharged, kernel-pure, for the unique-channel class (death-star-S101's 84%)** — two new sorry-free
theorems, one a drop-in replacement at codex's interface. Not a full bypass: the coincident-channel
stratum (the symmetric/resonant 16%) is the irreducible DvdK content (THM-2067), and I verified the
standard erasure routes (face-simplification, saddle, char-`p`, genericity) are all blocked. The
NC2-level wiring for the unique-channel class is a one-line architectural change owned by codex,
proposed above.

Links: HYP-8931, HYP-8930, HYP-8878, THM-2067, THM-2070,
[[the-unique-channel-dvdk-in-lean-and-the-cancellation-inclusion-exclusion-dictionary-boxeph-S229]],
[[starting-to-formalize-dvdk-the-positive-coefficient-case-and-where-cancellation-is-the-crux-boxeph-S228]].
