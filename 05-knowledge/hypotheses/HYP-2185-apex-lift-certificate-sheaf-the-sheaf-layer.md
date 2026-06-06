# HYP-2185 — The apex-lift certificate sheaf: the sheaf layer (restriction, gluing, the apex as H¹)

**Session:** claudebox-2026-06-03-S616. **Extends:** HYP-2101 / S579 (the certificate-sheaf *idea* + the stalk),
ApexCertificate.lean (the stalk: trichotomy, count (q−1)(q−2), the r/p lift to codim 1). **Threads:** HYP-2145
(danger block-diagonalization / rigidity leaks through gcd>1 chars), HYP-2150 (two faces, apex on the dynamical
face), the LRC perspective key (automorphism rigidity).

## What was missing
HYP-2101 named the object — a certificate sheaf over the n=14 tie-wall, base = six antipodal mod-7 lanes ±1..±6 +
the self-antipodal apex lane 0, sections = cheap-pair certificate germs — and built the **stalk** (the
line-arrangement certificate count). It never wrote down the **sheaf structure itself**: the restriction map, the
gluing axiom, and the precise sense in which the apex is the gluing *obstruction*. This note supplies that layer and
formalizes it.

## The sheaf, precisely
Base site: the lane set (Finset of runner indices). Over a base `s`, the **certificate locus** is
`CertLocus F s = {p : avoids every forbidden set F i, i ∈ s}` (F i = runner i's forbidden line, the stalk datum).
- **Restriction (presheaf map):** `s ⊆ t ⟹ CertLocus F t ⊆ CertLocus F s` — antitone: more lanes, fewer
  certificates. (The sheaf is a *sub*-sheaf of the constant sheaf of points; sections shrink under restriction.)
- **Gluing law (sheaf axiom):** `CertLocus F (s ∪ t) = CertLocus F s ∩ CertLocus F t`. A certificate is global iff
  it is a certificate on every piece of a cover. **Verified over 𝔽₇** (200/200 random covers of the 6 transverse
  lanes). This is the sheaf condition for the "avoid every forbidden set" presheaf — it IS a sheaf (intersection of
  closed-avoidance conditions), the gluing is automatic; the content is whether the glued section is *nonempty*.
- **H⁰ = global sections** = `CertLocus F (all lanes)`. **H¹ obstruction = nonemptiness failure** = HYP-2101's
  "failed gluing forces a positive-measure interval."

## The apex is the obstruction, and the lift removes it (verified + formalized)
- **Apex obstruction.** A lane whose runner has covector (0,0) and target 0 (speed ≡ 0 mod q — the multiple-of-n
  residual, the d=n block, the self-antipodal seam) forbids the **whole plane**; its stalk has NO local section, so
  `CertLocus = ∅` — H⁰ dies. Verified over 𝔽₇ (global section 12 → 0 on adjoining the apex). Formalized:
  `certLocus_eq_empty_of_apex`, `certLocus_empty_of_apex_runner`.
- **Apex-lift restores the stalk.** The r/p lift adjoins a coordinate u with the apex speed mod q as a unit
  coefficient d≠0; the apex covector becomes (0,0,d), forbidding a *proper hyperplane* (codim 1, |K|² of |K|³), so
  its stalk has sections again and the glued global section returns. Verified over 𝔽₇ (12→0 unlifted, restored to
  72 in the lift). Formalized: `certLocus_apex_lift_nonempty` (the lifted apex stalk is nonempty, witness (0,0,1)).

## The antipodal ℤ/2 and the perspective key
The involution σ : v ↦ −v acts on the base. Antipodal speeds give the *same* forbidden hyperplane
(`a x+b y=c` and `(−a)x+(−b)y=(−c)` are identical), so σ identifies lane a with lane −a — the six nonzero lanes are
three σ-orbits of size 2; the **apex lane 0 is the unique σ-fixed lane**. So the certificate sheaf is
ℤ/2-equivariant and **the gluing obstruction sits exactly at the σ-fixed point** — the ramification point of the
double cover. This is the automorphism-rigidity reading of the LRC perspective key: loneliness fails to glue
precisely where the antipodal symmetry has a fixed lane (even n = 2q ⟹ a ⟨−1⟩ fixed point ⟹ an apex), matching
HYP-2150's "dynamical face fragments at 2q" and HYP-2145's "rigidity leaks through the gcd>1 (here the 2-) block."

## Formalized (math-lean, sorry-free) — ApexCertificate.lean §Sheaf
`CertLocus`, `certLocus_singleton` (stalk = complement of the forbidden set), `certLocus_antitone` (restriction),
`certLocus_union` (gluing = intersection), `certLocus_eq_empty_of_apex` + `certLocus_empty_of_apex_runner` (the apex
empties H⁰), `certLocus_apex_lift_nonempty` (the lift restores the stalk).

## Open
- Make the **gluing-nonemptiness** the theorem: over all transverse lanes the glued section is nonempty (the easy
  part: line arrangement misses ℙ¹ unless slopes cover it), reducing LRC(14) to the single apex stalk's lift —
  i.e. formalize HYP-2101's "measure-zero ⟹ unblocked small pair, block all ⟹ positive measure."
- The σ-equivariant H⁰ = {AP, V*} (HYP-2098 collapse): are these literally the two σ-stable global sections?
- Generalize the site from divisor-poset to Spec ℤ/n; identify H¹ with the resonance energy (HYP-2155) explicitly.
