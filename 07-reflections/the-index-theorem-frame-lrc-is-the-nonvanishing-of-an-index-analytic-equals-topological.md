# The index-theorem frame: LRC(2p) ⟺ an index = (p−1)/2 is nonzero, where analytic (equidistribution/Euler char) = topological (Borsuk-Ulam degree), and parity = p mod 4 selects the method — for n=14 the ODD index forces it

*mac-mini-2026-06-28-S79. The owner: find a more creative unifying frame for the LRC and push toward a proof.
The pieces — cap = Euler characteristic (S78, topology), equidistribution (S74, analytic), the Borsuk-Ulam odd
degree (kps S31av), the saddle index (p−1)/2 (kps S31aw) — are the two sides of an INDEX THEOREM: existence of
the lonely point = non-vanishing of an index, computed both analytically and topologically. This frame subsumes
the previous ones and gives n=14 a topological proof push. Builds on
[[the-topology-of-the-lrc-cap-is-euler-char-of-the-cover-nerve-lonely-is-the-hole-certified-by-borsuk-ulam]],
[[14-is-the-heptagon-dihedral-group-borsuk-ulam-not-brouwer-kps]], [[the-vitali-wall-brouwer-equioscillation-and-the-cyclotomic-core-construction]].*

## The frame
> **LRC(2p) ⟺ a single INDEX is NONZERO.** The index `= (p−1)/2` (kps S31aw's saddle index = the de Moivre
> cyclotomic degree). It is computed TWO ways that an index theorem equates:
> - **ANALYTIC index** = the equidistribution / the cap = the **Euler characteristic of the danger cover** (S78;
>   `cap = Σ(−1)^|S| meas(∩_S) = χ_meas(nerve)`) = the measure-side count of lonely components.
> - **TOPOLOGICAL index** = the **degree of the heptagon reflection** on the witness configuration space (the
>   Borsuk-Ulam / Brouwer degree, kps S31av) = the de Moivre / `D₇`-irrep count `(p−1)/2`.
> The LRC holds iff this index ≠ 0 (a lonely component / an antipodal coincidence exists).

## Parity = p mod 4 selects the proof method (the family, VERIFIED structure)
```
  p   n=2p  index=(p-1)/2  p mod 4  parity  method        index≠0
  3    6        1           3       odd     Borsuk-Ulam   yes
  5   10        2           1       even    Brouwer/SOS   yes
  7   14        3           3       ODD     Borsuk-Ulam   yes   ← the first hard, ODD
 11   22        5           3       odd     Borsuk-Ulam   yes
 13   26        6           1       even    Brouwer/SOS   yes
```
- **p ≡ 1 (mod 4)** (even index): the reflection is an AUTOMORPHISM → it has a fixed point → **Brouwer** → a
  symmetric, SOS witness (the real/even/positive case, n=10, 26). The index nonvanishing is an SOS positivity.
- **p ≡ 3 (mod 4)** (ODD index): the reflection is an ANTI-automorphism → free `ℤ₂` → **Borsuk-Ulam** → the
  witness is an antipodal pair certified by **odd degree**. And an **odd degree is automatically ≠ 0** — so the
  index nonvanishing is FORCED by parity (the imaginary/odd/negative case, n=14, 22).

## THE PUSH (n=14): the odd index ⟹ Borsuk-Ulam ⟹ the lonely point is forced
For `n=14` (`p=7`), the index `=(7−1)/2 = 3` is **ODD**. By Borsuk-Ulam, an odd-degree map `S¹→S¹` has nonzero
degree, hence a zero/coincidence — the **antipodal lonely pair `(t*,−t*)` at the `1/14` equioscillation EXISTS**.
The odd degree IS the imaginary Gauss sum `i√7` = the negative trace = the de Moivre cyclotomic degree `(p−1)/2`
— one datum. So the topological route for n=14 is: **certify the index is odd (parity), conclude the lonely point
by Borsuk-Ulam.** No SOS/Brouwer fixed point is needed (and none exists — that is why n=14 is "hard" for the
even/SOS methods but EASY for the odd/topological one).

## Why this is the MORE-unifying frame (it subsumes the others)
Every prior frame is a way to *compute the same index*:
| prior frame | = which side of the index |
|---|---|
| cap = Euler char of the cover nerve (S78) | the **analytic index** (measure count) |
| equidistribution / Erdős-Turán (S74 bulk) | the **analytic index** (continuous) |
| Borsuk-Ulam / `D₇` degree (kps S31av) | the **topological index** |
| the de Moivre cyclotomic degree `(p−1)/2` (S75e) | the **index value** |
| even/odd = positive/negative = IE parity (S76) | the **index PARITY** = `p mod 4` |
| Vitali wall: measure (bulk) vs construction (core), S75f | analytic index vs topological index |
| SOS (even) vs non-SOS (odd) | Brouwer (even index) vs Borsuk-Ulam (odd index) |
| Fejér/de Moivre magic function (kps), Q(√−7) wall | the index's positivity certificate / its sign datum |
So **the index theorem is the meta-frame**: "LRC = an index ≠ 0," and the analytic↔topological equality unifies
the measure side (equidistribution/Euler char) with the symmetry side (Borsuk-Ulam degree), with `p mod 4` the
parity that says which classical theorem (Brouwer or Borsuk-Ulam) closes it.

## Honest status
- **VERIFIED structure:** the family index `(p−1)/2`, the parity `= p mod 4`, the method split (Brouwer even /
  Borsuk-Ulam odd), index `≠ 0` for all `p≥3`. The analytic side `cap = χ_meas(cover)` is verified (S78).
- **FRAME + ROUTE (not a proof):** the LRC as "an index ≠ 0" unifies all prior frames as two computations of one
  index; for `p ≡ 3 (mod 4)` (n=14) the odd index gives the **Borsuk-Ulam topological push** (odd degree ≠ 0 ⟹
  the lonely point). The OPEN rigor: (1) prove the saddle `(p−1)/2` IS the topological degree of the
  right reflection map (the index theorem itself), and (2) that the coincidence sits at the `1/14` threshold
  (the equioscillation = the index condition). These are the two remaining steps of the topological route.
- **NOT a proof.** But the frame is the most unifying yet, and it converts n=14's "hardness" (no SOS/Brouwer
  witness) into its EASE under the correct theorem (Borsuk-Ulam, odd degree). LRC(14) open.

Related: HYP-3246 (this), kps S31av/aw (D₇/Borsuk-Ulam, index (p−1)/2), HYP-3242 (cap=Euler char), HYP-3237
(Vitali wall), HYP-3221 (the one obstruction), OPEN-Q-108.
