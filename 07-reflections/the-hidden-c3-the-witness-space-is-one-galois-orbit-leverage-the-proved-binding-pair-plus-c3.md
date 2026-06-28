# The hidden C₃: the LRC(14) witness space is a single C₃-Galois orbit — the creative leap is to leverage the PROVED single binding pair (HYP-2909) + C₃ to generate the whole equioscillation

*mac-mini-2026-06-28-S83. The owner: see more of the underlying hidden structure and make a creative leap on
direction. The hidden structure: every "3" in the LRC(14) is ONE C₃, and the three witnesses are a single
C₃-Galois orbit. The leap: the proof should be C₃-equivariant, generating the full equioscillation from the
already-PROVED single binding pair. Builds on
[[reframing-the-last-two-details-resonant-reduces-to-off-grid-bulk-rigidity-is-a-finite-equioscillation-system]],
[[the-cap-is-a-totally-real-cyclotomic-quantity-and-n14s-two-cyclotomic-heads]], kps S255 (equioscillation/units), HYP-2909.*

## The hidden structure (CONFIRMED): every "3" is ONE C₃, and the witnesses are a single orbit
The three binding-pair witnesses `{1,13}, {3,11}, {5,9}` form a **single C₃-Galois orbit** — `×3 mod 14` cycles
them: `{1,13} → {3,11} → {5,9} → {1,13}` (verified). And this C₃ is the same object seen six ways:
```
  C₃  =  (ℤ/14)* / {±1}   (the 3 antipodal unit-pairs)
      =  (ℤ/7)*  / {±1}    (the 3 antipodal residue-pairs mod 7, via CRT 14=2·7)
      =  Gal( ℚ(cos 2π/7) / ℚ )   (the de Moivre cubic's Galois group)
      =  {2cos(2πj/7) : j=1,2,3}   (the 3 de Moivre angles, one Galois orbit)
      =  the index (p−1)/2 = 3    (kps S79/S256)
      =  #QR mod 7 = 3            (kps S256, the Gauss sum)
```
So **the LRC(14) is fundamentally a C₃-statement.** The witness space (the 6 units = 3 pairs) is a single
C₃-orbit; the index is `|C₃|`; the cap lives in the fixed field of C₃ (ℚ(cos 2π/7), S75e).

## The creative leap (the DIRECTION): a C₃-equivariant proof from the proved single binding pair
The whole remaining proof should be run **C₃-equivariantly on the witness space**, not config-by-config:
1. **One witness is already PROVED.** HYP-2909 (Lean-formalized, sorry-free): `M=1/14 ⟹ a binding pair
   `14 | (s_i + s_j)` at the optimum `t*` — i.e. ONE of the 3 pairs `{a, 14−a}` is forced.
2. **C₃ generates the other two.** The witness space is the single C₃-orbit `{1,13}→{3,11}→{5,9}`; the proved
   one pair, transported by `×3`, yields all three — **the full equioscillation at the 6 units (kps S255) is
   generated from the proved single binding pair + the C₃ symmetry.** No new binding-pair proof needed; the
   symmetry does the work.
3. **The cap is the C₃-trace** (hence rational — S75e's `disc 7²`, the Galois-invariance): `cap = Tr_{ℚ(cos2π/7)/ℚ}`
   of a de Moivre element, summed over the C₃-orbit. The rationality IS the C₃-invariance.
4. **The two last details become C₃-equivariant** (S82):
   - rigidity (a) = the C₃-equivariant finite system (the 3-orbit equioscillation has only AP/GW integer solutions);
   - equidistribution (b) = the C₃-symmetric off-grid bulk (the danger is C₃-symmetric on the unit grid, the bulk
     survives by the symmetry-averaged equidistribution).

> **The leap:** stop attacking LRC(14) on the 13-dimensional config space. Run it on the **C₃ witness space**:
> the proved binding pair (HYP-2909) is ONE point of the C₃-orbit; C₃ generates the equioscillation; the cap is
> the C₃-trace; the rigidity and equidistribution are the C₃-equivariant residuals. The proof's natural home is
> `ℚ(cos 2π/7)` with its Galois group C₃ acting on everything.

## Why this is a real direction-change
Previous work treated the witnesses (6 of them), the index (3), the de Moivre cubic (degree 3), the QR mod 7
(3 of them), and the binding pairs (3) as separate facts. They are ONE C₃. So:
- the index-theorem (S79) `index = (p−1)/2` is `|C₃|`;
- the equioscillation (kps S255) at 6 units is the C₃-orbit of 3 pairs;
- the cyclotomic cap (S75e) lives in the C₃-fixed field;
- the proved binding pair (HYP-2909) is ONE C₃-point.
Unifying them means the remaining proof has **one symmetry group (C₃) acting on one object (ℚ(cos 2π/7))**, and
the proved input (HYP-2909) is a single orbit-point — the rest is "spread it by C₃ + the totally-positive trace."

## Honest status
- **CONFIRMED:** the 3 binding-pair witnesses are a single C₃-Galois orbit (×3-cycled); C₃ = (ℤ/14)*/{±1} =
  Gal(ℚ(cos2π/7)) = de Moivre angles = index = #QR mod 7.
- **LEAP (direction, not a proof):** run the proof C₃-equivariantly — HYP-2909 (one pair, proved) + C₃ generates
  the equioscillation; the cap is the C₃-trace; the last two details are C₃-equivariant residuals. This is the
  natural home of the proof (ℚ(cos 2π/7), Galois C₃).
- **Residual (now C₃-framed):** (a) the C₃-equivariant finite system has only AP/GW (the rigidity); (b) the
  C₃-symmetric off-grid bulk has M ≥ 1/14 (the equidistribution). Both inherit the C₃ symmetry. NOT a proof;
  LRC(14) open — but the remaining work is now one symmetry (C₃) on one field (ℚ(cos 2π/7)).

Related: HYP-3257 (this), HYP-3255 (the reframe), HYP-2909 (the proved binding pair), S75e (cap in ℚ(cos2π/7)),
kps S255 (equioscillation/units), S256 (index=Gauss sum), OPEN-Q-108.
