# Two geometries, one algebra: the Cayley transform bridges the staircase and the circle — complement is a reflection in both, and the transitive spectrum is the roots of unity

*opus-2026-07-01-S22. The owner asked to find more algebra↔geometry links like R = complement / V = span-swap,
notice the pattern, and extend. The pattern generalizes sharply: every tournament algebra fact has geometric
incarnations in **two** spaces — the staircase tiling and the unit circle — and the Cayley transform is the
bridge. This ties the repo's tournament/tiling pillar to its LRC/runners pillar.*

## The pattern: one algebra, two geometries
The finding "R (anti-diagonal reflection) = complement" was a *geometric reflection = algebraic involution*.
The generalization: tournament algebra lives simultaneously in the **staircase** (the tiling model) and on the
**unit circle** (via the Cayley transform of the skew-adjacency matrix), and the same algebraic operation is a
reflection in **both**.

**The Cayley bridge.** A tournament `T` gives a skew-symmetric sign matrix `S = A − Aᵀ` (`S_ij = ±1`). Since
`S` is real skew, `I+S` is always invertible, and the **Cayley transform** `U = (I+S)⁻¹(I−S) ∈ SO(n)` is
orthogonal, so its eigenvalues lie **on the unit circle**. Those eigenvalue-points are a geometric fingerprint
of `T` — and they are exactly the "runners on a loop" of the LRC side.

## The dictionary (verified n≤5, spectra n≤9)
| Algebra (tournament) | Staircase geometry | Circle geometry (Cayley spectrum) |
|---|---|---|
| **complement** `Tᵒᵖ` (`S → −S`) | anti-diagonal reflection **R** (`x+y ↦ 2n+2−(x+y)`) | **complex conjugation** (angle → −angle; reflect across the real axis) |
| **self-complementary** (SC) | **grid-symmetric** tiling (R-fixed) | **conjugation-symmetric** spectrum (mirror-image about the real axis) |
| **transitive** (empty tiling) | the all-forward corner | **roots of unity** — n odd: the `n`-th roots `e^{2πik/n}`; n even: the primitive `2n`-th roots `e^{iπ(2k+1)/n}` |
| **rotational / circulant** | difference-striped tiling | exact roots of unity (the runner cloud) |
| **span** `d = x−y` | main-diagonal grade (span 1 = cut space / Ham path; span ≥2 = cycle space / tiles) | (eigenvalue spacing — open) |
| **3-cycles** `c₃` | grid triangles with cyclic orientation (`c₃ = C(n,3) − Σ C(sᵢ,2)`) | spectral moment (open) |

The first three rows are the striking ones:

- **Complement is a reflection in both geometries.** On the staircase it is R (reflection across the
  anti-diagonal `x+y=n+1`); on the circle it is complex conjugation (reflection across the real axis). *One
  involution, two mirrors.* (Verified `S → −S` ⇒ spectrum angles negate, n=3,4,5.)
- **SC is the fixed set in both.** Grid-symmetric tilings (R-fixed) ⟺ self-complementary ⟺ conjugation-symmetric
  spectra (real-axis-mirror-invariant point sets). (Verified n=3,4,5.)
- **The transitive tournament's spectrum is the roots of unity** — for odd `n` the `n`-th roots (with `1` at
  angle 0), for even `n` the primitive `2n`-th roots (no `1`). This is the same **odd/even parity split** that
  runs through the whole project, now as "is `1` in the spectrum." (Verified exactly n=3..9.)

## Why this is the bridge to the LRC pillar
The Cayley eigenvalues are points on the unit circle — **the LRC runners**. The transitive tournament maps to
the `n`-th roots of unity, which is exactly the AP runner configuration `{e^{2πik/n}}` at the lonely time
`t=1/n` (my heptagon/HYP-3802 result: the 7 units → roots of `z⁷=±1`). So:
- The **staircase** (tiling model, tournament combinatorics) and the **circle** (LRC runners, OPUC/Verblunsky)
  are two coordinate systems for the *same* objects, glued by the Cayley transform `S ↔ U`.
- **Complement** is the shared involution: anti-diagonal reflection on one side, conjugation on the other.
- The project's mantra "**Everything is the Triangle**" gains a dual, "**Everything is the Circle**," and the
  Cayley map is the dictionary between them.

## The odd/even parity, unified
The `n` odd/even split shows up in *both* geometries as the *same* phenomenon:
- **Circle:** odd `n` ⇒ the spectrum contains a real-axis point (`1`, angle 0); even `n` ⇒ it does not.
- **Algebra (S20/HYP-3811):** an SC tournament's anti-automorphism has exactly one **fixed vertex** iff `n` is
  odd. A conjugation-symmetric set of odd cardinality *must* have an odd number of real-axis points ⇒ ≥1 —
  and that real eigenvalue **is** the anti-automorphism's fixed vertex. (Theoretical corollary; the odd-vs-even
  fixed-point rule of HYP-3811 and the "is 1 in the spectrum" rule are the same fact seen in two geometries.)
So the `≡ 2 mod 4` cycle rule (HYP-3811, algebra), the roots-of-unity parity (circle), and the transpose-self
tiling parity (staircase) are three views of one thing.

## Extensions (open, proposed)
1. **Span on the circle.** Span `d=x−y` grades cut/cycle on the staircase; find its circle image (does flipping
   a span-`d` tile rotate/perturb specific eigenvalues?).
2. **`c₃` and the spectral moments.** `c₃ = C(n,3) − ΣC(sᵢ,2)` (grid triangles); express via `Tr` of powers of
   `U` or `S` — a third face (algebraic invariant ↔ grid count ↔ spectral moment).
3. **V (span-swap) on the circle.** V is staircase-geometric-only (S21); its circle image (if any) would say
   whether "swap span↔anti-diagonal" is spectrally meaningful.
4. **The order-4 rotation** `ρ = R∘V` (a 90° turn of the tile grid, `ρ⁴=id`) — is there an order-4 algebraic
   operation (a "√complement") matching the `i` of the Cayley/Cayley-Dickson tower?

## Status
- **Verified:** complement = conjugation of the Cayley spectrum (n=3,4,5); SC = conjugation-symmetric spectrum
  (n=3,4,5); transitive spectrum = roots of unity, odd `n` → `n`-th, even `n` → primitive `2n`-th (n=3..9);
  circulant → exact roots.
- **Pattern:** one algebra, two geometries (staircase + circle); complement is a reflection in both; the Cayley
  transform is the bridge; the odd/even parity is "is 1 in the spectrum" = the anti-automorphism fixed vertex.
- **Bridges** the tournament/tiling pillar (R=complement, grid-sym=SC, S18–S21) to the LRC pillar (runners =
  Cayley eigenvalues = roots of unity, HYP-3802/3795).
- **Honest:** the three headline rows are verified; the span/c₃/V/ρ circle-images are proposed, not computed.

Related: S18 (R = complement, staircase reflection), S21 (V = span-swap), HYP-3811 (SC = cycles ≡2 mod4 + one
fixed vertex iff odd n — = the real-axis eigenvalue), HYP-3802/3795 (heptagon: units → roots of unity, Cayley
duality, Verblunsky), "everything is the triangle" (the staircase side). HYP-3813 (this). Scripts:
04-computation/cayley_{alggeo_dictionary, transitive_roots_of_unity}_opus_20260701.py.
