# The Galois group of the apex prime: LRC(14) in ℚ(ζ₇), and exactly where its C₆-symmetry breaks

*kind-pasteur-2026-06-28-S257. The owner: find even more creative reframes, test them, and be inspired by
the concurrent work. This builds directly on mac-mini's hidden-C₃ (S83) and fuses it with my census split
(HYP-3258) and ℚ(√−7) floor (HYP-3254). The reframe: the entire LRC(14) is the arithmetic of the apex prime
7, encoded in ℚ(ζ₇) with Galois group C₆ — and the proof's two halves are the two subfields, with the hard
core located precisely where the C₆-symmetry breaks.*

## The whole runner set is two Galois orbits plus a fixed point

mac-mini found the witnesses are one C₃-orbit. The full picture, under the *complete* Galois group
`(ℤ/14)* ≅ C₆` acting on `ℤ/14` by multiplication, is cleaner still — **the 13 runners are exactly two
regular C₆-orbits plus the fixed apex prime:**

```
  {1,3,5,9,11,13}  — the units (ℤ/14)*        — ONE C₆-orbit — the BINDING skeleton
  {2,4,6,8,10,12}  — the evens                 — ONE C₆-orbit — the COVERING runners
  {7}              — the apex prime            — FIXED (7u ≡ 7) — covering, self-paired
```

`6 + 6 + 1 = 13`. Two copies of the regular representation of C₆, plus the trivial. The apex prime that
*names* the field is the unique fixed runner — it sits in the covering layer, not the skeleton.

## The two halves of the proof are the two subfields of ℚ(ζ₇)

`C₆ = C₂ × C₃`, where `C₂` is complex conjugation = the complement involution `a ↦ −a`, and `C₃` is
mac-mini's pair-cycle `×3 : {1,13}→{3,11}→{5,9}`. By the Galois correspondence, ℚ(ζ₇) has exactly two
intermediate subfields, and they carry the two halves of the proof:

| subfield | degree | fixed by | Galois grp | carries |
|---|---|---|---|---|
| `ℚ(cos 2π/7)` | 3, **real** | `C₂` (conj) | `C₃` | mac-mini's **cap = C₃-trace**, the **equioscillation** at the 3 pairs |
| `ℚ(√−7)` | 2, **imaginary** | `C₃` | `C₂` | my **floor** / the **Gauss sum** `i√7` / the `±` complement |

The cap is rational because it is a `C₃`-trace from the real cubic field; the floor's resonances split
by the `C₂`-arithmetic of the imaginary quadratic field. **The equioscillation and the decorrelation
floor are not two unrelated techniques — they are the two Galois-dual subfields of the apex prime.**

## The sharp dichotomy: the residue half is C₆-symmetric, the magnitude half breaks it

This is the load-bearing test (verified). Apply `×u` for any unit `u`:

- **Residue / equioscillation half — C₆-INVARIANT.** Multiplying the AP by any unit permutes `{1,…,13}`
  (the AP residues *are* the full C₆-invariant set), so `M = 1/14` is preserved at all six units. The
  binding-pair structure is C₆-equivariant: THM-568 (`14 ∣ D` at the optimum — PROVED) forces a binding
  pair at *each* unit optimum, and mac-mini's C₃ organizes the three into one orbit. **This half is
  rigorous.**
- **Magnitude / census half — C₆-BROKEN.** Among the evens orbit `{2,4,6,8,10,12}`, *only* 12 admits a
  tight magnitude-lift (`12→24` = GW); its orbit-mates `2,4,6,8,10` have **none**. The integer magnitudes,
  not the residues, select 12 uniquely. **This half is where the symmetry breaks — and that is the hard
  core.**

So the Galois symmetry does real work: it *closes* the residue half and *localizes* the entire remaining
difficulty to the symmetry-breaking in the magnitude layer — the unique doubling site, and the floor.

## The symmetry-breaking is q-specific (cross-N, verified)

The doubling site is not universal. For the AP at `n = 2q`, the only nontrivial tight single-replacement:

```
  n= 6 (q=3): 2 → 9      (exotic, not a doubling)
  n= 8 (q=4): 6 → 12     = (n−2) → 2(n−2)   ✓ doubling
  n=10 (q=5): NONE       — census = {AP} alone, fully rigid
  n=12 (q=6): NONE       — census = {AP} alone, fully rigid
  n=14 (q=7): 12 → 24    = (n−2) → 2(n−2)   ✓ doubling  (= GW)
```

So GW-type symmetry-breaking exists only for special `q` (here `q = 4, 7`), absent for `q = 5, 6`. The
census size is itself `q`-dependent: a single rigid point for `q = 5, 6`; two points for `q = 4, 7`. The
"order parameter" of the breaking — which even runner can double — is on or off depending on the apex
prime's arithmetic (the Jacobsthal / Goddyn–Wong gcd-window). This is why the magnitude half resists a
uniform argument: it is genuinely `q`-sensitive, while the residue half is `q`-uniform.

## What this buys toward rigor

A clean partition of the remaining work along the Galois decomposition:

- **Residue half (C₆-symmetric, `ℚ(cos 2π/7)`):** equioscillation at the six units, binding pairs forced
  by THM-568 at each unit optimum, organized by mac-mini's C₃. Rigorous.
- **Magnitude half (C₆-broken, the hard core):** (a) the census magnitude rigidity — among integer
  realizations, only AP and the unique `q`-gated doubling GW; (b) the decorrelation floor `R' ≥ c` in the
  covering branch. Both are symmetry-broken and `q`-specific; this is where the proof is not yet rigorous.

The reframe says: stop looking for a single uniform mechanism. The proof is **two-phase** — a
Galois-symmetric residue skeleton (closed) and a symmetry-broken magnitude residual (open, `q`-gated). The
residual is small and explicit (one doubling site, one floor constant), and it lives in the magnitude
layer that no congruence or Galois argument can see.

## The pointer beyond

If the symmetry-breaking is an order parameter gated by the apex prime's arithmetic, the natural question
is *which primes `q` switch it on*. The data (`q = 4, 7` on; `5, 6` off) is the seed of a criterion — the
GW doubling `(n−2) → 2(n−2)` survives iff `2(n−2)`'s danger arcs realign with the cover, a condition on
`q` mod small moduli (the gcd-window `[2,3]`). A clean such criterion would turn the census from a
case-by-case rigidity into a statement about the apex prime, completing the Galois picture: the residue
skeleton is `q`-uniform, and the magnitude break is governed by a single arithmetic switch on `q`.

— Related: [[lrc14-thread]], HYP-3411 (this reframe), [[the-hidden-c3-the-witness-space-is-one-galois-orbit-leverage-the-proved-binding-pair-plus-c3]]
(mac-mini's C₃), HYP-3258 (census split), HYP-3254 (ℚ(√−7) floor), HYP-3256, THM-568 (`14∣D`), THM-523;
`the-census-splits-rigid-unit-skeleton-and-one-doubling.md`, `lonely-runner-as-chebyshev-equioscillation.md`.
