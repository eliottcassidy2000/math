---
id: THM-1440
title: "REFLECTION = TORUS — THM-1375's named open question resolves AFFIRMATIVELY. The meridian monodromy of a loop around the Jelonek quartic and the lambda = -1 torus involution are the SAME transposition of the fibre. Mechanism: the fibre cubic is DEPRESSED — P(x) = L*x^3 + (4-3bc)*x - 2c with no x^2 term — so the three fibre x-coordinates always sum to ZERO; on the tau-fixed locus the fibre is therefore {0, r, -r}, the two roots that escape as L -> 0 are exactly the pair sigma = (x -> -x) swaps, and the square-root branch of x ~ +-sqrt(-(4-3bc)/L) transposes them. Also: independent re-derivation of the Jelonek polynomial L = 27a^2c^2 - 18abc + 16a + b^3c - b^2, cross-validating canon (-L|_{c=0} = b^2 - 16a); and TWO sheets are lost over {L=0}, not one, so the fibre drops 3 -> 1"
status: PROVED (three-line argument on the tau-fixed locus) + VERIFIED NUMERICALLY (4000-step continuation, permutation [1,0,2] both ways) + CROSS-VALIDATED against canon's L
author: opus-2026-07-20-S405
depends_on: [THM-1300 (the counterexample), THM-1330 (Jelonek = Zariski three-cuspidal quartic, referee-grade), THM-1350 (sigma/tau fixed-locus reduction), THM-1375 (the lattice; poses this question as NEXT (a))]
---

# THM-1440 — Reflection = torus

## 0. The question, from THM-1375

> **NEXT: (a) do the two reflections coincide? compute the monodromy of a loop around the
> Jelonek quartic and compare with the `λ = -1` action — if they agree the counterexample
> is 'reflection = torus'.**

**Answer: they agree.** Both are the transposition fixing the `σ`-fixed sheet and swapping
the other two. Proof in §3; the mechanism in §4 is more informative than the fact.

## 1. The fibre cubic (independently re-derived)

`F3 = c` is linear in `z` for `x ≠ 0`; substituting `z = (2x − 3x²y − c)/x³` into `F1 = a`,
`F2 = b` and eliminating `y` by resultant leaves a single cubic factor in `x`:

```
P(x)  =  L·x³  +  (4 − 3bc)·x  −  2c
L(a,b,c)  =  27a²c² − 18abc + 16a + b³c − b²
```

**Cross-validation with canon.** THM-1330/1375 record the Jelonek set as `{L = 0}` with
`−L|_{g=0} = b² − 16a`. Here `−L|_{c=0} = b² − 16a` — **exact match**, obtained by an
independent elimination. Sanity checks all pass: `det JF = −2`; `F(σp) = τ(Fp)` identically
with `σ = (−x,−y,z)`, `τ = (a,−b,−c)`; `F(0,0,z) = (z,0,0)`; and at `(a,b,c) = (1,0,0)`,
`P = 16x³ + 4x` with roots `0, ±i/2`, matching canon's `F⁻¹(1,0,0) = {(0,0,1)} ∪
{(±i/2, ±3i, −26)}`.

## 2. Two new structural facts

**(a) The cubic is DEPRESSED — there is no `x²` term, for every target `(a,b,c)`.**
Hence, by Vieta, **the three `x`-coordinates of every fibre sum to zero**:

```
x₁ + x₂ + x₃ ≡ 0        on every fibre of F
```

The fibre's `x`-centroid is identically zero. (This looks like it may *be* the repo's
"polynomial-centroid conjecture" for this map; flagged for whoever owns that thread.)

**(b) TWO sheets are lost over the Jelonek set, not one.** As `L → 0` the degree drops
from 3 to 1, so the fibre drops `3 → 1`. The surviving root is `x = 2c/(4 − 3bc)`.

**(c) `L` is `τ`-invariant** (`L(a,−b,−c) = L(a,b,c)`, verified symbolically), and
`P_τ(x) + P(−x) ≡ 0`. So **`σ` acts on the fibre exactly as `x ↦ −x` on the roots of `P`.**

## 3. The proof

Work on the `a`-axis `{b = c = 0}`, which lies inside the `τ`-fixed locus, so the meridian
loop and `σ` act on the *same* three roots — the comparison is direct, not up to conjugacy.

*Transversality.* `∂L/∂a = 54ac² − 18bc + 16 = 16 ≠ 0` at the origin, so `a = 0` is a
**smooth** point of the quartic and the `a`-axis crosses it transversally. Hence
`a = εe^{iθ}` is a genuine meridian.

*The fibre.* On `b = c = 0`, `P(x) = 16a·x³ + 4x = 4x(4a x² + 1)`, with roots

```
x = 0 ,      x = ± ½(−a)^{−1/2}.
```

*Monodromy.* Transporting `θ : 0 → 2π` sends `(−a)^{−1/2} ↦ e^{−iπ}(−a)^{−1/2} =
−(−a)^{−1/2}`. So the meridian **swaps the two non-zero roots and fixes `x = 0`.**

*σ.* By §2(c), `σ` is `x ↦ −x`, which **swaps the two non-zero roots and fixes `x = 0`.**

Identical permutations. ∎

**Numerical confirmation** (independent of the argument): 4000-step continuation of the
roots of `16a x³ + 4x` around `a = 0.3·e^{iθ}` returns the permutation `[1, 0, 2]`; `σ`
induces `[1, 0, 2]`. The fixed label carries the root `x = 0`.

## 4. The mechanism — why this is not a coincidence

The agreement is *forced by depression*. Because `P` has no `x²` term, the roots sum to
zero. Over the `τ`-fixed locus the `σ`-fixed root is `x = 0` (it is the point `(0,0,a)`,
where `F(0,0,z) = (z,0,0)`), so the remaining two must sum to zero: the fibre is
`{0, r, −r}`. Now:

- **`σ` acts as `x ↦ −x`**, so it must swap `r ↔ −r` and fix `0`;
- **as `L → 0`, the two escaping roots satisfy `Lx³ + (4−3bc)x ≈ 0`, i.e. `x ≈ ±√(−(4−3bc)/L)`**
  — a `±` pair, transposed by the square-root branch cut around `L = 0`.

Both operations are "negate the balanced pair, fix the odd one out", and there is only one
such pair because the centroid vanishes. **The `±` structure that makes `σ` an involution
is the same `±` structure that makes the escaping sheets a conjugate pair.** So
"reflection = torus" is a corollary of the depressed cubic, and the honest statement is
that the torus reflection and the asymptotic monodromy are two readings of one Vieta
relation.

## 5. Scope and what is NOT claimed

- The proof is on the `τ`-fixed locus `{b = c = 0}`. That is the right place to compare (it
  is where both act on a common fibre), but it does **not** compute the full monodromy
  representation `π₁(ℂ³ ∖ {L=0}) → S₃`. Canon records that image as **all of `S₃`**
  (THM-1375); a single transposition generates only `ℤ/2`. **The remaining work is the
  other generators** — loops around the three cusps, where Zariski's non-abelian
  `π₁(ℙ² ∖ C)` of order 12 lives.
- Nothing here weakens THM-1300/1330. This is structure *inside* the known counterexample.

## 6. Open, and now sharper

1. **The cusp monodromy.** Meridians give transpositions; the `S₃` image needs the cusps.
   Zariski's three-cuspidal quartic has `π₁(ℙ²∖C)` metacyclic of order 12 surjecting onto
   `S₃`. Computing the local monodromy at a cusp `(b²/12, b, 4/(3b))` (canon's Gröbner-exact
   cuspidal edge) would complete the representation.
2. **Does depression persist?** Is `Σxᵢ ≡ 0` special to this `F`, or does every degree-3
   Keller counterexample admit a coordinate in which the fibre polynomial is depressed? If
   the latter, "reflection = torus" is a *general* feature of degree-3 counterexamples, not
   a property of this one — that is the version worth chasing.
3. **Two sheets, not one.** The `3 → 1` drop (§2b) is a stronger degeneration than the
   generic `3 → 2`. Whether that is forced by the three-cuspidal (rather than nodal)
   Jelonek type is a concrete question about the rung between nodal and cuspidal that
   THM-1375 wants to insert below Smith.

## Verification

`04-computation/jelonek_monodromy_opus_S405.py` — symbolic elimination (sympy resultant),
the `τ`-invariance and `P_τ(x) = −P(−x)` identities, transversality, and the 4000-step
numerical continuation. Output in `05-knowledge/results/jelonek_monodromy_opus_S405.out`.
