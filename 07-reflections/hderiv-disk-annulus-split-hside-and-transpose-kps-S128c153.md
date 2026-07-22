# The disk/annulus (Wiener–Hopf) structure of the Weierstrass split — and why it trivializes the h-side

*kind-pasteur-2026-07-22-S128c153. Owner: work creatively on the frame factorization + h-side lemma;
pull often; prioritize mathematical reasoning and exploration above builds. These are the two last open
inputs to `hderiv` (the assembly, (c) degree lemma, my F=D_m leg, and boxeph's bridge are all done,
kernel-pure).*

## The one idea: `[x⁰]` is a ring homomorphism on the disk, singular on the annulus

The frame lives in `PowerSeries(LaurentSeries F) = F⸨x⸩⟦t⟧`, and `xCoeff0` reads the `x⁰` Hahn-coefficient
of each `t`-coefficient. The crucial fact, easy to miss:

> **`xCoeff0` is NOT multiplicative on Laurent series, but IS multiplicative on power series.**
> For `f,g ∈ F⟦x⟧` (disk, `x`-support `≥ 0`), `[x⁰](f·g) = [x⁰]f · [x⁰]g` (the constant term of a product
> of power series is the product of constant terms). For Laurent `f,g` this is false: `[x⁰](x·x⁻¹) = 1`
> but `[x⁰]x · [x⁰]x⁻¹ = 0`.

So on the **disk subring** `F⟦x⟧⟦t⟧ ↪ F⸨x⸩⟦t⟧`, `xCoeff0 = PowerSeries.map (constantCoeff)` is a **ring
homomorphism**, and it commutes with the `t`-derivative `∂_t` (they act on different variables).

### This makes the h-side a triviality — but ONLY because `h` is a disk series

The h-side lemma is `xCoeff0(h_t/h) = g_t/g` with `g := xCoeff0(h)`. **It is false for a general Laurent
unit** — verified computationally (`sympy`):

| `h` | `xCoeff0(h_t/h)` | `g_t/g` | equal? |
|---|---|---|---|
| `1 + t(x + x⁻¹)` (Laurent) | `−2t − 6t³ − …` | `0` | **NO** |
| `1 + t(x + x²)` (disk) | `0` | `0` | yes |
| `1 + t·x + t²·x³` (disk) | `0` | `0` | yes |

It holds **exactly because the Weierstrass unit `h ∈ (F⟦t⟧)⟦X⟧` is a genuine power series in `X`**
(non-negative `x`-support). On the disk, with `φ := map(constantCoeff)` a ring hom commuting with `∂_t`:

```
xCoeff0(logDeriv h) = φ(∂_t h · h⁻¹) = φ(∂_t h)·φ(h⁻¹)      [φ ring hom]
                    = ∂_t(φ h)·(φ h)⁻¹                       [φ commutes with ∂_t; preserves unit inverse]
                    = ∂_t(g)·g⁻¹ = g_t/g.                    [g = φ h = xCoeff0 h]
```

**The h-side is `logDeriv` commuting with a ring hom — the whole content is that `h` is on the disk.**

## The disk/annulus split is Weierstrass = Wiener–Hopf = Birkhoff

This is the structural payoff. `Φ = P·h` splits the Laurent world into two regimes, and `[x⁰]` (`= ∮ dx/x`,
the Cauchy/residue functional) sees each differently:

- **`P` = the annulus/"minus" part** (the small-root factor, monic of `x`-degree `M`). `P⁻¹` has *negative*
  `x`-powers; `[x⁰](P_t/P) = 0` is a **degree/pole count** — proved by the reversed-polynomial trick
  `P_t/P = y·Q(y)/P*(y)`, `P*(0)=1` (mac-mini's (c)). Here `[x⁰]` is *not* multiplicative; it is a residue.
- **`h` = the disk/"plus" part** (the unit, `h ≡ 1 mod t`, genuine `X`-power-series). `[x⁰]` is the
  *evaluation at `x=0`* — multiplicative, holomorphic. The h-side (a) is automatic.

So the two "hard" `hderiv` inputs are the **two halves of one Birkhoff factorization**, and each is hard/easy
for the *same* reason `[x⁰]` behaves the way it does on that half. My F=D_m leg is the third functional
`[x⁰](−R/Φ)` — the full Laurent inverse (geometric series), the "moment/trace" `∑ D_m tᵐ = log-det`.
The three legs are `[x⁰]` applied to `logDeriv` of the three objects `Φ` (trace), `P` (annulus/degree),
`h` (disk/value): `trace = degree + value`, i.e. `−(F−1)/t·t = 0 + t·g_t/g` collapses to `g_t = 0` under `F=1`.

## The transpose (frame factorization) — the exact target

`hderiv_of_frame` needs `PhiFrame = Pfr·hfr` with `hfr` a unit. The map is the **re-indexing isomorphism**
`T : (F⟦t⟧)⟦X⟧ ≅ MvPowerSeries (Fin 2) F ≅ (F⟦X⟧)⟦t⟧`, then `map(HahnSeries.ofPowerSeries)` into `F⸨x⸩⟦t⟧`.

Two cautions that shape the formalization:
1. **It is NOT an `eval₂`/substitution `X ↦ x`.** `x = single 1 1` is a *unit* in `F⸨x⸩`, not topologically
   nilpotent, so the power-series universal property does not apply. It must be the coefficient re-indexing
   `∑ a_{n,j} tⁿ Xʲ ↦ ∑ a_{n,j} xʲ tⁿ` (valid in `MvPowerSeries`, symmetric in the two variables).
2. **`T` lands `h` in the disk subring.** `h ∈ (F⟦t⟧)⟦X⟧` has `X`-support `≥ 0`, so `T(h) = map(ofPowerSeries)(H)`
   with `H := (swap) h ∈ F⟦X⟧⟦t⟧ = PowerSeries(PowerSeries F)`. **This is precisely the hypothesis the h-side
   needs.** So the transpose and the h-side compose: the transpose produces `hfr = map(ofPowerSeries) H`, and
   the h-side machinery (below) then gives (a) for free. `T(P)` is a Laurent *polynomial* factor (degree `M`).

So the division is clean: **transpose (death-star) delivers `hfr = map(ofPowerSeries) H`; the h-side is then
a corollary of "`map(ofPowerSeries)` and `xCoeff0` restrict to the disk ring hom `map(constantCoeff)`, which
commutes with `logDeriv`."** I am formalizing that corollary (`GMC2DvdKFrameHSide`), self-contained on the
disk subring, so it plugs into `hderiv_of_frame`'s `ha` the moment the transpose lands.

## Formal reduction of the h-side (what I am building)

`xCoeff0(logDeriv (map ofPowerSeries H)) = derivativeFun(g) * Ring.inverse g`, `g := xCoeff0(map ofPowerSeries H)`,
from three reusable facts:
- **(A1) `derivativeFun_map`**: `map ψ (derivativeFun f) = derivativeFun (map ψ f)` (any ring hom `ψ`).
- **(A2)** `map ψ (Ring.inverse u) = Ring.inverse (map ψ u)` for `u` a unit (ring homs preserve unit inverses)
  ⇒ `map ψ (logDeriv u) = logDeriv (map ψ u)` — **`logDeriv` commutes with ring homs**.
- **(B)** `xCoeff0 (map ofPowerSeries H) = map constantCoeff H` (`[x⁰]∘ofPowerSeries = constantCoeff`).

Then `xCoeff0(logDeriv(map ofPS H)) = map cc (logDeriv H) = logDeriv(map cc H) = logDeriv g = g_t/g`, applying
(A2) at `ψ = ofPowerSeries` and `ψ = constantCoeff` and (B) in the middle. All disk-side, no annulus, no Puiseux.

## Cross-links
`GMC2DvdKFrame` (death-star frame) · `GMC2DvdKHderiv.hderiv_of_frame` (the assembly consuming `ha`, `hfact`) ·
`GMC2DvdKFrameExtraction` (my F=D_m leg) · `GMC2DvdKFrameDegree` (mac-mini (c)) · THM-1550 · HYP-9016.
