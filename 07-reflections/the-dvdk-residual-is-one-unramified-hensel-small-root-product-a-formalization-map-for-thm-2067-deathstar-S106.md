# The DvdK residual = one unramified-Hensel small-root product: a formalization map for THM-2067

**death-star-2026-07-22-S106** (HYP-8935). Owner: finish the GMC(2) formalization by working the *remaining
reasoning* (not Lean builds); find creative ways to bypass DvdK or replace it with something easy to formalize;
mine past threads. **This is not a completed formalization.** It is (a) a precise dependency map showing the
GMC(2) DvdK residual is already closed on paper by codex THM-2067, whose formalization decomposes into four
Mathlib-ready pieces plus **one** genuinely valued-field gap; (b) an elementary re-derivation of the single
analytic input (THM-1550) with no residues/monodromy; and (c) a creative simplification of the gap — the one hard
object, the small-root product `Π(t)`, is computable by an **unramified Hensel lift**, replacing the ramified
Puiseux/Newton-polygon construction (the "≈4–9 person-months" item of the S95 roadmap) with standard Hensel
(already in Mathlib) plus one explicit degree-`M` base change. All four claims verified numerically
(`dvdk_residual_formalization_map_deathstar_S106.py`, 4 faces incl. both hard cancellation examples).

## 0. Where the dependency actually sits (mining S100/S101/S223/S226–S229 + THM-2067/2070)

GMC(2)/THM-2022 needs exactly one non-Mathlib fact (`GMC2DvdKInterface.DvdK1`): for the lowest face
`f_F = Σ_{i∈F} c_i u^{q_i}` (distinct integer charges straddling 0, `c_i ≠ 0`), some power has a nonzero constant
term, `∃ m0 ≥ 1: CT(f_F^{m0}) ≠ 0`. The confinement threads already peel most of this off:
- **Unique minimal balanced cycle** (S101/HYP-8878, boxeph S229): `CT(f^{m*})` is a single nonzero term —
  DvdK-free, coefficient-independent (84% of small supports). **Formalized kernel-pure.** boxeph S230 (HYP-8931,
  same day) then localizes the entire spine dependency to the single seed lemma
  `GMC2FaceSeed.exists_nonzero_lowest_face_seed` and gives a drop-in DvdK-axiom-free replacement for the
  unique-channel class — so DvdK now bites **only** on the coincident-channel stratum.
- **Positive coefficients** (S100, boxeph S228) and **two-charge/edge** (S100, boxeph S226/S227): elementary,
  **formalized kernel-pure.**
- The genuine residual is **≥2 coincident minimal cycles with signed/complex coefficients that cancel** — the
  "resonant" stratum. boxeph S230 pins down *why it is nonempty and irreducible*: the involution `u ↦ −1/u`
  (`f(−1/u) = −f(u)`, THM-2070) pairs balanced compositions, forcing **even** multiplicity at every mass, so
  symmetric faces (e.g. `{−2,−1,1,2}`) are **never** unique-channel; and boxeph blocked every elementary erasure
  route (face-simplification is impossible — THM-2070: *any* Laurent polynomial is some GMC lowest face; the
  saddle route was retracted; char-`p` is worse; genericity conflates feasibility with cancellation). boxeph's
  S223 coprime-interval/numerical-semigroup route is likewise *refuted* (THM-2070's `u^2+u+u^{-1}-u^{-2}` has
  support-return at every `m≥2` yet `CT(f^m)=0` for all odd `m`, first nonzero at `m=4`): support feasibility
  decides *which* balanced words exist, **not** whether the coefficients cancel.

So the residual is a real cancellation problem, and on paper it is **already solved** by **codex THM-2067**
(Galois orbit-product). The remaining task is *formalization*. This note maps that task.

## 1. THM-2067 in one paragraph, and its formalization decomposition

Shift the face to `Λ(u) = u^{-M} R(u)`, `R = r_0 + … + r_d X^d`, `r_0 r_d ≠ 0`, `d = M+N`, `M,N ≥ 1`
(`M` = negative extent). Then `CT(Λ^m) = [u^{Mm}] R(u)^m =: D_m` (Check A). Let `Φ(X) = X^M − t R(X)`, with roots
`Ω` in a splitting field of `ℂ(t)`. THM-2067:
1. **THM-1550 (Wiener–Hopf criterion):** `D_m = 0 ∀m ≥ 1  ⟺  Π(t) := ∏_{small roots} = c·t`, `c = (−1)^{N+d+1}r_0 ≠ 0`.
2. **Irreducibility:** `Φ` is irreducible over `ℂ(t)` — degree 1 in `t` + `gcd(X^M, R) = 1` (from `R(0) ≠ 0`) +
   Gauss. Hence `Gal` is **transitive** on `Ω`.
3. **Vieta:** `C_Φ := ∏_{Ω} α = (−1)^d r_0/r_d ∈ ℂ^*` (valuation 0).
4. **Orbit-product + valuation:** a proper `M`-subset `S_0` (the small roots) with `∏_{S_0} = Π ∈ ℂ(t)` forces,
   by transitivity/double-counting, `Π^r = C_Φ^η`. `t`-adic valuation: `r·val(ct) = η·val(C_Φ)` i.e.
   `r·1 = η·0 = 0`, so `r = 0` — contradiction. Therefore some `D_m ≠ 0`.

**Formalization status of each piece (this is the map):**

| piece | content | Mathlib status |
|---|---|---|
| orbit-product lemma (§1) | finite Galois + double-counting + `p^r = C^η` | **ready** (`IsGalois`, transitivity from `Irreducible`) |
| irreducibility of `X^M−tR(X)` | Gauss + `gcd(X^M,R)=1` | **ready** (`Polynomial.Monic.irreducible`, Gauss) |
| Vieta `∏Ω = (−1)^d r_0/r_d` | product of roots | **ready** (`Polynomial.prod_roots`/coeff) |
| `t`-adic valuation contradiction | order of vanishing at `t=0` on `ℂ(t)` | **ready** (`RatFunc`/`Valuation`) |
| **THM-1550: `∏_{small} = ct`** | Newton-slope factor of `Φ` over `ℂ((t))` | **the one gap** |

Everything except THM-1550 is ordinary finite algebra already in Mathlib. The valued-field content is confined to
one lemma.

## 2. Elementary re-derivation of THM-1550 (no residues, no monodromy) — verified

Because `1 − tΛ(u) = u^{−M}Φ(u)` and `Φ(X) = −t r_d ∏_{α∈Ω}(X−α)`, grouping roots by `|α|` gives the Wiener–Hopf
factorization
```
1 − tΛ(u) = const(t) · ∏_{small u_i}(1 − u_i/u) · ∏_{large a_j}(1 − u/a_j),
const(t) = (−t r_d)(−1)^N ∏_j a_j.
```
The first product has only `≤0` powers of `u`, the second only `≥0`, each with `u`-constant-term 1. Taking `CT_u`
of `log`, both root products die, leaving
```
CT_u log(1 − tΛ) = log const(t).
```
But `CT_u log(1 − tΛ) = −Σ_{m≥1} CT(Λ^m) t^m/m = −Σ_{m≥1} D_m t^m/m`. Hence
```
−Σ_{m≥1} D_m t^m/m = log const(t),   and via Vieta  const(t) = (−1)^{N+d+1} t r_0 / Π(t).
```
So `D_m = 0 ∀m ⟺ const(t) = 1 ⟺ Π(t) = (−1)^{N+d+1} r_0 t = c t`. This is THM-1550, obtained from a formal
`log` of a Laurent-series factorization — **no contour integration, no residue theorem, no analytic continuation.**
The only structure used is the small/large **slope factorization of `Φ` over `ℂ((t))`.**
*Verified (Check B):* `log const(t)` matches `−Σ D_m t^m/m` to machine precision on all four faces; Vieta `C_Φ`
and `Π(t)/(ct) → 1` confirmed (Check D). Corollary: `Π(t) = ct·exp(Σ D_m t^m/m) ∈ t·ℂ[[t]]` — **`Π` is
unramified** (integer powers of `t`), even though the individual small roots have valuation `1/M`.

## 3. The creative simplification: the small-root product is an UNRAMIFIED Hensel lift — verified

The one gap is "compute `Π(t) = ∏_{small roots}` and know it lies in `ℂ((t))`." The small roots have `t`-valuation
`1/M` (ramified), which is why the roadmap invoked Puiseux series / extending the valuation to
`AlgebraicClosure ℂ((t))` (the months-long item). **This is avoidable.** Substitute `X = sZ`, `t = s^M`:
```
Φ(sZ) = s^M Z^M − s^M R(sZ) = s^M · ψ(Z),   ψ(Z) = Z^M − R(sZ) ∈ ℂ[[s]][Z].
```
Modulo `s`, `ψ(Z) ≡ Z^M − r_0`, which is **separable** (`M` distinct `M`-th roots of `r_0 ≠ 0`). So by **ordinary
(unramified) Hensel over the complete local ring `ℂ[[s]]`** — no ramification, no Puiseux — `ψ` has a monic
degree-`M` factor `A(Z) ∈ ℂ[[s]][Z]` with `A ≡ Z^M − r_0 mod s`, and the small roots of `Φ` are `u_i = s·Z_i`
with `Z_i` the roots of `A`. Hence
```
Π(t) = ∏_i u_i = s^M ∏_i Z_i = t · (−1)^M A(0),   A(0) ∈ ℂ[[s]],  A(0) ≡ −r_0 mod s.
```
`D_m = 0 ∀m ⟺ A(0) = −r_0 constant (∈ ℂ) ⟺ Π = ct`. *Verified (Check C):* `Π(t)` from the unramified Hensel
factor matches `Π(t)` from direct root-finding on all four faces.

**What this buys.** The gap changes from
> *"construct the Puiseux/ramified splitting field and extend the `ℂ((t))` valuation to its algebraic closure"*
(the S95 ≈months item)
to
> *"Hensel-lift the separable factor `Z^M − r_0` of `Z^M − R(sZ)` over `ℂ[[s]]`, plus the fixed degree-`M`
> base change `t = s^M`."*

Hensel over a complete local ring is **in Mathlib**; `LaurentSeries`/`PowerSeries` and their valuation are in
Mathlib. The residual valued-field work is the **local–global bridge** — matching the `M` Hensel roots in
`ℂ((s))` to a proper `M`-subset of the abstract `ℂ(t)`-splitting-field `Ω` and confirming `∏ = ct ∈ ℂ(t)` — over
the *unramified* `ℂ((s))` (with the explicit `s^M = t`), rather than an abstract ramified closure. That is a
bounded, standard target, not an open research programme.

## 4. Honest scope

- **Not** a Lean proof of `DvdK1`, and **not** a claim that DvdK is fully bypassed: the residual (mixed-sign
  coincident-cycle cancellation) is genuine and is settled on paper by THM-2067, which this note takes as given.
- **Is** three things: (i) a decomposition showing 4 of the 5 THM-2067 ingredients are Mathlib-ready; (ii) an
  elementary, residue-free proof of the 5th's analytic input (THM-1550); (iii) a reduction of the one
  valued-field object from ramified Puiseux to unramified Hensel + a fixed base change. Net: "formalize the DvdK
  residual" moves from an open ≈months valuation-theory project to a mapped task whose only non-elementary piece
  is a standard Hensel factorization plus a local–global bridge.
- The confinement (S100/S101) still does the up-front work: on the 84%+ unique-cycle stratum none of this is
  needed (single-term `CT`), so THM-2067 + this map bite only on the resonant coincident-cycle set.
- No effective first-return bound (`m0` size) is produced or claimed — the Sturmfels/ESV effective problem is
  untouched, exactly as THM-2067 notes.

## 5. Suggested next Lean targets (for the formalizing agents)

1. The **orbit-product lemma** (§1) — self-contained finite Galois; likely a few hundred lines against
   `Mathlib.FieldTheory.Galois`. Discharges the abstract core independent of the gap.
2. **Irreducibility** of `X^M − tR(X)` over `ℂ(t)` via Gauss + `R(0) ≠ 0`.
3. The **unramified Hensel factor** `A(Z)` of `Z^M − R(sZ)` over `ℂ[[s]]` and `Π = t·(−1)^M A(0)`; then the
   local–global bridge to a proper `M`-subset of `Ω` with rational product `ct`.
   Piece 3 is the only one touching valued fields, and (§3) shows it is unramified-Hensel-shaped.

Cross-links: THM-2067 (Galois orbit-product, the paper closer — this note maps/simplifies its formalization),
THM-1550 (Wiener–Hopf criterion, re-derived elementarily here), THM-2070 (dihedral cancellation example),
THM-2022 (GMC2), S100/HYP-8877 (edge/positive confinement), S101/HYP-8878 (unique-cycle, the DvdK-free 84%),
boxeph S223/S226–S230 (coprime-interval refutation, formalized special cases, and S230's localization of the
dependency to one seed lemma + the involution characterization of the residual — this note maps the formalization
of exactly the THM-2067 residual S230 leaves open, and dovetails with S230's proposal to parameterize the descent
by the seed lemma), S95 (DvdK roadmap: the Puiseux/Monsky ≈months item this simplifies),
memory `nc2-gmc2-lean-formalization-state`,
`gmc-lrc-same-positivity-manoeuvre`. Script `04-computation/dvdk_residual_formalization_map_deathstar_S106.py`
(+ `05-knowledge/results/…S106.out`). HYP-8935.
