# THM-415 — Signed-LRC: prime `2n−1` ⟹ full sign-orbit, via cyclic homometry and a Galois-stable zero set

**Status:** PROVED (Galois argument; verified exhaustively `n=3..22`, `C=5..43`, and consistent with
all composite data). This **closes the forward direction of HYP-2270** (prime ⟹ orbit `=2^{n−2}`),
which THM-413 left open (THM-413 proved only the `3∣C` ⟹ degenerate direction).
**Source:** monad-explorer-2026-06-06-S705. Builds on THM-413 (order-3 silent flip), THM-401/403
(shells / modulus `C=2n−1`), and S699's signed-LRC theory (T1 gauge-invariance, T2 sign=cut).
Independent of, and complementary to, THM-414 (S704, additive-energy/Krawtchouk face).

---

## The reframing: sign-orbit collisions = CYCLIC HOMOMETRY

Setup as in THM-413: runners `V={1,…,n−1}`, modulus `C=2n−1`, a cut `ε∈{±1}^{n−1}` (up to global
swap, `2^{n−2}` of them), folded pair-clock `ρ(ε_i i − ε_j j)`, `ρ(s)=min(s mod C, C−s mod C)`.

**Key observation.** Map the signed runner to a point `u_i = ε_i·i ∈ ℤ/C`. Because
`{0} ∪ {±1,…,±(n−1)} = ℤ/C` (since `2(n−1)+1 = C`), the point set `S_ε = {ε_i i}` is a
**half-system selection**: an `(n−1)`-subset of `ℤ/C\{0}` taking exactly one of `{i, C−i}` per
magnitude. The folded pair-clock `ρ(u_i − u_j)` is exactly the **circular distance** between the
points. Hence:

> two cuts collide (same folded clock-multiset) **⟺ `S_ε`, `S_{ε'}` are HOMOMETRIC** (same difference
> multiset in `ℤ/C`) **⟺ they have the same Patterson power spectrum `|χ̂|²`.**

This is the classical cyclic-homometry / phase problem of crystallography, specialized to half-system
selections. **Verified** (`signed_lrc_homometry_s705.py`): the three partitions of the `2^{n−2}` cuts
— by folded clock-multiset, by `ℤ/C` difference multiset, and by `|DFT|²` — coincide exactly, all
`n≤13`.

## The closed form (the engine)

Write `ζ = e^{2πi/C}` and `f̂_ε(t) = Σ_i ζ^{t ε_i i}`. Since `cos` is even and `sin` is odd in the
sign,
```
   f̂_ε(t) = Σ_i cos(2π t i/C)  +  i · Σ_i ε_i sin(2π t i/C)  =  A(t) + i·Φ(ε)_t,
```
where `A(t)=Σ_i cos(2π t i/C)` is **cut-independent** and `Φ(ε)_t = Σ_{i=1}^{n−1} ε_i sin(2π t i/C)`
is the **signed sine sum**. Therefore
```
   |f̂_ε(t)|²  =  A(t)²  +  Φ(ε)_t² .
```
**Verified** exhaustively (`signed_lrc_sine_product_s705b.py`, all cuts, `n≤11`; equivalent form
`|f̂|²=|F|²−4σ_Pσ_N`). Immediately:

> **Collision criterion.** `ε, ε'` collide **⟺ `Φ(ε)_t² = Φ(ε')_t²` for every `t=1,…,(C−1)/2`**,
> i.e. `Φ(ε')_t = ±Φ(ε)_t` with an **independent sign at each frequency `t`**.

The map `Φ(ε) = M ε` is linear, `M[t,i] = sin(2π t i/C)`, a `(C−1)/2 × (C−1)/2` **discrete sine
transform**. `M` is invertible for every odd `C` (the half-range sines are orthogonal:
`Σ_t sin(2πti/C)sin(2πtj/C) = (C/4)δ_{ij}`; `|det M|` confirmed nonzero, `C≤29`,
`signed_lrc_collision_anatomy_s705c.py`).

---

## Theorem (prime ⟹ full orbit) — PROVED

> If `C = 2n−1` is **prime**, the only collision is the trivial global swap `ε' = −ε`. Hence the
> `AP_n` folded sign-orbit equals `2^{n−2}` exactly.

**Proof.** Suppose `ε ≠ ±ε'` collide. Put `δ = ε − ε'` and `s = ε + ε'`; for each `i` exactly one of
`δ_i, s_i` is `0` (the other `±2`), so `δ` is supported on the **disagreeing** coordinates
`Diff = {i : ε_i ≠ ε'_i} ≠ ∅` and `s` on the agreeing ones. By the collision criterion and
`Φ(ε)²−Φ(ε')² = Φ(δ)·Φ(s)` (linearity),
```
   Φ(δ)_t · Φ(s)_t = 0      for every t.                                  (★)
```
Let `h(x) = Σ_{i∈Diff} ε_i x^i ∈ ℤ[x]` (a nonzero polynomial; coefficients `±1`). Then
`Φ(δ)_t = 2·Im(h(ζ^t)) = (h(ζ^t) − h(ζ^{−t}))/i`, so
```
   Φ(δ)_t = 0  ⟺  h(ζ^t) = h(ζ^{−t})  ⟺  h(ζ^t) ∈ ℝ .
```
Apply any `σ_a ∈ Gal(ℚ(ζ_C)/ℚ) = (ℤ/C)^*` (`σ_a(ζ)=ζ^a`). Since `h` has integer coefficients,
`σ_a(h(ζ^t) − h(ζ^{−t})) = h(ζ^{at}) − h(ζ^{−at})`. Thus `Φ(δ)_t = 0 ⟹ Φ(δ)_{at} = 0` for all
`a∈(ℤ/C)^*`. The zero set `Z = {t ∈ ℤ/C : Φ(δ)_t = 0}` (extended antisymmetrically, `Φ(δ)_{−t} =
−Φ(δ)_t`, so `Z = −Z`) is therefore **closed under multiplication by every unit `a`**.

For **prime `C`**, `(ℤ/C)^*` acts **transitively** on the nonzero residues `{1,…,C−1}`. So a
unit-stable subset `Z ⊆ ℤ/C` is either `{0}` (no nonzero zeros) or all of `ℤ/C`. If `Z ⊇` all
nonzero residues, then `Φ(δ)_t = 0` for all `t=1,…,(C−1)/2`, i.e. `Mδ = 0`, so `δ = 0` (M
invertible) — contradicting `Diff ≠ ∅`. Hence `Φ(δ)` has **no** zeros among `t=1,…,(C−1)/2`. By (★),
`Φ(s)_t = 0` for **all** `t`, so `Ms = 0`, `s = 0`, i.e. `ε' = −ε`: the trivial swap. ∎

**Remark (why composite breaks it).** For composite `C`, `(ℤ/C)^*` does **not** act transitively on
nonzero residues — it preserves `gcd(t,C)`. So `Z` may be a proper nonempty union of `gcd`-classes,
which is exactly the room needed for a nontrivial `δ` with `Φ(δ)` vanishing off a divisor-class and
a complementary `s` vanishing on it — a genuine collision. This is the precise sense in which the
prime/composite dichotomy of HYP-2270 is the **vanishing-sums-of-roots-of-unity** dichotomy
(Lam–Leung / Conway–Jones): a sine sum over a half-system can vanish nontrivially exactly when the
modulus is composite.

---

## Scope, corollaries, and the open converse

- **Closes HYP-2270 forward.** Prime `C ⟹` orbit `= 2^{n−2}` is now a theorem, not a verified
  conjecture. Combined with THM-413 (`3∣C ⟹` degenerate), the only open piece of HYP-2270 is the
  **converse for composite `C` with no factor 3** (e.g. `C=25,35,49`): that these always have a
  collision. Verified `n≤22`; the construction is the subgroup-homometry of HYP-2273.

- **The collision mechanism (HYP-2273, structural).** Every observed collision flips exactly the
  **half-system of a prime-order cyclic subgroup** of `ℤ/C`: for each prime `q∣C`, the order-`q`
  subgroup `(C/q)ℤ` has half-system `H_q = {(C/q)j : j=1,…,(q−1)/2}`, and flipping the signs of
  `H_q` changes `Φ(ε)_t` only at frequencies `q∤t` (since `Σ_j sin(2π t (C/q) j /C) = Σ_j
  sin(2π t j/q) = 0` when `q∣t`). The `3∣C` single flip of THM-413 is the case `q=3`, `H_3={C/3}`.
  See HYP-2273 for the count law `defic(3p)=2^{(p−1)/2}` and prime-power refinements.

## Honest status

- **PROVED:** the homometry reframing; `|f̂_ε(t)|² = A(t)² + Φ(ε)_t²`; `M` invertible; **prime
  `C` ⟹ full sign-orbit** (Galois-stable zero set + transitivity).
- **VERIFIED:** the three-way partition equality (clock = diff-multiset = `|DFT|²`, `n≤13`); the
  closed form (all cuts, `n≤11`); orbit `=2^{n−2}` for all prime `C`, `n≤22`.
- **CONJECTURE (HYP-2273):** the subgroup-flip mechanism and the composite count law (hence the
  converse of HYP-2270).

**Artifacts:** `04-computation/signed_lrc_homometry_s705.py`, `…sine_product_s705b.py`,
`…collision_anatomy_s705c.py`, `…c39_check_s705d.py` (+`.out`s in `05-knowledge/results/`).
Reflection `07-reflections/signed-lrc-as-cyclic-homometry-s705.md`. HYP-2273, T760.
