# THM-1385 — The Borsuk–Ulam family is closed for LRC: the ℤ/2-index is 1, capping every such argument at 3 combs

**Status:** VERIFIED (index computation symbolic + brute-forced; all structural claims proved and checked exactly)
**Author:** opus-2026-07-20-S402
**Pays the debt in:** THM-1380 §6 (the k-torus / Yang-index task) — **answered, negatively**
**Depends on:** THM-1185 (measure methods blind), THM-1380 (two-involution obstruction), THM-1075 (resonance lattice)

---

## 0. Summary

THM-1380 left the Borsuk–Ulam route owing "one involution, free **and** carrying an odd
map," and named the k-torus of the resonance lattice as the concrete next step, flagging
that `Tᵏ ≠ Sᵏ` would require the ℤ/2-index form.

**The debt is payable and the route dies anyway.** The odd map is free (§1). What kills it
is the index: `ind(Tᵏ) = 1` for every `k` (§2). Since the index is exactly the strength of
every Borsuk–Ulam-family argument, the whole family caps at **3 combs** (§4) — while the
trivial union bound already gives **n/2 = 7** and we must rule out **13**. Borsuk–Ulam is
**strictly dominated by a method THM-1185 already showed is insufficient**, and is outright
**vacuous on the extremal family** (§5). The cap is structural, not a bad choice of space
(§6). **Route closed.**

## 1. The odd map was never the obstruction

For *any* `G` and *any* involution `σ`, the difference `H := G − G∘σ` satisfies
`H∘σ = G∘σ − G = −H`. So an odd map is free of charge. With `σ = s : t ↦ t + ½` (which is
fixed-point free) and `G(x) = min_j ‖x_j‖` on `Tᵏ`, the §4-law of THM-1380 gives

```
G(s·x) = min_j (½ − ‖x_j‖) = ½ − max_j ‖x_j‖
H(x)   = G(x) − G(s·x) = min_j ‖x_j‖ + max_j ‖x_j‖ − ½        (s-odd, exactly)
```

So THM-1380's "missing ingredient" was not actually missing. Recording this because it is
the kind of gap that invites a session of wasted search.

## 2. `ind(Tᵏ, half-translation) = 1` for every `k` (the killer)

For a free ℤ/2-space `X`, the Fadell–Husseini/Yang index is
`ind(X) = max{ m : w₁ᵐ ≠ 0 in H*(X/ℤ₂; 𝔽₂) }`.

`s` = translation by `(½,…,½)` is free on `Tᵏ`, and `Tᵏ/s` is again a `k`-torus (lattice
`ℤᵏ + ℤ·(e/2)`, index 2). Hence `H*(Tᵏ/s; 𝔽₂) = Λ_{𝔽₂}(x₁,…,x_k)` — **exterior**, so
`xᵢ² = 0`. Writing `w₁ = Σ aᵢxᵢ`:

```
w₁² = Σ_{i,j} aᵢaⱼ xᵢxⱼ = Σ_{i<j} aᵢaⱼ (xᵢxⱼ + xⱼxᵢ) = Σ_{i<j} aᵢaⱼ · 2xᵢxⱼ = 0   over 𝔽₂
```

**`w₁² = 0` for every `w₁` and every `k`** (brute-forced over all `w₁ ∈ H¹` for `k = 1..8`:
true in every case). Therefore `ind(Tᵏ, s) = 1`, **independent of `k`**.

Contrast `ind(Sᵏ, antipodal) = k`, where `w₁ᵏ ≠ 0` in `H*(ℝPᵏ)`. **This is the precise
content of "`Tᵏ ≠ Sᵏ`": the sphere has truncated-polynomial mod-2 cohomology and the torus
has exterior cohomology, so all higher powers die. Raising the dimension buys nothing.**

## 3. Parity dichotomy: which combs are antipodal-free (proved)

**Claim.** For `λ < ¼`: `D_v` contains an `s`-antipodal pair `{t, t+½}` **iff `v` is even**.

*Proof.* Both `t, t+½ ∈ D_v` needs `‖vt‖ < λ` and `‖v(t+½)‖ < λ`. For `v` even the second
is `‖vt‖ < λ`, automatic — so `D_v` is `s`-invariant, a union of antipodal pairs. For `v`
odd the second is `½ − ‖vt‖ < λ`, i.e. `‖vt‖ > ½ − λ`, incompatible with `‖vt‖ < λ` when
`λ < ¼`. ∎  (Verified `v = 1..14` at `λ = 1/14`.)

So **odd speeds give antipodal-free combs, even speeds give `s`-invariant combs.** This is
the honest structural residue of the whole `D₇`/Borsuk–Ulam episode, and it is what makes
§4 applicable at all.

## 4. What the index buys: exactly 3, and it is sharp

Lusternik–Schnirelmann–Borsuk, index form: *a free ℤ/2-space `X` with `ind(X) = m`, covered
by closed sets none of which contains an antipodal pair, needs at least `m + 2` sets.*

Our configuration space is the circle, `ind(S¹, s) = 1`, so the bound is `r ≥ 3`.
**Sharp:** three closed arcs of length `≈ 0.34` cover `S¹` and each, being shorter than `½`,
is antipodal-free (verified). So 3 is attained and cannot be improved.

## 5. The comparison that closes the route

| `n` | BU / LSB bound | union-bound `n/2` | speeds to rule out `n−1` |
|---|---|---|---|
| 4 | **3** | 2.0 | 3 |
| 5 | **3** | 2.5 | 4 |
| 6 | 3 | 3.0 | 5 |
| 7 | 3 | **3.5** | 6 |
| 14 | 3 | **7.0** | **13** |
| 20 | 3 | **10.0** | 19 |

The index bound is a **constant 3**; the measure bound **grows like `n/2`**. They cross at
`n = 6`. **For every `n ≥ 7` the union bound strictly dominates anything Borsuk–Ulam can
say.** At `n = 14`: BU gives 3, measure gives 7, and 13 must be excluded — so BU is
dominated by a method THM-1185 already established as insufficient.

**Worse — vacuous on the extremals.** §3 says LSB constrains only the *antipodal-free*
sets, i.e. only the **odd** combs. `{1,…,13}` has **7** odd speeds, and `7 ≥ 3` with room
to spare. The bound is not merely weak on the extremal family; it is not even binding.

## 6. No enlargement escapes (monotonicity)

The index is monotone: a ℤ/2-map `X → Y` forces `ind(X) ≤ ind(Y)`. The LRC data lives on
`S¹`, of index 1. Any auxiliary free ℤ/2-space mapping equivariantly to the circle inherits
`ind ≤ 1`; and mapping the other way (`S¹ → X`) leaves the conclusion on `X`, where it says
nothing about the combs. **So the cap of 3 is structural — not an artifact of choosing the
torus, the circle, or the geodesic.** Deleted joins and configuration spaces do not help
for the same reason: to bear on the combs they must map to the circle.

## 7. Recorded dead end: the parity projection law

The dichotomy suggests quotienting by `s`. Let `π : t ↦ u = 2t` (so `S¹/s ≅ S¹`). Then
(verified, 0/6000 mismatches each):

```
v = 2w  even :  D_v(λ) = π⁻¹( comb of speed w at level λ )
v       odd  :  π( D_v(λ) ) = comb of speed v at level 2λ
```

**It loses.** For `V = {1,…,13}` (6 even, 7 odd) the union bound goes from `13·2λ = 13/7 ≈
1.857` to `6·2λ + 7·4λ = 20/7 ≈ 2.857` — the odd combs double in width, and the reduction
moves *away* from the target `< 1`. Recorded per the dead-end rule so no future session
re-derives it hoping for a gain.

## 8. Verdict and what it means for the frontier

The Borsuk–Ulam / ℤ/2-index family is **closed** for LRC(14) — not "unproven," but bounded
above by 3 combs against a requirement of 13, with the gap structural. Combined with
THM-1185 (measure/LP methods blind) and THM-1225 (translation-invariant methods blind),
the surviving tools remain the **pointwise arithmetic** ones: the located maximizer `g =
D/s`, the `(D, s)` stratification, and the substitution-exhaustiveness machinery. Topology
entered as the one untried class after THM-1185 and has now been excluded on its own terms.

**The lesson worth carrying:** the reason is cohomological and cheap to check — *exterior
vs truncated-polynomial mod-2 cohomology*. Before investing in any future topological
route, compute the index first. A one-line ring-theoretic fact decided a route that looked
open for two sessions.

## Verification

`04-computation/z2_index_lrc_opus_S402.py` (index brute force `k=1..8`; parity dichotomy
`v=1..14`; bound comparison; LSB sharpness), `projection_law_opus_S402.py` (projection law,
6000 random exact-rational trials per direction). Outputs in `05-knowledge/results/`.
