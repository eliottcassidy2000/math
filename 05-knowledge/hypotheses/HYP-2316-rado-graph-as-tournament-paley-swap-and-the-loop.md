# HYP-2316 — The Rado graph as a tournament: the random/universal tournament, its self-converse "swap," and the Paley-7 avatar where arcs = cube roots of unity

**Session:** S638
**Status:** CONFIRMED (finite avatar formalized + verified; infinite claims are standard model theory, cited)
**Provenance forward:** math-lean `Math/Tournaments/PaleyRado.lean` (sorry-free)
**Prompt:** *"a loop of the input causes a swap of the two outputs. consider the Rado graph as a tournament."*

> **Convergence note.** opus took the *same* prompt in parallel (HYP-2282 / S699p) with a complementary
> reading: "loop ⟹ swap" = the **monodromy** of the FTA branched cover (loop the coefficients around
> the **discriminant** and two roots braid-swap = a Galois transposition; the discriminant = the
> worry-set), and Rado = the countable random tournament = the **generic fiber** (extension property =
> "always a witness"), with the worry-set fiber abelian/cyclotomic vs the generic fiber symmetric. That
> is the "two outputs = two roots" reading. **This session is the complementary "two outputs = out-set
> / in-set, swap = self-converse" reading, plus the part opus did not do: a sorry-free FORMALIZATION**
> of the finite avatar (the Paley-7 tournament) and the observation that its **arcs are the cube roots
> of unity**. The two readings agree: opus's symmetric monodromy of the generic fiber = `Aut(𝒯)`, whose
> orbit counts are the perspective numbers (§2). Cross-reference HYP-2282.

---

## 0. Reading the two clues

**"Consider the Rado graph as a tournament."** The Rado graph `R` is the unique countable graph with
the extension property (Fraïssé limit of finite graphs; almost surely the result of flipping a fair
coin for every pair). Its tournament analogue is the **countable random / universal homogeneous
tournament** `𝒯` — the Fraïssé limit of the class of finite tournaments:

- finite tournaments form a Fraïssé class (hereditary, joint embedding, **amalgamation**: amalgamate
  `B, C` over `A` by orienting all new `B\A`–`C\A` pairs `B → C`), so the limit `𝒯` exists and is the
  unique countable tournament with the **extension property**: *for all finite disjoint `U, V` there is
  a vertex `w` with `u → w` for all `u ∈ U` and `w → v` for all `v ∈ V`* (a vertex realizing any
  prescribed pattern of wins/losses);
- almost surely, orienting every pair of a countable vertex set by an independent fair coin yields a
  tournament isomorphic to `𝒯` (Borel–Cantelli on the extension axioms + back-and-forth). So `𝒯` *is*
  "the Rado graph as a tournament."

**"A loop of the input causes a swap of the two outputs."** Each vertex `v` has two "outputs": its
**out-set** `O(v) = {w : v → w}` and its **in-set** `I(v) = {w : w → v}`. The operation that **swaps
`O(v) ↔ I(v)` at every vertex simultaneously** is arc-**reversal** (the converse `𝒯*`). The clue says a
"loop" induces this swap — and indeed:

- `𝒯` is **self-converse**: reversal is an isomorphism of the Fraïssé class (a bijection on finite
  tournaments preserving embeddings), so `𝒯 ≅ 𝒯*`. The swap is an intrinsic (anti-)symmetry of the
  universal object.
- The "loop" is the project's order-2 involution **`σ : v ↦ -v`** (antipodal / complement / CM
  conjugation, HYP-2185): applying it once swaps the two outputs (reverses arcs); applying it twice
  returns the input (`σ² = id`). In the antipodal double-cover picture (HYP-2185), `σ` is the deck
  transformation — a non-trivial **loop** in the base lifts to a path that **swaps the two sheets**.

So the two clues are one statement: **the universal ("Rado") tournament carries the swap-involution
`σ` as its self-converseness, and `σ` is the loop that exchanges the two outputs.** This session pins
that down on the arc's magic prime `7`.

---

## 1. The finite avatar: the Paley tournament on `𝔽₇` (formalized)

The smallest nontrivial vertex-transitive, self-converse, "quasi-generic" tournament is the **Paley
(quadratic-residue) tournament** `QR_7`. It exists because **`7 ≡ 3 (mod 4)`** (and `7 = Φ₃(2) = N(3+ω)`
— the arc's prime, S628/S637). Orient `x → y` iff `x − y` is a nonzero square mod 7.

**The connection set IS the cube roots of unity.** The nonzero squares mod 7 are `QR = {1, 2, 4}`, and
these are **exactly the cube roots of unity** `μ₃(𝔽₇) = {x : x³ = 1}` — because `(7−1)/2 = 3` makes the
index-2 QR subgroup coincide with the order-3 cube-root group, and `2, 4` are the two *primitive* cube
roots (S637, `eval_point_two_cube_root_mod_seven`). So:

> **The Paley-7 tournament orients an arc exactly when the difference is a cube root of unity.** The
> `ω`-resonance the whole arc converges on literally *defines the edges* of the smallest random-
> tournament avatar.

**Why it is a tournament — the `σ`-disjointness.** Because `7 ≡ 3 (mod 4)`, `−1` is a **non**-residue,
so `QR` and `−QR = {3, 5, 6}` are **disjoint** and `𝔽₇* = QR ⊔ (−QR)`. That partition is *precisely*
"for `x ≠ y`, exactly one of `x → y`, `y → x`." (If `−1` were a residue — `q ≡ 1 (mod 4)` — you'd get
the Paley **graph**, undirected; the cube-root/`q ≡ 3` condition is what makes it a tournament.)

**The swap.** The antipode `σ : x ↦ −x` sends each arc to its reverse:
`σ(x) → σ(y) ⟺ (−x) − (−y) = −(x−y) ∈ QR ⟺ (x−y) ∈ −QR ⟺ y → x`. So **`σ` reverses every arc** — it is
an anti-automorphism, the exact realization of "a loop of the input swaps the two outputs," and the
finite shadow of `𝒯 ≅ 𝒯*`.

### Formalized (math-lean, sorry-free): `Math/Tournaments/PaleyRado.lean`
- `paleySet := {1,2,4}`, `paleyAdj x y := (x − y) ∈ paleySet`.
- `paleySet_eq_cube_roots : x ∈ paleySet ↔ x ^ 3 = 1` — **the arc set = `μ₃`**.
- `neg_one_not_residue : (-1 : ZMod 7) ∉ paleySet` — `7 ≡ 3 (mod 4)`.
- `paley_irrefl`, `paley_total`, `paley_asymm` — **it is a tournament** (`𝔽₇ = QR ⊔ −QR`).
- `paley_antipode_converse : paleyAdj (-x) (-y) ↔ paleyAdj y x` — **the loop `σ` swaps the outputs.**
- `antipode_involutive : -(-x) = x` — `σ² = id` (looping twice returns).

Verified numerically (`rado_tournament_paley_swap_s638.py`): `QR = μ₃ = {1,2,4}`; `−1` non-residue;
out-degrees all `3` (doubly-regular); the swap reverses all `49` arc-checks; the `(|U|+|V|≤2)`
**extension axiom holds on all 98 disjoint pairs** (Paley-7 already approximates the generic
extension property).

---

## 2. The perspective tie: orbits = iso-types, made exact

The universal tournament closes a long-standing thread (lrc-perspective-key: "perspectives" =
automorphism orbits). By **homogeneity + Ryll–Nardzewski**, the number of `Aut(𝒯)`-orbits on
`k`-element subsets **equals the number of isomorphism types of `k`-vertex tournaments** = `A000568(k)`
= `1, 1, 2, 4, 12, 56, …` (`k = 1,2,3,4,5,6`). The "perspective" counting the user keeps returning to
is, on the *generic* tournament, **exactly** the tournament-iso-type sequence — with **no finite
coincidence-breaking** (the breaks at `n = 5` in HYP-2130 are finite-size artifacts; the homogeneous
limit realizes the orbit↔type identity cleanly, that being the definition of `ℵ₀`-categoricity).

At `k = 3` the two orbits/types are **transitive vs the 3-cycle** — the 3-cycle being the cube-root
resonance (`σ³ = 1`, the same `ω` that orients Paley-7). The generic tournament is the object where
*every* finite tournament appears as an induced sub-pattern (universality) and `σ` (the swap) is a
symmetry of the whole — the perfect home for the partition-function `H` story: every finite `H`-value
is realized inside `𝒯`.

---

## 3. Why this matters for the arc

- **The swap `σ` is now an intrinsic symmetry of the universal object**, not just a per-configuration
  involution: the self-converseness of `𝒯` is the global form of the perspective key (HYP-2185's
  antipodal cover, the CM conjugation, the complement `τ`). "A loop swaps the two outputs" names the
  deck transformation.
- **Arcs = cube roots of unity** on `7`: the Paley-7 tournament fuses three arc threads — the magic
  prime `7 = Φ₃(2) = N(3+ω)` (S628/S637), the cube root `ω = μ₃ = QR`, and the `σ : v ↦ -v` swap — into
  one object. The chromatic-plane `7` (S637) and the tournament `7` are joined to the *Paley* `7`.
- **`q ≡ 3 (mod 4)` is the tournament condition** = `−1 ∉ QR` = "the swap has no fixed arc-type" — the
  same `2`-adic/antipodal seam (`⟨−1⟩` gaining/losing a fixed point) that runs through the LRC apex
  analysis (lrc-perspective-key, the `n` even ⟺ apex fixed point seam).
- **Universality of `H`**: every finite tournament (hence every achievable `H`-value, the
  strong-atom monoid) embeds in `𝒯`; the forbidden values `{7, 21}` are the *non*-embeddable
  `H`-spectrum gaps — a property of finite strong atoms, unaffected by the universal ambient.

## 4. Open / handoffs
- Formalize the **extension property** abstractly and that **reversal preserves it** (⟹ `𝒯 ≅ 𝒯*`) —
  needs Mathlib's `ModelTheory`/Fraïssé machinery (heavier; the finite avatar here is the down-payment).
- The **higher Paley primes** `q ≡ 3 (mod 4)`: `q = 11, 19, 23, …` — `q = 11` ties to the Moser-spindle
  field `ℚ(√−11)` (fleet HYP-2277); is there a Paley/Heegner bridge (the `q ≡ 3 (mod 4)` tournament
  primes vs the `4N−1` Heegner map)?
- `σ`-fixed structure: `σ : x ↦ -x` on `𝔽₇` fixes only `0` (the apex) — the antipodal cover's single
  ramification point, exactly HYP-2185. Connect the Paley apex to the LRC apex.
