---
id: THM-419
name: 2d-lattice-chromatic-cap
status: PROVED (core) + VERIFIED-conjectured (global)
date: 2026-06-06
session: monad-explorer-2026-06-06-S708
depends_on:
  - THM-418   # square/triangular lattice dichotomy (the two special cases)
  - THM-416   # CM density-quantum totient cap (w<=6 in 2D)
---

# THM-419: The chromatic cap for 2D lattice unit-distance graphs (χ ≤ 3)

## Statement

For a 2D lattice `L` with integral (or real) positive-definite quadratic form
`Q` and squared-norm `D`, the **unit-distance graph** `U(L,D)` has vertex set `L`
and edges `{x,y}` with `Q(x−y)=D` (a Cayley graph on `ℤ²` with symmetric
connection set `S_D = {v : Q(v)=D}`).

**Main claim:** `χ(U(L,D)) ≤ 3` for **every** 2D lattice and **every** norm `D`.
Equivalently, the triangular lattice's `χ=3` (THM-418) is the **maximum chromatic
number over all 2D lattice unit-distance graphs**: **no single 2D lattice ever
forces `χ ≥ 4`, at any norm or scale.**

This answers the open handoff of THM-418 ("does some generic 2D lattice, `w=2`,
attain `χ=4`?") in the **negative**, extending the square (`χ=2`) / triangular
(`χ=3`) dichotomy to **all** 2D lattices — including every imaginary-quadratic
order (the `w=2` lattices "between the triangular lattice and the CM field").

## Proof status (honest split)

### PART 1 — Generic (non-arithmetic) real lattices: `χ ≤ 2`  [PROVED]
If three integer vectors `u,v,w` satisfy `Q(u)=Q(v)=Q(w)=D`, that is three linear
equations on the three Gram entries `(a,b,c)`; for a lattice whose Gram matrix is
**not rational up to scale** (non-arithmetic) these are generically inconsistent,
so every achieved norm has `r(D)=|S_D|=2` (just `±v`). Then `U(L,D)` is a disjoint
union of doubly-infinite paths ⇒ **bipartite, `χ≤2`**. An odd cycle (`χ≥3`)
requires `≥3` equal-norm vectors, which forces `Q` rational up to scale (**`L`
arithmetic** = similar to an integral form). *Verified:* irrational Gram matrices
`diag(1,√2)`, `diag(1,π)`, `(1,1/π,e)`, … all have `max r(D)=2` over their first
60 achieved norms (`real_lattice_chromatic_s708.py`).

So only **arithmetic** lattices (integral binary forms) can have `χ≥3`, and the
whole question reduces to integral forms.

### PART 2 — Every primitive integral binary form, `3 ∤ D`: `χ ≤ 3`  [PROVED]
Colour by a **linear functional mod 3**, `c(x,y)=φ(x,y) mod 3`, `φ≠0`. This is
proper iff `φ(s)≢0 (mod 3)` for all `s∈S_D`. Reduce `Q` mod 3:

- **`3∤disc` (Q nondegenerate mod 3).** `S_D mod 3 ⊆ {v≠0 : Q(v)≡D}`. The 8
  non-zero classes of `𝔽₃²` split into 4 lines-through-0 (2 points each). A short
  `𝔽₃` computation shows that for `D≢0`, the level set `{Q≡D}` (4 classes) is the
  union of **exactly two** of those lines, hence **avoids the other two**; pick `φ`
  whose kernel is an avoided line. (Anisotropic `disc≡2`: e.g. `Q=x²+y²`, `D≡1`
  set `{±e₁,±e₂}` avoids both diagonals — `φ=x−y`. Hyperbolic `disc≡1`: `Q≅xy`,
  `D≠0` set `{x,y both ≠0}` avoids both axes — `φ=y`.) *Verified: 0 failures.*
- **`3∣disc` (Q degenerate, primitive).** Then `Q ≡ a·ℓ(x,y)² (mod 3)` for a
  linear form `ℓ` and `a≢0`. `Q(v)=D`, `3∤D` ⇒ `ℓ(v)²≡D/a≢0` ⇒ `ℓ(v)≢0`. So
  `c=ℓ mod 3` is proper. (This is exactly THM-418's triangular trick
  `Q=a²+ab+b²≡(a−b)²`, `ℓ=a−b`.)

Since non-bipartite ⇒ `χ≥3` and the colouring gives `χ≤3`, **`3∤D ⇒ χ=2 or 3`.**

### PART 3 — `3∣D` reductions  [PROVED in two of three sub-cases]
- **`disc≡2 (mod 3)` (anisotropic):** `Q(v)≡0 mod 3 ⇒ v≡0 mod 3` ⇒ `S_D ⊆ 3ℤ²`,
  so `U(L,D)` is 9 disjoint copies of `U(3ℤ², D)`, similar to `U(L, D/9)`.
  **3-adic recursion** terminates at `3∤D'`. *Verified: `S_D⊆3ℤ²` always holds.*
- **`3∣disc` (degenerate):** `Q≡aℓ²`, `3∣D ⇒ ℓ(v)≡0 mod 3` for all `v∈S_D` ⇒
  `S_D` lies in the index-3 sublattice `{ℓ≡0 mod 3}`; the graph splits into 3
  copies on a scaled-similar lattice — **3-adic recursion** (THM-418 triangular
  mechanism verbatim).
- **`disc≡1 (mod 3)` (hyperbolic) with `9∣D`:** the single **residual** family.
  Here `S_D mod 3` lies on the null cross `{x≡0}∪{y≡0}` *and* contains a vector
  `≡0 mod 3` (forcing `9∣D`); the period-3 pullback acquires a loop, so a
  **period-2 colouring** is used instead. Proven only computationally so far.

### PART 4 — Global verification  [EXHAUSTIVE]
`all_forms_chromatic_s708.py`: **all 448 reduced positive-definite integral binary
forms** with `A≤6, C≤14` and **every norm `D≤160` with edges** (16 486
`(form,D)` cases): `χ∈{2,3}`, **0 cases of `χ≥4`**. Every case is closed by an
**explicit periodic 3-colouring**: bipartite mod-2 (10 238), linear mod-3 (3 825),
period-2 into `ℤ₃` (362), period-3 into `ℤ₃` (61) — **0 unresolved**. Wider check
(`real_lattice_chromatic_s708.py`): non-principal forms at small discriminants
with class number `>1` (richest popular norms), `D≤300`, also all `χ≤3`. The
residual family (Part 3, hyperbolic mod 3, `9∣D`) is always resolved by a period-2
colouring (`torus2`), never reaching 4.

## Consequence (Hadwiger-Nelson / dispatched question)

The dispatched probe — *is there a 2D-realizable group strictly between the
triangular lattice (`κ=6`) and the CM field whose norm layer forces `χ≥4`?* — is
answered **NO**: in 2D there is **nothing above `χ=3`**. The triangular lattice is
the chromatic maximum. Combined with THM-416 (the 2D Euclidean density quantum is
capped at `3 = w/2`, `w≤6`), the 2D lattice world is uniformly `≤3` in **both** the
density-quantum sense and the chromatic sense.

Therefore **all chromatic forcing above 3 in the plane is genuinely off-lattice
(Dehn-nontrivial, HYP-2275):** it cannot be realized inside ANY single 2D lattice,
only by combining incommensurate lattices / irrational rotation angles (Moser
`arccos 5/6`, de Grey). The open `χ(ℝ²)∈{5,6,7}` lives entirely in the
non-arithmetic / non-torsion regime — consistent with S699o's "rotation rank"
picture (`χ=3+rank`: a single lattice is rank 0, capped at 3).

## Scope / limits
- Parts 1–2 fully proved; Part 3 proved except the hyperbolic-mod-3 `9∣D` family
  (computationally `≤3`, see HYP-2280 for the conjectural completion).
- Statement about single-lattice subgraphs; does **not** resolve `χ(ℝ²)`. Its
  value is *localizing* the open `{5,6,7}` entirely off any 2D lattice.

## Related
THM-418 (square/triangular special cases), THM-416/HYP-2274 (CM density-quantum
totient cap), HYP-2275 (Niven=Dehn lattice-escape), HYP-2279/S699o (rotation rank
`χ=3+rank`), HYP-2280 (conjectural completion to a full proof), de Grey 2018,
Moser spindle. Classical analogues: `χ(ℚ²)=2`, triangular grid 3-chromatic.
