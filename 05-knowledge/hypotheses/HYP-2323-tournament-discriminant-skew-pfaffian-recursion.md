# HYP-2323 — The tournament discriminant is the skew-adjacency det: 0 (odd n) / odd-square Pf² (even n); the Pfaffian is the √-discriminant, and the recursion is parity (n→n−1) + Pfaffian (n→n−2)

**Session:** S645
**Status:** CONFIRMED (odd-n vanishing formalized; spectrum + recursion verified exhaustively n≤6)
**Provenance forward:** math-lean `Math/Tournaments/SkewDiscriminant.lean` (sorry-free)
**Answers:** S643 open thread #3 (*what is the tournament discriminant?*) + the S644 falling-factorial model.

---

## 0. The tournament discriminant

A tournament `T` has a **skew-adjacency matrix** `M`: `Mᵢⱼ = +1` if `i→j`, `−1` if `j→i`, `0` on the
diagonal, so `Mᵀ = −M`. Define the **tournament discriminant** `:= det M`. It is a genuine isomorphism
invariant (`det(PMPᵀ) = det(P)² det M = det M`), and skew-symmetry forces the discriminant dichotomy:

| `n` | `det M` | square? | √-discriminant |
|---|---|---|---|
| **odd** | `0` (formalized) | `0 = 0²` | `Pf = 0` |
| **even** | `Pf(M)²` (classical) | perfect square | `Pf(M)` = the **Pfaffian** |

So the tournament discriminant is **always a perfect square** — exactly mirroring `Aut(T) ⊆ Aₙ` (S643)
and the falling-factorial discriminant `∏k!²` (S644). The **Pfaffian** `Pf(M)` is the tournament analogue
of the **Vandermonde / √discriminant**: it flips sign under odd relabelings and is fixed under even ones
(`Aut(T) ⊆ Aₙ`).

---

## 1. The spectrum (verified exhaustively, `tournament_discriminant_recursion_s645.py`)

```
  n | #iso (A000568) | discriminant det M        | Pfaffian (even n)
  1 |  1             | {0}                       | 0
  2 |  1             | {1}                       | {1}
  3 |  2             | {0}                       | 0
  4 |  4             | {1, 9}                    | {1, 3}
  5 | 12             | {0}                       | 0
  6 | 56             | {1, 9, 25, 49, 81}        | {1, 3, 5, 7, 9}
```

Three clean facts emerge:

1. **The Pfaffian is always ODD.** `n=6`: `Pf ∈ {1,3,5,7,9}` (every odd up to 9); `n=4`: `{1,3}`. So the
   tournament √-discriminant joins `H` (Rédei, odd) and `|Aut(T)|` (odd, S643) as an **odd invariant** —
   "everything about tournaments is odd." The discriminant `det = Pf²` is an **odd square**.
2. **Maximal discriminant `= 3^{n−2}`** (powers of 3): `n=2,4,6 → 1, 9, 81 = 3^0, 3^2, 3^4`; max Pfaffian
   `= 3^{(n−2)/2} = 1, 3, 9`. The extremal (most-regular / doubly-regular) tournament's discriminant is a
   **power of 3** — the cube-root prime (`Φ₃`, the arc's `ω`) governs the discriminant's extreme.
3. **Rank `M = n − [n odd]`.** The skew rank is `n` (even, `det ≠ 0`) or `n−1` (odd, `det = 0`):
   discriminant vanishes ⟺ rank-deficient-by-1 ⟺ `n` odd.

---

## 2. The recursion in n (the heart of the question)

The tournament discriminant changes with `n` by **two interleaved recursions** — exactly the repo's two
reduction modes (CLAUDE.md):

- **Mode A (`n → n−1`, add/remove a vertex): the PARITY FLIP.** Each vertex toggles `n`'s parity, so the
  discriminant alternates `0` (odd, rank-deficient) ↔ `Pf²` (even, full rank). The "discriminant on/off"
  is the vertex-insertion recursion.
- **Mode B (`n → n−2`, both legs / Cayley–Dickson descent): the PFAFFIAN RECURSION.** The √-discriminant
  expands `Pf(M) = Σⱼ (−1)ʲ M₁ⱼ · Pf(M_{1̂ĵ})` — the `(n−2)`-Pfaffian cofactor expansion. So the even-`n`
  discriminant is *built recursively from the `(n−2)` discriminants*, the slow time scale.

> **Formalized (math-lean, sorry-free): `Math/Tournaments/SkewDiscriminant.lean`** — `skew_odd_det_zero`
> (`Mᵀ = −M`, `Odd n` ⟹ `det M = 0`, via `det M = det Mᵀ = det(−M) = (−1)ⁿ det M = −det M`); `disc_two`
> (`n=2`: `det = 1 = Pf²`). The even-`n` `det = Pf²` is classical (Mathlib lacks the Pfaffian; flagged).

---

## 3. Other tournament invariants and their recursion (the broader picture)

| invariant | `n = 1..6` | recursion / law |
|---|---|---|
| #iso classes | `1,1,2,4,12,56` | **A000568** — the master count (Burnside on `Sₙ`, CLAUDE.md) |
| `max |Aut(T)|` | `1,1,3,3,5,9` | odd (S643); `Aut ⊆ Aₙ`; grows slowly |
| `H` (Ham-paths, Rédei) | all **odd** | deletion–contraction / tie-induction (S625/S633) |
| discriminant `det M` | `0` / odd-square | **parity (A) + Pfaffian (B)** — this session |
| `Pf(M)` (√-disc) | **odd**, max `3^{(n−2)/2}` | Pfaffian recursion `n→n−2` |
| `rank M` | `n − [n odd]` | parity (rank drops iff `n` odd) |

The unifying picture: **the odd/even parity of `n` is the master switch** — it flips the discriminant
(0 ↔ square), the rank (deficient ↔ full), and the Pfaffian (0 ↔ odd). The three "odd" invariants
(`H`, `|Aut|`, `Pf`) are the tournament's signature of living in the **alternating/sign-kernel** world
(S643/S644): everything is odd because the swap is forbidden.

---

## 4. Ties & new threads
- **Closes S643 thread #3:** the tournament discriminant = `det(skew-adj)`, always a square, with the
  Pfaffian as √-disc — the exact analogue of "disc is a square ⟺ Galois ⊆ Aₙ," and it holds *always*
  because `Aut(T) ⊆ Aₙ` always.
- **Ties S644:** the falling factorial `(x)_n` has discriminant `∏k!²` (a square) — the *polynomial*
  model; the tournament has discriminant `Pf²` (a square) — the *combinatorial* model; both are the
  "even/sign-kernel" side.
- **New: max discriminant = `3^{n−2}` (power of 3).** The extremal tournament's discriminant is governed
  by the cube-root prime — does the maximizer = the doubly-regular / Paley-type tournament, and is the
  `3`-power exact for all even `n`? (ties the Φ₃/forbidden-`7` cube-root line, S637/8/40/42/43.)
- **New: the Pfaffian's combinatorial meaning** — `Pf(M) = Σ` signed perfect matchings of `K_n` weighted
  by the orientation. Is the tournament Pfaffian a known invariant (signed 1-factor count)? Worth
  identifying.
- **Handoff:** formalize the even-`n` `det = Pf²` (needs a Pfaffian in Mathlib — a real gap to fill); and
  the rank law `rank M = n − [n odd]`.
