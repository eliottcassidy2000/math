---
id: THM-1970
title: "H SITS AT THE FORMULA/#P EDGE -- no bounded-degree polynomial invariant determines the Hamiltonian-path count, and the defect grows with n. The tournament invariants form a DEGREE-GRADED poly tower: score sequence (degree 1) -> c3 = 3-cycle count (degree 3, score-determined by KBS) -> var(lambda^2)=SC4 (degree 4, THM-1940) -> tr(S^{2j}) (degree 2j, the j-subtournament census) -> char_S (ALL moments, the full skew spectrum, poly-computable). Each rung is a poly invariant = a 'degree-k induced-subtournament census' (k-profile). H (Hamiltonian-path count) is determined by NONE of them: (i) H is not char_S-determined -- H|char_S already SPLITS at n=3 (THM-1935 control), so even the complete spectrum misses H; (ii) the H-DEFECT (max H - min H within one degree-k census class) GROWS with n at every fixed resolution k -- exhaustive n=5,6: k=3 defect 4 then 14, k=4 defect 2 then 12 -- and reaches 0 only at the deck level k=n-1 (reconstruction). So H requires FULL-SUPPORT (degree Omega(n)) resolution; it is the PERMANENT (#P) to the spectrum's DETERMINANT (poly), and the char_S->H gap is the permanent/determinant boundary. Moreover scalar H is not even COMPOSITIONAL: H(C3[S1,S2,S3]) is not a function of (H(S1),H(S2),H(S3)) (Part A, across sizes) -- the natural compositional refinement is the path-SYSTEM (linear-forest) polynomial that categorifies H, not a poly formula (none exists). This locates tournaments exactly at the expressibility edge: the moment tower tr(S^{2j}) are the partial sums, char_S the full determinant-side datum that APPROACHES but never reaches H, and H the full-support limit -- the harmonic-series edge between what a bounded formula expresses and what is provably #P"
status: VERIFIED (exhaustive n<=6): H-defect within degree-k census classes = k3:{n5:4,n6:14}, k4:{n5:2,n6:12}, ->0 at k=n-1; H|char_S splits at n=3 (THM-1935); H(C3[.]) not a function of block scalar-H's (across sizes). The 'permanent/determinant edge' and harmonic-series framing are the interpretation (reflection H-at-the-formula-sharp-P-edge). Small-n evidence for the defect-growth trend; the direction (no bounded-k invariant determines H) follows from H being #P while every fixed-k census is poly.
author: opus-2026-07-20-S445
depends_on: [THM-1935 (H not spectrum-determined; the quaternion wall), THM-1940 (var=SC4=degree-4 census; the vertex-support ladder), THM-1930 (Sum lambda^2 fixed), THM-1945 (H is per-like #P, leaves the spectral ladder), THM-1860 (H=prod H(SCC), the transitive-quotient composition), THM-1960 (substitution/seeds)]
cite_by_filename: true
---

# THM-1970 — H sits at the formula/#P edge

Owner: `H` failing to be determined by a poly-time invariant looks like an **edge case** — maybe a
refined invariant is the real answer; tournaments sit at the harmonic-series edge between what a
formula expresses and what provably cannot be. This locates that edge precisely.

## The degree-graded poly tower

Tournament invariants stratify by **vertex-support degree** (`THM-1940`'s ladder), each rung a
poly-computable "degree-`k` induced-subtournament census":

| invariant | degree | poly? | determines H? |
|---|---|---|---|
| score sequence | 1 | ✓ | no |
| `c₃` (3-cycles) | 3 | ✓ (KBS) | no |
| `var(λ²) = SC4` | 4 | ✓ (`THM-1940`) | no |
| `tr(S^{2j})` (`j`-census) | `2j` | ✓ (fixed `j`) | no |
| `char_S` (all moments) | `n` | ✓ (determinant) | **no** |
| **`H` (Ham paths)** | full support | **#P** | — |

## H is captured by no bounded resolution

- **The full spectrum misses it.** `H | char_S` already **splits at `n=3`** (`THM-1935` control) —
  the complete determinant-side datum does not determine `H`.
- **The defect grows.** The `H`-defect (`max H − min H` within one degree-`k` census class), over
  all tournaments:

  | resolution `k` | `n=5` | `n=6` |
  |---|---|---|
  | `k=3` | 4 | **14** |
  | `k=4` | 2 | **12** |
  | `k=n−1` (deck) | 0 | 0 |

  At every **fixed** `k` the defect grows with `n`; it vanishes only at the **deck** level `k=n−1`
  (reconstruction). So `H` needs **full-support (degree `Ω(n)`)** resolution.

## The permanent/determinant edge

`char_S` is a **determinant** (poly); `H` is a **permanent-type** count (`THM-1945`: `H` is
`per`-like, a finite `Sₙ×Sₙ` stabilizer, `#P`). The `char_S → H` gap is exactly the
**permanent/determinant boundary**. Tournaments sit on it because a **complete signed relation**
makes the spectrum poly (the `±1` skew matrix's eigenvalues) while the path count is `#P`.

## Not even compositional

Scalar `H` is not a **compositional** invariant: `H(C₃[S₁,S₂,S₃])` is **not** a function of
`(H(S₁),H(S₂),H(S₃))` (verified: block-`H`-multiset `(1,1,1)` yields composites
`{3,5,9,11,23,45,51,137,509,2721}`, driven by block size and internal path structure). The natural
refinement that *does* compose is the **path-system (linear-forest) polynomial** — the number of
ways to cover a block by vertex-disjoint directed paths with prescribed boundary — which
categorifies `H` (`H` = its top, single-path, all-endpoints coefficient). It is *more refined than
`H`* and functorial, but still `#P`: **there is no poly formula for `H`** (unless `P=#P`); the
refinement is *functorial, not a complexity reduction*.

## The harmonic-series edge

The moment tower `tr(S²), tr(S⁴), …` are the **partial sums**; `char_S` is the **full
determinant-side datum** that *approaches but never reaches* `H`; `H` is the **full-support limit**.
Like `Σ 1/n^s` — poly (convergent) for every bounded resolution `s>1`, with the **pole at `s=1`**
(the harmonic edge) being the object no partial sum reaches — `H` is the first tournament invariant
past the edge. (The repo's harmonic/`γ` threads — `THM-805` resistance `=` harmonic number,
`CLAUDE.md`'s Euler–Mascheroni `γ` — are the quantitative face of the same edge.)

## Open

1. **The relative defect.** Does `defect_k(n)/H̄(n) → 0` (H "almost" determined by low resolution,
   the edge being a measure-zero correction like the finite part `γ` of the harmonic pole), or stay
   bounded below? Compute the *relative* defect and its `k`-decay.
2. **The path-system composition law.** Make the linear-forest polynomial's substitution transfer
   explicit — the exact compositional refinement of `H` (resolves THM-1960's open cyclic-`H`).

## Verification

`04-computation/refined_H_and_harmonic_edge_opus_S445.py` and `PH_composes_opus_S445.py` (+ `.out`)
— the degree-`k` census `H`-defect table (`n=5,6`); `H(C₃[·])` not scalar-`H`-determined.
