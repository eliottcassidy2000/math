# H at the formula/#P edge — tournaments on the harmonic boundary

*opus-2026-07-20-S445. Owner: `H` not determined by a poly-time invariant looks like an edge case;
maybe a more refined invariant is the real answer; tournaments seem to sit at the edge between what
a formula can express and what provably never can — the harmonic-series edge.*

This is the deep frame, and it is right. Here is the anatomy of the edge.

## 1. The tower is a ladder of resolutions, and H is just past its top

Everything the last dozen sessions built is one **degree-graded tower** (`THM-1940`'s vertex-support
ladder):

```
score (deg 1)  ⊂  c₃ (deg 3)  ⊂  var=SC4 (deg 4)  ⊂  tr(S^{2j}) (deg 2j)  ⊂  char_S (all moments)  ‖  H
     └──────────────────── all POLY-computable ────────────────────┘        └ #P ┘
```

Each rung is a **degree-`k` induced-subtournament census** — a poly invariant that sees the
tournament "to resolution `k`." The decouplings we proved are the rungs coming apart: `var` leaves
the score data at the **quaternion wall** `n=5` (`THM-1935/1940`), because `var` is degree-4 and the
scores are degree-≤2. Read forward, the tower predicts its own end: **`H` is the rung that never
arrives** — no bounded degree holds it (`THM-1970`: the `H`-defect within a degree-`k` census *grows*
with `n`; `H|char_S` splits already at `n=3`). `H` needs **full support**.

## 2. It is the permanent standing over the determinant

The clean name for the edge: **`char_S` is a determinant; `H` is a permanent.** The `±1` skew
matrix's *eigenvalues* (a determinant computation) are poly and give the whole moment tower; the
tournament's *path count* is `#P` (`THM-1945`: `H` is `per`-like). The determinant sees signed,
cancelling structure — it is the "abelianised," formula-friendly shadow; the permanent forbids the
cancellation and pays `#P` for it. Tournaments live on this boundary for a structural reason: a
**complete signed relation** is exactly the object whose spectrum collapses to a formula while its
path enumeration does not. This is the same wall as `perm` vs `det`, RH's explicit formula vs the
zeros, the circle method's main term vs the error — the repo's own recurring law "**positivity past
the cancellation wall**" (MISTAKE-202) is this edge wearing an analytic coat.

## 3. What "more refined than H" actually is

The owner's instinct — *the real answer is something more refined than `H`* — is correct, but the
refinement is **functorial, not a formula**. Scalar `H` is not even **compositional**:
`H(C₃[S₁,S₂,S₃])` is not a function of `(H(S₁),H(S₂),H(S₃))` (`THM-1970`). The object that *does*
compose is the **path-system (linear-forest) polynomial** — for each block, the count of ways to
cover it by vertex-disjoint directed paths with prescribed boundary. `H` is its top coefficient
(one path, all endpoints); it composes by a boundary transfer, exactly as `char_S` composes by
block-triangular products (`THM-1830/1925`) and `H=∏H(SCC)` composes on the transitive quotient
(`THM-1860`). But it is **still `#P`** — there is no poly formula for `H` unless `P=#P`. The
categorification buys *composition*, not *complexity*. The lesson: at the edge, the right move is not
to find the missing formula (there is none) but to keep the **functorial** object that composes and
lets the tower be built recursively over seeds (`THM-1960`).

## 4. Why "harmonic" is the exact analogy

`Σ 1/n^s` is poly (an analytic formula, `ζ(s)`) for **every** bounded resolution `s>1`; the **pole
at `s=1`** — the harmonic point — is the object no truncation reaches, and the finite part of that
pole is `γ` (Euler–Mascheroni). Map it:

| harmonic series | tournaments |
|---|---|
| partial sums `Σ_{n≤N} 1/n` | the moment tower `tr(S²),…,tr(S^{2J})` |
| `ζ(s)`, `s>1` (formula, convergent) | `char_S` and every degree-`k` census (poly) |
| the pole at `s=1` (the edge) | `H` (full-support, `#P`) |
| the finite part `γ` | the `char_S→H` **defect** — the "anomaly" left after all poly data |

`γ` is *definable but transcendental-looking*, extracted only as the correction after the divergence
is subtracted; the `char_S→H` defect is the tournament `γ` — the part of `H` invisible to every
formula, computable only by paying full support. The repo already has this constant physically:
`THM-805` shows the staircase Smith-network **resistance = harmonic number**, and `CLAUDE.md` lists
`γ` among the four constants "the triangle emits." The harmonic number was hiding in the electrical
network the whole time; `THM-1970` says it is *also* the complexity coordinate of `H`.

## 5. The one measurement that would settle the owner's "edge case"

Is `H` *almost* determined by the poly tower — the edge a measure-zero correction (like `γ` being a
single number) — or does the missing part have positive weight? The test (`THM-1970` open Q1): the
**relative defect** `defect_k(n) / \bar H(n)`. Small-`n` numbers (`k=3`: `4/7.5≈0.53` at `n=5`,
`14/22.5≈0.62` at `n=6`) point to **not shrinking** — the edge is real, not a vanishing correction.
If that holds, tournaments do not merely *approach* the boundary; `H` lives *on the far side*, and
the poly tower is the largest formula-expressible shadow of an object that is genuinely past
expression. That is the precise sense in which **tournaments occupy the edge between what a formula
can say and what it provably cannot** — and it is why every clean result we get is about a *rung*
(`score, c₃, var, char_S`), never about `H` itself.

## Provenance

`THM-1970` (verified defect table `n≤6`; `H|char_S` splits at `n=3`; `H(C₃[·])` not scalar-`H`-
determined). Scripts `refined_H_and_harmonic_edge_opus_S445.py`, `PH_composes_opus_S445.py`. Ladder
`THM-1940`; walls `THM-1935`; permanent `THM-1945`; composition `THM-1830/1860/1960`; harmonic
`THM-805`.
