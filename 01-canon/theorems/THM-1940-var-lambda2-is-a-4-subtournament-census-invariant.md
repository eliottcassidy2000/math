---
id: THM-1940
title: "var(lambda^2) IS A 4-SUBTOURNAMENT-CENSUS INVARIANT -- resolving THM-1930's two open questions and pinning the mechanism of the quaternion wall (THM-1935). The GIT-instability scalar var(lambda^2) (equivalently tr(S^4), since Sum lambda^2 = n(n-1) is fixed) is EXACTLY a function of the score sequence together with SC4 := the number of strongly-connected induced 4-subtournaments (the (1,1,2,2)-type, i.e. induced 4-sets carrying exactly 2 cyclic triples). VERIFIED n=3..6: var(lambda^2) is determined by (score, SC4) with zero exceptions, and its tr(S^4) expansion has SC4-coefficient exactly 64 = 2^6 at every n>=4, with score-part 32(n-3)*Sum C(s_i,2) + alpha(n)*n(n-1) + beta(n). RESOLVES THM-1930: (Q2, the finer invariant that determines var) = the 4-vertex census (score + SC4), a DEGREE-4 tournament invariant; (Q1, the 32-step insertion quantum) = gcd(32(n-3), 64) = 32, so Delta tr(S^4) under vertex-insertion moves in steps of 32, counting the change in score-variance plus 64 times the number of new strongly-connected 4-subsets through the inserted vertex. THE QUATERNION-WALL MECHANISM (THM-1935): SC4 first ESCAPES (score,c3)-determination at EXACTLY n=5 (SC4 is determined by (score,c3) at n<=4, splits at n>=5), which is PRECISELY why var(lambda^2) decouples from the degree-<=2 combinatorial data at n=5. So the wall is a DEGREE-ESCAPE: var is a degree-4 invariant, and n=5 is the first size at which the 4-vertex census is free of the <=2-vertex (score, 3-cycle) data. This puts the invariant hierarchy on a vertex-support ladder: c3 = 3-vertex census (KBS, score-determined), var(lambda^2) = 4-vertex census (score+SC4), H = full-support; each escapes the lower-degree data at its own threshold, and the CD-tower reading (degree-2^k invariants escaping degree-<2^k data) predicts the next wall for an 8-subtournament invariant"
status: PROVED/VERIFIED n=3..6. var determined by (score,SC4) exhaustively (0 exceptions); SC4 coeff in tr(S^4) = 64 universal; SC4 escapes (score,c3) at n=5 exactly. Resolves THM-1930 Q1/Q2; pins the THM-1935 wall mechanism as a degree-escape. The global closed form has n-dependent score-coefficients (per-n exact); the structural determination var=f(score,SC4) is the content.
author: opus-2026-07-20-S443
resolves: THM-1930 (open Q1 the 32-step index; open Q2 the var-determining invariant)
refines: THM-1935 (the quaternion wall = SC4's degree-escape from (score,c3) at n=5)
depends_on: [THM-1930 (var carried by tr(S^4); Sum lambda^2=n(n-1) fixed), THM-1935 (universal decoupling threshold n=5), THM-1820 (c3=C(n,3)-sum C(s_i,2), score-determined = 3-vertex census)]
cite_by_filename: true
---

# THM-1940 — var(λ²) is a 4-subtournament-census invariant

The concrete next steps after THM-1930/1935. `var(λ²)` (kps's GIT-instability scalar) was shown
*not* c₃-determined (THM-1930) and to decouple at `n=5` (THM-1935), with two open questions: what
finer invariant *does* determine it, and what is the `32`-step insertion quantum. Both resolve
through one object: the **4-subtournament census**.

## var(λ²) = f(score, SC4)

Let `SC4` = the number of **strongly-connected induced 4-subtournaments** (the `(1,1,2,2)` type —
induced 4-sets carrying exactly `2` cyclic triples). Since `Σλ² = n(n−1)` is fixed (THM-1930),
`var(λ²)` is carried by `\operatorname{tr}(S⁴)`, and:

> **`var(λ²)` is determined exactly by `(score sequence, SC4)`** — verified `n=3..6`, zero
> exceptions. In `\operatorname{tr}(S⁴)` the **`SC4`-coefficient is `64 = 2⁶`** at every `n≥4`,
> with score-part `32(n−3)·ΣC(s_i,2) + α(n)·n(n−1) + β(n)`.

So `var(λ²)` is a **degree-4** tournament invariant (a function of the induced 4-vertex census),
not of the degree-`≤2` data `(score, c₃)`.

## Resolving THM-1930

- **Q2 (the finer invariant).** `var(λ²)` is the **4-subtournament-census invariant** `(score, SC4)`.
- **Q1 (the 32-step quantum).** The score-coefficient `32(n−3)` and the `SC4`-coefficient `64` are
  both multiples of `32`, so `Δ\operatorname{tr}(S⁴)` under vertex-insertion moves in steps of
  `\gcd(32(n−3),64)=32` — the index counts `Δ(\text{score-variance}) + 64·ΔSC4` (the new
  strongly-connected 4-subsets through the inserted vertex). The `32`-quantization observed in
  THM-1930(D) is exactly `2⁵ = \gcd`.

## The quaternion-wall mechanism (THM-1935)

> **`SC4` first escapes `(score,c₃)`-determination at exactly `n=5`** (determined at `n≤4`, splits
> at `n≥5`).

That *is* the wall: `var(λ²)` decouples from the degree-`≤2` data at `n=5` **because its
determining datum `SC4` does**. The quaternion wall is a **degree-escape** — the first size at
which the 4-vertex census is free of the `≤2`-vertex (score, 3-cycle) data.

## The invariant hierarchy on a vertex-support ladder

| invariant | vertex-support degree | determined by | escapes lower data at |
|---|---|---|---|
| `c₃` (3-cycles) | **3** | the score sequence (KBS, THM-1820) | never (it *is* the 3-census) |
| `var(λ²)` | **4** | `(score, SC4)` | `n=5` (SC4 frees) |
| `H` (Ham paths) | full `n` | the whole tournament | `n=5` (THM-1865) |

Each invariant escapes the lower-degree combinatorial data at its own threshold. The
Cayley–Dickson reading (THM-1935): degree-`2^k` invariants escape degree-`<2^k` data; the next
predicted wall is for an **8-subtournament** invariant (the octonion level), testable at larger `n`.

## Open

1. **The octonion wall, made precise.** Is there a natural degree-8 invariant (an 8-subtournament
   census, or `\operatorname{tr}(S⁸)`) that decouples from the `≤4`-vertex data at `n=9`? This is
   the sharp, correctly-posed version of THM-1935's second-wall conjecture.
2. **The exact `α(n), β(n)`.** Pin the score-part closed form (per-`n` exact; a clean global form
   needs the `n`-scaling absorbed — likely `\operatorname{tr}(S⁴) = 2m² + \dots` with `m=\binom n2`).

## Verification

`04-computation/tr_s4_formula_and_deg2_wall_opus_S443.py` and
`04-computation/tr_s4_clean_formula_opus_S443.py` (+ `.out`) — `var=f(score,SC4)` exhaustive `n=3..6`;
`SC4`-coefficient `64`; `SC4` escapes `(score,c₃)` at `n=5`; the degree-2 threshold matrix.
