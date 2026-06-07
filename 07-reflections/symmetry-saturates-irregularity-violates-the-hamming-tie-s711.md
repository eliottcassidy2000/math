---
source: monad-explorer-2026-06-07-S711 (deep-research; OPEN-Q-057 frontier)
status: REFLECTION on THM-432 (verified) + a cross-domain pattern (CONJECTURE-grade
  analogy, honestly flagged) tying the unit-distance 3N-crossover to the LRC worry-set
  via the recurring shape "the symmetric construction saturates the bound; only an
  irregular one violates it." Grounded in exact computation (THM-432) on the UD side;
  the LRC side cites repo results (opus-S699 signed-LRC AP/V* split, THM-425/420).
tags: [unit-distance, hamming-graph, cartesian-product, kissing-number, 3N-crossover,
  LRC, lonely-runner, worry-set, AP, V-star, cyclotomic-cayley, HYP-2170, meta-pattern,
  symmetry-saturates]
---

# Symmetry saturates, irregularity violates — the Hamming tie at n = 27

## The fact (exact, THM-432)

The unit-distance "tie at `n = 27`" — the one mysterious `+0` in THM-431's deficit
ladder `−6,−5,−4,−3,−2,+0,+1` — is not a numerical accident. It is the **Hamming
graph `H(3,3) = K₃ □ K₃ □ K₃`**, the 3-fold Minkowski sum of unit equilateral
triangles: `27` points, **6-regular**, exactly `81 = 3·27` unit distances (verified
exactly in `Q(√3)`). The `3³` everyone kept noticing *is* `K₃^□3`.

And the reason it **ties rather than beats** is a single sentence: *a Cartesian
product of triangles is forced 6-regular, and `6 = κ` is the kissing cap, so every
vertex hits the threshold and none can exceed it.* To get `u(N) > 3N` you must
produce a vertex of degree `≥ 7`; a product's average is dragged *below* `κ` by each
factor's low-degree boundary, so the product is structurally pinned at-or-below the
line. The first crossing `N* ∈ [25,28]` is therefore necessarily a **non-product,
irregular blob** (THM-432(C): the smallest product that beats is `N = 32`).

## The shape

Strip the geometry and a clean shape remains. There is a **threshold** (`3N`, set by
the kissing number `κ`). There is a **maximally symmetric construction** sitting on a
divisor-structured `N` (`27 = 3³`) that meets the threshold **with perfect
equality**, every vertex identical, because symmetry forces *exact* regularity at the
critical degree. And there is the **theorem that this symmetric object can never cross
the threshold** — equality is the most it can do. Crossing is reserved for objects
that **break the symmetry** and concentrate the violated quantity on a few sites.

> **Symmetry saturates the bound; only irregularity can violate it.**

The symmetric extremal is always the *boundary case*, never the *counterexample*.

## Where else this exact shape lives in the repo: the LRC worry-set

The project's other great "threshold" problem is the **Lonely Runner**: the
conjectured floor is `M ≥ 1/n`, and the **worry-set** is the set of speed
configurations achieving `M = 1/n` exactly (sitting *on* the threshold). A
counterexample to LRC would be a config with `M < 1/n` — i.e. one that *crosses*.

The repo's signed-LRC work (opus-S699; THM-425/THM-420) found that the worry-set
**splits into two species** at `n = 14`:

- **AP** `= {1,…,n−1}` — the arithmetic progression. Maximally symmetric / structured
  / "the lattice of speeds." It sits exactly on the floor (`M = 1/n`) and is
  **shell-partner-free** (no pair `a+b ≡ 0 mod (2n−1)`; max sum `2n−3 < 2n−1`).
- **V\*** `= {1,…,11,13,24}` (AP with `12` doubled to `24`) — *not* the AP, carries a
  **shell-partner** `(3,24)`, `3+24 = 27 = 2·14−1`. Also exactly on the floor.

Now overlay the unit-distance shape:

| role | unit distance (THM-432) | LRC worry-set (S699/THM-425) |
|---|---|---|
| threshold | `3N` (kissing `κ=6`) | floor `1/n` |
| symmetric construction *on* the line | `H(3,3)=K₃^□3`, 6-regular | **AP** `{1..n−1}`, shell-free |
| why it cannot cross | forced `deg ≡ κ` (regular) | structured / shell-free; the gauge "cut" finds no across-cut sum to exploit |
| the crosser must be… | irregular blob, some `deg ≥ 7` | non-AP, shell-partner-carrying (V\*-type), a *doubled* speed |
| arithmetic of the special `n` | `27 = 3³` (`K₃` cubed) | `n=14`: `2n−1 = 27 = 3³` (the prime-cube shell) |

The same number `27 = 3³` sits at the crux of **both** problems — as `K₃^□3` (a
product of three triangles, on the construction side of unit distance) and as the
`2n−1 = 27` shell modulus (the homometry/shell side of LRC at `n=14`). THM-427's
"two-tower witness group" already noted the LRC `27` as the `3-adic` shell tower; here
the *same* `27` appears as the unit-distance product `K₃^□3`. Whether this is one
structure seen twice or a coincidence of small numbers, I flag honestly as open — but
the **shape** (symmetric-saturator vs irregular-crosser) is identical, and that is the
durable content.

## What the analogy predicts (testable, honestly CONJECTURE)

If the shape is real, it makes a transferable prediction in each direction:

1. **UD → LRC.** The unit-distance crosser must *break product structure and
   concentrate degree*. The LRC analogue: a worry-set config that goes *below* `1/n`
   (a counterexample), if it exists, cannot be an AP or any direct-product /
   group-structured speed set — it must be a "blob," a speed set whose pair-sum
   (shell / additive-energy, S674) structure is irregular, with the deficit
   concentrated on a few doubled / shell-partner pairs (the V\*-doubling motif). This
   says: **search for LRC counterexamples among shell-partner-RICH, non-AP sets, never
   among APs or product configs** — the symmetric ones are provably-on-the-line by the
   saturation half of the pattern. (Consistent with: AP is a *good* config; the
   worry lives off-AP, S3/MISTAKE on "tight=consecutive-block.")

2. **LRC → UD.** The LRC worry-set splits *exactly at the prime-cube* `2n−1 = 27 = 3³`
   (the first `n` whose shell admits a doubled-speed partner, THM-427). The UD
   analogue would be: the product/symmetric UD construction is *tight with the global
   optimum* exactly at the prime-power `N = 3³ = 27` and starts losing to irregular
   blobs immediately after — which is **exactly what THM-432(C) measures** (product
   deficit `+0` at 27, then beaten by `1` at 28, the blob `H(3,3)+`?). So the LRC
   "split at `3³`" and the UD "product tight at `3³`" are the same event read on two
   problems: the `3³` object is the last place pure symmetry suffices.

## A concrete next probe (for the UD ceiling, OPEN-Q-057)

The blob `u(28) ≥ 85 = 81 + 4` beats `3·28 = 84` by exactly `1`, and `85 = e(H(3,3)) +
4`. This invites the question: **is the `n=28` crosser `H(3,3)` plus one extra point
adjacent to `4` of its vertices?** Generically a point in the plane has only 2 DOF, so
being unit-distant from `4` prescribed vertices is over-determined — but at the
*special* angles that maximise coincidences it may close. If so, `N* = 28` gets a
transparent "`H(3,3) + 1`" construction and the tie/break boundary becomes a single
added vertex. If not, the `28`-blob is genuinely un-product-like, sharpening the
"irregularity is essential" reading. Either outcome is durable. (Pure-product attempts
are futile below `N = 32` — THM-432(C).)

## The take-away

Two of the project's hardest thresholds — the unit-distance kissing line `3N` and the
Lonely-Runner floor `1/n` — share a single architecture: a **maximally symmetric,
divisor-structured construction saturates the line with perfect regularity and is
provably incapable of crossing it; crossing is the exclusive privilege of
symmetry-broken, irregular, degree/shell-concentrated objects.** The Hamming graph
`H(3,3) = K₃^□3` is the unit-distance avatar of this principle, and `27 = 3³` is where
both problems say, in their own language, *"this is the last rung pure symmetry can
reach."* Follow the irregular blob, not the lattice, to find the violation.
