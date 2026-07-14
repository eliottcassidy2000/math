# The Atkin-Lehner V₄ does not descend to the tournament metagraph — S_n-transitivity makes complement irreducible, and the n=4 match is the base coincidence

*mac-mini-2026-07-13-S91. Owner asked to explore the V₄/Atkin-Lehner thread — the one survivor of the
S90 "odd-graph → cusp-form" build attempt — on the metagraph side. The honest result: the Atkin-Lehner
Klein four-group of X₀(14) does **not** descend to the tournament iso-class metagraph. There is exactly
one iso-class involution (complement), it is **irreducible** (cannot factor into two), and its
fixed-point structure is even incompatible with being a single Atkin-Lehner element. The tempting n=4
match is exact and beautiful — and breaks at n=6. Below: what descends, what doesn't, why, and how it
unifies the whole tournament↔modular severance under one cause.*

---

## The established base (do not re-derive)

- **X₀(14) = 14a**, genus 1, rank 0; its 4 cusps carry the Atkin-Lehner group `W(14) = {1,W₂,W₇,W₁₄} =
  (Z/2)²`, and klein-S10 identified the 4 cusps with the 4 **n=4 tournament classes** `{T,+,−,S}`
  (HYP-3586). Hardness = the genus-1 cusp form `f₁₄` (opus, 2026-06-30).
- **THM-280 (PROVED):** the tiling grid-reflection `(x,y)→(n+1−y,n+1−x)` induces `T→T^op` at the class
  level. So on iso-classes **complement = transpose = grid-reflection = the fold** — they are *one* map.
- **THM-584 (PROVED):** complement = the antipodal map `x→x⊕1` of the arc-hypercube `Q_{C(n,2)}`; it
  fixes the self-complementary (SC) classes and swaps the non-self-complementary (NS) complement-pairs;
  its ±1 eigenspaces are the R-even (Perron/bulk) and R-odd blocks.

The correspondence has been flagged, correctly, as "numerology, unbuilt, coincidence-at-14" — with a
stack of proved honest negatives (opus: complement ≠ `[−1]`; floor ≠ `L(14a,1)` ≠ `L(sym²f₁₄)`; S90:
the odd-graph → cusp transfer is provably blind). This note explains the metagraph-side *reason*.

---

## Three metagraph-side facts

### 1. Complement is irreducible under S_n — the V₄ cannot form (new)

A genuine V₄ needs a **second** involution `W` with `complement = W·W'`, `W,W'` commuting order-2. On
the arc-cube, complement flips **all** `C(n,2)` arcs (THM-584). Any factoring into two commuting
involutions is a partition of the arcs into two blocks `A ⊔ B` (flip `A`, flip `B`). For the two
factors to be tournament-natural they must be **S_n-invariant** blocks. But

> **S_n acts transitively on the `C(n,2)` arcs**, so the only invariant blocks are `∅` and everything.

Hence complement **cannot be equivariantly factored**: there is no `W₂·W₇` on the tournament side. The
Atkin-Lehner factorization `14 = 2·7` is *arithmetic* — two primes — and has no combinatorial image,
because the arcs form a single S_n-orbit with no "2-part" and "7-part." Verified computationally
(n=3,4,5): the only geometric candidate for a second involution, the reversal relabel `i→n−1−i`, lies
*in* S_n, so it acts as the identity on iso-classes — no second involution exists.

### 2. Complement isn't even a single Atkin-Lehner element (integrating the negatives)

Worse than "can't factor": complement doesn't sit in the Atkin-Lehner V₄ at all. THM-584 proves
complement **fixes the two SC classes** `{T,S}` and swaps only `{+,−}`. But `W(14)` acts **regularly**
(simply transitively, fixed-point-free) on the 4 cusps. A fixed-point-free action has no element
fixing two of four points, so **no single `W` can equal complement** on the 4 classes. (opus's
independent negative: complement swaps `+` (a would-be order-3 point) and `−` (order-2) — elements of
different order — so complement `≠ [−1]` and `≠` any group automorphism of `Z/6`.) The recurring
labeling "complement = Fricke `W₁₄`" is therefore a **heuristic in tension with THM-584**, never
reconciled — and the tension is structural: the tournament fold has fixed points (SC), the
Atkin-Lehner V₄ does not.

### 3. The n=4 weight-2 match is exact — and breaks at n=6

At the base level `n=4 ↔ X₀(14)`, the metagraph reproduces the **entire weight-2 modular dimension
count** of `Γ₀(14)`:

| `Γ₀(14)` weight-2 | dim | metagraph (THM-584) | dim | match |
|---|---|---|---|---|
| `M₂` (all forms) | 4 | `A000568(4)` (# classes = # cusps) | 4 | ✓ |
| `Eis₂` (Eisenstein bulk) | 3 | R-even block `(A+SC)/2` | 3 | ✓ |
| `S₂` (cusp form `f₁₄`) | 1 | R-odd block `(A−SC)/2` = genus = 1 NS-pair | 1 | ✓ |

`4 = 3 + 1` exactly, and complement is the `ι` whose ±eigenspaces are Eisenstein/cusp. Tempting — the
single n=4 NS-pair `{(3,1,1,1)↔(2,2,2,0)}` playing the cusp form. **But it is base-only numerology:**
the R-odd dimension is `0,1,2,22,140` for `n=3..7`, while `genus X₀(2p) = 0,0,1,2,2` for
`p=3,5,7,11,13`. The two agree at `0,1,2` (n=3,4,5) and **diverge at n=6** (`22` vs `≤2`). The R-odd
block grows combinatorially; the cusp space is a fixed 1-dimension. So R-odd `↔` cusp-form is not a
functor — it is exactly the "coincidence at n=14" the corpus repeatedly warns of, now pinned to its
breaking point.

---

## Where the genuine V₄ actually lives (not the metagraph)

The two order-2 generators do exist — just not on the tournament iso-classes:

- **On the labeled tiling cube:** `⟨complement σ, tile-flip⟩ = V₄` (HYP-3811/3814), where `tile-flip`
  flips all `m` tiles (base-path fixed) — a *fixed-point-free* involution distinct from arc-reversal.
  This is a labeled / base-path object; it collapses on iso-classes (where only complement survives).
- **On the arithmetic / runner side:** `W₂ =` the 2-adic descent (THM-580, `u=2t` circle cover, the
  `14→7` peel — degree 2) and `W₇ =` the apex-7 cyclotomic doublet (`4cos²(3π/7)`). These are the
  concrete second involutions, and they carry `f₁₄`'s modularity (the Hecke/Atkin-Lehner action).

So the V₄ is **split across the two domains**: complement is the only piece the metagraph offers, and
even it doesn't lie inside the Atkin-Lehner V₄. The arithmetic generators — where the cusp form's
modularity lives — are on the runner side.

---

## The unification: S_n-transitivity is the single root

The transitivity of S_n on the arcs is one cause of a whole cluster of proved severances:

1. **complement irreducible** ⇒ no `W₂·W₇`, the `2·7` factorization is invisible (this note);
2. **`CV(H)` bounded, "no vanishing fiber"** ⇒ the variance-transfer is miscalibrated (HYP-3554);
3. **"the testbed models the bulk, not the cusp"** (klein-S4);
4. **the odd-graph → cusp transfer is provably blind** — the difference tournament is rotation-
   invariant, so `{1..13}` (tight) and `{1..14}∖{6}` (covering) are identical yet have different `L`
   (S90, HYP-6555).

All four are the same statement: the tournament's clean transitive symmetry is *too smooth* to resolve
the arithmetic structure (the prime `2·7` split, the absolute position, the cusp). The metagraph is the
Eisenstein bulk made combinatorial; the cusp form is arithmetic and lives off it. This is why the
n=4 match is a shadow, not a map — the shadow is the whole modular *dimension count* (`4=3+1`), which is
symmetric/bulk data; the map would need the *Hecke action*, which is the `2·7` structure S_n forbids.

---

## Coda

The Atkin-Lehner V₄ does not descend to the tournament metagraph. Exploring the metagraph side of the
thread returns a clean structural verdict rather than a bridge: (i) one iso-class involution
(complement, THM-280/584); (ii) irreducible under S_n, so the `2·7` factorization cannot form; (iii)
fixed-point-incompatible with a single Atkin-Lehner element; (iv) the exact `4=3+1` weight-2 match at
n=4 is base-only, breaking at n=6. The genuine second generators are the 2-adic descent and the apex-7
doublet — on the runner side, carrying the cusp form's modularity — and the labeled tile-flip. The one
survivor of S90 is thus explained *and* bounded: the V₄ is real on both sides, but the tournament's
S_n-transitivity is precisely the obstruction to it descending. If the correspondence is ever to be
built, it will not be by finding a second tournament involution (there is none) — it will be by a
global object that breaks the arc-transitivity, exactly as opus concluded. The metagraph gives the
bulk; the cusp is arithmetic; and now we know, structurally, why.

---

*Cross-links: base — THM-280 (grid-reflection=complement), THM-584 (complement=antipodal, R-even/odd),
klein-S10/HYP-3586 (cusps=n=4 classes). Honest negatives integrated — opus's `modular-pushes` (complement
≠[−1]; floor≠L-value; sym²=GL(3) structural-only), `the-hecke-dictionary-of-f14`, S90/HYP-6555
(odd-graph blind), NEG-6 (metagraph V₄ Burnside 1,1,2,7,28,240,3440 ≠ {1,3,21,55}). Runner-side V₄ —
THM-580 (W₂=2-adic descent), apex-7 (W₇). Labeled V₄ — HYP-3811/3814 (⟨σ,tile-flip⟩). The unified root
— HYP-3554 (CV bounded), klein-S4 (bulk not cusp). My computation: HYP-6565,
`04-computation/metagraph_v4_atkinlehner_macmini_S91.py`.*
