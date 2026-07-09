---
source: opus-2026-07-09-S182
status: EXTENDED/REDIRECTED the S181 forward lead. The density-floor resonance sum R (L=(6/7)^13+R) is a
  theta-sum over the resonance lattice whose MINIMAL vectors are the height-3 SCHUR TRIPLES a+b=c (norm 3),
  NOT the height-4 additive-energy relations a+b=c+d (norm 4). Confirmed: the density-floor deficit
  (6/7)^13 - L tracks the ORDER-3 additive energy E3=#{(a,b,c): a+b=c} (corr +0.79) better than the
  order-2 energy/doubling E2 (corr +0.65). CLEAN SEPARATOR: the sum-free odd set {1,3,..,25}=2*{1..13}-1
  has MAXIMAL E2=1469 (=the tight AP) but E3=0 (odd+odd=even, no Schur triple) and is LOOSE (deficit 0.019)
  -- resolving the S181 translation puzzle. The STRUCTURAL reason is a SYMMETRY MATCH: loneliness is
  dilation- but NOT translation-invariant; E3 has exactly that invariance group; E2 is ALSO
  translation-invariant, so E2 cannot tell the tight AP from its loose translate. The AP {1..13} uniquely
  maximizes E3 (=78) among nonzero speeds. REDIRECTS the forward lead: the a-priori density floor is a
  SCHUR-TRIPLE / SUM-FREE structure theorem (E3), NOT a Freiman/doubling (BSG, 3k-4) theorem (E2).
tags:
  - lrc14
  - schur-triples
  - additive-energy
  - sum-free
  - resonance-lattice
  - symmetry-matching
  - density-floor
---

# The resonance sum is Schur triples, not doubling — a symmetry match

**opus-2026-07-09-S182.** Owner: extend the S181 forward lead (make the density floor a-priori via additive
combinatorics), be inspired to explore.  The lead pointed at Freiman/doubling tools (BSG, 3k−4).  Following
the resonance sum to its leading order redirected it to a *different, sharper* invariant — and the
redirection has a clean symmetry reason that also resolves the S181 translation puzzle.

## The resonance sum is a Schur-triple count

The lonely measure decomposes `L(S) = (6/7)¹³ + R(S)`, with `R(S) = Σ_{n:n·v=0, n≠0} 𝒲̂(n)` a theta-sum
over the resonance lattice `Λ={n∈ℤ¹³ : Σnᵢvᵢ=0}`.  A theta-sum is dominated by the lattice's MINIMAL
vectors, and the arc-coefficient product `∏b₀(nᵢ)` has one factor per nonzero `nᵢ` with `|b₀|<1` — so
FEWER factors ⇒ larger term.  The minimal relations:

- **height 3**: `a+b=c`, i.e. `n=e_a+e_b−e_c` (`‖n‖₁=3`) — a **Schur triple**.  Three factors.
- **height 4**: `a+b=c+d`, i.e. `n=e_a+e_b−e_c−e_d` (`‖n‖₁=4`) — the **additive energy** relation. Four.

So the theta-sum is led by the Schur triples (norm 3 = the relation lattice's kissing vectors, mac-mini-S25),
not the additive-energy relations (norm 4).  Measured (`lrc14_resonance_is_schur_triples`), across
APs/near-APs/GAPs/Sidon/Fibonacci/random, the density-floor deficit `(6/7)¹³−L` correlates with the
order-3 energy `E₃=#{a+b=c}` at **+0.79**, beating the order-2 energy `E₂` at **+0.65**.

## The clean separator and the symmetry match

The decisive row is the **sum-free odd set `{1,3,…,25}` = `2·{1..13}−1`**:

| set | `E₂` (energy) | `E₃` (Schur) | `L` | verdict |
|---|---|---|---|---|
| AP `{1..13}` | 1469 (max) | **78 (max)** | 0 | tight |
| sum-free `{1,3,…,25}` | **1469 (max)** | **0** | 0.116 | LOOSE |

Same maximal `E₂`, opposite `E₃`, opposite looseness.  `E₂` cannot see the difference; `E₃` sees it exactly.
The reason is an invariance match:

> **`M(S)` is dilation-invariant (`M(cS)=M(S)`, THM-522) but NOT translation-invariant.  `E₃` (`a+b=c`) has
> exactly that group: dilation-invariant, but broken by `x↦x+t` (`a+b=c` becomes `a+b+t=c`).  `E₂` (`a+b=c+d`)
> is ALSO translation-invariant, so it is constant on the translation orbit `{1..13} → {1,3,…,25}` — it
> literally cannot distinguish the tight AP from its loose translate.**  A quantity can only be governed by
> an invariant with the SAME symmetry group.  Loneliness's group is dilation-only; so is `E₃`'s; `E₂`'s is
> strictly larger.  Therefore `E₃`, not `E₂`, is the additive-combinatorial governor of loneliness.

This resolves the S181 puzzle (why the max-energy translate `{1,3,…,25}` is loose): translation preserves
`E₂` but annihilates the Schur triples (`odd+odd=even ∉` odds), and loneliness follows `E₃`.  It also
sharpens S180: the AP is the tight extremal because it maximizes `E₃` (`=78` among nonzero speeds — the
unrestricted maximizer `{0,…,12}` scores 91 via free triples `a+0=a`, but `0` is not a legal speed, exactly
as `0` would trivially cover the origin).

## The redirect of the forward lead

The S181 lead was "BSG (energy ⇒ small-doubling subset) → Freiman 3k−4 (small doubling ⇒ short AP)."  Those
are `E₂`/doubling tools.  Since the density floor is governed by `E₃`, the a-priori route is instead a
**Schur-triple / sum-free structure theorem**:

1. **Extremal:** `E₃(S) ≤ E₃(interval) = C(k,2)` among nonzero `k`-sets, equality iff `S` is a dilated
   interval `c·{1..k}`.  (The interval maximizes Schur triples — the small elements' sums land back in the
   set; this is a clean, classical-flavored extremal problem, the additive-triple analogue of Freiman.)
2. **Theta-sum bound:** `|R| ≤ f(E₃)` with `f(E₃(AP))=(6/7)¹³` and `f` increasing ⇒ `|R|<(6/7)¹³` for every
   non-interval set ⇒ loose (density floor).

The `+0.79` (not `1.0`) correlation says step 2 is not a pure scalar law — the Schur triples' DIMENSION /
coherence still modulates (the near-AP's `E₃=71` aligned in 1-D gives deficit `0.102`; the 2-D GAP's
`E₃=66` spread over two directions gives only `0.044`).  So the honest object is the **Schur-triple
sublattice geometry** (count + alignment), with `E₃` the leading scalar — the order-3 refinement of S181's
"resonance geometry, not a scalar."  BSG/3k−4 remain relevant only insofar as small doubling forces the
Schur triples into a low-dimensional progression; the primary bound is on `E₃`.

## Ledger

- The density-floor resonance sum `R` is led by height-3 SCHUR TRIPLES `a+b=c` (the relation lattice's
  minimal vectors), so the deficit `(6/7)¹³−L` tracks `E₃=#{a+b=c}` (corr `+0.79`) over `E₂`/doubling
  (`+0.65`). Separator: sum-free `{1,3,…,25}` (max `E₂`, `E₃=0`, loose).
- STRUCTURAL reason = symmetry match: loneliness and `E₃` share the dilation-only invariance group; `E₂` is
  additionally translation-invariant, so it cannot govern loneliness (can't separate AP from its translate).
  Resolves S181's translation puzzle; sharpens S180 (AP maximizes `E₃=78` among nonzero speeds).
- REDIRECTS the S181/HYP-5682 lead: the a-priori density floor is a Schur-triple/sum-free extremal theorem
  (`E₃(S)≤E₃(interval)`, equality iff dilated interval) + a theta-sum bound, NOT a Freiman/doubling (BSG,
  3k−4) theorem. Dimension/coherence of the Schur sublattice is the second-order correction.
- Files: `lrc14_resonance_sum_vs_doubling_opus_S182` (+out), `lrc14_resonance_is_schur_triples_opus_S182`
  (+out), `lrc14_maxE3_extremal_opus_S182` (out). -> opus-S181 (resonance geometry)/S180 (AP extremal)/S179,
  THM-515 (theta-sum)/522 (invariance), mac-mini-S25 (relation-lattice kissing = additive triples),
  opus-2026-06-29 (AP maximizes higher energy E₃,E₄), HYP-5682 (redirected).
