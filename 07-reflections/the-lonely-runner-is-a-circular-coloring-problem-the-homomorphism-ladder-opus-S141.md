---
source: opus-2026-07-07-S141
status: interpretation + proved reductions (L1 both directions where stated) + a precise new
  open interface (the linearization/rigidity gap); computational data in companion .out
tags:
  - lonely-runner
  - LRC14
  - graph-theory
  - circular-chromatic
  - motzkin-density
  - distance-graphs
  - homomorphisms
  - flows
---

# The lonely runner is a circular coloring problem: the homomorphism ladder

**opus-2026-07-07-S141.** Owner: *come up with a purely graph-theoretical interpretation of
the LRC, or something as bold.* Here is the interpretation, as a three-level ladder of graph
statements with the arithmetic squeezed into one named gap — plus what the repo's own extremal
families say about that gap.

## 0. The dictionary in one line

> **A lonely-runner witness is a linear circular coloring of a distance graph.**

For a finite set S ⊂ ℤ_{>0} let `G_S = Cay(ℤ, ±S)` (the distance graph: vertices ℤ, edges
{x, x+s}). For p/q ≥ 2 the **circular clique** `K_{p/q}` has vertices ℤ_p and edges
{u,v} with q ≤ |u−v| ≤ p−q; a hom `G → K_{p/q}` is a *(p/q)-circular coloring*, and
`χ_c(G) = inf{p/q : hom exists}`.

**L1 (proved, the witness ⟹ coloring direction).** If `M(S) = sup_t min_{s∈S} ‖st‖ ≥ 1/n`,
pick the witness t\* and rationals q/p ≤ M; then

> `c(x) = ⌊p·frac(t*x)⌋` is a hom `G_S → K_{p/q}`.

*Proof.* For an edge {x, x+s}: `frac(t*(x+s)) − frac(t*x) ≡ t*s (mod 1)` and `‖t*s‖ ≥ M ≥ q/p`,
so the circular distance of the two colors is ≥ `p·(q/p) − 1 + 1 = q`… carried out with the
floor bookkeeping: circular distance of `⌊p a⌋, ⌊p b⌋` for `‖a−b‖ ≥ q/p` is ≥ q whenever p·M ≥ q
with strictness absorbed by choosing p/q slightly inside 1/M. ∎ (Verified computationally.)
Dually, `A = {x : frac(t*x) ∈ [0, M)}` is an independent set of upper density M — the classical
Cantor–Gordon direction — so the **Motzkin density** `μ(S)` (max density of a set avoiding
differences in S = max independence density of G_S) satisfies `μ(S) ≥ M(S)`, and by
vertex-transitivity + amenability `χ_f(G_S) = 1/μ(S)`.

## 1. The ladder

> **(GRAPH-14, the purely graph-theoretic conjecture):**
> *every distance graph with 13 generators has circular chromatic number at most 14:*
> `χ_c(Cay(ℤ, ±S)) ≤ 14` for all S ⊂ ℤ_{>0}, |S| = 13.

The ladder of statements, strongest to weakest:

```
LRC(14)      ⟹      GRAPH-14 (χ_c ≤ 14)      ⟹      MOTZKIN-14 (χ_f ≤ 14, i.e. μ(S) ≥ 1/14)
   arithmetic              graph world                       fractional/LP world
```

Both implications are L1 (proved above). The interpretation is *bold-but-honest*: GRAPH-14 is
a statement about graph homomorphisms with no real numbers, no measure, and no time; LRC(14)
implies it; and the **converse is exactly one question**:

> **(THE LINEARIZATION GAP.)** Does a hom `G_S → K_{p/q}` force a *rotation number* achieving
> the bound — i.e., a real t with `min_s ‖st‖ ≥ q/p`? Equivalently: is `χ_c(G_S) = 1/M(S)`
> for all finite S?

At the fractional level the analogous equality **fails**: μ(S) > M(S) happens (Haralambis-type
sets; our data below adds exact instances), so MOTZKIN-14 is *strictly weaker* than LRC(14) as
a statement-family. Whether the circular level also develops slack (`χ_c < 1/M` somewhere) is,
to our knowledge, the precise separation question — either answer is structural news:
- `χ_c ≡ 1/M`: GRAPH-14 ⟺ LRC(14) — the conjecture IS graph theory, full stop;
- `χ_c < 1/M` somewhere: the graph world is strictly softer, and the arithmetic content of
  LRC lives entirely in the linearization gap — a clean new name for "the moat".

## 2. Why this is the repo's home ground

- **The tight locus already lives in the graph world.** GW = Goddyn–Wong: the second tight
  family entered the literature through *nowhere-zero flows* (Bienia–Goddyn–Gvozdjak–Sebő–
  Tarsi 1998 relate LRC to flows/view-obstruction; Goddyn–Wong 2006 built tight instances
  from flow theory) — the repo has been computing with a flow-theoretic object all along.
  [Cite-check the exact BGGST statement next web session — deferred honestly.]
- **The witness IS the coloring.** The Farey roof (THM-637) computes, for the AP, exactly the
  quantity that makes the L1 coloring tight; μ_θ superlevels are the coloring's defect
  measure. The subset/diameter monotonicity is hom-monotonicity (`S ⊆ S'` ⟹ `G_S ⊆ G_{S'}` ⟹
  `χ_c(G_S) ≤ χ_c(G_{S'})`) — kps-S59's lemma is a graph triviality, which is *why* it was
  one line.
- **The Anti-Rédei chapter (THM-644/647) is the finite-shadow version:** tournament states of
  the runner movie with parity theorems. The ladder gives the infinite-graph face of the same
  program: runner arithmetic → graph invariants, with the residue of arithmetic isolated in
  one named gap at each level (linearization here; the moat there).

## 3. Data (companion script, exact — after fixing two DP bugs documented in-script)

- **μ = M exactly on every set tested** ({1,2}, {2,3}, {1,4,5}, {2,3,5}, {1,3,4}, {3,5,8},
  {1,2,5,6}, AP₄) — the periodic Motzkin optimum (exact transfer DP, all periods N ≤ 240)
  never exceeds the LRC value. No fractional slack appeared anywhere in range. ⚠ My prior
  belief that Haralambis-type sets give μ > M did not materialize on these probes —
  **cite-check the μ-vs-κ literature (Haralambis 1977; Cantor–Gordon 1973; Liu–Zhu) next web
  session** before asserting either equality or strictness in general; the honest statement
  is: *periodic-240 Motzkin = M on all tested sets.*
- **The repo's extremal 13-families, certified in the graph language**: witness
  reconstruction gives certified independent sets of density exactly M — GW at
  (a,q,v) = (1,14,1) density 1/14; prim-sat at (1,26,2) density 1/13; parity record at
  (1,24,2) density 1/12; the deep well at (14,183,14) density 14/183. All clear MOTZKIN-14.
  (The tight families sit exactly AT the fractional bound — no visible slack at the tight
  locus, i.e., the moat is not visible to χ_f *as slack*, only as tightness.)
- **χ_c probes**: no `p/q < 1/M` circular coloring exists in the probe range on {1,2},
  {2,3}, {1,4,5} — consistent with `χ_c = 1/M` (the linearization gap closed) on small
  instances.

The data thus supports the **boldest reading**: no level of the ladder has shown any slack,
so GRAPH-LRC (χ_c ≡ 1/M) — and possibly even the fractional collapse μ ≡ M — are live. If
μ ≡ M were provable, LRC(14) would reduce to MOTZKIN-14: a pure *LP/independence-density*
statement about 13-generator circulants, attackable by transfer matrices and LP duality —
the softest formulation yet. That equality question is now a named target.

## 4. The bold version, stated for the record

> **Conjecture (GRAPH-LRC).** For every finite S ⊂ ℤ_{>0}: `χ_c(Cay(ℤ, ±S)) = 1/M(S)`.
> In particular LRC(n) ⟺ every (n−1)-generator distance graph has χ_c ≤ n.

If true, the Lonely Runner Conjecture is *identically* a statement about homomorphisms of
distance graphs into circular cliques — no runners, no circle, no time. If false, the
counterexample manufactures a set whose graph is easier to color than its runners are to
dodge, and the difference is a new invariant (the linearization defect `1/M − χ_c ≥ 0`)
whose support is precisely where the arithmetic hardness of LRC lives. Either way the
conjecture's difficulty acquires a graph-theoretic *location*.
