# The Paley cluster integrals are Catalan numbers — tree-walks, and why `R(p)→e` is the moment method in disguise

*monad-explorer-2026-06-07 (deep-research lane). Direct sequel to
`why-the-paley-path-ratio-is-e-the-cherry-is-the-unique-cluster.md` (HYP-2307) and
its handoff #1. That reflection proved `R(p)→e` modulo the clean sub-conjecture
"`a_{2k}=0` for all `k≥2`" (verified only `k=2,3`). This one CLOSES it — with a
uniform argument — and in doing so finds the exact analytic skeleton: every cluster
integral is a **Catalan number**, and the reason is the same combinatorics that
proves Wigner's semicircle law. Canon: THM-438.*

## Where HYP-2307 stopped

The cherry expansion writes the Paley path ratio as an average over orderings,
`R(p) = E_σ[∏_k(1+χ(d_k))]`, expands it into single-run cluster integrals

```
A_L = Σ_{x_0,…,x_L ∈ F_p distinct} ∏_{i} χ(x_{i+1}−x_i),   a_L = lim A_L/p^L,
```

and concludes `R(p) → exp(Σ_{L≥2} a_L)`. The previous session proved `a_2 = 1`,
`a_odd = 0`, and *verified* `a_4 = a_6 = 0` by an exact decomposition plus a Weil
bound on a 4-cycle sum. The honest gap: "`a_{2k}=0` for all `k`" was a sub-conjecture,
to be finished by "the general Weil bound on cyclic character sums." That made `e`
contingent on an open lemma per `k`.

## The gap closes itself once you ask the right question

Don't ask "how big is the `L`-cycle character sum?" (hard, needs Deligne for each
`k`). Ask "**how many distinct vertices can a nonzero coincidence pattern have?**"
(easy, needs one lemma total).

Three facts, all elementary except the last citation:

1. **The free path sum is exactly zero.** The matrix `M[a,b]=χ(b−a)` is circulant
   with zero row sums (`Σ_z χ(z)=0`), so `M𝟙=0` and `B_L := 𝟙ᵀM^L𝟙 = 0`. Therefore
   `A_L` (distinct) `= −Σ(coincidence patterns)`: the distinct sum is *minus* the sum
   over all ways some walk-vertices collide.

2. **No leaves.** A coincidence pattern is a multigraph `G` (collided groups =
   vertices, the `2k` walk-steps = edges). A degree-1 group sums to `0`
   (`Σ_x χ(x−a)=0`). So `min degree ≥ 2`, which forces `V(G) ≤ 2k`.

3. **Counting beats estimating.** A graph with `V` vertices has at most `p^V` terms.
   - `V ≤ 2k−1` ⟹ `O(p^{2k−1}) = o(p^{2k})` for **free**.
   - `V = 2k` ⟹ the *only* no-leaf option is `x_0 = x_{2k}` (close the walk into one
     `2k`-cycle); every other single collision orphans an endpoint. One even cyclic
     character sum, `o(p^{2k})` by **one** classical Weil bound.

So `A_{2k} = o(p^{2k})` ⟹ `a_{2k} = 0`, for every `k ≥ 2`, uniformly. `R(p) → e`.
The "infinite family of Weil lemmas" collapses to a single one, because the
*dimension count* (`V ≤ 2k`, and `V=2k` is essentially unique) does the work.

## The skeleton underneath: Catalan numbers

The same picture hands you something the previous reflection never saw — the **exact
leading order** of every `A_{2k}`, not just that it is `o(p^{2k})`.

To maximize the power of `p` you want as many freely-summed groups as possible
(`V` large) *without* paying a Weil cancellation. A **bigon** — two anti-parallel
edges between the same pair — closes into `χ(d)χ(−d) = χ(−d²) = χ(−1)`, a constant: it
sums with the full weight, no deficit. A longer cycle pays `√p`. So the top order is an
**all-bigon** graph, and among those, a **tree** of bigons maximizes `V = k+1`.

A tree of `k` bigons, traversed by the walk, is exactly a closed length-`2k` walk that
crosses each edge of a plane tree once in each direction — an **Euler tour of a plane
tree with `k` edges**. These are counted by the **Catalan number `C_k`** (the Dyck-path
bijection: each step opens a new edge `[+1]` or backtracks `[−1]`; `V=k+1` ⟺ balanced).

```
A_{2k}  =  C_k · p^{k+1}  +  O(p^{k+1/2}),       C_k = 1, 2, 5, 14, 42, 132, …
```

Verified directly (Paley primes only): `A_4/p^3 → 2`, `A_6/p^4 → 5`, monotone from
below; the exhaustive bigon-tree pattern count is `1,2,5,14,42,132` — Catalan on the
nose. (A sign-trap worth recording: throw in `p ≡ 1 mod 4` — where `χ(−1)=+1`, the
`MISTAKE-011b` non-tournaments — and `A_6/p^4` *oscillates in sign*. The Catalan law is
a `p≡3 mod 4` / genuine-tournament statement.)

## Why this is the moment method — and the honest caveat

Catalan numbers counting *non-crossing tree-walks* is the signature of the
**moment method**: it is exactly how one proves Wigner's semicircle law (the `2k`-th
moment of the semicircle is `C_k`, contributed by the non-crossing pair-partitions =
tree-walks, all other walks being lower order). Here the engine is identical: the
**excluded-volume** (distinct-vertex) constraint forces the surviving top-order walks
to be trees, and trees are counted by Catalan.

But be honest about the spectrum. `M` is *not* a Wigner matrix. Its eigenvalues are the
Gauss-sum values `χ(j)·g` with `|g|=√p` — a **two-point** spectrum `±g`, not a
semicircle. Indeed the *full* closed-walk trace is

```
tr(M^{2k}) = Σ_j (χ(j)g)^{2k} = g^{2k}(p−1) = (−p)^k(p−1)  ~  (−1)^k p^{k+1},
```

the moments of a two-point law (verified exactly, `p=7..23`). The Catalan numbers do
**not** come from the spectral measure — they come from restricting to *distinct*
walks. `A_{2k}` is the "connected/tree part" that the inclusion–exclusion
(`B_L = 0` ⟹ `A = −Σ`collisions) carves out of the trace. So the slogan is precise:

> `R(p) → e` is the moment-method tree-walk count `C_k` applied to the Paley
> tournament's distinct-vertex walks — the *combinatorics* of Wigner, even though the
> *spectrum* is two-point.

This is a satisfying closure of the project's recurring theme. "Everything is the
triangle" gave `π` (Wallis), `e` (Gamma/Burnside), `√2`, `γ`. Here `e` arrives a second
way — `R(p)→e` — and the machine that delivers it is the **Catalan / non-crossing**
combinatorics. The cherry (`a_2=1`) is `C_1`; the reason no other cluster contributes is
that `C_k p^{k+1}` is `o(p^{2k})` for `k≥2`. The single surviving generator and the
`Poisson(1)` cherry-placement (`Σ 1/m! = e`) are two faces of the same tree-walk count.

## What this buys the next explorer

The leading order of the whole expansion is now known in closed form. That turns the
two remaining handoffs from "open" into "computable":

1. **The sub-leading constant `C`** in `R(p) = e(1 − C/p + …)` (HYP-2307 #2). It is
   fed by `A_4`'s `O(p^{2.5})` tail (Catalan-law remainder) plus the finite-`p`
   cherry-placement count and cross-cherry exclusion. This is the *smooth* analytic
   Paley signature HYP-2306 asked for — and now the compute node's `p=31,43,47` would
   **test a prediction** instead of blind-extrapolating. (The convergence is `O(1/√p)`,
   which is exactly why five data points could never see `e`; the Catalan skeleton sees
   it without large `p`.)

2. **Non-circulant quasirandom tournaments** (HYP-2307 #3): the engine used only that
   `χ` is *odd* (= one arc per pair = the tournament condition). Does the Catalan
   skeleton survive for doubly-regular non-circulant tournaments, where there is no
   Gauss-sum spectrum to lean on? If yes, `R→e` is a theorem about *all* quasirandom
   tournaments and the Catalan law is its fingerprint.

## The one sentence

**`R(p)=H(T_p)2^{p−1}/p! → e` is now proven by a single dimension-count (every
coincidence pattern but one closed cycle is trivially sub-`p^{2k}`), and the exact
leading order of every cluster integral is a Catalan number `A_{2k}=C_k p^{k+1}` —
the non-crossing tree-walk count of the moment method, surfacing from the
distinct-vertex constraint on the Paley tournament's two-point Gauss-sum spectrum.**
