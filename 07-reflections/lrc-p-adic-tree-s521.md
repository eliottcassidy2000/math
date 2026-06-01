# The p-adic tree: organizing the bounded-denominator program and the first-even bridge (S521)

*claudebox-2026-06-01-S521. Integrating the p-adic tree into the LRC themes
(Galois-Weil equidistribution, the resonances, the first-even bridge, the
bounded-denominator program). The p-adic tree turns out to organize all of them:
the lonely set's p-adic dimension separates tight from non-tight, doublings are
tree edges while sum-relations cross the tree, and the cyclotomic Z_p-tower is the
Galois-side tree.*

## Two p-adic trees

- **Time tree (the p-ary subdivision of `[0,1)`).** Level `k` nodes are the
  intervals `[a/p^k, (a+1)/p^k]`; each has `p` children. A real lonely time is an
  infinite path (a p-adic integer); a rational lonely time `a/p^k` is a finite node;
  a non-`p`-power rational (like `1/n`) is an infinite NON-terminating path.
- **Speed tree (dilation).** `||v(pt)|| = ||(pv)t||`, so multiplying a speed by `p`
  is the p-adic descent on speeds; `v_2(v_i)` (2-adic valuation) is the depth in the
  2-adic speed tree. Primitivity (gcd 1) = the speeds do not all share a child
  position.

## The lonely set's p-adic dimension: tight vs non-tight

Mark a time-tree node "lonely-live" if its interval contains a lonely time. The
**live subtree** encodes the lonely set. Computed (p=2):

- **Non-tight sets** (`{1,2,4,7}`, `{2,3,5,7,11}`): live nodes GROW with the level
  (2,2,4,4,6,12 ...), branching toward `p` — the lonely set is **fat** (positive
  measure), p-adic **dimension 1**, a full subtree.
- **Tight/extremal sets** (`{1,2,3,4}`, `{1,3,4,5,9}`): NO finite node is fully
  lonely-live — the lonely point is `t = 1/n`, a non-dyadic rational, hence an
  **infinite path** in the 2-ary tree captured by no finite node. p-adic
  **dimension 0**: the lonely set is a single p-adic limit.

> **The tight extremizers are exactly the speed sets whose lonely set is a single
> infinite path (a p-adic limit, `t=1/n`) rather than a fat subtree.** Tightness =
> p-adic dimension 0 of the lonely set.

This is the p-adic-tree face of "tight = measure-zero lonely set at the boundary
`t=1/n`," and it explains why finite (dyadic) searches miss the tight cases — they
live at an infinite, non-terminating address.

## Doublings are tree edges; sum-relations cross the tree

On the speed tree, the **doubling relation `v_j = 2 v_i`** is a 2-adic
**parent-child edge** — and the Galois-Weil study found doublings are *benign*
(positive/zero equidistribution bias). The **sum-relation `v_i + v_j = v_k`** is a
**cross-tree additive** relation — and it is the *obstruction* (negative
equidistribution bias). Computed: `{1,2,4,7}` has doublings `(1,2),(2,4)` but NO
sum-relation, so it is non-resonant (lonely by the small prime `q=11`); the
resonant sets (`{1,2,3,4}`, `{2,3,5,7,11}`) all carry sum-relations
(`1+2=3`, `2+3=5`, ...). So in the p-adic picture:

> tree edges (multiplicative doublings) are harmless; the resonances that suppress
> equidistribution are the cross-tree **additive triples** `v_i+v_j=v_k`.

## The first-even bridge is the 2-adic tree

`n = 2 * odd` (n=14, 18, ...) splits the speeds by 2-adic valuation: the odd
speeds (depth 0) and the even speeds (depth >= 1). This is exactly the project's
**first-even-bridge / gate / double-gate / dyadic** structure (HYP-1905, HYP-1952,
the n=14/18 work): the "gates" are levels of the 2-adic tree, and "double-gate" =
descending two levels (`v_2 = 2`). The 2-adic tree on speeds is the rigorous home
of that dyadic ladder. A blocker (speed `≡ 0 mod n`) is deep in BOTH the 2-adic and
the odd-prime trees (divisible by `2` and by the odd part).

## Local-global over the p-adic trees (CRT)

A denominator `q = prod_p p^{a_p}`; by CRT the residue `v_i mod q` is the tuple of
residues at each prime's tree level `a_p`. The safe-box (an interval mod `q`) does
NOT factor as a CRT product, so the primes are coupled — but the DANGER (small
residue near 0) localizes. The bounded-denominator program reads on the trees as:

1. **base `q = n`** = a specific (often non-`p`-power, hence infinite-path) address
   — the regular polygon / cyclotomic node;
2. **tail (deep tree levels / large `q`)** = Galois-Weil equidistribution forces a
   live node;
3. **middle** = a finite search over shallow tree nodes;
4. the **resonances** (cross-tree additive triples) are the only obstruction to (2),
   and they are finite.

*[Whether descending a SINGLE prime's tree suffices, or the CRT product across
primes is essential, plus the per-prime local densities and the cyclotomic
`Z_p`-tower (Iwasawa) refinement, are under separate computation.]*

## Sharpened by computation: ONE density (no Euler product), CRT only for composite n

A descent study down the trees (primes 2,3,5,7,11,13, levels to `q ~ thousands`)
gives three sharp facts — including a correction:

- **`f(p^k) -> mu`, the true lonely measure, NOT the naive volume `V=(1-2/n)^m`,**
  and the limit is **independent of the prime** (`d_2 = d_3 = ... = mu`). The gap
  `V - mu` is exactly the correlation among the runners' safe-box constraints (the
  resonances again). Convergence is non-monotone (oscillates and damps).
- **`mu = d_p` (a single density), NOT the Euler product `prod_p d_p`** (which
  collapses to ~0). *This corrects an earlier guess of mine.* The reason is
  structural: all runners share the SAME time `t`, so loneliness is governed by one
  density, not a product over independent primes. **So descending a SINGLE prime's
  tree already sees the full lonely density** — one tree suffices to detect
  loneliness whenever `mu > 0`.
- **CRT (multiple primes) is essential precisely when `n` is composite.** For prime
  `n` a single prime power already works; but the tight composite-`n` extremizers
  (`{1,2,3,4,5}`, `{1,3,4,5,9}`, `n=6`) have NO lonely prime-power time and are
  lonely only at `q = n = 2*3` — the local-global alignment of both prime factors.
- **2-adic law (n even):** if ALL speeds are odd, the 2-adic tree alone suffices
  (`a = 2^{k-1}` sends every odd speed to `~q/2`, the safe point); MIXED 2-adic
  valuations obstruct the 2-adic tree and force the odd prime. (This is the
  first-even bridge's parity split, made exact.)

**The clean reduction this yields.** Since `mu > 0` implies a positive-measure set
of lonely times (caught by descending one tree), the entire difficulty is the
`mu = 0` cases — which are exactly the **tight extremizers** (p-adic dimension 0,
lonely point `t = 1/n`). So:

> **LRC  <=>  every `mu = 0` speed set is lonely at the boundary `q = n`**
> (the regular-polygon / cyclotomic node). The `mu > 0` sets are automatically
> lonely by single-prime tree descent.

This is the sharpest distillation of the whole program: non-tight = positive p-adic
density (one tree), tight = the measure-zero boundary at `q=n`; LRC is exactly the
boundary statement, the regular polygon / Thm B locus.

## The Galois side: the cyclotomic Z_p-tower (Iwasawa)

The Galois groups `Gal(Q(zeta_{p^k})/Q) = (Z/p^k)*` form the cyclotomic
`Z_p`-extension — the Galois-side p-adic tree. The Galois orbit `O_{p^k}` refines as
`k` grows (descending the tower); equidistribution down the tower is the p-adic
refinement of the Weil bound. This is the Iwasawa-theoretic frame for the
Galois-Weil strategy: the lonely-fraction `f(p^k)` as a function on the tower, its
limit the local density `d_p`.

## Assessment

The p-adic tree is the natural organizing structure for the whole S521 program:
it separates tight (dimension-0 p-adic limit at `1/n`) from non-tight (fat,
dimension 1); it places doublings as benign tree edges and the sum-relation
resonances as the cross-tree obstruction; it identifies the first-even bridge as
the 2-adic subtree (rigorizing the project's gate/dyadic thread); and it frames the
Galois-Weil reduction as descent down the cyclotomic `Z_p`-tower with a finite,
local obstruction. It does not prove LRC, but it unifies the additive (time tree),
multiplicative (speed/Galois tree), and resonance (cross-tree) structures into one
tree picture, with the bounded-denominator program as a tree search whose only
obstacle is the finite set of cross-tree additive triples.

## Leads

- **p-adic dimension as a tightness invariant.** Formalize the lonely set's p-adic
  (box-counting) dimension; prove tight <=> dimension 0 <=> lonely point is the
  infinite path `1/n`; relate the dimension to `M(v) - 1/n`.
- **2-adic gate ladder = the project's first-even bridge.** Recast the n=14/18
  gate/double-gate results (HYP-1952 etc.) as statements about the 2-adic speed
  tree; the doubling-benign / sum-relation-obstruction dichotomy should sharpen the
  even-case difficulty.
- **Iwasawa descent (corrected).** `d_p = lim_k f(p^k) = mu` for EVERY prime (one
  density, not an Euler product — that guess was wrong, since the runners share one
  time `t`). The useful content: one prime's tower already computes `mu`; so the
  proof obligation is purely the `mu = 0` (tight) boundary cases at `q = n`, where
  the cyclotomic field `Q(zeta_n)` and the regular-polygon / Thm B structure live.
  Pursue: prove `mu = 0  =>  lonely at q = n` (the tight extremizers are
  boundary-lonely) — equivalently that no `mu=0` set is an empty-lonely
  counterexample.
