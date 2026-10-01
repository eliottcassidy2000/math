# Trying to prove Collatz connectivity from rigidity: four routes and where each stops, the transfer barrier (3n−1 and 5n+1 are rigid too), two-place density of the trivial component (every finite view of every integer occurs in it), the quantifier exchange that remains, and the minimal counterexample's two-sided trap

**Session:** opus, `collatz-functional-uniqueness-20261001` (S15 series, eighteenth note), 2026-10-01.
**Owner's directive:** "try proving collatz connectivity using the rigidity theorem".
**Inherits:**

* The [seventeenth note](collatz_functional_uniqueness_20261001.md): Theorem R, Theorem R_a, Proposition M (no periodic invariant), Proposition L (cycle half first-order, divergence half not), Theorem C (`Collatz ⟺ (N,T) ≅ Gamma_1`).
* THM-4517 (`01-canon/theorems/THM-4517-no-root-uniform-positive-density-of-collatz-basins.md`): the sheet cap, i.e. a sheet-blind cycle-uniform basin bound is `≤ 1/3`, and the trunk-entry points `a_j = (4^j − 1)/3`.
* HYP-9165: the 3n−1 basin densities `0.327, 0.325, 0.348`.
* The synthesis's requirement that a mechanism must see the 3n−1 sheet sign and the `log 3 < 2 log 2` drift together.

**Status:**

* **No proof of connectivity. Collatz OPEN.**
* PROVED (elementary): Proposition 1 (two-place density of `Gamma_1`), Corollary 1 (finite views), Proposition 2 (quantifier exchange), Proposition 3 (the minimal counterexample is a component minimum, with its forward and backward sieves), Theorem B (the transfer barrier).
* FINITE-EXACT: every table, with script
  `04-computation/experiments/collatz_connectivity_from_rigidity_20261001.py` and its output `.out` beside it, ending `ALL CHECKS PASSED`.
* Independent audit OWED; the seventeenth note's audit is in progress.

## 0. The answer in one paragraph

Rigidity says that the unlabelled Collatz graph remembers every integer:
the backward tree of `n` determines `n`. I tried to turn this into
connectivity along four routes, and each stops at the same place.

* **A second component would have to create a symmetry.** It need not; and
  rigid disconnected graphs exist.
* **Every backward type could be realised in the trivial component.** Every
  *finite* view of every integer is realised there (new, Proposition 1),
  but the *infinite* view is realised only by the integer itself.
* **A minimal counterexample would be trapped.** It is trapped, but the minus
  sheet has integers caught in exactly the same trap (5 and 17).
* **Density or ergodicity.** Finite places see one component even on the
  minus sheet, which has three components of densities `0.327, 0.325, 0.348`.

The decisive fact is Theorem B. Rigidity holds verbatim for `3n−1` on `N`
(three components) and for `5n+1` (three cycles, positive drift), so it is
blind to the sheet sign and to the drift: the two coordinates the synthesis
demands of any proof mechanism. Rigidity adds exactly one thing beyond
finite views, namely the integer itself, and that is what the conjecture is
about. So it cannot carry the proof. It *does* locate the remaining content
exactly: **Collatz ⟺ every positive integer is its own bounded 3-adic
approximant from the trivial component** (Proposition 2). This is a pure
statement about archimedean size, and on the minus sheet its failure is
visible: the approximants of 5 and 17 grow like `3^k`.

## 1. What rigidity puts on the table

* From the seventeenth note:
  * `Aut(N, T) = 1`.
  * Backward trees identify vertices prime to 3.
  * Components are pairwise non-isomorphic and each is rigid.
  * The trivial component `Gamma_1` is categorical.
  * `Collatz ⟺ (N, T) ≅ Gamma_1`.
* A proof of connectivity from these would have to produce, from a second
  component `K`, either a contradiction with rigidity or a contradiction with
  categoricity. The routes below are the natural ways to try.

## 2. Route 1: a second component forces a symmetry? (fails)

The candidate symmetries are of three kinds.

* **Swapping `K` with `Gamma_1`.** Impossible: the labels differ, and `K` has
  no branch three-cycle.
* **Shifting along a divergent ray.** Not an automorphism, since in-degrees
  change along the ray.
* **Exact versions of approximate symmetries.** These come from recurrence:
  pigeonhole gives `n_k ≡ n_k' (mod 3^J)` along any orbit. By Lemma 3 such
  coincidences break by backward depth about `5J`. A coincidence deep enough
  relative to size would force `n_k = n_k'`, which is a cycle, but pigeonhole
  delivers only coincidences shallower than size. That is THM-4476's
  thin-divergence mechanism, already exhausted to within `(log X)^0.514`.

Nothing turns a second component into a symmetry, and Theorem B exhibits
rigid graphs that do have several components.

## 3. Route 2: realise every backward type in the trivial component (stops at the quantifier exchange)

By rigidity, connectivity holds iff every `n` has some `m ∈ Gamma_1` with
`B(m) ≅ B(n)` (then `m = n`). The finite-depth version is a theorem.

**Proposition 1 (two-place density of `Gamma_1`).** For all `a, b ≥ 0` and
every residue `r mod 2^a 3^b`, infinitely many `n ≡ r (mod 2^a 3^b)` reach 1.
Explicitly, each such `n` reaches, after `a` shortcut steps, a trunk-entry
point `a_j = (4^j − 1)/3`.

*Proof.*

1. **Trunk points cover every class modulo `3^k`.** `4 = 1 + 3` has order
   `3^k` modulo `3^(k+1)`. So `j ↦ 4^j` is a bijection from `Z/3^k` onto
   `1 + 3Z/3^(k+1)`, and `j ↦ (4^j − 1)/3 mod 3^k` is a bijection of `Z/3^k`.
   Each `a_j` is odd with `3a_j + 1 = 4^j`, so `a_j ∈ Gamma_1`.
2. **Terras reduction.** Let `r_0 ∈ [1, 2^a]`, `r_0 ≡ r (mod 2^a)`. Every
   `n = r_0 + 2^a t` shares the first `a` shortcut steps of `r_0`, so
   `T_1^a(n) = A_0 + 3^p t`. Here `p` is the number of odd steps and
   `A_0 = T_1^a(r_0)`.
3. **Matching the 3-adic class.** `n ≡ r (mod 3^b)` iff `t ≡ t_0 (mod 3^b)`.
   So `T_1^a(n)` runs over the values `≥ L_0 = A_0 + 3^p t_0` in the class
   `L_0 mod 3^(p+b)`.
4. **Choice of `j`.** Choose `j` with `a_j ≡ L_0 (mod 3^(p+b))` and
   `a_j ≥ L_0`, and set `t = (a_j − A_0)/3^p`. Raising `j` by multiples of
   `3^(p+b)` gives infinitely many `n`. ∎

FINITE-EXACT (`.out`, section B): the construction was run on all `2520`
classes modulo `2^a 3^b` with `a ≤ 5`, `b ≤ 3`. In each, the constructed `n`
lies in the class and `T_1^a(n) = (4^j − 1)/3`, checked by iteration.

**Corollary 1 (finite views).** The radius-`R` view of a vertex is its forward
path of length `R` together with the depth-`R` backward trees along it.
Wherever it is acyclic, it depends only on the class of the vertex modulo
`2^(R+1) 3^R` (the parities and 3-adic labels of the iterates are affine in
`n`). Hence **every acyclic finite view of every vertex of `(N, T)` occurs
in `Gamma_1`**, at infinitely many vertices.

* A hypothetical divergent component is therefore locally indistinguishable
  from the trivial one.
* A hypothetical cycle of length `L` *is* visible, in views of radius `≥ L`.

This is the cycle-half / divergence-half split again, now as local
realisability.

**Proposition 2 (the quantifier exchange).** For `n ∈ N` put
`δ_k(n) = min{m ∈ Gamma_1 : m ≡ n (mod 2·3^k)}`, which is finite by
Proposition 1. Then

`n ∈ Gamma_1  ⟺  sup_k δ_k(n) < ∞  ⟺  δ_k(n) = n for all large k`.

*Proof.* If `δ_k(n) ≤ C` for all `k`, take `3^k > max(C, n)`: two positive
integers below `3^k` that agree mod `3^k` are equal. ∎

So connectivity is the exchange of "for every depth there is a witness in
`Gamma_1`" (a theorem) into "there is a witness for every depth". The
exchange holds exactly when the witnesses stay bounded, and boundedness is
archimedean.

On the minus sheet, where second components exist, the failure is visible.
Table from `.out`, section C (basins of `3n−1` on `[1, 3·10^6]`, densities
`0.3271 / 0.3244 / 0.3485`):

| `n` (component) | `δ_k(n)` for `k = 1, ..., 11` |
|---|---|
| 5 (`{5,7,10}`) | 11, 59, 59, 653, 1463, 1463, 4379, 39371, 39371, 118103, `> 3·10^6` |
| 17 (`{17,...}`) | 11, 53, 71, 827, 1961, 4391, 4391, 52505, 118115, 118115, 2480075 |
| 41 (`{17,...}`) | 11, 59, 95, 365, 3929, 13163, 13163, 13163, 275603, 944825, 1062923 |

`δ_k(n)/3^k` stays between about 2 and 14: the approximants grow like
`3^k`, the signature of a genuine second component. On the plus sheet
`δ_k(n) = n` for every `n` in any verified range, and the conjecture says
this persists.

## 4. Route 3: the minimal counterexample's two-sided trap (realised on the minus sheet)

**Proposition 3.** If Collatz fails and `n_0` is the least positive integer
outside `Gamma_1`, then `n_0` is the minimum of its component: every
forward iterate, every predecessor and every cousin of `n_0` exceeds `n_0`.

* **Forward (classical).** `n_0` is odd and `n_0 ≡ 3 (mod 4)`, and its orbit
  never drops below `n_0`. This confines `n_0` to a 2-adic set of measure
  zero (the stopping-time sieve): within `k` steps the surviving fraction is
  `≈ 2^((h*−1)k)`.
* **Backward.** `n_0` has no smaller predecessor. So `n_0 ≢ 2 (mod 3)`,
  because `(2n_0 − 1)/3` would be one, and `n_0 ≢ 4 (mod 9)`, because
  `(8n_0 − 5)/9` would be one. In general `n_0` survives the backward sieve:
  no predecessor of multiplier `2^e/3^j < 1`.
* **The backward sieve is weak.** The fraction of odd 3-adic classes prime to
  3 that survive depth `D` is
  `0.500, 0.333, 0.333, 0.315, 0.315, 0.315, 0.311, 0.305, 0.305, 0.305`
  for `D = 4, ..., 22` (`.out`, section D). Multiples of 3 always survive,
  since their backward tree is a bare ray. Unlike the forward sieve, the
  backward one appears to keep a positive proportion. This is consistent with
  sparse, frequently blocked branching, but it is only observed.

*Proof.* Predecessors and iterates of `n_0` lie in its component, and
everything below `n_0` lies in `Gamma_1`. The congruences are direct
computations. ∎

**The trap is not contradictory.** On the minus sheet `5` and `17` are
exactly such component minima (`.out`, section E). Any contradiction must
therefore use the sign of the carry. On the minus sheet, cycle words with
`2^A < 3^p` close up: the cycles through 1, 5 and 17 have
`(A, p) = (1, 1), (3, 2), (11, 7)`. On the plus sheet a cycle word must have
`2^A > 3^p`.

## 5. Route 4: density and ergodicity (fails on the minus sheet)

* Proposition M (no periodic invariant) and Proposition 1 say the finite
  places see one component.
* An ergodic argument ("an invariant set has density 0 or 1") would close
  the gap. It is false for this kind of graph: the minus sheet's three
  basins have densities near `1/3` each (HYP-9165; recomputed here on
  `[1, 3·10^6]`).
* THM-4517's sheet cap applies directly. A sheet-blind method cannot prove
  even a cycle-uniform basin density above the smallest minus-sheet basin.

## 6. Theorem B: the transfer barrier (PROVED)

**Theorem B.**

1. The conclusions of Theorem R hold for the graph of `3n − 1` on `N`. That
   graph is `(−N, T)` relabelled by `n ↦ −n`, and Theorem R(iii) covers `Z`.
   It has at least three components (cycles through 1, 5, 17).
2. The conclusions of Theorem R_a (`a = 5`) hold for `5n + 1` on `N`, which has
   at least three components (cycles through 1, 13, 17) and positive drift.
3. Proposition 1 and Corollary 1 hold on the minus sheet with the trunk
   points `(2^(2j+1) + 1)/3`, which cover every class modulo `3^k`.

Consequently no argument whose inputs are rigidity, finite views,
Collatz-typed local rules and absence of periodic invariants can prove
connectivity. All of these hold for `3n−1`, which is disconnected. The
missing input must distinguish `N` from `−N`.

FINITE-EXACT (`.out`, section A):

* `3n−1`: the depth-25 backward trees of the 1334 vertices `n ≤ 2000` prime
  to 3 are pairwise distinct; its cycles have least elements `1, 5, 17`.
* `5n+1`: the depth-33 trees of the 500 vertices `n ≤ 624` prime to 5 are
  pairwise distinct; its cycles have least elements `1, 13, 17`.

**Positivity is invisible at every finite depth.** For every `n > 0` and every
`D` there are negative integers with the same depth-`D` backward tree: take
any `m < 0` with `m ≡ n (mod 2·3^⌈D/2⌉)`. Positivity is a property of the
infinite tree, which rigidity guarantees determines it, and of no truncation.

In the synthesis's terms, rigidity is **sheet-blind** (it holds for `3n−1`)
and **drift-blind** (it holds for `5n+1`). Those are exactly the two
coordinates a proof mechanism must see together.

## 7. What a proof would have to add

The remaining content is a single archimedean statement in three equivalent
dresses:

* **Approximation.** For every `n`, `δ_k(n)` is bounded (Proposition 2).
* **Components.** 1 is the only positive integer that is the minimum of its
  component (Proposition 3).
* **Two places.** A component minimum must survive both the 2-adic forward
  sieve and the 3-adic backward sieve. The two sieves interlock only through
  size: a positive `n < 6^L` is determined by its residues modulo `2^L` and
  `3^L`. The conjecture says no integer above 1 survives both all the way
  down.

This last form is the synthesis's "2-versus-3 digit transversality",
reached from the graph side. Rigidity is the statement that the 3-adic side
is a complete invariant. It supplies no transversality, because the minus
sheet satisfies everything rigidity says and is not transversal: 5 and 17
survive both sieves.

## 8. Verdicts

| claim | status |
|---|---|
| Collatz connectivity from rigidity | **NOT PROVED**; blocked by Theorem B |
| Proposition 1 (`Gamma_1` dense in `Z_2 × Z_3`; trunk points cover `Z/3^k`) | PROVED + FINITE-EXACT (2520 classes) |
| Corollary 1 (every acyclic finite view occurs in `Gamma_1`) | PROVED |
| Proposition 2 (quantifier exchange: `Gamma_1` membership ⟺ bounded `δ_k`) | PROVED; minus-sheet growth FINITE-EXACT |
| Proposition 3 (minimal counterexample = component minimum; sieves) | PROVED; backward survival table FINITE-EXACT |
| Theorem B (rigidity holds for `3n−1`, `5n+1`; positivity invisible at finite depth) | PROVED + FINITE-EXACT |
| backward sieve keeps a positive proportion | OBSERVED (D ≤ 22) |
| Collatz | OPEN |

## 9. Directions

* **D62.** The limit of the backward-sieve survival fraction (observed
  `0.305` at depths 18–22): a branching random walk with killing on the
  3-adic labels; is the limit positive?
* **D63.** A sheet-sensitive tree invariant. Positivity is an infinite-depth
  property (§6). Is there a *uniform* depth-`D(n)` test, `D(n) ≍ log n`, that
  reads the sign from the tree? That would be the first rigidity-side
  quantity able to see the sheet.
* **D64.** Growth law of `δ_k(n)` on the minus sheet: observed `≍ 3^k`. Is it
  `≍ 3^k/dens(basin of 1)`? A plus-sheet analogue would bound how a
  hypothetical second component's approximants must grow.
* **D59 (continued).** Full profinite density of `Gamma_1` at primes `≥ 5`.
  These primes are passive for the graph, so Corollary 1 does not need them.
