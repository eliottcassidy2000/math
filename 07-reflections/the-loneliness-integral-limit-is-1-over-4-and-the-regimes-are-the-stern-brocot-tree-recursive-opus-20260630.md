# Chasing both threads recursively: (1) the LONELINESS INTEGRAL has the clean limit L(AP_n) = 1/4 + 3/(8n) + o(1/n) — the 1/4 is ∫₀¹‖c‖dc = the average nearest-integer distance (the TRIVIAL-loneliness floor, the t=0 witness), and the LRC's hard part is a vanishing 3/(8n) correction (the covmin difficulty has L¹-measure→0); (2) the THREE-DISTANCE REGIMES are the STERN-BROCOT TREE — the ordering complexity O=Φ(n−1) grows by the MODE-A recursion O(n)−O(n−1)=φ(n−1) (vertex insertion splits φ(n−1) regimes via Farey MEDIANTS), each regime's three gaps come from the endpoint denominators (q,q',q+q') = the local continued-fraction state, and the covering-min optimum t=1/n is the SIMPLEST leaf [0;n]

*opus-2026-06-30. Owner: chase both the L limit and the continued fraction, think recursively, connect to
prior work. Both landed: L→1/4 (a clean constant + 1/n correction), and the regimes ARE the Stern-Brocot tree
grown by the Mode-A recursion. The LRC's combinatorics is the Euclidean/Stern-Brocot recursion.*

> **⚠ CORRECTION (opus, same day):** the `1/n` **correction coefficient below is wrong** — I wrote
> `L = 1/4 + 3/(8n)`, but `3/8` was a `Qmax=2n` under-estimate. With adequate Qmax, `M_c(AP_n) = 1/n +
> c(n−2)/n` exactly (linear!), so **`L = 1/4 + 1/(2n)`** (coefficient `1/2`), the **spread regime is all of
> `[0,1/2)`**, and the binding is the **antipode `c*=1/2`** (NOT `c*=1/3`). See
> `CORRECTION-the-inhomogeneous-AP-LRC-is-exactly-linear-…`. **The `L→1/4` limit and the entire
> Stern-Brocot / Mode-A / three-distance content below STAND** — only the `3/n`-coefficient and the
> spread-binding `c*` are corrected.

## Thread 1 — the loneliness integral: L(AP_n) = 1/4 + 3/(8n) + o(1/n)
`L(S) = ∫₀¹ M_c(S) dc`. The **`t=0` witness** (all runners at `0`, observer at distance `‖c‖`) gives
`M_c ≥ ‖c‖`, so `L ≥ ∫₀¹‖c‖dc = 1/4`. Computed:
| n | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 24 |
|---|---|---|---|---|---|---|---|---|---|
| `L` | .313 | .297 | .288 | .281 | .277 | .274 | .271 | .269 | .266 |
| `(L−¼)·n` | .377 | .375 | .376 | .375 | .377 | .376 | .376 | .376 | .375 |
> **`(L−1/4)·n ≈ 0.375 = 3/8`** across all `n` ⇒ **`L(AP_n) = 1/4 + 3/(8n) + o(1/n)`**. So:
> - the **limit is `1/4`** = `∫‖c‖dc` = the **average nearest-integer distance** = the *trivial* loneliness
>   the observer gets for free at `t=0`. As the runners fill in (`n→∞`), the observer can do no better than
>   sit at its `t=0` distance to the lattice point `0`; the mean over observers is `1/4`.
> - the **LRC's hard part is a vanishing `3/(8n)` correction**: `∫(M_c − ‖c‖)dc → 0`. The covering-min
>   difficulty (where `M_c ≫ ‖c‖`) has **L¹-measure → 0** — it is a `1/n` ripple on the `1/4` floor. The
>   `c=0` worst case (`M_0 = 1/n`) is a measure-zero spike; the *bulk* of observer-space is trivial.
> (The `3` plausibly = the three-distance gap count; the `8 = 2³`. Exact derivation open.)

## Thread 2 — the regimes ARE the Stern-Brocot tree (the recursion)
The `O = Φ(n−1)` three-distance regimes (the orderings as `t` sweeps) are organized by the **Stern-Brocot /
Farey** tree, grown by the **Mode-A recursion**:
- **`O(n) − O(n−1) = φ(n−1)`** (verified n=4..12): inserting the new vertex (`n−1 → n`) creates exactly
  `φ(n−1)` **new regimes** — the new Farey fractions of denominator `n−1`, born as **mediants** between
  existing Farey neighbors. This is the project's **Mode A** (vertex insertion, the "fast time scale") realized
  as Stern-Brocot growth: `Φ(n−1) = Φ(n−2) + φ(n−1)`.
- **Each regime's three gaps come from the endpoint denominators.** For a regime `t ∈ (p/q, p'/q')` (Farey
  neighbors, `q p' − p q' = 1`), the mediant is `(p+p')/(q+q')` and the gaps are `{1/(q+q'), 1/q' …}` —
  e.g. `(1/5,1/4)→` mediant `2/9`, gaps `{1/9, 2/9}`. **The gap-pattern IS the local continued-fraction
  state**; the CF of `t` is the Stern-Brocot **address** of the regime. The CF does *not* change the regime
  COUNT (still `Φ(n−1)`) — it gives the **recursive tree structure** and the mediant refinement.
- **The covering-min optimum `t = 1/n` is the simplest leaf `[0; n]`** — the equally-spaced node, where all
  gaps equal `1/n` and `M = 1/n`. The escape lives at a Stern-Brocot vertex (the project's "escape = Farey
  mediant/convergent", now placed: it is the depth-`n` leaf `[0;n]`).

## Recursive synthesis (connecting prior work)
> **The LRC's combinatorics is the Stern-Brocot / Euclidean recursion.** Three threads of the project converge:
> - **Mode A** (`n→n−1`, vertex insertion, fast time scale) = the Stern-Brocot growth `+φ(n−1)` regimes.
> - **Three-distance** (construction gaps `{1,n,2n}`, mac-mini HYP-3702) = the per-regime gap law from the
>   endpoint denominators; the *same* Sós–Steinhaus engine.
> - **The covering-min escape = Farey mediant/convergent** (prior) = the leaf `[0;n]`, `t=1/n`.
> The Stern-Brocot tree is the recursive object underneath all of them — built by mediants (the Euclidean
> algorithm), addressed by continued fractions, counted by `Φ` (Farey length). The covering-min sits at the
> simplest leaf; the orderings are the nodes; the loneliness is the gap at the node.
## Prior-work connections (research sweep — confirms AND sharpens)
The Stern-Brocot frame is already deep in the repo; this result slots in and corrects one placement:
- **The covering-min is the BASE of the Stern-Brocot ray `[0; n−1, k]`** (`one-stern-brocot-ray-floor-
  covering-min-construction`, HYP-3732, HYP-3722): that ray runs `1/n` (k=1) → … → construction `n/Φ₆(n)`
  (k=n) → `1/(n−1)` (k→∞). My `t = 1/n = [0;n] = [0; n−1, 1]` is the **base (k=1)**. **KEY correction-
  consistency:** by SECOND-CORRECTION the true covering-min is the base `1/n`, **not** the construction
  `n/Φ₆(n)` (k=n) — the construction is a *non-extremal* point further up the ray (the refuted frame). The ray
  structure stands; the extremum is its leaf, not its interior.
- **The mediant rule `(p+p')/(q+q')` IS the functional skeleton's `b(+1)` step** (`the-functional-
  decomposition-divide-and-add-one`): the dyadic group `⟨a=÷2, b=+1⟩` *is* the Stern-Brocot/Farey tree. My
  per-regime mediant is the `b`-step; the regimes ARE the tree the skeleton generates — `Φ₆=2T_{n−1}+1` is
  `b` on the doubled triangle.
- **The Sylvester/apex path = the Φ₆ iteration** (`the-recursive-mindset-phi6-iteration-is-sylvester`,
  HYP-3724) — **answers my open question:** Sylvester `2,3,7,43` is the `Φ₆(x)=x²−x+1` path down the tree
  (apex primes `3,7 = Φ₆(2),Φ₆(3)`); Cayley-Dickson `n=2^k+1` is the dyadic spine. Distinct recursions, one
  tree.
- **The regimes = Farey geodesics** (`cuts-as-farey-geodesics`): the witness denominators `2,…,n` are
  consecutive Farey neighbors; cuts step along the tree. My regime boundaries (the Farey fractions) are those
  geodesic steps — and the tight witness `1/n` is killed by multiples of `n`, pushing deeper (the additive-
  mediant / multiplicative-divisibility duality).
- **Three-gap → covering-min** (HYP-3717, `how-the-route-sharpened…cap-kernel`): the `{1,n,2n}` three-gap
  gives the covering set; the cap kernel is a three-gap/Stern-Brocot function with breakpoints at CF
  convergents — the same Sós–Steinhaus engine as the regimes.

## Status
- **Thread 1 (verified):** `L(AP_n) = 1/4 + 3/(8n) + o(1/n)`; limit `1/4 = ∫‖c‖dc` (avg nearest-integer
  distance, the `t=0` trivial floor); LRC hard part is a `3/(8n)` L¹-vanishing correction. Exact `3/8` open.
- **Thread 2 (verified):** regimes = Stern-Brocot nodes; **Mode-A recursion `O(n)−O(n−1)=φ(n−1)`**; per-regime
  gaps from endpoint denominators (= local CF state); covmin `t=1/n` = simplest leaf `[0;n]`.
- **Recursive synthesis:** Mode A + three-distance + covmin-escape all live on the Stern-Brocot tree (mediants
  = Euclid = the skeleton's `b(+1)` step, addresses = CF, count = Φ). The LRC's combinatorics IS the
  Stern-Brocot recursion. **Confirmed against prior work:** covmin = base of the `[0;n−1,k]` ray; Sylvester =
  the Φ₆ path; regimes = Farey geodesics (the cuts).
- **Correction-consistency:** the covmin is the ray's BASE `1/n` (k=1), not the construction `n/Φ₆(n)` (k=n) —
  agrees with SECOND-CORRECTION; the prior "covmin = n/Φ₆(n)" placement was the refuted (non-extremal) frame.
- **Open:** the exact `3/8` coefficient (derive `∫(M_c−‖c‖)dc = 3/(8n)`); the Mode-B (`n→n−2`) analogue on the
  tree; whether the loneliness integral `L`'s `3/(8n)` ripple is the Stern-Brocot depth-`k(n)` (HYP-3732).

Related: the-LRC-for-the-AP-IS-the-three-distance-theorem (the engine), new-invariants-…ordering-complexity
(`O=Φ(n−1)`), the-ordering-complexity-hamiltonian-path-bridge (time vs arcs); Mode A/B, the covering-min escape
= Farey mediant, mac-mini HYP-3702 (three-distance), Φ₆/Sylvester, Cayley-Dickson; A002088 (Φ); OPEN-Q-108.
