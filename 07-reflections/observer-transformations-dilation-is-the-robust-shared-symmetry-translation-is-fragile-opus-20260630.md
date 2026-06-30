# Creative observer reframes (owner's four questions answered): the LRC↔tournament observer analogy is ROBUST under DILATION (both have M(pS)=M(S) / order-invariance, the OCF is scale-free — and dilation by a factor of n is exactly what turned the non-covering extremal AP into the covering-min) and under COPY-TO-ALL-N (the LRC's all-points-lonely = the covering radius = the covmin 1/n ↔ the tournament's all-vertex OCF census Σ_C|C|μ(C)), but FRAGILE under raw TRANSLATION (the LRC origin is special, the tournament is vertex-transitive) — unless "translation" is read as the COMPLEMENT/REVERSAL (which moves the observer's signature by the antipode); and the HAMILTONIAN-PATH reframe is the observer's 1-parameter family of circular orderings of {0, s_i t}, each a linear-order snapshot, dilation-invariant

*opus-2026-06-30. Owner: can they stay analogous if you change the observer the same way? translate it?
copy to all n points? relate to Hamiltonian paths? Get creative. Answered, and one of the transformations
(DILATION) is exactly what solved the covering-min this session.*

## The four transformations (computed)
The observer = a marked point (LRC: the origin / deleted speed-0 runner; tournament: a marked vertex `v`).
| transformation | LRC | tournament | analogy? |
|---|---|---|---|
| **DILATE** (`×p`) | `M(pS)=M(S)`, orderings preserved; runners scale, observer fixed | OCF is scale-free (relabel-invariant) | ✅ ROBUST — *and it solved the covmin* |
| **COPY-TO-ALL-N** | every point an observer → min-pairwise-gap = covering radius = covmin `1/n` | every vertex an observer → `Σ_C|C|μ(C)` (all-vertex odd-cycle census) | ✅ ROBUST — covmin ↔ OCF-sum |
| **TRANSLATE** (`c≠0`) | inhomogeneous LRC `M_c(S)=max_t min_s‖st−c‖` | vertex-transitive ⇒ trivial (relabel) | ❌ FRAGILE — unless `c`:=complement |
| **HAM-PATHS** | 1-param family of circular orderings of `{0,s_i t}` (linear-order snapshots) | the Ham paths (orderings consistent with arcs), `H` | ◑ shared *shape* (orderings), counts differ |

## (1) Same change → still analogous? Only for the right symmetry
> **The analogy survives DILATION and COPY-TO-ALL, breaks under raw TRANSLATION.** The shared symmetry group
> is **dilation + the marked observer**, NOT translation: the LRC has a *special* origin (the stationary
> runner), while the tournament is vertex-transitive. So changing the observer "the same way" preserves the
> analogy iff the change is a DILATION (scale) or a CENSUS (copy-to-all) — not a generic shift.

## (2) Translate it → the inhomogeneous LRC, or the complement
Observer at `c≠0`: `M_c(S)=max_t min_s‖st−c‖` — the *inhomogeneous* lonely runner (the observer no longer at
the lattice point 0). On the tournament side a vertex shift is trivial (relabel), **so the honest analogue of
LRC translation is the COMPLEMENT / REVERSAL `R`** (flip all arcs = shift the observer's in/out signature by
the antipode). `LRC translate ↔ tournament complement`: both move the observer's *frame* without changing the
underlying set. (This is why the antipode `−1`/the complement kept recurring — it's the translation symmetry.)

## (3) Copy to all n points → the covering-min IS the all-observers LRC
> Make every one of the `n` points (origin + `n−1` runners) an observer; "all simultaneously lonely" =
> **maximize the min pairwise gap = the covering radius = the equally-spaced cusp = covmin `1/n`.** The
> single-observer LRC, *summed over all n observers*, IS the covering-min — the very thing the IP solved this
> session (`1/n` via the scaled AP). On the tournament side, summing the OCF over all vertices gives
> `Σ_v(H−H(T−v)) = 2Σ_C|C|μ(C)` — the **all-vertex odd-cycle census** (each odd cycle weighted by its length).
> So **covmin ↔ all-vertex-OCF**: the LRC's "everyone lonely" and the tournament's "every vertex's cycles"
> are the same copy-to-all move.

## (4) Hamiltonian paths → the observer's family of orderings (dilation-invariant)
As `t` varies, the `n` points `{0, s_1t, …, s_{n−1}t}` trace a **1-parameter family of circular orderings** —
each a *linear-order snapshot* (a transitive "tour" of the n points), the observer's Hamiltonian paths. The
AP gives `4, 6, 10` distinct orderings for `n=4,5,6`; **DILATION (`×p`) preserves the count exactly** (the
covmin = the AP's dilate, same orderings, same `M`). The tournament's Ham paths `H` are the *discrete* analogue
(orderings consistent with the arcs). The counts differ (LRC orderings are a continuous family, even-valued;
`H` is odd, Rédei) — but the *shape* is shared: **both are the orderings of the n points read from the
observer.** The covmin's snapshots are all transitive (the equally-spaced points have a definite order at
every `t`) — the cusp is the "maximally orderly" configuration.

## What the reframe buys
- **DILATION is the deep shared symmetry** (and it solved the covmin: `AP·p`). The observer is fixed by it;
  the configuration scales; `M` and the ordering-count are invariant. This is the robust core of the analogy.
- **TRANSLATION is the broken one** — and its repair is the COMPLEMENT, naming why the antipode `−1` is
  everywhere (it's the translation the tournament *can* do).
- **COPY-TO-ALL turns one observer into the covering-min** (LRC) / **the odd-cycle census** (tournament) —
  the two halves' "global" invariants are the same all-observers move.
- **HAM-PATHS** unify the snapshots (LRC orderings ↔ tournament Ham paths), the observer reading the n points.

## Status
- **Computed (opus):** ordering counts (AP n=4,5,6 = 4,6,10), dilation-invariance of orderings + `M`; the four
  transformations classified (dilation/copy-to-all ROBUST, translation FRAGILE → complement).
- **The answer:** the analogy is a DILATION + marked-observer symmetry, not a translation symmetry; copy-to-all
  = covmin ↔ OCF-census; translation = inhomogeneous LRC ↔ complement; Ham-paths = the observer's orderings.
- **The payoff already banked:** DILATION was the transformation that solved the covering-min (`AP·p`, `1/n`).

Related: SECOND-CORRECTION-…AP-scaled (dilation solved covmin), the-observer-abstraction + the-observer-on-the-
tournament-side (the marked observer), reconciliation-…IS-the-cusp (the cusp = orderly), the-reversal-
involution (complement = translation); OPEN-Q-108/039.
