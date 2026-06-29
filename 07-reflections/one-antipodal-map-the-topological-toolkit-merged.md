# One antipodal map: Borsuk–Ulam, Ky Fan, Ham Sandwich, Kaczynski, all on the arc-hypercube

*mac-mini-2026-06-29-S12. Merging the topological-combinatorics toolkit into the per-level signed-cycle-index thread, prompted by arXiv:2512.09332 (El Sahili–El Zein, oriented Hamiltonian paths under arc deletion) and the owner's list (Borsuk–Ulam, Kaczynski, Ky Fan, Ham Sandwich). New: HYP-3544.*

## The map that ties the toolkit together

klein's THM-584 named the object the whole topological toolkit lives on: the **arc-hypercube `Q_d`** (`d = C(n,2)`, one bit per pair), with the **complement as its antipodal map** `x ↦ x ⊕ 1`. That single identification turns four classical theorems into four readings of one free `Z_2`-action:

- **Borsuk–Ulam** is *the* theorem about that antipodal map. On the LRC side it is the **witness**: lonely times come in antipodal pairs `{t*, −t*}`, the saddle index is `(p−1)/2` (odd for `p ≡ 3 mod 4`, so `n=14` is index 3 → Borsuk–Ulam, even → Brouwer/SOS), and the odd degree forces `M(S) ≥ 1/14`.
- **Ky Fan** is the *combinatorial* Borsuk–Ulam — it counts alternating simplices and returns an **odd** number. On the tournament side it is **Rédei/Forcade graded by type**.
- **Tucker** is the discrete `Z_2`-labeling (no complementary edge) — the same antipodal action, asking for a labeling, *unused* so far.
- **Ham Sandwich** is the measure-equipartition corollary — *unused*, but the natural tool for the R-odd cap obstruction `M_odd`.

And **Kaczynski** sits on the *other* (analytic, `R`-even) side: the sieve machinery (`Σ μ²/φ ∼ log x + 1/ζ(2)`, `Σ φ`) supplying the `1/ζ(2) = 6/π²` totient density that makes the floor positive. So the same `R = ` complement = antipodal map splits the whole program into an `R`-even analytic half (Kaczynski, Brouwer, SOS, the floor) and an `R`-odd topological half (Borsuk–Ulam, Ky Fan, the witness, the obstruction).

## What I grounded: Ky Fan is Rédei/Forcade, graded by type

OPEN-Q-059 has long asked to read Rédei as a *tournament* Ky Fan lemma. The bridge is the **type-hypercube** `{+,−}^{n-1}`: an oriented path's type `τ` is a corner, reversal and complement are its antipodal `Z_2`, and Ky Fan's odd alternating-count is the per-type Hamiltonian-path count. I verified the keystone fact (Forcade 1973): **`N_τ(T) mod 2` is tournament-independent** — for every type, across *all* tournaments, the parity is the same, equal to the transitive (descent-set permutation) count. Not one type out of `2^{n-1}` had variable parity at `n=4,5,6`. And at `n = 2^k` (n=4) *every* type is odd — Forcade's special case. The directed corner (Rédei) and the antidirected corner (El Sahili–Abi Aad, `≡ 2 mod 4`) are two corners of one graded Ky-Fan odd-count, and the **signed cycle index (klein's open HYP-3540) is its generating function** — the metagraph eigenvalue `d − 2k` (reversed-arc level) and the per-type parity are two gradings of the *same* `Z_2`-equivariant counting on `Q_d`.

That oriented-path counts have tournament-independent parity is exactly the kind of statement Ky Fan delivers (the count is a *topological* invariant of the antipodal action, blind to which tournament realizes it). Proving the tournament Ky Fan — that this parity *is* a Ky-Fan alternating count for arbitrary `T`, not just the transitive order — would make Rédei and Forcade **shadows of Borsuk–Ulam**, which is the entire point of OPEN-Q-059, and the natural home of the signed cycle index.

## The arXiv paper is the edge-stability of this count

El Sahili–El Zein (arXiv:2512.09332) opens with Rédei's odd count and proves: deleting **one arc** from a tournament of order `n ≥ 8` preserves every oriented Hamiltonian path, except two explicit special exceptions. In klein's encoding one arc deletion is **one `Q_d` edge** (a wiggly line). So the paper is the **edge-robustness of the graded Ky-Fan count along the arc-flip graph** — the same graph whose `S_n`-quotient is the metagraph and whose antipodal map is the complement. And its classical exception list contains the **Paley tournament of order 7** (Grünbaum's antidirected exception) — the apex prime `7` of LRC(14), special in the oriented-path theory and the runner theory at once. The three threads — Rédei/Forcade parity (Ky Fan), arc-deletion robustness (the paper, `Q_d` edges), and the metagraph antipodal spectrum (klein) — are one structure on `Q_d`.

## The two unused tools, and where they fit

- **Ham Sandwich** wants to bisect the danger-cover measure (or the R-odd `M_odd`, the cap obstruction) by the spectral eigenvector: the obstruction is a *signed* measure on the antipodally-symmetric danger set, and Ham Sandwich/Borsuk–Ulam says a symmetric perturbation cannot move it to zero — a possible existence-without-construction certificate for the cap bound, complementary to the descent route.
- **Tucker** wants a `Z_2`-labeling of `Q_d` (or the metagraph) with no complementary edge; its non-existence is the discrete Borsuk–Ulam, and a Tucker labeling encoding "which runner binds" would be the finite, checkable form of the saddle-index obstruction.

## The one-line merge

Everything in the owner's list is one free `Z_2`-action — the complement = antipodal map of the arc-hypercube. Its `R`-even half is analytic (Kaczynski, `1/ζ(2)`, the floor, Brouwer/SOS); its `R`-odd half is topological (Borsuk–Ulam, Ky Fan, the witness, the obstruction). Rédei/Forcade are the graded Ky-Fan odd-count (verified tournament-independent); the signed cycle index is its generating function; arc deletion is its edge-stability (the arXiv paper); and the apex prime `7` is where the tournament and runner instances coincide. The remaining moves are named: prove the tournament Ky Fan (OPEN-Q-059), and try Ham Sandwich on `M_odd`.

See [[the-one-involution-three-spectra]] (HYP-3543, klein THM-584), [[the-two-indices-redei-is-odd-lonely-is-even-half-tiling-is-the-quotient]] (THM-582), [[the-saddle-index-is-p-minus-1-over-2-borsuk-ulam-forcing-meets-the-vitali-core-kps]]. New: HYP-3544; OPEN-Q-059 (grounded). External: arXiv:2512.09332.
