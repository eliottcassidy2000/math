# The bottleneck was the constraints, not the problem

*klein-2026-07-01-S61. A reflection on HYP-3779 — brainstorming past the ILP, and the cut that reached the construction scale.*

Last session the exact set-cover ILP resolved the covering-min up to speeds `4n` and then died: at
larger speed bounds the solver's time limit kicked in and the answers turned to noise. I recorded that
honestly as a residual — a beater might hide among the large speeds, in `(4n, n(n-1)]`, where the ILP
couldn't look. Asked this session for methods more creative than the ILP, I started reaching for heavier
machinery — LP relaxations, then SDPs, then the whole zoo of relaxation hierarchies — as if the problem
needed a stronger hammer.

It didn't. The LP relaxation, the first thing I tried, was *too weak*: feasible for every `n`, an
integrality gap, no certificate. That was the wrong direction — loosening the problem to make it fast
throws away exactly the integrality that carries the answer. The thing that worked was the opposite
move: keep the ILP exact, but stop feeding it dead weight. The monolithic ILP carries about forty
thousand danger constraints — one per breakpoint `k/d` — and at `V = n(n-1)` that matrix is what chokes
the solver. But almost all of those constraints are slack at the optimum. Only a couple hundred witnesses
actually bind. So: solve a tiny ILP with just the size and divisibility constraints, look at the candidate
it returns, compute its deepest hole, and add *only that one witness* as a cut. Re-solve. Repeat. At
`n = 12`, after 208 cuts, the ILP went infeasible — proving no covering set with speeds up to `n(n-1)`
beats the construction. The residual I'd flagged is closed for `n = 12`. Same solver, same exactness, but
fed the two hundred constraints that matter instead of the forty thousand that don't.

The lesson is one I keep having to relearn in a different costume. Last session it was "the edge was
tool-resistant, not search-resistant" — the problem was fine, my heuristics had blind spots. This
session it is the sharper version: the *exact* tool was fine too; it was drowning in irrelevant
constraints. When a solver times out, the reflex is to conclude the instance is too big and reach for a
weaker, faster relaxation. But "too big" often means "carrying too much dead weight," not "intrinsically
hard." Cutting-plane and column generation exist precisely for this: the problem is small at the optimum;
you just have to discover *which* small problem it is, one binding constraint at a time. Benders knew
this in 1962. I rediscovered it after a wrong turn through the relaxation hierarchies.

There is a pleasing echo of the mathematics in the method. The whole covering-min story has been about
witnesses — the `q`-witness, the `k`-witness, the `(n+q)`-witness, the sheaf of local dangers. The
cutting-plane algorithm is that theory made operational: each cut *is* a witness, the deepest hole of the
current candidate, and the proof that no beater exists is literally a finite list of two hundred
witnesses that no covering set can simultaneously dodge. The LP dual of the set-cover was always "a
packing of pairwise-incompatible lonely witnesses"; the lazy-cut run computes that packing explicitly, on
demand. The algorithm and the number theory are the same object — a bookkeeping of which specific time
catches which specific escape — seen from the solver's side.

The multi-cut variant then made the point emphatically. Feed *every* lonely witness of each candidate
back, not just the deepest, and `n = 13` and `n = 14` — the target — close in **three rounds** each,
where the single-cut version took 208 for `n = 12`. Same exactness, same `V = n(n-1)`, same certificate
(a finite packing of a few thousand witnesses that no covering set can dodge), but three solves instead
of two hundred. The lesson compounds: not only "feed the solver only the binding constraints," but "when
you must guess which bind, guess generously" — a candidate's every hole is cheap to compute and each is a
valid cut, so there is no reason to trickle them in one at a time.

One honest edge remains. The bound is "speeds `<= n(n-1)`," the construction's own scale; a beater with a
still-larger speed is not excluded by this run, though the construction never needs one and HYP-3745
argues large speeds dig their own holes. That tail is a theory question, not a search one. Everything
below it is settled: the covering-min is the construction `n/Phi_6` for `n = 12, 13, 14`, rigorously, and
the proof is a finite list of witnesses. Don't reach for a bigger hammer when the problem is small at the
optimum. Find the constraints that bind, feed the solver those, and feed it all of them at once.
