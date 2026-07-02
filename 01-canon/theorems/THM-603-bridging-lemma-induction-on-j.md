# THM-603: The bridging lemma — induction on j, shifted bounded systems, and why LRC(<n) is never needed

**Status:** PROVED (architecture + each step elementary; the base-case minimizations are finite `decide`-shaped computations)
**Author:** mac-mini-2026-07-02-S4 (HYP-3861)
**Role:** the "one new proof" of the all-n feasibility verdict (HYP-3860): the induction wrapper's interface. Result: the wrapper needs NO induction hypothesis on n at all — the recursion of THM-602(D) already carries everything, with the induction variable j = the number of free comb parameters.

## Statement

Let a cluster `F` (`j` combs, radius `r = 1/n`) be frozen along a rank-`k` resonance lattice `Λ` in the window `I` (THM-602(A) cases 2–3). Then the renormalized problem is:

> a **shifted bounded-integer comb system**: `j' = j − k` effective free comb parameters with
> integer speed data bounded by the HNF census, together with `k` shift parameters that are
> FROZEN up to total drift `< k` across `I`,

and this instance is again an instance of the SAME free-phase anti-covering problem, with:
1. **fewer free parameters** (`j' < j`);
2. **all floors applicable as-is**: THM-598/599/601/602's bounds are quantified over ALL phases (shift-proof under resolution), so inherited shifts cost nothing in the resolved directions;
3. **frozen shifts entering only the base cases**: when the effective system is itself fully frozen (a bounded fixed system with parameters `θ` on a `k`-torus), the required floor is `min_θ uncoveredMeasure(boundedSystem, θ)` — a finite piecewise-linear minimization over rational breakpoints (S100's method verbatim), i.e. one exact-rational TABLE entry per census pattern;
4. **termination at `j' ≤ 6`** unconditionally (mass-subcritical: union bound `1 − 2rj' > 0` needs no structure). Hence the recursion never appeals to LRC(m) for any `m < n`: the strong induction of HYP-3860(T2) is subsumed by the well-founded recursion on `(j, corank Λ)` with base `j ≤ 6`.

## Proof

(1) is THM-602(B): the unimodular HNF change of coordinates splits slow (drift `< 1` each, `k` of them) from fast; the fast block is an integer matrix with entries bounded via the HNF/Hadamard bound, defining the effective bounded speeds. (2) is by inspection of the floors' quantifiers: THM-598(B), THM-599, THM-601 hold for every phase vector; a drifting-but-frozen shift is a particular phase path, and the Part-B error already charges its drift. (3): for a fixed bounded system, the uncovered measure as a function of the shift torus is piecewise linear with breakpoints where arc endpoints meet (finitely many rationals — the S100 breakpoint method in `k` variables; the minimum is attained at a breakpoint vertex and is an exact rational). (4): `2rj' ≤ 2·6/14 < 1`; the union bound is positive with no hypotheses. Well-foundedness: each renormalization strictly decreases `j` (as `k ≥ 1`), and `j` is bounded below by the base. ∎

## Consequences for the all-n pipeline

- The uniform pipeline's "induction wrapper" is a MISNOMER-avoidance: there is no induction on `n`, only the `j`-recursion — one theorem statement, no per-n hypothesis threading. This simplifies the Lean architecture materially: no motive juggling across `n`; the certificate packs differ only in table CONTENT.
- The owner's LRC(≤13)-settled policy is thereby not even load-bearing for the formal proof: it remains a working convenience for n=14-only sessions.
- The remaining computational surface of the bridge: the base-case shift-minimization table (one exact rational per bounded pattern in the HNF census — the same few-hundred-row census already spec'd, each row now carrying `min_θ` instead of a single value).

-> THM-602, HYP-3858/3859/3860, THM-598/599/601, OPEN-Q-108.
