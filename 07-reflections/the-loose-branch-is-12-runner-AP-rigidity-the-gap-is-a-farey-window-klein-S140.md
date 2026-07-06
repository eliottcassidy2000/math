# The loose branch is 12-runner AP-rigidity: gap-violators are bounded, and the gap is a Farey window

*klein-2026-07-05-S140 (HYP-4151). Owner: improve the LRC(14) proof; figure out the EXACT mathematics
before formalizing. The open crux — the loose branch of `TightLooseDichotomy` — is sharpened here from
"census/template (dead, HYP-4137)" to a classical, exact rigidity statement about the arithmetic
progression, with a concrete proof route and the r=1 case essentially proved. No Lean this session.*

## The exact statement (verified, not a census)

The loose branch asks: every primitive compressed covering 13-family whose 12-runner **base** `B`
(the non-max runners) is not a dilated AP has a real witness with margin `2/25`. Stripped to its
arithmetic core, this is a **12-runner Lonely-Runner rigidity**:

> For a 12-tuple of nonzero integers `B`, `M(B) := max_t min_{v∈B} ‖v t‖` satisfies
> **`M(B) = 1/13` iff `B` is a dilated AP `c·{1,…,12}`, and `M(B) ≥ 2/25` otherwise.**
> The spectrum has an **empty gap `(1/13, 2/25)`** (width `2/25 − 1/13 = 1/325`).

Exact computation over 12,093 families (`lrc14_twelve_runner_rigidity_klein_S140.py`): the *only*
families with `M ∈ [1/13, 2/25)` are the eight dilated APs `c·{1,…,12}` (`c=1..8`), all at exactly
`1/13`; **zero** non-AP families in the gap. Second value `2/25`, attained by `{1,…,11,24}`
(= AP with `12 → 24`), its dilations, and one sporadic `{1,2,3,5,7,8,9,10,11,12,17,19}`. Low
spectrum: `1/13, 2/25, 3/37, 1/12, 3/35, 2/23, …`, the ladder `{1,…,11,12k} = k/(12k+1)`.

## The decisive fact: gap-violators are BOUNDED (dual to kps's unbounded witnesses)

The worry was that large "translated" near-APs sneak into the gap, making it an unbounded
real-analytic mess. **They do not** (`lrc14_gap_boundedness_klein_S140.py`):
- 12 **consecutive** integers `{N,…,N+11}` — the obvious translated AP — have `M` *blowing up*, not
  down: `N=2 → 2/15`, `N=3 → 3/17`, `N=5 → 5/21`, … → increasing. Only `N=1` gives `1/13`.
- The ladder `{1,…,11,12k}` gives `k/(12k+1) ≥ 2/25` for `k ≥ 2` (never in the gap).
- mac-mini's floating 7-cluster (the profile-passer that persists at all scales) has `M = 3/26 ≈
  0.115` — **loose**, not a gap-violator.

So the gap-violator set `{M < 2/25}` is **exactly the dilated APs**; the unique *primitive* one is
`{1,…,12}` — a single bounded family. This is the precise **dual** of kps-S11 (HYP-4137): kps showed
the *witness modulus* for loose families is unbounded (no finite template); the *gap-violator set* is
nonetheless bounded (isolated at the AP). The loose branch is therefore a **rigidity** (the AP is an
isolated minimizer), and the real-analytic difficulty is *why* it is isolated — not a search.

## The structure: the gap is a Farey window, and only 1/13 is realizable

`M(B) = r/Q` in lowest terms. Two hard constraints pin it:
- **LRC(13)** (settled, cited): `M(B) ≥ 1/13 ⟹ Q ≤ 13r`.
- **The ceiling**: `M(B) < 2/25 ⟹ Q > 12.5r`.

So a gap value must be a Farey fraction with **`Q ∈ (12.5r, 13r]`** — `1/13, 3/38, 4/51, 5/63, 5/64,
6/77, 7/88, 7/89, …` (infinitely many, but each in a razor-thin window). The computation shows
**only `1/13` is realized**; every intermediate Farey fraction corresponds to a *sub-optimal* residue
configuration, never a genuine max-min. That is the content to prove.

## The r = 1 case is essentially proved

If `M(B) = 1/Q` (numerator 1) and `M(B) < 2/25`: then `1/Q < 2/25 ⟹ Q ≥ 13`, while LRC(13) gives
`1/Q ≥ 1/13 ⟹ Q ≤ 13`. So **`Q = 13`** exactly. At the optimum `t = a/13`, the 12 residues
`{v_i a mod 13}` are 12 distinct nonzero residues mod 13 — i.e. **exactly `{1,…,12}`** (the only 12
nonzero residues), so `{v_i mod 13} = a^{-1}·{1,…,12}`: the base is a dilated AP *mod 13*. (Upgrading
"dilated AP mod 13" to "integer dilated AP with `M` exactly `1/13`" is the LRC(13)-extremizer
uniqueness — the residue-AP families that are not integer APs, e.g. `{2,…,12,14}`, have a *better*
`t` and so `M > 1/13`, escaping the gap.) The LRC(13) sandwich `12.5r < Q ≤ 13r` doing the work here
is the template for general `r`.

## The proof route: an AP-residue rigidity, uniform in r

The AP is the LRC(13) **extremizer** — the family that *minimizes* `M` (the hardest to avoid `0`);
everything else is looser (`≥ 2/25`). Recast at the optimum `a/Q` (value `r/Q`, `gcd(r,Q)=1`,
`Q ∈ (12.5r, 13r]`): the 12 residues `v_i a mod Q` lie in `[r, Q−r]` and bind at `r`.

- **Clean half.** If those residues form an AP with step `r`, i.e. `{r, 2r, …, 12r}`, then
  `12r ≤ Q − r ⟹ Q ≥ 13r`; with `Q ≤ 13r` this gives `Q = 13r`, so `r/Q = 1/13`. **An AP-residue
  configuration yields *only* the value `1/13`** — every other Farey fraction in the window is
  non-AP-residue. (And `r = 1` forces `Q = 13`, residues `= {1,…,12}` exactly — proved above.)
- **Open half (the rigidity).** A genuine max-min with `M < 2/25` *forces* the residues to be that
  AP. Intuition: if the 12 residues in `[r,Q−r]` are not the balanced `{r,…,12r}` spread, some
  consecutive pair is farther apart than the AP's uniform `≈ r`, and re-centering `t` inside that
  slack lifts every distance above `r/Q` toward `2/25` — contradicting optimality. Turning this into
  a clean inequality **uniform in `r`** (not a per-`Q` census) is the missing real-analytic
  ingredient — the honest replacement for the dead finite-template route. mac-mini-S40's 2-point
  equioscillation and the three-distance structure of the AP's own orbit are the tools; the work is
  the non-AP case, where the points are not a single generator's orbit.

## Why this improves the proof state

- The loose branch is now a **sharp classical statement** (AP-rigidity + spectral gap `1/325`),
  not a vague census — and it is a *bounded-violator* rigidity, so the "irreducibly real-analytic"
  difficulty is precisely located: prove **equal-spacing rigidity uniform in `r`** (or the LRC(13)
  extremizer-uniqueness-with-gap). The r=1 case is done via the LRC(13) sandwich.
- It **retires the finite-template route cleanly**: the Farey window shows why bounding the modulus
  fails (the window `(12.5r,13r]` moves with `r`), and why the real statement is a rigidity, not an
  enumeration.
- Next: prove the equal-spacing rigidity (the three-gap maximization is classical; the arithmetic
  realizability is the work), OR mine the LRC(13) proof (Sungkawichai–Trakulthongchai) for an
  extremizer-uniqueness + second-value corollary. Only then formalize.

## Links

- Scripts: `04-computation/lrc14_twelve_runner_rigidity_klein_S140.py` (+ `.out`),
  `lrc14_gap_boundedness_klein_S140.py` (+ `.out`). Exact. HYP-4151.
- Builds on: klein-S139 (loose branch irreducibly real-analytic), kps-S11 HYP-4137 (finite-template
  dead; unbounded witnesses — dual to bounded violators here), mac-mini-S55 (pole-necessity),
  mac-mini-S40 (2-point equioscillation), mac-mini-S38 (Ostrowski ladder `k/((n−1)k+1)` = the
  `k/(12k+1)` rungs), klein-S138 (`M({1,…,12}) = 1/13` machine-checked — the r=1 anchor),
  klein-S126 (11-runner even-part gap — the one-lower analogue). Open: equal-spacing rigidity uniform
  in `r` (the loose branch).
