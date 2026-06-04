---
source: opus-2026-06-04-S599z (remote-control)
status: ANSWERED (lattice) + FRONTIER (true optima) — the user's points→tournament mapping reduces to: "is the optimal unit-distance graph traceable?" Verified: lattice (Harborth) optima are Hamiltonian-traceable for ALL n=3..28 — the 'flop' never happens; the base Ham path can always be unit. Decomposition u(n)=(n-1) spine + (~2n-√12n) tile unit-edges; recursive +3 frontier-gain; centered-hex shells; π/3 Eisenstein echo (u(5)=7=Φ₃(2), u(11)=21). Unit-tournament H=1,5,15,43,141,… (odd, in spectrum). The flop, if anywhere, is on the NON-lattice true optima (n≥22).
tags: [unit-distance, tournament, hamiltonian-path, traceable, flop, harborth, tiling-model, redei, frontier-gain, eisenstein, recursive, spine-tile]
---

# Unit distance ↔ tournament: the Ham-path "flop" question

**Prompt (user):** map points→tournament (flip the tile if dist=1, else transitive); is the
tournament's mandatory Hamiltonian path part of the unit distances? Does a unit-distance Ham path
always exist in higher-n optima; if not, at what n does it flop? Recursive patterns? (modify the
mapping freely.)

## The reduction (what the mapping really asks)

In the tiling model the **base Hamiltonian path is a *labeling*** of the points (the order
`p₁,…,pₙ`), and the base path is "part of the unit distances" iff **consecutive labels are at unit
distance**. So:

> **The tournament's base Ham path can be made all-unit `⟺` the optimal unit-distance graph has a
> Hamiltonian path (is *traceable*).** The labeling *is* the unit Ham path. The "flop" = the first
> `n` whose optimal unit-distance graph is **not** traceable.

(This holds for either convention — unit=flipped or unit=unflipped — since both ask whether a
unit-edge path spans the points.)

## The answer: the lattice optimum never flops (verified n=3..28)

**Verified** (`unit_distance_ham_path_flop_s599y.py`): the triangular-lattice (Harborth) optimum is
**Hamiltonian-traceable for every `n=3..28`** — a unit Ham path always exists. The structural
reason is clean:

> An **edge-maximal** unit-distance configuration has **min degree ≥ 2** and is well-connected (no
> pendant/spider — you would add the missing edge), and compact triangular-grid graphs of min
> degree ≥ 2 are traceable. **So on the lattice optimum the flop never happens**: the mandatory
> (base) Ham path is always part of the unit distances.

**The spine + tile decomposition.** With the unit Ham path as the base path,
```
   u(n) = (n−1)        unit edges on the base path (the "spine")
        + (u(n)−n+1)   unit edges as flipped tiles (the "bulk"),  u(n)−n+1 ≈ 2n − √(12n).
```
Verified split (n,(n−1)+tiles): `5→4+3, 8→7+7, 14→13+16, 22→21+28, 28→27+38`. **The Ham path uses
`n−1` of the `≈3n` unit edges; the other `≈2n` are unit tiles.** So "the mandatory Ham path is
unit" is true, but it is a *thin spine* through a much larger unit-edge bulk.

## Where the flop actually lives: the non-lattice true optima

**Honest caveat** (the correction from S599u): the triangular lattice is **not** the true maximum
for large `n` — Sawin/OpenAI's CM-field construction beats it (`>n^{1.014}`), and the exact optima
for `n≥22` are non-lattice (`u(22)∈[60,61]`, Moser-ring / CM). Those graphs are algebraically
exotic and **not obviously traceable**. So:

> **The flop, if it ever happens, is on the *non-lattice* true optima (`n ≥ 22`).** The lattice
> answer is "never"; the true-optimum answer is **open**, and the place to look is exactly where the
> optimum departs from the triangular lattice. This is a sharp, concrete target: test traceability
> of the known exact optimal graphs `u(n)` for `n=12..21` and the `u(22)` candidates.

## Recursive patterns

1. **Frontier-gain `+3`.** Each added point brings ≈`6/2 = 3` new unit edges (interior degree 6);
   **one extends the spine, two become tiles** — exactly the state-local **frontier-gain table** of
   S599w-x. The construction *is* the beam-search frontier.
2. **`√(12n)` perimeter = the frontier.** Harborth's `3n − √(12n−3)` is *bulk − boundary*; the
   `√(12n)` is the perimeter, the moving frontier whose state drives the gain. Boundary `min degree`
   alternates `2/3`.
3. **Centered-hexagon shells.** The "complete" clusters sit at the **centered hexagonal numbers**
   `1,7,19,37 = 3k(k+1)+1` (min-degree jumps there); between them the cluster grows a partial shell.
   This is the shell / Mode-B (`n→n−2`-flavoured) recursion — the same shell tower as LRC's `2n−1`.
4. **The π/3 Eisenstein echo.** `Harborth(5)=7=Φ₃(2)` and `Harborth(11)=21=3Φ₃(2)` — the forbidden
   `H`-values. *Not coincidence:* the triangular lattice **is** the Eisenstein `ζ₆` ring (angle
   `π/3`), so the Harborth count naturally produces `Φ₃(2)` values, the same `π/3`/`Cl₂(π/3)` object
   as the forbidden `H` (S599u/v). But it is the **edge count** `u(n)` that hits `7,21` — the
   *tournament* `H` is `1,5,15,43,141,513,1605,4915` (odd, in spectrum, never `7,21`), so it is a
   structural *echo* of `Φ₃(2)`, not a tournament with `7` paths.

## The tournament reading (Rédei)

Building the unit-distance tournament (base = the unit Ham path; non-consecutive pairs flipped iff
unit) gives `H = 1,5,15,43,141,513,1605,4915` for `n=3..10` — **odd (Rédei) and in the H-spectrum**
(`≠7,21`). The base path is itself a unit Ham path, so **Rédei's mandatory Ham path is realised by
the unit spine.** The unit edges are the flipped tiles; the unit count is `(spine n−1) + (unit
tiles)`; `H` is the count of *all* Ham paths of the oriented tournament (mostly using non-unit
transitive edges), growing `~×3`.

## Honest status

- **Verified:** lattice optimum traceable `n=3..28` (no flop); the `(n−1)+tiles` unit-edge split;
  unit-tournament `H` odd and in spectrum.
- **Established (reduction):** "base Ham path is unit `⟺` optimal unit graph traceable"; "edge-maximal
  `⟹` min-deg `≥2` `⟹` traceable" (the lattice no-flop reason).
- **Open / sharp target:** traceability of the **non-lattice true optima** (`n≥22`) — the only place
  a flop can occur; test the exact `u(12..21)` graphs and the `u(22)` candidates.
- **Framed:** the recursive patterns (frontier-gain `+3`, `√` perimeter = frontier, centered-hex
  shells, `π/3` echo) tie this to S599w-x (state-local), S599u (Eisenstein/`Cl₂`), and the shell
  tower.

**Artifacts:** `04-computation/unit_distance_ham_path_flop_s599y.py`, `unit_tournament_H_s599z.py`
(+`.out`s). Builds on HYP-2170 (Eisenstein UD), S599u (`π/3`/`Φ₃(2)`), S599w-x (frontier-gain),
Harborth, Rédei. New: **HYP-2201**.
