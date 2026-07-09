# The abstract object is a winding circle touching the medial axis — every ruler we invented was a shadow of its intrinsic contact lattice

**mac-mini-2026-07-09-S65.** Written to answer the owner's standing question directly: *what
abstract object are we fundamentally interacting with?* The answer below is not a metaphor; each
identification is a proved statement in the repo, assembled.

## 1. The object

Fix the 13 speeds `v = (v_1, …, v_13) ∈ Z^13_{>0}`. The object is the **closed winding circle**

    C_v : T → T^13,   τ ↦ (v_1 τ, …, v_13 τ)  (mod 1)

— a closed geodesic of the flat 13-torus with winding vector `v` (a higher-dimensional torus
knot). Because the speeds are integers, `C_v` is CLOSED (period 1): there is no ergodic limit to
take, no equidistribution to invoke — the entire orbit is present in one period. Its closure is
cut out exactly by the **relation lattice** `Λ = {m ∈ Z^13 : m·v = 0}` (rank 12): the curve is
`Ann(Λ)`, and the gap/tooth process is **Haar on that annihilator** — which is literally THM-639's
Theorem C on the runner side (`L_0(E)` the balanced lattice, complete invariant), i.e. the same
statement the tournament project proved for tilings.

Two dual structures stratify the ambient torus:

- the **walls** `W_l = {x_l = 0}` (the danger set of runner `l`), and
- the **medial axis** of the wall arrangement — the locus where the two nearest walls are
  equidistant: the sheets `{x_i = −x_j}` (mirror sheets) and `{x_i = ±1/2}` (peak sheets).

Then, with `m(τ) = min_l ‖v_l τ‖ = dist(C_v(τ), walls)`:

- **LRC(14)** = "`C_v` meets the closed box `[1/14, 13/14]^13`" = `max_τ dist(C_v(τ), walls) ≥ 1/14`.
- **THM-668 (this session):** the max of the wall-distance along the curve is attained where the
  curve crosses the **medial axis**, and the crossing times with sheet `(i,j)` are exactly
  `τ = p/(v_i+v_j)`: the intersection number of `C_v` with the mirror sheet `{x_i = −x_j}` is
  **`v_i + v_j`** — the pair-sum. The witness ALWAYS lives on a pair-sum ruler.
- **Schur-triple kill rule:** the sheet `(i,j)` is useless iff it lies inside a wall along the
  curve — iff `(v_i+v_j) | v_l` for some `l` (in particular a Schur triple `v_i+v_j = v_l`). This
  is the mechanism behind opus-S182's discovery that E3 (Schur triples), not additive energy,
  governs tightness: **each additive relation of height 3 deletes one contact sheet from the
  witness supply.** The AP `{1..13}` maximizes E3 = 78 → maximal deletion → the witness supply
  collapses to the single surviving sheet `q = 14` (its non-covering modulus), where the curve is
  TANGENT to the box corner (`σ = 0` exactly, verified): opus-S177's pinch/lemniscate node, in
  closed form.

So the tight extremal is not "a bad velocity set"; it is **the maximally coherent circle — the one
whose relation lattice is 1-dimensional and dilate-coherent (opus-S181), whose contact lattice is
maximally cannibalized by its own resonances, osculating the box corner.** Maximal coherence = the
AP = CM-like extra symmetry (the Φ6/Eisenstein thread, klein-S206: universal in n, only the
cushion `n²/Φ6(n) → 1` is n-dependent).

## 2. Every obstruction of the last month was a non-intrinsic shadow

Each "wall" the fleet hit was an artifact of projecting the circle to a quotient it does not
naturally live on:

| obstruction | what it actually was |
|---|---|
| drift `e_i φ/Vmax` (hembed refutations, mac-mini-S64) | evaluating the curve at ruler points `j/Vmax` instead of on itself; klein-S207: a discretisation artefact, gone at real τ |
| ruler points never lonely (klein-S207) | the observer coordinate collapses on its OWN ruler — `Vmax·(j/Vmax) ∈ Z` — the quotient kills the coordinate that defined it |
| grid-invisible pinches (MISTAKE-130, opus-S177) | sampling an open condition on a grid that provably cannot see its measure-zero boundary |
| "the witness is on a different ruler: 39 ∤ 91" (klein-S207) | the intrinsic rulers are the CONTACT moduli `v_i+v_j`, not the coordinate modulus `Vmax`; `39 | 24+54` — lawful (THM-668) |
| the Mertens wall / absolute residual bounds (kps-S98, klein-S202) | bounding the theta-sum over Λ absolutely when the signed/geometric structure of Λ (dimension + coherence, opus-S181) carries the content |
| "no bounded-q reduction" (mac-mini-S64 OtR) | true absolutely, but the RELATIVE bound `q ≤ 2Vmax` holds — the contact lattice is bounded by the curve's own winding |

The convergent prescription (kps-S112's continuum bridge, klein-S207's real-τ criterion, THM-668's
event form) is one sentence: **work on the circle itself, and when a finite check is needed, use
its intrinsic contact lattice — never an extrinsic grid.**

## 3. The two halves of the repo study the same exact sequence

The tournament project's deepest frame ("everything is the triangle") is the decomposition
`arcs = base-path (cut space) ⊕ wiggly tiles (cycle space)` on the hypercube of tilings, quotiented
by `S_n` (the metagraph). THM-639 ported it verbatim to runners: `steps = path`,
`balanced lattice L_0(E) = tiles`, `Haar on Ann(L_0) = the gap law`, palindrome/reversal = the
involution, and the AP = the "transitive tournament" of runner families (coarsest order-cell
complex, wall count `C(k+1,3)`).

So both halves interact with one abstract object:

> **A 1-parameter orbit in a compact abelian group, presented by an exact sequence
> `0 → relations (cycles/resonances) → ambient lattice → hierarchy (cut/path) → 0`,
> interrogated through finite quotient shadows, with all difficulty concentrated on the orbits of
> maximal coherence.**

- Tournament side: hierarchy = score/base-path; relations = cycle space; shadows = tilings/
  metagraph classes; maximal-coherence extremals = transitive & rotation/Paley (the "regularity is
  extremal" of kps-S13); duality = complementation.
- LRC side: hierarchy = the fast phase / observer (the Vmax-ruler direction); relations =
  `Λ = v^⊥` (Schur triples = its minimal vectors, mac-mini-S25's kissing vectors); shadows =
  rulers/grids; maximal-coherence extremal = the AP (Freiman minimal-sumset, opus-S180); duality =
  co-offset reversal `v ↔ e = Vmax − v` (THM-639's palindrome involution).
- Even the failure taxonomy matches: the tournament side's MISTAKE-033 (confusing complement
  tiling with `T^op`) and the LRC side's MISTAKE-130 (grid vs open set) are both "shadow mistaken
  for object".

The four constants of the triangle (√2, π, e, γ) and the LRC's Φ6/Eisenstein constants are the
same phenomenon at the two ends: **an orbit's contact with a symmetric arrangement produces CM
points at maximal symmetry** — lemniscate `Z[i]` on the tiling side (opus-S177), equianharmonic
`Z[ω]`/Φ6 on the runner side (klein-S206). The heptagon (7 = 14/2, non-Fermat, non-constructible)
is why the box threshold `1/7` resists coordinate constructions.

## 4. What the object prescribes (falsifiable, actionable)

1. **Dispatch on the contact lattice.** Any finite/native_decide leg should enumerate pair-sum
   events (exact integers, provably complete — THM-668 Part 2), not grids. This subsumes and
   sharpens monad-S1's "certified rational checks" recommendation.
2. **The a-priori bridge is a Schur-budget statement.** Covering ⟹ the constant rulers `q ∈
   {2..14}` all carry a zero residue; the conjecture is that the 91 pair-sum rulers cannot ALL be
   dead or band-empty. Since E3(S) ≤ E3(interval) = C(13,2) = 78 with equality iff dilated
   interval (opus-S182's target), and dilated intervals are non-covering at their own modulus, the
   covering hypothesis caps the kill budget strictly below saturation. **The open problem = show a
   sub-saturated kill budget leaves a live banded ruler.** This is a finite-geometry statement
   about residues, not an equidistribution statement.
3. **Expect the same object at every n.** Nothing in §1 uses 14 except the box size. klein-S206's
   cushion `n²/Φ6(n) → 1` says the coherent extremal presses the box corner harder as n grows —
   the object predicts LRC(n) difficulty is monotone in coherence-vs-cushion, with no new
   phenomena, only thinner margins.
4. **For the engineering mandate:** distance-to-wall maxima along integer winding circles =
   max-min scheduling on cyclic resources; THM-668's event enumeration is an exact O(Σ(v_i+v_j))
   algorithm for the max-min "just-in-time" window of k periodic processes — a clean library
   candidate alongside mod_rank (the same pair-sum events are the collision times of round-robin
   schedulers). Worth a spec in the engineering synthesis.

## 5. The one-sentence answer

**We are studying how much coherence a closed integer-winding circle can carry before it stops
meeting the middle of the torus — and every tool that worked (co-offsets, good periods, density
floors, pair-sum witnesses, Schur counts) is a coordinate on its contact with the medial axis of
the wall arrangement, while every tool that failed was a shadow of the circle on a ruler it never
lived on.**

→ THM-668, THM-639, THM-420; klein-S206/S207; kps-S112; opus-S177/S180/S181/S182; monad-S1;
mac-mini-S25/S62 (broken clock)/S64; MISTAKE-030/130; everything-is-the-triangle;
merged-metagraph-invariants.
