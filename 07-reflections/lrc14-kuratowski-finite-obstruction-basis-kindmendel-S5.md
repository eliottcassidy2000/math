# LRC(14) through the Kuratowski / Robertson–Seymour lens: a finite obstruction basis for the witness

*kind-mendel-2026-06-22-S5. Inspiration (owner): Kuratowski/Wagner ({K5,K3,3} = forbidden minors for
planarity), Robertson–Seymour (minor-closed ⇒ finite forbidden-minor set), the repo's {7,21} forbidden
H-values, tournament↔even-graph equinumerosity, the two simplifications (suppress degree-2 / contract
edge), and Beurling–Selberg/trig. Goal: many ideas, prove what we can, honest about the rest.*

## The mapping that makes the analogy bite

The repo ALREADY has a Robertson–Seymour story: **{7,21} are permanently forbidden H-values**
(THM-029/075/079), and the obstruction is literally **Kuratowski K5 = five mutually-overlapping odd
cycles** in the conflict graph Ω(T) (THM-519); "Ω(T) planar" is conjecturally characterized by a finite
forbidden set (OPEN-Q-106). Deletion–contraction `H = H\e + H/e` (THM-082) is the "contract edge"; the
Alcuin "+1" is the repo's example of a property that is *not* minor-closed under contraction (THM-520).
And **kps-S36 (HYP-2867) just proved the LRC GOOD event depends only on the complement-EVEN part of the
cluster** — so the event lives in even-graph land, where tournaments↔even-graphs are equinumerous.

This session I push the analogy onto my HYP-2864 (every covering 13-set has a bounded-denominator
lonely witness `τ=a/D`), which is the cleanest place a *finite obstruction set* could live.

## Idea 1 (TESTED, strong) — a FINITE certificate basis of denominators

A counterexample must be a *covering* set (THM-523). Over **602 primitive covering 13-sets** (scales to
10⁴, incl. the loosest known `{1..11,13,84}`):
- the **minimal witness denominator is always ≤ 41** (distribution concentrated in `D∈{15,…,25,27,41}`);
- **a basis of just THREE denominators `{83, 89, 21}` certifies ALL of them** (greedy set-cover): every
  covering set is lonely via `τ=a/D` for one `D` in the basis.

This is the Kuratowski/RS payoff in spirit: the hard case of LRC(14) is, empirically, a **finite residue
check** — "for each covering set, one of a fixed finite list of cheap rational witnesses works." The
basis isn't canonical (greedy picks permissive large primes), but its *existence at size 3* is the
content. (`04-computation/lrc14_forbidden_obstruction_kindmendel.py`.)

## Idea 2 (TESTED) — the Beurling–Selberg / character-sum count (why witnesses exist; why a basis)

Count witnesses by a character sum on `Z/D`:
> `N(S,D) = #{a coprime mod D : ‖s·a/D‖ ≥ 1/14 ∀s∈S} = Σ_{a coprime} ∏_{s} 1_safe(s a mod D)`.
Fourier-expanding `1_safe` gives **main term `(6/7)^13·φ(D)`** (all-frequencies-zero = independence),
with error a sum over **resonances `Σ_s k_s s ≡ 0 (mod D)`** — *the same lattice as the Node-3
spectrum-intersection*. Measured: for covering sets `N ≈ (6/7)^13 φ(D)` (ratio `N/main ≈ 0.5–3.4`,
fluctuating); **positive on average but individual `D` can dip to `N=0` via a resonance** — which is
*exactly why a finite basis, not a single `D`, is required*. The Beurling–Selberg minorant turns this
into a rigorous one-sided finite bound `N ≥ Σ_{|k|≤K}(…)`. (`..._character_sum_witness_kindmendel.py`.)

The tight **AP `{1..13}` (non-covering, M=1/14) has `N=0` for all `D∤14`** — confirming the clean split:
the analytically hard tight sets are non-covering (handled by the `q=14` witness `τ=1/14`), and the
covering sets are uniformly loose, so they admit bounded-`D` witnesses.

## Idea 3 (PROVED sub-fact) — why the minimal `D ≥ 15`: the apex prime 7 is the coarsest obstruction

Covering forces a multiple of `14=2·7`. At `D=14`, `‖(14k)·a/14‖ = ‖k a‖ = 0` for **every** multiplier
`a` — the forced multiple-of-14 runner sits *exactly on the observer* — so **`D=14` never certifies a
covering set** (verified). This proves the empirical floor `min-witness-D ≥ 15`: the witness must live
in the band just *above* the threshold denominator. (Mirrors why `14=2·7` composite defeats the
polynomial method: the apex prime 7 is the obstruction at the coarsest scale, resolved only by going
finer.) This is the cleanest fragment of an actual obstruction theorem.

## Idea 4 (conceptual, evocative) — {7,21}/K5 obstruction in the winding tournament

kps-S36: "LRC is a question about random round tournaments `T(x)`; observer lonely ⟺ marked vertex is a
source (THM-381); `H(T(x))` = Rédei count." The analogy: just as `H∈{7,21}` are forbidden because they'd
force a **K5 of overlapping odd cycles** (THM-519/029/079), perhaps *non-loneliness of a covering set*
forces `T(x)` (for all `x`) into a forbidden conflict-graph class — and the covering structure cannot
realize it. Untested (building `T(x)` and its Ω is heavy), but it's the most direct {7,21}↔LRC bridge:
**both loneliness and the H-spectrum are governed by the odd-cycle/overlap (conflict-graph) structure of
the same winding tournament.** Worth a dedicated session.

## Idea 5 (conceptual) — even-graph forbidden family; the n=7 odd holes

By HYP-2867, GOOD only sees the **complement-even part** of the cluster. The even-graph metagraph `E_n`
is **chordal (perfect) for n≤6 but develops ODD HOLES at n=7** (CLAUDE.md). Since `7` is the apex prime
and `n=7`/`q=7` is exactly LRC(14)'s threshold structure, the obstruction to a clean perfect-graph
argument appearing *precisely at n=7* is suggestive: **the LRC(14) hardness may be the same phenomenon as
`E_7` losing chordality** (odd holes = odd cycles = the same odd-cycle obstruction as K5/Rédei). A
concrete test: do the "bad" (low-witness, near-tight) covering clusters map to `E_7` odd-hole even-graphs?

## Idea 6 (honest negative) — loneliness is NOT minor-closed under runner-deletion

Deleting a runner *increases* M (fewer constraints), so "non-lonely at 1/14" is **not** preserved by
deletion; there are no smaller "non-lonely cores," so the naive forbidden-subset characterization fails.
The structure-preserving operation is *drop a runner AND relax 1/14→1/13* — which is exactly the
proven-LRC(≤13) induction (THM-525), hard at each step. The minor analogy must therefore live in a
*different* category (residue patterns mod D, or even-graphs), not in the speed-subset lattice. This is
the key reason the analogy is subtle, and why Ideas 1–2 (residues) are the productive home for it.

## More ideas (brainstormed; for future sessions)
- **CRT constructive witness**: build `a mod D` factor-by-factor (avoid danger mod each `p^e | D`) using
  covering's rich divisibility — a potentially *provable* route to Idea 1.
- **Second-moment / Paley–Zygmund** on `N(S,D)`: `Var = ` the resonance sum; `Var < E[N]^2 ⇒ N>0` for
  most `D` (the variance is the Node-3 spectrum again — unifies the count with decorrelation).
- **One-octave band**: minimal `D∈[15,28]` empirically — conjecture a witness always exists with
  `15≤D≤2·14`, a very clean bounded check.
- **Equinumerosity as a count**: `#tournaments=#even graphs (A000568)`; does `N(S,D)` admit an
  even-graph-counting closed form?
- **Beurling–Selberg majorant for the cap (Node 2)** dual to the minorant for the floor — both
  band-limited, both finite-frequency, the shared `√`-cancellation backbone (HYP-2862).

## Honest status

No node is *proved closed* this session. But the Kuratowski/RS lens produced two concrete, novel,
mutually-reinforcing results — a **size-3 finite certificate basis** (Idea 1) and its **character-sum
explanation with main term `(6/7)^13 φ(D)`** (Idea 2) — plus a genuinely *proved* fragment (the `D≥15`
apex-7 floor, Idea 3) and a sharp conjecture (one-octave band). The most promising path to an *elementary*
LRC(14): prove the finite-basis / bounded-`D` witness (HYP-2864/2872) via the CRT construction or the
character-sum main-term domination. The most beautiful (unproven) lead: the **n=7 odd-hole / K5 odd-cycle
obstruction** as the shared root of {7,21}, `E_7` non-chordality, and LRC(14) hardness.

→ HYP-2872 (new), HYP-2864/2847 (mine), HYP-2867 (kps complement-even), THM-029/075/079/519/520/523/525,
OPEN-Q-106/108.
