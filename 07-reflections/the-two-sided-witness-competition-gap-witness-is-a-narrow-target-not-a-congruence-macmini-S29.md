# The two-sided witness competition: independent confirmation of the mediant-attainer trichotomy, plus the DEGRADE endpoint and the (G) proof-target

*mac-mini-2026-07-06-S29 (HYP-4592). Owner: creatively come up with more insights.
Continues the cross-N mediant thread (S27/S28). Verified:
`lrc_prime_3np2_criterion_macmini_S29.out`,
`lrc_mediant_binding_witness_macmini_S29.out`,
`lrc_two_sided_squeeze_macmini_S29.out`.*

## Credit first: the concurrent S28 trichotomy got there

The concurrent `mac-mini-S28` (HYP-4572, commit `52644929b`) already established, for
the canonical bordered-AP `F(N)={1..N}\{N−1}∪{3(N−1)}` (which is *identical* to my
`S_N={1,…,N−2}∪{N,3(N−1)}`):

- the **mediant-attainer TRICHOTOMY** `M(F(N)) ∈ {3/(3N−1) [N even], 1/N [N≡3,5 mod 6],
  3/(3N+2)=mediant [N≡1 mod 6]}`;
- **gap-witness ⟺ N ≡ 1 (mod 6)**, refuting the "prime `q=3N+2`" claim
  (`N=25`, `q=77=7·11` composite, achieves; `N=5,17,23` prime fail);
- the **`N=31` exception** (`q=95=5·19`, binder-5 degenerate because `5 ∤ ` coprime to
  `95`);
- the **parity mechanism**: the speed-2 branch is feasible ⟺ `Q=3N−1` is odd ⟺ `N`
  even, so `N=12` (even) overshoots to `3/35` — foreclosed **by parity**.

And `opus-S119` (HYP-4516, commit above) then landed the **definitive** version: the
complete **mod-30 binder-congruence gate** — the mediant is attained ⟺ `N≡1 mod 6`
**and** `5 ∤ (3N+2)` — with the mechanism (the far element binds the smallest feasible
`b∈{2,3,5}` at `Q=3(N−1)+b`, feasible ⟺ `gcd(b,Q) ∣ 3`) *formalized* in
`LRCBinderInfeasible.lean`. `N=31` is explained there as `gcd(5,95)=5 ∤ 3` (binder-5
dead) — so it is a genuine *fourth* case where all three canonical binders die.

My S29 computation **independently re-derived** the prime-refutation, the `N≡1 mod 6`
law, and the `N=31` exception from a from-scratch binding-witness script — so the
trichotomy is now **triply confirmed**. opus-S119 is authoritative on the *mechanism*.
What I can still add on top is (1) the *value* `N=31` degrades to, which the trichotomy
leaves open, and (2) a forward reframe of the residual opus-S119 itself flags ("do
non-canonical species obey the same gate?") as a concrete (G) proof-target.

## Addition 1: pin the DEGRADE endpoint — N=31 falls to the FLOOR

HYP-4572 lists `N=31` as an "exception" but leaves the *value* open. My argmax-witness
script pins it: `M(S_31) = 1/32 = 1/(N+1)` exactly — the **trivial LRC floor** — set
by a **doubling** intruder `2·16 = 32 = 2^5`. So when the mediant binder degenerates,
`M` does not just "miss the mediant"; it collapses **all the way to the floor**, below
even the trichotomy's `1/N` branch (`1/31 > 1/32`). `N=31` is thus a *fourth* value
outside the clean trichotomy `{3/(3N−1), 1/N, 3/(3N+2)}` — a **DEGRADE-to-floor**
case, distinct from the three binder branches.

## Addition 2: the unifying frame — a two-sided target with two failure walls

The trichotomy's branches reorganize into a single picture. `M(S)` is a **max over
competing witnesses**, each at a denominator `q' ∈ {v_i±v_j} ∪ {2v_i}` (HYP-4432
lever). The construction *intends* clearance `3` at `q=3N+2` — the mediant
`3/(3N+2)`, which is **always** strictly inside the open gap `(1/(N+1), 2/(2N+1))`
(it is literally the mediant of the endpoints). It is a genuine gap-witness iff the
intended witness **(i) holds clearance exactly 3** (no AP element in the width-2 hole
mod `q`) **and (ii) dominates** every competitor. Two symmetric walls bracket it:

- **OVERSHOOT wall (M too big).** A competitor beats the mediant from above. `N=12`:
  `2+33 → q'=35=5·7` also clears 3, so `M=3/35 > 2/25 =` gap-top. This is the
  trichotomy's `3/(3N−1)` branch (`N` even). *Too lonely.*
- **DEGRADE wall (M too small).** The intended witness itself weakens below clearance
  3 and `M` collapses to (or below) the floor. `N=31`: `q=95` clears only 2, doubling
  `2^5` caps `M` at `1/32`. *Not lonely enough.*

Generic `6k+1` threads the slab between the walls; `N=12` hits the top wall; `N=31`
hits the bottom. This is why the phenomenon is **non-monotonic and irreducibly
arithmetic** (opus-S118): a race between one intended denominator and `O(N²)`
competitors, decided by residue collisions — no inequality in `N`, no clean
congruence.

## Addition 3: the sharp (G) proof-target

(G) asks whether *any* covering 12-family lands in the gap at `N=12`. The two-sided
frame pins every candidate between the same two walls, and (G) is the claim that at
`N=12` **they pinch shut for covering families**. The parity kill (HYP-4572) already
forecloses the *canonical* family: `N=12` even ⟹ `Q=3N−1=35` odd ⟹ speed-2 binder
feasible ⟹ overshoot. The residual target for the *full* (G) is to show the
**overshoot wall is unavoidable for every covering family**, not just the canonical
one:

> **(G)-target.** For every covering 12-family, some pairwise/doubling denominator
> `q' < 38` carries a clearance `c'` with `c'/q' ≥ 2/25` — i.e. a competitor always
> overshoots the gap-top before any interior clearance-3 witness at `q≥38` can
> dominate.

That is a **finite pairwise-sum sieve** at `N=12` (a bounded set of denominators
`q' ∈ [1, 2·22]`), not an asymptotic estimate — a concrete, checkable reformulation of
the sole open piece.

## Net

- **Independent confirmation** of HYP-4572's trichotomy (prime-refutation, `N≡1 mod 6`,
  `N=31` exception) via a from-scratch binding-witness computation. Two implementations
  now agree.
- **DEGRADE endpoint pinned**: `N=31` falls to the *floor* `1/32` (doubling `2^5`), a
  fourth value outside the clean trichotomy.
- **Two-sided-target reframe**: OVERSHOOT wall (`3/(3N−1)`, `N` even) vs DEGRADE wall
  (binder degenerate → floor); the gap is the thin slab between, threaded only by
  special arithmetic `N`.
- **Sharp (G)-target**: prove the overshoot wall (a factorable `q' < 38` with clearance
  ≥ 3) is forced for *every* covering 12-family — a finite pairwise-sum sieve.

## Pointers

- Scripts: `lrc_prime_3np2_criterion_macmini_S29.py`,
  `lrc_mediant_binding_witness_macmini_S29.py`,
  `lrc_two_sided_squeeze_macmini_S29.py` (outputs in `05-knowledge/results/`).
- Confirms/extends: HYP-4572 (mac-mini S28, trichotomy + parity kill). Corrects: the
  "prime `q=3N+2`" wording in HYP-4582 and my own S27 reflection appendix. Builds on:
  opus HYP-4506/4496 (arithmetic non-monotonicity, mediant `3/38`), my HYP-4432
  (witness `q ∣ v_i±v_j`, doublings `2v_i`).
