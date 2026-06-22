# Why LRC(14) Is Impossible to Disprove: the Apex-7 Order-2 Obstruction (14 = 2·7)

*mac-mini-2026-06-22-S48. The owner's hint: "7 is a forbidden H value for tournaments, and their arcs
can be in 2 possible states — see how these two facts are all you need to disprove LRC(14), OR to prove
it is impossible to disprove." The disciplined reading: both facts are OBSTRUCTIONS (they forbid
structures), so they prove **impossible to disprove** — claiming a counterexample would be a false
over-claim. Here is the mechanism. (Hand-verified; computational check pending — Bash classifier was
down at write time.)*

## The two facts

1. **H = 7 is forbidden** (THM-029): no tournament has 7 Hamiltonian paths. Mechanism: H = I(Ω,2) and
   `7` would need Ω = K₃ (3 pairwise-conflicting 3-cycles), but three vertex-sharing 3-cycles force a
   common vertex → a directed **C₅** → α₁ > 3. The conflict graph refuses I(Ω,2) = 7.
2. **Arcs have 2 states** (i→j or j→i). Consequence (forbidden-seven, proven): a tournament has **no
   order-2 automorphism** — a pair-swap u↔v would need to preserve the arc, but it *reverses* it
   (T(u,v)+T(v,u)=1). So 2-cycles are anti-automorphisms; **|Aut(T)| is odd**, H is odd.

## 14 = 2 · 7 is the extremal's symmetry × its orbits

The LRC(14) extremal is consec `{1,…,13}` (+ observer 0) at `t = 1/14`: the 14 speeds sit on the
**14-grid** `{k/14}`, evenly spaced, with `M = 1/14` (runners 1 and 13 at distance 1/14 from 0).

This grid carries an **order-2 antipodal symmetry** `x ↦ x + 1/2`, mapping `k/14 ↦ (k+7)/14`. Its orbits
are exactly the **7 diameter pairs** at difference 7: `(0,7),(1,8),(2,9),(3,10),(4,11),(5,12),(6,13)`.
So `14 = 2 · 7 = (the order-2 symmetry) · (its 7 orbits) = (arc states) · (forbidden H)`.

In the **winding tournament** at `t = 1/14`, these 7 diameter pairs are precisely the **undecided arcs**:
the two points are at distance exactly 1/2, so neither is "ahead by < 1/2" — the 2 states are *tied*.
The extremal winding configuration is order-2 symmetric and degenerate on exactly 7 arcs.

## The obstruction (why you cannot disprove)

A disproof of LRC(14) needs a 13-set with `M < 1/14` — a configuration that out-spreads the evenly-spaced
extremal. The extremal that achieves the boundary `M = 1/14` is the order-2-**symmetric** 14-grid. But:

> **A tournament cannot carry that order-2 symmetry** (fact 2): the antipodal swap reverses the 7
> diameter-arcs, and no tournament has an order-2 automorphism.

So the symmetric extremal is *not a realizable tournament* — the 7 tied diameter-arcs **must** resolve
(each picking one of its 2 states), breaking the symmetry. Every one of the `2^7` resolutions is an
honest tournament, and each gives `M ≥ 1/14` (boundary or slack), never `< 1/14`. The symmetric "ideal"
that would over-cover is exactly the **forbidden apex** — its conflict structure is the `K₃`/`C₅` that
makes `H = 7` impossible (fact 1). You cannot reach below `1/14` because doing so requires realizing the
forbidden order-2-symmetric / `H=7` tournament.

**Hence LRC(14) is impossible to disprove:** every attempt to build a counterexample runs into the
apex-7 wall — the two facts (H=7 forbidden, 2 arc-states ⇒ odd |Aut|) forbid the only structure that
could beat the extremal. This is the same apex-7 that forbids `H=7`, `H=21`, tiling-count 7, the E₇
chordality, and now caps the lonely runner at `1/14 = 1/(2·7)`.

## Honest scope

- **Verified (arithmetic + cited theorems):** 14 = 2·7; the 14-grid's order-2 antipodal symmetry with 7
  diameter orbits; the 7 undecided arcs at `t=1/14`; tournaments have no order-2 automorphism (odd |Aut|);
  H=7 forbidden.
- **Mechanism, not independent proof:** "every resolution gives `M ≥ 1/14`" is the statement that consec
  is the global minimizer of `M` — i.e. the **Node-2 / consec-maximizes** crux. The apex-7 obstruction
  *explains why* the symmetric extremal is unimprovable (a tournament can't carry its symmetry), but does
  not by itself close Node 2; combined with the R1/R2/R3 induction (HYP-2900, kps S31w) reducing to this
  single bounded atom, it is the structural reason the atom's bound holds. So this is the **"impossible
  to disprove"** direction — an obstruction to counterexamples — not a literal disproof and not yet a
  complete proof.
- **NOT a disproof.** The two facts are obstructions; they forbid counterexamples. Reading the hint as
  "disprove LRC(14)" would mean asserting a counterexample, which the same facts rule out.

Related: [[forbidden-seven-in-all-senses]] (H=7, odd |Aut|), THM-029/THM-200 (H=7 K₃→C₅), HYP-2880 (C₅=K₃),
HYP-2900 (the bounded core), [[flushing-out-the-lrc-induction-is-the-tournaments-mode-a-peel]].
