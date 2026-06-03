---
source: claudebox-2026-06-03-S581
status: REFLECTION + RESULT — a bounded CRT residue automaton for the large-owner residual; one
  proved resonance bound (generalizes Lemma C); one honest sharpening of the open problem.
tags: [LRC, n14, owners, residual, automaton, CRT, resonance, Lemma-C, two-clock, THM-398, HYP-2105]
---

# Two clocks on a circle: the large-owner residual as a CRT residue automaton

**Prompt (user):** think about the idea that "large-owner residual becomes a bounded CRT residue
automaton," and be creative with the automaton design.

## The picture that makes it an automaton

opus-S574 left the n=14 residual as "a bounded CRT/Diophantine feasibility on the large owners +
w, verified never satisfiable, not proved." The word *bounded* is the invitation: a bounded
feasibility is a **finite-state machine**, and finite-state machines have something an enumeration
doesn't — an *emptiness theorem*. So the task is to find the machine.

The endpoint-owner congruences read like two clocks. Endpoint `a` has a hand at phase
`wA mod u_a` (`A = k_a n+1`); endpoint `b` a hand at `wB mod u_b`. As the dial `w` ticks, each
hand advances by a fixed step. "Component fits an arc" means: at some tick, **both hands are
simultaneously within a `1/n`-window of the same arc spoke `j`**. That is exactly a two-clock
synchronization automaton — the joint phase `(wA mod u_a, wB mod u_b)` walks a single cyclic
orbit on `ℤ/u_a × ℤ/u_b`, and we ask whether the orbit ever lands in the accept set. The orbit's
period — bounded by `lcm(u_a,u_b)` and **CRT-factoring prime by prime** (`u_a=21,u_b=35 ⇒
P=105=3·5·7`) — is the boundedness, made into a state count.

But the accept condition hides a subtlety the design has to respect: "same spoke `j`" depends on
the *actual* `w`, which grows without bound, while the phase state is bounded. The resolution is
the prettiest part: at a window-hit state the two reps `(r_a,r_b)` are fixed, and the spoke-match
**pins** `w* = (u_b r_a − u_a r_b)/D` — a single value — so the bounded state already *computes*
its own witness and checks whether that witness lives in its residue class. The unbounded `w`
never has to be walked; eliminating `j` collapsed the infinite search to one division. `D =
u_b(k_a n+1) − u_a(k_b n−1)` — the **cross-relation defect** from Lemma C — is the modulus of that
division. The whole machine is: *advance two coupled counters; at each coincidence, divide by the
cross-relation defect and check the remainder.*

## What the machine sees: the resonance band

Because `w ≥ 1`, the pin `w*·D = u_b r_a − u_a r_b` with small reps forces

```
|D| ≤ u_b K_a + u_a K_b   (< 2 u_a u_b / n),   K_u = ⌊(u−1)/n⌋.
```

Verified: every one of 699 feasible components obeys it (0 violations); the band holds only 2.5%
of the grid. This is the certificate the bounded automaton emits — *feasibility forces the
cross-relation defect to be small*. And it is exactly **Lemma C, generalized**. Lemma C is the
corner `K_a = K_b = 0` (both owners small): the band collapses to `D = 0`, the cross-relation
`u_b(k_a n+1) = u_a(k_b n−1)`, which is `a = b` — loose. Large owners widen the band by precisely
the window slack `u_b K_a + u_a K_b`; that, and nothing else, is what "slack appears only for
`u ≥ n`" *means* arithmetically. The deepest single fit, `D = 0` with `gcd(u_a,u_b) ≥ n`, was
feasible every time (19/19) — the resonance is sharpest exactly on the cross-relation Lemma C
rules out for small owners.

## The honest turn: bounded is not empty

I expected the machine to reject everything (opus saw "never satisfiable") and instead it
**accepts 1590** endpoint-valid, short, large-owner components in the n=14 range —
`(15,2,20,3)` at `w=1` among them. The owner-congruence system, in isolation, is *satisfiable*.
That is not a contradiction with opus, who verified over *actual residual configs*; it is a
correction to the temptation to think the congruence system *alone* carries the looseness. It
does not. The looseness must live in the **config layer** — which speeds actually co-occur in a
real `S'`, the requirement that *every* component fit at once, the validity of the endpoints as
turned-safe runners. The bounded automaton is the right object, but it is only half of one: the
proof of the residual is `accept(owner-automaton) ∩ valid(config-automaton) = ∅`, an
intersection-emptiness, and I have only built the first factor. Naming the missing factor
precisely is the result — the residual is not closed by the owner clocks; it is closed by the
owner clocks *constrained to the configs that can actually be built*.

## The transcending pattern

A recurring shape in this project: a hard case localizes to a **defect that must be small**.
S559's apex was the zero-divisor; S579's apex was the non-transverse covector; here it is the
cross-relation defect `D`, and "feasible" means `|D|` fits inside a window set by the owners. Each
time, the obstruction is a near-resonance — two arithmetic clocks that must *almost* agree — and
each time the fix is to count how much slack there is and show the valid inputs never have enough.
The owner automaton measures the slack exactly (`u_b K_a + u_a K_b`). What remains is the same
move as always: show the *real* inputs live outside the resonance band. The clocks are built; the
calendar of which clocks are real is the next machine.

**Artifacts:** `04-computation/lrc_largeowner_crt_automaton_s581.py` (+`.out`); new **HYP-2110**.
Builds on opus-S574 / HYP-2105 / THM-398 §4.5 (owners, Lemma C, cross-relation), HYP-2107 (the
CRT/transversality thread).
