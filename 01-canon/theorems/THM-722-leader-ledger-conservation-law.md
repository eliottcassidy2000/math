---
id: THM-722
title: The Leader-Ledger Conservation Law (with the chain partition, the handoff parity law, and the climb lemma LEM-025)
status: PROVED (elementary; exact computational verification on 50+ families, 04-computation/lrc14_leader_ledger_boxeph_S21.py)
source: boxeph-2026-07-12-S21
related: [THM-668, THM-373, THM-381, HYP-6240, HYP-6250, HYP-6280, LEM-025, MISTAKE-141]
---

# THM-722 — The Leader-Ledger Conservation Law

**Setting.** Let S = {v_1 < … < v_n} be distinct positive integers, **not all odd** (so
M(S) < 1/2; every covering family qualifies since it has an even element). On the time
circle R/Z define

- the **signed phase** p_i(t) = v_i t − round(v_i t) ∈ (−1/2, 1/2],
- the **distance-to-loneliness** f(t) = min_i |p_i(t)| (so M(S) = max_t f(t)),
- the **leader** λ(t) = argmin_i |p_i(t)| (ties resolved lexicographically; ties occupy a
  finite set of t),
- the **leader phase** φ(t) = p_{λ(t)}(t), so |φ(t)| = f(t) ≤ M(S) < 1/2.

**(A) Structure of φ.** φ is piecewise linear with slopes v_{λ(t)} > 0. Every leader
change happens at a **pair-ruler point**: a *difference event* ((v_j − v_i)t ∈ Z, where
p_i = p_j; φ is continuous there) or a *sum event* ((v_i + v_j)t ∈ Z, where p_i = −p_j;
φ jumps). At a jump the transition is always a FALL: φ(t⁻) = +x → φ(t⁺) = −x with
x = f(t) ≥ 0. Call a sum event with x > 0 where the leadership actually changes a
**sum-handoff**; let H⁺ = #sum-handoffs per period and x_h = f(t_h) its **depth**.

**(B) The conservation law.**

    ∫₀¹ v_{λ(t)} dt  =  2 · Σ_{sum-handoffs h} x_h .

**(C) The chain partition.** The H⁺ sum-handoffs cut the circle into H⁺ **chains**. On
each chain φ rises strictly (through difference-handoffs, continuously) from −x_in to
+x_out, crossing 0 exactly once — at a **leader landing** (a time k/v where the leader v
passes its own lattice point). Per chain,

    ∫_chain v_λ dt  =  x_in + x_out  ≤  2 M(S) ,

and within a chain the leader speed is unimodal: difference-handoffs pass leadership
slow→fast while φ < 0 and fast→slow while φ > 0 (the chain "climbs the speed ladder,
lands, descends").

**(D) Consequences.**
1. **M(S) = max_h x_h** — the maximum of f is attained at a sum-handoff, so the witness
   denominator divides a pair sum v_i + v_j (Kravitz's reduction / THM-668-mac-mini's
   pair-sum ruler, rederived dynamically).
2. **M(S) ≥ (∫₀¹ v_λ dt)/(2 H⁺)** — the maximum depth is at least the average depth. The
   **ledger efficiency** η(S) = (∫v_λ)/(2 H⁺ M) ∈ (0, 1] measures equioscillation: η = 1
   iff every handoff has depth exactly M (attained by {1,2}; AP {1..13} has η ≈ 0.785).
3. **H⁺ = #leader landings** (one zero-crossing per chain): of the Σ_i v_i lattice
   passages per period, exactly H⁺ are made while leading.

**(E) The parity law.** If S contains an even element then **H⁺ is EVEN**. (The
involution ι: t ↦ −t maps chains to chains reversing orientation; its fixed points 0 and
1/2 each lie in the INTERIOR of a chain — at t = 0 all runners land so the leader crosses
zero there; at t = 1/2 the even runners land and the slowest even leads — so exactly two
chains are ι-fixed and the rest pair up. H⁺ ≡ 0 (mod 2).) This is the ledger home of the
ι-pairing of witnesses (klein-S270's even-witness-count constraint): the deepest handoff
at +t* pairs with its mirror at −t*.

## Proof

Everything follows from three observations about the finite piecewise-linear system.

**(1) Leader changes are pair-ruler events.** The functions |p_i| are piecewise linear;
a change of argmin at t₀ forces a tie |p_i(t₀)| = |p_j(t₀)| between the old and the new
leader, i.e. p_i = p_j (so (v_j − v_i)t₀ ∈ Z) or p_i = −p_j (so (v_i + v_j)t₀ ∈ Z). At a
multi-way tie all tied phases have the same modulus x, hence lie in {+x, −x}, and the
old→new transition is still of one of the two types. The leader's own phase never wraps:
|φ| = f ≤ M < 1/2 (this is where "not all odd" is used: an all-odd family has f(1/2) = 1/2).
So on a maximal interval of constant leader, φ = p_λ is linear with slope v_λ, and the
only possible discontinuities of φ are at leader changes with p_old = −p_new ≠ 0.

**(2) Sum-jumps fall.** Suppose at t₀ leadership passes i → j with p_i(t₀) = −p_j(t₀) =
−x < 0 (an upward jump). Just before t₀, |p_i| = x + v_i δ and |p_j| = x − v_j δ < |p_i|
(δ = t₀ − t > 0), contradicting i's leadership. So every jump is +x → −x: DOWN by 2x.
The same one-sided comparison shows difference-handoffs pass slow→fast when the common
phase is negative and fast→slow when positive (part C's unimodality): below zero the
faster runner's |p| shrinks faster, above zero the slower runner's |p| grows slower.

**(3) The ledger balances.** φ has period 1, is piecewise linear with positive slopes
v_λ, and its only discontinuities are downward jumps of size 2x_h. Over one period the
net change is zero, so total rise = total fall:

    ∫₀¹ φ'(t) dt = ∫₀¹ v_λ(t) dt   and   total fall = Σ_h 2x_h .            ∎ (B)

Between consecutive downward jumps φ rises strictly from −x_in to +x_out, hence crosses 0
exactly once; at the crossing p_λ = 0, i.e. the leader sits on its own lattice point —
a leader landing. Integrating the slope across one chain gives x_in + x_out, and each
x ≤ max f = M. ∎ (C, D)

For (E): ι(t) = −t satisfies p_i(−t) = −p_i(t), so f and the leader are ι-invariant and
φ(−t) = −φ(t); ι maps sum-handoffs to sum-handoffs with the same depth and maps the chain
interval (h, h′) to the chain (ι h′, ι h). A chain is ι-fixed iff it contains a fixed
point of ι (0 or 1/2) in its interior. t = 0 is interior to a chain: every runner lands
at 0, f(0) = 0, and no positive-depth jump occurs there (the leader v_min crosses zero
continuously). If S has an even element, the same holds at t = 1/2: the even runners land
(f(1/2) = 0), the slowest even runner leads through it, and 1/2 is not a positive-depth
handoff. Two ι-fixed chains; all others in disjoint pairs; hence H⁺ even. ∎ (E)

For (D1): f is continuous, piecewise linear with slopes ±v; its local maxima occur where
the active slope turns + → −, i.e. exactly at the downward jumps of φ (a difference-
handoff keeps φ continuous and rising). So max f = max_h x_h. ∎

## LEM-025 — the climb lemma (stopping-time form of the far-element floor)

**Statement.** Let S be distinct positive integers, v = min S, f = max S, B = max(S∖{f}),
q = v + f. Then for every integer k with 1 ≤ k ≤ q/(B+v),

    f(k/q) = v·k/q   exactly;   hence   M(S) ≥ v·⌊q/(B+v)⌋/q .

**Proof.** At t = k/q: for u ∈ S∖{f}, uk ≥ vk and uk ≤ Bk ≤ q − vk (since (B+v)k ≤ q),
so uk mod q = uk ∈ [vk, q − vk] and ||u k/q|| ≥ vk/q, with equality at u = v. For f:
fk = (q − v)k ≡ −vk (mod q), so ||f k/q|| = ||v k/q|| = vk/q (vk ≤ q/2 since B ≥ v). ∎

**The four extremals are all tight instances at the (min,max) ruler** (verified exactly):

| family | v, B, q | ⌊q/(B+v)⌋ | bound | true M |
|---|---|---|---|---|
| AP {1..13} | 1, 12, 14 | 1 | 1/14 | **1/14** ✓ tight |
| ladder {1..12, 13k} | 1, 12, 13k+1 | k | k/(13k+1) | **k/(13k+1)** ✓ tight |
| deep well {1..12, 182} | 1, 12, 183 | 14 | 14/183 | **14/183** ✓ tight |
| compressed 2·{1..12}∪{13} | 2, 22, 26 | 1 | 1/13 | **1/13** ✓ tight |

In ledger language: the (v_min, v_max) sum-ruler carries a staircase of handoffs at
t = k/q with depths climbing linearly (vk/q) until the **stopping time** k* = ⌊q/(B+v)⌋,
where the second-largest runner's lander cuts the climb (for the deep well: the runner 12
lander at 183 − 12k < k ⟺ k > 183/13, so k* = 14 — this is mac-mini cont.56's "omit the
distance-1 lander" as a stopping time, and the +1 of Φ₆ = 13·14 + 1 is the ruler offset).
Verified: the deep well's first 14 sum-handoffs are exactly t = k/183, pair (1,182),
depth k/183, k = 1..14.

**Two-line closure of the covering {1..12}∪{f} stratum** (statement known from mac-mini
cont.55/56; this proof is citation-free): covering forces 13 | f and 14 | f, so 182 | f,
q = f + 1 ≡ 1 (mod 13), ⌊q/13⌋ = f/13, and LEM-025 gives M ≥ (f/13)/(f+1) ≥ 14/183 for
f ≥ 182 — no LRC(≤13) citation, no exact-M computation. (If 13 ∤ f the family is
non-covering and t = 1/13 clears it.)

## Scope and honest limits

- The conservation law is an identity, not an inequality engine by itself: the crux
  (covering ⟹ M ≥ 14/183) is NOT implied. The bound (D2) is weak when many handoffs are
  shallow.
- The conjecture "covering ⟹ max-chain-mass ≥ 28/183" is **REFUTED** (this session): the
  deep well itself has max chain mass 2639/17751 < 28/183 = 2716/17751 — its witness
  chain does not equioscillate (ratio 0.9716). Recorded in HYP-6280.
- LEM-025 only sees the (min, max) ruler; families whose witness lives on another
  pair-sum ruler (e.g. GW {1..11,13,24}, witness 3/14 on the (1,13) ruler) are not tight
  for it. The bound degrades to v/(B+v) as f → ∞ — the classical far-element shape.

*Files: 04-computation/lrc14_leader_ledger_boxeph_S21.py,
05-knowledge/results/lrc14_leader_ledger_boxeph_S21.out. Reflection:
07-reflections/the-leader-ledger-the-metric-lives-on-the-walls-boxeph-S21.md.*
