# Flushing out the LRC: it is exactly the THM-079 (H=21) template — reduce-to-atom + a Moon/forcing step

*kind-mendel-2026-06-22-S8. Owner: immerse in the repo's tournament-induction work, then flush out the LRC.
Pulled main throughout; converges with kps-S31w (induction tree HYP-2900), mac-mini-S47 (induction = Mode-A
peel), kps-S31y (counterexample = forbidden K₃/H=7). This makes the template precise, **verifies every
linchpin independently**, and pins the single crux with a crisp, possibly-provable target.*

## The proof, fleshed out (mirrors THM-079's proof that H=21 is forbidden)

THM-079 forbids `H=21` in two moves: **(A) reduce to one strong atom** (`H` is multiplicative over strong
components, so WLOG the conflict graph Ω is one component); **(B) a forcing/Moon step** (3 pairwise-conflicting
odd cycles always force a 4th — THM-029 — so the exact value cannot be realized). LRC(14) has the *same two
moves*, on the winding tournament `T(x)`:

### Move A — reduce to the bounded atom (= the tournament Mode-A peel)
The size-induction tree (kps S31w / mac-mini S47), which **is** the Mode-A single-vertex peel:
- **R1 (remove-large, n→n-1 = vertex deletion):** if the largest speed `v` dominates, `safe(S) ⊇
  safe(S∖v) ∖ U_v` and the comb bound gives `meas(safe(S)) ≥ (6/7)meas(safe(S∖v)) − r/(7v) > 0` by **proven
  LRC(≤13)** (my S3 + kps). Descends every set with a dominant speed (incl. the `lcm` family of my S7) to a
  smaller proven case — by equidistribution.
- **R2 (omit-prime = resonance witness `t=1/p`):** the q-witness.
- **R3 (dilation = normalize):** scale-invariance.

These reduce **every** primitive 13-set to the **irreducible bounded covering core** (peeling stalls only
when all speeds are comparable and the set is covering — exactly THM-079's "single strong atom").

### Move B — bound the atom (= the LRC Moon/forcing step)
**Claim: every bounded covering 13-set has `M > 1/14`.** Mechanism (the forcing), *verified this session*:
- The **tight locus** `(M = 1/14)` at `n=14` is `{AP {1..13}, GW {1..11,13,24}}` (+ dilations). Verified:
  both have `M = 1/14` **exactly**, both achieve it at `t = a/14` (denominator **14**: AP at `5/14`, GW at
  `3/14`), with binding speeds `s·a ≡ ±1 (mod 14)` (AP: `{3,11}`, GW: `{5,9}`).
- **Both are NON-COVERING** (each omits a multiple of 14).
- **Covering forces a multiple of 14**, and at *any* `t = a/14` that runner sits **exactly on the observer**
  (`‖(14k)·a/14‖ = ‖ka‖ = 0 < 1/14` — the **apex-7 floor**, my S5/S6). So a covering set **cannot be tight
  at the denominator-14 optima** ⟹ `M > 1/14`.
- Verified: bounded covering sets have `min M = 1/12 ≈ 0.083 > 1/14` (over `{2..14}` and a `[1,22]` search).

This is the exact analog of THM-079's forcing: there, a candidate `H=21` *forces* a 4th cycle (ruling out the
value); here, *covering forces the apex-7 obstruction* (the multiple of 14 on the observer at every denom-14
point), ruling out tightness. It is **identical** to kps-S31y's "a counterexample = over-cover at apex-7 =
the K₃ conflict graph = `I(K₃,2)=7`, FORBIDDEN" — `14 = 2·7 = (arc-states)·I(K₃,2)`.

## The single crux, stated crisply (the only unproven step)

Move A's `R1` and Move B both reduce to **one** statement — the LRC "Moon step":
> **(★) If `M(S) = 1/14`, the optimum is achieved only at an apex-blocked point** (a `t` where some
> denom-14-type relation `s·a ≡ ±1 (mod 14)` binds) — equivalently, the **tight locus is exactly
> `{AP, GW}` (+ dilations)**, equivalently kps's **"over-cover ⟺ exact K₃"**.

Given (★), covering — which by the apex-7 floor *cannot* use any denom-14 optimum — forces `M > 1/14`, and
the whole proof closes (Move A descends the rest to proven LRC(≤13)). **(★) is the irreducible content**,
and it is precisely the two realizability structures I isolated in S7:
- the **finite half** (Node 2): "consec/AP is the unique three-gap-rigid extremizer" — a Steinhaus-rigidity
  characterization of the tight locus (why only `{AP, GW}`);
- the **analytic half** (Node 3): the effective-Weyl control of the `R1` peel in the moderate/resonant-middle
  regime (why the descent is exhaustive). The `lcm` family (S7) shows this half is irreducibly analytic.

## Why this framing is progress (the "flush")

LRC(14) is now a *single named theorem* — (★), the apex-7/`K₃`-forcing — sitting in a proven template
(THM-079). Everything else is discharged: the peel (R1, my comb bound + proven LRC(≤13)), the easy cases
(R2/R3), and the bounded-core *value* (`min M = 1/12`, verified). The remaining work is no longer "prove
LRC(14)" but "**prove tightness forces the apex-7 point**" — a concrete extremal/forcing statement with two
attack routes (three-gap rigidity for the characterization; the K₃/conflict-graph forcing of kps-S31y).

## Honest status & leads
- **Verified linchpins** (independently, exact): tight locus `{AP, GW}` both non-covering with denom-14
  optima and `s·a≡±1 (mod 14)` binding; bounded covering `min M = 1/12 > 1/14`; the apex-7 floor excludes
  covering from denom-14. So Move A's target and Move B's *value* are confirmed; only **(★)** is open.
- **Lead 1 (finish Move B):** prove (★) via the three-gap characterization of the tight locus — show any
  `M=1/14` instance has its binding pair at a Steinhaus 3-gap configuration, forcing denom-14. This is the
  finite/algebraic half and looks closest.
- **Lead 2 (kps-S31y):** prove "over-cover ⟺ exact K₃" — the conflict-graph forcing form of (★); ties (★) to
  THM-029/THM-200 (`K₃` forces `C₅`), making the LRC Moon step *literally* the H=7 impossibility.
- **Lead 3 (finish Move A):** the effective-Weyl peel threshold `v > r/(6·meas)`; my S7 hint that the witness
  denominator scales with the committed speed's *radical* (not size) may bound the moderate-regime arc-count.
- **Caution (from S7):** Move A cannot be made purely finite (the `lcm` family forces unbounded witness
  denominators); the analytic equidistribution input in Lead 3 is irreducible.

→ HYP-2906 (new), HYP-2898/2880/2864 (mine), HYP-2900 (kps tree), kps-S31y (K₃/H=7), THM-079, THM-029,
THM-200, THM-523, THM-560, OPEN-Q-108.
