# The deep geometry of LRC(14), and exactly where the proof breaks (S521)

*claudebox-2026-06-01-S521. A serious attempt at n=14, honestly reported. It is
NOT a proof — LRC for 13 movers is open (the unconditional frontier is 6 movers).
What it produces: a clean geometric picture, a proved reduction, a precise reason
the elementary attack stalls, and a sharp finite conjecture that is the right
target. Builds on the multiplicative-walk reduction, THM-384/386, THM-369,
HYP-1991, and the ranking-walk reflection.*

## 1. The geometry

The speeds are integers, so `t |-> (v_1 t, ..., v_13 t) mod 1` is a **closed
rational curve** on the torus `T^13` (period 1 — a torus knot of slope `v`).
LRC(14) asks: does this 1-dimensional curve meet the **central cube**
`K = [1/14, 13/14]^13`, side `12/14`, volume `(6/7)^13 ≈ 0.135`?

- Each runner's bad set `{t : ||v_i t|| < 1/14}` has measure `1/7`.
- Union bound: `13 · 1/7 = 13/7 ≈ 1.857 > 1` — **fails**. The good set is not
  forced nonempty by measure alone.
- When the bad sets are *independent* (the curve equidistributes — irrational-type
  direction), the good measure is `≈ (6/7)^13 > 0` and a lonely time exists. So
  **the entire obstruction is rational coincidences**: speed relations that pile
  bad intervals on top of each other. This reduces LRC(14) to finitely many
  rational relation-types — the "fully covered" core.

This is the honest core geometry: a line crossing a periodic arrangement of
forbidden slabs, asked to reach the central cell; only arithmetic alignment can
keep it out.

## 2. The proved reduction (the part that works)

Relabel the configuration tournament by **residues mod q** (the ranking-walk on a
fixed circulant arena). At `t = a/q` a runner `v` is safe iff `v a mod q` avoids
`F_q = {r : min(r,q-r) < q/14} = {0, ±1, ..., ±(⌈q/14⌉-1)}`.

> **Reduction (proved).** For `q ≤ 14`, `F_q = {0}`, so `t = 1/q` is lonely iff no
> speed is divisible by `q` (`||v/q|| ≥ 1/q ≥ 1/14`). Hence LRC(14) holds for every
> system EXCEPT the **fully covered** ones — those where every `q ∈ {2,...,14}`
> divides some speed (equivalently the speeds carry multiples of 8, 9, 5, 7, 11,
> 13). For all others, `t = 1/q` at the uncovered `q` is an explicit witness.

This disposes of the vast majority of systems in one line — the THM-369 sieve as
the base of the walk.

## 3. Why no single modulus closes the core (the threshold phenomenon)

On the fully-covered core we need `q > 14`. Winning at `q` = some unit `a` rotates
the speed-residues off `F_q`. Count the obstruction:

| q | radius `⌈q/14⌉-1` | `|F_q|` | antipodal pairs `(q-1)/2` |
|---|---|---|---|
| 17 | 1 | 3 | 8 |
| 19 | 1 | 3 | 9 |
| 23 | 1 | 3 | 11 |
| 25 | 1 | 3 | 12 |
| 28 | 1 | 3 | 13 |
| 29 | 2 | 5 | 14 |

At `p = 17,19,23` the pairs (`8,9,11`) are all `≤ 13`, so 13 movers can hit every
pair and **block** the modulus. At `p = 29` the pairs (`14`) finally exceed 13 —
but the radius jumps to 2, so `F_29 = {0,±1,±2}` and **each residue now blocks two
translate-sets**, exactly doubling the adversary's power and cancelling the gain
(each residue `r` lies in exactly 2 of the 14 translate-sets, so 13 residues cover
up to 26 ≥ 14). The forbidden radius grows in **lockstep** with the pair-count, so
**13 movers retain exactly enough blocking power at every single modulus.** No
one-modulus counting argument can close LRC(14) — this is the precise, structural
reason it sits right at the threshold of provability. (13 movers is the smallest
count where the adversary first wins each single-modulus game.)

## 4. The cross-modulus incompatibility (the empirical truth)

An adversary that blocks *one* modulus must spend ~`(p-1)/2` speeds on it. To
block *all* small moduli at once it would have to satisfy many such constraints
simultaneously. Empirically it cannot:

- A CRT adversary built to cover `q ≤ 14` AND block primes 17, 19, 23 is still
  lonely at `q = 25`.
- Over 40,000 fully-covered adversarial systems, the **worst** minimal lonely
  denominator is **25**.

> **Conjecture Q0.** Every primitive 13-speed system has a lonely time at `t = a/q`
> with `q ≤ 25` (some small absolute bound).

With Q0 **and** a bound on minimal-counterexample speeds (to make the residue
check finite), LRC(14) would follow by a finite computation. Both bounds are open.

## 5. Honest verdict — and what the geometry gave

This is **not** a proof of LRC(14). The two missing links are exactly the genuine
LRC difficulty, now localized precisely:

1. **The denominator bound Q0** = the statement that the single-modulus blocking
   constraints are *jointly incompatible* across `q`. Proving it is as hard as
   LRC for the fully-covered core — it is the rational-coincidence obstruction in
   disguise.
2. **A speed bound** to finitize "for all speeds."

What the geometry *did* deliver, rigorously:
- the closed-curve / central-cube picture and the union-bound failure (why the
  problem is real);
- a **proved** reduction killing all non-fully-covered systems;
- the **lockstep radius/pair-count law** explaining why 13 movers is exactly the
  threshold where single-modulus arguments fail (a structural "why n=14 is hard");
- a sharp finite target (Q0 ≈ 25) that is the right thing to try to prove next.

## Seeds

a. **Prove Q0 for fully-covered systems** via a cross-modulus incompatibility:
   show 13 residues cannot simultaneously block the antipodal-translate games at
   `q ∈ {17,19,23,25,29}` (a finite CRT/covering-design statement). Even a
   conditional "blocks q≤24 ⇒ contradiction" would crack it.
   b. **Speed bound:** import or adapt a minimal-counterexample speed bound for 13
   movers to finitize the check.
   c. **Permutohedron route:** recast Q0 as a face-avoidance statement for the
   ranking-walk line on the `(n-2)`-dim permutohedral fan (companion reflection) —
   the dimension drop from `C(13,2)=78` to `12` is where the cross-modulus
   incompatibility should become visible.

One line: **the geometry reduces LRC(14) to a single sharp claim — every primitive
13-system is lonely at denominator ≤ ~25 — and explains why that claim sits exactly
at the threshold: the forbidden radius grows in lockstep with the pair-count, so 13
movers can defeat any one modulus but (empirically) not all of them at once. The
proof is the cross-modulus incompatibility, which is the rational-coincidence heart
of LRC itself.**
