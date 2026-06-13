---
source: claude-2026-06-03-S603
status: framework + proved identities + a Vitali-wall resolution of the finer-invariant question
tags: [LRC, depth-polynomial, independence-polynomial, OCF, Helly, Vitali-wall, correlation, additive-resonance, recursive-log]
---

# The LRC depth polynomial, and why the finer invariant lives only at top order

The prompt gave the map I needed:

> Helly number = "how many orders of overlap you must keep"; Vitali wall =
> "finite moments cannot decide `p_0`"; Collatz/two-block = "correlated residue
> where density is blind"; OCF/partition functions = the sibling world where
> independence polynomials already play the depth-distribution role.

All four land on one object: the **depth generating polynomial**
`P(z) = Σ_k p_k z^k = ∫₀¹ z^{depth(t)} dt`. This is the LRC **independence
polynomial** (the OCF sibling). Everything below is a property of `P`.

## The identity that ties the threads

`depth(t) = Σ_i 1_{F_i}(t)`, so `z^{depth} = ∏_i(1+(z-1)1_{F_i})` and

```
P(z) = Σ_{S⊆[n]} (z-1)^{|S|} m_S = Σ_j M_j (z-1)^j,
      M_j = Σ_{|S|=j} m_S = E[C(depth,j)]  (total j-fold OVERLAP = Helly data),
p_0  = P(0) = Σ_j (-1)^j M_j               (inclusion–exclusion, verified exact).
```

So **tight ⇔ `0` is a root of `P`**, and the `M_j` are literally "the orders of
overlap." The **Helly number** is how far up the `M_j` ladder you must read.

## Finding 1 — tightness is a correlation effect (density is blind)

If the arcs were independent, `depth ~ Binomial(n, 2/(n+1))` and
`p_0^indep = (1-2δ)^n = ((n-1)/(n+1))^n → e^{-2} ≈ 0.135 > 0` for **every** `n`.
**Density alone never produces `p_0 = 0`.** A tight config evacuates the depth-0
stratum purely by **arithmetic correlation**: the high-order overlaps
`M_3, M_4, …` run `1.3×–9×` above the independent prediction. This is exactly the
Collatz/two-block phenomenon — the *correlated residue that density cannot see*.

## Finding 2 — the Helly number is full

The partial alternating sums `S_h = Σ_{j≤h}(-1)^j M_j` reach `0` **only at
`h = n`** for tight configs. You cannot truncate; every order of overlap is load-
bearing. (For non-tight configs `S_h` settles to `p_0 > 0`.)

## Finding 3 — the Vitali wall *is* the answer to "find the finer invariant"

I went looking for a finer invariant separating sporadic tights from the many
non-tight additive chains. The data says there is **provably none at finite
order**:

```
(1,3,4,7) tight : S = +1.000 −0.600 +0.281 −0.057  0.000
(1,2,3,5) NOT   : S = +1.000 −0.600 +0.347 −0.027 +0.053
```

They agree at orders 0–1 and stay close through the middle, **separating only at
the top order**. So no bounded set of overlap orders — equivalently no finite set
of moments — can decide tightness. The finer invariant lives *only* in the
degree-`n` cancellation `P(0)`. This is the **Vitali wall** in the cleanest
possible form: a structural obstruction, not a failure of search. It explains why
the earlier negative results had to happen — Lucas is tight only at `n=4`, Paley
QR only at `p=11` — because **any bounded-complexity family (Lucas, Paley, a
moment threshold) is exactly the kind of finite-order signature the wall forbids.**

## Finding 4 — the arithmetic source: additive resonance

The high-order correlation has a concrete origin. A triple `v_k = v_i + v_j`
**co-resonates** at the origin: where `v_i t, v_j t ≈ 0`, also
`v_k t = (v_i+v_j)t ≈ 0`, so three forbidden arcs pile up. Measured at `δ=1/7`:
`m_{1,3,4} = 0.071` vs independence `0.023`. This is why **additive-chain is
necessary** — you need these high-order resonances to supply the correlation
that cancels `p_0^indep`. It is not sufficient because small speeds overlap
heavily regardless, so a chain can over- or under-shoot the exact cancellation.

## Finding 5 — `P(z)` is an independence polynomial but not real-rooted

Unlike claw-free independence polynomials (Chudnovsky–Seymour), `P(z)` carries
complex roots. Tight ⇔ `min|root| = 0`; the non-zero roots do **not** classify
sporadic-vs-AP. No simple algebraic separator — again exactly what the Vitali
wall predicts.

## What this resolves, and the next object

The "finer invariant" question is answered in the negative, and that negative is
the content: the collapse family is the zero set of the top-degree cancellation
`P(0)=0`, a condition the additive chain makes *possible* (necessary high-order
resonance) but not *certain*. So the productive next object is **not another
moment** but the **resonance graph** (which triples `v_k=v_i+v_j` fire) and the
**witness lattice** (S553's half-division points) — the combinatorics of *which*
resonances align, which is where a genuine classification, if any, must live.

## Open / next

1. **Resonance graph → tightness.** Build the graph on `V` with an edge/triangle
   for each additive relation; is tightness a property of this graph plus the
   speeds' sizes? (A graph invariant would be the real "finer invariant," and it
   is allowed to be unbounded-order.)
2. **Witness lattice.** Do the sporadic tights' measure-zero witnesses sit on a
   common lattice (half-division points `j/(2n)` per S553)? A shared witness
   lattice would be a classification handle the Vitali wall does *not* forbid
   (it is not a finite-moment condition).
3. **OCF transfer.** `P(z)` is an independence polynomial; import a genuine OCF /
   independence-polynomial identity (deletion–contraction on the resonance graph)
   to compute `p_0` structurally rather than by inclusion–exclusion.
