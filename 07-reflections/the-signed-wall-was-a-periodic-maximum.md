# The signed wall was a periodic maximum

**Source:** mac-mini-2026-06-21-S6. Dispatch: push hard on the LRC proof, generate creative angles,
push hypotheses live, integrate the other agents in real time. Canon: THM-563, HYP-2787. Built on
the team's diagnosis (HYP-2784: the absolute bound is 125× lossy, the wall is signed cancellation;
HYP-2785/2786: the one-far error is a Dedekind/Abel endpoint object) and kps's HYP-2788 (the regime
dichotomy reducing multi-far to single-far).

## A wall made of estimates

For many sessions the single-far deviation `Δ_w = p0(B∪{w}) − Φ(B)` resisted. The natural bound —
put absolute values on every Fourier mode — gave `(6/49)V/w`, and the team proved (HYP-2784) it is
intrinsically `~125×` too big at the binding `w`: the signed value `~0.007` against an absolute
budget `~0.9`. So the bound *had* to keep the sign, and "signed cancellation" became the wall. The
next instinct was the right kind of object — Dedekind sums, the Abel endpoint identity — but it came
wrapped in the expectation that bounding it would need *reciprocity*, an *equidistribution tail*, an
*estimate*. The wall was assumed to be analytic.

## It was periodic

Write the deviation exactly:
```
Δ_w · w = Σ_{j=1}^{6} Σ_{endpoints t of A_j} ± S_j(frac(w·t)),
```
`A_j = {x : B misses exactly sector j}`, `S_j` the centered sawtooth. This is a generalized
Dedekind sum — that part the team had. The one thing that dissolves the wall is a triviality once
seen: **the arcs `A_j` depend only on the base `B`, never on the far runner `w`.** The far runner
enters *only* through `frac(w·t)`, and the endpoints `t = k/(7e)` are fixed rationals. So `w·Δ_w`
is a finite sum of sawtooths evaluated at `frac(w·t)` — and that is a **periodic function of the
integer `w`**, with period `7·lcm(B)`.

A periodic function has no tail to estimate. Its supremum is a maximum over one period — a finite,
exact, rational computation. There is nothing to bound: you compute it. For the binding bases,

```
k=8:  period 420,  max(w·Δ_w) = 1        (at w=175)
k=9:  period 2940, max = 43/49 ≈ 0.878   (at w=1659)
k=10: period 5880, max = 1007/980 ≈ 1.03 (at w=231)
```

all below `15·margin_k` — so `Δ_w ≤ max/w ≤ max/15 < margin` for every `w ≥ 15`. The single-far
case closes for all far runners at once, with no asymptotics, no Koksma constant, no reciprocity.
The "conditionally convergent signed sum" was never a sum to be summed; it was a periodic sequence
to be maximized.

## The lesson

When an obstruction is described as "signed cancellation in a conditionally convergent sum," the
reflex is to find the cancellation analytically — a clever majorant, a reciprocity law, a tail
estimate. But conditional convergence is a property of how you chose to *write* the object, not of
the object. Here the object — the exact deviation — is a finite combinatorial sum whose only
`w`-dependence is through fractional parts of `w` times *fixed* rationals. Fixed rationals force
periodicity, and periodicity replaces every estimate with a maximum. The wall was analytic only in
appearance; underneath it was a finite list of rational numbers, and the largest one is `1`.

The same shape has recurred all session: kps's `L7` two-clock discrepancy turned out to be
**residue-only** (a `7×7` table, HYP-2739), not a Koksma constant; the multi-far gap turned out to
**reduce to single-far** by a regime dichotomy (HYP-2788), not to need a genuine joint estimate. Three
times the "hard analytic input" was a finite object hiding behind an analytic description. The far
runner, which looked like it carried the difficulty off to infinity, only ever moved the base's fixed
arcs around a circle — and a point moving around a circle comes back.
