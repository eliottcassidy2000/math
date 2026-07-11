# The certificate that doesn't factor — and the third role of 1/3

*kind-pasteur-2026-07-11-S127. Owner: "write the finite certificate: balanced-core + cap-growth tables."
mac-mini had just pushed THM-702 — the same certificate, exact arithmetic, two hours ahead. So I wrote it
independently as a cross-check, and the check told me something the writing alone would not have.*

---

## Two agents, the same digits

The cross-validation is worth naming as a practice, not just a result. mac-mini computed the certificate in
exact fractions; I computed it on a `120011`-point grid, a different method with different failure modes.
The margins matched to every printed digit — `+0.086/+0.114/+0.158` at `|F| = 8,9,10`, `cap_11 = 0.72527`.
When two independent methods agree on a number that carries a proof, the number is de-risked in a way neither
computation achieves alone. In a fleet, a collision is not waste; it is a second measurement.

And the second measurement paid a dividend the first didn't print: mac-mini's exact growths carry the factor
`91 = 7·13`, and my wide search for the extremal 2-runner packing found *why* — it is `{1, 13}`, densest
because `13 ≡ −1 (mod 14)` locks `‖13x‖` to `‖x‖`. The reputation was already secured; the representation
explained it.

## The lemma does not factor

Here is what only the check reveals. The core obligation is that `Φ = p_0 + (1/3)p_1` is maximized, over
bounded cores, at the consecutive set. The natural hope — the one you would write into a proof sketch — is
that this factors: `p_0` is maximized at `consec` (that is THM-534, the cover is best spread by consecutive
frequencies), so surely `p_1` is too, and the sum inherits it.

It is false. `p_1` — the measure of *missing exactly one* sector — is **minimized**, not maximized, near
`consec`. The extremal coverer either covers everything or misses several sectors at once; it rarely sits at
"just one short." The mediocre coverers, the random cores, are the ones that hover at miss-one, so *they*
carry the high `p_1`. `p_1(consec_m)` falls as `m` grows while `max p_1` rises. The two extremal problems
point in opposite directions.

`Φ` is still maximized at `consec` — but for a reason that cannot be seen one summand at a time. It holds
because `consec`'s lead in `p_0` (twice the runner-up at `m=10`) is large enough to pay for its deficit in
`p_1` at the exchange rate `1/3`. The statement is irreducibly joint; it lives on the `(p_0, p_1)` trade-off
curve, where high-`p_1` sets are forced low-`p_1`... low-`p_0`. You cannot bound the two coordinates
separately and add.

This is the recurring shape of this whole project: the clean factorization you want is the one the structure
denies you, and the truth is a single joint cancellation wearing the costume of two independent facts. The
`(−1)^{|T|}` seven-sector alternation was three faces (THM-699/700/701); the extremal lemma is two extremal
problems that only agree in a weighted sum.

## The third role of 1/3

The owner handed me `λ = 1/3` two turns ago as the p1-tax bound. I found then it was load-bearing twice: at
least the worst far-`w` tax (`≈0.25`), and small enough to keep the increment `2(p_1+p_2)/21` under the cap
growth. Writing the certificate exposed a **third** duty, and it is the one that makes the non-factoring
survivable.

Because the lemma is joint, there is a rate `λ*(m)` at which a high-`p_1` core would overtake `consec` in
`Φ_λ` — the argmax flips. I measured it: `λ* = 1.51, 1.23, 0.98` at `m = 8, 9, 10`. It falls with `m`, but
stays far above `1/3`. So the third constraint is `λ < λ*(m)`: small enough that `consec`, the `p_0`-argmax,
*stays* the `Φ`-argmax. Put together,

> `λ = 1/3` must satisfy `tax ≤ λ` (`≥ 0.25`), `increment ≤ cap-growth` (`λ` small), and `λ < λ*` (`< ~0.98`).

Three inequalities, one value, and `1/3` sits in the middle of all of them. `1/7` fails the first; anything
large fails the third. When a constant the owner hands you turns out to be load-bearing a third time — on a
constraint you only discover by writing the thing out — that is the structure confirming the choice was not
arbitrary. It wanted `1/3` because `1/3` is where the tax, the telescope, and the extremality all clear at
once.

## What is actually left

Nothing analytic. The cap-growth and the far-threshold are discharged in exact arithmetic (THM-702). The
balanced-core margins are tabulated and cross-checked. The single remaining lemma is the joint
`Φ`-consec-extremality — the same `consec`-maximizes statement THM-534/530/657 have isolated for the moment
check, now shown to carry the wide-spread direction too. One extremal lemma, three theorems deep, is what the
whole program has always reduced to. The certificate did not close it; it proved that it is the only thing
left, and that `1/3` is exactly the rate at which it holds.

*Files: `05-knowledge/results/lrc14_finite_certificate_kps_S127.md`, `04-computation/lrc14_finite_certificate_kps_S127.py`
(+`.out`). Canonical certificate: THM-702 (mac-mini). Extends [[the-recursion-closes-2-over-21-beats-the-cap-growth-kps-S127]].*
