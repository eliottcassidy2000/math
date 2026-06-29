# The LRC cap is Claim A without the parity

*mac-mini-2026-06-29-S7. A niche reframe found while looking for small creative connections: the lonely measure is the project's own Claim A, run on the even (x=−1) side instead of the odd (x=2) side.*

## The shared skeleton

Claim A — the project's central deletion–contraction — says `H(T) − H(T−v) = 2 Σ_{C odd cycle ∋ v} μ(C)`: adding a vertex bumps the Hamiltonian-path count by twice a sum over the cycles through it. I went looking for the Lonely-Runner analog and it is exact (verified n=6,8,10,14):

> `cap(S\{s}) − cap(S) = meas(lonely(S\{s}) ∩ D_s) = 2 · Σ_teeth μ_LRC(s,k)`.

Adding a speed bumps the lonely measure *down* by its **conditional danger** — the part of `s`'s danger comb where everything else is already safe — and that conditional danger is *twice* a sum over the `s` **teeth** of `D_s`. The teeth are the cyclic positions `k/s` where speed `s` wraps the circle: they are the "odd cycles through `s`." The factor 2 is the same complement `Z₂` that runs through the whole project — `lonely ∩ D_s` is invariant under `t↦−t`, so it is twice its half-circle part. Peel speed by speed from `cap(∅)=1` and you get the measure-valued OCF: `cap(S) = 1 − Σ_s (conditional danger of s)`, and the floor is exactly `Σ(conditional dangers) < 1` — a conditional union bound, strictly tighter than the hopeless naive `|S|/n > 1`.

## Where the twins part ways

So the LRC cap and Rédei's `H` are built on the *same* deletion–contraction skeleton, with the *same* factor 2 and the *same* "sum over cyclic through-pieces." Yet one is famous and elementary (Rédei) and the other is a Fields-medal-grade open problem. The reframe makes the reason precise, and it is exactly the two-index split (THM-582):

- **Tournament side (the x=2 evaluation, odd).** `H` of the empty start is `1` — *odd* — and every deletion–contraction step adds `2 Σ μ(C)`, which is *even*. Oddness is preserved by induction. Rédei is a **parity** theorem; the measure of the increment never matters, only that it is even.
- **LRC side (the x=−1 evaluation, even).** `cap(∅) = 1`, and every step *subtracts* `2 Σ μ_LRC` — but `cap` is a real number in `[0,1]`, not a parity. There is no "stays odd"; there is only "stays positive," and that is a genuine **measure inequality**.

The LRC cap is *Claim A without the parity*. The skeleton transfers perfectly; the proof does not, because the tournament invariant is integer-valued (so oddness can do the work) while the runner invariant is measure-valued (so you actually have to bound something). This is the cleanest small statement of why the two problems — which the project keeps sensing are "the same" — are genuinely the odd and even halves of one structure: same deletion–contraction, different category of conclusion.

## A sibling route to the descent

The chain rule hands over an alternative target for the covering floor: order the peel **odd speeds first, then even** — the same ordering the 2-adic descent (THM-580) uses — and bound the per-tooth conditional dangers `μ_LRC(s,k)`. Each is a circular-arc overlap measure (lonely-so-far ∩ one tooth), a clean Diophantine object. If `Σ < 1` falls out of bounding those, the floor closes through the Claim-A skeleton rather than through the SOS — a sibling of the descent, in the same even category, and possibly easier to make uniform because each term is a single-arc measure.

The meta-lesson, again: the project's invariants are evaluations of one polynomial at two points, `x=2` (odd, parity, done) and `x=−1` (even, measure, open). Every "is the LRC like Rédei?" question is answered by *yes, at the other evaluation point* — and the work that remains is always the measure bound the parity argument got to skip.

See [[the-two-indices-redei-is-odd-lonely-is-even-half-tiling-is-the-quotient]] (THM-582), [[the-topology-of-the-lrc-cap-is-euler-char-of-the-cover-nerve-lonely-is-the-hole-certified-by-borsuk-ulam]] (HYP-3242). Claim A: THM-070. The LRC twin: HYP-3537.
