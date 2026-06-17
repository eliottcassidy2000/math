# Why every LRC-14 positivity method fails: a signed integral with no floor

**Session:** mac-mini-2026-06-16-S2
**Complementary to:** kind-pasteur THM-522 / HYP-2561 / T833 (the *constructive* route — scale-invariance + quantization + compactness, inf L = 1/1260) and MISTAKE-075. This reflection is about the *negative* space: why four independent positivity attacks all died, and what their shared cause says.

## The collision, honestly

I spent this session attacking inf_S L(S) > 0 for LRC-14 while kind-pasteur-S7 concurrently solved the headline. We both re-walked the same restricted family — the "interior-drop cores {1..13}\{j} ∪ {14m}" with a *multiple-of-14* stranger — and I pinned its exact infimum, **7/858 − 1/7³ = meas(Lonely({1..13}\{6})) − 1/343 = 1543/294294**, with the resonant stranger 98 = 2·7² removing *exactly* 1/7³ (a clean 7-adic fact). But that family is the artifact MISTAKE-075 names: the true extremizers are *sporadic-tight perturbations* of the maximally-tight AP {1..13} (which has L = 0 exactly), and the real infimum is **1/1260 ≈ 0.000794**, six times smaller, at {1,…,11,13,36}. I independently verified 1/1260, L({1..13})=0, and scale-invariance, and a broader search (single perturbations to w≤600, doubles to w≤140) found *nothing* with 0 < L < 1/1260 — corroborating HYP-2561. So the exact-value frontier is kind-pasteur's; my contribution is the corroboration and the map of what does *not* work.

## The four closed doors

An eight-angle assault (LLL, 7-adic split, Selberg–Beurling, Abel summation, OCF-bridge, large-speed decoupling, exact-rational, finite-reduction) produced four clean dead-ends:

1. **Lovász Local Lemma / Shearer.** The 13 danger events d_i = {‖v_iτ‖ ≤ 1/14} are *positively* correlated — pairwise overlaps reach 1/14, which is 3.5× the independent value. Shearer's independent-set polynomial is *negative* (Z(G;−1/7) = −36/343 < 0); p = 1/7 sits past the Shearer radius ≈ 0.117. LLL wants negative dependence; the LRC danger events have the opposite sign.

2. **Selberg–Beurling nonnegative minorant.** No nonzero nonnegative trig polynomial minorizes the indicator of a closed band; the only valid finite-degree pointwise minorant is the zero function. This is the exact *dual* of Bedert's Riesz-product wall — both extremal-function routes die on the indicator's sharp edge.

3. **Abel / Cesàro / Euler resummation of the level masses Λ_k.** The signed Λ_k do not decay geometrically (ratios 1.4–3.7) and the unsigned per-level mass diverges for k ≥ 3. No elementary summation method brackets L from below; a per-level Pólya–Vinogradov bound Gᵏ is simply false.

4. **The OCF / tournament-partition-function bridge.** I had hoped L(S) = Σ_{t∈Λ} ∏h(t_i) might *be* an independence polynomial I(Ω, 2) = H, transporting the project's home-turf machinery. It is not. The decisive difference: **H = Σ 2^{|U|} is termwise positive — it has a structural floor of 1 — whereas L is a signed sum whose level masses grow.** The "two sevens" (LRC's band 14 = 2·7, the OCF's forbidden phantom Φ₃(2) = 7) are coincidental: one is half a runner count, the other a cyclotomic value at the fugacity.

## The one obstruction underneath all four

Every one of these failures has the same root, and naming it is the point of this reflection:

> **L(S) is an archimedean *signed* singular integral with no termwise floor, evaluated on *positively correlated* events.**

- *Signed, no floor* kills the minorant routes (Selberg, OCF-bridge): there is no positive quantity to bound below termwise, because positivity of L is *manufactured entirely by cross-level (−1)^{|T|} cancellation*. This is THM-503(4)'s "archimedean singular integral, not an Euler product" seen from four new directions at once.
- *Positive correlation* kills the probabilistic routes (LLL/Shearer): the danger events cluster rather than repel, so locality lemmas built for near-independence point the wrong way.
- *Growing level masses* kill the resummation route: the series converges only conditionally, by alternation, not by decay.

This is why kind-pasteur's route is the one that lives. Scale-invariance and quantization don't try to *bound a signed sum* at all — they sidestep positivity entirely, reducing inf L > 0 to a *compactness* statement (the L-minimizers have bounded lcm) plus the rational quantization floor L ≥ 1/(14·lcm S). The lesson is almost a slogan: **when positivity is an emergent cancellation rather than a termwise fact, stop trying to certify it termwise — quantize the values and compactify the domain instead.**

## The meta-pattern worth keeping

Twice this week the same human dispatch landed on two machines at once (the Fibonacci/Heron theme → kind-pasteur's boat T830 vs my octal lens; LRC-14 → kind-pasteur's THM-522 vs this). Both times the collision was *productive*: the second machine, arriving at the already-mapped headline, was forced into the negative space — the dead-ends, the restricted-family artifacts, the "why the obvious thing fails." That negative space is not waste. A proof program is only as trustworthy as its list of *ruled-out* alternatives, and four adversarially-tested closed doors are exactly what tells the next session that the quantization+compactness route is not one option among many but the only one left standing. The infimum 1/1260 is kind-pasteur's; the knowledge that nothing termwise can ever reach it is the thing I can hand forward.
