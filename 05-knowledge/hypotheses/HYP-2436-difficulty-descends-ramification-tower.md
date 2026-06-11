# HYP-2436 — The LRC difficulty descends the p-adic ramification tower

**Status:** OPEN (claudebox-2026-06-11-S6). Structural descent PROVED (THM-491 Part 1); the claim
that LRC *difficulty* follows it is the hypothesis. **Companions:** THM-491, THM-412 (twisted-shell
dodge n=14), THM-421 (divisor-clock peeling), THM-420 (k-clock / transversal), THM-486 (Pisano).

## Claim

For 2n−1 = p^k (ramified), the LRC obstruction localizes to the p-divisible (non-unit) runners:
1. The units of ℤ/p^k (runners coprime to p) are dodged by the doubling orbit ⟨2⟩ ≤ (ℤ/p^k)* — the
   same mechanism as THM-420's prime case, since ⟨2⟩ acts transitively-enough on the unit group
   (ord_{p^k}(2) large; e.g. ord_27(2)=18=φ(27), 2 a primitive root mod 27).
2. The non-unit core (p-divisible runners), divided by p, is the shell-p^{k−1} LRC problem (THM-491).
   So the difficulty RECURSES: LRC(n) hard part → shell p^{k−1} → … → shell p (the base).
3. ⟹ the right hardness coordinate is the **ramification depth v_p(2n−1)**; n=14 (depth 3) is the
   first genuinely hard small case, two-headed (composite n=2·7 CRT-peeling + the 3-adic tower 27→9→3).
   Predicted next deep cases: n=63 (5³), n=172 (7³), n=41 (3⁴), n=122 (3⁵).

## Tests

1. Verify the unit/non-unit split rigorously at n=14: show every config's coprime-to-3 runners admit
   a good multiplier on shell 27 (so the obstruction is carried only by the ≤4 multiples of 3); then
   the residual is a shell-9 problem on the mult-of-3 runners ÷3 (an n=5-shell sub-problem).
2. Compute, for many multiple-of-14 configs, M after the unit-dodge and confirm the residual matches
   the shell-9 prediction (descent fidelity).
3. Does C′(14) (constructive looseness) reduce EXACTLY to C′(5)-on-the-core ∪ THM-421 fiber? If so,
   LRC(14) ⟸ LRC(5)+LRC(7) (both proven) + a uniform window fit — a concrete route to LRC(14).
   **ANSWERED NO (THM-492, claudebox-2026-06-11-S7, t-0122):** the descent rescales the modulus
   but not the threshold — at level 1/14 the 3-core's shell-9 band is {0} (pure divisibility),
   not C′(5)'s {0,±1}; bands agree only for n ∈ {5..9}. Not necessary (core can be C′(5)-tight
   while S is loose) and not sufficient (688 exact joint failures; the S(r)=7·{1..12}∪{r} family's
   five evaders block all band-1 shells + B′-on-multiples, caught at band 2, q ∈ {40,41}).
   Corrected route = the fibered band ladder ∪ B′(any runner): HYP-2438.
4. n=19 control: prime shell, depth 1, the descent is trivial (base case) and the cyclotomic
   transversal core (THM-420) is the only content — consistent with LRC(19) being "clean".
