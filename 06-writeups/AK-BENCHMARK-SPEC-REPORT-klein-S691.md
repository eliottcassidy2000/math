# Bug + divergence report for the arithmetic-Kakeya certificate benchmark
(prepared klein-S691, 2026-07-28 — self-contained; suitable to forward upstream)

Concerns the problem "find X ⊂ Z², an X-constructible graph G and forcing
pair (R,T) with score ≤ 1.675" (constructible graphs / forcing pairs /
score = (m+|R|)/(n−|T|)).

## Finding 1 — the verifiable set-up's rule (1), as written, is unsound

Rule (1) reads: if f_i(a_1,…,a_i) = x and the first i coordinates of e_1
are (a_1,…,a_i) and the first i coordinates of e_2 are (a_1,…,a_i+1),
then the function [e_1 ↦ x, e_2 ↦ −x] may be added to R. **The suffixes
of e_1, e_2 are unconstrained.** Differencing two such functions sharing
e_2 yields within-layer transporters x·(δ_{(p,s)} − δ_{(p,s')}), which the
glued-copies picture never produces. This admits certificates of score
→ 1, e.g. (machine-verified against the literal rules; forcing succeeds):

```
score=14/9 m=10 r=4 n=9 t=0
[(0, 0), (0, 1), (1, 0), (1, 1)]
[3, 3]
[{(1,): (1, 0), (2,): (0, 1)}, {(1, 1): (1, 1), (1, 2): (1, 1), (3, 1): (1, 1), (3, 2): (1, 1)}]
[]
[{(1, 1): (1, 0)}, {(1, 1): (0, 1)}, {(2, 1): (0, 1)}, {(3, 1): (1, 0)}]
```

and the family dims [D, n₀] (alternating axis-1 labels, z-paths in the two
end copies, one seed per copy) has score 1 + 1/D + 1/n₀ − 1/(Dn₀) → 1.
Since score s is claimed to prove AK(s), the literal rules "prove" the
full arithmetic Kakeya conjecture; hence the intended rule must require
e_1, e_2 to AGREE ON ALL COORDINATES except the i-th (equal suffixes).
If the deployed verifier implements the literal text, the benchmark is
winnable by the certificate above; if it implements equal suffixes, the
problem text should say so.

## Finding 2 — the two set-ups are not equivalent as written (both ways)

- The INTUITIVE definition glues k+1 copies of ONE graph H, with per-step
  AND per-vertex (i.e. per-suffix) choices: identify OR edge with a
  per-vertex label. Its label data is (step, suffix)-indexed; inner
  structure must repeat across copies (same H); identifications exist.
- The VERIFIABLE format's f_i depends on the full PREFIX (inner labels may
  vary across outer positions — not obviously same-H-realizable), labels
  are suffix-UNIFORM within a step, and there is no identification device.
  ((0,0) ∈ X is needed only so the R-functions are total; it does not
  encode identification, since n(G) = d_1···d_k is the raw product.)

So the label algebras are exactly complementary — (prefix+step) vs
(step+suffix) — and the equivalence claimed in the problem statement
needs nontrivial compilation lemmas in both directions. Empirically (exact
searches, this session): the equal-suffix verifiable game's best found
score is 13/7 ≈ 1.857; allowing per-suffix labels (intuitive-style) yields
7/4 exactly (three witnesses, two shapes, n ≤ 9); allowing identifications
yields 12/7 ≈ 1.714 (witnesses available on request). If those games'
infima genuinely differ, the two published set-ups define different
benchmarks; at minimum, small certificates do not transfer.

## Suggested fixes

1. State rule (1) with equal suffixes explicitly (or fix the verifier).
2. Either add an identification device to the verifiable format or
   publish the compilation argument that removes identifications
   score-preservingly.
3. Clarify whether f_i's prefix-dependence is intended (it exceeds
   same-H gluing) — if yes, the intuitive definition's step-2 should say
   copies may differ; if no, restrict f_i to (e_i, ·).

Verification engine and all witnesses: `04-computation/ak_forcing_engine.py`
and companions, this repository, session klein-S691.
