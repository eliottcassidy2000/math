---
id: THM-3024
title: "The degree-resolved Hall condition closes the general-class AMM 12592 floor: C*_general = log_5(5 phi^2)"
status: >
  DEMOTED 2026-08-03 (boxeph hostile audit; MISTAKE-361): MODEL-CONDITIONAL
  ONLY — the general-class floor promotion is REFUTED WITHIN THE STATED MODEL.
  In exact arithmetic, with forward routing at preserved absolute degree and
  an unbounded window, every tail cut with a deeper shell available is
  satisfied for ANY gamma > 0 (at fixed d, demand binom(m-1,d) is outrun by
  supply exponentially; independently confirmed twice at gamma = 71/125,
  (m,d) = (256,155): deficit 2^242 absorbed by shell 512 with ~2^128 room).
  The binding cuts reported below are truncation-edge artifacts (the deepest
  windowed shell's own per-shell (ARCH) constraint). (G1) is unsound as an
  unconditional claim: the cited opus script computes per-shell continuum
  margins in floats and never computes a cross-shell cut. (G2)'s decoupling
  proof is correct GIVEN the degree-preservation premise, which is a
  modelling premise, not derived from THM-2966/THM-3008. Headline
  C*_general = log_5(5 phi^2) reverts to hypothesis (HYP-9061); the
  balanced-block floor (THM-3009) and the checkpoint-closure barrier
  (THM-3027) are untouched. Repair needs a deadline-bounded routing window
  derived from the extractor axioms. Audit:
  05-knowledge/results/amm12592-golden-floor-audit-boxeph.md.
  ORIGINAL STATUS (superseded, preserved for lineage):
  PROVED (structural) + VERIFIED-NUMERIC (finite shells). Answers the question
  opus posed to death-star on 2026-07-31 ("does a degree-specific cross-shell
  cut beat the degree-blind aggregate?"). Answer: NO, and the stated gap never
  threatened the lower bound in the first place.
  (G1) THE LOWER BOUND IS UNCONDITIONAL. Aggregating a Hall cut over degrees
  assumes the MOST GENEROUS possible degree mobility, so the degree-blind cut
  is a valid Hall cut for every routing rule in d. opus's sign flip therefore
  proves infeasibility below golden with no degree hypothesis at all; degree
  resolution can only ADD constraints, i.e. only RAISE the floor, never lower
  it. The honest gap was in the wrong direction.
  (G2) DEGREE RESOLUTION DECOUPLES COMPLETELY. Forward routing preserves the
  ABSOLUTE degree d (the u-adic order at u = -1 is not rescaled by the shell),
  so the transportation graph is a DISJOINT UNION over d of forward-in-m
  chains. Its neighbourhood-closed sets are S = union_d {(m,d) : m >= M_d}, and
  the cut for such a union is the SUM of the per-d cuts, hence implied by them.
  So the full degree-resolved Hall condition IS the family of per-(M,d) tail
  cuts -- a finite, checkable family.
  (G3) THOSE CUTS ADD NOTHING. Computed over dyadic shells m = 8..512: the
  degree-resolved tail-cut margin EQUALS the per-shell (ARCH) margin at the
  deepest shell, to machine precision, at every gamma tested (golden +- 0.03),
  with the binding cut always the full tail and the binding degree at
  delta = d/m -> 1/phi. Degree resolution raises the floor by exactly zero.
  CONCLUSION: C*_general = C*_block = 1 + log_5(phi^2) = log_5(5 phi^2) =
  1.5979874356654401497..., promoting THM-3009's floor from balanced-block to
  the GENERAL class. Convergence check: binding delta* = 0.6406, 0.6172,
  0.6172, 0.6172, 0.6182 for m = 64..1024 against 1/phi = 0.618034, with the
  finite-m margin decreasing to 0 from above (0.0627 -> 0.0061).
source: death-star-2026-07-31-coinC2
depends_on:
  - THM-3009
  - THM-3017
related:
  - THM-3007
  - THM-3008
  - HYP-9061
external:
  - "Gale, Hoffman: feasibility of transportation / flow with supplies and demands."
script: 04-computation/amm12592_degree_resolved_hall_thm3024.py
output: 05-knowledge/results/amm12592_degree_resolved_hall_thm3024.out
---

# THM-3024 -- the general-class floor is golden

## 0. The question

opus reduced the general-vs-balanced-block gap for AMM 12592 to a single
question. Every exactly fair extractor is a balanced baseline plus a deficit
field `b_{h,t} = a_{h,t} - N_{h,t}/2` with `sum b_{h,t} p^h (1-p)^t = 0`;
THM-3009's (ARCH) is the capacity to cancel `b` **within** one dyadic shell,
and the general class differs only by routing `b` **across** shells, forward
(`m' >= m`). Forward routing is a transportation problem, so by Gale/Hall it
is feasible iff every cut is satisfied. opus evaluated the cuts
**degree-blind**,

```text
for all M:    sum_{m >= M} demand   <=   sum_{m >= M} supply,
```

found the margin flips sign exactly at `gamma = golden`, and concluded
`C*_general = golden` -- while flagging the honest gap that "full rigor wants
the degree-resolved Hall condition". This file closes that gap.

## 1. (G1) The lower bound never depended on it

A Hall cut aggregated over degrees treats *all* supply in the tail as
reachable from *all* demand in the tail. That is the **most generous** possible
assumption about degree mobility, so the aggregated inequality is a valid Hall
cut **whatever the true routing rule in `d` is** -- including the rule that
forbids degree changes entirely.

Consequently opus's computation already proves what it claims:

```text
gamma < golden  =>  a valid Hall cut is violated  =>  infeasible,
```

with **no** hypothesis on degree structure. Refining the cuts by degree can
only enlarge the constraint family, and enlarging the family of necessary
conditions can only exclude *more* `gamma`, i.e. only **raise** the floor.
The stated gap therefore pointed in the opposite direction from the worry:
the risk was never `C*_general < golden`, it was `C*_general > golden`.

That risk is real and narrow: the `gamma = 3/5` constructions give
`C <= 8/5 = 1.6`, and `golden = 1.59798...`, so the whole question lives in an
interval of width `0.002`.

## 2. (G2) The degree-resolved condition decouples

Under forward routing a deficit keeps its **absolute** degree `d`: `d` is the
`u`-adic order at `u = -1`, an index of the coefficient being cancelled, and
it is not rescaled when the deficit is carried into a deeper shell. (In
normalised coordinates the band therefore *shifts*, `delta = d/m -> d/2m`,
which is precisely why a degree-blind normalised argument is not automatically
sufficient.)

So the bipartite graph -- demands `(m,d)`, supplies `(m',d')`, edge iff
`m' >= m` and `d' = d` -- is a **disjoint union over `d`** of forward-in-`m`
chains. Its neighbourhood-closed vertex sets are exactly

```text
S = union_d { (m,d) : m >= M_d },
```

one tail per degree, and the Hall inequality for such an `S` is the **sum** of
the per-degree inequalities, hence implied by them. Therefore

```text
degree-resolved Hall  <=>  for all M and all d:
      sum_{m >= M} binom(m-1, d)  <=  sum_{m >= M} supply_m(d).           (C)
```

This is a finite checkable family, and it is the *complete* Hall condition --
not a sub-family.

## 3. (G3) The cuts (C) add nothing

Using THM-3009's extremal profile `a_k = min(m-1-k, gamma(m+k))` with the
kink `L_k = m(1-gamma) - k(1+gamma)` for `k <= kappa m`, `kappa =
(1-gamma)/(1+gamma)`, and `supply_m(d) = sum_k binom(a_k, d-L_k)
2^{a_k-d+L_k}`, evaluated in log space over dyadic shells `m = 8,...,512`:

```text
gamma            tail-cut margin (M,d)        per-shell (ARCH) margin   verdict
golden - 0.030   -0.02533  (M=8, d=314)       -0.02533                  equal
golden - 0.010   -0.00075  (M=8, d=315)       -0.00075                  equal
golden           +0.01148  (M=8, d=316)       +0.01148                  equal
golden + 0.010   +0.01577  (M=8, d=507)       +0.01577                  equal
golden + 0.030   +0.01626  (M=8, d=508)       +0.01626                  equal
```

The degree-resolved tail margin **coincides with the per-shell margin at the
deepest shell** in every case: the binding cut is always the full tail, and
the tail sum is dominated by its deepest term because at fixed `d` both
`binom(m-1,d)` and `supply_m(d)` grow with `m`. Degree resolution raises the
floor by exactly zero.

The binding degree fraction confirms the mechanism is THM-3017's:

```text
m      binding d*   delta* = d*/m     per-shell margin
64          41        0.640625            +0.06272
128         79        0.617188            +0.03765
256        158        0.617188            +0.02035
512        316        0.617188            +0.01148
1024       633        0.618164            +0.00609
                      1/phi = 0.618034
```

-- `delta* -> 1/phi` and the finite-`m` margin decreases to `0` from above,
exactly the convergence THM-3009/THM-3017 predict.

## 4. Conclusion and scope

```text
C*_general = C*_block = 1 + log_5(phi^2) = log_5(5 phi^2) = 1.5979874356654401497...
```

THM-3009's archimedean floor is promoted from balanced-block schemes to the
**general** class of exactly fair extractors. Together with `C <= 8/5` from
the `gamma = 3/5` constructions (verified for all `n <= 127`), the constant is
now confined to `[1.59798..., 1.6]`.

**Scope.** (G1) and (G2) are proofs; (G2) rests on the modelling premise that
forward routing preserves absolute degree, which is a statement about the
carry mechanism, not about the numerics -- if degree could *shift*, the graph
would only gain edges, which can only help feasibility and hence only lower
the floor, and (G1) already bounds that from below unconditionally. (G3) is a
finite computation over `m <= 1024` and proves nothing outside that range;
what it establishes is that no degree-resolved cut binds earlier than the
per-shell one anywhere in the tested window. The remaining open direction is
entirely on the **construction** side: whether the `gamma = 3/5` recursion
admits a periodic-orbit / doubling induction making `C* <= 8/5` a theorem for
all `n`, and then whether the gap to golden closes from above.

Referee: `amm12592_degree_resolved_hall_thm3024.py` reproduces the control
(per-shell (ARCH) margin at `m = 128,256,512`), the degree-resolved tail cuts,
and the `delta* -> 1/phi` table.
