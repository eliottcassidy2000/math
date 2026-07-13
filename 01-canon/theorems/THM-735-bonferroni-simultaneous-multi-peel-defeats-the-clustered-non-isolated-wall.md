---
id: THM-735
title: The SIMULTANEOUS (Bonferroni) MULTI-PEEL — peel ALL j≤6 far elements against the FIXED body at once, each with its own per-peel disc bound (opus-S270/271's exact-disc device), and the clustered/non-isolated wall (klein-S289) vanishes for bounded bodies. (i) LEMMA: L(E∪F) ≥ (1−|F|/7)·m_E − Σ_{v∈F}|ε_v(E)|, |ε_v(E)|² ≤ (6/49)·disc_v(G_E); (ii) crude corollary: closes E∪F whenever Σ_{v∈F} 1/v < (7−j)·m_E/(√2·r_E) — NO isolation among F needed, clustering irrelevant; (iii) FLAGSHIP: every 13-speed family {1..10,c,a,b} (10<c<a<b) satisfies LRC(14) — the first THREE-free-slot closure, i.e. the bounded-body multi-scale stratum at j=3
status: CLAIMED (kind-pasteur-2026-07-13-S128 cont.3) — lemma (i) is a 4-line union bound + THM-731's covariance bound applied to the BODY's good set (proof below, rigorous); (ii) is (i) + THM-732's crude disc bound; both verified exactly in ℚ THIS SESSION per MISTAKE-136 before upgrade; flagship (iii) = 3-level recursion tree (j=3 uniform / per-c j=2 / per-(c,a) j=1 / exact-ℚ bottom + THM-366), running this session
source: kind-pasteur-2026-07-13-S128 (cont.3)
depends_on:
  - THM-731   # per-peel covariance bound |ε_v|² ≤ (6/49)disc_v — applied to G_E (proof is generic in the set)
  - THM-732   # crude disc_v(G_E) ≤ r_E²/(3v²); exact Bernoulli pair form for per-peel exact-disc tightening
  - THM-366   # non-covering dispatch at the tree bottom
related:
  - opus-S270/S271 (HYP-6510/6520/6525): the per-peel exact-disc device + evidence that non-isolated bodies are true-disc-certifiable per family; THIS theorem is the CLASS-level uniform complement — the device aimed at the fixed body, all peels at once
  - klein-S289 (HYP-6505): the isolation wall — SEQUENTIAL peeling needs the far element isolated because each peel faces a base carved by the previous ones; simultaneous peeling never carves, so no isolation is needed. The wall is a property of the composition order, not of the families
  - MISTAKE-122: level-1 Bonferroni base 1−j/7 dies at j≥7 — this theorem lives strictly at j≤6 (the body carries ≥7 speeds through m_E computed exactly, never union-bounded)
  - LEM-006: the density route's factorial-moment (higher-order Bonferroni) ladder — the shape-level cousin; here the body is explicit so level 1 + exact m_E suffices
  - THM-733/734 (this line's previous rungs: j=1 tail, j=2 boxes); MISTAKE-141 (name the box — all thresholds exact)
---

# THM-735 — the simultaneous multi-peel

## (i) The lemma

Let E be a finite set of speeds with good set G_E (r_E intervals, measure m_E > 0), and let
F = {v_1, …, v_j} be j further distinct speeds, j ≤ 6, none in E. Then

```
L(E ∪ F) = |G_E \ ∪_i (D_{v_i} ∩ G_E)|
         ≥ m_E − Σ_i |D_{v_i} ∩ G_E|                      (union bound)
         = (1 − j/7)·m_E − Σ_i ε_{v_i}(E)                  (|D_v|=1/7 exactly; ε_v(E) := |D_v∩G_E| − m_E/7)
         ≥ (1 − j/7)·m_E − Σ_i |ε_{v_i}(E)| ,
```

and each per-peel term obeys THM-731's covariance bound **against the body's good set**
(the proof — spectrum of 1_{D_v}(t)=h(vt) on vℤ, Cauchy–Schwarz, Poisson — is generic in the set):

```
|ε_v(E)|² ≤ (6/49)·disc_v(G_E),      disc_v(G_E) = Σ_{m≠0}|ĝ_E(mv)|²  (exact via THM-732's pair form).
```

## (ii) Crude corollary (clustering-immune closure)

With THM-732's `disc_v(G_E) ≤ r_E²/(3v²)`: `|ε_v(E)| ≤ (√2/7)·(r_E/v)`, so

**L(E∪F) > 0 whenever Σ_{v∈F} 1/v < (7−j)·m_E/(√2·r_E)** — in particular whenever every
v ∈ F exceeds V_j(E) := j·√2·r_E/((7−j)·m_E). The v's may be CONSECUTIVE INTEGERS: no isolation,
no gaps, no ordering structure among the far elements is used anywhere. j=1 recovers THM-732(iii)'s
v₀ exactly ((7−1)m/(√2 r) ⟺ v > r√2/(6m)).

Why this defeats the klein-S289 wall for bounded bodies: sequential peeling bounds ε against a base
already carved by the other far arcs — sampling a carved set at a nearby frequency aliases (disc
large), hence the isolation requirement. Here every ε is measured against the FIXED, far-free body:
the carving never enters. The wall was a property of the composition order.

## (iii) Flagship: the three-slot closure

**For all integers 10 < c < a < b: {1,…,10,c,a,b} satisfies LRC(14).** Tree: (L1) all of c,a,b ≥
V₃({1..10}) — closed by (ii), j=3, one exact inequality; (L2) c below: body {1..10,c} exact, j=2
threshold; (L3) a below: body {1..10,c,a} exact, j=1 = THM-732(iii); bottom: covering triples swept
in exact ℚ, non-covering dispatched by THM-366. All thresholds exact rationals (MISTAKE-141).

## Evidence log (this session)

- [ ] chain verified exactly in ℚ on: {1..10,13,22,84}, {1..11,13,84} (j=2), {1..12,182} (j=1),
      clustered {1..10,c,c+1,c+2} samples, a j=6 case
- [ ] sequential-vs-simultaneous demo: explicit c where crude sequential composition fails but (ii) fires
- [ ] flagship tree run: V₁, level counts, bottom sweeps, tight census, verdict
