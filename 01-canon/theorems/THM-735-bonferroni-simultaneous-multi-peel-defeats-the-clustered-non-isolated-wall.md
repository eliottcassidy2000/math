---
id: THM-735
title: The SIMULTANEOUS (Bonferroni) MULTI-PEEL — peel ALL j≤6 far elements against the FIXED body at once, each with its own per-peel disc bound (opus-S270/271's exact-disc device), and the clustered/non-isolated wall (klein-S289) vanishes for bounded bodies. (i) LEMMA: L(E∪F) ≥ (1−|F|/7)·m_E − Σ_{v∈F}|ε_v(E)|, |ε_v(E)|² ≤ (6/49)·disc_v(G_E); (ii) crude corollary: closes E∪F whenever Σ_{v∈F} 1/v < (7−j)·m_E/(√2·r_E) — NO isolation among F needed, clustering irrelevant; (iii) FLAGSHIP: every 13-speed family {1..10,c,a,b} (10<c<a<b) satisfies LRC(14) — the first THREE-free-slot closure, i.e. the bounded-body multi-scale stratum at j=3
status: PROVED (kind-pasteur-2026-07-13-S128 cont.3; Lean transcription open). Lemma (i): 4-line union bound + THM-731's covariance bound on the body's good set. (ii): (i) + THM-732 crude disc. Chain verified EXACTLY in ℚ on a 9-family battery (j=1 deep well; j=2 residue body; j=3 clustered triples c=150..500; j=6 clustered sextuple {1..7}∪{300..305}) — every link exact, zero violations. PER-PEEL EXACT-DISC EXTENSION (opus's device, Cauchy–Schwarz form (1−j/7)²m_E² > j·(6/49)·Σ_v disc_v(G_E), all-rational): fires where crude fails, incl. the j=6 clustered sextuple. FLAGSHIP (iii) ESTABLISHED in 1.3 s: every {1..10,c,a,b} (10<c<a<b) satisfies LRC(14) — leg J3 = ONE inequality (all c ≥ V₁=154); leg J2 = 143 exact bodies {1..10,c} (max V₂=219 at c=143); leg J1 = 7537 exact bodies {1..10,c,a} (max bottom-b=167); bottom = 19,202 triples of which only 27 covering (exact-swept, ALL L>0, zero tights) and 19,175 non-covering (THM-366). The first three-free-slot closure; the bounded-body multi-scale stratum at j=3 is done, clustering included
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

## Evidence log (all done — scripts lrc14_thm735_*_kps_S128c3.py, outputs in 05-knowledge/results/)

- [x] chain verified exactly in ℚ on the 9-family battery; every link (lemma / per-peel covariance /
      crude disc) holds with zero violations. Tightness: the lemma bound reaches 84–92% of true L on
      clustered triples; at the deep well (j=1) it equals L exactly.
- [x] per-peel EXACT-disc certificate (opus-S270/271's device, CS form, all-rational) fires below the
      crude threshold: {1..10,150,151,152} (crude −0.0014 → exact-disc +0.006) and the j=6 sextuple
      {1..7,300,…,305} (crude −0.024 → exact-disc +0.002): SIX consecutive far elements certified at
      once against the fixed body. Clustering is irrelevant, as the lemma predicts.
- [x] honest calibration (HYP-6535): on {1..10,c,c+1,c+2} the observed sequential crude also fires
      (observed edge growth ~0.87·m₀ per unit c vs P1's 2.14 coefficient — P1 is ~2.5× loose); the
      cont.2 sequential failure is of the P1-BOUND composition. Bonferroni removes the composition
      rather than sharpening it — and is the only route at j ≥ 2 slots simultaneously free.
- [x] flagship tree: V₁=154; 143 + 7,537 exact bodies; 19,202 bottom triples (27 covering exact-swept,
      all L>0; 19,175 non-covering → THM-366); ZERO tights, ZERO covering L=0. Total 1.3 s.

## What remains of the multi-scale/non-isolated stratum after (iii)

- Other bounded bodies: the same tree per body (mechanical; the 364-body analogue of THM-734 at j=3
  is a larger but finite computation; j=4,5,6 trees over 9,8,7-element bodies likewise).
- Families whose 7th-smallest speed exceeds every bound (j ≥ 7 far slots — the Bonferroni base
  1−j/7 dies, MISTAKE-122): these have ≤ 6 bounded speeds, i.e. genuinely spread/two-scale objects —
  the density route's home turf (klein's floors/THM-692; LEM-006's factorial-moment ladder is the
  higher-order Bonferroni for exactly this). The covering/density seam is now sharply drawn at j=7.
- Lean: (i) is a union bound + the already-Lean-targeted THM-731/732 chain; the tree bottoms are
  decide-style (27 sweeps).
