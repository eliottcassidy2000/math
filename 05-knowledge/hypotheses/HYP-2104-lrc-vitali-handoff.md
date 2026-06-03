---
id: HYP-2104
status: CONNECTION (unifies S551o Vitali wall with THM-398) + Criterion B' = Vitali-covering iff (PROVED), quantified 72-97%; alignment residual OPEN
source: opus-2026-06-03-S573
related:
  - THM-398
  - HYP-2103
  - HYP-2102
  - THM-369
---

# HYP-2104: the Vitali handoff equation is n|v

**UNIFICATION:** S551o's "Vitali wall" (LRC = positive-measure bulk by MEASURE + measure-zero core by CONSTRUCTION; Vitali marks the handoff) IS THM-398's split. The handoff S551o left abstract is the equation **n|v**: n-does-not-divide-any-v => construction side (t=1/n, measure-blind, reaches the measure-0 worry-core which has no multiple of n by S564); n|v => measure side (C': positive measure).
**TWO VITALIS:** the Vitali SET (S551o) = the boundary (proves measure can't see the n-point core); the Vitali COVERING lemma / Lebesgue density = the TOOL on the measure side. It APPLIES to n|v because D_v (v=nw) is a periodic bounded-eccentricity arc family (centres k/(nw), radius 1/(n^2 w), gaps (n-2)/(n^2 w)>0) = a genuine Vitali cover; the worry-core is n isolated points {k/n} = no cover structure = measure-blind. So the multiple of n is exactly what makes the question Vitali-answerable.
**CRITERION B' = THE VITALI-COVERING IFF (PROVED):** G(S\{v}) is a finite union of intervals; D_v's gaps are open nonempty => an E-interval lies in D_v iff it fits in one arc (len<=2/(nv)). So S loose <== some component of G(S\{v}) longer than 2/(nv); S tight(measure-0) ==> every component short AND arc-aligned. Residual = the arc-ALIGNMENT case.
**QUANTIFIED (lrc_vitali_covering_residual_s573.py):** B' proves looseness for 72.4%(n=6),78.7%(8),88.9%(10),91.5%(12),**96.8%(14)** of mult-of-n configs -- proved fraction GROWS toward the frontier. All-short residual is NEVER tight (0 across n=6..14): periodic arcs fail to align to cover the short intervals.
**RESIDUAL (open, <=3% at n=14):** Diophantine alignment -- the interval centres of G(S\{v}) (rationals with S\{v}-denominators) would all have to be 1/(n^2 w)-close to the 1/(nw)-lattice; ripe for Erdos-Turan/three-distance.

**See:** `07-reflections/lrc-the-vitali-handoff-equation-is-multiple-of-n-s573.md`, `01-canon/theorems/THM-398-...md`, `04-computation/lrc_vitali_covering_residual_s573.py` (+.out); S551o (Vitali wall), HYP-2103/2102, THM-369, S564.
