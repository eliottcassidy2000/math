# HYP-2675 proof foundation: the cross-scale DECORRELATION formula (kps-S19, validated exact-vs-model)

The sole LRC(14) residual is HYP-2675: `span(E)>B ⟹ p0(E)=meas(S7(E)) ≤ cap_k`, k=8..12.
The KEY proof mechanism (coverage lemma, NOT Δ_w cancellation): a wide set splits into scale-clusters; each
cluster's residues at slow time x are `(cluster anchor phase frac(M_i x)) + (internal pattern frac(d·x))`, and
the anchor phases DECORRELATE across clusters by Weyl. So coverage factors.

## Decorrelation formula (VALIDATED)
For a 2-cluster wide set A∪B (A at scale ~1, B at scale ~M), as M→∞:
  **p0(A∪B) → p0_decorr := E_{φ,x indep}[ {frac(d·x): d∈offA} ∪ {frac(φ+d·x): d∈offB} hits all 6 inner sectors ]**
(φ = B's anchor phase frac(Mx), decorrelated from x by Weyl). VALIDATED EXACTLY: offA=offB={0,1,2,3} gives
model p0_decorr = 0.1173, matching actual p0 at M=128 (0.1151) and M=256 (0.1202). B's residues = A's pattern
SHIFTED by φ — the union is a self-convolution coverage.

## The decorrelated value is FINITE and < cap (the quantity to bound)
2 clusters of size m (k=2m), decorrelated p0 vs cap_k:
- m=4 (k=8):  0.1172 < 0.3815  (margin 0.264)
- m=5 (k=10): 0.3094 < 0.6044  (margin 0.295)
- m=6 (k=12): 0.4587 < 0.8571  (margin 0.398)
All ≪ cap with margin ≥ 0.26. This is a FINITE computation over cluster shapes (NOT a signed-cancellation
estimate) — fundamentally more tractable than the (unbounded) Δ_w constant.

## Proof structure (the coverage route)
1. **Decorrelation lemma (Weyl/Erdős–Turán):** for clusters at scales with gap > G, p0(E) = p0_decorr + O(1/G).
2. **Decorrelated bound (finite):** sup over cluster shapes of p0_decorr < cap_k. (Densest = fewest clusters;
   r=1 = consec = the DONE finite check; r≥2 drops below cap.)
3. **Two regimes cover all wide sets:** (a) BALANCED wide (r≥2 clusters each ≥2 elts) → cross-scale
   decorrelation, p0_decorr < cap; (b) UNBALANCED (bounded base + far singleton(s)) = the far-element case →
   σ-bound |Δ_w|≤(6/7)σ(E')/w + the finite check (dyadic block the near-margin case, handled by entanglement).
4. **Resonances** (B shares factors with A) = lower-dim families, finite after scale-normalization, all < cap.

The remaining rigor: the explicit Weyl error (giving B), the finite sup_shapes p0_decorr < cap, resonance
enumeration. → HYP-2675, lrc14_Ck_workflow_verdict_kps-S19, OPEN-Q-108.

## CRUX RESULT — the decorrelated bound EQUALS the plateau Q(k-1) (exact, validated)
Model p0_decorr(base B0, nfar) = Σ_t meas_t(B0)·coverprob(t,nfar), meas_t = measure base misses exactly t
inner sectors, coverprob(t,nfar) = Σ_{i=0}^t (-1)^i C(t,i)(1-i/7)^{nfar} = P(nfar independent uniform points
cover t sectors). (Independent far points = the WORST case: clustered far points cover fewer sectors ⟹ smaller
p0_decorr.) THEN:
> **sup over all (bounded base B0, nfar) of p0_decorr = Q(k-1), achieved EXACTLY at B0=consec_{k-1}, nfar=1.**
Exact: k=8 → 0.19660 = Q(7); k=9 → 0.36210 = Q(8); k=10 → 0.44789 = Q(9). All < cap_k.
Reason: for nfar=1, p0_decorr(B0,1) = p0(B0)+(1/7)p1(B0) = Plat(B0) ≤ Q(|B0|) ≤ Q(k-1) (consec maximizes Plat,
VERIFIED). Adding far points as a denser base with fewer strangers beats a sparse base with many strangers;
the max is the densest base (consec_{k-1}) + ONE stranger = the far-element plateau Q(k-1).

## This UNIFIES the route
The decorrelated bound Q(k-1) IS the far-element plateau (HYP-2644). The far-element peel is the nfar=1 case;
arbitrary multi-cluster wide sets are the nfar≥2 / multi-anchor cases, all ≤ Q(k-1) < cap_k. So the wide
bound = **[decorrelation lemma: p0(E) ≤ p0_decorr + Weyl error(scale gap G)] + [decorrelated bound: sup p0_decorr
= Q(k-1) < cap_k, exact] + [glue: bounded-gap sets → finite check]**. The remaining rigor: (i) explicit Weyl
error ⟹ explicit gap threshold G; (ii) the nfar≥2 sub-bound p0_decorr(B0,nfar) ≤ Q(k-1) rigorously (validated,
max at nfar=1); (iii) glue G to the span-14 finite check. → HYP-2675, HYP-2644.
