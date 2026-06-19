# Wide-spread bound grounding: L_y^iid, margins, the multi-scale relation structure (kps-2026-06-19-S9)

Target: prove **L_y(E) ≤ cap_k for all E** (k=8..12) ⟹ LRC(14) (with k≤7 + upstream glue).
THM-534: meas(S7(E)) ≤ L_y(E) [PROVED]; L_y(consec) ≤ cap [VERIFIED]. Suffices: consec maximizes L_y.

## The iid (large-spread) limit and margins (EXACT, k=8)
- `cap_8 = 2243/5880 ≈ 0.3815`.
- `L_y(consec_8) = 2633/7350 ≈ 0.3582` — the MAXIMUM. Margin cap−consec = **2243/5880 − 2633/7350 ≈ 0.0233** (thin).
- `L_y^iid(8) = 40573/823543 ≈ 0.0493` (the relation-free / large-spread limit; S_r^iid = C(6,r)(1−r/7)^{k−1}).
- relation-free E (e.g. {0,7,53,311,1009,4999,131,9001}): L_y ≈ 0.049 = L_y^iid ✓.
- resonant w≡0 mod7 ({0..6,49}, {0,7,14,21,28,35,42,1}): L_y ≈ 0.21–0.23 — intermediate, SAFE (≪ cap, < consec).
- perforated AP {0,2,3,4,5,6,7,9}: L_y ≈ 0.151.

So **L_y ranges from ≈0.049 (iid) up to the consec max ≈0.358**, all ≤ cap_8 ≈ 0.3815. The correction
`L_y(E) − L_y^iid` is maximized at consec (≈0.309), and the whole game is bounding it ≤ cap − L_y^iid ≈ 0.332.

## The structural subtlety (why "spread > B" is the wrong threshold)
- WIDE-spread + RELATION-FREE ⟹ L_y ≈ L_y^iid ≪ cap (the absolute Weyl bound suffices — sparse lattice).
- BUT a wide-spread E can carry a PERSISTENT short relation from a sub-AP: e.g. `{0,1,2,3,4,5,6,N}`
  (the HYP-2608 "stranger-pull-in") has the relation 1+2−3=0 (λ₁=√3) for EVERY N, so its correction does
  NOT vanish as N→∞ — L_y oscillates (resonance dips at N≡0 mod 7 via the 7-vanishing ĉ_T(7m)=0) toward
  L_y({0..6} + one equidistributed point), still < cap.
- So the right notion is the RELATION LATTICE Λ°(E)=e^⊥ (HYP-2606), not metric spread. The correction is
  a SIGNED sum over Λ° dominated by SHORT vectors; the short vectors come from low-dim sub-AP structure.
  The proof must (a) for shapes whose ONLY short relations are "long" (large covolume sub-lattice), use
  the absolute Weyl bound (correction small); (b) for shapes with persistent short relations (sub-AP),
  reduce to the lower-k sub-problem + the equidistributed remainder (Weyl on the free directions),
  bounding the limit < cap; (c) the genuinely-dense AP and near-AP shapes are the bounded-spread finite
  check (where signed cancellation matters, HYP-2606's 5× — but it's a FINITE family).
- RESONANCES (w≡0 mod7) are HELPED by the 7-vanishing (ĉ_T(7m)=0 kills that offset's contribution), so
  they DIP, not spike — consistent (L_y(resonant)≈0.23 < consec 0.358).

## NET (the cleanest framing of the wide-spread bound)
L_y(E) = L_y^iid(k) + Σ_{0≠n∈Λ°(E), all n_j≢0 mod 7} K_y(n) [signed; 7-vanishing restricts support].
The correction is maximized by the densest short-relation lattice = the AP (HYP-2606 F2: AP=min covolume).
Closing it = (i) the AP/near-AP finite check (bounded covolume class) + (ii) an absolute-or-signed Weyl
bound for the rest (large-covolume directions), with the sub-AP persistence handled by induction on k.
The thin 0.023 margin at k=8 means the bound must be near-sharp. → THM-534, HYP-2606/2607/2608.
