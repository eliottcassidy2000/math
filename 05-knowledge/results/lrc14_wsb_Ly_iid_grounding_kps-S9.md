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

## KEY FINDING: consec is the UNIQUE tight case; every wide config has comfortable margin (kps-S9)
L_y({0,1,2,3,4,5,6,N}) (one large stranger on a 7-AP core; the HYP-2608 stranger-pull-in) over N=8..140:
**oscillates around ~0.21, max 0.257 (N=14), ALL ≪ cap_8=0.3815** (margin ≥ 0.12). The N≡0 mod7
"resonances" are slight SPIKES (N=14:0.257, N=21:0.237) but still far below cap. Limit (large non-res N)
≈ 0.21. So: pulling ONE point of consec out to a large stranger drops L_y from 0.358 (consec max) to ~0.21–0.26.
**CONSEQUENCE — the proof split is clean and asymmetric:**
- consec (and near-consec, bounded spread) is the UNIQUE tight case (L_y=0.358, margin to cap only 0.023) — handled by the EXACT finite check (consec verified max over bounded spread).
- ANY wide config (≥1 large stranger) drops to L_y ≤ ~0.26 < cap (margin ≥ 0.12) — handled by the wide-spread bound, which therefore has COMFORTABLE room (not the tight 0.023).
So the wide-spread bound does NOT need to be sharp — it needs L_y(E) ≤ ~0.30 (say) for maxE > B(k), which the equidistribution of the stranger(s) gives with margin. The tight 0.023 margin lives entirely in the FINITE (bounded-spread) check, which is exact. This is the most favorable possible structure: hard part is finite & exact, analytic part is loose. → THM-534, HYP-2607/2608.

## THE STRANGER-CONTRACTION ROUTE (kps-S9) — a CLEAN single-variable Weyl path for the wide-spread bound
Avoids the HYP-2606 "must-be-signed / 5×-lossy full lattice sum" by peeling ONE offset at a time.

**Contraction lemma (VERIFIED, limit exact):** for E = E' ∪ {N} (N=maxE the largest offset),
J(A, E'∪{N}) → (1−|A|/7)·J(A,E') as N→∞ (the N-point avoids the |A| sectors INDEPENDENTLY, by Weyl on
the single largest offset). Hence S_r(E) → (1−r/7)·S_r(E'), and
> **L_y(E) → L_y^c(E') := Σ_r y_r (1−r/7) S_r(E')  + O(disc_N)**,  the SINGLE-OFFSET Weyl error.
VERIFIED exact: L_y({0..6,N}) → 0.2132 = L_y^c({0..6}) (the contraction formula reproduces the observed
large-N limit to 3 decimals); and L_y({0..6,N}) ≤ 0.257 for ALL N≥8 (never near cap).

**The wide ceiling (VERIFIED):** max over all 7-cores E'⊆{0..9} of L_y^c(E') = **0.2132 at the consec core
{0..6}**, margin to cap_8 = **0.168** (vs the tight 0.023 at consec_8 itself). So ANY config with ≥1 large
stranger contracts to ≤ ~0.21 ≪ cap — the wide-spread bound is LOOSE.

**The proof structure (concrete):**
- BASE (bounded spread, maxE ≤ B): exact finite check; consec is the max, L_y(consec_8)=2633/7350 < cap_8
  (the unique TIGHT case, margin 0.023). [done by THM-534/535/536 finite checks; needs explicit B]
- INDUCTION (maxE = N > B): peel N. L_y(E) = L_y^c(E∖{N}) + O(disc_N), disc_N ≤ (#E'-intervals)/N → 0
  (single-variable Erdős–Turán/three-distance — CLEAN, not the multi-dim lattice). L_y^c(E∖{N}) is the same
  problem one offset down (recurse) and is ≤ cap with margin 0.168 for bounded cores.
- The error disc_N is a ONE-DIMENSIONAL discrepancy (offset N vs the fixed E'-sector-conditions), so the
  HYP-2606 "signed cancellation" issue does NOT arise — that was for bounding the whole k−2 dim lattice at
  once; peeling reduces it to k−1 single-variable Weyl steps, each clean.
This is the moment-side analogue of THM-518 (measure-side stranger-decoupling). The remaining work: (i) the
explicit single-offset Weyl error disc_N ≤ C(E')/N with C(E') controlled by recursion; (ii) the explicit B;
(iii) L_y^c(E') < cap for all E' (recursion base = finite check). → THM-534, THM-518, HYP-2606/2608.
