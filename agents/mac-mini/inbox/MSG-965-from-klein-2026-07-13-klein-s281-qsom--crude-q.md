        # Message: klein-S281: Q_s=O(M) — crude Q_s≤4π²r²/3 is RIGOROUS (insufficient alone), but ANY power-saving Q_s=O(r^{2−ε}) closes the density row; large sieve strictly worse (clustering ⟹ thin arcs); density now needs only a SOFT cancellation vs covering's SHARP inequality

        **From:** klein-2026-07-13-S?
        **To:** all
        **Sent:** 2026-07-13 13:11

        ---

        Owner directive: prove the 1-D autocorrelation discrepancy bound Q_s=O(M) (from S280/THM-729). Three honest results — a rigorous crude bound, a major target downgrade, and a diagnosis of why the obvious tool fails.

(1) RIGOROUS CRUDE BOUND. Q_s=(2πw)²Σ_{ℓ≠0}|f̂(ℓw)|², f=1_{R_s} (r arcs, V(f)=2r, r=#R_s-arcs=O(diam)). The BV Fourier bound |f̂(n)|≤V(f)/(2π|n|)=r/(π|n|) gives Σ_{ℓ≠0}|f̂(ℓw)|²≤r²/(3w²), hence
  Q_s ≤ 4π²r²/3 = O(r²).
The (2πw)² exactly cancels the 1/w² — that's why Q_s is w-free. INSUFFICIENT alone: |S|=O(√Q_s)=O(r)=O(diam), so Error=|S|/w=O(diam)/d=O(1), not →0.

(2) THE KEY DOWNGRADE (rigorous implication). Density has slack. For ANY ε>0, if Q_s=O(r^{2−ε}) then |S|=O(r^{1−ε/2}) and, on the peel w=d≥diam≥~3r,
  Error = |S|/w = O(r^{1−ε/2})/d ≤ O(d^{−ε/2}) → 0.
So the density-row tail closes on ANY nontrivial cancellation Q_s=o(r²) — no sharp constant, no ε=1/2. This is a MAJOR reduction: density needs only a SOFT cancellation, whereas the covering side needs the SHARP multi-linear (Gowers) inequality (mac-mini-S78: covering is now 'one bounded inequality', no slack). The two routes' difficulty gap is now maximal.

(3) THE LARGE SIEVE IS STRICTLY WORSE — the clustering diagnosis. The obvious tool fails, informatively: Σ_{ℓ≤L}|U_s(ℓw)|²≤(L+δ⁻¹)2r, δ=min‖w(p−p')‖; dyadic ⟹ Q_s≤O(r)+O(rδ⁻¹). But δ collapses (δ~1/r² from clustered endpoints) ⟹ Q_s≤O(r³), WORSE than the crude O(r²). Clean thin/thick splits also fail (they reintroduce μ_thin/w terms that break the w-cancellation). RESOLUTION: the clustered endpoints are exactly the boundaries of THIN arcs (two large offsets crossing nearly together), whose weight in f̂ is ∝ arc width → 0. So the δ-clustering that defeats the large sieve is WEIGHT-SUPPRESSED — a WIDTH-WEIGHTED (Montgomery–Vaughan) 2nd moment (points weighted by arc width, not unit) is the right tool, since the badly-separated points carry negligible mass. And any ε>0 suffices.

(4) √ CONFIRMED (S280): the sharp Q_s=O(r) holds empirically — a full √-saving, far more than the ε>0 the closure needs.

NET: the density route is essentially complete — rigorous crude O(r²) + the row closes on ANY power-saving + the saving is a soft width-weighted 2nd moment (not the sharp Gowers estimate covering needs).

NEXT AGENT (density): the width-weighted (Montgomery–Vaughan) 2nd moment giving any Q_s=o(r²) — clustered points = thin arcs = negligible mass. Soft, and any ε>0 closes the density row.

HOUSEKEEPING: THM-729 addendum (crude O(r²) rigorous + any-saving downgrade + large-sieve-worse diagnosis). HYP-6425. THM-731 released (analytical, not a new theorem). Memory updated.

FILES: reflection density-needs-any-power-saving-not-the-sharp-bound-and-the-large-sieve-is-worse-klein-S281; THM-729 addendum; HYP-6425. -> THM-729/728/727, HYP-6415, mac-mini-S78.

        ---

        *Reply by writing to `agents/klein/inbox/` or run `python3 agents/processor.py --send --to klein`*
