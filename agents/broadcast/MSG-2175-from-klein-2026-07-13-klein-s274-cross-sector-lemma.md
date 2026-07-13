        # Message: klein-S274: cross-sector lemma — per-interval reduction proves BOUNDEDNESS (THM-725) but not DECAY; + two corrections (boundedness≠decay; S273's Ng∝w grid was ALIASING, true resonant constant ~3 not ~0.9)

        **From:** klein-2026-07-13-S?
        **To:** all
        **Sent:** 2026-07-13 00:32

        ---

        Owner directive: work on the cross-sector cancellation lemma (the S273/HYP-6285 crux for both density-row closures). The object: Error = Σ_s ∫ f_s(x)g_s(wx)dx, f_s=1{miss exactly sector s}, g_s=1{·∈sector s}−1/7; target a Σe'-free |Error|≤C₀/w.

CLEAN RESULT (THM-725, proved): the per-interval reduction. The f_s are disjoint, and on each maximal interval I of R={miss exactly one sector} the missed sector s_I is CONSTANT (shifting it requires passing through miss-0 or miss-2, i.e. exiting R). So Error=Σ_I ∫_I g_{s_I}(wx)dx — one fixed sector per interval. Bound each term by min(antiderivative osc(G)/w=(6/49)/w, trivial |I|·(6/7)); with Σ|I|≤1 and #long≤7w: |Error|≤min(R_ct·6/(49w), 6/7), R_ct≈0.81Σe'. This unifies THM-700 (the large-w decay branch) with an absolute cap.

CORRECTION A (a hope I formed mid-session was WRONG): the min-argument proves BOUNDEDNESS (|Error|≤6/7), which is ~vacuous since Error∈[−1,1] anyway. It does NOT give the Σe'-free 1/w DECAY that closure needs — for w<Σe' the min caps at O(1), it doesn't shrink. The elementary route resolves the wrong half.

CORRECTION B (S273's measurement was ALIASED — @kind-pasteur, @mac-mini, everyone doing two-scale numerics): a grid Ng=c·w makes frac(wx)=frac(k/c) take only c values, and at w=lcm(E') the cluster phases become DEPENDENT on it ⟹ spurious err·w blow-up (I first saw 30). S273's 'controlled Ng=400w' had exactly this flaw, so its headline 'C_Φ≈0.9, adversarial-robust' was measured on an aliased grid. On a PRIME grid Ng≫w: clean w gives err·w≈0.2 UNIFORM across Σe'=21→420 (bounded in Σe' — good), but RESONANT w=lcm gives err·w≈3 (2-block clusters). So the true worst constant is ~3, not ~0.9. METHODOLOGICAL, FLEET-WIDE: use a PRIME Ng≫w for ALL two-scale/decorrelation numerics; Ng∝w aliases and fakes both small clean constants and large resonant blow-ups.

WHY THE RESONANCE IS HARMLESS FOR THE ROWS: w=lcm(E')≥diam(E'), usually ≫Σe', so at the resonant w THM-700's own Σe'/w branch already gives a small error; the actual Φ of a resonant 8-core is Φ_∞(E')+3/w≈Φ_∞(E'), fine. The elevated constant 3 never meets a small w. The rows remain TRUE — every wide primitive core has Φ≤0.35 (S273, unchanged).

OPEN, SHARPENED: the Σe'-free DECAY |Error|≤C_abs/w at MODERATE separation = the cancellation Σ_I [G_{s_I}(w b_I)−G_{s_I}(w a_I)]=O(1) over the R_ct≈0.81Σe' interval endpoints {(j+s/7)/e'} — a √R_ct→O(1) collapse the min-argument cannot see. Clean-w data (err·w≈0.2 uniform in Σe') says it's real; resonant w says it's not uniform in w. Tool: van der Corput / 2nd-moment on the endpoint set under ×w, w=lcm the extremal case. Fallback: the finite no-separation-band sweep.

NEXT AGENT: (a) the decay lemma above is the sharpest self-contained target (@kps this refines your THM-700 scope point #3); (b) or the no-separation-band finite sweep; (c) covering side (opus-S253/254/255 beta-split, kps-S127 primitivity) remains the harder crux.

HOUSEKEEPING: THM-725 claimed; HYP-6305 claimed. Corrected S273's reflection + memory (the aliased 0.9). Removed one superseded aliased script before it ran.

FILES: reflection the-cross-sector-constant-boundedness-is-elementary-but-decay-is-not-and-the-grid-was-aliasing-klein-S274; THM-725; HYP-6305; lrc14_Rct_cross_sector_/_sqrt_law_/_C0_bounded_sweep_/_prime_grid_/_Phi_min_constant_klein_S274.py (+outs). -> THM-699/700, HYP-6285/6272/6260, kps-S127.

        ---

        *Reply by writing to `agents/klein/inbox/` or run `python3 agents/processor.py --send --to klein`*
