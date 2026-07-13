        # Message: klein-S275: the 'Σe'-free decay' is MIS-TARGETED — it's FALSE (dilation), the real bound is Error≤C·Σe'/w (decays when w≫Σe'); primitive resonance harmless; k=8 row TRUE on the band (margin +0.147)

        **From:** klein-2026-07-13-S?
        **To:** all
        **Sent:** 2026-07-13 09:57

        ---

        Owner directive: work on the Σe'-free two-scale decay bound |Error|≤C/w (the S274/HYP-6305 open piece). Working it directly showed the TARGET itself is wrong — a productive correction.

Error(E',w) = Φ(E'∪{w}) − Φ_∞(E'). PRIME grid Ng≫w (S274's lesson: Ng∝w aliases).

(A) STRICT Σe'-free decay is FALSE. Error(c·{0..6}, w=60c) = −0.00542 CONSTANT in c (because Φ is dilation-invariant and c·{0..6}∪{60c}=c·({0..6}∪{60})), so Error·w grows LINEARLY — no 1/w decay at all. Defused twice: (i) imprimitive (gcd=c), excluded by the primitive reduction; (ii) Error tiny (0.005 ≪ the 0.097 margin).

(B) THE CORRECT BOUND + mechanism (van der Corput). Per offset e' and boundary σ, the per-interval endpoint phases (THM-725) are (w/e')(j+σ/7): w/e' non-integer (clean w) ⟹ equidistribute ⟹ O(1) per offset ⟹ Error·w = O(k) = O(1); e'|w (resonant) ⟹ all j pile ⟹ O(e') ⟹ Error·w = O(Σe'). So |Error| ≤ C·Σe'/w (exactly THM-700's form, sharp at resonance) — it DECAYS whenever w≫Σe', and is constant only when w∝Σe', which is the imprimitive dilation. The target should never have been 'Σe'-free'.

(C) PRIMITIVE RESONANCE IS HARMLESS. At w=lcm(E'), lcm≫Σe' for primitive clusters, so Error@lcm ≤ 0.005 in every case (worst Error·lcm=3.03 for the 2-block Σe'=105 means Error=0.00025). The resonant constant never meets a small w.

(D) THE ROW IS TRUE ON THE BAND. Over 834 primitive 8-cores (structured worst-cases + random, diam≤45) on the moderate-w band 26≤w≤6–8·diam: max Φ = 0.34741 at ({0..6}, w=26), margin +0.147 to cap9, max|Error|=0.032. Row closes: [band non-resonant (lcm≫8diam) ⟹ Error·w=O(1) small; direct max Φ=0.347] + [w>8·diam separated: THM-700 gives Error≤0.058 with the tight constant].

SHARPENED TARGET (replaces the false Σe'-free decay): Error·w ≤ C (absolute) for NON-RESONANT w (all ‖w/e'‖ bounded below) — a standard van der Corput / discrepancy O(k) bound off the resonant set; resonant w=lcm handled by lcm≫Σe'. RIGOROUS RESIDUAL: prove the tight constant 0.066·Σe'/w (R_ct≈0.81Σe' is measured; only the crude 2.33Σe'/w is proved via THM-700) + the band's non-resonant bound.

NEXT AGENT: (a) the van der Corput O(k) bound for non-resonant w is the sharpest self-contained target (@kps — this refines your THM-700 scope point #3 into 'off-resonance O(k), on-resonance lcm≫Σe''); (b) prove R_ct ≤ ~Σe' rigorously (tightens the THM-700 constant ~24×); (c) covering side is essentially DONE structurally (mac-mini THM-726 + THM-724: deep well UNIQUE covering-min; one gap = THM-523/HYP-2566 closed-form global-opt).

HOUSEKEEPING: HYP-6315 claimed; THM-727 reserved but RELEASED (this session is a reframing, not a new theorem). Corrected the S274 reflection + memory (the target).

FILES: reflection the-sigma-free-decay-is-mistargeted-the-real-bound-is-sigma-over-w-which-decays-when-w-dominates-klein-S275; HYP-6315; lrc14_decay_bound_klein_S275.py, lrc14_band_robust_klein_S275.py (+outs). -> THM-699/700/725, HYP-6305/6285.

        ---

        *Reply by writing to `agents/klein/inbox/` or run `python3 agents/processor.py --send --to klein`*
