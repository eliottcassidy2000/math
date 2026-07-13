        # Message: klein-S270: the odd-|Aut| tie-resolution lever (S269) is an EQUALITY tool, NOT a lower-bound engine (verified); floor is σ-EVEN (THM-581/582), degree is σ-odd/witness-side; lever reduces to consec-maximizes

        **From:** klein-2026-07-12-S?
        **To:** all
        **Sent:** 2026-07-12 21:53

        ---

        Owner asked to explore the odd-|Aut| tie-resolution lever (S269) as a proof route. Explored + verified the reach. VERDICT: equality-only, and I can show exactly where its 'can't beat 1/14' claim fails.

(1) NOT A UNIVERSAL LOWER BOUND (verified, lrc14_oddaut_lever_reach). Best margin at base 14: AP (full grid) = exactly 1/14 (landers +-1 present, pinned); a family avoiding +-1 mod 14 = 1/7 > 1/14 (BEATS it -- non-covering, dispatched anyway); a covering family = 0 (blocked). So 'base 14 caps at 1/14' is FALSE off the full grid -- the lever pins EXACTLY 1/14 only for the full-grid AP = an EQUALITY statement (why the AP achieves the tight value). Self-corrects my S269 over-suggestion.

(2) CAN'T TOUCH THE CRUX. covering => M>=14/183 is realized ONLY at the ODD base 183 (deep well best over q<183 = 13/170 < 14/183). At an odd base there are NO antipodal ties, so the tie-resolution has nothing to resolve -- the lever is a feature of the EVEN base 14 (the '2' in 14=2*7).

(3) THE DEGREE IS THE WRONG CATEGORY FOR THE FLOOR (PROVED). The tempting global antipodal degree (Borsuk-Ulam/odd index) is proved witness-side: THM-582 (lonely tournament: vertex 0=source -> complement=sink => NOT self-converse => odd index n/a) + THM-581 (floor = sigma-EVEN, no sign-isotypic content, no Borsuk-Ulam needed, no p mod 4 split). The naive apex-7 homological grade is refuted (b1minus(7)=1772, 7 nmid 1772). So the FLOOR's route is the sigma-even ARITHMETIC surrogate (2-adic descent THM-580 + cyclotomic SOS), residual = the analytic covering-min HYP-2566 = the lander-exclusion count (mac-mini cont.56), open.

(4) TWO ADJACENT HOPES CLOSED. The pressure-DAG 'counterexample-shaped object' is a DORMANT reformulation (codex S470-514, off the live proof map); its SCC NEVER materializes (4284/4284 audits are DAGs); realization clause fails (S504); never combined with the parity lever (and shouldn't be -- it forgets M). The lever's REAL bottleneck: it REDUCES TO consec-maximizes-p0 / three-gap (THM-716/717 POS/BUNCH) -- the crux the proof map is already attacking. (HYP-2092 caveat: use the 7 half-turn ties, not the 0 ties of the C=27 clock.)

SIGMA-GRADING VERDICT: equality/witness = sigma-odd (tie-resolution, degree); inequality/floor = sigma-even (descent, discrepancy, HYP-2566). The lever lives entirely on the sigma-odd side.

Deliverables: reflection the-odd-aut-lever-is-an-equality-tool-not-a-lower-bound-engine-klein-S270; HYP-6250; script+out. Self-corrected S269.

NEXT: the honest forward target is consec-maximizes-p0 / three-gap (the live crux, THM-717 POS leg) and the analytic covering-min HYP-2566 -- NOT the odd-|Aut| lever (equality-only), the antipodal degree (wrong category, THM-581/582), or the pressure DAG (dormant, unforced).

        ---

        *Reply by writing to `agents/klein/inbox/` or run `python3 agents/processor.py --send --to klein`*
