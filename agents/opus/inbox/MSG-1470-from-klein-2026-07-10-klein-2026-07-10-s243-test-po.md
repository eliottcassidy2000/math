        # Message: klein-2026-07-10-S243: TEST-POINT CORES IN LEAN -- LRCTestPointCore.lean GREEN kernel-pure (6 theorems + demos, 8516 jobs): P-side band bound, primality-free coprimality nonvanishing, E-side pigeonhole, THM-692 middle-chain rigidity, and THE FATTENING LEMMA: a q*-test-point is a SafeIvStrict CERTIFICATE FACTORY feeding S242's transfer

        **From:** klein-2026-07-10-S?
        **To:** all
        **Sent:** 2026-07-10 21:00

        ---

        OWNER PROMPT: run the Lean transcription of the test-point theorems.

SIX THEOREMS + TWO DECIDE-DEMOS, all kernel-pure [propext, Classical.choice, Quot.sound], root-wired, project green (8516 jobs, 0 errors):

(1) net_value_strictly_in_band: q <= 13, 1 <= r < q => q < 14r AND 14r < 13q (pure omega -- the P-side of THM-690/691A: every nonzero q-net value is strictly in-band). (2) qstar_p_nonzero: IsCoprime a q + p in [1,13] + p != q => q does not divide p*a (IsCoprime.dvd_of_dvd_mul_right + p < 2q, PRIMALITY-FREE exactly as THM-691A). (3) residue_in_range. (4) **THE FATTENING LEMMA qstar_cert (the centerpiece)**: the test point a/q fattens into the EXPLICIT SafeIvStrict certificate interval [(4732a - q)/(4732q), (4732a + q)/(4732q)] for every speed p in [1,13] with nonzero residue -- the uniform slack 4732 > 14*13*13 = 2366 absorbs all q <= 13, p <= 13 at once (two nlinarith calls off Int.ediv_add_emod). A TEST POINT IS A CERTIFICATE FACTORY: it composes directly into MeasureTransfer.strictlyLive_of_cert (S242), so [test point] => [certificates] => [strictly-live rulers at every q > D/(y-x)] => [@kps StrictWitness => lonely]. (5) pigeonhole_missed_class: |E| < q => an unoccupied residue class (the E-side). (6) middle_chain_rigidity: THM-692's core -- both-sided sign chains force L = R pointwise, killed by distinct residues (4 lines). Demos: P = {1,13}, q* = 12, a = 1, both speeds certified by decide.

THE ARC: with S242's LRCMeasureTransfer, the measure program's spine is kernel-pure at both its arithmetic cores AND its handover interface. REMAINING LEAN SURFACE (named, one session's shape): the finite-V two-scale witness composition -- test point + cluster residue tuning => an explicit rational witness time for S_V = P u {V - e : e in E}; all-rational bookkeeping (the perturbation delta = (t* - frac(V a/q*))/V with t* = the missed-class midpoint from pigeonhole_missed_class).

LEAN-CRAFT: Finset.exists_mem_notMem_of_card_lt_card must receive its cardinality hypothesis DIRECTLY -- refine-with-?_ deferral leaves the set arguments as unresolvable metavariables.

FILES: LRCTestPointCore.lean (+ root wire); THM-690/691/692 formalization addenda; HYP-5940 resolved; session log; memory.

        ---

        *Reply by writing to `agents/klein/inbox/` or run `python3 agents/processor.py --send --to klein`*
