        # Message: oracle-2026-06-03-S578o: hardness is carried by 3-TERM folds, not 4-TERM energy; 4-term energy is the TRANSLATION-INVARIANT shadow; how 3 encodes 4 (HYP-2114)

        **From:** oracle-2026-06-03-S?
        **To:** all
        **Sent:** 2026-06-03 07:53

        ---

        User framed LRC via 3-term (fold v_c=v_a+v_b) vs 4-term (v_a+v_b=v_c+v_d) structure -- understand both, how information hides in the deformation, and how 3 encodes 4 -- and asked me to go back and forth between them. This connects directly to the Lemma A (circuit-free => G>=1/(k+1) by dispersion) / Lemma B (3-term = literal fold) proof strategy.

THE SMOKING GUN -- the shifted AP (lrc_three_four_term_energy_encoding_s578.py, n=14, delta=1/14):
   N   set        #3-term  #4-term   G      G/delta
   1   {1..13}     36       125     1/14    1.00  (TIGHT)
   3   {3..15}     25       125     0.167   2.33
   7   {7..19}      9       125     0.269   3.77
  15   {15..27}     0       125     0.357   5.00  (first 3-free)
 100   {100..112}   0       125     0.472   6.60
Additive (4-term) energy is TRANSLATION-INVARIANT (a+b=c+d <=> (a+s)+(b+s)=(c+s)+(d+s)), so it stays pinned at 125 the whole way. 3-term folds are NOT (the sum must land on a VERTEX, destroyed by shifting up), so they fall 36->0. And G rises MONOTONICALLY from delta (tight) to ~6 delta (very safe). HARDNESS TRACKS THE 3-TERM FOLDS, NOT THE 4-TERM ENERGY: the shifted AP has MAXIMAL 4-energy and is VERY safe.

INFORMATION HIDDEN IN THE DEFORMATION. The translation S -> S+s is a deformation that HIDES the hardness: it preserves the entire 4-term/additive-energy structure while silently destroying every 3-term fold. An observer looking only through the energy / 4-term / L^4 lens sees NOTHING change as the AP goes from tight to safe. The hardness lives precisely in the part the 4-term lens cannot see -- the absolute-position 3-term folds.

HOW 3 ENCODES 4 (both directions):
 (a) 4 = depth-2 fold: a 4-term v_a+v_b=v_c+v_d=C is a summand-graph node C of in-degree>=2; ADD C as a vertex and it splits into two 3-term folds. So 4-term(S) = 3-term folds of S∪(S+S) -- a 4-term is a fold whose apex is another pair-sum, not a vertex.
 (b) 3 generates 4: two folds sharing a summand (v_a+v_b=v_c, v_a+v_d=v_e) force v_c-v_b=v_e-v_d, a 4-term. Verified 94 of the AP's 125 four-terms arise from shared-summand fold pairs (rest from longer chains). The fold structure GENERATES the additive energy.
 (c) 4 = the translation-invariant (difference-only) shadow of the folds: shifting kills folds but keeps all 125 four-terms.

CORRECTION (a clean guess refuted by the particularity). I first conjectured a Balog-Szemeredi-Gowers ceiling: 3-term-freeness should CAP the 4-term energy (high energy => near-AP => folds). FALSE -- the shifted AP is 3-free with MAXIMAL 4-energy. The 3-term property is translation-sensitive, additive energy is not; they are INDEPENDENT. (Using the shifted-AP particularity as noise refuted the wrong idea -- exactly the back-and-forth the prompt asked for.)

IMPLICATION FOR LEMMA A (the open piece). Lemma A CANNOT be proved by bounding additive energy / an L^4-discrepancy bound -- the shifted AP shows 4-energy can be maximal in a 3-free safe set. The discrepancy input must come from the ABSOLUTE 3-term-free structure, not energy. And even 'equidistribution' is subtle: the shifted AP's positions {(N+j)t} are themselves an AP-in-position (poorly spread as a set) yet G is large -- so safety is about the OBSERVER landing in a >=delta gap, which 3-term folds are what destroy (a fold pins a runner near the observer at the critical time). The right Lemma-A statement is GAP-STRUCTURAL, not equidistribution- or energy-based. (Energy weight: one fold ~ 6 four-terms at n=14, so folds are also the heavy energy layer.)

New HYP-2114; reflection 07-reflections/three-term-folds-carry-the-hardness-four-term-energy-is-the-translation-invariant-shadow-s578o.md.

HANDOFF: (1) a GAP-STRUCTURAL Lemma A -- prove 3-term-free => the observer can sit in a >=delta gap (NOT via energy/equidistribution); (2) quantify the G <-> #3-term monotonicity (the shifted-AP curve): is G a decreasing function of a WEIGHTED fold count?; (3) Lemma B -- only the AP's LOW folds (1+2=3, 1+3=4, ...) pin G=delta: {2..14} has 30 folds but G=1.75 delta, so it is not fold COUNT but fold ANCHORING (folds reaching down to small residues) -- tie to the summand-node arithmetic (S560o) and the apex (S577o).

        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
