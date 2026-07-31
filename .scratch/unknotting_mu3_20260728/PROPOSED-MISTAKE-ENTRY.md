## MISTAKE-ID-TBD — `u(K)<=1` was promoted to `u(K)=1`

- **Status:** DEMONSTRATED; repair and regression prepared under
  `.scratch/unknotting_mu3_20260728/proposed/`; canonical promotion pending.
- **Affected artifacts:**
  - `04-computation/unknot1_decider.py`;
  - `05-knowledge/results/UNKNOT1-DECIDER-STATUS-AND-DESIGN-klein-S691.md`.
- **What:** `TRUE_CERTIFIED` was emitted when a visible crossing change
  reduced to the empty diagram after the greedy input simplifier stalled.
  The change certificate proves only `u(K)<=1`.  Greedy failure is not a
  certificate that the input is nontrivial, so exact equality was
  unjustified.
- **Minimal demonstrated witness in this audit:** the six-crossing PD

  ```text
  [[1,11,2,10],[6,10,7,9],[3,8,4,9],
   [11,5,12,4],[7,2,8,3],[5,1,6,12]]
  ```

  represents the unknot: it is obtained from the one-crossing unknot by
  `R1+,R2+,R2+,R3,R3`.  The input R1/R2 reducer stalls.  Changing crossing
  4 reduces by `R2,R1,R2,R1`, so the original program returns
  `TRUE_CERTIFIED` although `u(K)=0`.
- **First failed implication:** “the current diagram does not simplify
  greedily” implies “the represented knot is nontrivial.”
- **Destroyed coordinate:** the distinction between failure of one
  simplification strategy and topological nontriviality of the knot.
- **Strongest survivor:** the emitted crossing/move sequence is a valid
  certificate of `u(K)<=1`.  It certifies `u(K)=1` when combined with any
  independent certificate that `K` is nontrivial.
- **Repair:** require both (i) an unknot certificate after one crossing
  change and (ii) an independent input-nontriviality certificate.  The
  minimal stdlib repair uses `det(K)!=1` or `sigma(K)!=0`; otherwise it
  returns `UNKNOWN` while retaining the `u<=1` certificate.  A complete
  exact unknot oracle on the input is the preferred future replacement.
- **Regression:** add the six-crossing PD; require repaired verdict
  `UNKNOWN`, a retained upper-bound certificate, and identical normal /
  `python3 -O` results.  The proposed suite passes all 13 original and 4
  new hostile checks.
- **Reproduction:**

  ```bash
  python3 .scratch/unknotting_mu3_20260728/unknot1_true_soundness_witness.py
  python3 .scratch/unknotting_mu3_20260728/proposed/unknot1_decider.py --test
  python3 -O .scratch/unknotting_mu3_20260728/proposed/unknot1_decider.py --test
  ```
