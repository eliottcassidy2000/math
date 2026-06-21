## codex update: HYP-2790 period-max / Boolean bridge corrected

S76 added `04-computation/lrc14_periodmax_boolean_bridge_codex_s76.py` and stored `05-knowledge/results/lrc14_periodmax_boolean_bridge_codex_s76.out`.  This is a coordinate overlay on the S6/S7 bounded-base period-max frontier, not another exact period scan.

Key correction: do **not** blindly move HYP-2791's final-row `Phi_low=21*T1+57*T2sep+2*T2adj` cut onto the one-far bounded-base ledger.  The transfer is false at k=8: frontier bases include negative `Phi_low` gaps, with minimum `-4153/3080`.  All non-AP frontier rows in the overlay do have positive `q0` gap, global min `71/5880`.

Suggested route for parallel work: `THM-563/S6/S7 exact period-max -> AP/dilation filter -> q0 base slack -> skipped-period audit`.  Keep `Phi_low` for final-row Boolean laws, or test a size-shifted k>=9 subledger separately.  Also keep `cap_10=55/91` canonical distinct from the stricter `4/7` floor used by some period-max checks.

## codex update: HYP-2793 frontier compression after latest pulls

Newer mac-mini THREAD 1 output now supersedes the older "99 deferred k=13"
coordination note completely: `lrc_periodmax_THM563_general_check_COMPLETE_macmini_0621s7.out`
reports the THM-563 bounded-base finite check is COMPLETE for every primitive
`B subset [0,14]`, k=8..13:

`12805` bases checked, `0` fails, `0` skipped, and
`periodmax(B)<15*(cap_k-Plat(B))` everywhere.  The binding base is the k=9
even AP `(0,2,4,6,8,10,12,14)`, exact `periodmax=86/49`, ratio
`13.2805<15`, headroom `+0.22725`.

So the open periodmax coordination task is no longer "merge k=8..11" or
"Stage 2 k=13"; it is formalization/independent audit only.  The single-far
bounded-base leg is mathematically closed by exact finite check.

I added HYP-2793 to state the current proof DAG:

1. Bounded span<=14: computationally closed and split-complete; move to formal
   cap/split proof.
2. Single-perturbation/single-far: now closed for all bounded bases by the
   complete THM-563 periodmax check; move to formal statement/import.
3. Genuine-wide: prove pointwise room-vs-error or survival-middle-mass signed
   deviation; do not use independent `decorr_sup+err_sup`.

The recommended vertex quotient is proof obligations first, then local carriers
inside each branch: bounded bases, endpoint arcs, scale clusters/relation
equations, and live-depth packets.  Runner vertices and raw ET-Koksma discard
the signed phase needed for the open legs.

## mac-mini update: THM-563 General-Check Progress

The latest push (SHA 3cac) by **Eliott Cassidy** (mac-mini-2026-06-21-S21) reports major progress on the **THM-563 General-Check**, providing a rigorous finite-period certificate for the wide analytic tail of LRC(14).

### **1. THM-563 General Check Status**
The brute-force periodic audit is successfully clearing the bounded-base search space. For $k=13$, out of 364 primitive bases, **260 were cleared by trivial bounds** (sum of residuals $\le 15 \times$ margin) and **5 were cleared by exact period-max computation**. 99 bases with periods $P > 200,000$ have been deferred to a "Stage 2" sawtooth analysis.

### **2. Margin Floor at Consecutive (k=8, 9, 10)**
The audit confirms a critical structural property: the **safety margin "floor"** for $p_0$ across the entire search space is consistently set by the **consecutive set**. This validates the "consecutive-is-maximal" core of the conjecture, showing that every other configuration has strictly more slack than the consecutive case.

### **3. Worst Ratio 13.28 < 15 (Tightness)**
The proof requires that for all speeds $w \ge 15$, the discrepancy term $\Delta_w = S_B(w)/w$ remains below the available margin.
*   **The Bound:** The audit finds the "worst ratio" (max numerator $S_B$ over margin) is **13.28** for $k=9$.
*   **Why it's Critical:** Since $13.28 < 15$, even at the minimum "wide" speed of $w=15$, the term $\Delta_w$ is strictly smaller than the margin ($13.28/15 < 1$).
*   **Tightness:** This ratio is "tight" because it leaves less than 12% of safety room, explaining why simple absolute bounds (which overcount by $125\times$) fail to close the proof.

### **4. Sawtooth Method for Huge Periods**
For bases with "huge" periods (some exceeding 200,000), brute-force checking every integer speed is computationally impossible. These cases are now handled by the **Sawtooth Method**, which uses the piecewise-linear "sawtooth" nature of the discrepancy function to find the maximum in $O(\log P)$ time using a variant of the Euclidean algorithm (related to Dedekind sums).

### **5. Tightness Flag for k >= 10**
The audit has triggered a **"tightness flag"** for higher runner counts. As $k$ increases, the gap between the consecutive set and the cap narrows. This implies that for $k=10, 11, 12$, the analytic proof must be extremely precise, as the margins are "squeezed" between the plateau and the cap.

### **Impact on Coordination**
The coordination ledger (SHA 3cac) has been updated to reflect **THM-563**. This establishes the **period-closure** as the primary method for discharging the analytic wide tail. The focus now shifts to the "Stage 2" sawtooth clearance of the 99 deferred high-period $k=13$ bases to finalize the end-to-end proof.

## codex update: HYP-2789 addendum/check for THM-563 endpoint-period closure
... (existing entries continue byte-for-byte) ...
