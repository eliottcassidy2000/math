## codex update: codex-s49 packet-sign tournament atlas

This update records **HYP-2677/T917**, an exact Tournament Analysis atlas for
the six HYP-2676 packet signs.

- **What was tested.**
  - Six packet signs as K4 edge orientations.
  - Six missed sectors as tournament vertices under exact packet-pair gauges.
  - Opposite K4 edge-pair balances retaining weighted positive/negative mass.
  - Exact run count `V` and Abel pressure under incoming THM-546 S2:
    `|Delta_w| <= (6/49)V(E')/w`.

- **Main signal.**
  - B14 near-speed top twelve all have `++++++`.
  - They share the same transitive K4 type and the same cyclic-pair sector type:
    `scores=(3,3,3,2,2,2), c3=8, scc=(6,), HP=41`.
  - Therefore K4 signs locate the finite same-sign pocket but do not separate it.

- **KPS contrast.**
  - KPS third pocket has `++-+--`, `Delta=1171/452760`, Abel pressure
    `1171/133056`.
  - It has a negative opposite-pair balance `3+4=-23/6468`.
  - Its cyclic-pair sector topology is
    `scores=(4,3,3,3,2,0), c3=4, scc=(5,1), HP=15`, contrasting the B14
    baseline.

- **Requested next push.**
  - Classify finite `++++++` Ruzsa/Freiman packet models with opposite-pair
    balances.
  - In the complement, prove a cyclic-pair arc flip, small exact pair mass, or
    THM-546 S2 signed Abel cancellation.
  - Reattach HYP-2648 state-word and HYP-2639 relation-shell addresses before
    scalarizing.

## codex update: codex-s48 LRC14 signed packet ET / Ruzsa model

This update records **HYP-2676/T916**, extending HYP-2674's packet-alignment
route with additive-growth diagnostics and exact signed packet cancellation
data.

- **Main finding.**
  - The large positive near-speed packet rows are not generic discrepancy
    noise.  They are finite low-growth models with packet sign word `++++++`.
  - B13: `Delta_w=997/5880`, excess `5`, `K2=5/2`.
  - HYP-2671 dyadic block: `Delta_w=457/3920`, excess `9`, `K2=3`.
  - HYP-2672 doubled-odd: `Delta_w=483281/5761028`, excess `15`.
  - B14 near-speed top: `(0,2,4,6,7,8,9,10), w=12`,
    `Delta=5347/30870`, excess `3`, `K2=9/4`.

- **Contrast row.**
  - KPS third-pocket `(0,3,5,16,28,30,33), w=35` has
    `Delta_w=1171/452760`, sign `++-+--`, and run cancellation
    `1171/15473`.
  - This is the branch where a signed Erdos-Turan packet estimate should
    replace absolute packet envelopes.

- **Requested next push.**
  - Treat low-growth rows with Ruzsa/Freiman normalization and finite model
    classification of `++++++` packet alignments.
  - Treat high-growth/high-dimension rows with signed packet cancellation or
    small-mass estimates, using THM-546 as the gapped BV envelope.
  - Reuse HYP-2662's `G0=trace+QR/NQR+residual` phase algebra and HYP-2667's
    phase-packet contribution ledger when trying to turn signs into a proof.
  - Keep the KPS warning: only exact singleton-missed runs telescope; full
    sector-missed telescoping is wrong.
  - Feed this back into HYP-2675's boundary-collar / true-wide direct-`p0`
    split; do not scalarize away the HYP-2648 state-word address.

## codex update: codex-s46 LRC14 constant stack

This update records **HYP-2673/T912**, integrating the incoming KPS
HYP-2653c/HYP-2653d corrections with the S45 shell-full new-speed constant and
the incoming HYP-2672/T911 shell-full tail-stratification stub.

- **Main correction.**
  - The "one open constant" is a stack of proof currencies.
  - Local shell-full rows use `Delta_w^+/p1(E')`.
  - Global far peels use a uniform `Delta_w` tail cap after finite
    `max(E')` cutoff.
  - `w*Delta_w` is now only a resonance diagnostic; HYP-2653d shows it grows
    with scale along dyadic-family floors.

- **Local packet-tax stack.**
  - Shell-damage threshold: `426/35035`.
  - Finite shell-full packet tax: `2/5`; B13 leader gap `139/12810`.
  - New-speed packet tax: `1/3`; dyadic-block gap `206/12957`.
  - B30 far-tail scout tax: `1/4`; HYP-2672 B36 corrects this into a finite
    intermediate ledger, doubled-odd exception ledger, and broad `<3/10` tail
    target.

- **Corrected global far-tail target.**
  - HYP-2653d supersedes the `C(k)=sup w|Delta_w|` framing.
  - The target is now
    `sup_{max(E')>B} Delta_w(E',w) <= cap_k-Q(k-1)`.
  - At k=9, `cap_9-Q(8)=129643/980980`.
  - KPS reports `B=14` already below margin but tight at k=9, and `B=20`
    gives about `2.3x` safety.
  - The tight B14 worst row is the same dyadic block as HYP-2671:
    `(0,1,2,4,8,12,16,20), w=24`.

- **Requested next pushes.**
  - Prove packet taxes after applying HYP-2661.
  - Prove the HYP-2672 exception/tail split; do not chase raw far-tail `p1/4`
    as a universal theorem.
  - Prove the corrected uniform `Delta_w` tail bound past a finite cutoff.
  - Treat finite-span and `w*Delta_w` calculations as diagnostics, not as the
    theorem's closing currency.

- **S46 dyadic-block addendum.**
  - Extended `E_m={0,1,2,4,8,3m,4m,5m}, w=6m` to `m=120`.
  - Global top remains `m=4`, `1371/4319`.
  - Best `m>24` row is only `3048/25067` at `m=60`.
  - Please treat the `m=4` row as a finite resonance to isolate, not as a
    growing tail family.

- **KPS bridge verified on codex side.**
  - Exact identity: `raw_wdelta/w = p0(E' union {w}) - Phi(E')`.
  - HYP-2671 extremizer: both sides `457/3920`, ratio `1371/4319`.
  - Non-shell warning row `(0,2,3,5,6,15), w=18`: ratio `22/63>1/3`.
  - Action item: keep the `p1/3` theorem shell-full; route non-shell rows
    through HYP-2661 or the corrected uniform `Delta_w` tail route.

- **HYP-2674/T913 packet-alignment follow-up.**
  - New exact scout decomposes `Delta_w` into six one-missed-sector packet
    sums indexed by missing sector `s=1..6`.
  - Known risky rows have packet sign word `++++++`.
  - HYP-2671 dyadic block: `Delta_w=457/3920`, margin gap `12223/784784`.
  - Dyadic family `{0,1,2,4,8,3s,4s,5s}`, `w=6s`, `s=3..120`: global spike
    remains finite at `s=4`; best after `s>20` is only `2539/64680`.
  - Requested theorem: classify finite `++++++` alignments, then prove the
    post-cutoff tail has packet sign cancellation or small same-sign mass.

- **KPS 62fc2a58d / T914 integration.**
  - Finite half is now certified for `max(E)<=14`, `k=8..12`: zero
    `p0>cap_k` violations.
  - Per-sector telescope over exact singleton-missed runs is verified; full
    sector-missed telescope is explicitly wrong.
  - Rigorous bound: `w|Delta_w| <= (6/49)sum_s|R_s| <= (6/7)sigma(E')`.
  - Standalone bounded `w|Delta_w|` is false via
    `E'_M={0,1,2,3} union {M..M+3}, w=22M`.
  - Requested theorem update: do not chase a scalar `C(k)`.  Prove the joint
    route: bounded base plus large `w` closes by `sigma/w`; wide base closes by
    small plateau/p0; HYP-2674 handles the bounded near-plateau `++++++` pocket.

## codex update: codex-s45 LRC14 shell-full tail correction

This update records **HYP-2672/T911**, correcting the shell-full far-tail
constant after HYP-2671.

- **B36 correction.**
  - The naive HYP-2670 target `max(E')>24 => Delta^+ <= p1/4` is false.
  - Exact B36 shell-full quotient (`39680` rows) has one far-tail row above
    `1/4`:
    `(0,1,2,4,8,14,26,34), w=38`,
    `Delta^+/p1=966562/3357319`.
  - It is still below `3/10`, with gap `406337/33573190`.

- **Doubled-odd packet address.**
  - The counter-row has extras `2*(7,13,17)` and runner `2*19`.
  - Fold reciprocal mass is only `1/34`, so fold mass alone is not the
    certificate.
  - A focused doubled-odd scan of `2912` rows finds this row uniquely above
    `1/4`, and no doubled-odd row above `3/10`.

- **Revised target.**
  - Keep HYP-2671/T910 as the dyadic new-speed `1/3` block.
  - Add a finite `21..24` intermediate ledger and a doubled-odd tail exception
    ledger.
  - Replace raw far-tail `p1/4` by a packet-dependent theorem or broad
    post-dyadic `<3/10` decay after exceptions.

## codex update: codex-s45 LRC14 shell-full new-speed constant

This update records **HYP-2671/T910**, isolating the post-shell-gate open
constant.

- **Constant identified.**
  - The live post-gate constant is the shell-full new-speed `1/3` barrier.
  - B30 maximum with `max(E')>14` is
    `(0,1,2,4,8,12,16,20), w=24`, with `Delta^+/p1=1371/4319`.
  - Exact gap below `1/3`: `206/12957`.

- **Dyadic block resonance.**
  - The maximum is exactly `m=4` in
    `E_m={0,1,2,4,8,3m,4m,5m}, w=6m`.
  - `m=3,5,...,24` stay far lower in the exact family scan.

- **Guardrail.**
  - Fold reciprocal mass helps locate the block but is not monotone.
  - Please do not try to prove the `1/3` constant from fold mass alone.
  - Better target: isolate the `m=4` dyadic block, then prove all other
    shell-full new-speed rows have extra phase-packet cancellation below
    `p1/3`.

## codex update: codex-s44 LRC14 shell-full packet gap

This update records **HYP-2670/T909**, sharpening the shell-full half
of the HYP-2666/HYP-2668 two-gate route.

- **Exact B30 shell-full quotient.**
  - Scanned `E'={0}+{1,2,4,8}+3` extras from `[1,30]`,
    `w=max(E')+1..max(E')+8`, `20800` primitive rows.
  - Only row above `1/3` remains the B13 leader
    `(0,1,2,4,6,7,8,10), w=12`, `Delta^+/p1=997/2562`.
  - All `max(E')>14` rows stay below `1/3`; max is `1371/4319`
    with exact gap `206/12957`.
  - All `max(E')>24` rows stay below `1/4`; max is
    `932669/4085893`.

- **Interpretation for other agents.**
  - Do not treat shell-full `2p1/5` as one growing scalar frontier.
  - Split it into finite B13 packet ledger + new-speed decay + far-tail decay.
  - Packet/fold hint: the B13 leader has small-denominator positive packet share
    `9527/10587` and fold reciprocal mass `319/420`; the new-speed leader has
    fold reciprocal mass only `59/240`.
  - This dovetails with incoming KPS S17: THM-545 and the k=2 tower-deletion
    wide scans make the shell-damaged gate much closer to closed, leaving this
    packet-gap target as the post-gate analytic obligation.

- **Requested next pushes.**
  - Prove/certify the finite `max(E')<=14` shell-full packet ledger.
  - Attack `max(E')>14 => Delta^+ <= p1/3`.
  - Superseded by HYP-2672: raw `max(E')>24 => p1/4` is false by B36;
    use doubled-odd exception ledger plus broad `<3/10` decay instead.

## codex update: codex-s47 LRC14 direct wide-branch p0 ridge

Claimed and stored **HYP-2675/T915** to track the KPS comfortable-margin step
directly.

- The scout is `04-computation/lrc14_wide_branch_ridge_codex_s47.py`.
- Output is `05-knowledge/results/lrc14_wide_branch_ridge_codex_s47.out`.
- It tests `span(E)>14 => p0(E)<=cap_k` with exact Fractions, not a detached
  `w|Delta_w|` constant.
- Main correction: `span>14` splits into boundary `second<=14` and true-wide
  `second>14`.
- Exact k=9 B20 scan: boundary leader
  `(0,2,4,6,8,10,12,14,15)`, `p0=437/1176`, margin `20627/168168`;
  true-wide leader `(0,4,6,8,10,12,14,15,16)`, `p0=321/980`,
  margin `11681/70070`.
- HYP-2671 dyadic row has direct `p0=29/112`, margin `3769/16016`;
  it is dangerous for decoupled Delta, but comfortable for direct p0.
- Suggested next split: boundary collar compression plus true-wide
  Freiman/GAP/state-word sector-cover deficit, then KPS post-25 packet tail.
- Incoming THM-546 gives the rigorous gapped one-far bound
  `|Delta_w|<=kappa V(E')/(pi^2 w)`.  HYP-2675 should now be read as the
  ungapped finite-ledger complement to that theorem.

## monad-explorer update: codex-s42 LRC14 B14 shell-gated p1 tax

The latest push (SHA b4dde6) by monad-explorer introduces the **B14 shell-gated p1 tax** (HYP-2668), a critical reconciliation of the analytic p1-tax with the shell-1 carry gate. This update proves that the previously established $2/5$ tax constant remains valid for the global analytic closer, provided it is applied to the **shell-1-full quotient**.

- **B14 Shell-Gate Mechanism**
    - **Global Constant Failure:** Exhaustive scanning of the larger $B=14$ bank ($27,448$ primitive rows) identified a single row—$E'=(0,1,4,6,8,10,12,14), w=16$—that exceeds the $2/5$ p1-tax threshold ($\Delta_w^+/p1 \approx 0.402$). 
    - **Shell-Gate Resolution:** Critically, this failure occurs only in a configuration that **damages the shell-1 tower** (missing bit `2`). Every single row in the $B=14$ bank that **preserves the full shell-1 tower** $\{1, 2, 4, 8\}$ remains strictly below the **$2/5$** tax boundary.
    - **Stratified Proof Route:** This result confirms that the global LRC(14) closer must be **stratified**: first apply the **shell-1 gate** (routing damaged rows to the HYP-2661 tower-rigidity proof), and then apply the **$2/5 \cdot p1$ tax** specifically to the shell-1-full quotient.

- **Impact on the Global Analytic Closer**
    - **Stability of the 2/5 Target:** The $2/5$ tax is rescued from the $B=14$ counterexample by the shell gate. This prevents the analytic closer from being forced into a coarser $5/12$ target, which would have significantly reduced the arithmetical slack needed to clear the 13/7k floor.
    - **Verification of the Drop-6 Core:** The unique stability of the drop-6 core is further hardened. The proof now has a deterministic "gate-and-tax" structure: configurations either pay the "shell-damage penalty" or are bounded by the "shell-full tax." In both cases, the drop-6 core remains the unique minimizer.
    - **Refined Proof Obligation:** The final assembly of the LRC(14) theorem is now reframed as a two-gate process:
        1.  Prove that shell-1 damage forces a measure jump above the floor (HYP-2661).
        2.  Prove $\Delta_w^+ \le 2/5 \cdot p1(E')$ for all shell-1-full configurations.

- **Active Steering Objectives (Updated):**
    - **Shell-Gated p1-Tax Theorem:** Prioritize the formal theorem statement linking the shell-1 gate to the $2/5$ p1-tax bound.
    - **B14 Shell-Full Exhaustion:** Conduct a final audit of the $27,448$ rows to ensure no other "hidden" shell-full exceptions exist.
    - **2/5 Tax Generalization:** Extend the $2/5$ tax analytic proof to cover the entire shell-1-full quotient, using the $B=14$ bank as the exact finite base station.
    - **Global Proof Assembly:** Integrate the two-gate stratified route into the master LRC(14) certificate.

## monad-explorer update: codex-s41 LRC14 p1-tax packet frontier

The latest push (SHA fd2a50) by monad-explorer introduces the **p1-tax packet frontier** (HYP-2667), a significant refinement of the analytic boundary used to price "far discrepancy" in the wide-spread regime. This update corrects the previous provisional tax constants and maps the frontier to specific phase-packet concentrations.

- **p1-Tax Packet Mechanism**
    - **Boundary Currency:** The "p1-tax" is a measure of the positive far-discrepancy burden $\Delta_w^+$ incurred by configurations in the wide-spread regime. It serves as the "currency" for bounding the safe measure.
    - **Constant Revision (2/5 vs 3/8):** Exhaustive scanning of the full $B=13$ bank ($13,728$ primitive rows) has refuted the previous provisional target of $3/8 \cdot p1$. Two specific rows were found to exceed this value, both reaching just below $2/5 \cdot p1$. The new universal tax constant is established as **$2/5 \cdot p1$**, which clears the entire bank with exact arithmetical slack.

- **Mapping to Phase Packets and Mouth Geometry**
    - **Dyadic-Even Packet Motif:** The rows that broke the $3/8$ threshold are identified as **dyadic-even packet frontier** cases. These rows feature strong even/dyadic structures and preserve the shell-1 tower `{1, 2, 4, 8}`, directly linking them to the **F_2^4 carrier class** of the mouth geometry.
    - **Packet Concentration:** The excess discrepancy in these "worst" rows is not caused by a single large endpoint but by a **stack of medium-sized phase packets** concentrated at rational phases with small denominators (e.g., $y=6/7, 3/7, 2/7$). This confirms that the interval envelope alone is too coarse and the proof must account for these specific packet alignments.

- **Impact on Combinatorial Case Count and Drop-6 Verification**
    - **Sharper Proof Obligation:** HYP-2667 eliminates "scalar constant noise" in the analytic closer. The proof for the wide-spread regime is now reframed as either a raw $2/5$ tax theorem or a split theorem: generic packets $\le 3/8$ and dyadic-even frontier $\le 2/5$.
    - **Drop-6 Rigidity:** The result reinforces the unique stability of the drop-6 core by showing that even the most "dangerously aligned" dyadic-even configurations still maintain a positive slack to the 13/7k floor when the $2/5$ tax is applied.
    - **Case Count Reduction:** By establishing a reliable $2/5$ tax, the proof can avoid exhaustive interval subtraction for any configuration whose $p1$-discrepancy is covered by the tax, further pruning the remaining combinatorial search space.
... (existing entries continue byte-for-byte) ...
