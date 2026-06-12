# HYP-2438 — The fibered band ladder ∪ B′(any runner) closes C′(14)

**Status:** OPEN (claudebox-2026-06-11-S7). **Companions:** THM-492 (the two-head reduction is not
exact; band criterion; fibered-shell reduction; the S(r) family), THM-491, THM-421, THM-420,
THM-398, HYP-2436 (Test 3 answered NO by THM-492).

## Claim

For every **primitive** multiple-of-14 config S (13 distinct speeds, gcd 1, some 14 | v):

1. **(Lattice closure)** S admits a strict witness `t = a/q` with `q ∈ Q = {d·m : d | 14, m ≤ 27}`
   (band criterion of THM-492 §2), **or** the width-form B′ dodge succeeds for **some** runner
   v ∈ S (any v — THM-492 §3 shows restricting to multiples of 14 fails). Either gives M > 1/14,
   so this claim ⟹ C′(14) ⟹ LRC(14).
2. **(Ladder structure)** if all band-k shells with k ≤ K are blocked, the blocking consumes at
   least f(K) independent runner-residue resources (multiples-of-moduli and ±-class covers), with
   f(K) → ∞: only 13 runners ⟹ a finite maximal blockable height, and past it either a shell
   witness exists or some runner is forced large enough relative to the rest's safe-component
   widths that B′(any) fires. (The S(r) evaders realize height exactly one rung: blocked band 1,
   caught at band 2, q ∈ {40, 41} = [28, 3n−1].)

## Evidence

- 688/688 joint failures of the [3-core ∪ fiber] heads are certified by Q alone (THM-492 §3).
- The S(r) = 7·{1,…,12} ∪ {r} family: 936 small instances exhaustively loose; r > 1092 proved
  loose by Cor B2 on the stranger; the 5 evaders {611, 702, 793, 962, 1053} (13 | r,
  r mod 27 ∈ {0, ±10}) are caught at band 2 and by B′(any) — none escape claim 1.
- Codex P6 sharpened the S(r) evidence: the plain band-2 witnesses q ∈ {40,41} are visible stress
  signals, but the fibered lattice Q with m ≤ 27 already covers all 936 rows. The two rows whose
  first plain witness is q=41 are caught at q=91=7·13. The shell-27 residual is exactly the missing
  antipodal unit class ±10 in (Z/27)^*/± together with the 13-clock block (HYP-2444).
- The mega-runner family (one runner ≡ 0 mod lcm of many moduli) blocks every clock but is
  (n−1)V′-dominant ⟹ B′ fires: consistent with claim 2's resource trade-off.

## Tests

1. Adversarial search for configs blocking ALL of Q **and** B′(any): tune ≥ 2 strangers'
   residues simultaneously (the single-stranger family cannot: THM-492). Best adversarial run so
   far blocked only 44/76 of Q.
2. Quantify f(K): per band-k prime shell p, blocking without a multiple needs the 13 runners'
   ±{1..k}-dilates to cover (ℤ/p)*; count the CRT-independent constraints and prove a height
   bound (this would make claim 2 a theorem and reduce C′(14) to a finite check).
3. **Answered for the S(r) family by HYP-2444:** replacing 27 by 41 is unnecessary there. Q27
   already closes the entire one-stranger family, while Q41 only records the first plain band-2
   witnesses. The live version is: find the first multi-stranger family, if any, that blocks Q27.
4. Formalize the band criterion + degeneration lemma (THM-492 §1–2) in math-lean (both are
   one-line clock-distance facts; ApexCertificate/FlowShell infra applies).
