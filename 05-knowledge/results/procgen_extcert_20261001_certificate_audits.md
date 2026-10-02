# External certificate audits (collatz-procgen-20260922, 2026-10-01/02)

**Owner request (2026-10-01):** "download the LRC(14) and other certificate audits".

**Status:** a living record. Sections marked PENDING are being worked by the lrcaudit lane and will be updated here.

## Downloads

All downloads are in `scratch/lrc_certs/`, which is not committed. Each was MD5-checked against its Zenodo record.

| Package | Repo meaning | Content |
|---|---|---|
| Zenodo 22066772, Allikvere, "Fourteen lonely runners" | **repo LRC(14)**: 14 runners, 13 speeds, gap 1/14 | 111 prime gates; 174.6 MB package, paper, code guide |
| Zenodo 22667683, Allikvere, "Fifteen lonely runners" | **repo LRC(15)**: 15 runners, 14 speeds, gap 1/15 | 71 gate archives (4.86 GB), manuscript source, code and logs, the author's audit summary |
| Zenodo 21739363, Allikvere, "The edge multiset dimension of hypercubes" | arXiv:2608.09983 | certificates and result files |
| `results/result_10` from github.com/vzsky/13-lonely-runners at commit 755b116 | the LRC(9) input of both LRC packages | the Sungkawichai-Trakulthongchai successive-lifting log |

The two LRC papers form arXiv:2609.02604 ("Fourteen and fifteen lonely runners"). Their convention is LRC(k) = k moving runners, so the paper's LRC(13) is the repo's LRC(14).

**Policy.** Only data was read from the packages. None of the authors' code was run, including their audit scripts. Every check below uses the orchestrator's own code.

## 1. Edge multiset dimension (arXiv:2608.09983): VERIFIED

Script: `04-computation/experiments/procgen_extcert_20261001_edim_paper_check.py`, with the C engine `procgen_extcert_20261001_edim_q5.c`. Output: `procgen_extcert_20261001_edim_paper_check.out`.

**Infinite cases.**
- `edim_m(Q_d) = infinity` for d = 2, 3, 4, by checking every subset.
- For d = 5: none of the 698,635 Aut(Q_5)-orbit representatives with |S| <= 16 resolves the edges. The representative system was checked two ways: against the Burnside count, and by orbit sizes summing to C(32, a). Complements cover |S| >= 16.

**The paper's intermediate Q_5 data, recomputed exactly.**
- R4 = 14,887,680.
- 3,056,640 directional survivors.
- Survivor orbits per size: 12, 38, 73, 30, 153, 184, 153, 30, 73, 38, 12.
- The 796 archived representatives meet every survivor orbit exactly once.
- All 796 archived collision witnesses are genuine.

**Finite cases.** Every archived resolving set resolves:
- the Q_6 descent chain from 29 down to 15, where the last set is THM-4525's;
- the Q_7..Q_10 sets of sizes 63, 115, 246 and 492.

**Rational bounds.** The 40 rational bounds U_d < 1 (11 <= d <= 50) are fractions below 1. They are not needed after THM-4534.

**Packaging slips in the README.** Neither is mathematical.
- The size-15 Q_6 set is in `q6_min_result.txt`, not in `q6_quick.txt`, which holds a 29-set.
- The 796 orbits are Aut(Q_5)-orbits, not Aut(Q_5) x complement orbits. Under Aut(Q_5) x complement there are 402.

Recorded in THM-4525 (`external_audit`) and THM-4534.

## 2. Fifteen lonely runners (repo LRC(15))

### 2a. Analytic half: VERIFIED

Script: `procgen_extcert_20261001_lrc15_flagbound_check.py`, output `.out`. This covers Section 3, the lattice flag bound.

- **Constants (exact).**
  - `K_14 = prod c_i` matches the paper's rational. Its minimum ratio is 847/513 > 4/3.
  - `A_14 = 206083383792625.44...`.
  - The threshold is `14 log(A_14/28) = 414.7793645567`.
  - The target is 414.7793645567 - log 360360 = 401.98450574593444.
  - Table 1 for n = 13 and 15 is reproduced: 341.03 and 497.03.
- **Lemma 3.6 (product minimization), by brute force.** All C(25,13) = 5,200,300 active sets were enumerated. Of them, 4,096 are feasible vertices, and the minimum is attained exactly at t = c. 4,096 = 2^12 is the number of block structures, consistent with the paper's vertex description.
- **Lemmas 3.1 and 3.2 (ellipsoid inside the zonotope; exact covolume 1/(S F_n)).** Checked on 300 random primitive speed vectors with n = 3..7.
- **Lemma 3.7 (shape bound).** The global minimum is H = 798.3079 > 784, at (1, t0, ..., t0) with t0 = 0.078207. This agrees with the paper's one-variable reduction. The displayed bound 796.4 holds exactly.
- **Prime mass.**
  - The 71 gate primes are distinct primes > 15, from 89 to 569.
  - Their sum is `sum log p = 408.8233173785861 > 401.9845`, a margin of 6.84.
- **The shift lemma is load-bearing.** The 52 strict `J(14,p) = {}` gates alone give only 310.0034. The 19 gates closed with the few-exception shift lemma (p = 89..241 except 239) are therefore essential.

### 2b. Read by the orchestrator

- **Lemma 2.1, prime-divisibility criterion, including (ii).** Correct.
- **Lemma 4.4, few-exception shift.** Correct and elementary.
  - The shifts `t0 + q/d` fix the d-divisible block.
  - Each exceptional speed visits a full d-grid.
  - An open arc of length 2/15 contains at most ceil(2d/15) grid points.
  - So `e ceil(2d/15) < d` leaves a good shift.
  - It assumes LRC(m) for m <= 12.
- **The rest of Section 3.** Read; no gap found.
  - The prefix inequalities via Giri-Kravitz's Lemma 3.3. The induction and its angle-minimization step were checked.
  - The KZ ratio 3/4.
  - The forced divisor lcm(2..15) = 360360.
  - Corollary 3.10.

### 2c. Gate archives: integrity and internal consistency VERIFIED

Script: `procgen_extcert_20261001_lrc15_archive_check.py`, which streams every archive. Output `.out`; 568 s.

- **Integrity.**
  - All 71 archives match `SHA256SUMS.txt`.
  - All 21,305,138 manifest entries match the archive contents, by count and by an order-independent SHA-256 digest.
  - The three tool hashes are identical across the 71 gates: `bgk14` cc275840..., `cascade_k14` 807b250e..., `tight_lift_15` ef8bea4d....
- **Gate summaries.** Every SUMMARY.json has k = 14, p equal to the file's prime, GATE_CLOSED and smax = 12.
  - The variant is "decomposition" exactly when no cover on at most 12 classes exists.
  - The claim types match the paper's table: 52 strict gates and 19 shift gates.
- **Kill records.** There are 21,174,429 persistent orbits in all, one kill record each.
  - Every strict gate has the single orbit (1, 2, ..., 14) with counts (3, 5, 15, 7, 0, 0), as in Section 4.3.
  - In every record of every shift gate, improper_after_neargcd = 0.
  - Before the shift lemma, up to 17,152 lifts per orbit remain improper after the gcd clause (p = 131).

**Limit.** The kill records store counts, not the lifts. The shift-lemma applications and the covering searches can therefore only be confirmed by recomputation; see 2e.

### 2d. Dependency chain

The proof uses the following as inputs:
- LRC(m) for m <= 12:
  - Rosenfeld: the paper's LRC(7), i.e. 8 runners.
  - Trakulthongchai: the paper's LRC(8) and LRC(9).
  - Sungkawichai-Trakulthongchai: the paper's LRC(10..12).
- LRC(13) in the paper's convention, which is the repo's LRC(14) from the fourteen-runner package (Section 3).

**The LRC(9) link.** The paper reports a bug in the intersection step of the published LRC(9) program. It relies instead on the ST26 log `results/result_10`. The orchestrator fetched and checked that log:
- 51 primes in [19, 317]; the absent ones are 23, 29, 31, 37, 41, 43, 47, 61;
- every prime closes with final lifted set S.size() = 0;
- `sum log p = 257.8701 > log B_9 = 9(8 log 45 - log 9) = 254.3047`, a margin of 3.57.

### 2e. PENDING (lrcaudit lane, round 2)

- An independent recomputation of `J(9,p) = {}` for the 51 primes.
- A sampled independent recomputation of the level-15 kills, including the shift-lemma applications.
- A full independent recomputation of the cheapest strict gate, p = 239.
- A reading of Lemmas 4.1-4.3.

**Verdict so far.** The analytic reduction, the prime mass and the certificate bookkeeping hold. The computational content of the gates is not yet independently reproduced.

## 3. Fourteen lonely runners (repo LRC(14)): PENDING

The lrcaudit lane, round 1, is checking:
- the mathematical chain;
- integrity;
- the mass sum;
- an independent recomputation of small gates.

**Checked so far by the orchestrator** (`procgen_extcert_20261001_lrc15_flagbound_check.py`, section F):
- **The gate table.** The package's gate table (`paper.tex`, Table 6.1) lists 111 distinct primes:
  - 5 small gates;
  - the 47 consecutive primes 199..479;
  - a tail of 59 primes, 487..877.

  The block sums 24.9221, 272.6869 and 383.9202 are reproduced, as is the total 681.5292.
- **The mass against the MSS threshold.** Without the extra gate p = 877, `sum log p = 674.7527 > log B_13 = 13(12 log 91 - log 13) = 670.3497`, the MSS threshold used by that paper.
- **The alternative route via the newer flag bound.** The fifteen-runner paper's Theorem 3.8 also re-derives the 14-runner threshold, n = 13 in Table 1: 341.03 against 670.35. By its Remark 3.9, the first 59 gates (through p = 523, with `sum log p = 341.1747 > 341.0320`) would already suffice. The orchestrator reproduced both numbers.
- **The n = 13 shape bound.** It needs R_13 < 101/200, i.e. H > (2600/101)^2 = 662.680. The multistart minimum is 664.926, so it holds but with a thin margin (0.34%).

**The package's own scope note.** The master summary (`LRC13_MASTER_SUMMARY.json`, written in Estonian) is candid about how independent its checks are: "the only truly independent recomputation is the l2-spot (about 500 rows); the rest is an internal consistency check of the same binary's artefacts". This is exactly the gap the lrcaudit lane's independent recomputation is meant to address.
