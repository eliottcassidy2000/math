# Independent audit: S69 ninja pairs / prehistory of the square-sum numbers (THM-4850, THM-4851, HYP-9447)

Referee: independent adversarial audit (opus, 2026-10-31), worktree at commit 918129c522. Everything here was
recomputed with code written for this audit (`ss_enum.c`, `ss_terminal_big.c`, Python checkers). The session's C
enumerator was compiled only to replay the session's own replay script.

**VERDICT: SOUND WITH CORRECTIONS.** Every number in THM-4850, THM-4851 and HYP-9447 was reproduced independently, and
the proofs are logically complete. The novelty claims of THM-4850 are wrong, though. The bases 42..52, the recursion
from m = 42, the base block [42, 2081], the thresholds [42; 24, ..., 24] and the terminal cycles for N >= 2082 are
already in the repository: `session124_square_prefix_threshold42_20261031`, committed nine hours earlier. The headline
"the square-sum numbers start at 42" also holds only for session124's normalisation. With Gerbicz's original "nice
pair" (same radix 49, same digits 24..72) the pairs start at 41. With the repository's own free nice pairs they start
at 15. The 22/33 and 24-25 readings are selective, and they vanish against the obvious baselines.

## Findings

### MAJOR-1. The threshold-42 numeration and replay are prior art (in the repository) and are claimed as new

- **Where.**
  - THM-4850 lines 45-46: "The threshold 42, the gate, and the start of the numeration at 42 (the base block extended
    from 53 down to 42 and replayed) are new".
  - Note lines 215-217: "New to the repository: the exact ninja threshold 42 ...; the replay of the digit operators on
    the bases 42..52 (the numeration starts at 42)".
  - Note line 125 and THM-4850 line 39 cite only session124's T_l = [53; 24, ..., 24].
  - THM-4851 (c) line 29, note line 147 and HYP-9447 lines 28-29 cite N >= 2621 as the large-N input. HYP-9447 also
    declares "odd N in [63, 2619]" the open window.
- **Evidence.** `05-knowledge/results/session124_square_prefix_threshold42_20261031/` (commit 727bfdf464,
  2026-10-31 00:44 -0600; S69's commits are 09:58 and 10:18 the same day) contains all of the following:
  - README lines 34-35: `SMALL_BASES.json` gives all eleven weak B pairs for 42 <= n <= 52. The audit re-checked
    them with its own checker (`audit_analysis.out`, section I).
  - README lines 39-48: the bases are used "only through 2081", with recursion "for every n >= 2082". The paired
    operators apply at m = 42 ("Strong induction constructs B(n) for every n >= 42"). This is the base block
    [42, 2081].
  - INDEPENDENT_AUDIT.md lines 16-21: "Thus the universal proof is valid at 42".
  - README lines 16-22: S_l = (85*49^l - 1)/2 = 42, 2082, 102042, 5000082, i.e. the numerals [42; 24, ..., 24]
    (checked in `audit_analysis.out` H).
  - README lines 18-19 and 59: for every N >= S_1 = 2082 there is a Hamiltonian cycle beginning 1,3,6,19,45,151 and
    ending 8 (N even) or 15 (N odd), i.e. a terminal cycle.
  - The packets S69 itself lists as context (note lines 23-27) consume these bases:
    - `session124_mixed_radix_square_rays_20261031/TERNARY_SCHEDULE_THEOREM.md` line 53: "Use the certified weak
      bases 42..2081"; `core.py` line 8; `bases.py` line 40;
    - `square_sum_spectrum_overnight_20261031b/README.md` lines 13 and 31.
  - The threshold42 README makes no minimality claim (line 68). So the following are genuinely new: the exhaustive
    non-existence for n <= 41 (the sharpness of 42), the uniqueness of the gate path, and the counts.
- **Correction.**
  1. **THM-4850 Novelty**, replace lines 45-46 with: "Weak B pairs for every n >= 42 are prior work: the bases
     42..52, the validity of session124's paired operators at m = 42, the base block [42, 2081] (cutoff 2081) and the
     thresholds S_l = (85*49^l - 1)/2 = [42; 24, ..., 24] are all in
     `session124_square_prefix_threshold42_20261031` (independently audited there). The 539-output replay
     re-confirms them with new input pairs. New here are three things:
     - the exhaustive non-existence for n <= 41, so 42 is sharp;
     - the uniqueness of the gate path;
     - the exact path and pair counts for n <= 53."
  2. **THM-4850 (3) and note line 125**: add "the threshold42 extension's S_l = (85*49^l - 1)/2 are the numerals
     [42; 24, ..., 24]".
  3. **Note §6**: move "the replay ... (the numeration starts at 42)" to "CITED and reproduced (threshold42)".
     Reword "the exact ninja threshold 42" as "the sharpness of 42 (no weak B pair for n <= 41)".
  4. **THM-4851 (c), note line 147, HYP-9447 lines 28-29**: "N >= 2082: the threshold42 terminal theorem (session124's
     core gives N >= 2621)". The search window 2082..2620 becomes a redundant cross-check, and the open window before
     the pre-registration was odd N in [63, 2081].
  5. Add threshold42 and `session124_gerbicz_reconstruction_20261031` to the context list in note §0.

### MAJOR-2. "Ninja pair" is most plausibly Gerbicz's "nice pair", and the threshold depends on the normalisation; the headline presents one normalisation as "the" start

- **Where.**
  - Note title and line 1; note line 34 (reading); note line 40 ("appears nowhere else in the repository or in a web
    search").
  - Ladder table line 164.
  - THM-4850 title and heading line 14 ("the ninja numbers start at 42").
  - Results INDEX line 3.
- **Evidence.**
  - **Gerbicz's own term and recursion.** "Nice pair" is Gerbicz's term. In mersenneforum "The Square-Sum problem"
    post #22 (2018-01-21), seq0 and seq1 of lengths n and n+1 form a nice pair if every common term sits at positions
    of equal parity. The same post gives the construction "every 49*n+res, where res=24..72" from a base table
    n = 41..2032, with cycles normalised a[1] = 1, a[n] = 8 (n even), 3 (n odd).
  - **The same term in the repository.**
    - The B-base packet S69 cites is titled "...square-sum nice pairs" (`square_sum_B_base_interval_20261019.md`,
      line 1). Its line 5 defines a nice pair with no prefix or endpoint condition; the 1,3 / 8 / 10 "state" is an
      extra condition.
    - `session124_gerbicz_reconstruction_20261031/GERBICZ_RECONSTRUCTION.md` lines 43, 91 and 97 reproduce
      Gerbicz's definition and base [41, 2032].
  - **Thresholds under each normalisation.** Recomputed exhaustively by the audit (`audit_analysis.out` B, B2):

    | normalisation (equal positional parity on [1, n] in all cases) | pairs exist exactly for |
    |---|---|
    | free nice pair (B-base packet line 5) | n in {15, 16, 25, 26, 27} u [29, oo): exhaustive to 41 (28 is a hole), then every weak B pair is one |
    | Gerbicz 2018 (start 1; odd order ends 3, even order ends 8) | n >= 41: exhaustive for 14 <= n <= 50 (witness at 41 in the .out), then Gerbicz's audited base 41..2032 and recursion |
    | session124 weak B (start 1,3; even order ends 8, odd order ends 10) | n >= 42 (THM-4850) |

  - So Gerbicz's original base-49 numeration, with radix 49 and digits 24..72, starts at 41 with base block
    [41, 2032]. That is the minimum possible for his normalisation. Its numerals [41; 24, ..., 24] =
    (83*49^l - 1)/2 = 41, 2033, ... begin exactly at his 2033 = 49*41 + 24. The odd "gate" member ending at 10 is an
    artefact of session124's prefix-protected state, whose odd member does not close (10 + 1 is not a square).
- **Correction.**
  1. **Note title.** "Weak B pairs start at 42 (Gerbicz's nice pairs at 41): thresholds of square-sum pair
     numerations, the prehistory before them, and the base-49 numeration".
  2. **THM-4850 title.** Keep "weak B pairs exist exactly for n >= 42". Qualify the numeration clause as "(for
     session124's normalisation; Gerbicz's nice pairs, with the same radix and digits, start at 41, base block
     [41, 2032])".
  3. **Note line 34.** Replace with: "Reading (INTERPRETATION): the owner's 'ninja pair' is most plausibly Gerbicz's
     'nice pair'. Three normalisations are in use, with thresholds 15 (free), 41 (Gerbicz) and 42 (session124's weak
     B); this note studies the last."
  4. **Note line 40.** Replace with: "The term 'ninja pair' appears neither in the repository nor on the web. The
     attested term is Gerbicz's 'nice pair'."
  5. **Ladder line 164.** "42: weak B pairs for all larger n (Gerbicz-normalised nice pairs from 41)".
  6. **Note §6.** Add Gerbicz's nice pairs and his numeration (radix 49, res 24..72, base 41..2032) as the origin of
     the construction.

### MINOR-1. The "3 with 22 and 33" and "24-25" statistics are selective and have no baseline

- **Where.**
  - Note lines 69-72 and 175-183.
  - Results INDEX line 3: "The first ninja pairs leave 3 through 22 or 33 ...; the consecutive edge 24-25 ... lies in
    most members".
  - The driver asserts the 22/33 pattern only for n in (42, 43) (`ninja_square_sum_number_system_20261031.py`,
    lines 186-187).
- **Evidence.** `audit_analysis.out` G compares the members of all weak B pairs with all paths of the same prefix and
  end:

  | n | members | 22/33 among members | 22/33, baseline | 24-25 among members | 24-25, baseline |
  |---|---|---|---|---|---|
  | 42 | 20 | 80% | 65% | 85% | 79% |
  | 43 | 11 | 100% | 62% | 100% | 72% |
  | 44 | 120 | 57% (mode 13: 41) | 56% | 82% | 69% |
  | 45 | 5 | 80% | 39% | 100% | 77% |
  | 46 | 16 | 6% | 37% | 94% | 72% |
  | 47 | 24 | 50% | 32% | 21% | 73% |
  | 48 | 1994 | 25% | 29% | 81% | 73% |
  | 49 | 179 | 59% | 38% | 79% | 72% |

  - The note picks n = 42, 43, 45 and skips n = 44, where 13 is the most frequent third vertex. 22 and 33 are already
    the most common third vertices of all such paths.
  - The 24-25 edge is in about 70-80% of all paths, because 25 has only the neighbours 11, 24 and 39 there.
  - "3's two largest neighbours" (line 70) is false once the order reaches 46 (3 + 46 = 49). The Longer at n = 45
    already has order 46.
- **Correction.** Replace both passages with the table above. Type them as NUMEROLOGY / weak observation, delete
  "almost always", "most early ninja members" and "two largest neighbours", and update INDEX line 3.

### MINOR-2. A wrong variable in the selector formula, and a loose description of the operators

- **Where.** Note line 115 and THM-4850 line 42: "the switch is at 49 + (-1)^(n+1) c <= d" (in the note, n = 49m + d
  is the target).
- **Evidence.**
  - Session124 README line 39 uses the input index: the extra input vertex m+1 sits at parity m+1, so block c uses the
    Longer iff 49 + (-1)^(m+1) c <= d. (-1)^(n+1) differs from this whenever d is odd.
  - "Low digits from the Shorter and high digits from the Longer" is also inaccurate. Every operator draws each of
    its 49 offset blocks from one member. The paired tables show an output of order 49m + e draws exactly e - 24
    blocks from the Longer (0/1 at d = 24, 48/49 at d = 72, both parities).
- **Correction.** Write "(-1)^(m+1)", with m the input index (session124's sign convention). Replace the low/high
  sentence with: "an output of order 49m + e takes exactly e - 24 of its 49 blocks from the Longer (those c with
  49 + (-1)^(m+1) c <= e) and the rest from the Shorter. Parity compatibility makes the two members interchangeable
  on [1, m]." Also write "98 paired operators (49 per input parity)" for "49 digit operators" (note line 31,
  THM-4850 title).

### MINOR-3. "Verified independently" is overstated, and no certificates are stored

- **Where.**
  - THM-4851 line 28: "verified independently".
  - HYP-9447 line 45: "the program's independent checker".
  - Note line 144: "two independent programs agree".
- **Evidence.**
  - `verify()` in `terminal_cycles_large_20261031.c` (lines 67-73) is a function inside the search program. The 2558
    cycles were not saved.
  - No `.out` records the second program's run on N <= 62; only the `exists` run of the session enumerator is in the
    driver output.
  - The audit has now re-found all 2558 N in 63..2620 with a different program (randomised Warnsdorff with forcing).
    Each cycle was re-verified in Python (`terminal_big_63_2620_receipt.json`; deterministic, certificate SHA-256
    97a7de77...). The audit's two programs also reproduce the complete searches at 30..33 and 39.
- **Correction.**
  - THM-4851 line 28: "verified by the search program's built-in checker; re-found by a different program and
    re-verified in Python by the S69 audit".
  - HYP-9447 line 45: "built-in checker".
  - Note line 144: "two programs of this session (the second run is unrecorded); reproduced by the audit's two
    programs".

### MINOR-4. The replay output was not recorded

- **Where.** THM-4850 (3) and the status line claim "539 outputs verified". The results INDEX line 3 points to the
  replay's `.out`.
- **Evidence.** At 918129c522, `04-computation/experiments/ninja_numeration_replay_20261031.out` is an untracked,
  empty (0-byte) file. The audit re-ran the session's script unchanged with the session's enumerator: "ALL 539
  outputs verified (checks=6340180)" in 4 min 40 s (`replay_session_script_with_session_enum.out`). With the audit's
  enumerator the result is the same (`replay_session_script_with_own_enum.out`). The audit's own replay checked
  15,435 outputs, from 315 input pairs at m = 42..52 including the threshold42 bases, with two checkers
  (`replay_own_pairs.out`).
- **Correction.** Commit the replay output, or cite the audit's files.

### MINOR-5. A stale cross-reference

- **Where.** `terminal_cycles_large_20261031.c` line 1 says "(HYP-9438)".
- **Evidence.** HYP-9438 is an unrelated hypothesis about binary matroids and square classes. The search program
  belongs to HYP-9447.
- **Correction.** "(HYP-9447)".

### NOTES

- **N1. The 32 -> 33 step.** It is the first record of the 2026-10-30 certificate
  `square_sum_spectrum_overnight_20261030/square_bridge.json.gz` (README line 23):
  `{'insert': 15, 'reverse': [[10, 16]]}`. Replayed, it gives the unique 33-cycle (`audit_analysis.out` D). Its steps
  33..42 are all non-direct, which is consistent with (d). The prose description is new; the fact is not. Cite it in
  THM-4851 (b) and note §6.
- **N2. The A071983 range.** THM-4851 (a) line 17 claims agreement with OEIS for N <= 42. The driver asserts
  A071983 only to N = 34 (lines 50-51), and note line 130 says 14..34. The audit confirms paths = A071983 for
  15..42 and cycles = A071984 for 32..42 against the OEIS lists fetched today. Harmonise to "<= 42 (audit)".
- **N3. The "<= 39" table row.** Note line 91 ("0 or 0 compatible") is garbled. Exact data: for n <= 32 both counts
  are 0; for 33 <= n <= 39 exactly one of the two members has no path. Also fill in the path pairs for n = 50..53:
  1090, 1351, 106235, 470586.
- **N4. "The threshold is set by the odd member".** At n = 40 and 41 the gate exists but is parity-incompatible
  with all 11 and all 139 partners. Odd-member scarcity gives n >= 40; incompatibility then gives 42. Note lines
  109-110 and THM-4850 lines 31-32 should say so.
- **N5. The Zeckendorf section.**
  - "radix phi^2-ish" (line 190) should be "radix phi": F_k grows like phi^k.
  - {1,3,8} is the triple inside Fermat's quadruple {1,3,8,120}.
  - 1,3,8 = F_2,F_4,F_6 is the k = 1 case of Hoggatt-Bergum's Diophantine triples {F_2k, F_2k+2, F_2k+4}
    (Fibonacci Quart. 15 (1977); the quadruple conjecture was proved by Dujella, 1999). That is the one genuine
    Fibonacci law behind "1, 3, 8". It concerns products + 1, not the sums of Q_N.
  - Berlekamp-Guy, "Fibonacci Plays Billiards" (arXiv:2002.03705), is the classical place where Fibonacci-sum and
    square-sum chains meet.
  - The numerology typing is otherwise correct, and all the quoted decompositions check (`audit_analysis.out` J).
- **N6. "Adding a node is the odometer".** This is arithmetic on the indices. B(n) and B(n+1) are not related by
  inserting a vertex; for d < 72 they are images of the same input B(m). "So" in note line 51 implies that uniqueness
  follows from the replay. It does not: uniqueness is the residue argument, and the replay gives constructivity.
- **N7. Premature status.** Note line 12 says "Independently audited (section 8)" while §8 says "(pending)". Update
  after this audit.
- **N8. HYP-9447 pre-registration.**
  - Commit 07c2eac4ba (09:58:01) contains P1-P3, the command and both C files.
  - Both `.out` files first appear in 918129c522 (10:18:53).
  - The HYP diff adds only the status line and the results section.
  - The `.c` diff is MAXV 2048 -> 4096 only. That change is needed for N >= 2048 in the window run and has no effect
    on 63..1000.
  - Git cannot date the run itself, and the `.out` has no timestamps. The predictions concern certificates, so a
    found cycle proves existence whatever the timing.
- **N9. The definition is faithful.** The weak B definition (note lines 80-84, THM lines 16-20) matches session124
  README line 37 and `core.py` `check_pair`. The proof structure of THM-4850 (1) is complete:
  - exhaustive search for 42..53;
  - the base for 53..2620;
  - for n >= 2621, m = floor((n-24)/49) >= 53.

  The even member closes by 8 + 1 = 9, so THM-4851 (c) for even N >= 42 follows. Session124 README lines 11-19 cover
  both parities for N >= 2621 (P_1 begins 1,3).

## Recomputed numbers

| quantity | session | audit (own code) | status |
|---|---|---|---|
| Hamiltonian paths, N = 14..42 (reversal identified) | = A071983 (asserted to 34) | 0,1,1,1,0,0,0,0,0,3,0,10,12,35,52,19,20,349,392,669,4041,17175,12960,14026,11889,29123,39550,219968,553694 = OEIS 15..42 | agree |
| Hamiltonian cycles, N = 32..42 | 1,1,11,57,31,20,25,50,64,464,1062 | same; none below 32 | agree |
| weak B pairs | none n <= 41; all 42..53 | same, exhaustive 2..53 | agree |
| Shorter / Longer paths, n = 40..53 | note table | identical | agree |
| compatible Shorter, n = 42..53 | 10,4,45,3,5,11,859,108,153,230,10710,32242 | identical | agree |
| compatible Longer, n = 42..53 | not given | 10,7,75,2,11,13,1135,71,182,318,12840,43753 | new data |
| path pairs, n = 42..53 | 20,13,221,3,11,17,3384,276, -, -, -, - | 20,13,221,3,11,17,3384,276,1090,1351,106235,470586 | agree / filled in |
| odd-order 1,3 -> 10 (odd N <= 41) | only N = 41, unique | same path, unique; 136 at N = 43 | agree |
| even-order 1,3 -> 8, first order | 34 | 34 | agree |
| Q_32 / Q_33 cycles, edge difference | as stated | identical; 33's neighbours 3,16,31 pairwise non-adjacent; both 8-1-15 | agree |
| terminal cycles, N = 34..46 | existence (counts to 42) | 7,3,4,1,2,0,11,17,139,369,946,2699,3864 | agree |
| terminal existence, N <= 62 | {34..38} u [40,62] | same; complete searches at 30..33 and 39 by two own programs | agree |
| Q_39 | 50 cycles, none terminal | 50, none terminal (none through 8-1-3 either) | agree |
| terminal, N = 63..2620 | all found (cycles not stored) | all 2558 re-found and verified in Python | agree |
| insertion-cost minima, 32 -> 42 | 2,3,2,3,3,4,4,3,2,2 | same; all histograms identical | agree |
| operator replay | 539 pass | 539 (session script + session exe); 539 (session script + audit exe); 15,435 (audit, two checkers) | agree |
| numeration | unique for b in [42,2081] | unique to 400,000; also [41,2032] and [53,2620] | agree |
| Gerbicz-normalised nice pairs | not considered | exist iff n >= 41 (exhaustive 14..50) | new |
| free nice pairs | not considered | {15,16,25,26,27} u [29, oo) | new |

## Prior-art log

- **Repository.**
  - `session124_square_prefix_threshold42_20261031` (727bfdf464, 2026-10-31 00:44): weak B bases 42..52, recursion
    from 42, S_l = [42; 24...], cycles for N >= 2082. **MAJOR-1.**
  - `session124_mixed_radix_square_rays_20261031` (4338480c37, 01:21) consumes weak bases 42..2081.
  - `square_sum_spectrum_overnight_20261031b` and later packets consume threshold42.
  - `session124_gerbicz_reconstruction_20261031` (cb81e1d9d6, 01:31): Gerbicz's nice pairs, base 41..2032.
    **MAJOR-2.**
  - `square_sum_B_base_interval_20261019` lines 1 and 5: the term "nice pair" and its free definition.
  - `square_sum_spectrum_overnight_20261030` (b52068c6ae, 2026-10-30 23:29): the 32 -> 33 bridge step (N1).
  - `session32_square_monster_20261023/cycle-threshold/README.md` lines 18-22: Q_32 has a unique cycle.
  - `collatz_mod6_20260922_w6_square_sum_hamiltonicity.md` §3: counts to 36, cycles to 40, including a(33) = 1.
  - No earlier statement was found of: the non-existence of weak B pairs below 42, the gate path, the terminal
    spectrum or the hole at 39, or the insertion-cost table. Search terms: "8-1-3", "15-1-3", "hole at 39", "B41",
    "below 42", "smallest weak", plus `weak B`/`nice pair` over 05-knowledge/results and 01-canon.
- **OEIS (fetched 2026-10-31).**
  - A071983: paths, reversals identified; data 15..58 match.
  - A071984: loops up to rotation and reversal; a(32) = a(33) = 1; Dobbelaere's comment covers n <= 31.
  - A090461: Gerbicz's 2018-01-21 comment (paths for k >= 25, cycles for k >= 32).
  - A272259: lexicographically first loops; 1 sits between 8 and 15 at n = 32 and 33.
  - A search for the terminal counts 7,5,4,2,2,0,11,38,139 found no sequence.
- **Web.**
  - "ninja pair": no hits.
  - Mersenneforum "The Square-Sum problem" (showthread.php?p=477787), post #22, 2018-01-21: the definition of a
    nice pair, 49*n+res with res = 24..72, base n = 41..2032, and the 1 / 8 / 3 normalisation.
  - Berlekamp-Guy, arXiv:2002.03705: Fibonacci-sum and square-sum chains; links Gerbicz's code.
  - Hoggatt-Bergum 1977 and Dujella 1999: the Fibonacci Diophantine triples, of which 1,3,8 is the first.
  - No literature statement was found of the weak B threshold 42, the gate, the terminal spectrum or the hole at 39.

## How to rerun

These commands run from this folder (Git Bash). On Linux, drop `-Wl,--stack,...` and use `ulimit -s unlimited`.

```text
gcc -O2 -o ss_enum.exe ss_enum.c
gcc -O2 -Wl,--stack,268435456 -o ss_terminal_big.exe ss_terminal_big.c
./ss_enum.exe pairs 2 47 B  > pairs_B_2_47.out          # 1 s
./ss_enum.exe pairs 48 53 B > pairs_B_48_53.out         # 1 min
./ss_enum.exe pairs 14 45 G > pairs_G_14_45.out; ./ss_enum.exe pairs 46 50 G > pairs_G_46_50.out   # 30 s
./ss_enum.exe pairs 14 36 F > pairs_F_14_36.out; ./ss_enum.exe pairs 37 41 F > pairs_F_37_41.out   # 15 s
./ss_terminal_big.exe 63 2620 > SCRATCH/terminal_big_63_2620.txt                                  # 2 s, 14.9 MB, deterministic
python verify_terminal_big.py SCRATCH/terminal_big_63_2620.txt 63 2620 terminal_big_63_2620_receipt.json
python audit_analysis.py ./ss_enum.exe ./ss_terminal_big.exe > audit_analysis.out                   # 25 s, all assertions
python replay_own_pairs.py ./ss_enum.exe 40 > replay_own_pairs.out                                  # 80 s
python ../ninja_numeration_replay_20261031.py ./ss_enum.exe > replay_session_script_with_own_enum.out
gcc -O2 -o session_enum_for_replay.exe ../ninja_pairs_square_sum_20261031.c   # only as replay input
python ../ninja_numeration_replay_20261031.py ./session_enum_for_replay.exe > replay_session_script_with_session_enum.out  # 5 min
rm -f *.exe
```

The executables are not committed. The 14.9 MB terminal certificate is kept out of the repository and is regenerated
byte-for-byte by the third command; its SHA-256 is in the receipt.
