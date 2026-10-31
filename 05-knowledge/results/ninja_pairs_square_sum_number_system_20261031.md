# Weak B pairs start sharply at 42 (Gerbicz's nice pairs at 41): thresholds of the square-sum pair numerations, the prehistory before them, and the base-49 numeration

opus-2026-10-31-S69, 2026-10-31.

**Status.**
- **PROVED:** THM-4850 and THM-4851, by exhaustive computation for small orders plus cited certified constructions.
- **FINITE-EXACT:** all counts.
- **CITED:** Gerbicz 2018; OEIS A071983, A071984; session124 and threshold42; the 2026-10-30 square bridge certificate.
- **HYP-9447:** pre-registered, P1-P3 CONFIRMED.
- **INTERPRETATION** for the owner's terminology; **NUMEROLOGY / weak observation** for the 22/33, 24-25 and
  Zeckendorf remarks.
- Independently audited, **SOUND WITH CORRECTIONS**, all applied (section 8; MISTAKE-694).

## 0. Directive and how the owner's terms were read

Owner, 2026-10-31: "start a new session where you consider a number system that starts at the 32 node square-sum
graph and considers that the start of the natural numbers, this means that all natural numbers can be reached by
adding nodes, but hidden before the start of the numbers is some interesting structure which dictates their future
behavior. note how the Shorter and Longer of a ninja pair correspond with consecutive integers, one must be odd and the
other even, and study 1 vs 3 and 8 as a triple, feel how this relates with zeckendorf decomposition related ideas.
think {15,16,17,23,25} 32 and 35. also 3 with 33 and 22, ninja pairs and 49,24,72".

The square-sum graph Q_N has vertices 1..N, with u ~ v when u + v is a perfect square.

**The construction behind all large square-sum cycles is Gerbicz's** (mersenneforum "The Square-Sum problem", post
#22, 2018-01-21).
- A **nice pair** is a Hamiltonian path on [1, n] and one on [1, n+1] whose common vertices sit at positions of equal
  parity.
- From a base table n = 41..2032, every 49n + res, res = 24..72, is reached.

The repository's 2026-10-19..31 lanes rebuilt and extended it:
- `square_sum_B_base_interval_20261019` ("square-sum nice pairs", base 53..2620);
- `session124_gerbicz_reconstruction_20261031`;
- `session124_recursive_square_prefix_20261031`: the **weak B** normalisation (start 1, 3; even-order member ends at
  8, odd-order member at 10) and 98 paired operators;
- `session124_square_prefix_threshold42_20261031`: weak B bases 42..52, base block [42, 2081], and terminal cycles for
  N >= 2082;
- `session124_mixed_radix_square_rays_20261031`;
- the `square_sum_spectrum_overnight_*` packets.

**Reading** (INTERPRETATION).
- The owner's "ninja pair" is most plausibly Gerbicz's "nice pair". The term "ninja pair" appears neither in the
  repository nor on the web.
- Shorter and Longer are its two members, of consecutive orders n and n+1, one odd and one even.
- "49, 24, 72" is the radix and residue range.
- "1 vs 3 and 8" is the closing triple 8-1-3 of the even member (1+3 = 4, 1+8 = 9).

Three normalisations are in use, and this note studies the last:
- free nice pairs (threshold 15);
- Gerbicz-normalised pairs (threshold 41);
- weak B pairs (threshold 42).

## 1. Results first

1. **Thresholds of the pair numerations** (THM-4850; exhaustive here and in the audit).
   - Weak B pairs exist exactly for n >= 42. Sharpness is new: none for n <= 41. Existence from 42 on is prior work
     (threshold42).
   - Gerbicz-normalised nice pairs exist exactly for n >= 41.
   - Free nice pairs exist exactly for n in {15, 16, 25, 26, 27} u [29, oo).
   - **The weak B gate** is the odd member ending at 10. Among odd orders <= 41 it exists only at order 41, as a single
     path:
     1,3,6,30,19,17,32,4,21,15,34,2,7,18,31,33,16,9,40,41,8,28,36,13,23,26,38,11,5,20,29,35,14,22,27,37,12,24,25,39,10.
     It is parity-incompatible with every partner at n = 40 and 41.
2. **The base-49 numeration.** Every n >= 42 is uniquely n = 49^k b + sum d_i 49^i with b in [42, 2081] and
   d_i in [24, 72] (the owner's 49, 24, 72).
   - Session124's 98 paired operators make this constructive. They were replayed here on new pairs at m = 42..52:
     539 outputs, all passing session124's own validator, which re-confirms threshold42.
   - In Gerbicz's normalisation the same numeration starts at 41, with base block [41, 2032].
3. **The prehistory** (THM-4851).
   - Paths exist for {1, 15, 16, 17, 23} u [25, oo). These are the owner's {15,16,17,23,25} (CITED).
   - Cycles exist from 32 (CITED). Q_32 and Q_33 have unique cycles, one segment reversal plus an insertion apart (as
     in the 10-30 bridge certificate). Both close 8-1-15.
   - **Terminal cycles** (8-1-3 even / 15-1-3 odd) exist **exactly for N in {34..38} u [40, oo)**. The owner's 35 is
     the first odd one, and 39 is a hole: Q_39 has 50 cycles, none terminal.
   - Adding a node from 32 to 42 always deletes at least two old edges.
4. **3, 22, 33 and 24-25** (weak observation; section 4). Against all paths with the same prefix and end, the members
   of the first weak B pairs leave 3 through 22 or 33 more often at some n and less often at others. The 24-25 edge
   (sum 49) is common in all such paths. Neither is a law.
5. **Zeckendorf** (section 5): a structural comparison, a classical Fibonacci fact behind 1, 3, 8, and a NUMEROLOGY
   table.

## 2. THM-4850 in detail

The exhaustive search (`ninja_pairs_square_sum_20261031.c`, mode `pairs`) enumerates every Longer, collects its parity
masks on [1, n], and then enumerates every Shorter.

| n | Shorter paths (order n) | Longer paths (order n+1) | Shorter paths with a compatible Longer | path pairs |
|---|---|---|---|---|
| <= 32 | 0 | 0 | 0 | 0 |
| 33..39 | one of the two members has no path at all | | 0 | 0 |
| 40 | 11 (end 8) | 1 (end 10: the gate) | 0 | 0 |
| 41 | 1 (the gate) | 139 (end 8) | 0 | 0 |
| 42 | 139 | 136 | 10 | 20 |
| 43 | 136 | 946 | 4 | 13 |
| 44 | 946 | 1129 | 45 | 221 |
| 45 | 1129 | 3864 | 3 | 3 |
| 46 | 3864 | 6083 | 5 | 11 |
| 47 | 6083 | 8212 | 11 | 17 |
| 48 | 8212 | 19615 | 859 | 3384 |
| 49 | 19615 | 8591 | 108 | 276 |
| 50 | 8591 | 24777 | 153 | 1090 |
| 51 | 24777 | 153594 | 230 | 1351 |
| 52 | 153594 | 302246 | 10710 | 106235 |
| 53 | 302246 | 481425 | 32242 | 470586 |

The threshold has two causes:
- odd-member scarcity forces n >= 40;
- parity incompatibility at 40 and 41 forces 42.

The input pairs at 42..52 also pass session124's own `check_pair`.

**Why pairs of consecutive orders.**
- The operator maps the vertex v in position j to 49v + (-1)^j c, with c in [-24, 24]. A square sum a + b = s^2 becomes
  49 s^2 = (7s)^2, so the radix must be a perfect square.
- An output of order 49m + e takes exactly e - 24 of its 49 offset blocks from the Longer (those c with
  49 + (-1)^(m+1) c <= e) and the rest from the Shorter.
- Parity compatibility makes the two members interchangeable on [1, m]. So one input pair serves all 49 residues.

**The numeration.**
- Uniqueness of n = 49^k b + sum d_i 49^i (b in [42, 2081], d_i in [24, 72]) is the residue argument. The replay
  supplies constructivity.
- Examples: 2082 = [42; 24]; 2130 = [42; 72]; 2131 = [43; 24]; 128453 = [53; 24, 24]. Threshold42's
  S_l = (85 * 49^l - 1)/2 are [42; 24, ..., 24]; session124's T_l = (107 * 49^l - 1)/2 are [53; 24, ..., 24].
- The odometer n -> n+1 is index arithmetic. For d < 72, B(n) and B(n+1) are images of the same input B(m); they are
  not related by inserting a vertex.

## 3. THM-4851 in detail: the prehistory

- **(a)** The path counts equal OEIS A071983 for 15 <= N <= 42 and the cycle counts equal A071984 for 32 <= N <= 42:
  1, 1, 11, 57, 31, 20, 25, 50, 64, 464, 1062 (audit; the driver asserts A071983 to 34). Existence for all N >= 25
  (paths) and N >= 32 (cycles) is Gerbicz's.
- **(b)** Q_32's unique cycle is 1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,29,20,16,9,27,22,14,2,23,26,10,15.
  - Q_33's is the same with 31,18,7,29,20 reversed and 33 inserted: delete 5-31 and 20-16; add 5-20, 31-33 and 33-16.
  - The step is no direct insertion, because 33's neighbours 3, 16 and 31 are pairwise non-adjacent in the 32-cycle.
  - The 2026-10-30 bridge certificate records it as {'insert': 15, 'reverse': [[10, 16]]}.
  - Both cycles contain 8-1-15, both of session124's closing endpoints. So neither is terminal.
- **(c) Terminal spectrum.** Q_N has a cycle through 8-1-3 (N even) or 15-1-3 (N odd) iff
  N in {34, 35, 36, 37, 38} u [40, oo).
  - N <= 62: exhaustive; complete searches give the absences.
  - 63..2620: every N found by a search whose built-in checker verified each cycle. 63..1000 was pre-registered as
    HYP-9447. The audit re-found all 2558 with a different program and re-verified them in Python.
  - N >= 2082: the threshold42 terminal theorem.
  - Every even N >= 42 also follows from THM-4850.
- **(d) Insertion costs** 32 -> 42: the minima are 2, 3, 2, 3, 3, 4, 4, 3, 2, 2. No step is a direct insertion.

**The ladder of starts.**

| threshold | event |
|---|---|
| 15 | first Hamiltonian path; free nice pairs begin |
| 25 | paths for all larger N |
| 32 | first Hamiltonian cycle (unique, as is the 33-cycle) |
| 34 / 35 | first terminal cycles, even and odd |
| 39 | the only later hole |
| 40 | terminal cycles for all larger N |
| 41 | Gerbicz-normalised nice pairs begin; the weak-B gate path (unique) |
| 42 | weak B pairs for all larger n |

The owner's "start of the natural numbers" at 32 is the cycle threshold. The pair numerations start at 41 and 42.
Between them lie the two unique cycles, the terminal onset at 34/35, the hole at 39 and the gate at 41.

## 4. The owner's remaining numbers (weak observations)

- **1, 3, 8.** The even weak-B member is a path 1,3,...,8, closed by 8 + 1 = 9 into a terminal cycle through 8-1-3.
  Odd terminal cycles close 15-1-3. The vertices 3, 8, 15 are 1's neighbours k^2 - 1.
- **3 with 22 and 33; the edge 24-25.** The audit compared the members of all weak B pairs with all paths of the same
  prefix and end (`audit_analysis.out` G):

| n | members | 22/33 after 1,3 (members) | 22/33 baseline | 24-25 (members) | 24-25 baseline |
|---|---|---|---|---|---|
| 42 | 20 | 80% | 65% | 85% | 79% |
| 43 | 11 | 100% | 62% | 100% | 72% |
| 44 | 120 | 57% | 56% | 82% | 69% |
| 45 | 5 | 80% | 39% | 100% | 77% |
| 46 | 16 | 6% | 37% | 94% | 72% |
| 47 | 24 | 50% | 32% | 21% | 73% |
| 48 | 1994 | 25% | 29% | 81% | 73% |
| 49 | 179 | 59% | 38% | 79% | 72% |

  - 22 and 33 (sums 25 and 36) are already the most common third vertices among all such paths at 42-43. The members
    exceed the baseline at some n and fall below it at others.
  - The 24-25 edge is common because 25 has only the neighbours 11, 24 and 39 there.
  - Weak observation only.
- **49, 24, 72.** The radix and residue range (section 2).

## 5. Zeckendorf

- **Structural comparison** (INTERPRETATION).
  - Zeckendorf's theorem writes every integer uniquely as a sum of non-consecutive Fibonacci numbers. That is a
    positional system of radix phi (F_k grows like phi^k).
  - The square-sum pair numeration is positional with radix 49 and a complete digit set [24, 72].
  - The square-sum radix must be a perfect square (scaling keeps sums square), so a Pisot radix such as phi cannot
    occur for this mechanism.
  - Session124's mixed-radix rays (radices 25 and 49) are the closest square-sum analogue of a mixed numeration.
- **The one Fibonacci law behind "1, 3, 8"** (CITED).
  - 1, 3, 8 = F_2, F_4, F_6 is the k = 1 case of Hoggatt-Bergum's Diophantine triples {F_2k, F_2k+2, F_2k+4}, whose
    pairwise products + 1 are squares (Fibonacci Quart. 15 (1977)). {1,3,8} lies inside Fermat's quadruple
    {1,3,8,120}, and Dujella (1999) proved the Hoggatt-Bergum quadruple extension.
  - This concerns products, not the sums of Q_N.
  - Berlekamp-Guy, "Fibonacci plays billiards" (arXiv:2002.03705), is where Fibonacci-sum and square-sum chains meet
    classically.
- **NUMEROLOGY:**
  - 32 = 21+8+3 and 33 = 21+8+3+1 = F_9 - 1; 34 = F_9; 35 = 34+1; 39 = 34+5; 41 = 34+5+2; 42 = 34+8;
  - 15 = 13+2, 16 = 13+3, 17 = 13+3+1, 23 = 21+2, 25 = 21+3+1.

## 6. Prior art and novelty

- **CITED and reproduced:**
  - Gerbicz's nice pairs and numeration;
  - OEIS A071983, A071984;
  - session124 (weak B pairs, operators, base 53..2620);
  - threshold42 (bases 42..52, base block [42, 2081], terminal cycles for N >= 2082; the replay here re-confirms it);
  - the 2026-10-30 bridge certificate (the 32 -> 33 step).
- **New to the repository:**
  - the sharpness of 42 (no weak B pair for n <= 41);
  - the uniqueness of the gate path;
  - the exact counts for n <= 53;
  - the normalisation table (the free and Gerbicz thresholds, from the audit);
  - the exact terminal spectrum with the onset at 34/35 and the hole at 39;
  - the insertion-cost table.

## 7. Reproduction

```text
gcc -O2 -o enum 04-computation/experiments/ninja_pairs_square_sum_20261031.c
gcc -O2 -o term 04-computation/experiments/terminal_cycles_large_20261031.c
python3 04-computation/experiments/ninja_square_sum_number_system_20261031.py ./enum   # about 11 min; all assertions
./term 63 1000 1 1000000000     # HYP-9447 (pre-registered): terminal_cycles_large_20261031.out
./term 1001 2620 1 1000000000   # terminal_cycles_large_20261031_window.out
python3 04-computation/experiments/ninja_numeration_replay_20261031.py ./enum          # ninja_numeration_replay_20261031.out
```

The audit's independent programs and receipts are in `04-computation/experiments/ninja_square_sum_20261031_audit_S69/`.

## 8. Audit record

Independent audit `ninja_square_sum_20261031_audit_S69/AUDIT.md`: **SOUND WITH CORRECTIONS**, all applied.

- **What reproduced.** Every count, the gate, the unique cycles, the terminal spectrum (with the 2558 cycles in
  63..2620 re-found by a different program and re-verified in Python), the insertion costs, and the replay (539, plus
  15,435 outputs with the audit's own pairs).
- **MAJOR.**
  - (1) Threshold-42 existence, the base block and the replay were prior work: session124_square_prefix_threshold42,
    the same morning. They are now credited, and only the sharpness, the gate and the counts are new.
  - (2) "Ninja pair" is most plausibly Gerbicz's "nice pair"; the threshold depends on the normalisation (15 / 41 / 42).
- **MINOR.**
  - Selective 22/33 and 24-25 statistics (the baseline table above).
  - (-1)^(n+1) corrected to (-1)^(m+1), and the operator description corrected.
  - "Verified independently" overstated.
  - The replay output was missing (now recorded).
  - A stale HYP reference in the search program.
- **NOTES.**
  - The 10-30 bridge certificate.
  - The OEIS range.
  - The table rows.
  - Two causes of the threshold.
  - Zeckendorf: radix phi; Hoggatt-Bergum; Berlekamp-Guy.
  - The odometer is index arithmetic.
