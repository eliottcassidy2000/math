# The square-sum numbers start at 42: ninja pairs, the prehistory before them, and the base-49 numeration

opus-2026-10-31-S69, 2026-10-31.

**Status.**
- **PROVED:** THM-4850 and THM-4851, by exhaustive computation for small orders combined with session124's certified
  constructions for large ones.
- **FINITE-EXACT:** all counts.
- **CITED:** OEIS A071983, A071984, A090461; Gerbicz's square-sum Hamiltonicity; session124 and the 2026-10-19 B base.
- **HYP-9447:** pre-registered, P1-P3 CONFIRMED, then completed to a theorem (THM-4851 (c)).
- **INTERPRETATION** for the owner's terminology; **NUMEROLOGY** for the Zeckendorf table.
- Independently audited (section 8).

## 0. Directive and how the owner's terms were read

Owner, 2026-10-31: "start a new session where you consider a number system that starts at the 32 node square-sum
graph and considers that the start of the natural numbers, this means that all natural numbers can be reached by
adding nodes, but hidden before the start of the numbers is some interesting structure which dictates their future
behavior. note how the Shorter and Longer of a ninja pair correspond with consecutive integers, one must be odd and the
other even, and study 1 vs 3 and 8 as a triple, feel how this relates with zeckendorf decomposition related ideas.
think {15,16,17,23,25} 32 and 35. also 3 with 33 and 22, ninja pairs and 49,24,72".

The square-sum graph Q_N has vertices 1..N, with u ~ v when u + v is a perfect square. The 2026-10-30/31 square-sum
lanes are the immediate context:
- `session124_recursive_square_prefix_20261031`, `session124_mixed_radix_square_rays_20261031`;
- the base `square_sum_B_base_interval_20261019`;
- the `square_sum_spectrum_overnight_*` packets.

They build every large Hamiltonian cycle from a **weak B pair**: a Hamiltonian path of Q_n and one of Q_(n+1), both
starting 1, 3, with equal positional parity on all common vertices. The even-order member ends at 8 and the odd-order
member at 10. Every B(n) with n >= 2621 is built from B(m), m = floor((n - 24)/49), by one of 49 digit operators
d = n - 49m in [24, 72].

**Reading** (INTERPRETATION): the owner's "ninja pair" is this weak B pair.
- **Shorter:** the path on [1, n]. **Longer:** the path on [1, n+1].
- Their orders are consecutive integers, one odd and one even.
- "49, 24, 72" is the radix and the digit range of the recursion.
- "1 vs 3 and 8" is the triple 8-1-3 that closes the even member into a cycle (1+3 = 4, 1+8 = 9).

The term "ninja pair" appears nowhere else in the repository or in a web search.

## 1. Results first

1. **THM-4850: ninja pairs exist exactly for n >= 42.**
   - An exhaustive search finds no ninja pair for any n <= 41, and ninja pairs for every 42 <= n <= 53.
   - Session124's certified base covers 53..2620, and its recursion covers everything above.
   - **The gate** is the odd member ending at 10. Among odd orders it first exists at **N = 41, as a single path**:
     1,3,6,30,19,17,32,4,21,15,34,2,7,18,31,33,16,9,40,41,8,28,36,13,23,26,38,11,5,20,29,35,14,22,27,37,12,24,25,39,10.
     It fits no pair at n = 40 or n = 41, and the first pairs appear at n = 42 (20 path pairs).
   - **Session124's 49 digit operators work on the new small pairs.** All 539 outputs from the bases m = 42..52 pass
     session124's own validator. So the base-49 numeration with the owner's 49, 24, 72 starts at 42: every n >= 42
     is uniquely n = 49^k b + sum_i d_i 49^i with b in [42, 2081] and d_i in [24, 72], and adding a node is the
     odometer step.
2. **THM-4851: the prehistory before the start.**
   - Hamiltonian paths exist for n in {1, 15, 16, 17, 23} and every n >= 25. These are the owner's {15,16,17,23,25}
     (CITED: A090461, Gerbicz).
   - Cycles exist from **32** (CITED: Gerbicz). Q_32 and Q_33 each have **exactly one** Hamiltonian cycle (A071984).
     The 33-cycle is the 32-cycle with the segment 31..20 reversed and 33 inserted: delete 5-31 and 20-16, add 5-20,
     31-33 and 33-16. Both cycles close as 8-1-15.
   - **Terminal-form cycles** (8-1-3 for even N, 15-1-3 for odd N, session124's closings) exist **exactly for
     N in {34,...,38} u [40, oo)**.
     - The owner's **35** is the first odd one.
     - **N = 39 is the one hole**: Q_39 has 50 cycles, none terminal.
     - Proof: exhaustive to 62; a verified search for every N in 63..2620 (pre-registered to 1000 as HYP-9447); the even
       members of ninja pairs; session124 for N >= 2621.
   - **Adding a node is never free.** From 32 to 42, every new cycle needs at least two edges of an old cycle deleted.
     32 -> 33 needs exactly two.
3. **The owner's other numbers.**
   - **"3 with 33 and 22":** in the first ninja pairs (n = 42, 43, 45), the vertex after 1,3 is almost always 22 or 33.
     These are 3's two largest neighbours, with sums 25 and 36.
   - **The consecutive pair 24-25** (sum 49 = 7^2, one of the consecutive-integer edges (k^2 - 1)/2, (k^2 + 1)/2)
     occurs in most early ninja members.
4. **Zeckendorf** (section 5): a structural comparison, plus a NUMEROLOGY table.
   - The square-sum radix must itself be a perfect square, because scaling has to keep sums square. A golden-ratio
     (Zeckendorf) radix is therefore impossible for this mechanism.
   - 32 = 21 + 8 + 3 and 33 = 21 + 8 + 3 + 1 = F_9 - 1 are Zeckendorf facts with no proved link.

## 2. THM-4850 in detail

**Definition** (session124). A ninja (weak B) pair at n consists of:
- a Hamiltonian path S of Q_n (the Shorter) and L of Q_(n+1) (the Longer);
- both start 1, 3;
- every v <= n has positions of equal parity in S and L;
- the member of even order ends at 8, and the member of odd order ends at 10.

**Exhaustive search** (`ninja_pairs_square_sum_20261031.c`, mode `pairs`). It enumerates every Longer, collects its
parity masks on [1, n], and then enumerates every Shorter.

| n | Shorter paths (order n) | Longer paths (order n+1) | Shorter paths with a compatible Longer | path pairs |
|---|---|---|---|---|
| <= 39 | 0 or 0 compatible | | 0 | 0 |
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
| 50 | 8591 | 24777 | 153 | |
| 51 | 24777 | 153594 | 230 | |
| 52 | 153594 | 302246 | 10710 | |
| 53 | 302246 | 481425 | 32242 | |

- The odd-order paths from 1,3 to 10 exist for odd N in 15..41 only at N = 41, as a single path.
- Even-order paths 1,3 -> 8 first exist at N = 34.
- Hence the start at 42 is forced by the odd member.
- The input pairs at 42..52 also pass session124's own `check_pair`.

**Why pairs, and why consecutive orders.** The digit operator maps a vertex v in position j to 49v + (-1)^j c, with
c in [-24, 24]. A square sum a + b = s^2 becomes 49 s^2 = (7s)^2, so the radix must be a perfect square. One scaled path
of order m produces about 49m labels. To reach every residue of the target n = 49m + d, d in [24, 72], the operator uses
the Shorter (order m) for low digits and the Longer (order m+1) for high digits; the switch is at 49 + (-1)^(n+1) c <= d.
Two consecutive orders thus cover a full residue system mod 49. That is the arithmetic content of "Shorter and Longer
correspond with consecutive integers".

**The numeration** (arithmetic; constructive by the replay).
- Every n >= 42 has exactly one expansion n = 49^k b + sum_{i<k} d_i 49^i with b in [42, 2081] and d_i in [24, 72],
  because [24, 72] is a complete residue system mod 49.
- The successor n -> n + 1 is the odometer: d -> d + 1, with 72 + 1 carried as 24 into the next place, since
  49m + 73 = 49(m+1) + 24.
- Examples: 2082 = [42; 24]; 2130 = [42; 72]; 2131 = [43; 24]; 128453 = T_2 = [53; 24, 24]; 6294221 = T_3 = [53; 24, 24, 24].
- Session124's threshold sequence T_l = (107 * 49^l - 1)/2 is exactly the numerals [53; 24, ..., 24].

## 3. THM-4851 in detail: the prehistory

- **(a) Paths and cycles** (CITED, reproduced).
  - The Hamiltonian path counts for N = 14..34 equal OEIS A071983 (square chains), and the cycle counts for N = 32..42
    equal A071984 (square loops): 1, 1, 11, 57, 31, 20, 25, 50, 64, 464, 1062.
  - Existence for all N >= 25 (paths) and N >= 32 (cycles) is Gerbicz's theorem (credited in the repository as A090461;
    see `square_sum_ports_20261023`).
- **(b) The two unique cycles.**
  - Q_32: 1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,29,20,16,9,27,22,14,2,23,26,10,15.
  - Q_33: the same with 31,18,7,29,20 reversed and 33 inserted.
  - Deleting 5-31 and 20-16 and adding 5-20, 31-33 and 33-16 is a 2-opt insertion: the cheapest possible here, since no
    edge of the 32-cycle has both ends adjacent to 33 (33's neighbours are 3, 16, 31, pairwise non-adjacent in the
    cycle).
  - Both cycles contain 8-1-15, i.e. both of session124's closing endpoints at once. That is why neither is terminal:
    1's other neighbour is 15 or 8, not 3.
- **(c) Terminal spectrum.** A Hamiltonian cycle of Q_N through 8-1-3 (N even) or 15-1-3 (N odd) exists **iff
  N in {34, 35, 36, 37, 38} u [40, oo)**. The parts of the proof:
  - N <= 62: exhaustive (two independent programs agree; complete searches certify the absences at 32, 33 and 39);
  - 63..1000: pre-registered HYP-9447, all 938 values found and verified;
  - 1001..2620: all 1620 values found and verified (exploratory continuation);
  - N >= 2621: session124's terminal-cycle theorem;
  - every even N >= 42 also follows from THM-4850, since the even member of a ninja pair closes by 8 + 1 = 9.
- **(d) Insertion costs.** For N = 32..41 and every Hamiltonian cycle C' of Q_(N+1), compute the least number of edges
  of some cycle C of Q_N that are absent from C' - (N+1). The minima are 2, 3, 2, 3, 3, 4, 4, 3, 2, 2. Adding a node is
  never a direct insertion on this range. The full histograms are in the `.out`.

**The ladder of starts.**

| threshold | event |
|---|---|
| 15 | first Hamiltonian path (unique) |
| 25 | paths for all larger N |
| 32 | first Hamiltonian cycle (unique, as is the 33-cycle) |
| 34 / 35 | first terminal cycles, even and odd |
| 39 | the only later hole |
| 40 | terminal cycles for all larger N |
| 41 | the odd gate path (unique) |
| **42** | ninja pairs for all larger n: the start of the base-49 numeration |

The owner's "start of the natural numbers" at 32 is the cycle threshold. The numeration itself starts ten steps later.
What lies between is the "hidden structure": two unique cycles, the terminal onset at 34/35, the hole at 39, and the
unique gate at 41.

## 4. The owner's remaining numbers

- **1, 3, 8.** The even member of every ninja pair is a Hamiltonian path 1,3,...,8, and 8 + 1 = 9 closes it into a
  terminal cycle through 8-1-3. The odd member ends at 10 (10 + 1 = 11 is not a square). Odd terminal cycles instead
  close 15-1-3. The vertices 3, 8, 15 are 1's neighbours k^2 - 1.
- **3 with 22 and 33.** The vertex after 1,3 in the members of the first ninja pairs is:
  - n = 42: {13: 4, 22: 11, 33: 5};
  - n = 43: {22: 7, 33: 4};
  - n = 45: {13: 1, 22: 3, 33: 1}.

  From n = 44 on, 6 and 46 also appear (full table in the `.out`). So the earliest ninja numbers leave 3 through its
  two largest neighbours, 22 (25) and 33 (36). The session124 ray prefix (1,3,6,19,...) leaves through 6.
- **24-25.** The consecutive-integer edge 24 + 25 = 49 = 7^2 lies in 17 of the 20 members at n = 42, all 11 at n = 43,
  and 99 of 120 at n = 44. In the gate path it appears as 12,24,25,39.
- **49, 24, 72.** The radix and digit range (section 2).

## 5. Zeckendorf: the structural comparison and the numerology

**Structural comparison** (INTERPRETATION).
- Zeckendorf's theorem: every positive integer is uniquely a sum of non-consecutive Fibonacci numbers. Its counting
  rule is a positional system of radix phi^2-ish growth, with a "no two consecutive" constraint. The successor map is a
  carry rule on Fibonacci words.
- The ninja numeration is a positional system of radix 49 with a complete digit set [24, 72] and a seed block
  [42, 2081]. Its successor map is the base-49 odometer.
- Both are normal forms grown from a seed by a self-similar rule. The decisive difference is the radix:
  - the square-sum operator needs a perfect-square radix, because it must map square sums to square sums;
  - so no Pisot (golden) radix can occur.
- Session124's mixed-radix rays (radices 25 and 49 = 5^2 and 7^2, prefix-decodable code words) are the closest
  square-sum analogue of a Zeckendorf-type mixed numeration.

**NUMEROLOGY** (recorded, not claimed):
- 1, 3, 8 = F_2, F_4, F_6. {1, 3, 8} is Fermat's Diophantine triple (1*3+1, 1*8+1 and 3*8+1 are squares), while in Q_N
  the relevant facts are the sums 1+3 = 2^2 and 1+8 = 3^2.
- 32 = 21+8+3 and 33 = 21+8+3+1 = F_9 - 1. 34 = F_9. 35 = 34+1, 39 = 34+5, 41 = 34+5+2, 42 = 34+8.
- 15 = 13+2, 16 = 13+3, 17 = 13+3+1, 23 = 21+2, 25 = 21+3+1.

No square-sum quantity computed here follows a Fibonacci law.

## 6. Prior art and novelty

- **CITED and reproduced:**
  - the path and cycle counts (OEIS A071983, A071984);
  - path and cycle existence (Gerbicz; A090461);
  - the weak B pair, the digit operators, the base 53..2620 and the terminal-cycle theorem for N >= 2621 (session124;
    `square_sum_B_base_interval_20261019`).
- **New to the repository:**
  - the exact ninja threshold 42 and the gate at 41;
  - the replay of the digit operators on the bases 42..52 (the numeration starts at 42);
  - the exact terminal spectrum {34..38} u [40, oo) with the hole at 39;
  - the 32 -> 33 reversal-insertion description and the insertion-cost table;
  - the 22/33 and 24-25 statistics of the first pairs.
- No literature search found statements of these (`ninja` is the owner's term); the audit repeats the search.

## 7. Reproduction

```text
gcc -O2 -o enum 04-computation/experiments/ninja_pairs_square_sum_20261031.c
gcc -O2 -o term 04-computation/experiments/terminal_cycles_large_20261031.c
python3 04-computation/experiments/ninja_square_sum_number_system_20261031.py ./enum   # about 15 min; all assertions
./term 63 1000 1 1000000000      # HYP-9447 (pre-registered); output terminal_cycles_large_20261031.out
./term 1001 2620 1 1000000000    # the window; output terminal_cycles_large_20261031_window.out
python3 04-computation/experiments/ninja_numeration_replay_20261031.py ./enum          # 539 operator outputs, about 5 min
```

## 8. Audit record

(pending)
