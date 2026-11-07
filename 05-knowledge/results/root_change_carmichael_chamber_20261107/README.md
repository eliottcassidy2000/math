# The Carmichael chamber of root-change words; two root changes to 10^9; the C = 1 spectrum below 10^9

kind-pasteur-2026-11-07-S14. **THM-4904** (PROVED), **HYP-9494** (CONJECTURE with FINITE-EXACT evidence), and the
pre-registered census extending HYP-9480 (D). Chan-Pak completeness and Carmichael's totient conjecture remain OPEN.

Owner prompt: merge in Smith numbers and repunit primes in all bases, investigate prior work on the Carmichael
conjecture and how 14, 26, 34, 38, ... scale against the repo's sequences, while working on the two-root-change
family.

**Concurrent same-prompt work, credited and not duplicated.** Session218 finished the Smith, repunit and nontotient
parts of that prompt one hour before this packet:
- the all-base Smith classification;
- the binary repunit multiplier obstruction;
- the Mersenne-Wagstaff to nontotient bridge;
- 14, 26, 34, 38 as even nontotients versus the Camion exponents, which diverge at 50;
- dyadic repair of nontotients.

See `../session218_smith_repunit_totient_20261107/`. This packet adds the Carmichael literature map, the word-level
Carmichael theorem and the two-root-change work.

## Results

1. **[CARMICHAEL_MAP.md](CARMICHAEL_MAP.md): Carmichael's totient conjecture in depth.**
   - The literature, read in Ford's own 1998 text: the bounds 10^37, 10^400, 10^10,000, 10^10,000,000 and
     10^(10^10); the Carmichael-Klee lemma; Pomerance's criterion; Ford's dichotomy; Sierpinski's conjecture,
     settled by Ford and Konyagin and by Ford; the 2025 special cases.
   - Every repository Carmichael-fiber result, with its fiber notion and exceptional set.
   - A mechanism dictionary between totients and words.
   - A short scaling check of 14, 26, 34, 38 against the root-change miss lists, which finds no fit.
2. **THM-4904 ([PROOF.md](PROOF.md), sections 1-4), PROVED.** For the word multiplicity M_K(N), the number of
   O/Z/J words with at most K root changes and value N:
   - **Twin decomposition:** M_K = [N=1] + 2 A_(K-1) + B_K, the analogue of phi(x) = phi(2x).
   - **Parity law:** M_K(N) = [N=1] + Fix_K(N) (mod 2), by Session199's adjoint.
   - **Chamber:** a singleton word is dagger-fixed, uses exactly K root changes, its value lies outside F_<=K-1, and
     it neither begins nor ends with J.
   - **Growth lemma:** U >= 2 c0^(#J+1) 2^(#O/Z), with c0 = 0.44257.
   - **Density zero:** singleton values number O_K(sqrt(X) (log X)^floor(K/2)), sharp for K = 0. This is the opposite
     of Ford's "none or a positive proportion" for phi.
3. **Two root changes to 10^9 (pre-registered, PROOF section 5).**
   - P1 CONFIRMED: exactly the 140 misses, the largest 3,780,967.
   - P3 CONFIRMED: no word-singleton above 1. This is Carmichael's conjecture for two root changes, through 10^9.
   - P2 and P4 FAILED: both predicted extremes of a heavy-tailed statistic.
4. **Root-change depth and the uniquely Hamiltonian spectrum (PROOF section 5.1), FINITE-EXACT.**
   - Every odd N in (2647, 10^9] is H of a strong tournament with exactly one Hamiltonian cycle, built with at most
     seven root changes.
   - Below 10^9, the C = 1 spectrum misses exactly Session197's 29 values. This extends the previous bound of 50,001.
   - Depth 7 is attained only by 1877 and 1897. Their literal 14-vertex witnesses are checked by subset DP.
5. **Exploratory, HYP-9494 (A)-(E).**
   - The singleton hierarchy: about 0.49 sqrt(X) for K = 0, sporadic for K = 1, none for K = 2.
   - A 3-adic local factor of the two-root-change intensity: mean M_2 = 206, 332 and 457 for v3(N) = 0, 1 and 2.
   - Near-misses are never divisible by 3.

## Files

| File | Role |
| --- | --- |
| `PREREGISTRATION.md` | P1-P4, pushed in 6d0e88f0e8 before the 10^9 output was read |
| `rc_census.c` | the census, byte-identical to the pre-registered version (OpenMP; coverage bitmaps per level, 8-bit saturating multiplicity) |
| `rc_1e9.txt`, `PREREG_OUTCOMES.txt`, `evaluate_prereg.py` | raw 10^9 output and its evaluation |
| `rc_depth.c`, `depth_K7_1e6.txt`, `depth_K3_4e6.txt` | up to seven root changes, full miss lists |
| `verify_chamber.c`, `VERIFY_CHAMBER_2p20.txt` | independent word-by-word checks of the adjoint (7.17 M words), (phys) and the growth lemma (76.3 M states), and THM-4904 (1)-(3) for K = 0, 1, 2 to 2^20; its window histograms equal `rc_census`'s exactly |
| `literal_witnesses.py`, `LITERAL_WITNESSES.txt` | exhaustive minimal-depth searches, literal C = 1 tournaments for the depth-3, -4, -5 and -7 values, and no word for the 29 C1 misses |
| `explore/` | near-miss list to 2^28, residue and 3-adic means to 2^24, singleton windows for K = 0, 1 to 2^26, odd-multiplicity counts |
| `CARMICHAEL_MAP.md` | literature and repository map |

## Reproduce

    gcc -O2 -fopenmp -o rc_census rc_census.c && ./rc_census 1000000000 2 512 > rc_1e9.txt   # about 31 min on 12 threads
    python3 evaluate_prereg.py rc_1e9.txt
    gcc -O2 -fopenmp -o rc_depth rc_depth.c && ./rc_depth 1000000 7 512 1 && ./rc_depth 4000000 3 512 1
    gcc -O2 -o verify_chamber verify_chamber.c -lm && ./verify_chamber 1048576 14
    python3 literal_witnesses.py

The `explore/*.c` files are `rc_census.c` with only their output section changed (16-bit counters in `rc_mod3.c`,
`rc_v3.c` and `rc_par.c`).

## Line endings

The packet's `.gitattributes` pins LF. Outputs were produced on Windows (MinGW gcc 12-thread OpenMP).
