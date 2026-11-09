# First creators at order 11, and how every tournament is generated from them

kind-pasteur-2026-11-08-S17. **THM-4942** (PROVED / FINITE-EXACT / KNOWN / ANALOGY typed) and **HYP-9532**. The
independent audit is PASS WITH CORRECTIONS ([audit/AUDIT.md](audit/AUDIT.md); MISTAKE-737). It reproduced every
number at orders 3-11 with its own code. Pre-registrations: [PREREGISTRATION.md](PREREGISTRATION.md) (R1-R9, af0112b596) and
[PREREGISTRATION-2.md](PREREGISTRATION-2.md) (R10-R11, c201ec7b34). Proofs and tables are in [PROOF.md](PROOF.md).
Chan-Pak remains OPEN.

Owner prompt (2026-11-08): run the order-11 test. For every odd number besides 7 and 21, keep track of the structure
that first created it as the order grows. See how all tournaments can be generated from this minimal set if you get
creative enough with the operations. Relate the objects used in the proofs that 7 and 21 never appear to fundamental
concepts and constructions.

## Results

1. **The engine.**
   - H(T + v) is a quadratic form in the attachment pattern of v, with split-path coefficients. This is KNOWN:
     HYP-2879, THM-4097 and THM-4114 (PROOF section 1).
   - One Gray-code walk gives every extension of a class.
   - Run on all 9,733,056 order-10 classes (10^10 extensions), it reproduces klein's S(11) exactly.
2. **Order 11 (pre-registered R1-R5, all CONFIRMED).**
   - Each of the 36,248 values new at order 11 is created by a one-vertex extension of an order-10 first creator. So
     63 is still the only orphan value through order 11.
   - No hypo analogue exists at order 11.
   - Every order-10 first creator is fertile, with at least 564 of its 1,024 extensions first creators.
   - At orders 3-10 the only first creators with fewer than half their extensions fertile are the two creators of 63.
   - Class level, settled by the independent audit (not pre-registered): every order-11 first creator has a
     first-creator deletion. There are 0 orphans among 506,912,096 candidates. So through order 11, every first
     creator except 63's two lies on a chain of first creators back to the 3-cycle.
3. **Generating every tournament (orders 3-11).** Every class is one of the following:
   - a first creator: 88.8% of all classes at order 10;
   - a one-vertex extension of a first creator;
   - obtained from strictly smaller tournaments by an ordinal sum, a substitution, or an ordinal sum with at most
     two arcs reversed.

   Each operation has an explicit H-calculus. The pre-registered "pendant" refinement failed (R7, R8): the right
   invariant is the arc-strong connectivity lambda.
4. **Connectivity floors.**
   - Let f_k(n) be the least H over n-tournaments with lambda >= k. f_1 is the Moon-Busch floor, the object that
     excludes 7 and 21.
   - Through order 10, f_3(n) > max H(n - 1). Hence every tournament whose scores all lie in [3, n - 4] is a first
     creator. (THM-4811 (B): for n <= 13, lambda >= 3 iff min degree >= 3.)
   - At order 11 this fails (R10, R11 CONFIRMED, by exhaustive search): f_3(11) = 11,475 < 15,745 = max H(10).
     Four of the seven 3-arc-strong values below 15,745 already occur at order 10: 11,475, 15,213, 15,273 and
     15,561.
   - Yet the five-operation grammar survives at order 11. All 7 such non-first-creator classes are substitutions.
     The floor is C3[R5, R5, 1]: the 3-cycle, the regular first creator R5 of 15 (twice), and a singleton.
5. **Fibers.**
   - Strong realizations of a value form a finite antichain (proved).
   - Among all values in S(11), so every odd h <= 80,451, the values born at their Busch bound are exactly 3, 5, 9,
     11, 13, 17, 19, 23, 35, 39 and 63. The 3,635 values of [80,453, 95,093] outside S(11) are conditional on
     HYP-9481 (a).
     The orphan values 3 and 63 are its smallest and largest members.
6. **The 7/21 proof objects** are the generation operations:
   - condensation is SUM;
   - insertion is EXT;
   - Moon's theorem gives the first-creator chains;
   - the floor becomes the graded floors.

   Their analogues are typed as ANALOGY: unique factorization, insertion sort, non-cut vertices, the Minkowski and
   Hermite bounds, and the hard-core gas (PROOF section 5).

## Files

| File | Role |
| --- | --- |
| `dump10.c` | all order-10 classes with H, deletion counts and strongness (from S16's `cls_9.txt`; 331 MB binary, not committed) |
| `ext.c` | the extension engine: attachment form, Gray walk, deletion forms, hypo and fertility counts, brute-force check mode |
| `ext_wit.c` | witness variant: one (first-creator parent, pattern) per new value |
| `analyze11.py`, `order11_summary.txt` | evaluation of R1-R5, with the per-slice logs |
| `gen.c`, `census_3_10.txt` | the generation census at orders 3-10. Lambda, pendant and cut-vertex diagnostics were added after the R7/R8 evaluation |
| `lam.c`, `floors_4_10.txt`; `mindeg.py`, `mindeg10.c`, `mindeg_4_10.txt` | connectivity floors and degree tables |
| `lam11.c`, `lam11_summary.txt`; `lam11_all.c`, `lambda3_witnesses_order11.txt` | the exhaustive order-11 search of 3-arc-strong tournaments with H <= 15,745 (R10, R11), and all 97 witnesses (11 classes) |
| `witnesses_order11.txt` | one witness per value new at order 11: an order-10 first creator and a pattern (R1) |
| `lineage.py`, `lineage_le10.txt` | the value-level genealogy ledger for all 7,661 values with n(h) <= 10: first creators, how many have a first-creator parent, and the parent values; 63 is the only ORPHAN line |
| `merge.py` | slice merging for the validation runs |
| `graft_check.py` | checks the one-back-arc graft identity H = H(A)H(B) - F_A(a)G_B(b) + U_A(a)U'_B(b) against direct DP |
| `audit/` | the independent audit |

## Reproduce

    gcc -O2 -o dump10 dump10.c && ./dump10 cls_9.txt classes10.bin                   # about 13 min; cls_9.txt from S16 classgen.c
    gcc -O2 -o ext ext.c && ./ext 10 bin classes10.bin <lo> <hi> S9.txt S10.txt o11_k 2  # 11 slices, about 20 min
    python3 analyze11.py
    gcc -O2 -o gen gen.c && ./gen 10 bin classes10.bin S8.txt S9.txt 0
    gcc -O2 -o lam lam.c && ./lam 10 bin classes10.bin S9.txt
    gcc -O2 -o lam11 lam11.c && ./lam11 10 bin classes10.bin <lo> <hi> S10.txt 15745 l11_k
