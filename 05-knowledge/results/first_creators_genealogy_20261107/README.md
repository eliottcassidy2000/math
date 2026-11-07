# First creators of odd Hamiltonian-path counts: converse fibers, the value colours, and genealogy

kind-pasteur-2026-11-07-S16. **THM-4922** (PROVED / FINITE-EXACT through order 10) and **HYP-9512**. The independent
audit is PASS WITH CORRECTIONS ([audit/AUDIT.md](audit/AUDIT.md); MISTAKE-735). Pre-registration
[PREREGISTRATION.md](PREREGISTRATION.md) was pushed in 9e07033ca4 before the order-10 run. Chan-Pak remains OPEN.

Owner prompt: merge the hypotraceable and hypohamiltonian work, and follow, for every odd number except 7 and 21, the
structure that first creates it as the order grows. The conjecture is that this structure is a single isomorphism class
or exactly two converse classes. The odd values would then split into two kinds, and larger structures might be
generated from smaller ones.

**Concurrent same-prompt work, credited.** Session224 (`../session224_first_occurrence_20261107/`, pushed first) has
the refutation at h = 9, the (s, p) birth ledger through order 8, and the strong-component structure of minimal
realizations.

## Results

1. **The dichotomy is false.**
   - The smallest failure is h = 9: two self-converse classes at order 5 (proved first by Session224).
   - h = 23 is the smallest value with two converse pairs.
   - The dichotomy holds for 69 of 1,520 values through order 9 and 55 of 6,143 at order 10, mostly at small values and
     near each order's maximum.
   - Fibers reach 257 classes through order 9 (h = 1503) and 2,458 at order 10 (h = 7799).
2. **The value colours.**
   - A value is blue if some first creator is self-converse, and black if not.
   - At even least orders every value 3 (mod 4) is black (THM-791). Elsewhere 82% are blue in aggregate (0.72, 0.66,
     0.84, 0.825 at orders 7-10).
   - **Proved (audit):** at every odd order the self-converse classes split exactly evenly between H = 1 and 3 (mod 4).
3. **Genealogy (the owner's generation idea, tested through order 10).**
   - At orders 4-10, every first creator except the two first creators of 63 (the Fibonacci fan F_6 and its converse,
     order 8) contains a vertex whose deletion is a first creator one order down.
   - So, among values h >= 5 with n(h) <= 10, 63 is the only one not grown from smaller first creators.
   - Proved above n max H(n - 2), from THM-4889's H(T) <= sum over v of H(T - v).
4. **Hypohamiltonian merge.**
   - No non-first-creator of order 3-10 has every vertex deletion a first creator, mirroring "no hypohamiltonian or
     hypotraceable tournaments".
   - The fully non-critical share of first creators rises: 24%, 39%, 50%, 62.8% at orders 7-10.
5. **Pre-registered order 10.**
   - Q1 (no new orphan), Q2 (no hypo analogue) and Q3 (full non-criticality above 50.4%): CONFIRMED blind.
   - Q4 (the dichotomy rare): CONFIRMED, but pre-determined by klein's fiber sizes.
   - Q5 (mixed fibers the majority): FAILED, forced by THM-791.

## Files

| File | Role |
| --- | --- |
| `classgen.c` | class generator to order 9, with H, c3, strong, SC, \|Aut\| and vertex-deletion classes; validated against A000568 and klein's census at orders 3-9 |
| `order10.c`, `order10.out` | the order-10 census and its aggregated first-creator statistics (lines "V h fiber #SC #pairs #fully-non-critical") |
| `analyze.py`, `hypo.py`, `split.py` | the order <= 9 analyses |
| `audit/` | independent audit: its own enumerator (`aud_enum.c`), the proof of the odd-order balance (`aud_sc_mod4.py`, `aud_sc_orbits.py`) and its outputs |

## Reproduce

    gcc -O2 -o classgen classgen.c && ./classgen 9 cls        # 37 s; writes cls_1..9.txt (30 MB at n = 9)
    python3 analyze.py 9 && python3 hypo.py && python3 split.py
    gcc -O2 -o order10 order10.c && ./order10 cls_8.txt cls_9.txt > order10.out   # 12.5 min
