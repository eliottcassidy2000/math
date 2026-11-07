# First creators of odd Hamiltonian-path counts: converse fibers, the two kinds, and genealogy

kind-pasteur-2026-11-07-S16. **THM-4922** (PROVED / FINITE-EXACT through order 10) and **HYP-9512**. Pre-registration
[PREREGISTRATION.md](PREREGISTRATION.md) was pushed in 9e07033ca4 before the order-10 run. Chan-Pak remains OPEN.

Owner prompt: merge the hypotraceable and hypohamiltonian work, and follow, for every odd number except 7 and 21, the
structure that first creates it as the order grows. The conjecture is that this structure is a single isomorphism class
or exactly two converse classes. The odd values would then split into two kinds, and larger structures might be
generated from smaller ones.

## Results

1. **The dichotomy is false.** It holds for 69 of 1,520 values through order 9 and 55 of 6,143 at order 10, only at
   small values and near each order's maximum. The smallest failure is h = 23: two converse pairs at order 6. Fibers
   reach 2,458 classes.
2. **The two kinds.**
   - Blue values have a self-converse first creator; black values have only converse pairs.
   - At even least orders, every value 3 (mod 4) is black (THM-791). Elsewhere about 83% are blue.
   - New finite observation: at odd orders the self-converse classes split exactly evenly between H = 1 and 3 (mod 4).
3. **Genealogy (the owner's generation idea, made exact).**
   - Through order 10, every first creator except the two first creators of 63 (order 8) contains a vertex whose
     deletion is a first creator one order down. So 63 is the only value not grown from smaller first creators.
   - Proved: this always holds above n max H(n - 2), by the endpoint bound H(T) <= sum over v of H(T - v).
4. **Hypohamiltonian merge.**
   - No non-first-creator has every vertex deletion a first creator, mirroring "no hypohamiltonian tournaments"
     (THM-4870).
   - The fully critical share of first creators rises: 24%, 39%, 50%, 62.8% at orders 7-10.
5. **Pre-registered order 10.** Q1 (no new orphan), Q2 (no hypo analogue), Q3 (full criticality above 50.4%) and Q4
   (the dichotomy rare) CONFIRMED. Q5 (mixed fibers the majority) FAILED, explained by THM-791.

## Files

| File | Role |
| --- | --- |
| `classgen.c` | class generator to order 9, with H, c3, strong, SC, \|Aut\| and vertex-deletion classes. Validated: A000568 and klein's census at orders 3-9 |
| `order10.c`, `order10.out` | the order-10 census and its aggregated first-creator statistics (lines "V h fiber #SC #pairs #fully-critical") |
| `analyze.py`, `hypo.py`, `split.py` | the order <= 9 analyses |
| `PROOF.md` | statements, tables, proofs |

## Reproduce

    gcc -O2 -o classgen classgen.c && ./classgen 9 cls        # 37 s; writes cls_1..9.txt (30 MB at n = 9)
    python3 analyze.py 9 && python3 hypo.py && python3 split.py
    gcc -O2 -o order10 order10.c && ./order10 cls_8.txt cls_9.txt > order10.out   # 12.5 min
