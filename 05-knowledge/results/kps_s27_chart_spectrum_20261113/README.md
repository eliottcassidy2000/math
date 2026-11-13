# The address chart cannot select every target; a map of where selection stands

kind-pasteur-2026-11-13-S27 · audit PASS WITH CORRECTIONS ([audit](audit/AUDIT.md); all applied; MISTAKE-770) · [THM-5152](../../../01-canon/theorems/THM-5152-the-3-adic-chart-class-omits-31-and-eighteen-more-counts-at-every-order.md) · [HYP-9679](../../hypotheses/HYP-9679-the-3-adic-chart-class-omits-infinitely-many-primes-including-every-mersenne-prime.md) · **Chan–Pak remains open.**

The owner asked for creative work toward target-dependent selection, or a bypass. The prompt pointed at Session
312's chart: digit-reversed addresses give positive Plücker coordinates, and the leading 3-adic digits recover
ancestry and orientation. It also pointed at the Zeckendorf tricolour.

Session 313 works the same prompt in parallel and holds the Zeckendorf selector, the physical switches and the mod-3
bypass. This packet does two complementary things.

## 1. What the chart can reach (THM-5152)

Session 312 shows the chart's tournaments are exactly the vertex switches of the grammar of ordinal joins and cyclic
C3 joins. Without switching, the grammar's spectrum is Session 268's and is thin. With switching, the class is rich,
but it is not complete.

- Every grammar class through order 11 (7,048 classes, every switch) was computed.
- Three facts make this decisive below 375: Busch's minimum for strong tournaments, heredity of the class, and
  multiplicativity of H. The audit added explicit membership checks of ordinal-sum cores for composite values.
- The chart class misses at **every** order (corrected after the audit, which found 215, 267, 305):

    7, 21, 31, 35, 53, 113, 127, 149, 173, 181, 203, 215, 217, 251, 263, 267, 277, 287, 305, 341, 343.

- These are 11 primes, 5 values 7p, and 215, 267, 305, 341, 343. The chart class is not closed under ordinal sums,
  so a product of chart values can still be missed; 125 is reached only at order 12. Each, except 7 and 21, is an ordinary count;
  for example 31 is attained at order 6.
- A selector that stays inside the address chart therefore cannot prove Chan–Pak. HYP-9679 conjectures it misses
  infinitely many primes.
- A curiosity: the Mersenne primes 7, 31 and 127 are all omitted.

This was pre-registered (P1, 961f30e87b) and CONFIRMED. The audit (PASS WITH CORRECTIONS, two majors, MISTAKE-770)
reproduced the computation and supplied the core checks. Replay: `python verify_s27.py` (pure Python, orders <= 8, passes under `-O`) and `scripts/`.

## 2. The selection map ([SELECTION_MAP](SELECTION_MAP.md))

Built from three repository sweeps.

- **Reduction.** The prime/7p reduction is THM-4094 §6.
- **Near-complete families.** The best is THM-4875 (D), covering every odd N <= 2^32 except 7 and 21, with explicit
  digit functions and a fixed bank. The others are the O/Z/J grammar, the 2-adic lines and the matrix orbits.
- **Obstructions.** Each selector mechanism met one: fixed banks, the mod-96 determinant law, site budgets,
  absorption radii, the Collatz descent sieve, and now the chart.
- **Single-lemma reductions closest to a proof.** The ternary gate; two root changes above 3,780,967; finitely
  generated matrix orbits; two-step line supply.
- **The Zeckendorf tricolour** (your September 26 diagonal colouring of the triangle) has no proved link to path
  counts. Its useful feature is 2-4 typed representations per integer, which is exactly the flexibility selection
  lacks. inverse65 shows why golden digits do not fit the binary/tribonacci gadgets.
- **One route the repository has never tried: escalation / finite truth sets** in the style of Bhargava's 15/290
  theorems. The nearest existing mechanism is the Busch-bound finiteness used in THM-5152.
