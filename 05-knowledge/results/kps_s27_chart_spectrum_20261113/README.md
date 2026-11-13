# The address chart cannot select every target; a map of where selection stands

kind-pasteur-2026-11-13-S27 · [THM-5152](../../../01-canon/theorems/THM-5152-the-3-adic-chart-class-omits-31-and-sixteen-more-counts-at-every-order.md) · [HYP-9679](../../hypotheses/HYP-9679-the-3-adic-chart-class-omits-infinitely-many-primes-including-every-mersenne-prime.md) · **Chan–Pak remains open.**

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
  multiplicativity of H.
- The chart class misses at **every** order:

    7, 21, 31, 35, 53, 113, 127, 149, 173, 181, 203, 217, 251, 263, 277, 287, 341, 343.

- Apart from 341 = 11·31 and 343 = 7³, these are primes and 7p values. Each, except 7 and 21, is an ordinary count;
  for example 31 is attained at order 6.
- A selector that stays inside the address chart therefore cannot prove Chan–Pak. HYP-9679 conjectures it misses
  infinitely many primes.
- A curiosity: the Mersenne primes 7, 31 and 127 are all omitted.

This was pre-registered (P1, 961f30e87b) and CONFIRMED. The order-10 data was reproduced independently by this
session's sweep. Replay: `python verify_s27.py` (pure Python, orders <= 8, passes under `-O`) and `scripts/`.

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
