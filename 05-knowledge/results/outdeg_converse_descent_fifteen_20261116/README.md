# The out-degree-1 law proved; the converse digit 0 ⟹ near-join; the dead-end descent; fifteen and twenty-five

Session klein-2026-11-16-outdeg (2026-11-16), on the owner's directive (the 15 moves on 5 labels, the 15 supersingular
primes, Monster and Baby Monster, the first 25 square-sum nodes in chunks of five, Fermat and Ogg primes, the first 25
primes; constructions that orient what they lose so that it recurs to the base case; "let the target choose the core or
admit a composite right source"; "try to prove the converse: digit 0 implies a near-join count") and the explicit request
"Prove the out-degree-1 law for digit-0 vertices of strong tournaments". Canon: THM-5264 (new); THM-5263 upgraded;
HYP-9453 updated.

## Status

- **PROVED** ([OUTDEG.md](OUTDEG.md), THM-5264). The two-path lemma (a strong tournament is partitioned by two paths, one
  from one of two given vertices and one ending at the other); the **out-degree-1 law** at every order (a digit-0 vertex
  has out-degree ≤ 1; exactly 1 in a strong tournament); the complete characterisation (digit 0 ⟺ sink, or sole exit
  whose port is outside the terminal component of `S − u`); **uniqueness** (a strong tournament of order ≥ 4 has at most
  one digit-0 vertex). Replayed: 58,650 vertices, 0 mismatches; 352,462 pairs for the lemma; uniqueness to order 9.
- **PROVED / FINITE-EXACT** ([CONVERSE.md](CONVERSE.md)). Dichotomy: a strong digit-0 tournament is a literal near-join
  (`k ≥ 1`) when its base has a sink (279 of 341 at orders 4–8), else a dead-end append over a strong right source `D`
  with `H = h(A)h(D) + B_{t_0}(A)P_D(2)`; all 62 such counts (orders ≤ 8) are near-join values from bases ≤ 7; **digit 0 ⟹
  near-join count for every strong tournament with `H ≤ 22,227`** (via THM-5263 and HYP-9453's census); the converse at
  all orders is exactly `{strong digit-0 counts} = {near-join values}`; OPEN above 22,227 in the strong-`D` case. Closure
  `{counts}·{near-join values} ⊆ {near-join values}` and the sink trick `(2^k+1)v` PROVED. Composite right sources are
  coprime linear readers `(h(D), P_D(2))`, but below 125 no value needs one.
- **FINITE-EXACT** ([DESCENT.md](DESCENT.md)). The dead-end descent `δ(S) = S − z` is a one-outward-edge system on strong
  tournaments with `H(S) = h(δS) + ε(S)`, escapes recurring to a sterile base; sterile strong counts of orders ≤ 9 = every
  odd in `[11, 199]` except 17, 21, 27; three kinds of counts below 125: sterile-only `11, 13, 47, 67` (the near-join
  exceptions), fertile-only `3, 5, 9, 17, 27` (permanent by Busch's floors 3, 5, 9, 15, 25, 45, 75, 125), mixed otherwise.
- **Typed** ([FIFTEEN.md](FIFTEEN.md)). 15 paths of `RT_5` vs 15 double transpositions: ANALOGY (three `Z_5`-orbits of
  five each); supersingular 15 induce a tree in the prime square-sum graph, the other ten below 100 one edge and all three
  isolated primes: NUMEROLOGY (THM-4860's verdict stands); the first 25 square-sum nodes in rows, residues, antidiagonal
  matchings, the unique path of 1..15 in chunks of five, Hamiltonian path counts 1,1,1,0,0,0,0,0,3,0,10 for `n = 15..25`;
  primes-vs-integers 25-sets: ANALOGY only.
- Pre-registration: [PREREGISTERED.md](PREREGISTERED.md) (14 predictions; P-4, P-5 and one clause each of P-6, P-10
  refuted, the rest confirmed).

## Files

- `scripts/`: `nearjoin_check.py` (helper block, copied from the previous packet; the scripts below `exec` its helpers),
  `twopath.py`, `converse_descent.py`, `converse_census.py`, `closure_check.py`, `squaresum25.py`, `fifteen.py`,
  `hampath_count.c`, `sterile9.c`.
- `receipts/`: `twopath.out`, `converse_descent.out`, `converse_census.out`, `closure_check.out`, `squaresum25.out`,
  `fifteen.out` (each with a byte-identical `.optimized.out` from `python3 -O`), `hampath_count.out`, `sterile9.out`.
- `prediction_outdeg.txt`, `SHA256SUMS`.

Reproduce: `cd scripts && python3 twopath.py` (minutes), `python3 converse_descent.py` (1 min), `python3 converse_census.py`
(5 min), `python3 closure_check.py`, `python3 squaresum25.py`, `python3 fifteen.py`, `cc -O2 hampath_count.c && ./a.out`,
`cc -O2 sterile9.c -o sterile9 && for n in 3 4 5 6 7 8 9; do gentourng -q $n | ./sterile9 $n; done` (gentourng from nauty).
