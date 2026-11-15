# Chan–Pak: does the pair spectrum (H, P₂) cover every line? Yes for a fixed kernel through order 10, and across the junctions to 1.2–1.6 million

Session klein-2026-11-15-lines. **Chan–Pak remains OPEN.**

Owner prompt (2026-11-15): "test whether the pair set (H, P₂) covers every line at the largest orders and other
strategies to prove chan-pak".

[Wave 35](../chanpak_waves_20261113/README.md#wave-35-positivity-confinement--no-kernel-gives-cancellation-so-prime-targets-have-finitely-many-port-pairs)
reduced prime targets to *pair-spectrum fullness*: the pairs (H(A), P₂(A)) should meet every line αx + βy = p of
one fixed kernel. This packet tests that directly ([LINES.md](LINES.md)).

## Findings

1. **For a fixed kernel the question is one-dimensional, and it is already true at order 10 (FINITE-EXACT).**
   Each pair lies on one line, so covering all lines means the single statistic αH + βP₂ is onto the admissible N.
   - The census of all 9,733,056 order-10 classes gives the following. With the one-vertex reader, 3H + P₂ meets
     **every** line N in [36,835, 241,521], a window of ratio 6.6, with at least 4 pairs per bulk line. With the
     two-vertex reader, 5H + 3P₂ meets every line in [188,269, 623,231].
   - So what is needed is σ > 1 (pairs per line), not the σ > 2 of full two-dimensional occupancy. That holds from
     order 9–10.
2. **Odd readers see only every other order for two thirds of the targets (PROVED from THM-5172 (C)).**
   - With an odd-order reader, odd-order cores give only N ≡ 0 (mod 3).
   - So the classes N ≢ 0 (mod 3) need windows two orders apart to overlap, a scale jump of about 37.
   - Even readers, or alternating readers, remove this.
3. **The junctions close (WITNESSED).** Tilted sampling of orders 11 and 12 (THM-5127 (B)'s ensemble) gives
   explicit witnesses, recounted independently.
   - Order 11 meets every K₁ line N ≡ 3 (mod 6) in [47,385, 1,294,989].
   - Order 12 saturates the hard even-order junction above 241,521: its empty lines fall 505 → 7 as 4.8M targeted
     samples arrive.
4. **Junctions are reducible, hence partly exact (NEW).**
   - *Mechanism.* Near-transitive low tails are mostly ordinal sums, and (H, F₂) multiplies under A ⇒ B. So the
     reducible pair sets of orders 11–14 are computed **exactly** from the census.
   - *Exact gains.* K₁ then misses no odd multiple of 3 up to 1.4·10⁶, and TT₂ closes to 818,949 without sampling.
   - *Combined.* Exact and witnessed data together cover every odd N in **[36,835, 1,238,097]** with the single
     kernel K₁, and **[188,269, 1,597,157]** with TT₂.
5. **A reduction to finitely many 1D statements (PROVED reduction).**
   - Padding a core by TT_k turns K₁ into the derived statistic H + 2^k F₂.
   - So single-kernel coverage of all large N follows from bulk coverage of the statistics H + 2^k F₂, for k ≤ 4,
     at every order. That is finitely many one-dimensional tilted local limit statements per order, strictly weaker
     than the 2D band fullness of HYP-9692.
   - All residues mod 3, 5, 7, 9, 16 and 25 occur, so there is no local obstruction; there is a mild mod-3 bias
     near the transitive end.
6. **New conjecture [HYP-9772](../../hypotheses/HYP-9772-single-kernel-conjecture-every-odd-n-from-36835-is-h-plus-f2-of-some-tournament-the-double-cone-count-and-fixed-kernel-line-coverage-is-one-dimensional.md).**
   - *Statement.* Every odd N ≥ 36,835 is H(A) + F₂(A), the Hamiltonian-path count of the double cone a → A → c → a.
   - *Record.* Certified to 1,238,097.
   - *Relation to Chan–Pak.* It implies Chan–Pak beyond the 2³² census.

## Replay

1. **Census.** For n in 2..10, run `gentourng -q n | ./paircensus n | sort -u > pairs_n.txt`. Order 10 takes about
   2 minutes split eight ways with `gentourng -q 10 r/8`.
2. **Samples.** Run `cc -O3 -o tiltsample tiltsample.c`, then regenerate the samples with the parameters in
   [receipts.txt](receipts.txt) (`run_junction12.sh` for the targeted run).
3. **Analyses.** Run `python3 linecov.py`, `multiplicity.py`, `junction11.py`, `reducible.py`, `union_k1.py`,
   `union_ext.py` and `equidist.py 12`.
4. **Witness checks.** `python3 verify_witnesses.py n q seed index ...` rebuilds any witness by replaying the
   sampler's random stream in Python and recounts H and P₂ with its own DP.

Credits:
- the reduction to lines is wave 35's, and the cyclic-join identity is THM-5222 / THM-5172 (C);
- the tilted ensemble is THM-5127 (B);
- wave 25 introduced the ridge band and the reader and padding windows; the derived-kernel chain is their mirror
  image on the core side.
