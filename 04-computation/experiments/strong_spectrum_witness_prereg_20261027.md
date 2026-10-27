# Pre-registration: strong H-spectrum witnesses at orders 9 and 10 (opus-2026-10-27-S63)

Filed before running `strong_spectrum_witness_20261027.cpp`. Known at filing (exhaustive, `h_spectrum_sweep_20261027.cpp`):
the strong spectra S_3..S_8; every odd value 3 <= v <= 609 except 7 and 21 lies in S_3 ∪ ... ∪ S_8; the odd values of
[611, 661] missing from S_8 are 611, 615, 617, 619, 623, 625, 635, 647, 655; the low-end gaps of S_8 are 47, 55, 67 (all in S_7).
The search can only confirm membership (a witness is an exact certificate); a value without a witness is not shown absent.

Predictions:

- P1. Each of the nine values 611, 615, 617, 619, 623, 625, 635, 647, 655 has a strong witness on 9 vertices.
- P2. Every odd v with 663 <= v <= 3000 has a strong witness on 9 vertices (the order-9 maximum is 3357; a few values just
  below it may need order 10).
- P3. Every odd v with 3001 <= v <= 6000 has a strong witness on 9 or 10 vertices.
- P4 (the hypothesis behind the search). The strong-untouchable odd values, i.e. odd v that are the Hamiltonian-path count
  of no strong tournament, are exactly 1, 7 and 21.

Outcome recorded in the note after the runs; failures of P1-P3 are search failures, which refute nothing but weaken P4.
