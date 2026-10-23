# Session38 electrical results

The full order-eight complement-line network has exact 2-primary critical
group `(Z/2)^3358 + (Z/4)^4 + Z/256`. Its spanning-tree count has 2-adic
valuation 3374. This is obtained from the full 2^21-port wiring ledger and
certified dyadic lifts, not a determinant approximation or precision guess.
It disproves continuation of the full-network exponent-four pattern seen
through order seven; no all-order exponent bound is proposed.

`PROOF.md` contains the all-size lift tower, the exact rank-two switch
criterion, a degree-preserving counterfeit invisible to the whole dyadic
tower, and the distinction between electrical constraints and authenticated
port acceptance. The same file gives the finite order-eight interpretation
and explains the lower-bit corrections used in the computation.

Main reproducible files:

- `check_electrical.py`, `receipt.json`: independent fixed-path censuses
  through order six, small Smith checks, and all 42 eligible order-five
  switches. No inherited producer imports.
- `n8_electrical_receipt.json`: full-network first-half-form projection.
- `check_n8_lift.py`, `n8_lift_receipt.json`: five certified mod-four lifts
  and the full second half-form.
- `lift_final_radical.py`, `n8_final_radical_receipt.json`: certified mod-eight
  lift of the final binary radical.
- `complete_n8_tower.py`, `n8_complete_tower_receipt.json`: every remaining
  lift through modulus 256, all integer Laplacian images, and termination.

The direct certificate audit can ignore how lifts were found: reconstruct
each vector, multiply by the full integer Laplacian, check the recorded
image, and compute the one-dimensional energy. The independent parent audit
uses its own endpoint reconstruction and a different grounding for q1.

Local default ledger location is the sibling `session38-wiring` directory.
`--base PATH` is accepted for the lift producer and tower completion;
`check_electrical.py --n8 PATH` accepts a ledger path. `lift_final_radical.py` also accepts `--base PATH`. Normal and optimized runs
have matching arithmetic outputs. These scripts do not mutate Git.
