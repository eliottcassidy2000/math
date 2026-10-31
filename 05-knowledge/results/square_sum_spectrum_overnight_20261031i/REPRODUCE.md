# Reproduction

Run from the repository root. Python 3 standard library only; the inherited packet-f atlas and packet-h carrier implementation must be present at the pinned versions in `MANIFEST.json`.

```bash
python3 05-knowledge/results/square_sum_spectrum_overnight_20261031i/one_copy_ports.py
python3 05-knowledge/results/square_sum_spectrum_overnight_20261031i/audit_integer_clock.py
python3 05-knowledge/results/square_sum_spectrum_overnight_20261031i/elliptic_ports.py
python3 05-knowledge/results/square_sum_spectrum_overnight_20261031i/audit_projective.py
python3 05-knowledge/results/square_sum_spectrum_overnight_20261031i/ratio_obstructions.py
python3 05-knowledge/results/square_sum_spectrum_overnight_20261031i/weighted_join.py
python3 05-knowledge/results/square_sum_spectrum_overnight_20261031i/verify_packet.py
python3 agents/check_docs.py
git diff --check
```

Each computational command prints the corresponding retained JSON. Every script was also run with `python3 -O`; the normal and optimized outputs are byte-identical. Validation uses explicit checks, not removable assertions. `VERIFICATION.json` records receipt hashes.

The independent paths are rational versus homogeneous-integer doubling, rational-affine versus integer-projective group arithmetic, chart formulas versus literal vertex pairs, and a symbolic weighted product versus the sixteen explicit edge-set products. The complete six-free-bit fibre in the first native square-class-one join additionally recovers both source cycles after undoing the declared scalings. General cycle completeness follows from the inherited input-carrier theorem and the two-edge-cut inverse proof, not from sampling these finite fibres.

No script is an unrestricted tournament enumeration or a proof of initial-interval Hamiltonicity. The delayed collision control uses arithmetic source sets without an asserted cycle structure.
