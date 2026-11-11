# Reproduce Session269

Python3 and SymPy are needed. No network is used by the research programs.
From this directory:

```bash
python3 reproduce.py
```

This runs each script normally and with `-O`, requiring byte-identical stdout.
Correctness tests use explicit exceptions, so optimized mode does not remove
checks. The final runner report is `FINAL_CHECKS.json`. The scripts' existing
JSON/text receipts are small, scoped mathematical controls, not global censuses.

Individual commands, run in the indicated subdirectory:

- `jacobian`: `python3 verify.py` checks symbolic flat connections, the explicit
  rational-certificate compiler, denominator locations and a punctured-surface
  hostile. It does NOT test the inherited analytic comparison theorem.
- `hadwiger`: `python3 replay.py` checks all91 pairs of the14-point core
  suspension and exact colouring channels, including the balanced triangle.
  `python3 equality.py` checks all21 pairs of the7-point equality gadget and
  all fixed-port colourings at q=3,4. The K4 uncolourability control is an
  abstract control, not a claim of unit-distance plane realization.
- `chanpak`: `python3 reproduce.py` runs32,764 residue targets for k3–15,
  58 literal graph controls and80 high-precision arithmetic targets. Literal
  subset DP and structural checks are distinct from the larger arithmetic-only
  tests. No independent rigidity census is claimed.
- `bridges`: `python3 check.py` runs the stated two-state positive/negative
  controls through60 transitions, including signed cancellation.

A concrete tournament witness:

```bash
cd chanpak
python3 reproduce.py --target 63 --bits 8 --graph
```

This prints adjacency rows and the ACTUAL integer H. The target parameter
means congruence modulo2^bits, not equality with63.

All-row, all-angle and arbitrary-precision conclusions are supported by proofs
in the corresponding notes, not extrapolated from these finite checks.
