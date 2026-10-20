# Exact leaf reduction for the paired-black quotient

Publication update, 2026-10-20: This is a preserved precursor. The later [complete deficit-two center classification](../deficit_two_centers_20261020/README.md) proves D(z)=z^4. The unknown-center and OPEN language below records the earlier boundary, not the current status; the leaf reduction remains valid.

This packet upgrades the former constructed NS leaf family to an all-order
classification and reduces every mixed leaf to a rigidly wrapped strong
deficit-two center. It does not claim all bridges are leaf edges.

- `BLACK_QUOTIENT_LEAF_CLASSIFICATION.md`: exact iff classification,
  generating series in the unknown center alphabet, sharp degree-three
  boundary, order-nine certificate, and bridge scope.
- `CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md`: universal cycle-cut bounds;
  excludes deficit-two centers at odd n>=7 unless 3 divides n and at
  all n divisible by eight.
- `check_deficit_two.py`, `verification.json`, `verification.out`:
  exact subset-DP test of every SC class at n=3,5,7,9.
- `audit/AUDIT.md` and its independent scripts/results: approved proof
  audit and complete independent SC9 verification.

The 2,752-class n=9 test needs only five automorphism checks, all on
nonstrong pure-blue candidates. All 2,661 strong classes fail earlier
necessary filters. No new deficit-two center exists at n=9.

Status: structural theorems PROVED from stated inherited inputs;
computation FINITE-EXACT; the full strong center alphabet OPEN.
The only known center is S4; its uniqueness beyond order nine is not proved.

The source SC catalog is published in the sibling
`odd_blue_spectrum_hollow_cubes_20261020` packet. From this directory:

    python3 check_deficit_two.py --source ../odd_blue_spectrum_hollow_cubes_20261020 \
      --output verification.json > verification.out

Run the independent audit from this directory:

    python3 audit/independent_audit.py --source ../odd_blue_spectrum_hollow_cubes_20261020
    python3 audit/check_closing_prefixes.py

The completed review checks cycle-cut orbit counting, closure of the
prefix witnesses, the safe direction of automorphism-order filters, and
the distinction between leaf enumeration and bridge classification.
