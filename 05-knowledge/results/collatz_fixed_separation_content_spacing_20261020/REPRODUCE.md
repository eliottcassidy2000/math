# Reproduction and evidence boundary

All scripts use Python's standard library and explicit runtime checks. Run from this packet directory:

    python3 verify_fixed_separation.py > receipt.replay.json
    python3 -O verify_fixed_separation.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json
    python3 certify_spacing_constants.py > spacing_constants.replay.json
    python3 -O certify_spacing_constants.py > spacing_constants.optimized.replay.json
    cmp spacing_constants.json spacing_constants.replay.json
    cmp spacing_constants.json spacing_constants.optimized.replay.json

The byte-identical imported absolute-gap audit runs locally:

    python3 inherited-absolute-gap/audit/audit_endpoint_large_range.py > absolute_gap.replay.json
    python3 -O inherited-absolute-gap/audit/audit_endpoint_large_range.py > absolute_gap.optimized.replay.json
    cmp inherited_absolute_gap.replay.json absolute_gap.replay.json
    cmp inherited_absolute_gap.optimized.json absolute_gap.optimized.replay.json

That script reads the adjacent `inherited-absolute-gap/large-range/` producer files, validates their source hash and matching normal/optimized receipts, and independently reconstructs all eighteen convergents using log3=(3log2+log(9/8))/2 and 112 positive-series terms. It checks the numerical specialization and consequences of the cited Matveev theorem; it does not reprove that theorem.

The principal note supplies the all-length proof. The physical tests are finite controls, not the proof of the unbounded claim and not a complete enumeration below the fixed-h cutoff. `SHA256SUMS` covers all publication files except itself; verify it before making replay outputs in the same directory.
