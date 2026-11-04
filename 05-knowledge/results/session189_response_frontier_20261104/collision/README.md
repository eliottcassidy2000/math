# Session189 internal quartic distance tomography

This packet proves an exact transitive-module reconstruction theorem, pending independent review. It does not claim a full-response collision or arbitrary-prime closure.

- `PROOF.md`: all-order distance polynomial, coefficientwise positivity, fixed-ambient transfer, physical-order reconstruction up to converse for module order at least five, and the corresponding necessary restriction on original-module full-response collisions.
- `bank.cpp`: literal Hamilton subset-DP producer for complete named two-port response polynomials. `bank.json` retains every coefficient, including every named cubic and quartic zero, plus all supporting small-support actual query values, for 44 physical tournaments of old orders5,7,8.
- `verify.py`: independent literal run enumeration, exact endpoint-capacity DP, literal Hamilton-word polynomial expansion, formula checks, reconstruction census, and arbitrary-weight actual-query checks.
- `distance_polynomials.json`: every polynomial for orders2 through7, compared with literal marked-cut enumeration.
- `order_census_5.json`, `order_census_6.json`, `order_census_7.part0.json`, `order_census_7.part1.json`: complete physical-order censuses and every shared quartic used in reconstruction. These are shared-quartic records, not complete-response banks; the full-response banks are separately listed above.
- `ambient_controls.json`: exact joint exterior weights and every internal shared-name quartic for every rooted exterior pattern with at most three exterior names, in two transitive-module orders per pattern.
- `receipt.json`, `receipt.optimized.json`: identical normal/optimized verification receipts.

Reproduce with existing compiler and Python; no installation is needed:

    g++ -O3 -std=c++17 bank.cpp -o /tmp/s189_collision_bank
    /tmp/s189_collision_bank bank.regenerated.json
    cmp bank.json bank.regenerated.json
    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json
    cmp receipt.json normal.json

The initial draft's coefficient prefactor was corrected from (r+1)! to r+1 after two independent literal checks exposed double-counting of already-labeled run slots at n=4,r=2. The retained proof and data use the corrected factor. Coefficientwise positivity and the strict constant-term ambient transfer are unaffected. No initial-draft conclusion was promoted before this check.

The unbounded statements rest on the bijections and inequalities, not an extrapolation from these finite controls. Section5's use of all-cover restoration deserves separate scope review.
