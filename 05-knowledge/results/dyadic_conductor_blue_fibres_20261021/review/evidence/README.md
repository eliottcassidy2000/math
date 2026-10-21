# Current conductor audit

FINAL PASS. This self-contained package contains only the approved current dyadic-conductor results and their current verification.

- Producer: 28,754 checks, normal and optimized output identical.
- Independent checker: 46,719 checks, normal and optimized output identical.
- INDEPENDENT_AUDIT.md: mathematical proof, source, and blue-fibre interface audit.
- FINAL_INPUT_PINS.json and current_inputs/: all active independent-checker dependencies.
- PUBLIC_APPROVAL.json: approved files and exact hashes.

Run from this directory:

    python3 -B independent_checker.py > /tmp/idoneal-independent.json
    python3 -B -O independent_checker.py > /tmp/idoneal-independent-O.json
    cmp /tmp/idoneal-independent.json /tmp/idoneal-independent-O.json
    cmp /tmp/idoneal-independent.json independent.normal.json
    python3 -B current_inputs/verify_conductor_tower.py > /tmp/idoneal-producer.json
    python3 -B -O current_inputs/verify_conductor_tower.py > /tmp/idoneal-producer-O.json
    cmp /tmp/idoneal-producer.json /tmp/idoneal-producer-O.json
    cmp /tmp/idoneal-producer.json producer.final.normal.json
    sha256sum -c SHA256SUMS

Keep FINAL_INPUT_PINS.json beside independent_checker.py; it selects the included current_inputs directory. No source files outside this package are needed at runtime.
