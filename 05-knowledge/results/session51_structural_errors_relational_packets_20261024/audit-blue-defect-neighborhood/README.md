# Independent exact blue-neighborhood audit

Read `REPORT.md` for the PASS, full-arc independent method, and scope.
Every blue edge incident to the 24 positive-error order-ten self-converse
nodes is classified, including all boundary target classes. No zero/zero or
black edges are reconstructed.

    python3 /path/to/packet/independent.py > /tmp/neighborhood.json
    python3 -O /path/to/packet/independent.py > /tmp/neighborhood-O.json
    cmp /tmp/neighborhood.json /tmp/neighborhood-O.json
    python3 /path/to/packet/validate.py /path/to/session51-blue-defect-neighborhood
    python3 /path/to/packet/verify_manifest.py

The independent checker exhausts all within-stable-cell permutations, up to
14400 per target, rather than the producer's at-most-25 rotations. Its output
retains independent class IDs, complete matrices and memberships, all target
coordinate witnesses, and the exact correspondence to producer labels.
