# Observer profile algebra packet

Read THEOREM.md first. It distinguishes all-size linear failure from finite individual separation and the order-six quadratic repair.

Reproduce, using only the Python standard library:

    python analyze.py 6
    python profile_collision.py
    python direct_check.py
    python family_check.py
    python quadratic_check.py

Repeat each with python -O. Normal and optimized outputs in this packet are identical.

- walsh_certificate.json: every relevant signed graph orbit through n=6.
- collision_certificate.json: Walsh-based order-six class/profile separation.
- direct_certificate.json: independent literal 34-by-78 profile matrix, rational rank32, explicit probability-mixture ambiguity.
- family_certificate.json: exact polynomial integration controls, supplementing the all-size proof.
- quadratic_certificate.json: exact rank34 integer minor using two squared coordinates.
- edge_actions.py: local pure helper dependency; no external packet is needed at runtime.
- SOURCE_PINS.json: inherited mathematical source hashes.

Independent audits are maintained separately by their auditors. No code here fetches the network, modifies the repository, or depends on an unrecorded package.
