# H533 is reachable with two fixed ports at exactly three cycles

The [positive proof](FEEDBACK533_C3.md) supplies a ten-vertex graph with H=533, C=3 and feedback vertex number2, including its literal fixed-port word. The [independent review](audit/REPORT.md) passes without correction.

Together with the separately published [all-order negative theorem](https://github.com/eliottcassidy2000/math/blob/231bda1d4b2f57506bfcd6a2b2b3decb5f5f089d/session197/lowcycle/feedback533/README.md), this proves that three is the minimum Hamilton-cycle count among strong H533 tournaments of feedback number at most two. No minimum-order, complete C3-classification, or unrestricted fixed-port spectrum claim is made.

The negative proof is copied unchanged for the verifier's hash binding. Its complete census/audit reproduction remains in the linked separate checkpoint.

Run from this directory, using Python3 standard libraries:

    python3 verify_feedback533_c3.py
    python3 -O verify_feedback533_c3.py
    python3 audit/check.py
    python3 -O audit/check.py

The source verifier rewrites its adjacent certificate with identical bytes; use a copy if preserving a pristine checkout matters. The independent checker does not write files and must reproduce audit/independent.json. FILES.sha256 binds every included byte. No navigation or index is changed.
