# Unconditional cubic profile positivity

`PROOF.md` gives an exact positive sum-of-squares identity for the equal-profile
cyclic-quotient gap functional. It needs no log-concavity hypothesis. Together
with the independently derived exact query formula, it gives strict separation
for arbitrary three child tournaments with the same ordered-cover profile.

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

The checker is standalone and standard-library only. It verifies formal cubic
coefficient dictionaries through cutoff 20, all 97,648 nonzero sequences of
lengths 1–7 over 0–4, seven zero controls, and 1,500 signed identity controls.
It includes the essential upper boundary term at j=m+1. These checks support
the all-order algebraic proof; they do not replace it.

All six files here belong to the packet: PROOF.md, README.md, verify.py,
receipt.json, receipt_optimized.json, SHA256SUMS.
