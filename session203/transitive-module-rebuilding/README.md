# Physical doubling-plus-one rebuilding families

**PROVED, independently audited in the stated scope.** Three explicit unbounded families of tournament rewrites satisfy H(child)=2H(parent)+1. Two replace one vertex, keep the same order, and take a C1 parent to a C3 child, where C counts directed Hamiltonian cycles modulo rotation.

- Family B: parent H=20*2^k+5 for every k>=0; replacement H=40*2^k+11.
- Family C: parent H=27*2^k-1 for k>=1; replacement H=54*2^k-1. At k=0 the separate values are 25 and 51.
- Family A retains every old vertex, reverses exactly one old arc, and adds one vertex. For k>=1 its parent count is 8*3^k+17*2^k+1 and its child count is twice that plus 1. At k=0 the values are 25 and 51.

All final parents and children are strong. The same-order replacements explicitly delete an old vertex after it becomes a universal sink. Family B is not closed under iteration: its child count is six larger than the next parent count.

These are sparse physical construction families, not a threshold theorem, a complete residue class, or an induction covering arbitrary odd integers. No novelty of the scalar counts is asserted. The all-order identity follows from an exact module-run expansion and a finite core coefficient census, not extrapolation from bounded graph tests.

[PROOF.md](PROOF.md) specifies all literal graphs and masks and proves the module expansion. [audit/REVIEW.md](audit/REVIEW.md) records independent allowed/mandatory-gap enumeration and literal path/cycle DFS on 48 sample graphs, checking over 1.1 million path prefixes. The producer independently uses gap-polynomial multiplication and subset DP. Both implementations agree in normal and optimized Python.

## Reproduction

    python3 verify.py > /tmp/module_producer.json
    python3 -O verify.py > /tmp/module_producer_optimized.json
    cmp /tmp/module_producer.json NORMAL.json
    cmp /tmp/module_producer_optimized.json OPTIMIZED.json
    python3 audit/check.py > /tmp/module_audit.json
    python3 -O audit/check.py > /tmp/module_audit_optimized.json
    cmp /tmp/module_audit.json audit/NORMAL.json
    cmp /tmp/module_audit_optimized.json audit/OPTIMIZED.json
    sha256sum -c SHA256SUMS

RELOCATION.json records exact release-local replays. The original reviewed proof and status-only promotion receipt are preserved.

## Dependencies and scope

The proof is self-contained: it uses elementary directed paths, universal-sink deletion, and surjections counted by r!S(k,r). There is no dependency on an exhaustive tournament bank, a numerical order cutoff, or Busch's lower bound. Both scripts use only the Python standard library.

The inherited motivation is the separately published all-order induced H25-in-H51 obstruction and its finite repair, at ../induced25-51-rebuilding/. That theorem is not needed to prove these three module identities. Its no-deletion ancestry restriction must not be applied to the same-order vertex replacements here. No minimal old-arc edit count is asserted when deletion is allowed.
