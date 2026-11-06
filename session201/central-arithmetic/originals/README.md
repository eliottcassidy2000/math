# Bounded central-word arithmetic

Unreviewed proof candidate and explicitly OPEN all-target lemma.

`PROOF.md` isolates a finite-menu two-index arithmetic statement sufficient
for the desired positive odd path-count coverage. Its compiler uses actual
C1 graphs. It does not turn an inverse cutoff into an existence proof.

Reproduce the all-index matrix certificates, complete small physical
inventory, progression inverse, hostile length-seven controls, and literal
graph counts:

    python verify.py --cap 5001
    python -O verify.py --cap 5001

Reproduce the larger finite arithmetic scan:

    python verify.py --header central_forms.h
    g++ -O3 -std=c++17 inventory.cpp -o /tmp/central-inventory
    /tmp/central-inventory 10000001

The generated header contains no empirical constants. Its matrices follow
exactly from the three operations in the proof. The C++ program uses exact
64-bit integer arithmetic in its documented cap range.

Dependencies preserve the reviewed physical adjoint and the older forty-form
arithmetic comparison. The 27 existing repair graphs are recounted literally,
so their inclusion does not rely on an unreviewed numerical claim.
