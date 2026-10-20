Complete merged metagraph realizability

Publication update, 2026-10-20: This is a preserved precursor. The later all-order Hamilton-path relocation theorem proves necessity as well as sufficiency of the local forest criterion, retaining every coordinate, complement, colour and exact class-count condition. See ../hamilton_path_relocation_connectivity_20261020/README.txt and LOCAL_FOREST_ADDENDUM.txt in that sibling directory. The finite-through-n8 wording below describes this earlier packet, not the current theorem boundary.

Read REALIZABILITY.txt first. It contains the full plaintext proof and exact
claim boundaries.

Main results:
1. An all-order necessary-and-sufficient coordinate/transport certificate,
   with a no-merger argument based on the exact Burnside node count.
2. An all-order sufficient local relocation/reflection forest certificate.
   Full local scans certify exactness through n=8, without a tournament atlas.
3. An n=5 false wiring that is a fully coherent Boolean cube partition and
   preserves the earlier endpoint/score-shell and signed-charge constraints.
4. For explicitly listed binary weighted inputs: polynomial-size certificates,
   NP membership, and reduction to four-coloured simple graph isomorphism.
   This does not give polynomial-time recognition. See COMPLEXITY_ADDENDUM.txt.

Evidence:
- COMPLEXITY_ADDENDUM.txt (separate encoding-sensitive complexity proof)
- verify_realizability.py, verify.out, audit.json
- certificate_n3.json, certificate_n4.json, certificate_n5.json
- local_move_classifier.cpp, local_move_classifier.out
- independent_local_check.py, independent_local_check.out
- independent_n8_scan.cpp, independent_n8_scan.out
- local_assignment_n3.txt through local_assignment_n7.txt

The first verifier uses only Python standard library and SymPy. The local
scan uses a C++17 compiler. The Python direct-arc audit independently checks
its partition through n=7, and a second direct-arc C++ implementation
independently repeats the complete n=8 scan. Exact n=8 local connectivity is finite evidence, not
an all-order connectivity theorem.
