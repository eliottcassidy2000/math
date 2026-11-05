# A physical three-cover gate with an exact restricted inverse

The reviewed `THREE_COVER_CYCLIC_GATE.md` derives the positive bilinear kernel for C3[T,U,Tr2] directly from maximal colored runs. It specializes the pinned S195 colored-cover algebra and carefully distinguishes unordered from ordered covers.

Pairing U=T^op gives a positive quadratic form, a strong self-converse graph, necessary prime/7p residue conditions, and a sharp input-order bound attained by transitive inputs. This yields a complete finite inverse within the paired gate, and a separate order bound gives completeness for its two-input version. Neither is a theorem of unrestricted spectrum coverage.

The explicit16-vertex witness has H=37050461=7*5292923. Its full rows, anti-involution and trial-divisor certificate are in `three_cover_cyclic_witness.json`. Both the source and an independent walk-inclusion-exclusion checker directly count the actual output graph. The independent PASS report and replay receipts are included in `three-cover-review/`.

Run `python3 verify_three_cover_cyclic_gate.py` and its `-O` form. The command `--target M` runs the complete paired inverse and can be expensive. Compile `three-cover-review/check.cpp` as C++17 for the independent audit. Previously published helper dependencies and the pinned S195 source are bound in `THREE_COVER_DEPENDENCIES.json`; `THREE_COVER_MANIFEST.json` binds every new payload file.

The later complete larger catalog and general Gram/forest theorem are separate and are not part of this release. This packet itself makes no minimality or infinite-primality claim.
