# Physical cyclic Gram kernels and rooted forests

For any fixed nonempty tournament V, the actual counts H(C3[T,U,V]) are a positive-definite bilinear form on the first |V|+1 path-cover coordinates. The reviewed proof gives its exact tridiagonal coefficient matrix, positive weighted-square decomposition, total nonnegativity, rooted-forest determinant, physical Cauchy-Schwarz inequality and exact rank. The one-vertex third block excludes every paired output 7p; the two-vertex block permits the separately reviewed 37050461 witness.

These are finite kernel statements. They do not assert total nonnegativity of arbitrary whole-module production matrices, nor of physical Gram arrays in arbitrary input order. Converse symmetry requires a reversing map for V. The separate content/gcd investigation is not included.

Run `python3 verify_cyclic_gram_forests.py` and its `-O` counterpart from this directory. The only computational import is the previously published `verify.py`. The verbatim independent PASS report, standalone C++ physical checker, Python incidence/forest checker and their receipts are in `gram-forest-review/`. Reproduction commands are in that report. The independent checks include 1331 small physical triples, 150 actual Gram entries, 64 larger-third-block controls, 4312 squared incidence minors and 18792 finite TN minors.

`GRAM_FOREST_DEPENDENCIES.json` binds the named inherited theorem and existing companion files. `GRAM_FOREST_MANIFEST.json` binds every new payload file. The proof and independent audit are preserved byte for byte.
