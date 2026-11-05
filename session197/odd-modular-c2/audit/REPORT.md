# Independent audit: odd-modulus C1 and C2 residue universality

**PASS, no correction to the frozen packet.** This is separate residue surjectivity for each class modulo every odd positive integer. It is not exact-value coverage, a cofinite theorem, a mixed even-modulus theorem or a small-order construction bound.

## 1. Snapshot and inherited inputs

Source: `session197/odd-modular-c2/`, with its unchanged eight-file manifest.

- `PROOF.md`: `1e011326969c4dc8a3f4f26ed11f673f451e51b153ae714a3fbf48760e03cb4b`.
- `verify.py`: `d1674999ecf81f9ff28ef280168bc17cef3488c67de2e31549a091050d4b71d3`.
- `core.py`: `dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042`.
- Both source receipts: `11e016df34e0b3acf27b7517668ec797427ffb70000f3a54d8c19a39f45e8a4f`.
- `WITNESSES.json`: `ed8164af79ec5842cbc134ef16f0dfb7b9d46c1be2d3fb4fe5fa356ef859d70b`.
- Copied ROOTED dependency: `ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e`.
- Copied BILINEAR_C2_INVERSE dependency: `5eab9deb2871819d8a1b5684c4987f2f4da807adc55063a30675d531c5ad3813`.

The physical gadget transfer and two-source C2 graft are inherited from those reviewed inputs. This review checks the new exact matrix argument, modular-to-positive-word realization and its physical use. The original files were preserved; replay used a copied subtree.

## 2. All-order proof review

The three large generators have determinant 16, and the singleton generator determinant 1. They are invertible modulo every odd m. A formal inverse in a finite matrix-group word can therefore be replaced by a nonnegative positive-generator power modulo that same full m. This is a finite terminating construction even without a useful order bound. It never applies a negative graph operation or asserts that a signed modular target vector itself is a physical source.

The exact central commutator is T12(3/4). The displayed commutator/conjugation chain produces all six elementary directions from T12(1), including the factor -4 that must be inverted away from 2. The rejected M5/M7 shortcut is indeed false; its first-row terms are nonzero. The corrected argument does not rely on it.

For m=3^a q, the chosen exponent makes the first relative shear identity modulo 3^a and T12(1) modulo q. All subsequent conjugations and commutators retain that identity on the first factor. Raising the -4 shear by its inverse exponent modulo q gives the required relative T32(1); the remaining directions then follow. This step uses only denominators that are units modulo the odd full modulus. Cases a=0 and q=1 are correctly treated as vacuous factors.

The two seed rays are actual strong start-one sources. Their exact readouts are 113+116*2^j and 1419+1472*2^j for C2, and 1+2^j and 9+10*2^j for C1. The elementary valuation proof gives the full order of 2 modulo each 3^a. Each pair of rays therefore covers the claimed complementary excluded residue classes, whose union is all residues. Positive j can always be chosen, preserving the strong-source boundary.

The explicit elementary seed reductions and target map have the right multiplication order. In particular the target vector (2r+3,35r+53,0) is primitive, has scalar readout r, and is the first column of the stated determinant-one product. Only its residue is requested. The relative operations adjust the q factor while preserving the already chosen entire 3-adic triple. This proves the desired CRT conclusion without assuming transitivity over a composite ring.

Every compiled positive word retains strongness and the unique-start property. For C2, the final M5 and the N graft with the nontrivial transitive two-vertex source give exactly 53D-3E+63U. For C1 no graft is needed: strongness supplies a Hamilton cycle, while cutting cycles at the unique-start root proves there can be only one. The proof therefore supplies a graph in each separate class, not merely in their union.

## 3. Independent numerical and physical controls

`check.py` imports no source implementation. Gadget matrices are derived from a direct binary-position sum and explicit reversed-complement response, matrix inverses use cofactors/adjugates rather than the source's Gaussian elimination, and modular orbits use only positive physical generators.

The checker verifies:

- The exact rational commutator and all six exact elementary directions, including the false-shortcut hostile.
- 122 relative CRT moduli through odd 255, with all 732 direction projections checked on both factors.
- Complete positive-generator orbits equal to every primitive vector modulo 5, 7, 11, 13, 25 and 35: respectively 124, 342, 1,330, 2,196, 15,500 and 42,408 vectors, 61,900 total.
- 28 explicit positive-power inverse controls, 204 exact seed reductions, 401 signed target-vector maps and 2,186 seed-ray controls through 3^7.

The graph builder independently puts the gadget apex first, starts from a literal triangle and constructs all seam/graft arcs directly. The separate `walks.cpp` counter uses inclusion-exclusion over closed and open walks in induced vertex subsets, not either source Hamilton-path DP. It confirms seven C1 source graphs and seven C2 final grafts through order 15. Their readouts include 345, 637, 4,363, 4,303, 8,659, 8,221 and 9,235 for C2. Optimized and assertion/undefined-behavior-sanitized C++ runs agree. The independent graphs and counts are finite controls, not premises of the arbitrary-modulus proof.

Source normal/-O replays reproduce both the frozen receipt and exact witness file. Independent normal/-O Python receipts agree as well.

## 4. Reproduction and hashes

From this audit directory:

```
python3 source/verify.py > source-normal.json
python3 -O source/verify.py > source-optimized.json
python3 check.py > independent-normal.json
python3 -O check.py > independent-optimized.json
g++ -std=c++20 -O2 walks.cpp -o /tmp/odd-modular-walks
/tmp/odd-modular-walks graphs.txt > walks-normal.txt
g++ -std=c++20 -O1 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined walks.cpp -o /tmp/odd-modular-walks-debug
/tmp/odd-modular-walks-debug graphs.txt > walks-debug.txt
```

The copied source writes its own VALIDATION and WITNESSES files. `check.py` writes only `graphs.txt`; the C++ checker reads that file. All dependencies are standard libraries and the complete copied source subtree.

- `check.py`: `c59ff85b055c75cd595aa12c95ce8d62ec69bfbb4dcd49a8d9a0639085392c0d`.
- `walks.cpp`: `7b1fc1015fbd8e09ab49465b33a48e4fe9a5800e43abf55d15afbffa08e92b59`.
- `graphs.txt`: `f0724e6a45d02ab8ad6c54a36e350a41f72d3245523dd91cd4efa47b9ce03b30`.
- Both independent Python receipts: `b9d791da6e81f7df17532e5e1f5eca9650081f1137603edb5982f50b58364ceb`.
- Both independent walk receipts: `9bc4b1f497eb9c288c8367b80cc0c281fdab94010bc80d5153a418bc670b84dd`.
- Both source replay receipts: `11e016df34e0b3acf27b7517668ec797427ffb70000f3a54d8c19a39f45e8a4f`.

No unresolved issue remains within the stated odd-modulus scope. Later proposed 3-adic group-lifting or dyadic extensions are not included in this binding.
