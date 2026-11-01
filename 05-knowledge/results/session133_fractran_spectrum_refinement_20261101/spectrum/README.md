# Session 133: small-H atlas and sharp local obstruction certificates

**Correction to checkpoint 5811bb4:** the initial H15 packing-profile hostile
compared trailing-zero padding, not distinct polynomials. It is replaced by
a normalized H9 control. See [CORRECTION.md](CORRECTION.md). All inventories
and obstruction results are unchanged.

## Results and status

* **PROVED, from cited classical ingredients:** if `r` is the smallest positive
  integer with `3^r>B`, every induced-minimal obstruction to `H<=B` has at
  most `3r` vertices. The bound is sharp for every integer `B>=1`, attained
  by the ordered sum of `r` directed triangles.
* **FINITE-EXACT with proved all-order cutoff:** `H<=21` has exactly **56**
  minimal induced obstructions: 29 strong order-six, 13 strong order-seven,
  13 nonstrong order-eight, and one nonstrong order-nine class.
* The **whole obstruction bases**, not merely their sizes, agree at
  thresholds 5/7 and 19/21. In general such basis equality is equivalent
  to no attainable H between the two thresholds.
* The first complete terms of the obstruction-count sequence at odd thresholds
  `1,3,...,21` are `1,2,7,7,11,16,26,48,53,56,56`.
* A fresh standalone strong atlas through H21 supplies verified atoms for
  tournament sequence compilers, including self-converse counts.

The exclusions of exact H7 and H21 are **inherited**, not claimed new. The
conjecture that every other positive odd number occurs remains **OPEN**.
No literature-priority claim is made for the new-to-this-packet consequence.

## Read and reproduce

`OBSTRUCTION_PROOF.md` gives the general sharp theorem, complete bounded
classification, source citations, and mathematical boundaries. `ATLAS_PROOF.md`
records the independent small atlas and inherited results. All tournament
adjacencies and deletion counts occur in `basis21.json`; `basis7.json` is the
smaller boundary control. `basis_sequence.json` contains the complete terms.

From any working directory:

```
python3 /path/to/this/package/run_checks.py
```

Only the Python standard library is used. The harness regenerates every
threshold in normal and optimized mode, checks frozen outputs, and runs an
independent fixed-Hamilton-cycle cover plus literal permutation checks.
`check_basis_independent.py` imports no producer. Its 434 literal checks cover
every published B21 obstruction and every one-vertex deletion. The independent
strong cover has 512+16,384 presentations. Mutation controls reject an altered
H count and nonminimal singleton padding.

Degree-canonical codes in the basis differ intentionally from the full
permutation-minimum codes in the base atlas; actual adjacency is supplied in
both. Do not identify their numeric codes across conventions.
