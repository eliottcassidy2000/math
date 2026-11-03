# Independent order-ten endpoint-block lower-bound audit

For each r, the first 2r old vertices form a fixed endpoint set E. The banked exterior contexts all touch every vertex of E. An independent endpoint-pair/path exact-cover census produces the following block catalogues and finite-field ranks:

| r | Columns | Rows | Rank mod 101 and 103 | det mod 101 | det mod 103 |
|---|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 4 | 4 |
| 2 | 6 | 6 | 6 | 88 | 10 |
| 3 | 60 | 60 | 60 | 12 | 47 |
| 4 | 735 | 679 | 679 | rectangular | rectangular |
| 5 | 945 | 945 | 945 | 21 | 18 |

For each block with more than one row, deletion or duplication reduces rank by exactly one over both fields. Size-one deletion yields rank zero; duplication is inapplicable. The r=4 result is a **679-row lower-bound certificate**, not an upper bound on the possible rank or an intrinsic obstruction. The source search stopped at 679 after its specified search; this says nothing conclusive about unsearched contexts.

For full-support columns a context touching E is zero unless the column's endpoint set contains E: any E vertex outside that set has old degree two, giving selected degree at least three. Independent rows in each representative block have a nonzero square minor; relabelling supplies the same size minor for each E, with nonzero sign gauges. Order these sets by inclusion. The selected rows and columns form a block-triangular minor with nonzero diagonal blocks, proving at least

`C(10,2) + 6 C(10,4) + 60 C(10,6) + 679 C(10,8) + 945 C(10,10) = 45,405`

independent no-isolate directions. The order-nine result and its lifting yield all 299,701 signatures with an isolate. Hence the all-fixed-source-oblivious-attachment dimension is at least 345,106 and at most the established pooled-parity bound 347,626. For normalized laws the corresponding bounds are 345,105 and 347,625. There is no S/B-only or nonlinear claim.

More generally, let d=(1,0,1,3,9,40,120,546,2233,9612). For n≥10, the certified lower bound is `sum_{k=0}^9 d_k C(n,k) + 45405 C(n,10)`, by lifting each observed support of size at most ten. This inherits the order-nine certificate, the lifting proof, and Session145's coefficient and upper-bound theorems. Exterior coefficient rows are rational linear combinations of actual fixed-cut evaluations, rather than single evaluations.

The verifier enumerates endpoint-pair paths over ordered interior subsets and constructs no-isolate forests by an exact-cover recursion. It independently checks all representative masks and signs, context format, and saturation, then uses exact finite-field elimination. Hostile inputs removing all edges from a selected endpoint and introducing loops, old-old edges, reversed duplicate edges, bad bounds, truncation or trailing data are rejected. Python normal and -O hostile receipts match.

Reproduce:

```sh
cat block5.contexts.part00 block5.contexts.part01 > block5.contexts
g++ -O3 -std=c++17 -o verify10 verify10.cpp
for r in 1 2 3 4 5; do ./verify10 "$r" "block$r.contexts" "block$r.reps"; done
python3 -B hostiles.py > hostiles.normal.log
python3 -B -O hostiles.py > hostiles.optimized.log
cmp hostiles.normal.log hostiles.optimized.log
sha256sum -c SHA256SUMS
```

If the 109,962-byte block-5 bank is distributed in ordered parts, concatenate the parts to `block5.contexts` and check the pinned SHA-256 before replay.
