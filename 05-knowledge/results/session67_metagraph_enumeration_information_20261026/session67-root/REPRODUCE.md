# Reproduction and artifacts

Run from this directory (C++17 compiler and Python3 standard library only):

```
g++ -std=c++17 -O2 -Wall -Wextra -pedantic enumerate.cpp -o /tmp/session67-enumerate
/tmp/session67-enumerate 8 > /tmp/session67-census.json
cmp output/census.json /tmp/session67-census.json
python3 verify_exports.py
python3 -O verify_exports.py
```

To regenerate all exported tables, create a new directory OUT, then run:

```
/tmp/session67-enumerate 8 OUT 8 > OUT/census.json
```

The repository keeps lines_8.tsv compressed using Python gzip.compress(data,mtime=0).
The input to that compression is the complete 438374-row coloured endpoint/multiplicity ledger, not a sampled graph. Node tables give every ordinary class at each order1..8. Output is deterministic for the exact source. Gzip SHA256 is recorded in export_check.json; decompressed bytes can be compared without relying on compression-library version. The first row in every TSV is a header.

Verified runs: O2 and O3+DNDEBUG census bytes agree at all orders through8. There are no C++ assert checks to remove: explicit exceptions protect invariants in both builds. Python export verification uses explicit RuntimeError and yields byte-identical normal/-O output. No network, class-atlas, Burnside sequence or graph-isomorphism software is used in the enumerator.

The independent local auditor separately checks move guards against complete arc relations, classes against vertex-permutation oracles, all artifact ledgers, and the all-order theorem/certificate proof. Audit scope is reported separately from production finite enumeration. Tables are not proofs beyond their listed orders.
