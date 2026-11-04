# Reproduce the incoming deletion check

This packet separates a reproducible incoming finite census from independent checks on its three reported collision pairs. Read PROOF.md for exact observer distinctions and scope.

The exact source files are source.c, source.out, HYP-9455.md. source_metadata.json records their Git blobs and the targeted later-source check. incoming_preregistration.md is a separately identified later preregistration, not a theorem used here.

Full replay with standard installed tools:

    mkdir -p data
    gcc -O3 -std=c11 source.c -o /tmp/s189_deletion_census
    /tmp/s189_deletion_census 9 data > reproduced.out
    python3 prepare.py
    g++ -O3 -std=c++17 bank.cpp -o /tmp/s189_deletion_bank
    /tmp/s189_deletion_bank bank.json < graphs.txt
    g++ -O3 -std=c++17 literal_cards.cpp -o /tmp/s189_deletion_cards
    /tmp/s189_deletion_cards < graphs.txt > literal_cards.json
    python3 prepare_independent.py
    g++ -O3 -std=c++17 predecessor.cpp -o /tmp/s189_deletion_predecessor
    /tmp/s189_deletion_predecessor < independent_queries.txt > independent_values.txt
    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json
    cmp normal.json receipt.json

The large full census exports, binaries, and generated independent query input are intentionally excluded from RELEASE_FILES.json. prepare.py reconstructs the exact six source rows and both converse choices from the export. All twelve full response vectors, independent actual counts, and independent deleted-card records are included. No installation is required. Counts through eleven augmented vertices fit in signed 64-bit integers.

The source reproduction is not an independent proof of the complete census or its ranks. The six-class exclusions have independent whole/card canonicalizations and independently recounted complete named response data. Normal and optimized Python receipts agree.
