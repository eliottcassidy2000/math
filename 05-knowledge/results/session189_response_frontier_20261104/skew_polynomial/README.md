# Reproduce grouped skew certificates

PROOF.md states the general sign-definite excursion-polynomial criterion, its source-root/ambient restoration, four exact prime-nine examples, and scope boundaries. No software installation is needed.

The complete fixed inputs are included. Rebuild the two independent large polynomial computations:

    g++ -O3 -std=c++17 weighted.cpp -o /tmp/s189_skew_weighted
    /tmp/s189_skew_weighted < selected.txt > weighted.new.json
    cmp weighted.new.json selected.out
    g++ -O3 -std=c++17 runword.cpp -o /tmp/s189_skew_runword
    /tmp/s189_skew_runword < runword.txt > runword.new.json
    cmp runword.new.json runword.out

The physical weighted-subset engine peaks around600MB at17 vertices. All selected polynomial coefficients are exact. Rebuild the actual physical query controls:

    python3 physical.py
    g++ -O3 -std=c++17 count_physical.cpp -o /tmp/s189_skew_physical
    /tmp/s189_skew_physical < physical_input.txt > physical_values.txt
    /tmp/s189_skew_physical < compression_input.txt > compression_output.txt

Rebuild all nine complete small named response banks and their independent counts:

    g++ -O3 -std=c++17 bank.cpp -o /tmp/s189_skew_bank
    /tmp/s189_skew_bank small_bank.json < small_graphs.txt
    /tmp/s189_skew_weighted < small_weighted.txt > small_weighted.out
    g++ -O3 -std=c++17 predecessor.cpp -o /tmp/s189_skew_predecessor
    /tmp/s189_skew_predecessor < small_queries.txt > small_values.txt
    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json
    cmp normal.json receipt.json

The verifier independently checks quotient primeness, nonregular degree sequences, all endpoint/repeated-run metadata, absence of the old ordered witness, the complete positive coefficient certificates, all physical masks and graph expansions, the literal contrasts and ambient weights, all compression controls, and every complete small-bank coefficient. Normal/optimized runs give identical receipts.

RELEASE_FILES.json is the explicit publication allowlist. Exploratory candidate scans, the C1 probe, binaries, caches and intermediate logs are excluded. selected.json fixes the four quotients and their physical contracted cores, so reproducing these certificates does not require the larger incoming census.
