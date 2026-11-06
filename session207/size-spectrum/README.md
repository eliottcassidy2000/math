# Exact tournament order spectra, session 207

Read PROOF.md for the exact n<=8 supports, general SCC and vertex-extension
recurrences, prior attribution and limitations. REVIEW.md records independent
mathematical review. census.json contains full and strong supports and witnesses.

From this directory, reproduce the main certificate:

    g++ -O3 -std=c++17 verify.cpp -o /tmp/session207-spectrum
    /tmp/session207-spectrum > /tmp/session207-census.json
    cmp census.json /tmp/session207-census.json
    python3 check_witnesses.py

Reproduce independent audit outputs without changing archived files:

    g++ -O3 -std=c++17 audit/independent.cpp -o /tmp/session207-independent
    mkdir -p /tmp/session207-audit/session207-joint-audit
    (cd /tmp/session207-audit && /tmp/session207-independent > independent.out)
    cmp audit/independent.out /tmp/session207-audit/independent.out
    python3 audit/check_hostile.py > /tmp/session207-hostile.out
    cmp audit/hostile.out /tmp/session207-hostile.out

The independent program writes full and strong supports to spectrum1.txt through
spectrum8.txt under its session207-joint-audit output directory; compare these
sets with the corresponding census.json fields. SHA256SUMS freezes release files.
