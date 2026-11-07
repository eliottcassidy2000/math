# Hamiltonian-cycle spectra, by order and by value

The all-order cycle-count spectrum is complete: every nonnegative integer
occurs. Every positive value has a strong realization of logarithmic order,
and every value at least two already occurs in the sharp feedback vertex
class two. This packet recovers a binary endpoint supply and a forced-entry
bijection, complementing the Session212 binary-comb packet's gate-free
construction and marked doubling grammar. No global path-spectrum
completeness or literature-priority claim is made.

Read `PROOF.md` for the mathematical statements and their inherited sources.
Exact complete spectra through order nine are freshly recomputed here.
Order ten is explicitly **inherited**, extracted from the pinned THM-4889
full census, with an independently counted maximum witness.

## Reproduce

Run from this directory:

    python3 generate_classes.py --output tournaments8.txt
    g++ -std=c++17 -O3 -Wall -Wextra -Werror census.cpp -o /tmp/cycle-spectrum
    for n in 3 4 5 6 7 8; do /tmp/cycle-spectrum --fixedcycle "$n" > "fixedcycle$n.json"; done
    /tmp/cycle-spectrum --extend8 tournaments8.txt > spectrum9.json
    python3 verify.py > /tmp/cycle-normal.json
    python3 -O verify.py > /tmp/cycle-optimized.json
    cmp /tmp/cycle-normal.json /tmp/cycle-optimized.json
    cmp /tmp/cycle-normal.json receipt-normal.json

To independently replay inherited-data extraction, point to a repository
containing the pinned commit, not merely an exported snapshot:

    python3 extract_inherited.py --repo /path/to/math > /tmp/inherited-receipt.json
    cmp /tmp/inherited-receipt.json inherited-receipt.json

The source files are read by `git show`, without changing that checkout.
The extractor checks every source against its committed SHA256SUMS.
Stored inherited summaries suffice for all other reproduction commands.

A debug C++ build with `-O0 -fsanitize=undefined` was also run against the
complete order-eight fixed-cycle and order-nine extension universes, and
both outputs byte-matched the `-O3` outputs. There are no external software
or data prerequisites besides Python 3 and a C++17 compiler; git and the
pinned source commit are needed only for optional inherited-data replay.

The raw representative/presentation multiplicities of `census.cpp` must
not be called isomorphism-class counts. The `inherited_summary.json`
multiplicities are class counts imported from THM-4889. These different
universes agree on supports, not on frequencies.
