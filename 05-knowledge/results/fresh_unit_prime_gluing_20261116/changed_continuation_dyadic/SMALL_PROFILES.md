# Exact low-deletion boundary of the two nine-count windows

This is an all-order boundary in the specified deletion layers, not an all-deletion impossibility. Notation and the positive bands are in [README.md](README.md).

## Exhaustive small-interface boundary

Borrow the deletion enumeration and weighted root contact identity
proved in the [deletion-45 replay](../ternary_deposit/SECOND_ROUND.md). The same
component bank covers all six deletion counts. Factor \(g\) as the
ordinal product of nonsingleton strong component counts, and distribute
singleton slots up to the Busch ceiling. A nonzero signed continuation
forces terminal root. A preceding ordinal prefix contributes a factor
\(a\mid g\) and \(a\mid H\), so enumerate all divisors \(a\) by scanning
terminal deletion count \(g/a\), then scale all fields. The divisor 45
is immaterial because no candidate H here is an odd multiple of 45.
Every tested divisor is realizable by prefix factors, and the program
counts all masks into each representative irrespective of automorphisms.

The replay yields 1175 terminal profiles, maximum four of nine
counts in one profile. The only aligned three-consecutive block
of the prescribed groups (indices 0–2, 3–5, 6–8) is indices 6–8
at q=163 and profile (1107,1122,15). No nonzero-signed profile
covers all nine in either window. This does not claim the exclusion
of every aligned triple for signed-zero continuations.

### Arbitrary signed-zero continuations cannot cover all nine either

For a prime target \(H\) greater than \(g\), a root with small
\(g\) is inside the sole nonsingleton strong component and all
ordinal predecessors and successors are singletons. In a strong
rooted tournament its \(E-O\) may or may not vanish; to safely screen
for signed-zero outputs include all strong terminal source sums
\(U=E+O\), whatever their signed difference. Adding \(k\ge1\)
ordinal sink singletons gives \(E=O=2^{k-1}U\) and preserves \(g\),
while adding no sink leaves signed-zero \(E=O=U/2\) when \(E=O\)
already. Thus any common signed-zero continuation has a common
odd part of \(U\) among the strong component possibilities.

For \(q=163\), the prime counts 1471, 1481, 1483 already give empty
intersection at each of the six \(g\)'s. For \(q=79\), prime counts
719 and 727 leave only \(g=3\) (odd part 1137) or \(g=9\) (845,
1051). But with these \(g\)'s any \(H\) coprime to 3 also forces
all prefix/suffix components to be singletons. To see this note
that the root component contains the root, and the product of
all other component counts divides both \(g\) and \(H\).
Then 713,715,719,721,725,727 together leave neither possibility.

Replay:
    
    g++ -O3 -std=c++17 scan_small_profiles.cpp -o /tmp/gscan
    /tmp/gscan >profiles.txt 2>witnesses.txt
    python check_profiles.py
