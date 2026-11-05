# Independent full-census audit

From this audit directory run:

    g++ -O3 -std=c++17 independent_census.cpp -o /tmp/source_census
    /tmp/source_census 50001 > /tmp/source_small.json
    /tmp/source_census 1000001 > /tmp/source_full.json
    cmp /tmp/source_small.json CAP50001.json
    cmp /tmp/source_full.json CAP1000001.json

The full independent calculation was also rebuilt with -O2 -DNDEBUG and produced identical bytes. This code uses exceptions, not assert. Its maximum declared cap is1,000,001, and the20-bit coordinate packing is checked on every input.

The exact universe is every positive-coordinate physical source reachable from the edge with at most three J letters and G<=cap. For each cost layer it closes O/Z from the J successors of the previous minimum-cost layer. A source reached earlier has a smaller or equal cost, so keeping only its earliest layer loses no optimum. Every retained source lies on e=0, p=0, or f=0; a canonical two-bit face tag and three20-bit coordinates encode it injectively. All transitions are actual source operations, and all accepted states stay within the monotone G cutoff. The singleton scalar3 is added separately at cost0.

This proves completeness for the declared finite universe. It does not extrapolate the stable gap list to higher targets. The all-order compiler, graph rebase and density arguments were separately reviewed from the proof.
