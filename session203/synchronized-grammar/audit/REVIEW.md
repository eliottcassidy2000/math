# Independent audit of eight-test synchronization

**PASS: the stated all-order equivalence is proved.** Existence of a physical seed pair satisfying it, and occupancy of any resulting common-base orbit, remain open.

Reviewed proof SHA256: `589addeec077c4793176221992d98732d8a00a164b9d6c2c8ac5ea1b52175248`. Producer verifier SHA256: `e001cbc101f3b242fa795186ecdbe438143bddfb001702fd323e88171b322423`. All six frozen producer-manifest entries and their manifest are preserved under `reviewed-inputs/`.

## Exact hypotheses and conclusion

The two inputs are actual nontrivial rooted tournaments with positive starts. Their residues are fixed as r=1 and r=2, and their admissible shifts j are fixed, with s=r+3j. Their adjusted bases agree initially. Every continuation is the same chronological O/Z/J word applied to both inputs, using the physical root change under J.

Under these hypotheses, the following are equivalent:

- Equality of adjusted bases after every common finite word.
- Equality after Z, J, OZ, OJ, ZZ, ZJ, JJ and JOZ.
- Equality of h,e,L,R,d,C,b,T, and s_1=s_2+2.

The scope does not include changing the shifts during continuation, choosing different words for the two sources, or dropping initial equality. The nontrivial-input hypothesis excludes the real singleton Z exception. Since all moves add a vertex, later inputs stay nontrivial.

## Algebraic proof review

Put G=h+L. The exact identity

    X=G-(s+2r)/3

follows from B=h+L-s and s=r+3j. The published nontrivial-source transfers preserve s. Therefore, after subtracting initial adjusted-base equality, synchronization along a word w is exactly

    (l M_w-l) delta=0,  l=(1,0,0,1,0,0,0,0,0).

The chronological convention is correct: the first physical letter acts first, hence is the rightmost matrix in the product. Independent basis-vector evaluation reproduces every displayed entry of the eight-row matrix.

After removing the zero s column and preserving column order (h,e,L,R,d,C,b,T), the determinant is exactly +32. Exact rational elimination confirms rank eight. Consequently the full nine-column kernel is precisely the span of the start-coordinate unit vector. This is a characteristic-zero statement; no invertibility modulo two is asserted or needed.

Thus the eight tests force equality of every non-s coordinate. Initial base equality then gives s_1+2=s_2+4, hence the claimed start difference. Conversely, none of the eight other transfer coordinates depends on s, and each transfer fixes the start-coordinate unit vector. Their equality and the start difference are therefore preserved by arbitrary words. This proves the all-order converse, not merely a finite-word observation.

For width eight, s_1=1+3j_1 and s_2=2+3j_2 give j_1=j_2+1. The admissible pairs are exactly (j_1,j_2)=(1,0),...,(8,7), with starts (4,2),(7,5),...,(25,23). The allowed numerical starts do not establish that any corresponding physical pair exists.

A nonzero difference in any non-s coordinate is detected by at least one of the listed words, given initial equality. This follows from the nonsingular minor. The report does not claim a minimal number of tests for the physical source class.

## Published physical dependency

The frozen nine-coordinate transfer proof has SHA256 `8e78951f9dd15decff0d23cfb3ebe64f358b3c9d4561b269b527a1e9aa902168`. It is the published [rooted-cover closure proof](https://github.com/eliottcassidy2000/math/blob/5067f5a9e43237ab3013399b3636dff92096d886/session203/rooted-cover-closure/PROOF.md), whose exact commit, path and hash are retained in the input provenance.

It supplies actual graph operations, all-input path-cover meanings, the singleton exception, and fixed-start transfers for nontrivial inputs. The synchronization argument uses only its nontrivial linear transfers. The all-word theorem is valid at arbitrary source order; it does not rely on a graph census or on the reported failure of any finite search to find a seed pair.

The constructive-consequence paragraph separately inherits the actual four-core readout H=3B+2s. Its reviewed published source is frozen as `READOUT_DEPENDENCY.md`, SHA256 `211623af96e90d516fbce7598450642878f81fc2892446dee5173b55fb41b099`, from [the published criterion dependency](https://github.com/eliottcassidy2000/math/blob/e5570ebddf67175234d0504af71da481747380ec/session203/staggered-profile-criterion/dependencies/FREE_ADDITIVE_PROOF.md). This ingredient is not needed for the eight-test equivalence itself.

The final frozen proof removed the ancillary discovery search. That search, its census and the reported 907,420 states are not part of this audited packet or its approved conclusions.

## Independent checks

`check.py` imports no symbolic algebra package and no producer code. It obtains each observation row by applying the transfers to the nine standard basis vectors, then uses exact Fraction elimination for the determinant and rank. It also checks:

- All 1,093 words of lengths 0 through 6 annihilate the start-only difference and preserve the difference 2 in starts.
- Every individual non-s coordinate has a separator among the eight listed words, using synthetic initially synchronized deltas. These are explicitly algebraic controls, not claimed physical sources.
- All eight admissible width-eight shift pairs and no others.

`check_physical.cpp` uses the independent append-at-end path recurrence and exact cover convolutions previously authored for the paired-carrier audit, not the producer's symbolic calculation. It checks every marked nontrivial labeled tournament through order 4: 284 inputs, including 166 with positive starts. All eight actual continuations are constructed and recounted, for 2,272 literal output graphs through order 7. Every full nine-coordinate transfer and every displayed observation row matches. The singleton Z start change 1->0 is a separate hostile control.

Normal/optimized Python receipts agree, and optimized/debug physical C++ receipts agree. The producer's SymPy verifier was also replayed normally and under Python -O; both outputs match its frozen receipt exactly. No validation gate depends on removable assertions.

## Reproduction

From this audit directory:

    python3 check.py > NORMAL.json
    python3 -O check.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json
    g++ -O3 -std=c++17 check_physical.cpp -o physical
    ./physical > PHYSICAL_NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 check_physical.cpp -o physical_debug
    ./physical_debug > PHYSICAL_DEBUG.json
    cmp PHYSICAL_NORMAL.json PHYSICAL_DEBUG.json
    python3 reviewed-inputs/verify.py > PRODUCER_NORMAL.json
    python3 -O reviewed-inputs/verify.py > PRODUCER_OPTIMIZED.json
    cmp PRODUCER_NORMAL.json PRODUCER_OPTIMIZED.json
    cmp PRODUCER_NORMAL.json reviewed-inputs/NORMAL.json

Independent symbolic source SHA256: `5ae2bee1d3b6d79ccfb310eb6ba7c6f7d70b53c1cd62dfe07b16f212bced45dc`.

Independent symbolic receipt SHA256: `43c0fcce599b0f53919e46bb94f943257d09416475196dbefbbe12de6caac35c`.

Independent physical source SHA256: `5e08d6a58d41d3a6bbe043ab82622dd746b48e8f99be601fbcd07314d16d5d2b`.

Independent physical receipt SHA256: `ef364e9cdea4223b7e2c5adf7a024a725d232d33e46c31dd274db9d55c949e0a`.

The manifest freezes reviewed inputs, independent sources and receipts, excluding local binaries. This is a proved synchronization criterion, not a physical seed-pair construction, an all-target coverage theorem or a threshold result.
