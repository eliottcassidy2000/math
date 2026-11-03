# Three different interpretations of square line sums

## 1. Nine square entries and eight equal line sums that are themselves squares

This interpretation is impossible for a nonzero 3 by 3 magic square of integer squares. In any 3 by 3 magic square with common line sum T and center entry c, add the middle row, middle column, and both diagonals. The resulting sum is 4T. Each noncentral cell appears once and the center four times, so it is also 3T + 3c. Hence T = 3c. If c = e², then T = 3e². For nonzero integer e, the exponent of 3 in 3e² is odd; therefore T cannot be an integer square. This is an elementary obstruction, not a result about the unresolved distinct-square magic-square problem with unrestricted T.

## 2. Nine distinct positive square entries, with eight square line sums allowed to differ

An explicit certificate is the grid

    2655²  2976²   320²
    2624²  2520²  1665²
    1440²   895²  3624²

Equivalently:

    7049025  8856576    102400
    6885376  6350400   2772225
    2073600   801025  13133376

Each of the three rows and each of the three columns sums to 16008001 = 4001². The main diagonal sums to 26532801 = 5151². The other diagonal sums to 8526400 = 2920². All nine entries are positive and distinct. Thus this weaker interpretation is solved affirmatively by direct arithmetic. It is NOT a magic square and is NOT a solution to the nine-distinct-squares magic-square problem or any prize requiring equal line sums.

### Construction and search scope

For integers a,b,c,d put N=a²+b²+c²+d² and

    M = [[a²+b²-c²-d², 2(bc-ad),       2(bd+ac)],
         [2(bc+ad),     a²-b²+c²-d², 2(cd-ab)],
         [2(bd-ac),     2(cd+ab),     a²-b²-c²+d²]].

The quaternion identity MMᵀ=N²I ensures that squaring the nine entries gives square row and column sums N². Search the six permutations selecting one entry in every row and column. If two such transversal sums are squares and their permutations agree in exactly one position, row and column permutations can place them on the two diagonals.

search.cpp enumerated primitive ordered nonnegative quaternions 1 <= a <= 120, 0 <= d <= c <= b <= a, stopping at the first witness. It rejected zero or repeated absolute matrix entries and tested the six squared transversal sums exactly (floating-point square root is only an initial estimate, followed by integer correction and integer equality). It found (a,b,c,d)=(48,32,23,12), N=4001, with transpose-selection permutations (0,2,1) and (1,2,0). Permuting columns to (0,2,1) and taking absolute values gives the displayed root grid.

The early-stopped run visited 233046 primitive ordered quaternions, of which 160271 passed the distinct nonzero-entry test, and encountered 1229 square transversals before stopping. This is not exhaustion of the a<=120 box, not a proof of smallest entries, and not a general nonexistence search. No global classification or originality claim is made. The witness is self-contained and requires no reliance on literature status.

Verification:

    python3 certificate.py
    g++ -O3 -std=c++17 search.cpp -o search
    ./search 120

## 3. Arbitrary integer entries and a common square line sum

If entries need not themselves be squares, scale the usual Lo Shu square by 15:

    120   15   90
     45   75  105
     60  135   30

Its eight line sums are all 225=15² and its entries are distinct positive integers. Its entries are not squares, so it also does not meet the nine-distinct-square-entry requirement.
