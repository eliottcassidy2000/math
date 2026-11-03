# Literal 6325 return, then the saved 18979 child

**FINITE-EXACT.** Start from `digit_cycle_ternary_20261103/return2107.json`,
not a new representative. The thirteen old names are

    [1,2,4,5,6,7,8,13,14,21,22,31,32].

The actual root is 31, spare is 6, and the root's unique in-port is 21.
The initial profile is (6325,2107,2107,6325), golden phase (6,3,1).

The one/two-start parent grammars were scanned through order twelve, using
an exact 64-bit C++ endpoint DP. There is one target in this bounded bank,
with parent order twelve and two starting paths. Its newborn has profile
**(6325,891,893,6327)** and order thirteen. A separate Python DP checks
the positive parent and child counts. No classification at larger orders
or among the exceptional ready families is claimed.

Fix both the actual root and its actual unique in-port, and compare all
11! identifications of this target. The minimum arc discrepancy is **ten**.
The branch-and-bound search is repeated without its degree lower bound,
giving the same minimum. It is a minimum in this one target class, not over
all ready sources or all native programs.

The degree bound is elementary. For each unmatched source vertex i and
candidate target j, let z_ij count discrepancies to already assigned names,
and let d_i,e_j be their outdegrees inside the remaining sets. Any completion
has future cost at least sum_i z_i,p(i) + (1/2)sum_i |d_i-e_p(i)|, since an
unmatched reversed arc changes at most two degrees. Taking independent row
minima of 2z+|d-e| is a valid lower bound on twice that cost. The second
search prunes only by accumulated discrepancies, checking this shortcut.

An ordering of the ten repairs keeps every completed arc-repair endpoint
strong. Each reversal has a supplied Hamilton path containing its directed
arc with the second endpoint nonterminal. The one-helper, two-native-move
compiler then reverses exactly that arc. All actual paths and intermediate
rows are saved. The cost is:

* twenty native complements;
* ten helper births and ten helper deletions;
* zero old deletions;
* peak order fourteen, returned order thirteen.

Every old name survives. Neither intermediate H nor phase is asserted constant.
Temporary helper stages are not asserted strong. An independent endpoint
algorithm recounts the returned profile. The retained root's golden observer
is (6325,57732,57759), phase **(4,2,2)**. The full eighteen-polynomial
root/spare observer and all five Q4 presentations are retained.

The two-birth reset at the same old root gives

    18979 = 3*6325+4.

Its new names are 43 and 44, its next root is 43, and spare 6 remains.
The fifteen-vertex child has profile (18979,6325,6325,18979), golden observer
(18979,634781,634808) and phase **(5,6,3)**. The source rows, actual old names,
and all Q4 phases are in `return6325.json`.

[O6](PROOFS.md) proves why merely cloning this repaired root is not a closed
return family: the scalar carry grows as 1342(2^m-2), and every clone mark
has ready defect at least 27 once m>=2. The next construction must repair
both coefficient data and readiness, or change the source family.
