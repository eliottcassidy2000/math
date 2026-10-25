# The first ordered birth-color ambiguity is at36

**FINITE-EXACT.** This supplements, without changing, [COLOR_WORD_BOUNDARY.md](COLOR_WORD_BOUNDARY.md). The order80 monochromatic hostile remains useful, but a completed exhaustive extension finds a much smaller failure.

Use the same precisely defined birth-color extension and orient every path from its smaller endpoint. The complete path search gives:

|order|Hamiltonian paths|endpoint plus ordered-color-word labels|largest label fibre|
|---|---:|---:|---:|
|33|669|669|1|
|34|4041|4041|1|
|35|17175|17175|1|
|36|12960|12957|2|

Combined with the main census through32 (and endpoint uniqueness at every earlier traceable order), this makes36 the first order where actual endpoints plus the entire ordered three-color word fail to identify the path. There are exactly three doubleton fibres at36. This is a statement about this coloring and this observer, not an abstract impossibility of small edge palettes.

One exact pair has endpoints6 and32:

`6,10,26,23,2,34,30,19,17,8,1,15,21,28,36,13,3,33,31,18,7,29,35,14,22,27,9,16,20,5,11,25,24,12,4,32`

`6,10,26,23,13,36,28,21,15,1,8,17,19,30,34,2,14,35,29,7,18,31,33,3,22,27,9,16,20,5,11,25,24,12,4,32`

Both are permutations of1..36, every adjacent sum is a square, and the ordered birth-color words agree entry by entry. Their actual edge sets differ. This is stronger than the earlier histogram collision: retaining the complete order of colors still loses which numerical endpoint carries each incidence.

The fast certificate [check_color_first_failure.py](check_color_first_failure.py) directly verifies these facts and emits the symmetric difference. It also checks the recorded census bookkeeping; it does not pretend that reading a count proves exhaustiveness. Exhaustiveness is reproduced separately by [enumerate_color_first_failure.py](enumerate_color_first_failure.py), using the primary DFS and writing [COLOR_FIRST_FAILURE_CENSUS.json](COLOR_FIRST_FAILURE_CENSUS.json). The complete replay took about four minutes in the research environment. Timings printed during replay are not mathematical data.

The finite birth-existence witnesses for33..36 come from the already frozen `birth_witnesses_25_80.json`, so the coloring of new edges does not assume an unverified infinite theorem. The independent session audit can check the pair immediately and the first-failure claim by its separate census.

The all-size survivors remain actual incidence reconstruction with endpoint and cotree coordinates, or the full numerical square-sum word with a marked starting vertex. A smaller unlabelled palette needs a proof of what endpoint ownership it retains; its apparent success through35 is not that proof.
