# Exact composite-tree census through order 15

Input: aggregate prime cycle index recovered internally; no composite representative census beyond the independent oracle through six. FINITE-EXACT recurrence outputs for n>6.

Reproduce the underlying joint counts and depth rows with `python3 check_marked_composites.py 15`. This emits exact JSON; the tables below display its marginals.

## Root families: primes, strong composites, and ordered sums

The singleton is separated. For n>=2, the three displayed columns partition all ordinary tournament classes. These exact coefficients follow from S=1−1/B and the recovered prime index; they are not a separate class census.

| n | modular primes | strong composites | nonstrong ordered sums | total |
|---:|---:|---:|---:|---:|
| 2 | 0 | 0 | 1 | 1 |
| 3 | 1 | 0 | 1 | 2 |
| 4 | 0 | 1 | 3 | 4 |
| 5 | 3 | 3 | 6 | 12 |
| 6 | 15 | 20 | 21 | 56 |
| 7 | 197 | 156 | 103 | 456 |
| 8 | 4008 | 2000 | 872 | 6880 |
| 9 | 136655 | 41478 | 13403 | 191536 |
| 10 | 7906161 | 1449788 | 377107 | 9733056 |
| 11 | 796876429 | 87588161 | 19288658 | 903753248 |
| 12 | 142980435791 | 9329713944 | 1798161433 | 154108311168 |
| 13 | 46453867277461 | 1780914985832 | 307332423619 | 48542114686912 |
| 14 | 27687212017612715 | 617279770545341 | 96931930964248 | 28401423719122304 |
| 15 | 30571668704134867685 | 392578842568016044 | 56754613652283119 | 31021002160355166848 |

## Prime-node count distribution

Each row lists counts for b=0,1,...,floor((n−1)/2).

| n | counts by number b of prime nodes |
|---|---|
| 1 | 1 |
| 2 | 1 |
| 3 | 1, 1 |
| 4 | 1, 3 |
| 5 | 1, 10, 1 |
| 6 | 1, 47, 8 |
| 7 | 1, 405, 48, 2 |
| 8 | 1, 6511, 345, 23 |
| 9 | 1, 187668, 3667, 195, 5 |
| 10 | 1, 9666765, 64465, 1756, 69 |
| 11 | 1, 901751307, 1980348, 20850, 731, 11 |
| 12 | 1, 153998625848, 109294727, 382649, 7735, 208 |
| 13 | 1, 48531130448192, 10972006600, 12130083, 99350, 2658, 28 |
| 14 | 1, 28399406711848458, 2016304961705, 700423488, 1856046, 31970, 636 |
| 15 | 1, 31020320111776590293, 681974061371154, 74457767736, 58989912, 438202, 9476, 74 |

## Linear-node count distribution

Each row lists counts for l=0,1,... through the last nonzero coefficient. These are canonical linear nodes, not SCC count or a binary parenthesization.

| n | counts by number l of linear nodes |
|---|---|
| 1 | 1 |
| 2 | 0, 1 |
| 3 | 1, 1 |
| 4 | 0, 4 |
| 5 | 4, 5, 3 |
| 6 | 15, 31, 7, 3 |
| 7 | 213, 155, 78, 8, 2 |
| 8 | 4109, 2085, 562, 119, 5 |
| 9 | 138298, 43458, 8574, 1101, 105 |
| 10 | 7943040, 1569721, 199357, 19653, 1228, 57 |
| 11 | 798257904, 97044988, 7899626, 522464, 27436, 806, 24 |
| 12 | 143067771384, 10483059038, 533352249, 23240267, 863693, 24216, 316, 5 |
| 13 | 46463443937593, 2014271189238, 62612043115, 1742252014, 44313819, 937457, 13617, 59 |
| 14 | 27689071668963044, 699113086965720, 13010064149039, 225079354503, 3761708268, 57304598, 672342, 4790 |
| 15 | 30572319159590447457, 443771449998227635, 4859925607164970, 51076721664076, 542745244484, 5641111362, 50994844, 310958, 1062 |

## Exact maximum prime-depth distribution

Each row lists counts for h=0,1,...,floor((n−1)/2).

| n | counts by exact prime-depth h |
|---|---|
| 1 | 1 |
| 2 | 1 |
| 3 | 1, 1 |
| 4 | 1, 3 |
| 5 | 1, 10, 1 |
| 6 | 1, 48, 7 |
| 7 | 1, 410, 44, 1 |
| 8 | 1, 6534, 334, 11 |
| 9 | 1, 187794, 3645, 95, 1 |
| 10 | 1, 9667810, 64373, 857, 15 |
| 11 | 1, 901766547, 1976911, 9626, 162, 1 |
| 12 | 1, 153999034411, 109116915, 158146, 1676, 19 |
| 13 | 1, 48531150653310, 10959611006, 4402732, 19617, 245, 1 |
| 14 | 1, 28399408557280212, 2014935032768, 226499742, 306703, 2855, 23 |
| 15 | 1, 31020320423535238151, 681714888987901, 21923159574, 7745707, 35169, 344, 1 |

The full joint (n,b,l) histogram is in the JSON emitted by `python3 check_marked_composites.py 15`. The checker internally verifies the full marked power-sum cycle index, not only these scalar marginals.
