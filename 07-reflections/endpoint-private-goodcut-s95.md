# Endpoint Private Good-Cut S95

The endpoint-transfer witness search found a triangular rank mechanism for
unmerged tournaments. The good-cut probe says this mechanism is not randomly
distributed across the tiling cube.

Every private odd child column found through `6->7` is good-cut pure:

```text
unmerged private columns:     [2, 4, 10, 29, 133]
unmerged width-zero columns:  [2, 4, 10, 29, 133]

merged private columns:       [2, 2, 8, 11, 74]
merged width-zero columns:    [2, 2, 8, 11, 74]
```

The merged case has a clean formal spine. A private odd column has odd column
sum, so the endpoint boundary law puts it in the odd child boundary. For
complement-merged tournaments, the odd child boundary is exactly the set of
self-complementary nodes. Thus every merged private pivot is forced to be SC.
The computation confirms this directly:

```text
merged private kind spectra:
2->3 {'SC': 2}
3->4 {'SC': 2}
4->5 {'SC': 8}
5->6 {'SC': 11}
6->7 {'SC': 74}
```

The surprising part initially looked like Morse purity. The follow-up
good-cut census closed that piece: THM-354 proves that good-cut count is
`n - scc(T)`, hence every tournament class is good-cut pure. The remaining
endpoint-specific signal is which SCC-defect levels provide private pivots.

## Pattern

At `6->7`, the private columns by good-cut bucket are:

```text
unmerged: {0: 1, 2: 5, 3: 4, 4: 19, 5: 47, 6: 57}
merged:   {0: 1, 2: 1, 4: 7, 6: 65}
```

The unmerged triangular witnesses spread across heights but still remain
height-pure class by class. The merged witnesses collapse harder into maximal
good-cut height, as if complement merging keeps only the self-complementary
barrier-crossing pivots.

Comparing against all SC child nodes gives a sharper constraint. In the tested
range, decomposable SC nodes are never the missing private pivots. At child
`n=7`, the non-private SC nodes by good-cut bucket are `{6: 14}`: all are
strongly connected.

## Working Picture

Endpoint insertion chooses a new endpoint signature. Good-cut count measures
strong-component defect. A private pivot seems to be a child whose
endpoint-deletion shadow has exactly one odd parent and whose SCC defect is
often very small, i.e. whose good-cut count is high.

So the rank proof might not require a global ordering on all classes. It may
require a recursive selection of Morse-pure pivots:

```text
parent class
  -> choose endpoint signature with maximal distinguishing cut profile
  -> obtain a child class with controlled SCC defect
  -> use it as a triangular pivot
```

Merged rows without private pivots may be the cases where two or more
self-complementary pivots with the same SCC defect collide after complement
identification.
Full rank survives because the collision matrix is still odd-Smith, but the
obvious triangular certificate disappears.

## Next Experiments

1. For each unmerged parent, choose the highest private good-cut pivot and test
   whether that choice is injective by a canonical invariant.
2. Compare private pivot heights with `H`, `F`, automorphism tax, and the
   bridge-exposure state `Q_T`.
3. Measure whether endpoint insertion has a bias toward decreasing SCC count
   by the largest possible amount.
4. Classify the strongly connected SC child nodes that fail to be private.
