# Endpoint-590 nine-for-nine exchange packet

**FINITE-EXACT relative to the fixed THM-4313 carrier, the labelled
30-speed rank-8/rank-9 mask model, and THM-4314's frozen endpoint-590 row
layer. This is neither a physical entry nor a proof of LRC(14).**

## Baseline and complete response boundary

The unchanged carrier has 3,925 masks (3,818 rank eight and 107 rank nine),
retains all 421 protected joint masks, and has ordered FNV
`a0d08a38c10bdab7`. Independent `-O2` and `-O3 -DNDEBUG` scans of all

```text
13 * binom(30,9) = 185,992,950
```

row-body cases on THM-4314's endpoint-590 layer are byte-identical. There are
exactly 100 failures, all at `(210,590)`, with ordered failure FNV
`8d19cba1e86e53b5`. The complete active absent rank-8/rank-9 response universe
induces 14,368 distinct nonempty signatures on these 100 bodies.

The nine selected additions, all rank nine, are

```text
20490236
22045017
29224016
0a439108
0b220096
0120403f
12844116
10686016
084a6016
```

Their ordered FNV is `d1cf49e4b811b958`. Direct arithmetic gives 136 response
incidences and hit distribution `1:69, 2:26, 3:5`, hence covers every one of
the 100 failures.

A complete depth-eight search over the lossless quotient of 1,165
inclusion-maximal signatures visits 7,163,197 nodes and returns `UNSAT` in
byte-identical `-O2` and `-O3` runs. The dominance reductions, pivot branches,
sum bound, and denominator-three integer dual bound are all replacement-safe
necessary conditions. Thus, within the complete active rank-8/rank-9 response
universe,

```text
minimum number of additions covering the frozen 100 failures = 9.
```

## Deletion quotient and simultaneous exchange

On the fixed carrier augmented by those nine additions, the exact protected
singleton-deletion quotient over the 480 inherited rows plus 13 endpoint-590
rows has:

```text
private nonjoint obligations    1,600   FNV e862a639d9536826
protected old masks               425   FNV 470279b13b453834
singleton-safe old masks         3,079   FNV 89be292ce50c2831
```

These are singleton-deletion statements only; they do not by themselves imply
that several safe masks may be deleted simultaneously. The chosen nine old
rank-eight nonjoint deletions are

```text
06021829
23222801
12444083
20827018
29c04082
02916180
13c00881
070c4840
2380408a
```

Their ordered FNV is `3546eb56552b4cde`; none belongs to the protected
421-mask joint deck. A separate direct simultaneous replay, not an inference
from the singleton quotient, reconstructs the exchanged 3,925-mask carrier:

```text
rank eight                 3,809
rank nine                    116
joint masks retained         421
carrier FNV     eeae5518d84ccac5
rows                         493
row FNV         1fef91ec25d074e5
body tests         7,053,424,950
failures                        0
pair-ledger FNV 1092fd57a8581a34
```

The full 493-row pair ledger and empty failure ledger are byte-identical under
`-O2` and `-O3 -DNDEBUG`. This certifies only the fixed simultaneous nine-add,
nine-delete exchange on the stated 493-row target.

## Typed frontier consumer

The typed consumer removes exactly THM-4314's 13 endpoint-590 rows from its
20,547-row residual and adds them to its 2,100-row typed union. Normal Python
and `python -O` outputs are byte-identical:

```text
new typed union       2,113   FNV c806cce6b836fdff
new residual         20,534   FNV 11285b5a49f4150d
next endpoint           589
rows at endpoint 589      28   FNV 5d9429c9f9971322
```

The exact next-frontier SHA-256 is
`005306a28c756a93862fd1745414fc058032f4d99e32086c079c29607bfed0c6`.

## Scope boundary

Every statement above is finite and exact for the frozen files and model.
The response optimum concerns only covering the 100 endpoint-590 failures.
The protected quotient distinguishes arbitrary single deletion from the
separately replayed selected simultaneous deletion. The packet proves no
terminating descent, physical lonely-runner entry, or LRC(14).
