# Session234: Hamiltonian-path counts 7 and 21 beyond tournaments

**Status: REFUTED (the proposed matrix condition); FINITE-EXACT (complete
loopless semicomplete census through order six); structural classification
under independent audit at this checkpoint.** Date: 2026-11-09.

The conjecture is false. A binary adjacency matrix can have exactly seven
or twenty-one directed Hamiltonian paths while its diagonal is entirely
zero and no mirrored off-diagonal pair consists of two zeros.

The distinction is **tournament versus semicomplete digraph**. A tournament
has exactly one arc between each pair of different vertices. A semicomplete
digraph has at least one, so it may contain reciprocal pairs. The classical
tournament omissions 7 and 21 do not survive that extra freedom. Self-loops
are irrelevant to vertex-distinct Hamiltonian paths.

## Smallest counterexamples

Rows of a four-vertex matrix with exactly seven paths:

```text
0111
0011
1001
1000
```

On vertices 0,1,2,3 its paths are exactly
`0123, 1203, 1230, 1302, 2013, 2301, 3012`.

Rows of a five-vertex matrix with exactly twenty-one paths:

```text
01111
00110
00011
10001
11000
```

These orders are minimal even among arbitrary digraphs: at most `n!`
directed Hamiltonian paths exist; three vertices allow at most six, and a
four-vertex graph missing any one directed nonloop arc has at most
`24-6=18`, whereas the complete bidirected four-vertex graph has 24.

## Complete finite catalogue

| order | labelled semicomplete matrices | H=7 | H=21 | H=7 classes | H=21 classes |
|---:|---:|---:|---:|---:|---:|
|1|1|0|0|0|0|
|2|3|0|0|0|0|
|3|27|0|0|0|0|
|4|729|72|0|3|0|
|5|59,049|720|1,760|6|16|
|6|14,348,907|6,480|32,640|9|48|

Exactly three strong seven-path cores occur, all of order four. There are
sixteen strong twenty-one-path cores of order five and sixteen of order
six. The latter include two classes with only **one** reciprocal pair.
Matrices, paths, reciprocal pairs, tournament resolutions and orbit sizes
are saved, not merely aggregate counts.

- [Finite census report and reproduction commands](census/REPORT.md)
- [Order 1–5 catalogue](census/results.json)
- [Order-six catalogue](census/classes6.json)
- [The two strong three-path cores](census/strong_h3.json)
- [Structural proof draft](STRUCTURE_DRAFT.md), independently being checked

**Important hostile control:** a directed triangle with one reverse arc
added still has three Hamiltonian paths. The added arc is unused. Thus
strong connectivity does not make arc addition strictly count-increasing.
Both three-path cores must be retained in a twenty-one-path product.

## Scope and next checkpoint

The structural reduction is through strongly connected components, whose
counts multiply. Every strong semicomplete digraph of order at least three
contains a spanning strong tournament. The inherited strong-tournament
floors therefore restrict a seven-path core to order four and a
twenty-one-path core to orders five or six. The next checkpoint will finish
the independent proof/census audit, all-order generation and exact counting
formulas. No completeness theorem for the remaining tournament spectrum is
claimed: the Chan–Pak conjecture remains open.

Inherited reference:
[THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
The similarly numbered elliptic-Jacobian file is unrelated.
