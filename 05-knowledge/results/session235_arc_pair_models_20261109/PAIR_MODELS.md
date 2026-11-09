# The pair-state dictionary

**Status: PROVED (elementary representation facts and expectations).** This
file distinguishes a binary *encoding* from the adjacency matrix decoded
from it. All off-diagonal statements refer to labelled vertices. Diagonal
bits do not affect vertex-distinct Hamiltonian paths.

## 1. Four ordinary pair states and the tournament boundary

For i<j write the raw pair as (a,b)=(A[i,j],A[j,i]). Ordinary adjacency
semantics give 00=empty, 10=i→j, 01=j→i, 11=both directions.

- A general digraph allows all four states.
- A semicomplete digraph allows 10,01,11, but not 00.
- An oriented graph allows 00,10,01, but not 11.
- A tournament allows exactly 10,01.

Thus “no reciprocal arcs” is not enough to mean tournament: it also needs
“no missing pair.” Every oriented graph has tournament completions, exactly
2^m of them when m unordered pairs are missing. Conversely every oriented
graph is obtained from a tournament by deleting arcs. There is no choice
of orientation to resolve at an actual missing pair until a completion is
specified.

## 2. Three explicit ways to finish the proposed 11=empty convention

The user specifies 11=empty; that alone does not specify 00. We analyze
three coherent extensions without claiming the user chose one.

| raw state | ordinary | D: deletion alphabet | C: cancellation | S: state swap |
|---|---|---|---|---|
|00|empty|forbidden|empty|both directions|
|10|i→j|i→j|i→j|i→j|
|01|j→i|j→i|j→i|j→i|
|11|both directions|empty|empty|empty|

**D (ternary deletion encoding)** restricts raw matrices to no off-diagonal
00 pairs. It is a bijection between raw semicomplete matrices and ordinary
oriented graphs. The word “semicomplete” describes the raw ordinary reading,
not the decoded graph. The two graph readings need not share their path
count, strong components or Hamiltonicity.

**C (both equal states empty)** applies to every raw binary matrix. The
ordinary adjacency of its oriented output is

`B[i,j]=A[i,j](1-A[j,i])`, for i≠j.

Equivalently, with `K=A-A^T`, `B=(|K|+K)/2`, entrywise. It forgets whether an
empty pair was written 00 or 11. An oriented graph with m empty pairs has
exactly 2^m loopless raw representations; counting encodings is therefore
not the same as counting oriented graphs. In particular K alone cannot
recover ordinary adjacency or ordinary H.

**S (swap the two symmetric states)** is the loopless involution

`tau(A)=J-I-A^T`.

It swaps 00 with 11 and fixes 10 and 01, retaining all information. It
also preserves K. On the no-00 domain, S and C coincide with D; on general
raw matrices they do not. S exchanges semicomplete and oriented graphs,
while leaving every tournament fixed. This is complement followed by
arc reversal, not ordinary complementation alone.

No interpretation silently replaces 11=empty by an undirected edge. An
empty pair supplies neither direction for a path.

## 3. Counting encodings and random-model expectations

Let N=binomial(n,2), n>=1. The labelled loopless counts are 4^N general
adjacencies, 3^N oriented graphs, 3^N semicomplete graphs, and 2^N tournaments.
The deletion alphabet is another encoding of the 3^N oriented objects;
the redundant cancellation alphabet is a 4^N set mapping onto them.

For uniform independent raw pair states:

- D, uniformly among 10,01,11: E H(decoded)=n!/3^(n-1).
- The same raw arrays, ordinary reciprocal reading: E H=n!(2/3)^(n-1).
- C, uniformly among all four raw states: E H(decoded)=n!/4^(n-1).
- S, uniformly among all four states: the *entire* H distribution equals
  that of ordinary uniform binary adjacency, because S is a bijection.

Proof of the expectations: each prescribed permutation requires n-1
arcs from distinct unordered pairs, so its validity probability is the
product of n-1 identical marginal probabilities. Sum over n! permutations.
The ratio 2^(n-1) for the first two models is a ratio of expectations, not
an identity for any individual matrix.

The four-state cancellation distribution instead weights each oriented
graph by 2^m. This sampling bias is a concrete consequence of the 00
ambiguity, not a cosmetic representation choice. Independent diagonal
bits would multiply every encoding count by 2^n, but leave H unchanged.

## 4. A boundary for the parity bridge

Classical complement parity gives H(A)=H(tau(A)) modulo two. Therefore
ordinary and 11-empty readings have equal parity on the no-00 domain.
The parity statement does **not** extend to cancellation C on arbitrary
raw inputs. The matrix

```text
0100
1001
1000
0000
```

has only path 2013, hence H=1. Cancellation removes its reciprocal pair
{0,1} and leaves no Hamiltonian path, hence H=0. This is the first failed
implication in confusing C with S when 00 states also occur. The repair is
to retain the whole four-state table or restrict explicitly to no-00 input.

The complement parity theorem is classical, not a novelty claim. The
session's theory file gives a self-contained proof and the path-cover
identity behind it.

## Verification

`python3 models/verify_models.py` independently exhausts all 4,165 loopless
four-state matrices of orders one through four, checks the swap involution,
parity, skew-part preservation, cancellation fiber sizes, exact expected
counts, and the hostile. It separately checks all ternary arrays through
four. Output is frozen in `models/validation.json`.
