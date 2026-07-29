# Rank-selective alternating-pair census on 73 marked suffixes

**Status: FINITE-EXACT SCOPED APPLICATION OF THM-2897.  NO UNIFORM
SEVEN-BODY CLOSURE.**

For the four roots in THM-2895, retain the globally ranked first-apex order
and the literal excluded prefix on every one of the `73` suffix branches.
On a branch carrier of mass `h`, let `q_5` be the fifth largest allowed
singleton coverage and `B_2` the exact global allowed pair-union cap.
THM-2897 proves that

```text
q_5+2B_2<h
```

excludes a five-cover.

The locked census gives

```text
branches                                      73
scalar top-five closed                        48
rank-selective q5+2B2 closed                  50
adaptive union closed                         51
new closures beyond scalar                     3
scalar closures missed by q5+2B2               1
adaptive residual                             22
finite pair-union evaluations               1388
```

The three new strict margins are

```text
E=(2,8,9,10,11,13,14), rank 5, a=24:
  15798539/3546663120;

E=(1,3,9,10,11,12,14), rank 4, a=19:
  5904653/1485764280;

E=(1,3,9,10,11,12,14), rank 6, a=46:
  52421/23207184.
```

The hostile minimum is

```text
-18659449/446185740
```

at the first branch of `E=(2,8,9,10,11,13,14)`.  Thus the partition cap is
neither a replacement for the scalar test nor a replacement for THM-2895
parity descent.  Its role is a cheap monotone activation layer before
enumerating an H4 pair complex.

The theorem-bearing object is the marked suffix.  Replacing `q_5` by the
unmarked global rank, or deleting future gate labels, would not prove the
branch.

```text
04-computation/lrc14_j6_rank_selective_alternating_pair_73_thm2897.py
SHA-256 9f0453010dcd92b1258ec949f266f7ab03c0c0dae3fc6ba65bdcfe3b085ab922

05-knowledge/results/lrc14_j6_rank_selective_alternating_pair_73_thm2897.out
SHA-256 2ed7dadef1488de134e1e7c018b4b9d4490d022e01da23e1029d8e26d7bc2a8f
```

Ordinary and optimized replays are byte-identical, there are no Python
`assert` statements, and the canonical ledger digest is

```text
97e4f5c0b8194f98b568245e88517e2c239875fb0512bb1cb8b40e0e6f90161b.
```
