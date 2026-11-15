# Six matched arms at h=15, and a real twelve-count fold

The [fixed-state inversion](INVERSE_STATE_FIBRE.md) rules out closing the first eight prescribed counts at h=15. It does not forbid better windows at different absolute targets. A second search over the same strong-component alphabet finds six adjacent odds with exactly matching continuation data at h=15. Their fold gives twelve actual adjacent counts at infinitely many heights.

## A new six-arm physical bank

Write TT marks as literal 0/1 words in an ordinal sequence. Let F, B, C have outgoing rows (36,1,2,7,15,30), (8,1,3,6), (4,1,2). Write F_m for the block with an apex outgoing bit mask m, and likewise for B,C. The apex in every row spans the entire ordinal sequence.

| Spine blocks in order | Marked count N | h | S+h |
|---|---:|---:|---:|
| F_29, 100 | 979 | 15 | 1203 |
| 1, B_6, 1001, C_7 | 981 | 15 | 8·1203 |
| F_16, 1100 | 983 | 15 | 1203 |
| 101, F_48 | 985 | 15 | 1203 |
| C_1, B_6, 100 | 987 | 15 | 1203 |
| 1, C_3, B_2, 0 | 989 | 15 | 1203 |

The replay checks all six starting states by literal path DP. Each first strong component has a mark, so the apex reaches the whole spine. For any g≥4 append g ones except on the N=981 arm, where append g−3; then append a common 0. The previous [neutral alignment law](../chanpak_carry_admission_20261115/MATCHED_WINDOWS.md) shows that S+h became A=1203·2^g before that 0. The final-zero guard returns to the apex from every vertex, and the six new strong seeds have precisely
```
N=979+A+{0,2,4,6,8,10},    (h,S,d)=(15,2(A−15),0).
```
Thus every canonical common binary suffix keeps all six adjacent. At g=4 the strong window starts at 20227. This is a six-wide matched bank at h=15, different from the preceding order-nine-derived bank at h=45. No claim of maximality for h=15 is made.

## Folding doubles that width, but not the matched interface

For the [binary suffix lift](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md), put F(t)=979+A+S t+15p(t), S=2(A−15). The adjacent-address law in the [previous fold](../chanpak_tail_return_20261115/TAIL_FOLDS.md) gives
```
F(2u+1)−F(2u)=S−15(u+(u mod 2)).
```
Since 1203≡3 mod15, the nearest signed fold difference is c=6,12,−6,−12 at g≡0,1,2,3 mod4. Explicitly take u=(S−c)/15, a positive even integer. At g≡1 or 3 the two six-count windows start exactly twelve apart, so their odd values concatenate without holes. They give twelve consecutive odd counts at every such preparation height, with literal tails at addresses 2u and 2u+1. At the even phases the six-step overlap gives nine consecutive odds instead.

The two tails finish at different starts, as the terminal-word observer predicts. The twelve are actual path counts rather than a matched-width-12 recursion. This is a concrete threshold distinction: a larger local matched bank focuses into a longer absolute interval, but closing that interval still costs an additional matching operation.

Replay: `python replay_six_twelve.py` checks the six literal base graphs and the exact fold arithmetic. These are constructive witnesses for every declared clock, not asymptotic estimates.
