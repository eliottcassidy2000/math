# Exact binary folds enlarge a four-count window to eight

Moving to the next binary **address** is not the same operation as appending another physical vertex. Consecutive addresses correspond to two different literal painted tails of the same anchored source. In the [suffix decoder](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md), the exact output is

```
F(a)=N_0+Sa+h p(a),
p(0)=0,    p(2u+1)=p(u),    p(2u)=p(u)+u+(u mod 2).
```

Consequently the graph of the tail function has an explicit alternating fold:

```
F(2u+1)−F(2u)=S−h[u+(u mod 2)].
```

For an odd address a=2^r(2v+1)−1 with r≥1 trailing ones, the other step is

```
F(a+1)−F(a)=S+h[(2^r−1)(2v+1)−v+1−(v mod 2)] ≥ S+2h.
```

The second identity strips the r trailing ones from a, while on a+1 the first doubling of the odd root contributes the sole adjacency bonus. Thus odd-to-even addresses always increase, whereas even-to-odd addresses eventually decrease; the two signs need not reflect any irregularity of the rule. For even S>0, the even-to-odd step is exactly ±2 precisely when M=(S∓2)/h is a nonnegative even integer and u is M or M−1 (discard the latter when M=0). In particular, positive S has a unit step on an even-to-odd address pair iff S≡±2 mod 2h. The condition is exact and follows from the two possible preimages of the round-up operation u↦u+(u mod 2).

## Physical enlargement at infinitely many preparation heights

Use the [four matched seeds](../chanpak_carry_admission_20261115/MATCHED_WINDOWS.md) at any preparation height j≥3. They have the same h=15 and S=2(17·2^j−15). Their initial counts are 27+17·2^j+{0,2,4,6}. At j≡3 mod 4, S≡2 mod 30. Put

```
M=(S−2)/15 = (34·2^j−32)/15,    a=2M.
```

Here M is positive even. The exact fold gives F(a+1)=F(a)+2 separately for each of the four matched seeds. Hence the two translated four-count windows overlap in three counts and their union is the entire five-count window

```
27+17·2^j+Sa+15p(a)+{0,2,4,6,8}.
```

At j=3, a=32, this is **8387, 8389, 8391, 8393, 8395**. The replays below count both literal tails independently. All the final tails are nonempty and end in mark zero. The two address choices do not have the same final start readout, so the union is **not** a width-five matched continuation state. It enlarges the actually attained interval and gives an explicit scale-by-scale position, a different goal from closing the interface for further growth.

## Other phases give six and eight

For the same seed bank, the nearest possible signed even-to-odd difference is determined by S modulo 30. Each row below uses the unique smallest feasible |S−15M| over nonnegative even M. An odd-to-even step is always too large. Overlaying two four-count windows at shift c joins them into one consecutive odd window exactly when 0≤c≤8; the new width is 4+c/2. This gives the following sharp widths **for one neighboring-address pair** at each preparation phase.

| j mod 4 | S mod 30 | nearest signed difference | maximum width by this trick |
|---:|---:|---:|---:|
| 3 | 2 | +2 | 5 |
| 0 | 4 | +4 | 6 |
| 1 | 8 | +8 | 8 |
| 2 | 16 | −14 | 4 (disjoint windows) |

For c=4 or 8, the exact choice is again M=(S−c)/15, a=2M. At j=4 this gives six consecutive odds from **36151** to **36161**; at j=5 it gives eight from **150401** to **150415**. Both forms hold at every preparation height in their corresponding phase, by the exact continuation and fold laws. These are independent constructions from the previously found slope-17 width-eight lane, and still do not match the full right interface across the two addresses.

There is also a literal period hiding in these fold addresses. Write j=4k+3. At k=0 we have a=32. For k≥1,

```
a=64 u_k,    u_k=(136·16^(k−1)−1)/15.
```

Thus the ordinary binary digits of a are `1001 (0001)^(k−1) 000000` at these later heights. The exact recursion is p(u_1)=4 and p(u_(k+1))=p(u_k)+7u_k+1, followed by p(a)=p(u_k)+63u_k+1. For c=4 and j=4k≥4 the hexadecimal digits of a are `4 (8)^(k−1) 4`; for c=8 and j=4k+1≥5 they are `8c` at k=1 and `9 (1)^(k−2) 0c` at k≥2. These are immediate from a=(68·2^j−60−2c)/15. This illustrates how a periodic address word, exact accumulation and a local fold can act in one construction. It does not imply arbitrary placement of the enlarged window or a universal count theorem.

Replay: `python replay_tail_folds.py` verifies the local formulas across 4096 addresses, all three periodic expressions at successive heights, and all eight literal first-height glues using subset-endpoint DP.
