# Threads, typed: the self-converse peak at 5 and the MRB terms; primes just above 10^k; 15, 63 and the holes; 756·2^s − 33

Owner's pasted block, item by item. Wave 54 §2–3 of `chanpak_waves_20261113` already typed 63 (Zsigmondy), 756 =
4·H(P_7) = 12·63 and the Ogg primes; wave 56 §3 the Mandelbrot count `2^{p−1} − 1`. Nothing below re-opens those.

## 1. The ratio of self-converse classes to converse pairs peaks at order 5 ([receipts/sc_ratio.out](receipts/sc_ratio.out))

Canonical labelling (nauty `labelg`) of every class and of its converse, orders 3–9:

| n | classes | self-converse | converse pairs | SC/pairs | strong | SC strong | pairs strong | ratio |
|---|---|---|---|---|---|---|---|---|
| 3 | 2 | 2 | 0 | ∞ | 1 | 1 | 0 | ∞ |
| 4 | 4 | 2 | 1 | 2 | 1 | 1 | 0 | ∞ |
| 5 | 12 | 8 | 2 | **4** | 6 | 6 | 0 | ∞ |
| 6 | 56 | 12 | 22 | 0.545 | 35 | 9 | 13 | 0.692 |
| 7 | 456 | 88 | 184 | 0.478 | 353 | 79 | 137 | 0.577 |
| 8 | 6880 | 176 | 3352 | 0.0525 | 6008 | 162 | 2923 | 0.0554 |
| 9 | 191536 | 2752 | 94392 | 0.0292 | 178133 | 2661 | 87736 | 0.0303 |

The self-converse counts are A002785 (2, 2, 8, 12, 88, 176, 2752); the ratio over all classes is finite from `n = 4` and
**maximal at `n = 5`** (four self-converse classes per converse pair), then falls by a factor 7 at 6 and 9 at 8. Every
strong tournament of order ≤ 5 is self-converse (no converse pair before order 6). Q-5's predicted values were wrong
(recorded); the location of the peak is as the owner remembered. THM-4922 holds the odd-order balance of the
self-converse classes between `H ≡ 1` and `3 (mod 4)`; HYP-2080 the Burnside count.

**MRB.** The MRB constant `Σ_{k≥1} (−1)^k (k^{1/k} − 1) = 0.18785964…` has terms `k^{1/k} − 1 = 0, 0.414, 0.442, 0.414,
0.380, 0.348, …`, maximal at `k = 3` because `x^{1/x}` peaks at `x = e`. Both sequences are unimodal with a small-argument
peak — the owner's microcosm (rise), inflection (peak), macrocosm (decay) — but the mechanisms differ: the tournament
ratio decays superexponentially (self-converse classes are a square-root-sized fibre of the class count, HYP-2080), the MRB
terms like `(ln k)/k`. **Type: ANALOGY** (shape only).

## 2. Primes just above 10^k ([receipts/primes_near.out](receipts/primes_near.out))

| k | least prime > 10^k | least twin pair > 10^k | Hardy–Littlewood pairs expected in (10^k, 10^k + 100] |
|---|---|---|---|
| 3 | +9 | (+19, +21) | 2.77 |
| 4 | +7 | **(+7, +9)** | 1.56 |
| 5 | +3 | (+151, +153) | 1.00 |
| 6 | +3 | (+37, +39) | 0.69 |
| 7 | +19 | (+139, +141) | 0.51 |
| 8 | +7 | (+37, +39) | 0.39 |
| 9 | +7 | **(+7, +9)** | 0.31 |
| 10 | +19 | (+277, +279) | 0.25 |
| 12 | +39 | (+61, +63) | 0.17 |
| 15 | +37 | (+2371, +2373) | 0.11 |

The owner's pairs: `10007, 10009` and `10^9 + 7, 10^9 + 9` are the only `k ≤ 15` with the least twin pair at offsets
`(7, 9)` (Q-6, second clause confirmed); `30011, 30013` are twin primes at `3·10^4 + (11, 13)`; `30021, 30027 = 3·(10007,
10009)` — tripling turns the `+2` pair into a `+6` pair. A twin pair within 100 of `10^k` occurs for 9 of the 15 values of
`k` (Q-6's "at most 5" refuted), in line with the Hardy–Littlewood expectation summed over `k`. The trend is the
Cramér-scale fact that the next prime after `N` is at `N + O(log² N)`; the pattern of offsets `7, 7, 7` at `k = 4, 8, 9`
is **NUMEROLOGY**.

**"Close enough approximations that select the integer."** In the count monoid the analogue is wave 56's target equation
`p = h + 2^k B`: the base `h` approximates `p` modulo `2^k` and the response `B = (p − h)/2^k` is then forced — the 2-adic
approximation selects the integer, exactly as `nextprime` selects the integer above `10^k`. Tripling `(10007, 10009) ↦
(30021, 30027)` is the closure `3·NJ ⊆ NJ` of THM-5264 (E); recovering a `+2` pair at the tripled level is THM-5260's
forcing problem, which tripling itself cannot solve. **Type: ANALOGY.**

## 3. 15, 63, the holes, and the periods whose component count contains both

For a prime period `p`, the Mandelbrot set has `ν_2(p) = 2^{p−1} − 1` hyperbolic components (`1, 3, 15, 63, 1023, 4095` for
`p = 2, 3, 5, 7, 11, 13`; wave 56 §3). `63 = 2^6 − 1` is Zsigmondy's exception (wave 54 §2). The trend behind "63 contains
both holes" is a theorem:

**Theorem 3 (arithmetic).** For a prime `p ≥ 5`: `7 | ν_2(p)` iff `p ≡ 1 (mod 3)`; `9 | ν_2(p)`, `21 | ν_2(p)` and
`63 | ν_2(p)` are each equivalent to `p ≡ 1 (mod 6)`. (The orders of 2 modulo 7, 9, 21, 63 are 3, 6, 6, 6.) So the
periods whose component count is divisible by both holes are exactly `p ≡ 1 (mod 6)`: `7, 13, 19, 31, 37, 43, 61, 67, 73,
79, 97` below 100 (Q-7 checked). `15 = ν_2(5) = H(RT_5)` is the only value shared with the regular tournament on five
vertices; `3, 15, 63, 1023 = 3·11·31, 4095 = 35·117` are all Hamiltonian-path counts; the Fermat quotients `ν_2(p)/p` are
`3, 9, 93, 315, 3855, 13797`. **Type: arithmetic (Theorem 3) / NUMEROLOGY (the link to the holes).**

## 4. `756·2^s − 33 = 3·(63·2^{s+2} − 11)` and the window at `s = 3` ([receipts/nj9.out](receipts/nj9.out))

Wave 54 read the family as `189·2^{s+2} − 33` (dyadic multiples of the Paley count `H(P_7) = 189`); the owner's reading
`3·(63·2^{s+2} − 11)` is the same numbers. Set `H_s = 63·2^{s+2} − 11`: 241, 493, 997, 2005, 4021, 8053, 16117, 32245,
64501, 129013. Every one is a near-join value, all from bases of order ≤ 8:

| `H_s` | base order | `h` | `B` | `k` |
|---|---|---|---|---|
| 241 | 6 | 17 | 28 | 3 |
| 493 | 5 | 13 | 15 | 5 |
| 997 | 7 | 101 | 112 | 3 |
| 2005 | 7 | 69 | 121 | 4 |
| 4021 | 7 | 117 | 244 | 4 |
| 8053 | 8 | 117 | 496 | 4 |
| 16117 | 8 | 373 | 492 | 5 |
| 32245 | 8 | 373 | 498 | 6 |
| 64501 | 8 | 245 | 251 | 8 |
| 129013 | 8 | 245 | 503 | 8 |

The family is **not a near-join line**: `H_{s+1} = 2H_s + 11`, while a fixed pair `(h, B)` gives `h + 2^{k+1}B = 2(h + 2^k B) − h`,
so a line through the family would need `h = −11`. The representations above change base at every step (`h = 17, 13,
101, 69, 117, 117, 373, 373, 245, 245`) — the dyadic approximation selects the integer anew each time. The window at
`s = 3`, `6015 = 3·2005`, `6017`, `6019`, is fully realised: `6015 = 111 + 2^4·369`, `6017 = 65 + 2^5·186`,
`6019 = 163 + 2^4·366` (bases of orders 7, 7, 8); `6015` is also the closure `3·2005`. Among one-vertex-core values
`h(D) + P_D(2)` (strong `D ≤ 8`), `241 = 17 + 224` occurs and 493, 997, 2005 do not (Q-9). **Type: arithmetic.**
