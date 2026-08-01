---
id: THM-3007
title: "Dyadic block rigidity for critical-run fair extractors"
status: >
  PROVED + VERIFIED-EXACT. Call [N, N+l) a BLOCK if a fair-coin extractor
  decides every critical value in [N, N+l) by flip N+l, and BALANCED if its
  terminal set (b^N y, y != b^l, both b) has every fixed-composition class of
  even size -- the condition making the shell/exactness lemma apply block by
  block. The class-size generating function is (1+u^N)(1+u)^l - 1 - u^(N+l),
  so a balanced block exists iff (1+u^N)(1+u)^l = 1 + u^(N+l) over F_2, and
  that holds IFF N = 2^a and N+l = 2^g with g > a. Hence every balanced block
  is a dyadic interval [2^a, 2^g): the dyadic shell structure of AMM 12592's
  classical solutions is FORCED, not a design choice. Corollaries: l >= N, so
  no balanced block has ratio (N+l)/N below 2; block boundaries at least
  double; and the block containing n is [2^a, 2^g) with 2^a <= n < 2^g.
  This generalizes THM-2160 S6.1 (the l = N case) to the full two-parameter
  family and closes the "shell of ratio below two" route of HYP-9061 sec 2b.
  IMPORTANT SCOPE: it constrains the BLOCK GEOMETRY, not the deadline slope;
  see THM-3008, which builds balanced block rules of slope well below 2.
source: opus-2026-07-31-amm12592-writeup
depends_on: []
related:
  - THM-2160  # S6.1: composition-exact h-shell extraction iff h is dyadic
  - THM-2225
  - THM-3032  # sharpened half-tail extractor
  - THM-3008  # within-block deadline profiles; slope below 2
  - HYP-9061  # the minimal linear deadline C*
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_dyadic_block_rigidity_thm3005.py
output: 05-knowledge/results/amm12592_dyadic_block_rigidity_thm3005.out
---

# THM-3007 -- balanced blocks are exactly the dyadic intervals

Bits are independent with `P(0)=p`, `P(1)=q=1-p`, `0<p<1` unknown; `n` is
the length of the maximal constant initial run.

**Definition.** A *block* `[N, N+l)` is a rule that decides every critical
value `n` in `[N, N+l)` by flip `N+l`. Its terminal set is

```text
B(N,l) = {0^N y : y in {0,1}^l, y != 0^l}
       u {1^N y : y in {0,1}^l, y != 1^l}.
```

The block is *balanced* if every fixed-composition class of `B(N,l)` has even
size, so that the class can be split evenly between heads and tails. By the
exactness lemma (THM-2225 sec 3), a scheme whose blocks are all balanced is
exactly fair for every `p`.

## 1. The criterion

Counting by total weight,

```text
sum_k |class k| u^k
  = [(1+u)^l - 1] + u^N[(1+u)^l - u^l]
  = (1+u^N)(1+u)^l - 1 - u^(N+l).                          (1)
```

Hence

```text
[N, N+l) is balanceable
  <=>  (1+u^N)(1+u)^l = 1 + u^(N+l)   in F_2[u].            (2)
```

## 2. Classification

**Theorem.** (2) holds if and only if `N = 2^a` and `N + l = 2^g` for some
`g > a >= 0`.

*Sufficiency.* If `N = 2^a` then `1+u^N = (1+u)^N` over `F_2`, so the left
side of (2) is `(1+u)^(N+l)`; if `N+l = 2^g` this is `1 + u^(2^g)`. QED

*Necessity.* First `N <= l`: if `N > l`, the coefficient of `u^l` on the left
of (2) is `binom(l,l) = 1` (the factor `u^N` contributes only in degrees
`>= N > l`), while on the right it is `0` because `0 < l < N+l`.

So assume `N <= l`. For `1 <= k <= N-1` the left side has coefficient
`binom(l,k)` and the right side `0`, so

```text
binom(l,k) = 0 mod 2        for 1 <= k <= N-1.               (3)
```

At `k = N` the left side is `binom(l,N) + 1` and the right side is `0`
(as `N < N+l`), so `binom(l,N)` is odd. Combining with (3),

```text
N = min{k >= 1 : binom(l,k) odd} = 2^(v_2(l))               (4)
```

by Lucas's theorem. In particular `N` is a power of two and `N | l`. Then
`1+u^N = (1+u)^N` over `F_2` and (2) reduces to
`(1+u)^(N+l) = 1 + u^(N+l)`, which holds iff `N+l` is a power of two. QED

## 3. Corollaries

1. **Forced dyadic shells.** Every balanced block scheme partitions the
   critical values into consecutive dyadic intervals
   `[2^(a_0), 2^(a_1)), [2^(a_1), 2^(a_2)), ...` with `0 = a_0 < a_1 < ...`.
   The standard shells `[2^r, 2^(r+1))` of THM-2160/2225/2996 are the finest
   admissible choice; any other scheme merely coarsens them.
2. **No sub-ratio-two shell.** `l >= N`, so `(N+l)/N >= 2`. HYP-9061 sec 2b's
   "ratio `rho < 2`" shells do not exist in balanced form at all; its corner
   analysis, which assumed such a shell and asked whether the two forced
   corner deficits could be routed to meet, is moot for balanced blocks. Any
   `rho < 2` mechanism must break block balance (cross-block compensation).
3. **Geometric boundaries.** `N_(j+1) >= 2 N_j`, so at most `log_2 n + 1`
   blocks lie below `n`.
4. **THM-2160 S6.1 is the diagonal.** That result is the case `l = N`:
   `[N, 2N)` is balanced iff `N` is a power of two. (2) is the full
   two-parameter statement.

## 4. Scope -- what this does NOT say

It constrains the block *geometry*, not the *deadline*. Inside a block
`[m, 2m)` a rule may decide many critical values long before flip `2m`, and
THM-3008 exhibits balanced block rules whose worst within-block deadline
ratio is `3/2, 14/9, 25/16, 11/7` for `m = 4, 8, 16, 32` -- all far below the
block ratio `2`. So "balanced blocks have ratio `>= 2`" must NOT be read as
"balanced schemes have deadline slope `>= 2`"; that inference is false.

## 5. Referee

```bash
python3 04-computation/amm12592_dyadic_block_rigidity_thm3005.py
```

Checks, for all `1 <= N, l <= 299`, that the three predicates agree --
"every class of `B(N,l)` has even size", the `F_2` criterion (2), and
"`N` and `N+l` are powers of two" -- with a direct enumeration of `B(N,l)`
as a fourth check for `l <= 12`. Zero mismatches; no balanced block with
`l < N`; no balanced block with non-dyadic `N`. QED.
