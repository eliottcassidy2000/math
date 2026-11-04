# Frozen reproduction record

2026-11-03. **FINITE-EXACT** producer replay; the unbounded proof in `PROOF.md`
remains **RESERVED** until independently audited. No finite experiment is used
as the proof of its arbitrary-block classification.

Commands, from the worktree root:

```text
python3 04-computation/multiport_grammar_20261103/closure.py
python3 -O 04-computation/multiport_grammar_20261103/closure.py
```

Both modes pass **213,973 always-active gates** and produce byte-identical
844-byte outputs, saved as `normal.out` and `optimized.out`.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| `04-computation/multiport_grammar_20261103/closure.py` | 11,305 | `3878a30e9e20ec1ae29caa5a523f0bd1df62a312c458c6e2014d11223cd6c962` |
| `normal.out` | 844 | `d26f2e4f531305d54134a5762c736b134d3ddd0e94a3cff695a1706e451b53c3` |
| `optimized.out` | 844 | `d26f2e4f531305d54134a5762c736b134d3ddd0e94a3cff695a1706e451b53c3` |

The declared labeled universe contains 1,099 tournaments of orders 1–5 and
33,866 actual tournament/port states. Of these, 27,475 meet the terminal-port
strongness guard. Every actual new block of orders 1–3, including its empty port
set, is retained in the raw extension bank. The program makes 4,066,300 raw
append/prepend slope comparisons. It finds exactly 775 first appends and 311
first prepends in that finite universe, and rejects 127,062 guarded candidate
second extensions. Direct ordered-cover enumeration separately checks 1,380
instances of the full seven-coordinate transfer law.

Finite append offsets are `{-42,-20,-16,-12,-10,-4}`. This is deliberately **not**
an unbounded sign theorem. The prepend digit is analytically negative for every
size, by `f=h`, the actual marked cut injection, and the single-port joint zero.

The code contains three distinct scope controls:

1. an actual first tripling step, of digit `-4`, proving nonvacuity;
2. actual strong uniform `+2` and `+4` constructions, whose differing adjacency
   or old-port data identify the exact forbidden overgeneralization;
3. a nonnegative formal signature permitting two steps when the even-`p`
   realization sidecar is discarded.

The independent subset/end-point DP does not call the ordered permutation
enumerator. Physical port masks and old-state embeddings are checked directly;
no quotient based on a matching marginal is used. No native operation schedule
is claimed by these static adjacency controls.
