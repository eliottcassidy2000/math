# H6 closed-danger replay validation

Date: 2026-07-15  
Theorem context: THM-857  
Arithmetic: exact integer/rational; no floating-point branch decisions

## Result

The independent closed-danger replay completed all `924` deletion roots.  The
canonical combined census is

```text
nodes       924,83881,8906315,559202706,12671505,53812,21
dead        0,0,0,555565824,12638291,53792,0
full        0,0,0,495797163,0,0,0
early       0,0,0,59768661,12638291,53792,20
edges       580918240
covers      1
loose       20
```

The sole cover is root `286`, with missing labels
`1,3,5,7,9,11` and ordered lifts
`1:14,3:16,5:18,7:20,9:22,11:24`.  Its root certificate is
`cb9d82c4580b7d9253914412df1d25a723539baef232802cf3da646e6649a4d5`.

The canonical all-root manifest is
`480ca266dba54d7dfa76a68baecd56b26f3294095d36786ab04c0dfd5b135dee`.

## Exact provenance hashes

```text
replay source   585c50601ef3c392a2afb020242115a7215d233878639c80d355984a8a16bd27
combiner source dd383b53015bea5e7aa76c1109ca4f382ca411611be726910951d8c4bf6b07e7
stored output   6d7316ceb1fe15987bb527f75ff7938ceef60a118facdbf697315da336a393be
shard executable 5d45f8c240ede29c7d4b2fc3163f446291ce98a4ed459cfc802ef7dc49c61bb5
```

The shard executable was compiled by Apple clang `17.0.0` for
`arm64-apple-darwin24.3.0`.  A fresh rebuild from the hashed source has the
same symbols and full disassembly and byte-identical `__text`, `__const`, and
`__cstring` sections.  Their hashes are respectively
`9daa73720dd8c261af39e8ff04354e951411b87068c36ce0ed48c93065d07e18`,
`c0b485c2aa41ddb0fd9070f872204ac5bd2dd028cba24421680f803f333fdd89`,
and `b2113bf15d05274b34c2de1e75f99fc366dcffa853f6a445d45219cac4e9c3bd`.
The whole Mach-O rebuild differs because the linker regenerates UUID/signature
metadata; no code, constant, or string section differs.

## Reproduction and aggregate validation

The strict optimized build command was

```bash
c++ -std=c++20 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow \
  04-computation/lrc13_hamming_six_closed_danger_union_replay_codex_S10.cpp \
  -o /tmp/lrc_h6_closed_O3_final
```

To balance measured deletion-root cost, the exact contiguous shard ranges
were

```text
0 214
214 93
307 131
438 70
508 65
573 87
660 100
760 164
```

Each range was run as

```bash
env LC_ALL=C LANG=C /tmp/lrc_h6_closed_O3_final \
  --root-start START --root-limit LIMIT \
  > /tmp/h6-independent-final/shardN.out \
  2> /tmp/h6-independent-final/shardN.err
```

All eight processes exited zero and emitted one `PASS` line.  The stored
transcript was then produced by

```bash
python3 -B \
  04-computation/lrc13_hamming_six_closed_danger_union_replay_combine_codex_S10.py \
  /tmp/h6-independent-final/shard{0..7}.out \
  > 05-knowledge/results/lrc13_hamming_six_closed_danger_union_replay_codex_S10.out
```

Because execution was sharded, the C++ program's monolithic
`check_all_root_aggregate` branch was not invoked.  This is not silently
treated as a monolithic run.  The separate combiner mechanically checked:

1. every shard's declared contiguous range, local totals, local manifest,
   cover-row count, and unique `PASS` line;
2. each root `0..923` exactly once, in the lexicographic order of the
   `C(12,6)` deletion sets;
3. every frozen global counter above and the literal sole-cover row; and
4. the canonical SHA-256 stream obtained from all `924` root rows in root
   order.

Only after all checks passed did it emit the `930`-line, `258721`-byte stored
output.

## Optimization and sanitizer diagnostics

Roots `0` and `286` were rerun under strict `-O3`, `-O0`, and
`-O1 -fsanitize=address,undefined -fno-omit-frame-pointer` builds.  Their
stdout was byte-identical across all three builds.  With
`ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1`, both sanitizer
runs exited zero without diagnostics.  Leak detection is unavailable in this
macOS ASan runtime and was not claimed.

```text
root 0 stdout SHA-256   d28656fd853f74f17f5fc453a72bfd0c07aa3968a535ec1e627d96c21d822533
root 286 stdout SHA-256 14b0f247d765ba955ab870181624c606ece571dc2d75dadbe08321d93dd582d0
root 462 O3 stdout SHA-256 6bdcb4c51661cacb0e5a852c00a77ba91989e2d9d24a9c04a179db16a3dbcd35
```

Roots `0`, `286`, and `462` also carry hard-coded diagnostic census
assertions in the replay source.

## Tournament and assumption audit

This certificate keeps literal open components rather than quotienting them
to a tournament.  Candidate lifts, gap components, wall crossings, and
cap-proof obligations were all considered as alternate vertex sets; no
tested pairwise orientation preserves the exact cover predicate.  The
root-level missing-label tournament therefore remains planning telemetry in
the primary verifier.  Importing it into this replay would weaken algorithmic
independence without changing the recursion.  The replay's literal component
carrier destroys no endpoint or remaining-comb information.
