# THM-4277 reproduction

Run from any directory inside the current repository checkout, keeping the
following variables in the same shell for all commands below:

```bash
REPO="${REPO:-$(git rev-parse --show-toplevel)}"
SRC="$REPO/04-computation/lrc14_two_dimensional_outsider_rectangle_thm4277"
OUT="$REPO/05-knowledge/results/lrc14_two_dimensional_outsider_rectangle_thm4277"
test -f "$SRC/primary_endpoint.cpp"
test -f "$OUT/primary_endpoint_O3.out"
```

For an unpacked promotion packet rather than a Git checkout, set `REPO`
explicitly to the absolute path of its `repo-ready` directory. The two C++
programs depend only on the C++20 standard library and
`lrc_rectangle_common.hpp`. The proof-graph postprocessor additionally reads
the frozen post-THM-4266 residual from the checkout supplied by `REPO`; it does
not edit that checkout.

All executable correctness gates use `require`/exceptions and remain active
under `-DNDEBUG`. There are no load-bearing C or C++ `assert` calls. Wall-clock
lines are removed before byte comparison and canonical hashing.

## 1. Primary endpoint and independent literal-wall O3 replays

```bash
c++ -std=c++20 -O3 -pthread "$SRC/primary_endpoint.cpp" \
  -o /tmp/thm4277-primary-O3
c++ -std=c++20 -O3 -pthread "$SRC/independent_literal.cpp" \
  -o /tmp/thm4277-literal-O3

/tmp/thm4277-primary-O3 > /tmp/thm4277-primary-O3.raw
/tmp/thm4277-literal-O3 > /tmp/thm4277-literal-O3.raw
sed '/^SECONDS /d' /tmp/thm4277-primary-O3.raw > /tmp/thm4277-primary-O3.out
sed '/^SECONDS /d' /tmp/thm4277-literal-O3.raw > /tmp/thm4277-literal-O3.out
cmp /tmp/thm4277-primary-O3.out "$OUT/primary_endpoint_O3.out"
cmp /tmp/thm4277-literal-O3.out "$OUT/independent_literal_O3.out"
python3 -B "$SRC/verify_transcripts.py" \
  "$OUT/primary_endpoint_O3.out" "$OUT/independent_literal_O3.out"
```

The primary sorts all `binom(30,8)` masks and integrates primitive endpoint
arcs. The independent engine uses a bounded heap for candidate selection,
literal `lcm(D,14q,14r)` grids, direct interval-list intersection, and a
recursive body generator. The fixed-pool cell/repair-geometry builder and
basic mask/FNV utilities are intentionally shared.

## 2. NDEBUG replay

```bash
c++ -std=c++20 -O3 -DNDEBUG -pthread "$SRC/primary_endpoint.cpp" \
  -o /tmp/thm4277-primary-O3-NDEBUG
c++ -std=c++20 -O3 -DNDEBUG -pthread "$SRC/independent_literal.cpp" \
  -o /tmp/thm4277-literal-O3-NDEBUG

/tmp/thm4277-primary-O3-NDEBUG > /tmp/thm4277-primary-O3-NDEBUG.raw
/tmp/thm4277-literal-O3-NDEBUG > /tmp/thm4277-literal-O3-NDEBUG.raw
sed '/^SECONDS /d' /tmp/thm4277-primary-O3-NDEBUG.raw \
  > /tmp/thm4277-primary-O3-NDEBUG.out
sed '/^SECONDS /d' /tmp/thm4277-literal-O3-NDEBUG.raw \
  > /tmp/thm4277-literal-O3-NDEBUG.out
cmp /tmp/thm4277-primary-O3-NDEBUG.out "$OUT/primary_endpoint_O3.out"
cmp /tmp/thm4277-literal-O3-NDEBUG.out "$OUT/independent_literal_O3.out"
```

## 3. Undefined-behaviour sanitizer replay

```bash
c++ -std=c++20 -O1 -g -DNDEBUG -pthread -fsanitize=undefined \
  -fno-omit-frame-pointer "$SRC/primary_endpoint.cpp" \
  -o /tmp/thm4277-primary-ubsan
c++ -std=c++20 -O1 -g -DNDEBUG -pthread -fsanitize=undefined \
  -fno-omit-frame-pointer "$SRC/independent_literal.cpp" \
  -o /tmp/thm4277-literal-ubsan

UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  /tmp/thm4277-primary-ubsan \
  > /tmp/thm4277-primary-ubsan.raw 2> /tmp/thm4277-primary-ubsan.err
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  /tmp/thm4277-literal-ubsan \
  > /tmp/thm4277-literal-ubsan.raw 2> /tmp/thm4277-literal-ubsan.err
sed '/^SECONDS /d' /tmp/thm4277-primary-ubsan.raw \
  > /tmp/thm4277-primary-ubsan.out
sed '/^SECONDS /d' /tmp/thm4277-literal-ubsan.raw \
  > /tmp/thm4277-literal-ubsan.out
cmp /tmp/thm4277-primary-ubsan.out "$OUT/primary_endpoint_O3.out"
cmp /tmp/thm4277-literal-ubsan.out "$OUT/independent_literal_O3.out"
test ! -s /tmp/thm4277-primary-ubsan.err
test ! -s /tmp/thm4277-literal-ubsan.err
```

## 4. Detached 256-bit comparator audit

```bash
python3 -B "$SRC/u256_comparator_audit.py" > /tmp/thm4277-u256.out
python3 -O -B "$SRC/u256_comparator_audit.py" \
  > /tmp/thm4277-u256-opt.out
cmp /tmp/thm4277-u256.out "$OUT/u256_comparator_audit.out"
cmp /tmp/thm4277-u256-opt.out "$OUT/u256_comparator_audit.out"
```

This exact Python-bigint path mirrors the eight-limb schoolbook product and
checks 81 boundary products plus 100,000 deterministic random comparisons.

## 5. Current proof graph after canonical THM-4276

```bash
python3 -B "$SRC/postprocess_current.py" --repo "$REPO" \
  > /tmp/thm4277-postprocess.out
python3 -O -B "$SRC/postprocess_current.py" --repo "$REPO" \
  > /tmp/thm4277-postprocess-opt.out
cmp /tmp/thm4277-postprocess.out "$OUT/postprocess_current.out"
cmp /tmp/thm4277-postprocess-opt.out "$OUT/postprocess_current.out"
```

The postprocessor hard-checks every inherited ledger through THM-4276, the
full rectangle ratio/gcd census, the exact 2,419-row novel deletion, the
172,322-row final residual, and the two-edge endpoint-670 top layer.
