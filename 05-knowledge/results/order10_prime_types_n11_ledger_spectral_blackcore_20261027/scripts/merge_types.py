"""Merge the eight order-10 TYPE part files into one (sum counts per signature)."""
import collections, re, sys, glob
agg = collections.Counter()
for fn in sorted(glob.glob(sys.argv[1] if len(sys.argv) > 1 else 'prime_types_n10_part*.txt')):
    for line in open(fn):
        m = re.match(r'TYPE count=(\d+) (.*)$', line.strip())
        if m: agg[m.group(2)] += int(m.group(1))
for k, v in sorted(agg.items(), key=lambda kv: (-kv[1], kv[0])): print(f"TYPE count={v} {k}")
print(f"# total primes {sum(agg.values())}, self-converse {sum(v for k, v in agg.items() if 'SC=1' in k)}, types {len(agg)}", file=sys.stderr)
