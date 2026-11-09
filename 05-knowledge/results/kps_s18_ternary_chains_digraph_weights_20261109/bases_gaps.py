#!/usr/bin/env python3
"""HYP-9530 (B) evidence: port chains with blocks of order <= 5 in bases b = 4..8.

For each base b this script:
- builds the block matrices (chain_blocks.py);
- runs chain_cover up to alpha-sum M (value 1 + b*M);
- adds every alpha-sum <= 80 that is not a known gap (S86's exact tables, HYP-9530 (B));
- closes the set under order joins: (1 + b m)(1 + b m') = 1 + b(m + m' + b m m').
It prints the alpha-sums in (80, M] that are still not attained. Only these can be further gaps of base b below M.

Usage: bases_gaps.py CHAIN_COVER_BINARY TMPDIR M
"""
import os
import subprocess
import sys

KNOWN_GAPS_LE_80 = {4: [3, 14, 27], 5: [3, 10, 27], 6: [3, 10, 14, 21, 46], 7: [3, 10, 14, 29], 8: [3, 35, 52]}


def main():
    prog, tmp, M = sys.argv[1], sys.argv[2], int(sys.argv[3])
    here = os.path.dirname(os.path.abspath(__file__))
    for b in range(4, 9):
        mats = os.path.join(tmp, f"m{b}k5.txt")
        subprocess.run([sys.executable, os.path.join(here, "chain_blocks.py"), str(b), "5", mats],
                       check=True, capture_output=True)
        out = subprocess.run([prog, mats, str(b), str(1 + b * M)], capture_output=True, text=True).stdout
        miss = [(int(t) - 1) // b for t in out.split("\n")[0].split(":")[1].split()]
        # chain_cover prints at most 200 misses: refuse to proceed if the list was truncated (audit m5)
        count = int(out.split("count")[1].split()[0])
        assert count == len(miss), f"base {b}: {count} misses but only {len(miss)} printed"
        att = set(range(M + 1)) - set(miss)
        att |= set(m for m in range(81) if m not in KNOWN_GAPS_LE_80[b])
        cl = set(att)
        changed = True
        while changed:
            changed = False
            for x in sorted(cl):
                if x == 0:
                    continue
                if 2 * x + b * x * x > M:
                    break
                for y in sorted(cl):
                    if y < x:
                        continue
                    z = x + y + b * x * y
                    if z > M:
                        break
                    if z not in cl:
                        cl.add(z)
                        changed = True
        rem = [m for m in range(81, M + 1) if m not in cl]
        print(f"base {b}: chain misses (alpha-sums <= {M}) {miss}")
        print(f"base {b}: known gaps <= 80 {KNOWN_GAPS_LE_80[b]}; undecided alpha-sums in (80, {M}]: {rem}")


if __name__ == "__main__":
    main()
