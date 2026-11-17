#!/usr/bin/env python3
"""residues.py (M-2): near-join misses (bases <= 9, nj9.bits up to 2^23) by residue class: mod 64, mod 8, mod 3, mod 7,
in the dyadic windows [2^19, 2^21) and [2^21, 2^23); the ranks of the classes 7 and 21 (mod 64)."""
import numpy as np, sys
bits = np.frombuffer(open(sys.argv[1] if len(sys.argv) > 1 else 'nj9.bits', 'rb').read(), dtype=np.uint8)
def member(v): return (bits[v >> 3] >> (v & 7)) & 1
for lo, hi in ((1 << 19, 1 << 21), (1 << 21, 1 << 23)):
    odd = np.arange(lo | 1, hi, 2, dtype=np.int64); m = member(odd) == 0; base = m.mean()
    print(f"window [{lo}, {hi}): odd numbers {len(odd)}, miss rate {base:.4f}")
    for mod in (64, 8, 3, 7):
        res = odd % mod; rates = {r: m[res == r].mean() / base for r in sorted(set(res.tolist()))}
        top = sorted(rates.items(), key=lambda kv: -kv[1])[:6]
        line = f"   mod {mod}: relative miss rate by class, top {top}"
        if mod == 64: line += f"; class 7: {rates[7]:.2f}, class 21: {rates[21]:.2f}; rank of 7: {sorted(rates, key=lambda r: -rates[r]).index(7) + 1}, rank of 21: {sorted(rates, key=lambda r: -rates[r]).index(21) + 1} of 32"
        if mod == 8: line += f"; 7 mod 8: {rates[7]:.2f} vs 1 mod 8: {rates[1]:.2f}"
        print(line)
