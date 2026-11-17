#!/usr/bin/env python3
"""residues10.py: the 1,197 near-join misses of bases <= 10 below 2^25 (nj10.bits) all lie in [2^24, 2^25); their
residue classes mod 64, 8, 3, 7 against the odd numbers of that window; the smallest ten."""
import numpy as np
bits = np.frombuffer(open('nj10.bits', 'rb').read(), dtype=np.uint8)
lo, hi = 1 << 24, 1 << 25; odd = np.arange(lo | 1, hi, 2, dtype=np.int64); miss = odd[((bits[odd >> 3] >> (odd & 7)) & 1) == 0]
print(f"window [2^24, 2^25): odd numbers {len(odd)}, misses {len(miss)} (rate {len(miss)/len(odd):.6f}); smallest ten: {miss[:10].tolist()}")
for mod in (64, 8, 3, 7):
    r = np.bincount((miss % mod).astype(int), minlength=mod); exp = len(miss) / (mod // 2 if mod % 2 == 0 else mod)
    cls = [(c, int(r[c]), round(int(r[c]) / exp, 2)) for c in range(mod) if r[c]]
    line = f"   mod {mod}: (class, misses, relative): {sorted(cls, key=lambda t: -t[1])[:8]}"
    if mod == 64: line += f"; class 7: {int(r[7])}, class 21: {int(r[21])} (expected {exp:.1f} each)"
    print(line)
