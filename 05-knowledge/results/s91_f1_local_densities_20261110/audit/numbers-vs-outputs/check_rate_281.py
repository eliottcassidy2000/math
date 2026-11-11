#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: max_k rho(T_k)/2 at the prime modulus 281 (HYP-9660 (B): 0.8130), full operator."""
import numpy as np
from verify_misc import twisted_matrix

m = 281
best, arg = 0.0, 0
for k in range(1, (m - 1) // 2 + 1):
    r = float(np.max(np.abs(np.linalg.eigvals(twisted_matrix(m, k))))) / 2
    if r > best:
        best, arg = r, k
print(f"m = {m}: max_k rho/2 = {best:.6f} at k = +-{arg}")
