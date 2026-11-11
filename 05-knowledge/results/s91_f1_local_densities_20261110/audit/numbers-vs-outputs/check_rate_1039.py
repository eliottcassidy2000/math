#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: rho(T_k)/2 at m = 1039, k = 64 (the maximiser reported in bridges/b3_f1_data.out),
from the full 2m x 2m operator of THM-5128 (B); and a few other k for comparison (not a proof of maximality)."""
import numpy as np
from verify_misc import twisted_matrix

m = 1039
for k in (64, 264, 304, 96, 1, 2):
    r = float(np.max(np.abs(np.linalg.eigvals(twisted_matrix(m, k))))) / 2
    print(f"m = {m}, k = {k}: rho/2 = {r:.6f}")
