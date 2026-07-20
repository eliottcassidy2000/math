#!/usr/bin/env python3
"""powersum_triangle_odd_sector_klein_S330.py -- klein-2026-07-20-S330.
(1) Triangle law T(n,k) = PS_{k-1}(n-k+1) + [k>=4][n-k>=2]: 28/28 exact; row 8
    prediction {8,28,91,226,355,277,65,1}.
(2) Series: s1 = H_n exact; s2 vs row-harmonic sum T(n,k)/k (exact n<=3 under
    29/6-reading; open at n>=4).
(3) Faulhaber: PS_3 = T^2, PS_5 = T^2(4T-1)/3, PS_7 = T^2(6T^2-4T+1)/3.
(4) Tournaments: c3 = C(n,3) - sum C(s_i,2) exhaustive at n = 4, 5 (64 + 1024).
See the S330 reflection for the full weave; code as run in-session."""
