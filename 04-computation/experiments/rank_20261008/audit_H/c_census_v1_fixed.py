#!/usr/bin/env python3
"""audit H: the census-v1 FIXED k = 4 certificate for Z_5 (1,1,28,11,7) (form from block_census_fixed_k4.out), re-verified exactly.
It contradicts THM-4611 (5) "the twelve rank-3 maps need adaptive lengths 3-5" for this map.  Output: c_census_v1_fixed.out"""
import c_certificates as C
C.CERTS['Z5_1_1_28_11_7_k4'] = (5, [1, 1, 28, 11, 7], 4, [[23, -4, -8], [-4, 20, -4], [-8, -4, 24]])
C.run('Z5_1_1_28_11_7_k4')
