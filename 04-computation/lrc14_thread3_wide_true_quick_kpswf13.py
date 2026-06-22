#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""C_wide on the TRUE GOOD (maxgap>1/7), bounded dilation range (kind-pasteur S37, Thread 3). Flushed."""
import sys
sys.path.insert(0, "04-computation")
from math import pi
from lrc14_spectrum_intersection_sum_kpswf12 import meas, intersect, safe_set, fourier_num_of_arcs
from lrc14_true_maxgap_good_kpswf13 import good_true

def chat(arcs, n):
    if n == 0:
        return complex(float(meas(arcs)), 0.0)
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)

def spec_true(P, E, good, N=2000):
    gp = safe_set(P)
    return sum(2.0 * (chat(gp, n) * chat(good, n).conjugate()).real for n in range(1, N + 1))

def Rt(P, E, good):
    gp = safe_set(P); base = meas(gp) * meas(good)
    return None if base == 0 else float(meas(intersect(gp, good)) / base)

def p(*a):
    print(*a, flush=True)

cases = [([1, 3, 4, 5], list(range(9)), "k9"),
         ([1, 2, 3], list(range(10)), "k10|P|3"),
         ([2, 3, 5], list(range(10)), "coprimeP")]
allMD = []
for P, E0, lab in cases:
    p(f"CASE {lab} P={P}")
    sp = {}; rr = {}
    for d in range(1, 25):
        E = [d * e for e in E0]
        g = good_true(E)
        sp[d] = abs(spec_true(P, E, g, N=2000))
        rr[d] = Rt(P, E, g)
        p(f"   d={d:3d} span={max(E):4d}  |SPEC_true|={sp[d]:.5f}  R_true={rr[d]:.5f}")
    for D in [1, 2, 3, 5, 8, 13]:
        dom = [d for d in sp if d >= D]
        am = max(dom, key=lambda d: sp[d]); M = sp[am]
        allMD.append(M * D)
        p(f"     env D={D:3d}: M(D)={M:.5f}  M(D)*D={M*D:.5f}")
p(f"CWIDE_TRUE = sup_D M(D)*D = {max(allMD):.5f}")
p("DONE")
