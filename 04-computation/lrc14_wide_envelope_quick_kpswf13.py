#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick wide sup-envelope C_wide pin (kind-pasteur S36, Thread 3). Flushed per line."""
import sys
sys.path.insert(0, "04-computation")
from math import pi
from lrc14_spectrum_intersection_sum_kpswf12 import (
    meas, complement, cover_set, safe_set, fourier_num_of_arcs,
)

def chat(arcs, n):
    if n == 0:
        return complex(float(meas(arcs)), 0.0)
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)

def spec(P, E, N=2500):
    gp = safe_set(P); covc = complement(cover_set(E))
    return sum(2.0 * (chat(gp, n) * chat(covc, n).conjugate()).real for n in range(1, N + 1))

def p(*a):
    print(*a, flush=True)

cases = [([1, 3, 4, 5], list(range(9)), "argmin k9"),
         ([1, 2, 3], list(range(10)), "k10 |P|=3"),
         ([2, 3, 5], list(range(10)), "coprime-P")]
allMD = []
for P, E0, lab in cases:
    p(f"CASE {lab}  P={P}")
    sp = {}
    for d in range(1, 49):
        sp[d] = abs(spec(P, [d * e for e in E0]))
    for D in [1, 2, 3, 5, 8, 13, 21, 34]:
        dom = [d for d in sp if d >= D]
        am = max(dom, key=lambda d: sp[d]); M = sp[am]
        allMD.append(M * D)
        p(f"   D={D:3d}  M(D)=sup_{{d>=D}}|SPEC|={M:.5f}   M(D)*D={M*D:.5f}   argmax_d={am}")
p(f"CWIDE_envelope_sup = sup_D M(D)*D = {max(allMD):.5f}")
p("=> |SPEC(dE)| <= CWIDE/d  =>  R'(dE) >= 1 - CWIDE/(meas(GP)(1-p0)*d); wide eta(V*) ~ CWIDE*span0/V* -> 0.")
p("DONE")
